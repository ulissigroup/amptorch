import sys
from ase import Atoms
from ase.calculators.lj import LennardJones as LJ
from ase.calculators.emt import EMT
from amptorch.gaussian import SNN_Gaussian
from amptorch.data_preprocess import AtomsDataset, collate_amp
from amptorch.core import AMPTorch
from amptorch.model import CustomMSELoss, BPNN
from amptorch.delta_models.morse import morse_potential
import numpy as np
import torch
from torch import optim
from skorch import NeuralNetRegressor
from skorch.callbacks import Checkpoint, EpochScoring
from amptorch.skorch_model.utils import forces_score, target_extractor, energy_score
import ase


def test_skorch_delta():
    from amptorch.skorch_model import AMP

    cp = Checkpoint(monitor="valid_loss_best", fn_prefix="valid_best_")

    distances = np.linspace(2, 5, 10)
    label = "skorch_example"
    images = []
    energies = []
    forces = []
    for l in distances:
        image = Atoms(
            "CuCO",
            [
                (-l * np.sin(0.65), l * np.cos(0.65), 0),
                (0, 0, 0),
                (l * np.sin(0.65), l * np.cos(0.65), 0),
            ],
        )
        image.set_cell([10, 10, 10])
        image.wrap(pbc=True)
        image.set_calculator(EMT())
        images.append(image)
        energies.append(image.get_potential_energy())
        forces.append(image.get_forces())

    image = Atoms("CuC", [(-1, 1, 0), (1, 1, 0)])
    image.set_cell([10, 10, 10])
    image.wrap(pbc=True)
    image.set_calculator(EMT())
    images.append(image)
    energies.append(image.get_potential_energy())
    forces.append(image.get_forces())

    energies = np.array(energies)
    forces = np.concatenate(np.array(forces))

    Gs = {}
    Gs["G2_etas"] = np.logspace(np.log10(0.05), np.log10(5.0), num=2)
    Gs["G2_rs_s"] = [0] * 2
    Gs["G4_etas"] = [0.005]
    Gs["G4_zetas"] = [1.0]
    Gs["G4_gammas"] = [+1.0, -1]
    Gs["cutoff"] = 6.5

    morse_model = morse_potential(images, Gs["cutoff"], label)
    morse_energies, morse_forces, num_atoms = morse_model.morse_pred(images)
    morse_data = [morse_energies, morse_forces, num_atoms, morse_model]

    forcetraining = True
    training_data = AtomsDataset(
        images,
        SNN_Gaussian,
        Gs,
        forcetraining=forcetraining,
        label=label,
        cores=1,
        delta_data=morse_data,
    )
    batch_size = len(training_data)
    unique_atoms = training_data.elements
    fp_length = training_data.fp_length
    device = "cpu"

    net = NeuralNetRegressor(
        module=BPNN(
            unique_atoms, [fp_length, 2, 2], device, forcetraining=forcetraining
        ),
        criterion=CustomMSELoss,
        criterion__force_coefficient=0.3,
        optimizer=torch.optim.LBFGS,
        optimizer__line_search_fn="strong_wolfe",
        lr=1e-2,
        batch_size=batch_size,
        max_epochs=100,
        iterator_train__collate_fn=collate_amp,
        iterator_train__shuffle=False,
        iterator_valid__collate_fn=collate_amp,
        device=device,
        train_split=0,
        verbose=0,
        callbacks=[
            EpochScoring(
                forces_score,
                on_train=True,
                use_caching=True,
                target_extractor=target_extractor,
            ),
            EpochScoring(
                energy_score,
                on_train=True,
                use_caching=True,
                target_extractor=target_extractor,
            ),
        ],
    )
    calc = AMP(training_data, net, "test")
    calc.train(overwrite=True)
    num_of_atoms = 3
    calculated_energies = np.array(
        [calc.get_potential_energy(image) for idx, image in enumerate(images)]
    )
    energy_rmse = np.sqrt(
        (((calculated_energies - energies) / num_of_atoms) ** 2).sum() / len(images)
    )
    last_energy_score = net.history[-1]["energy_score"]
    assert round(energy_rmse, 4) == round(
        last_energy_score, 4
    ), "Energy errors incorrect!"
    last_forces_score = net.history[-1]["forces_score"]

    calculated_forces = np.concatenate(
        np.array([calc.get_forces(image) for image in images])
    )
    force_rmse = np.sqrt(
        (((calculated_forces - forces)) ** 2).sum() / (3 * num_of_atoms * len(images))
    )
    assert round(force_rmse, 4) == round(
        last_forces_score, 4
    ), "Force errors incorrect!"
