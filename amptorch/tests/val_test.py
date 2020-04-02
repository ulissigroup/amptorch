import torch
from torch import optim
from skorch import NeuralNetRegressor
from skorch.dataset import CVSplit
from skorch.callbacks import Checkpoint, EpochScoring
from amptorch.gaussian import SNN_Gaussian
from amptorch.skorch_model import AMP as AMP_skorch
from amptorch.skorch_model.utils import target_extractor, energy_score, forces_score
from amptorch.model import FullNN, CustomMSELoss
from amptorch.data_preprocess import AtomsDataset, collate_amp
import numpy as np
from ase import Atoms
from ase.calculators.emt import EMT
from ase.io import read

def test_skorch_val():
    distances = np.linspace(2, 5, 100)
    label = "example"
    images = []
    energies = []
    forces = []
    for l in distances:
        image = Atoms(
            "CuCCu",
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

    energies = np.array(energies)

    Gs = {}
    Gs["G2_etas"] = np.logspace(np.log10(0.05), np.log10(5.0), num=2)
    Gs["G2_rs_s"] = [0] * 2
    Gs["G4_etas"] = [0.005]
    Gs["G4_zetas"] = [1.0]
    Gs["G4_gammas"] = [+1.0, -1]
    Gs["cutoff"] = 6.5

    forcetraining = True
    training_data = AtomsDataset(
        images,
        SNN_Gaussian,
        Gs,
        forcetraining=forcetraining,
        label=label,
        cores=1,
        delta_data=None,
    )
    batch_size = len(training_data)
    unique_atoms = training_data.elements
    fp_length = training_data.fp_length
    device = "cpu"

    net = NeuralNetRegressor(
        module=FullNN(
            unique_atoms, [fp_length, 2, 2], device, forcetraining=forcetraining
        ),
        criterion=CustomMSELoss,
        criterion__force_coefficient=0.3,
        optimizer=torch.optim.LBFGS,
        optimizer__line_search_fn="strong_wolfe",
        lr=1e-2,
        batch_size=10,
        max_epochs=20,
        iterator_train__collate_fn=collate_amp,
        iterator_train__shuffle=True,
        iterator_valid__collate_fn=collate_amp,
        device=device,
        train_split=CVSplit(0.1, random_state=1),
        verbose=0,
        callbacks=[
            EpochScoring(
                forces_score,
                on_train=False,
                use_caching=True,
                target_extractor=target_extractor,
            ),
            EpochScoring(
                energy_score,
                on_train=False,
                use_caching=True,
                target_extractor=target_extractor,
            ),
        ],
    )
    val_indices = [80, 84, 33, 81, 93, 17, 36, 82, 69, 65]
    val_images = [images[idx] for idx in val_indices]
    val_energies = energies[val_indices]
    val_forces = np.concatenate(np.array([forces[idx] for idx in val_indices]))
    calc = AMP_skorch(training_data, net, "test")
    calc.train(overwrite=True)
    num_of_atoms = 3

    last_energy_score = net.history[-1]["energy_score"]
    last_forces_score = net.history[-1]["forces_score"]

    calculated_energies = np.array(
        [calc.get_potential_energy(image) for image in val_images]
    )
    energy_rmse = np.sqrt(
        (((calculated_energies - val_energies) / num_of_atoms) ** 2).sum()
        / len(val_images)
    )
    assert round(energy_rmse, 5) == round(
        last_energy_score, 5
    ), "Validation energy errors incorrect!"
    calculated_forces = np.concatenate(
        np.array([calc.get_forces(image) for image in val_images])
    )
    force_rmse = np.sqrt(
        (((calculated_forces - val_forces)) ** 2).sum()
        / (3 * num_of_atoms * len(val_images))
    )
    assert round(force_rmse, 5) == round(
        last_forces_score, 5
    ), "Validation force errors incorrect!"

def test_energy_only_skorch_val():
    distances = np.linspace(2, 5, 100)
    label = "example"
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

    energies = np.array(energies)

    Gs = {}
    Gs["G2_etas"] = np.logspace(np.log10(0.05), np.log10(5.0), num=2)
    Gs["G2_rs_s"] = [0] * 2
    Gs["G4_etas"] = [0.005]
    Gs["G4_zetas"] = [1.0]
    Gs["G4_gammas"] = [+1.0, -1]
    Gs["cutoff"] = 6.5

    forcetraining = False
    training_data = AtomsDataset(
        images,
        SNN_Gaussian,
        Gs,
        forcetraining=forcetraining,
        label=label,
        cores=1,
        delta_data=None,
    )
    batch_size = len(training_data)
    unique_atoms = training_data.elements
    fp_length = training_data.fp_length
    device = "cpu"

    net = NeuralNetRegressor(
        module=FullNN(
            unique_atoms, [fp_length, 2, 2], device, forcetraining=forcetraining
        ),
        criterion=CustomMSELoss,
        criterion__force_coefficient=0,
        optimizer=torch.optim.LBFGS,
        optimizer__line_search_fn="strong_wolfe",
        lr=1e-2,
        batch_size=5,
        max_epochs=20,
        iterator_train__collate_fn=collate_amp,
        iterator_train__shuffle=True,
        iterator_valid__collate_fn=collate_amp,
        device=device,
        train_split=CVSplit(cv=0.1, random_state=1),
        verbose=0,
        callbacks=[
            EpochScoring(
                energy_score,
                on_train=False,
                use_caching=True,
                target_extractor=target_extractor,
            ),
        ],
    )
    val_indices = [80, 84, 33, 81, 93, 17, 36, 82, 69, 65]
    val_images = [images[idx] for idx in val_indices]
    val_energies = energies[val_indices]
    calc = AMP_skorch(training_data, net, "test")
    calc.train(overwrite=True)
    num_of_atoms = 3

    last_energy_score = net.history[-1]["energy_score"]

    calculated_energies = np.array(
        [calc.get_potential_energy(image) for image in val_images]
    )
    energy_rmse = np.sqrt(
        (((calculated_energies - val_energies) / num_of_atoms) ** 2).sum()
        / len(val_images)
    )
    assert round(energy_rmse, 5) == round(
        last_energy_score, 5
    ), "Validation energy errors incorrect!"
