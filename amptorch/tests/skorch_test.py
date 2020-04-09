import ase
import sys
import torch
import time
from torch.nn import MSELoss
from skorch import NeuralNetRegressor
from skorch.dataset import CVSplit
from skorch.callbacks import Checkpoint, EpochScoring
from skorch.callbacks.lr_scheduler import LRScheduler
import skorch.callbacks.base
from amptorch.gaussian import SNN_Gaussian
from amptorch.skorch_model import AMP
from amptorch.skorch_model.utils import (
    target_extractor,
    energy_score,
    forces_score,
)
from amptorch.model import FullNN, CustomMSELoss
from amptorch.data_preprocess import (
    AtomsDataset,
    factorize_data,
    collate_amp,
    TestDataset,
)
from torch.utils.data import DataLoader
from torch.nn import init
from skorch.utils import to_numpy
import numpy as np
from sklearn.pipeline import Pipeline
from ase import Atoms
from ase.calculators.emt import EMT
from ase.io import read


def test_skorch():
    distances = np.linspace(2, 5, 100)
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

    energies = np.array(energies)
    forces = np.concatenate(np.array(forces))
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
        lr=1,
        batch_size=100,
        max_epochs=150,
        iterator_train__collate_fn=collate_amp,
        iterator_train__shuffle=True,
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
        [calc.get_potential_energy(image) for image in images]
    )
    energy_rmse = np.sqrt(
        (((calculated_energies - energies) / num_of_atoms) ** 2).sum() / len(images)
    )

    calculated_forces = np.concatenate(
        np.array([calc.get_forces(image) for image in images])
    )
    force_rmse = np.sqrt(
        (((calculated_forces - forces)) ** 2).sum() / (3 * num_of_atoms * len(images))
    )
    l1_force = np.sum(np.abs(calculated_forces - forces) / num_of_atoms, 1)
    idx = 0
    force_loss_image = np.zeros((len(calculated_energies), 1))
    for i in range(len(calculated_energies)):
        force_loss_image[i] = np.sum(l1_force[idx: idx + 3])
        idx += 3
    force_loss_image /= 3

    reported_energy_score = net.history[-1]["energy_score"]
    reported_forces_score = net.history[-1]["forces_score"]
    assert force_rmse <= 0.005, "Force training convergence not met!"
    assert energy_rmse <= 0.005, "Energy training convergence not met!"
    assert round(reported_energy_score, 4) == round(
        energy_rmse, 4
    ), "Shuffled reported energy scores incorrect!"
    assert round(reported_forces_score, 4) == round(
        force_rmse, 4
    ), "Shuffled reported forces score incorrect!"

def test_e_only_skorch():
    distances = np.linspace(2, 5, 100)
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

    energies = np.array(energies)
    forces = np.concatenate(np.array(forces))
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
        batch_size=batch_size,
        max_epochs=100,
        iterator_train__collate_fn=collate_amp,
        iterator_train__shuffle=True,
        iterator_valid__collate_fn=collate_amp,
        device=device,
        train_split=0,
        verbose=0,
        callbacks=[
            EpochScoring(
                energy_score,
                on_train=True,
                use_caching=True,
                target_extractor=target_extractor,
            )
        ],
    )
    calc = AMP(training_data, net, "test")
    calc.train(overwrite=True)
    num_of_atoms = 3
    calculated_energies = np.array(
        [calc.get_potential_energy(image) for image in images]
    )
    energy_rmse = np.sqrt(
        (((calculated_energies - energies) / num_of_atoms) ** 2).sum() / len(images)
    )

    reported_energy_score = net.history[-1]["energy_score"]
    assert energy_rmse <= 0.005, "Energy training convergence not met!"
    assert round(energy_rmse, 4) == round(
        reported_energy_score, 4
    ), "Shuffled energy only energy score incorrect!"
