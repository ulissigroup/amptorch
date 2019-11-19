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
from amptorch.gaussian import Gaussian
from amptorch.skorch import AMP
from amptorch.model import FullNN, CustomLoss
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


def target_extractor(y):
    return (
        (to_numpy(y[0]), to_numpy(y[1]))
        if len(y) == 2
        else (to_numpy(y[0]), to_numpy(y[1]), to_numpy(y[2]))
    )


def energy_score(net, X, y):
    mse_loss = MSELoss(reduction="sum")
    energy_pred = net.infer()[0]
    device = energy_pred.device
    energy_target = torch.tensor(y[0]).to(device)
    num_atoms = torch.tensor(y[1]).to(device)
    dataset_size = len(energy_pred)
    sd_scaling = float(X[0][-1][0])
    mean_scaling = float(X[0][-1][1])
    raw_preds = (energy_pred * sd_scaling) + mean_scaling
    raw_preds_per_atom = torch.div(raw_preds, num_atoms)
    raw_targets = (energy_target * sd_scaling) + mean_scaling
    target_per_atom = torch.div(raw_targets, num_atoms)
    energy_loss = mse_loss(raw_preds_per_atom, target_per_atom)
    energy_loss /= dataset_size
    energy_rmse = torch.sqrt(energy_loss)
    return energy_rmse


def forces_score(net, X, y):
    mse_loss = MSELoss(reduction="sum")
    sd_scaling = float(X[0][-1][0])
    force_pred = net.infer()[1] * sd_scaling
    device = force_pred.device
    num_atoms = torch.tensor(y[1]).to(device)
    force_target = torch.tensor(y[-1], device=device)
    dataset_size = len(num_atoms)
    raw_force_target = force_target * sd_scaling
    num_atoms_force = torch.cat([idx.repeat(int(idx)) for idx in num_atoms])
    num_atoms_force = torch.sqrt(num_atoms_force).reshape(len(num_atoms_force), 1)
    force_pred_per_atom = torch.div(force_pred, num_atoms_force)
    force_targets_per_atom = torch.div(raw_force_target, num_atoms_force)
    force_mse = mse_loss(force_pred_per_atom, force_targets_per_atom)
    force_mse /= 3 * dataset_size
    force_rmse = torch.sqrt(force_mse)
    return force_rmse


class train_end_load_best_valid_loss(skorch.callbacks.base.Callback):
    def on_train_end(self, net, X, y):
        net.load_params("valid_best_params.pt")


def test_skorch():
    LR_schedule = LRScheduler("CosineAnnealingLR", T_max=5)
    cp = Checkpoint(monitor="valid_loss_best", fn_prefix="valid_best_")
    load_best_valid_loss = train_end_load_best_valid_loss()

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
        Gaussian,
        Gs,
        forcetraining=forcetraining,
        label=label,
        cores=1,
        lj_data=None,
    )
    batch_size = len(training_data)
    unique_atoms = training_data.elements
    fp_length = training_data.fp_length
    device = "cpu"

    net = NeuralNetRegressor(
        module=FullNN(
            unique_atoms, [fp_length, 2, 2], device, forcetraining=forcetraining
        ),
        criterion=CustomLoss,
        criterion__force_coefficient=0.3,
        optimizer=torch.optim.LBFGS,
        optimizer__line_search_fn="strong_wolfe",
        lr=1e-2,
        batch_size=batch_size,
        max_epochs=60,
        iterator_train__collate_fn=collate_amp,
        iterator_train__shuffle=False,
        device=device,
        train_split=0,
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
    assert force_rmse <= 0.005, "Force training convergence not met!"
    assert energy_rmse <= 0.005, "Energy training convergence not met!"
