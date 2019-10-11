import ase
import sys
import torch
from torch.nn import MSELoss
from skorch import NeuralNetRegressor
from skorch.dataset import CVSplit
from skorch.callbacks import Checkpoint, EpochScoring
import skorch.callbacks.base
from amp.descriptor.gaussian import Gaussian
from amptorch.skorch.model_skorch import FullNN, CustomLoss
from amptorch.skorch.skorch_data import AtomsDataset, factorize_data, collate_amp
from skorch.utils import to_numpy


def target_extractor(y):
    return(to_numpy(y[0]), to_numpy(y[1])) if len(y) == 2 else (to_numpy(y[0]),
            to_numpy(y[1]), to_numpy(y[2]))

def energy_score(net, X, y):
    mse_loss = MSELoss(reduction="sum")
    energy_pred = net.infer()[0]
    device = energy_pred.device
    energy_target = torch.tensor(y[0]).to(device)
    num_atoms = torch.tensor(y[1]).to(device)
    dataset_size = len(energy_pred)
    sd_scaling = scalings[0]
    mean_scaling = scalings[1]
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
    sd_scaling = scalings[0]
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
    force_mse /= (3*dataset_size)
    force_rmse = torch.sqrt(force_mse)
    return force_rmse

forcetraining = True
data = AtomsDataset("../../datasets/water/water.extxyz", descriptor=Gaussian(),
        cores=1, forcetraining=forcetraining)
scalings = data.scalings
unique_atoms = data.unique_atoms
fp_length = data.fp_length

device = 'cpu'
# device = 'cuda:0'

net = NeuralNetRegressor(
    module=FullNN(unique_atoms, [fp_length, 5, 5], device,
        forcetraining=forcetraining),
    criterion=CustomLoss,
    criterion__force_coefficient=0.3,
    optimizer=torch.optim.LBFGS,
    lr=1,
    batch_size=400,
    max_epochs=50,
    iterator_train__collate_fn=collate_amp,
    iterator_valid__collate_fn=collate_amp,
    device=device,
    train_split=None,
    callbacks=[EpochScoring(forces_score, on_train=True, use_caching=True,
        target_extractor=target_extractor), EpochScoring(energy_score,
            on_train=True, use_caching=True, target_extractor=target_extractor)],
)

net.fit(data, None)
