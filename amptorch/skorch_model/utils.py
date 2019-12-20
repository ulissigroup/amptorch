import sys
from skorch.utils import to_numpy
import torch
from torch.nn import MSELoss
import numpy as np


def target_extractor(y):
    return (
        (to_numpy(y[0]), to_numpy(y[1]))
        if len(y) == 2
        else (to_numpy(y[0]), to_numpy(y[1]), to_numpy(y[2]))
    )

def energy_score(net, X, y):
    mse_loss = MSELoss(reduction="sum")
    energy_pred, _ = net.forward(X)
    device = energy_pred.device
    if not hasattr(X, "scalings"):
        X = X.dataset
    num_atoms = torch.FloatTensor(np.concatenate(y[1::3])).reshape(-1, 1).to(device)
    energy_target = torch.tensor(np.concatenate(y[0::3])).to(device).reshape(-1, 1)
    if X.scaling_scheme is not "log":
        sd_scaling = X.scalings[0]
        mean_scaling = X.scalings[1]
        dataset_size = len(energy_pred)
        raw_preds = (energy_pred * sd_scaling) + mean_scaling
        raw_preds_per_atom = torch.div(raw_preds, num_atoms)
        raw_targets = (energy_target * sd_scaling) + mean_scaling
        target_per_atom = torch.div(raw_targets, num_atoms)
        energy_loss = mse_loss(raw_preds_per_atom, target_per_atom)
        energy_loss /= dataset_size
        energy_rmse = torch.sqrt(energy_loss)
    else:
        dataset_size = len(energy_pred)
        raw_preds = torch.exp(energy_pred) - 1
        raw_preds_per_atom = torch.div(raw_preds, num_atoms)
        raw_targets = torch.exp(energy_target) - 1
        target_per_atom = torch.div(raw_targets, num_atoms)
        energy_loss = mse_loss(raw_preds_per_atom, target_per_atom)
        energy_loss /= dataset_size
        energy_rmse = torch.sqrt(energy_loss)
    return energy_rmse

def forces_score(net, X, y):
    mse_loss = MSELoss(reduction="sum")
    outputs = net.forward(X)
    energy, forces = outputs
    if forces.nelement() == 0:
        raise Exception('Force training disabled. Disable force scoring!')
    device = forces.device
    if not hasattr(X, "scalings"):
        X = X.dataset
    num_atoms = torch.FloatTensor(np.concatenate(y[1::3])).reshape(-1, 1).to(device)
    force_target = torch.tensor(np.concatenate(y[2::3])).to(device)
    energy_target = torch.tensor(np.concatenate(y[0::3])).to(device).reshape(-1, 1)
    if X.scaling_scheme is not "log":
        sd_scaling = X.scalings[0]
        force_pred = forces * sd_scaling
        device = force_pred.device
        dataset_size = len(num_atoms)
        raw_force_target = force_target * sd_scaling
        num_atoms_force = torch.cat([idx.repeat(int(idx)) for idx in num_atoms])
        num_atoms_force = torch.sqrt(num_atoms_force).reshape(len(num_atoms_force), 1)
        force_pred_per_atom = torch.div(force_pred, num_atoms_force)
        force_targets_per_atom = torch.div(raw_force_target, num_atoms_force)
        force_mse = mse_loss(force_pred_per_atom, force_targets_per_atom)
        force_mse /= 3 * dataset_size
        force_rmse = torch.sqrt(force_mse)
    else:
        energy = torch.exp(energy) - 1
        energy_target = torch.exp(energy_target) - 1
        energy = torch.cat(
            [k.repeat(int(num_atoms[idx])) for idx, k in enumerate(energy)]
        ).reshape(-1, 1)
        energy_target = torch.cat(
            [k.repeat(int(num_atoms[idx])) for idx, k in enumerate(energy_target)]
        ).reshape(-1, 1)
        force_pred = forces * (energy + 1)
        device = force_pred.device
        dataset_size = len(num_atoms)
        raw_force_target = force_target * (energy_target + 1)
        num_atoms_force = torch.cat([idx.repeat(int(idx)) for idx in num_atoms])
        num_atoms_force = torch.sqrt(num_atoms_force).reshape(len(num_atoms_force), 1)
        force_pred_per_atom = torch.div(force_pred, num_atoms_force)
        force_targets_per_atom = torch.div(raw_force_target, num_atoms_force)
        force_mse = mse_loss(force_pred_per_atom, force_targets_per_atom)
        force_mse /= 3 * dataset_size
        force_rmse = torch.sqrt(force_mse)
    return force_rmse
