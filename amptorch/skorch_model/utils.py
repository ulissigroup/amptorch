import sys
from skorch.utils import to_numpy
import torch
from torch.nn import MSELoss, L1Loss
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
    return energy_rmse


def forces_score(net, X, y):
    mse_loss = MSELoss(reduction="sum")
    outputs = net.forward(X)
    energy, forces = outputs
    if forces.nelement() == 0:
        raise Exception("Force training disabled. Disable force scoring!")
    device = forces.device
    if not hasattr(X, "scalings"):
        X = X.dataset
    num_atoms = torch.FloatTensor(np.concatenate(y[1::3])).reshape(-1, 1).to(device)
    force_target = torch.tensor(np.concatenate(y[2::3])).to(device)
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
    return force_rmse


def energy_mad(net, X, y):
    l1_loss = L1Loss(reduction="none")
    energy_pred, _ = net.forward(X)
    device = energy_pred.device
    if not hasattr(X, "scalings"):
        X = X.dataset
    num_atoms = torch.FloatTensor(np.concatenate(y[1::3])).reshape(-1, 1).to(device)
    energy_target = torch.tensor(np.concatenate(y[0::3])).to(device).reshape(-1, 1)
    if X.scaling_scheme is not "log":
        sd_scaling = X.scalings[0]
        mean_scaling = X.scalings[1]
        raw_preds = (energy_pred * sd_scaling) + mean_scaling
        raw_preds_per_atom = torch.div(raw_preds, num_atoms)
        raw_targets = (energy_target * sd_scaling) + mean_scaling
        target_per_atom = torch.div(raw_targets, num_atoms)
        energy_loss = l1_loss(raw_preds_per_atom, target_per_atom)
        energy_mad_loss = torch.median(energy_loss)
    else:
        raw_preds = torch.exp(energy_pred) - 1
        raw_preds_per_atom = torch.div(raw_preds, num_atoms)
        raw_targets = torch.exp(energy_target) - 1
        target_per_atom = torch.div(raw_targets, num_atoms)
        energy_loss = l1_loss(raw_preds_per_atom, target_per_atom)
        energy_mad_loss = torch.median(energy_loss)
    return energy_mad_loss


def forces_mad(net, X, y):
    l1_loss = L1Loss(reduction="none")
    outputs = net.forward(X)
    energy, forces = outputs
    if forces.nelement() == 0:
        raise Exception("Force training disabled. Disable force scoring!")
    device = forces.device
    if not hasattr(X, "scalings"):
        X = X.dataset
    num_atoms = torch.FloatTensor(np.concatenate(y[1::3])).reshape(-1, 1).to(device)
    dataset_size = len(num_atoms)
    force_target = torch.tensor(np.concatenate(y[2::3])).to(device)
    energy_target = torch.tensor(np.concatenate(y[0::3])).to(device).reshape(-1, 1)
    if X.scaling_scheme is not "log":
        sd_scaling = X.scalings[0]
        force_pred = forces * sd_scaling
        device = force_pred.device
        raw_force_target = force_target * sd_scaling
        num_atoms_force = torch.cat(
            [idx.repeat(int(idx)) for idx in num_atoms]
        ).reshape(-1, 1)
        force_pred_per_atom = torch.div(force_pred, num_atoms_force)
        force_targets_per_atom = torch.div(raw_force_target, num_atoms_force)
        force_loss = l1_loss(force_pred_per_atom, force_targets_per_atom)
        force_loss = torch.sum(l1_loss(force_pred_per_atom, force_targets_per_atom), 1)
        sum_idx = 0
        force_loss_per_image = torch.zeros(dataset_size, 1)
        for idx, atoms in enumerate(num_atoms):
            atoms = int(atoms)
            force_loss_per_image[idx] = torch.sum(force_loss[sum_idx : sum_idx + atoms])
            sum_idx += atoms
        force_loss_per_image /= 3
        force_mad_loss = torch.median(force_loss_per_image)
    else:
        # TODO median force loss
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
        raw_force_target = force_target * (energy_target + 1)
        num_atoms_force = torch.cat([idx.repeat(int(idx)) for idx in num_atoms])
        force_pred_per_atom = torch.div(force_pred, num_atoms_force)
        force_targets_per_atom = torch.div(raw_force_target, num_atoms_force)
        force_mse = l1_loss(force_pred_per_atom, force_targets_per_atom)
        force_mse /= 3
        force_mad_loss = torch.median(force_loss)
    return force_mad_loss


def make_force_header(log):
    header = "%5s %12s %12s %12s %7s"
    log(header % ("Epoch", "EnergyRMSE", "ForceRMSE", "TrainLoss", "Dur"))
    log(header % ("=" * 5, "=" * 12, "=" * 12, "=" * 12, "=" * 7))


def make_energy_header(log):
    header = "%5s %12s %12s %7s"
    log(header % ("Epoch", "EnergyRMSE", "TrainLoss", "Dur"))
    log(header % ("=" * 5, "=" * 12, "=" * 12, "=" * 7))


def make_val_force_header(log):
    header = "%5s %12s %12s %12s %12s %7s"
    log(header % ("Epoch", "EnergyRMSE", "ForceRMSE", "TrainLoss", "ValidLoss", "Dur"))
    log(header % ("=" * 5, "=" * 12, "=" * 12, "=" * 12, "=" * 12, "=" * 7))


def make_val_energy_header(log):
    header = "%5s %12s %12s %12s %7s"
    log(header % ("Epoch", "EnergyRMSE", "TrainLoss", "ValidLoss", "Dur"))
    log(header % ("=" * 5, "=" * 12, "=" * 12, "=" * 12, "=" * 7))


def log_results(model, log):
    log("Training initiated...")
    if model.train_split != 0:
        if model.criterion__force_coefficient != 0:
            make_val_force_header(log)
            for epoch, ermse, frmse, tloss, vloss, dur in model.history[
                :,
                (
                    "epoch",
                    "energy_score",
                    "forces_score",
                    "train_loss",
                    "valid_loss",
                    "dur",
                ),
            ]:
                log(
                    "%5i %12.4f %12.4f %12.4f %12.4f %7.4f"
                    % (epoch, ermse, frmse, tloss, vloss, dur)
                )
        else:
            make_val_energy_header(log)
            for epoch, ermse, tloss, vloss, dur in model.history[
                :, ("epoch", "energy_score", "train_loss", "valid_loss", "dur")
            ]:
                log(
                    "%5i %12.4f %12.4f %12.4f %7.4f" % (epoch, ermse, tloss, vloss, dur)
                )
    else:
        if model.criterion__force_coefficient != 0:
            make_force_header(log)
            for epoch, ermse, frmse, tloss, dur in model.history[
                :, ("epoch", "energy_score", "forces_score", "train_loss", "dur")
            ]:
                log(
                    "%5i %12.4f %12.4f %12.4f %7.4f" % (epoch, ermse, frmse, tloss, dur)
                )
        else:
            make_energy_header(log)
            for epoch, ermse, tloss, dur in model.history[
                :, ("epoch", "energy_score", "train_loss", "dur")
            ]:
                log("%5i %12.4f %12.4f %7.4f" % (epoch, ermse, tloss, dur))
    log("...Training Complete!\n")
