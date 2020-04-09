import sys
import skorch
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
    dataset_size = len(energy_pred)
    energy_targets_per_atom = torch.tensor(np.concatenate(y[0::3])).to(device).reshape(-1, 1)
    energy_preds_per_atom = torch.div(energy_pred, num_atoms)
    energy_loss = mse_loss(energy_preds_per_atom, energy_targets_per_atom)
    energy_loss /= dataset_size
    energy_rmse = torch.sqrt(energy_loss)
    return energy_rmse

def forces_score(net, X, y):
    mse_loss = MSELoss(reduction="sum")
    _, force_pred = net.forward(X)
    if force_pred.nelement() == 0:
        raise Exception("Force training disabled. Disable force scoring!")
    device = force_pred.device
    if not hasattr(X, "scalings"):
        X = X.dataset
    num_atoms = torch.FloatTensor(np.concatenate(y[1::3])).reshape(-1, 1).to(device)
    force_targets_per_atom = torch.tensor(np.concatenate(y[2::3])).to(device)
    device = force_pred.device
    dataset_size = len(num_atoms)
    num_atoms_extended = torch.cat([idx.repeat(int(idx)) for idx in num_atoms])
    num_atoms_extended = torch.sqrt(num_atoms_extended).reshape(-1, 1)
    force_pred_per_atom = torch.div(force_pred, num_atoms_extended)
    force_targets_per_atom = force_targets_per_atom*num_atoms_extended
    force_mse = mse_loss(force_pred_per_atom, force_targets_per_atom)
    force_mse /= 3 * dataset_size
    force_rmse = torch.sqrt(force_mse)
    return force_rmse

class train_end_load_best_loss(skorch.callbacks.base.Callback):
    def __init__(self, filename):
        self.filename = filename

    def on_train_end(self, net, X, y):
        net.load_params("./results/checkpoints/{}_params.pt".format(self.filename))

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


