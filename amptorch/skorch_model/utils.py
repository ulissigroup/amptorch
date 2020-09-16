import skorch
import torch
from skorch.utils import to_numpy
from torch.nn import MSELoss


def target_extractor(y):
    return (
        (to_numpy(y[0]), to_numpy(y[1]))
        if len(y) == 2
        else (to_numpy(y[0]), to_numpy(y[1]), to_numpy(y[2]))
    )


def to_tensor(X, device, accept_sparse=False):
    return X


def energy_score(net, X, y):
    mse_loss = MSELoss()
    energy_pred, _ = net.forward(X)
    energy_target = torch.FloatTensor(y[0])
    energy_loss = mse_loss(energy_pred, energy_target)
    return energy_loss


def forces_score(net, X, y):
    mse_loss = MSELoss()
    _, force_pred = net.forward(X)
    force_target = torch.FloatTensor(y[1])
    force_loss = mse_loss(force_pred, force_target)
    return force_loss


class train_end_load_best_loss(skorch.callbacks.base.Callback):
    def __init__(self, filename):
        self.filename = filename

    def on_train_end(self, net, X, y):
        net.load_params("./checkpoints/{}_params.pt".format(self.filename))


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
