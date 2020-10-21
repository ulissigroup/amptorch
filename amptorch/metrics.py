import numpy as np
import torch
from skorch.callbacks import EpochScoring, Checkpoint
from amptorch.utils import target_extractor
from torch.nn import L1Loss, MSELoss

def mae_energy_score(net, X, y):
    mae_loss = L1Loss()
    energy_pred, _ = net.forward(X)
    if isinstance(X, torch.utils.data.Subset):
        X = X.dataset
    energy_pred = X.target_scaler.denorm(energy_pred, pred="energy")
    energy_target = X.target_scaler.denorm(
        torch.FloatTensor(np.concatenate(y[::2])), pred="energy"
    )
    energy_loss = mae_loss(energy_pred, energy_target)

    return energy_loss.item()

def mae_forces_score(net, X, y):
    mae_loss = L1Loss()
    _, force_pred = net.forward(X)
    if isinstance(X, torch.utils.data.Subset):
        X = X.dataset
    force_pred = X.target_scaler.denorm(force_pred, pred="forces")
    force_target = X.target_scaler.denorm(
        torch.FloatTensor(np.concatenate(y[1::2])), pred="forces"
    )
    force_loss = mae_loss(force_pred, force_target)

    return force_loss.item()

def mse_energy_score(net, X, y):
    mse_loss = MSELoss()
    energy_pred, _ = net.forward(X)
    if isinstance(X, torch.utils.data.Subset):
        X = X.dataset
    energy_pred = X.target_scaler.denorm(energy_pred, pred="energy")
    energy_target = X.target_scaler.denorm(
        torch.FloatTensor(np.concatenate(y[::2])), pred="energy"
    )
    energy_loss = mse_loss(energy_pred, energy_target)

    return energy_loss.item()

def mse_forces_score(net, X, y):
    mse_loss = MSELoss()
    _, force_pred = net.forward(X)
    if isinstance(X, torch.utils.data.Subset):
        X = X.dataset
    force_pred = X.target_scaler.denorm(force_pred, pred="forces")
    force_target = X.target_scaler.denorm(
        torch.FloatTensor(np.concatenate(y[1::2])), pred="forces"
    )
    force_loss = mse_loss(force_pred, force_target)

    return force_loss.item()

def evaluator(
        val_split,
        metric,
        identifier,
        forcetraining,
    ):

    callbacks = []
    isval = val_split != 0
    if isval:
        cp_on = "val"
    else:
        cp_on = "train"

    if metric == "mae":
        energy_score = mae_energy_score
        forces_score = mae_forces_score
    elif metric == "mse":
        energy_score = mse_energy_score
        forces_score = mse_forces_score
    else:
        raise NotImplementedError(f"{metric} metric not available!")

    callbacks.append(
        EpochScoring(
            energy_score,
            on_train=True,
            use_caching=True,
            name="train_energy_{}".format(metric),
            target_extractor=target_extractor,
        )
    )
    if isval:
        callbacks.append(
            EpochScoring(
                energy_score,
                on_train=False,
                use_caching=True,
                name="val_energy_{}".format(metric),
                target_extractor=target_extractor,
            )
        )
    callbacks.append(
        Checkpoint(
            monitor="{}_energy_{}_best".format(cp_on, metric),
            fn_prefix="checkpoints/{}/".format(identifier),
        )
    )
    if forcetraining:
        callbacks.append(
            EpochScoring(
                forces_score,
                on_train=True,
                use_caching=True,
                name="train_forces_{}".format(metric),
                target_extractor=target_extractor,
            )
        )
        if isval:
            callbacks.append(
                EpochScoring(
                    forces_score,
                    on_train=True,
                    use_caching=True,
                    name="val_forces_{}".format(metric),
                    target_extractor=target_extractor,
                )
            )
        callbacks.append(
            Checkpoint(
                monitor="{}_forces_{}_best".format(cp_on, metric),
                fn_prefix="checkpoints/{}/".format(identifier),
            )
        )
    return callbacks
