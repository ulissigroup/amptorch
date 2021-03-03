import numpy as np
import torch
from skorch.callbacks import Checkpoint, EpochScoring
from torch.nn import L1Loss, MSELoss

from amptorch.utils import target_extractor


def mae_energy_score(net, X, y):
    mae_loss = L1Loss()
    energy_pred, _ = net.forward(X)
    if isinstance(X, torch.utils.data.Subset):
        X = X.dataset
    energy_pred = X.target_scaler.denorm(energy_pred, pred="energy")
    energy_target = X.target_scaler.denorm(
        torch.FloatTensor(np.concatenate([i[0] for i in y])), pred="energy"
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
        torch.FloatTensor(np.concatenate([i[1] for i in y])), pred="forces"
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
        torch.FloatTensor(np.concatenate([i[0] for i in y])), pred="energy"
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
        torch.FloatTensor(np.concatenate([i[1] for i in y])), pred="forces"
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
        MemEffEpochScoring(
            energy_score,
            on_train=True,
            use_caching=True,
            name="train_energy_{}".format(metric),
            target_extractor=target_extractor,
        )
    )
    if isval:
        callbacks.append(
            MemEffEpochScoring(
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
            MemEffEpochScoring(
                forces_score,
                on_train=True,
                use_caching=True,
                name="train_forces_{}".format(metric),
                target_extractor=target_extractor,
            )
        )
        if isval:
            callbacks.append(
                MemEffEpochScoring(
                    forces_score,
                    on_train=False,
                    use_caching=True,
                    name="val_forces_{}".format(metric),
                    target_extractor=target_extractor,
                )
            )

    # callbacks.append(
    #     Checkpoint(
    #         monitor="{}_loss_best".format(cp_on),
    #         fn_prefix="checkpoints/{}/".format(identifier),
    #     )
    # )
    return callbacks


def to_cpu(X):
    if isinstance(X, (tuple, list)):
        return type(X)(to_cpu(x) for x in X)
    return X.detach().to("cpu")


class MemEffEpochScoring(EpochScoring):
    def __init__(
        self,
        scoring,
        lower_is_better=True,
        on_train=False,
        name=None,
        target_extractor=None,
        use_caching=True,
    ):
        super().__init__(
            scoring, lower_is_better, on_train, name, target_extractor, use_caching
        )

    def on_batch_end(self, net, y, y_pred, training, **kwargs):
        if not self.use_caching or training != self.on_train:
            return

        if y is not None:
            self.y_trues_.append(to_cpu(y))
        self.y_preds_.append(to_cpu(y_pred))
