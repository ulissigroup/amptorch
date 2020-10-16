import numpy as np
import skorch
import torch
from skorch.utils import to_numpy
from torch.nn import MSELoss


def target_extractor(y):
    # TODO better remove need for numpy call here before GPU support
    if len(y) == 2:
        return (to_numpy(y[0]), to_numpy(y[1]))
    elif len(y) == 1:
        return (to_numpy(y[0]), None)


def to_tensor(X, device, accept_sparse=False):
    return X


def energy_score(net, X, y):
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


def forces_score(net, X, y):
    mse_loss = MSELoss()

    
    _, force_pred = net.forward(X)
    print("=====================================")
    print(force_pred)
    if isinstance(X, torch.utils.data.Subset):
        X = X.dataset
    force_pred = X.target_scaler.denorm(force_pred, pred="forces")
    force_target = X.target_scaler.denorm(
        torch.FloatTensor(np.concatenate(y[1::2])), pred="forces"
    )
    force_loss = mse_loss(force_pred, force_target)

    return force_loss.item()


class train_end_load_best_loss(skorch.callbacks.base.Callback):
    def __init__(self, filename):
        self.filename = filename

    def on_train_end(self, net, X, y):
        net.load_params("./checkpoints/{}/params.pt".format(self.filename))
