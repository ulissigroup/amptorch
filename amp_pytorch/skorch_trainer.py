import sys
import torch
from torch.nn import MSELoss
from skorch import NeuralNetRegressor
from skorch.dataset import CVSplit
from skorch.callbacks import Checkpoint, EpochScoring
import skorch.callbacks.base
from amp.descriptor.gaussian import Gaussian
from amp_pytorch.model_skorch import FullNN
from amp_pytorch.model_skorch import CustomLoss
from amp_pytorch.skorch_data import AtomsDataset, factorize_data, collate_amp
from skorch.utils import to_numpy


def target_extractor(y):
    return(to_numpy(y[0]), to_numpy(y[1])) if len(y) == 2 else (to_numpy(y[0]),
            to_numpy(y[1]), to_numpy(y[2]))

def energy_score(net, X, y):
    mse_loss = MSELoss(reduction="sum")
    energy_pred = net.infer()[0]
    energy_target = torch.tensor(y[0])
    num_atoms = torch.tensor(y[1])
    dataset_size = len(energy_pred)
    sd_scaling = scalings[0]
    mean_scaling = scalings[1]
    raw_preds = (energy_pred * sd_scaling) + mean_scaling
    raw_preds_per_atom = torch.div(raw_preds, num_atoms)
    raw_targets = (energy_target * sd_scaling) + mean_scaling
    target_per_atom = torch.div(raw_targets, num_atoms)
    energy_loss = mse_loss(raw_preds_per_atom, target_per_atom)
    energy_loss /= dataset_size
    return energy_loss


forcetraining = False
data = AtomsDataset("../datasets/water/water.extxyz", descriptor=Gaussian(),
        cores=1, forcetraining=forcetraining)
scalings = data.scalings
unique_atoms = data.unique_atoms

device = 'cpu'

net = NeuralNetRegressor(
    module=FullNN(unique_atoms, [20, 20, 20], device,
        forcetraining=forcetraining),
    criterion=CustomLoss,
    criterion__force_coefficient=0,
    optimizer=torch.optim.LBFGS,
    lr=1,
    batch_size=400,
    max_epochs=100,
    iterator_train__collate_fn=collate_amp,
    iterator_valid__collate_fn=collate_amp,
    device=device,
    train_split=None,
    callbacks=[EpochScoring(energy_score, on_train=True, use_caching=True,
        target_extractor=target_extractor)],
)

net.fit(data, None)
