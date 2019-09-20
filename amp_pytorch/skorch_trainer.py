import sys
import torch
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
    return(to_numpy(y[0]), to_numpy(y[1]))

def custom_score(net, X, y):
    y_pred = net.predict(X)
    y_target = y[0]
    num_atoms = y[1]
    sys.exit()

    return 3

forcetraining = False
data = AtomsDataset("../datasets/water/water.extxyz", descriptor=Gaussian(),
        cores=1, forcetraining=forcetraining)
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
    # callbacks=[EpochScoring(custom_score, on_train=True, use_caching=True,
        # target_extractor=target_extractor)],
)

net.fit(data, None)
