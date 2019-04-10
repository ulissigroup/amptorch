import sys
import torch.nn as nn
from skorch import NeuralNetRegressor
from torch.utils.data import DataLoader
from amp.descriptor.gaussian import Gaussian
from NN_model import FullNN
from data_preprocess import AtomsDataset, factorize_data, collate_amp

training_data = AtomsDataset("../datasets/water.extxyz", descriptor=Gaussian())
unique_atoms, _, _, _ = factorize_data(training_data)
size = len(training_data)

atoms_dataloader = DataLoader(
    training_data, size, collate_fn=collate_amp, shuffle=False
)

net = NeuralNetRegressor(
    module=FullNN,
    criterion=nn.MSELoss,
    iterator_train_collate_fn=collate_amp,
)

net.fit(training_data, None)
