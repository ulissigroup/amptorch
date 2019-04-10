import sys
# from skorch import NeuralNetRegressor
from torch.utils.data import DataLoader
from amp.descriptor.gaussian import Gaussian
from NN_model import FullNN
from data_preprocess import AtomsDataset, factorize_data, collate_amp

k=AtomsDataset('../datasets/water.extxyz', descriptor=Gaussian())
