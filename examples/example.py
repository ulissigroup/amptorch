"""An example of how to utilize the package to train a model utilizing a LBFGS
optimizer"""

import numpy as np
from ase import Atoms
from ase.calculators.emt import EMT
import torch
import torch.nn as nn
import torch.optim as optim
from amp_pytorch import core

# locate training images
# IMAGES = "../datasets/water.extxyz"
IMAGES = [Atoms(symbols='PdOPd',
                        pbc=np.array([True, False, False], dtype=bool),
                        calculator=EMT(),
                        cell=np.array(
                            [[2.,  0.,  0.],
                             [0.,  2.,  0.],
                             [0.,  0.,  2.]]),
                        positions=np.array(
                            [[0.5,  1., 0.5],
                             [1.,  0.5,  1.],
                             [1.5,  1.5,  1.5]])),
                  Atoms(symbols='PdO',
                        pbc=np.array([True, True, False], dtype=bool),
                        calculator=EMT(),
                        cell=np.array(
                            [[2.,  0.,  0.],
                             [0.,  2.,  0.],
                                [0.,  0.,  2.]]),
                        positions=np.array(
                            [[0.5,  1., 0.5],
                             [1.,  0.5,  1.]])),
                  Atoms(symbols='Cu',
                        pbc=np.array([True, True, False], dtype=bool),
                        calculator=EMT(),
                        cell=np.array(
                            [[1.8,  0.,  0.],
                             [0.,  1.8,  0.],
                                [0.,  0.,  1.8]]),
                        positions=np.array(
                            [[0., 0., 0.]]))]

# specify whether a GPU is to be utilized
DEVICE = "cpu"
# DEVICE = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

# NN architectures across different atoms are identical with the first index
# representing the number of layers, and the second number representing the
# number of nodes in each hidden layer. i.e. [3,3] = 3 layers (2 hidden layers,
# 1 output layer) and 3 nodes in each hidden layer.
ARCHITECTURE = [3, 3]

# define model
MODEL = core.AMPtorch(IMAGES, DEVICE, structure=ARCHITECTURE, val_frac=0)

# define training parameters
CRITERION = nn.MSELoss()
OPTIMIZER = optim.LBFGS
RMSE_CRITERIA = 2e-3
LR = 0.8

# train the model
TRAINED_MODEL = MODEL.train(CRITERION, OPTIMIZER, LR, RMSE_CRITERIA)
# plotting
# MODEL.parity_plot(TRAINED_MODEL)
# MODEL.plot_residuals(TRAINED_MODEL)
