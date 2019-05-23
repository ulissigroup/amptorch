"""An example of how to utilize the package to train on energies and forces"""

import sys
import torch
from ase import Atoms
from ase.calculators.emt import EMT
import ase
import torch.nn as nn
import torch.optim as optim
from amp_pytorch.NN_model import CustomLoss
from amp_pytorch import AMPCalc
from amp_pytorch.core import AMPtorch
from amp.descriptor.gaussian import Gaussian

# define training images
IMAGES = "../datasets/water.extxyz"
# IMAGES = "../datasets/reaxff_data/15.traj"
images = ase.io.read(IMAGES, ":")
IMAGES = []
for i in range(300):
    IMAGES.append(images[i])

# specify whether a GPU is to be utilized
DEVICE = "cpu"
# DEVICE = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

# NN architectures across different atoms are identical with the first index
# representing the number of layers, and the second number representing the
# number of nodes in each hidden layer. i.e. [3,5] = 3 layers (2 hidden layers,
# 1 output layer) and 5 nodes in each hidden layer.

# define model
calc = AMPCalc(model=AMPtorch(IMAGES, descriptor=Gaussian()))

# train the model
TRAINED_MODEL = calc.train()
# calc.train(CRITERION, OPTIMIZER, lr=LR, rmse_criteria=RMSE_CRITERIA)
# plotting
calc.model.parity_plot(TRAINED_MODEL)
calc.model.parity_plot_forces(TRAINED_MODEL)
# MODEL.plot_residuals(TRAINED_MODEL)
