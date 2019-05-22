"""An example of how to utilize the package to train on energies and forces"""

import torch
from ase import Atoms
from ase.calculators.emt import EMT
import ase
import torch.nn as nn
import torch.optim as optim
from amp_pytorch import core
from amp_pytorch.NN_model import ForceLossFunction

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
ARCHITECTURE = [3, 5]

# define model
MODEL = core.AMPtorch(IMAGES, DEVICE, batch_size=None, structure=ARCHITECTURE,
        val_frac=0)

# define training parameters
CRITERION = nn.MSELoss(reduction='sum')
OPTIMIZER = optim.LBFGS
RMSE_CRITERIA = 2e-2
LR = 1

# train the model
TRAINED_MODEL = MODEL.train(CRITERION, OPTIMIZER, lr=LR, rmse_criteria=RMSE_CRITERIA)
# plotting
MODEL.parity_plot(TRAINED_MODEL)
MODEL.parity_plot_forces(TRAINED_MODEL)
# MODEL.plot_residuals(TRAINED_MODEL)
