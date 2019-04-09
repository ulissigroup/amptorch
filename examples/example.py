"""An example of how to utilize the package to train a model utilizing a LBFGS
optimizer"""

import torch
import torch.nn as nn
import torch.optim as optim
from amp_pytorch import core

# locate training images
IMAGES = "../datasets/water.extxyz"

# specify whether a GPU is to be utilized
DEVICE = "cpu"
# DEVICE = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

# define model
MODEL = core.AMPtorch(IMAGES, DEVICE, val_frac=0)

# define training parameters
CRITERION = nn.MSELoss()
OPTIMIZER = optim.LBFGS
RMSE_CRITERIA = 1e-3
LR = 0.8

# train the model
TRAINED_MODEL = MODEL.train(CRITERION, OPTIMIZER, LR, RMSE_CRITERIA)
# plotting
MODEL.parity_plot(TRAINED_MODEL)
MODEL.plot_residuals(TRAINED_MODEL)
