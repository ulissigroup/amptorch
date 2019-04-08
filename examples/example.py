"""An example of how to utilize the package to train a model utilizing a LBFGS
optimizer"""

import torch
import torch.nn as nn
import torch.optim as optim
from amp_pytorch import core

# locate training images
IMAGES = "../datasets/defect-trajectory.extxyz"

# specify whether a GPU is to be utilized
# DEVICE = "cpu"
DEVICE = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

# define model
MODEL = core.AMPtorch(IMAGES, DEVICE, val_frac=0)

# define training parameters
CRITERION = nn.MSELoss()
OPTIMIZER = optim.LBFGS
RMSE_CRITERIA = 2e-3
LR = 1

# train the model
MODEL.train(CRITERION, OPTIMIZER, LR, RMSE_CRITERIA)
