"""An example of how to utilize the package to train a model utilizing GD based
optimizers - including SGD and Adam
"""

import torch
import torch.nn as nn
import torch.optim as optim
from amp_pytorch import core_GD

# locate training images
IMAGES = "../datasets/water.extxyz"

# specify whether a GPU is to be utilized
DEVICE = "cpu"
# DEVICE = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

# define model
MODEL = core_GD.AMPtorch(IMAGES, DEVICE, val_frac=0, batch=None)

# define training parameters
CRITERION = nn.MSELoss()
OPTIMIZER = optim.Adam
RMSE_CRITERIA = 5e-2
LR = 0.01

# train the model
TRAINED_MODEL = MODEL.train(
    criterion=CRITERION, optimizer_ft=OPTIMIZER, lr=LR,
    rmse_criteria=RMSE_CRITERIA
)
# plotting
MODEL.parity_plot(TRAINED_MODEL)
MODEL.plot_residuals(TRAINED_MODEL)
