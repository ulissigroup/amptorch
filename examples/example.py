from amp_pytorch import core
import torch.nn as nn
import torch.optim as optim

images = "../datasets/water.extxyz"
device = "cpu"

model = core.AMPtorch(images, device)

criterion = nn.MSELoss()
optimizer = optim.LBFGS
rmse_criteria = 2e-3

model.train(criterion, optimizer, rmse_criteria)
