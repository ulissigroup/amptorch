"""An example of how to utilize the package to train on energies and forces"""

import sys
from ase import Atoms
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

# define model
calc = AMPCalc(model=AMPtorch(IMAGES, descriptor=Gaussian()))
"""model parameters can be modified directly here:

   calc.model.device: Define whether training is to be performed on a cpu or
   GPU. default: 'cpu'

   calc.model.batch_size: Number of images in a training batch, where 'None'
   represents the entire dataset as one batch. default: None

   calc.model.structure: First index represents the number of
   layers, including output layer, and the second index represents the number
   of nodes in each hidden layer. i.e. [3,5] = 3 layers (2 hidden layers, 1
   output layer) and 5 nodes in each hidden layer. default: [3,5]

   calc.model.val_frac: Proportion of dataset to be used as a validation test
   set. default: 0

   calc.model.lossfunction: Define the loss function, and in turn, whether
   force training is on or off.
   default:CustomLoss(force_coefficient=0)

   calc.model.optimizer: Define the training optimizer to be utilized.
   default: optim.LBFGS

   calc.model.scheduler: Specify whether a learning rate decay scheme is to be
   utilized. default:None

   calc.model.lr: Define the model learning rate. default: 1

   calc.model.convergence: Define the training convergence criteria.
   default: {'energy':0.02, 'force': 0.02}

"""
calc.model.lossfunction = CustomLoss(force_coefficient=0.4)
calc.model.convergence = {"energy": 0.02, "force": 0.02}
# train the model
calc.train(overwrite=True)
# predictions = calc.get_potential_energy([images[0], images[1], images[2]])

# plotting
calc.model.parity_plot()
calc.model.parity_plot_forces()
calc.model.plot_residuals('energy')
