"""An example of how to utilize the package to train on energies and forces"""

from ase import Atoms
import ase
import sys
import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
from amp_pytorch.NN_model import CustomLoss
from amp_pytorch import AMP
from amp_pytorch.core import AMPModel
from amp.descriptor.gaussian import Gaussian

# define training images
IMAGES = "../datasets/water.extxyz"
images = ase.io.read(IMAGES, ":")
IMAGES = []
for i in range(300):
    IMAGES.append(images[i])

# define the number of threads to parallelize across
torch.set_num_threads(2)
# define calculator, model, and descriptor
# turn force training on by defining a force coefficient>0
calc = AMP(model=AMPModel(IMAGES, descriptor=Gaussian(), force_coefficient=0))

# define the convergence criteria
calc.model.convergence = {"energy": 0.002, "force": 0.02}

# train the model
calc.train(overwrite=True)
# predictions
energy_predictions = np.concatenate([calc.get_potential_energy(image) for image in images[:10]])
forces_predictions = np.concatenate([calc.get_forces(image) for image in images[:10]])

# plotting
# calc.train() needs to be run whenever plotting occurs.
calc.model.parity_plot("energy")
calc.model.plot_residuals("energy")
