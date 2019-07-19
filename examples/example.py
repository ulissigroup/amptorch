"""An example of how to utilize the package to train on energies and forces"""

from ase import Atoms
import ase
import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
from amp_pytorch.NN_model import CustomLoss
from amp_pytorch import AMP
from amp_pytorch.core import AMPModel
from amp.descriptor.gaussian import Gaussian
from ase.visualize import view

# define training images
# IMAGES = "../datasets/CO_Pt.traj"
IMAGES = "../datasets/water.extxyz"
images = ase.io.read(IMAGES, ":")
IMAGES = []
for i in range(100):
    IMAGES.append(images[i])

# define the number of threads to parallelize training across
torch.set_num_threads(1)
# define calculator, model, and descriptor
# turn force training on by defining a force coefficient>0
# define the number of cores to parallelize across for fingerprint calculations
calc = AMP(model=AMPModel(IMAGES, descriptor=Gaussian(), cores=1,
    force_coefficient=0.3))

# define the convergence criteria
calc.model.val_frac = 0.3
# calc.model.convergence = {"energy": 0.02, "force": 0.02}

# train the model
calc.train(overwrite=True)
# plotting
# calc.train() needs to be run whenever plotting occurs.
# calc.model.parity_plot("energy")
# calc.model.plot_residuals("energy")

# predictions
# energy_predictions = np.concatenate([calc.get_potential_energy(image) for image in images[:10]])
# forces_predictions = np.concatenate([calc.get_forces(image) for image in images[:10]])

