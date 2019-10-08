"""An example of how to utilize the package to train on energies and forces"""

import sys
from ase import Atoms
import ase
import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
from amp_pytorch.NN_model import CustomLoss, LogCoshLoss
from amp_pytorch import AMP
from amp_pytorch.core import AMPModel
from amp.descriptor.gaussian import Gaussian
from ase.visualize import view
from amp_simple_nn.convert import make_amp_descriptors_simple_nn
from ase.build import molecule
from ase.calculators.emt import EMT

# define training images
# IMAGES = "../datasets/water/water.extxyz"
IMAGES = "../datasets/COCu/COCu_pbc_300K.traj"
images = ase.io.read(IMAGES, ":")
IMAGES = []
for i in range(2):
    IMAGES.append(images[i])

G2_etas = [0.005]
G2_rs_s = [0] * 4
G4_etas = [0.005]
G4_zetas = [1., 4.]
G2_gammas = [1., -1]
cutoff = 6

# define the number of threads to parallelize training across
torch.set_num_threads(1)
# define calculator, model, and descriptor
# turn force training on by defining a force coefficient>0
# define the number of cores to parallelize across for fingerprint calculations
calc = AMP(model=AMPModel(IMAGES, descriptor=Gaussian(), cores=1,
    force_coefficient=0.3))

# calc.model.val_frac = 0.3
calc.model.convergence = {"energy": 0.005, "force": 0.02}
# calc.model.fine_tune = "results/trained_models/amptorch.pt"
calc.model.structure = [5, 5]
calc.model.batch_size = 100
calc.lr = 1
# calc.criterion = LogCoshLoss
calc.criterion = CustomLoss

# train the model
calc.train(overwrite=True)
# plotting
# calc.train() needs to be run whenever plotting occurs.
# calc.model.parity_plot("energy")
# calc.model.plot_residuals("energy")

# predictions
# energy_predictions = np.concatenate([calc.get_potential_energy(image) for image in images[:10]])
# forces_predictions = np.concatenate([calc.get_forces(image) for image in images[:10]])

