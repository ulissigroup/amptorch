"""An example of how to utilize the package to train on energies and forces"""

import sys
from ase import Atoms
import ase
import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
from amptorch.NN_model import CustomLoss, LogCoshLoss
from amptorch import AMP
from amptorch.core import AMPTorch
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
for i in range(100):
    IMAGES.append(images[i])

# define symmetry functions to be used
GSF = {}
GSF["G2_etas"] = np.logspace(np.log10(0.05), np.log10(5.0), num=4)
GSF["G2_rs_s"] = [0] * 4
GSF["G4_etas"] = [0.005]
GSF["G4_zetas"] = [1.0, 4.0]
GSF["G4_gammas"] = [1.0, -1]
GSF["cutoff"] = 6.5

# define the number of threads to parallelize training across
torch.set_num_threads(1)
calc = AMP(
    model=AMPTorch(
        IMAGES,
        descriptor=Gaussian,
        Gs=GSF,
        cores=1,
        force_coefficient=0.3,
        lj_data=None,
    )
)
# define model settings
calc.model.val_frac = 0.2
calc.model.convergence = {"energy": 0.005, "force": 0.02}
calc.model.structure = [5, 5]
calc.lr = 1
calc.criterion = CustomLoss

# train the model
calc.train(overwrite=True)
