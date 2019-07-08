"""An example of how to utilize the package to train on energies and forces"""

import sys
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
from amp_pytorch.lj_model import lj_optim

# define training images
IMAGES = "../datasets/water.extxyz"
images = ase.io.read(IMAGES, ":")
IMAGES = []
for i in range(100):
    IMAGES.append(images[i])

# lj optimization
p0 = [2.886, 0.044, -2700, 3.50, 0.060, -2700]
p0 = [1, 1, -2700, 1, 1, -2700]
params_dict = {"H": [], "O": []}
cutoff = 6.5
lj_model = lj_optim(IMAGES, p0, params_dict, cutoff)
fitted_params = lj_model.fit()
lj_fitted_data = lj_model.lj_pred(fitted_params)
# lj_model.parity(lj_fitted_data[0], lj_fitted_data[1])

# define the number of threads to parallelize training across
torch.set_num_threads(1)
# define calculator, model, and descriptor
# turn force training on by defining a force coefficient>0
# define the number of cores to parallelize across for fingerprint calculations
calc = AMP(
    model=AMPModel(
        IMAGES,
        descriptor=Gaussian(cutoff=cutoff),
        cores=1,
        force_coefficient=0,
        lj_data=lj_fitted_data
    )
)

# define the convergence criteria
calc.model.convergence = {"energy": 0.002, "force": 0.002}

# train the model
calc.train(overwrite=True)
# plotting
# calc.train() needs to be run whenever plotting occurs.
calc.model.parity_plot("energy")
calc.model.plot_residuals("energy")

# predictions
# energy_predictions = np.concatenate(
    # [calc.get_potential_energy(image) for image in images[:10]]
# )
# forces_predictions = np.concatenate([calc.get_forces(image) for image in images[:10]])
