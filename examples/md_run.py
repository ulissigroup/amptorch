"""An example of how to utilize the package to train on energies and forces"""

import sys
import time
from ase import Atoms
import ase
import ase.db
import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
from amp_pytorch.NN_model import CustomLoss
from amp_pytorch import AMP
from amp_pytorch.core import AMPModel
from amp.descriptor.gaussian import Gaussian
from amp_pytorch.lj_model import lj_optim
from ase.visualize import view
import ase.io
from ase import units
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md import VelocityVerlet, Langevin
from ase.calculators.emt import EMT
import matplotlib.pyplot as plt


# define training images
# IMAGES = "../datasets/COPt/COPt_conT_2x2.traj"
# IMAGES = "../datasets/CuPt/CuPt_conT.traj"
# IMAGES = "../datasets/COCu/COCu.traj"
IMAGES = "../datasets/COPt/COPt_conT.traj"
images = ase.io.read(IMAGES, ":")
IMAGES = []
for i in range(100):
    IMAGES.append(images[i])

# lj optimization
eV_kcalmol = 0.043372093
# p0 = [3.495, 0.005 * eV_kcalmol, 1, 2.754, 0.080 * eV_kcalmol, 1]
# p0 = [4.61454387, 2.69573286e-4, 9.41039131e-1, 2.93605221, 0.00459339, -0.0189848]
# p0 = [3.851, 0.105 * eV_kcalmol, 0, 3.500, 0.060 * eV_kcalmol, 0, 3.495, 0.005 *
        # eV_kcalmol, 0]
p0 = [3.851, 0.105 * eV_kcalmol, 0, 3.500, 0.060 * eV_kcalmol, 0, 2.754, 0.080 *
        eV_kcalmol, 0]
params_dict = {"C": [], "O": [], "Pt": []}
# params_dict = {"Cu": [], "Pt": []}
cutoff = 5.876798323827276
lj_model = lj_optim(IMAGES, p0, params_dict, cutoff)
fitted_params = lj_model.fit()
# fitted_offsets = lj_model.fit()
# p0[2::3] = fitted_offsets
# fitted_params = p0
lj_energies, lj_forces, num_atoms = lj_model.lj_pred(IMAGES, fitted_params, params_dict)
lj_data = [lj_energies, lj_forces, num_atoms, fitted_params, params_dict, lj_model]
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
        force_coefficient=0.03,
        lj_data=lj_data,
    ),
    label="amptorch_COPt.pt",
)

# device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
# calc.model.convergence = {"energy": 0.02, "force": 0.02}
# calc.model.optimizer = optim.Adam
# calc.model.batch_size = 20
calc.model.lr = 0.1
# calc.model.device = device

# train the model
calc.train(overwrite=True)

def generate_data(count, filename, cons_t=False):
    """Generates test or training data with a simple MD simulation."""
    traj = ase.io.Trajectory(filename, "w")
    slab = images[0]
    slab.set_calculator(calc)
    slab.get_potential_energy()
    traj.write(slab)
    MaxwellBoltzmannDistribution(slab, 300.0 * units.kB)
    if cons_t is True:
        dyn = Langevin(slab, 5 * units.fs, 300.0 * units.kB, 0.002)
    else:
        dyn = VelocityVerlet(slab, dt=1.0 * units.fs)
    for step in range(count):
        dyn.run(50)
        traj.write(slab)


generate_data(100, "MD_results/MLMD_LJ_COPt_conT.traj", True)
