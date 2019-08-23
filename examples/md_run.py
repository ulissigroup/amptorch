"""An example of how to utilize the package to train on energies and forces"""

import sys
import time
from ase import Atoms
import ase
import ase.db
import torch
import random
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

def ml_lj(IMAGES, filename, count, const_t=False, lj=False):
    # lj optimization
    lj_data = None
    cutoff = 5.876798323827276
    if lj:
        eV_kcalmol = 0.043372093
        p0 = [
            0.005,
            0.105 * eV_kcalmol,
            0,
            0.005,
            0.060 * eV_kcalmol,
            0,
            0.005,
            0.005 * eV_kcalmol,
            0,
        ]
        params_dict = {"C": [], "O": [], "Cu": []}
        lj_model = lj_optim(IMAGES, p0, params_dict, cutoff)
        fitted_params = lj_model.fit(method="L-BFGS-B")
        # fitted_params = lj_model.fit()
        lj_energies, lj_forces, num_atoms = lj_model.lj_pred(IMAGES, fitted_params, params_dict)
        lj_data = [lj_energies, lj_forces, num_atoms, fitted_params, params_dict, lj_model]

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
            force_coefficient=0.3,
            lj_data=lj_data,
        ),
        label=filename+".pt",
    )

    # calc.model.device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    calc.model.convergence = {"energy": 0.02, "force": 0.1}
    calc.model.lr = 1
    # calc.model.val_frac = 0.2
    calc.model.structure = [10, 10, 10]

    # train the model
    calc.train(overwrite=True)
    md_run(IMAGES, count, calc, filename, const_t)

def md_run(images, count, calc, filename, cons_t=False):
    """Generates test or training data with a simple MD simulation."""
    traj = ase.io.Trajectory("MD_results/"+filename+".traj", "w")
    slab = images[0]
    slab.set_calculator(calc)
    slab.get_potential_energy()
    traj.write(slab)
    MaxwellBoltzmannDistribution(slab, 300.0 * units.kB)
    if cons_t is True:
        dyn = Langevin(slab, 1 * units.fs, 300.0 * units.kB, 0.002)
    else:
        dyn = VelocityVerlet(slab, dt=1.0 * units.fs)
    for step in range(count):
        dyn.run(50)
        traj.write(slab)

# define training images
IMAGES0 = "../datasets/COCu/COCu.traj"
images0 = ase.io.read(IMAGES0, ":")
data0 = []
for i in range(100):
    data0.append(images0[i])
images1 = ase.io.read("MD_results/MLMD_LJ_COCu.traj", ":")
points = random.sample(range(100), 20)
# points = np.arange(10)
for point in points:
    image = images1[point]
    image.set_calculator(EMT())
    data0.append(image)


ml_lj(data0, "MLMD_LJ_COCu_sample3", count=100, const_t=False, lj=True)
# ml_lj(data0, "MLMD_COCu_new", count=100, const_t=False, lj=False)
