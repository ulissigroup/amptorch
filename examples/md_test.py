"""An example of how to utilize the package to train on energies and forces"""

import sys
import time
from ase import Atoms
import ase
import ase.db
import torch
import torch.nn as nn
import torch.optim as optim
import scipy.stats as stats
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
import matplotlib.pyplot as plt


# define training images
# IMAGES = "../datasets/COCu/COCu.traj"
# IMAGES = "../datasets/COCu/COCu_pbc.traj"
images = ase.io.read(IMAGES, ":")
IMAGES = []
for i in range(100):
    IMAGES.append(images[i])

# lj optimization
eV_kcalmol = 0.043372093
p0 = [
    3.851,
    0.105 * eV_kcalmol,
    1e-5,
    3.5,
    0.060 * eV_kcalmol,
    1e-5,
    3.495,
    0.005 * eV_kcalmol,
    1e-5,
]
# p0 = [2.362, 0.056 * eV_kcalmol, 0]
params_dict = {"C": [], "O": [], "Cu": []}
# params_dict = {"Cu": [], "Pt": []}
# params_dict = {"C": [], "Cu": []}
# params_dict = {"He": []}
cutoff = 5.876798323827276
# cutoff = 6.5
# lj_model = lj_optim(IMAGES, p0, params_dict, cutoff)
# fitted_params = lj_model.fit()
# fitted_offsets = lj_model.fit()
# p0[2::3] = fitted_offsets
# fitted_params = p0
# lj_energies, lj_forces, num_atoms = lj_model.lj_pred(IMAGES, fitted_params, params_dict)
# lj_data = [lj_energies, lj_forces, num_atoms, fitted_params, params_dict, lj_model]
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
        force_coefficient=0.3,
        # lj_data=lj_data,
    ),
    label="test.pt",
)

# device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
# calc.model.convergence = {"energy": 0.02, "force": 0.02}
# calc.model.optimizer = optim.Adam
# calc.model.batch_size = 20
calc.model.lr = 0.1
# calc.model.device = device

# train the model
calc.train(overwrite=True)
calc.model.parity_plot()
calc.model.parity_plot("forces")
sys.exit()


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


generate_data(20, "test_MLLJ.traj")
sys.exit()


# generate_data(100, "MD_results/MLMD_LJ_COPt.traj", True)

images0 = ase.io.read("../datasets/CuPt/CuPt_conT.traj", ":")
images1 = ase.io.read("MD_results/MLMD_CuPt_consT.traj", ":")
images2 = ase.io.read("MD_results/MLMD_LJ_CuPt.traj", ":")
forces0 = np.array([np.amax(np.abs(image.get_forces())) for image in images0]).reshape(
    -1, 1
)
forces1 = np.array([np.amax(np.abs(image.get_forces())) for image in images1]).reshape(
    -1, 1
)
forces2 = np.array([np.amax(np.abs(image.get_forces())) for image in images2]).reshape(
    -1, 1
)

sd0 = np.std(forces0)
avg0 = np.mean(forces0)
sd1 = np.std(forces1)
avg1 = np.mean(forces1)
sd2 = np.std(forces2)
avg2 = np.mean(forces2)
forces0 = np.sort(forces0, 0)
forces1 = np.sort(forces1, 0)
forces2 = np.sort(forces2, 0)

# forces0 = forces0[forces0 < avg0 + 2 * sd0]
# forces1 = forces1[forces1 < avg1 + 2 * sd1]
# forces2 = forces2[forces2 < avg2 + 2 * sd2]


def force_dist(forces0, forces1, forces2, data):
    plt.plot(forces0, stats.norm.pdf(forces0, avg0, sd0), label="EMT")
    plt.plot(forces1, stats.norm.pdf(forces1, avg1, sd1), label="ML")
    plt.plot(forces2, stats.norm.pdf(forces2, avg2, sd2), label="ML-LJ")
    # plt.hist(forces0, bins=100, alpha=0.5, label="EMT")
    # plt.hist(forces1, bins=100, alpha=0.5, label="ML")
    # plt.hist(forces2, bins=100, alpha=0.5, label="ML-LJ")
    plt.legend(loc="upper right")
    plt.savefig("MD_results/" + data + "_force_dist.pdf")
    plt.show()


force_dist(forces0, forces1, forces2, "ALL_CuPt_gauss")
