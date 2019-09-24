"""An example of how to utilize the package to train on energies and forces"""

import sys
import os
import time
import json
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
from ase.calculators.emt import EMT
from amp.descriptor.gaussian import Gaussian
from amp_pytorch.lj_model import lj_optim
import matplotlib.pyplot as plt


def gen_data(images, count):
    IMAGES = []
    for i in range(count):
        IMAGES.append(images[i])
    return IMAGES


def trainer(images, count, optimizer=False):
    data = gen_data(images, count)
    eV_kcalmol = 0.043372093
    if optimizer is True:
        p0 = [
            3.851,
            0.105 * eV_kcalmol,
            0,
            3.500,
            0.060 * eV_kcalmol,
            0,
            3.495,
            0.005 * eV_kcalmol,
            0,
        ]
        params_dict = {"C": [], "O": [], "Cu": []}
        lj_model = lj_optim(data, p0, params_dict, cutoff)
        fitted_params = lj_model.fit()
        lj_energies, lj_forces, num_atoms = lj_model.lj_pred(
            data, fitted_params, params_dict
        )
        lj_data = [
            lj_energies,
            lj_forces,
            num_atoms,
            fitted_params,
            params_dict,
            lj_model,
        ]

        calc = AMP(
            model=AMPModel(
                data,
                descriptor=Gaussian(cutoff=cutoff),
                cores=1,
                force_coefficient=0.03,
                lj_data=lj_data,
            )
        )
    else:
        calc = AMP(
            model=AMPModel(
                data,
                descriptor=Gaussian(cutoff=cutoff),
                cores=1,
                force_coefficient=0.03,
            )
        )
    # calc.model.convergence = {"energy": 0.02, "force": 0.10}
    calc.model.lr = 0.1
    calc.train(overwrite=True)
    energy_predictions = np.concatenate(
        [calc.get_potential_energy(image) for image in IMAGES]
    ).reshape(-1, 1)
    force_predictions = np.concatenate([calc.get_forces(image) for image in IMAGES])
    MSE_e = ((energy_predictions - energy_targets) ** 2).mean()
    MSE_f = ((force_predictions - force_targets) ** 2).mean()
    return MSE_e, MSE_f


def create_learning_plot(data, data_lj, datatype):
    fig = plt.figure(figsize=(7.0, 7.0))
    ax = fig.add_subplot(111)
    ax.plot(data_size, data, "b-", lw=0.3, label='ML')
    ax.plot(data_size, data_lj, "r-", lw=0.3, label='ML-LJ')
    ax.set_xlabel("# of Training Data")
    ax.set_ylabel("MSE")
    ax.set_title("MSE vs. # Data")
    fig.savefig("results/" + datatype + "_mse.png")
    plt.show()


if not os.path.exists("MD_results"):
    os.mkdir("MD_results")

images = ase.io.read("../datasets/COCu/COCu.traj", ":")
IMAGES = []
for i in range(100):
    IMAGES.append(images[i])

energy_targets = np.array([image.get_potential_energy() for image in IMAGES]).reshape(
    -1, 1
)
force_targets = np.array([np.amax(np.abs(image.get_forces())) for image in IMAGES])
force_targets = np.concatenate([image.get_forces() for image in IMAGES])

data_size = np.arange(10, 110, 10)
cutoff = 5.876798323827276
torch.set_num_threads(1)
