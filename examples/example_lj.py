"""An example of how to utilize the package to train on energies and forces"""

import os
from ase import Atoms
import ase
import ase.db
import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import matplotlib.pyplot as plt
from .NN_model import CustomLoss
from . import AMP
from .core import AMPModel
from ase.calculators.emt import EMT
from amp.descriptor.gaussian import Gaussian
from .lj_model import lj_optim


def gen_data(images, count):
    IMAGES = []
    for i in range(count):
        IMAGES.append(images[i])
    return IMAGES


def trainer(images, count, optimizer=False):
    data = gen_data(images, count)
    eV_kcalmol = 0.043372093

    GSF = {}
    GSF["G2_etas"] = np.logspace(np.log10(0.05), np.log10(5.0), num=4)
    GSF["G2_rs_s"] = [0] * 4
    GSF["G4_etas"] = [0.005]
    GSF["G4_zetas"] = [1.0, 4.0]
    GSF["G4_gammas"] = [1.0, -1]
    GSF["cutoff"] = 6.5

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
        lj_model = lj_optim(data, p0, params_dict, GSF["cutoff"])
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
                descriptor=Gaussian,
                Gs=GSF,
                cores=1,
                force_coefficient=0.03,
                lj_data=lj_data,
            )
        )
    else:
        calc = AMP(
            model=AMPModel(
                data,
                descriptor=Gaussian
                Gs=GSF,
                cores=1,
                force_coefficient=0.03,
            )
        )
    # calc.model.convergence = {"energy": 0.02, "force": 0.10}
    calc.model.lr = 0.1
    calc.train(overwrite=True)
