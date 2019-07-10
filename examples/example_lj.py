"""An example of how to utilize the package to train on energies and forces"""

import sys
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

# define training images
# IMAGES = "../datasets/water.extxyz"
IMAGES = "../datasets/H2_Pt.traj"
images = ase.io.read(IMAGES, ":")


def gen_data(images, count):
    IMAGES = []
    for i in range(count):
        IMAGES.append(images[i])
    return IMAGES


data_size = list(range(10, 1000, 50))
cutoff = 5
torch.set_num_threads(1)


def trainer(count, optimizer=False):
    data = gen_data(images, count)
    if optimizer is True:
        # p0 = [2.886, 0.044, -2700, 3.50 , 0.06, -2700]
        p0 = [0.044, 2.886, 0, 0.060 , 3.50, 0]
        params_dict = {"Pt": [], "H": []}
        lj_model = lj_optim(data, p0, params_dict, cutoff)
        fitted_params = lj_model.fit()
        lj_fitted_data = lj_model.lj_pred(fitted_params)

        calc = AMP(
            model=AMPModel(
                data,
                descriptor=Gaussian(cutoff=cutoff),
                cores=1,
                force_coefficient=0.3,
                lj_data=lj_fitted_data
            )
        )
    else:
        calc = AMP(
            model=AMPModel(
                data, descriptor=Gaussian(cutoff=cutoff), cores=1, force_coefficient=0.3
            )
        )
    calc.model.convergence = {"energy": 0.02, "force": 0.02}
    calc.train(overwrite=True)


for size in data_size:
    trainer(size, False)
    trainer(size, True)
