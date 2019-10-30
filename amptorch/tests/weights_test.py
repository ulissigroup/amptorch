"""An example of how to utilize the package to train on energies and forces"""

import sys
from ase import Atoms
import ase
import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
from amptorch.NN_model import CustomLoss, WeightedCustomLoss
from amptorch import AMP
from amptorch.core import AMPTorch
from amptorch.utils import get_hash
from pickle import load
from ase.calculators.emt import EMT
from ase.build import molecule
# from amp.descriptor.gaussian import Gaussian
from amptorch.gaussian import Gaussian


def test_weighting():
    # define training images
    N2 = molecule('N2')
    images = []
    for i in np.linspace(-.1,1,10):
        image = N2.copy()
        image.set_cell([10,10,10])
        image.center()
        image[1].position[2] += i
        image.set_calculator(EMT())
        image.get_forces()
        images.append(image)

    # define symmetry functions to be used
    GSF = {}
    GSF["G2_etas"] = [0.005, 1]
    GSF["G2_rs_s"] = [0, 0]
    GSF["G4_etas"] = [0.005]
    GSF["G4_zetas"] = [1.0, 4.0]
    GSF["G4_gammas"] = [1.0, -1]
    GSF["cutoff"] = 4

    weights_dict = {}
    for image in images:
        weights_dict[get_hash(image, GSF)] = [1., 1.]

    # define the number of threads to parallelize training across
    torch.set_num_threads(6)

    # declare the calculator and corresponding model to be used
    calc = AMP(
        model=AMPTorch(
            images,
            descriptor=Gaussian,
            Gs=GSF,
            cores=6,
            device='cpu',
            force_coefficient=2,
            criterion=WeightedCustomLoss,
            lj_data=None,
            label='example',
            save_logs=True,
            weights_dict=weights_dict,
        )
    )
    # define model settings
    calc.model.structure = [3, 5]
    calc.model.val_frac = 0.2
    calc.model.convergence = {'energy': 0.02, 'force': 0.02}
    calc.model.epochs = 10
    calc.lr = 1
    calc.criterion = CustomLoss
    calc.optimizer = optim.LBFGS

    # train the model
    calc.train(overwrite=True)
