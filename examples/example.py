"""An example of how to utilize the package to train on energies and forces"""

import sys
from ase import Atoms
import torch
import numpy as np
import torch.nn as nn
import torch.optim as optim
import numpy as np
import ase
from amptorch.NN_model import CustomLoss
from amptorch import AMP
from amptorch.core import AMPTorch
#from amp.descriptor.gaussian import Gaussian
from amptorch.gaussian import Gaussian

# define training images
IMAGES = "../datasets/water/water.extxyz"
#IMAGES = "/home/bcomer3/data/simple_NN/nn_trains/psi4/fixed_md.traj"
images = ase.io.read(IMAGES, ":")
IMAGES = []
for i in range(360):
    IMAGES.append(images[i])

# define symmetry functions to be used
GSF = {}
GSF["G2_etas"] = np.logspace(np.log10(0.05), np.log10(80.0), num=4)
GSF["G2_rs_s"] = [0] * 4
GSF["G4_etas"] = [0.005]
GSF["G4_zetas"] = [1.0, 4.0]
GSF["G4_gammas"] = [1.0, -1]
GSF["cutoff"] = 6.5

# define the number of threads to parallelize training across
torch.set_num_threads(1)

label = 'example'
cores=1
# declare the calculator and corresponding model to be used
calc = AMP(
    model=AMPTorch(
        images,
        descriptor=Gaussian,
        Gs=GSF,
        force_coefficient=0.3,
        label=label,
        db_path='/storage/home/hhive1/bcomer3/data/mass_amps/fps/',
        save_logs=True,
    )
)
# define model settings
calc.model.device = "cpu"
calc.model.structure = [3, 10]
calc.model.val_frac = 0
calc.model.convergence = {
    "energy": 0.02,
    "force": 0.02,
    "epochs": 5,
    "early_stop": False,
}
calc.model.loader_params = {
        "batch_size": None,
        "shuffle": False,
        "num_workers": 0}
calc.model.criterion = CustomLoss
calc.model.optimizer = optim.LBFGS
calc.model.lr = 1e-2
calc.model.fine_tune = None

# train the model
calc.train(overwrite=True)
parity_plot(calc, images, data="energy", label=label)
parity_plot(calc, images, data="forces", label=label)
