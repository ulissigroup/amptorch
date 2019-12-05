"""An example of how to utilize the package to train on energies and forces
using the physics coupled ML-LJ model"""

import torch
import torch.nn as nn
import torch.optim as optim
from ase import Atoms
from ase.calculators.emt import EMT
import numpy as np
from amptorch.model import CustomLoss
from amptorch import AMP
from amptorch.core import AMPTorch
from amptorch.analysis import parity_plot
from amptorch.gaussian import SNN_Gaussian
from amptorch.lj_model import lj_optim

# define training images
distances = np.linspace(2, 5, 100)
label = "example"
images = []
for l in distances:
    image = Atoms(
        "CuCO",
        [
            (-l * np.sin(0.65), l * np.cos(0.65), 0),
            (0, 0, 0),
            (l * np.sin(0.65), l * np.cos(0.65), 0),
        ],
    )
    image.set_cell([10, 10, 10])
    image.wrap(pbc=True)
    image.set_calculator(EMT())
    images.append(image)

# Define symmetry functions to be used
Gs = {}
Gs["G2_etas"] = np.logspace(np.log10(0.05), np.log10(5.0), num=4)
Gs["G2_rs_s"] = [0] * 4
Gs["G4_etas"] = [0.005]
Gs["G4_zetas"] = [1.0]
Gs["G4_gammas"] = [+1.0, -1]
Gs["cutoff"] = 6.5

#############LJ Optimization##############
p0 = [
    1.33905162,
    0.12290683,
    6.41914719,
    0.64021468,
    0.08010004,
    8.26082762,
    2.29284676,
    0.29639983,
    0.08071821,
]
params_dict = {"C": [], "O": [], "Cu": []}
lj_model = lj_optim(images, p0, params_dict, Gs["cutoff"], label)
fitted_params = lj_model.fit()
lj_energies, lj_forces, num_atoms = lj_model.lj_pred(images, fitted_params, params_dict)
lj_data = [lj_energies, lj_forces, num_atoms, fitted_params, params_dict, lj_model]

############ML Model######################

# define the number of threads to parallelize training across
torch.set_num_threads(1)

# declare the calculator and corresponding model to be used
calc = AMP(
    model=AMPTorch(
        images,
        descriptor=Gaussian,
        Gs=Gs,
        cores=1,
        force_coefficient=0.3,
        label=label,
        save_logs=True,
        lj_data=lj_data,
    )
)
# define model settings
calc.model.device = "cpu"
calc.model.structure = [2, 2]
calc.model.val_frac = 0
calc.model.convergence = {
    "energy": 0.02,
    "force": 0.02,
    "epochs": 1e10,
    "early_stop": False,
}
calc.model.loader_params = {"batch_size": None, "shuffle": False, "num_workers": 0}
calc.model.criterion = CustomLoss
calc.model.optimizer = optim.LBFGS
calc.model.lr = 1e-2
calc.model.fine_tune = None

# train the model
calc.train(overwrite=True)
parity_plot(calc, images, data="energy", label=label)
parity_plot(calc, images, data="forces", label=label)
