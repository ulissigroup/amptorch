"""An example of how to utilize the package to train on energies and forces"""

import sys, copy, random, time
from ase import Atoms
import ase
import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
from amp_pytorch.NN_model import CustomLoss, LogCoshLoss
from amp_pytorch import AMP
from amp_pytorch.core import AMPModel
from amp.descriptor.gaussian import Gaussian, make_symmetry_functions
from amp_pytorch.lj_model import lj_optim
import ase.io
from ase import units
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md import VelocityVerlet, Langevin
from ase.calculators.emt import EMT
from ase.visualize import view
import os


def ml_lj(IMAGES, filename, count, temp, GSF, dir="MD_results/", const_t=False,
        lj=False, fine_tune=None, activation_fn='l2amp'):
    if not os.path.exists(dir):
        os.mkdir(dir)
    # lj optimization
    lj_data = None
    cutoff = GSF['cutoff']
    if lj:
        p0 = [2.08022879, 1.89536258e-11, -5.47894512e-3, -2.10675310,
                1.94943321e-3, 7.18881277e-3, 2.29444069, 1.00651095e-1,
                3.624737e-2]
        params_dict = {"C": [], "O": [], "Cu": []}
        lj_model = lj_optim(IMAGES, p0, params_dict, cutoff)
        # fitted_params = lj_model.fit(method="L-BFGS-B")
        fitted_params = lj_model.fit()
        # fitted_params = p0
        lj_energies, lj_forces, num_atoms = lj_model.lj_pred(
            IMAGES, fitted_params, params_dict
        )
        lj_data = [
            lj_energies,
            lj_forces,
            num_atoms,
            fitted_params,
            params_dict,
            lj_model,
        ]
    # define the number of threads to parallelize training across
    torch.set_num_threads(1)
    calc = AMP(
        model=AMPModel(
            IMAGES,
            descriptor=Gaussian,
            Gs=GSF,
            cores=1,
            force_coefficient=0.3,
            lj_data=lj_data,
            label=filename,
            save_logs=True,
        ),
    )

    calc.model.convergence = {"energy": 0.002, "force": 0.02}
    calc.model.lr = 1e-2
    # calc.model.fine_tune = fine_tune
    # calc.model.optimizer = optim.SGD
    if activation_fn == 'l2amp':
        calc.model.criterion = CustomLoss
    elif activation_fn == 'logcosh':
        calc.model.criterion = LogCoshLoss
    # calc.model.val_frac = 0.2
    calc.model.structure = [20, 20, 20]

    # train the model
    calc.train(overwrite=True)
    md_run(IMAGES, count, calc, filename, dir, temp, const_t)

def md_run(images, count, calc, filename, dir, temp, cons_t=False):
    """Generates test or training data with a simple MD simulation."""
    traj = ase.io.Trajectory("".join([dir, filename, ".traj"]), "w")
    slab = copy.deepcopy(images[0])
    slab.set_calculator(calc)
    slab.get_forces()
    traj.write(slab)
    if cons_t is True:
        dyn = Langevin(slab, 1.0 * units.fs, temp * units.kB, 0.002)
    else:
        dyn = VelocityVerlet(slab, dt=1.0 * units.fs)
    for step in range(1):
        dyn.run(5)
        traj.write(slab)


'''Runs multiple simulations of ML and LJ models and saves corresponding
trajectory files'''
def multiple_runs(images, filename, dir, num_images, num_iters, temp,
        activation_fn, GSF):
    data = []
    for idx in range(num_images):
        data.append(images[idx])
    for i in range(num_iters):
        lj_name = filename + "_LJ_%s" % str(i+1)
        ml_name = filename + "_%s" % str(i+1)
        ml_lj(data, ml_name, count=num_images, dir=dir, temp=temp,
                GSF=GSF, const_t=True, lj=False, activation_fn=activation_fn)
        # ml_lj(data, lj_name, count=num_images, dir=dir, temp=temp,
                # GSF=GSF, const_t=True, lj=True, activation_fn=activation_fn)

GSF = {}
GSF['G2_etas'] = np.logspace(np.log10(0.05), np.log10(5.0), num=8)
GSF['G2_rs_s'] = [0] * 4
GSF['G4_etas'] = np.array([0.005])
GSF['G4_zetas'] = np.array([1.0, 4.0])
GSF['G4_gammas'] = np.array([1, -1])
GSF['cutoff'] = 4.5

# define training images
images0 = ase.io.read("../../datasets/COCu/COCu_pbc_300K.traj", ":")
# images_aimd = ase.io.read("../datasets/COCu/COCu_pbc_aimd_300K/1.OUTCAR", ":")
# images_LJ = ase.io.read("MD_results/COCu/pbc_300K/val_cl2/MLMD_COCu_pbc_300K_cl2_LJ_1.traj", ":")
# images_ML = ase.io.read("MD_results/COCu/pbc_300K/val_cl2/MLMD_COCu_pbc_300K_cl2_1.traj",  ":")

multiple_runs(images0, filename="MLMD_COCu_pbc_300K_l2amp",
        dir="./", num_images=10, num_iters=1, temp=300, activation_fn='l2amp', GSF=GSF)
