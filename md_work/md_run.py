"""An example of how to utilize the package to train on energies and forces"""

import sys
import time
import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
from amptorch.NN_model import CustomLoss, TanhLoss
from amptorch import AMP
from amptorch.core import AMPTorch
from amptorch.gaussian import Gaussian
from amptorch.lj_model import lj_optim
from amptorch.utils import Logger
from amptorch.analysis import parity_plot
import ase.io
from ase import units
from ase.md import VelocityVerlet, Langevin
from ase.calculators.emt import EMT
from ase.visualize import view
import os
import multiprocessing


def ml_lj(
    IMAGES,
    filename,
    count,
    temp,
    GSF,
    dir="MD_results/",
    const_t=False,
    lj=False,
    fine_tune=None,
    loss_fn="l2amp",
    save_logs=True,
):
    if not os.path.exists(dir):
        os.mkdir(dir)
    # lj optimization
    lj_data = None
    cutoff = GSF["cutoff"]
    eV_kcalmol = 0.043372093
    if lj:
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
        lj_model = lj_optim(IMAGES, p0, params_dict, cutoff, filename)
        # fitted_params = lj_model.fit()
        fitted_params = p0
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
        model=AMPTorch(
            IMAGES,
            descriptor=Gaussian,
            Gs=GSF,
            force_coefficient=0.3,
            lj_data=lj_data,
            label=filename,
            save_logs=save_logs,
        )
    )

    calc.model.lr = 1e-2
    if loss_fn == "l2amp":
        calc.model.criterion = CustomLoss
    elif loss_fn == "tanh":
        calc.model.criterion = TanhLoss
    calc.model.convergence = {
        "energy": 0.02,
        "force": 0.02,
        "epochs": 1000,
        "early_stop": True,
    }
    calc.model.loader_params = {"batch_size": 20, "shuffle": False, "num_workers": 0}
    calc.model.val_frac = 0.2
    calc.model.structure = [30, 30]
    calc.model.optimizer = optim.Adam
    calc.model.scheduler = optim.lr_scheduler.CosineAnnealingLR

    # train the model
    calc.train(overwrite=True)
    parity_plot(calc, IMAGES, filename, data="forces")
    parity_plot(calc, IMAGES, filename, data="energy")
    md_run(IMAGES, count, calc, filename, dir, temp, const_t)


def md_run(images, count, calc, filename, dir, temp, cons_t=False):
    """Generates test or training data with a simple MD simulation."""
    log = Logger("results/logs/" + filename + ".txt")
    traj = ase.io.Trajectory("".join([dir, filename, ".traj"]), "w")
    slab = images[0].copy()
    slab.set_calculator(calc)
    slab.get_forces()
    traj.write(slab)
    if cons_t is True:
        dyn = Langevin(slab, 1.0 * units.fs, temp * units.kB, 0.002)
    else:
        dyn = VelocityVerlet(slab, dt=1.0 * units.fs)
    time_start = time.time()
    for step in range(count):
        dyn.run(20)
        traj.write(slab)
    time_elapsed = time.time() - time_start
    log("MD Simulation Dynamics: %s" % dyn)
    log("MD Simulation Time: %s \n" % time_elapsed)


"""Runs multiple simulations of ML and LJ models and saves corresponding
trajectory files"""


def multiple_runs(
    images, filename, dir, num_images, num_iters, temp, loss_fn, GSF, save_logs
):
    data = []
    for idx in range(num_images):
        data.append(images[idx])
    for i in range(num_iters):
        lj_name = filename + "-LJ-%s" % str(i + 1)
        ml_name = filename + "-%s" % str(i + 1)
        # ml_lj(
            # data,
            # ml_name,
            # count=num_images,
            # dir=dir,
            # temp=temp,
            # GSF=GSF,
            # const_t=True,
            # lj=False,
            # loss_fn=loss_fn,
            # save_logs=save_logs,
        # )
        ml_lj(
            data,
            lj_name,
            count=num_images,
            dir=dir,
            temp=temp,
            GSF=GSF,
            const_t=True,
            lj=True,
            loss_fn=loss_fn,
            save_logs=save_logs,
        )


GSF = {}
GSF["G2_etas"] = np.logspace(np.log10(0.05), np.log10(5.0), num=4)
GSF["G2_rs_s"] = [0] * 4
GSF["G4_etas"] = [0.005]
GSF["G4_zetas"] = [1.0, 4.0]
GSF["G4_gammas"] = np.array([+1.0, -1.0])
GSF["cutoff"] = 5.876798323827276

# define training images
images0 = ase.io.read("../../datasets/COCu/COCu_pbc_300K.traj", ":")

jobs = []
p1 = multiprocessing.Process(target=multiple_runs, args=(
    images0,
    "MLMD_COCu_pbc_300K_l2amp_30x30_reg",
    "MD_results/COCu/pbc_300K/l2amp/",
    100,
    2,
    300,
    "l2amp",
    GSF,
    True,
))
jobs.append(p1)
p1.start()

p2 = multiprocessing.Process(target=multiple_runs, args=(
    images0,
    "MLMD_COCu_pbc_300K_tanh_30x30_reg",
    "MD_results/COCu/pbc_300K/tanh/",
    100,
    2,
    300,
    "tanh",
    GSF,
    True,
))
jobs.append(p2)
p2.start()
