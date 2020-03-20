"""An example of how to utilize the package to train on energies and forces"""

from ase import Atoms
import sys
import time
import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
from amptorch.model import CustomLoss, TanhLoss
from amptorch import AMP
from amptorch.core import AMPTorch
from amptorch.gaussian import SNN_Gaussian
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

    if lj:
        p0 = [
            1.0,
            6.3535,
            0,
            1.0808,
            8.5357,
            0,
            2.1717,
            3.7575,
            0,
            12
        ]
        params_dict = {"C": [], "O": [], "Cu": []}
        lj_model = lj_optim(IMAGES, p0, params_dict, cutoff, filename)
        fitted_params = lj_model.fit()
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
            descriptor=SNN_Gaussian,
            Gs=GSF,
            force_coefficient=0.04,
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
        "epochs": 500,
        "early_stop": False,
    }
    calc.model.loader_params = {"batch_size": 10, "shuffle": False, "num_workers": 0}
    calc.model.val_frac = 0.2
    calc.model.structure = [3, 5]
    calc.model.optimizer = optim.Adam
    calc.model.scheduler = optim.lr_scheduler.CosineAnnealingLR

    # train the model
    calc.train(overwrite=True)
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

def multiple_runs(images, filename, dir, num_iters, temp, loss_fn, GSF, save_logs):
    for i in range(num_iters):
        lj_name = filename + "-LJ-%s" % str(i + 1)
        ml_name = filename + "-%s" % str(i + 1)
        steps = len(images)
        ml_lj(
            images,
            lj_name,
            count=steps,
            dir=dir,
            temp=temp,
            GSF=GSF,
            const_t=True,
            lj=True,
            loss_fn=loss_fn,
            save_logs=save_logs,
        )

        ml_lj(
            images,
            ml_name,
            count=steps,
            dir=dir,
            temp=temp,
            GSF=GSF,
            const_t=True,
            lj=False,
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
images0 = ase.io.read("../../datasets/COCu/COCu_pbc_300K.traj", ":100")

multiple_runs(
    images=images0,
    filename="test",
    dir="./",
    num_iters=2,
    temp=300,
    loss_fn="l2amp",
    GSF=GSF,
    save_logs=True,
)
