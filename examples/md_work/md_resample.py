"""An example of how to utilize the package to train on energies and forces"""
import sys
import copy
import random
import time
from ase import Atoms
import ase
import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
from amptorch.NN_model import CustomLoss, TanhLoss, MSLELoss
from amptorch import AMP
from amptorch.core import AMPTorch
from amptorch.gaussian import Gaussian, make_symmetry_functions
from amptorch.lj_model import lj_optim
from amptorch.analysis import parity_plot
from amptorch.utils import Logger
import ase.io
from ase import units
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
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
    resample=None,
):
    if not os.path.exists(dir):
        os.mkdir(dir)
    # lj optimization
    lj_data = None
    cutoff = GSF["cutoff"]
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
        lj_model = lj_optim(
            filename, IMAGES, p0, params_dict, cutoff, forcetraining=True
        )
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

    # calc.model.device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    calc.model.convergence = {
        "energy": 0.02,
        "force": 0.02,
        "epochs": 500,
        "early_stop": False,
    }
    calc.model.loader_params = {"batch_size": None, "shuffle": False, "num_workers": 0}
    calc.model.lr = 1
    if loss_fn == "tanh":
        calc.model.criterion = TanhLoss
    elif loss_fn == "l2amp":
        calc.model.criterion = CustomLoss
    calc.model.val_frac = 0.2
    calc.model.resample = resample
    calc.model.structure = [2, 2]

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


"""Runs multiple simulations of resampled LJ models and saves corresponding
trajectory files"""


def multiple_samples(
    images,
    sample_images,
    filename,
    dir,
    num_images,
    sample_points,
    num_iters,
    temp,
    lj,
    GSF,
    loss_fn,
    fine_tune=None,
    save_logs=True,
    resample=None,
):
    data = [images[idx] for idx in range(num_images)]
    for idx in sample_points:
        sample_images[idx].set_calculator(EMT())
        data.append(sample_images[idx])
    for i in range(num_iters):
        name = filename + "_%s_resample_%s" % (num_samples, str(i + 1))
        if lj:
            name = filename + "_LJ_%s_resample_%s" % (num_samples, str(i + 1))
        ml_lj(
            data,
            name,
            count=num_images,
            dir=dir,
            temp=temp,
            GSF=GSF,
            const_t=True,
            lj=lj,
            loss_fn=loss_fn,
            fine_tune=fine_tune,
            save_logs=save_logs,
            resample=sample_points,
        )


# define training images
images0 = ase.io.read("../../datasets/COCu/COCu_pbc_300K.traj", ":")
images_ML_l2 = ase.io.read(
        "MD_results/COCu/pbc_300K/l2amp/paper/MLMD_COCu_pbc_300K_l2amp-2.traj", ":"
)
images_LJ_l2 = ase.io.read(
    "MD_results/COCu/pbc_300K/l2amp/paper/MLMD_COCu_pbc_300K_l2amp-LJ-2.traj", ":"
)
images_ML_tanh = ase.io.read(
        "MD_results/COCu/pbc_300K/tanh/paper/MLMD_COCu_pbc_300K_tanh-1.traj", ":"
)
images_LJ_tanh = ase.io.read(
    "MD_results/COCu/pbc_300K/tanh/paper/MLMD_COCu_pbc_300K_tanh-LJ-1.traj", ":"
)

GSF = {}
GSF["G2_etas"] = np.logspace(np.log10(0.05), np.log10(5.0), num=4)
GSF["G2_rs_s"] = [0] * 4
GSF["G4_etas"] = np.array([0.005])
GSF["G4_zetas"] = np.array([1.0])
GSF["G4_gammas"] = np.array([1, -1])
GSF["cutoff"] = 5.876798323827276

num_images = 100
num_samples = 10
sample_points = random.sample(range(1, num_images), num_samples)
jobs = []
# p1 = multiprocessing.Process(target=multiple_samples, args=(
        # images0,
        # images_ML_l2,
        # "MLMD_COCu_pbc_300K_l2amp",
        # "MD_results/COCu/pbc_300K/l2amp/paper/",
        # 100,
        # sample_points,
        # 2,
        # 300,
        # False,
        # GSF,
        # "l2amp",
        # None,
        # True,
        # None
    # ))
# jobs.append(p1)
# p1.start()
# p2 = multiprocessing.Process(target=multiple_samples, args=(
        # images0,
        # images_LJ_l2,
        # "MLMD_COCu_pbc_300K_l2amp",
        # "MD_results/COCu/pbc_300K/l2amp/paper/",
        # 100,
        # sample_points,
        # 2,
        # 300,
        # True,
        # GSF,
        # "l2amp",
        # None,
        # True,
        # None
    # ))
# jobs.append(p2)
# p2.start()
p3 = multiprocessing.Process(target=multiple_samples, args=(
        images0,
        images_ML_tanh,
        "MLMD_COCu_pbc_300K_tanh2",
        "MD_results/COCu/pbc_300K/tanh/paper/",
        100,
        sample_points,
        2,
        300,
        False,
        GSF,
        "tanh",
        None,
        True,
        None
    ))
jobs.append(p3)
p3.start()
p4 = multiprocessing.Process(target=multiple_samples, args=(
        images0,
        images_LJ_tanh,
        "MLMD_COCu_pbc_300K_tanh2",
        "MD_results/COCu/pbc_300K/tanh/paper/",
        100,
        sample_points,
        2,
        300,
        True,
        GSF,
        "tanh",
        None,
        True,
        None
    ))
jobs.append(p4)
p4.start()
