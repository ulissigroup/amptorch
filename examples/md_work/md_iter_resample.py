import os
import random
import time
from ase import Atoms
import ase
import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
from amptorch.NN_model import CustomLoss, TanhLoss
from amptorch import AMP
from amptorch.core import AMPTorch
from amptorch.gaussian import Gaussian, make_symmetry_functions
from amptorch.lj_model import lj_optim
from amptorch.utils import Logger
from amptorch.analysis import parity_plot
import ase.io
from ase import units
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md import VelocityVerlet, Langevin
from ase.calculators.emt import EMT
from ase.visualize import view
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
        lj_model = lj_optim(filename, IMAGES, p0, params_dict, cutoff)
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

    calc.model.convergence = {
        "energy": 0.02,
        "force": 0.02,
        "epochs": 500,
        "early_stop": False,
    }
    calc.model.loader_params = {"batch_size": None, "shuffle": False, "num_workers": 0}
    calc.model.lr = 1e-2
    if loss_fn == "l2amp":
        calc.model.criterion = CustomLoss
    elif loss_fn == "tanh":
        calc.model.criterion = TanhLoss
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


def sampler(
    images,
    sample_images,
    filename,
    dir,
    num_images,
    num_samples,
    iteration,
    temp,
    lj,
    GSF,
    loss_fn,
    fine_tune=None,
    save_logs=True,
):
    data = [images[idx] for idx in range(num_images)]
    sample_points = random.sample(range(1, num_images), num_samples)
    for idx in sample_points:
        sample_images[idx].set_calculator(EMT())
        data.append(sample_images[idx])
    name = filename + "_%s_iter_%s" % (num_samples, str(iteration + 1))
    if lj:
        name = filename + "_LJ_%s_iter_%s" % (num_samples, str(iteration + 1))
    ml_lj(
        data,
        name,
        count=num_images,
        dir=dir,
        temp=temp,
        GSF=GSF,
        const_t=True,
        lj=lj,
        fine_tune=fine_tune,
        loss_fn=loss_fn,
        save_logs=save_logs,
        resample=sample_points,
    )


def iterative_sampler(
    images,
    sample_images,
    num_images,
    sample,
    filename,
    dir,
    iter,
    lj,
    temp,
    GSF,
    loss_fn,
    save_logs,
):
    for i in range(iter):
        sampler(
            images,
            sample_images,
            filename,
            dir=dir,
            num_images=num_images,
            num_samples=sample,
            iteration=i,
            temp=temp,
            GSF=GSF,
            lj=lj,
            loss_fn=loss_fn,
            fine_tune=None,
            save_logs=save_logs,
        )
        if lj:
            sample_images = ase.io.read(
                dir + filename + "_LJ_%s_iter_%s.traj" % (sample, str(i + 1)), ":"
            )
        else:
            sample_images = ase.io.read(
                dir + filename + "_%s_iter_%s.traj" % (sample, str(i + 1)), ":"
            )


# define training images
images0 = ase.io.read("../../datasets/COCu/COCu_pbc_300K.traj", ":")
images_LJ_tanh = ase.io.read(
    "MD_results/COCu/pbc_300K/tanh/paper/MLMD_COCu_pbc_300K_tanh-LJ-1.traj", ":"
)
images_ML_tanh = ase.io.read(
    "MD_results/COCu/pbc_300K/tanh/paper/MLMD_COCu_pbc_300K_tanh-1.traj", ":"
)
images_LJ_l2 = ase.io.read(
    "MD_results/COCu/pbc_300K/l2amp/paper/MLMD_COCu_pbc_300K_l2amp-LJ-2.traj", ":"
)
images_ML_l2 = ase.io.read(
    "MD_results/COCu/pbc_300K/l2amp/paper/MLMD_COCu_pbc_300K_l2amp-2.traj", ":"
)

Gs = {}
Gs["G2_etas"] = np.logspace(np.log10(0.05), np.log10(5.0), num=4)
Gs["G2_rs_s"] = [0] * 4
Gs["G4_etas"] = np.array([0.005])
Gs["G4_zetas"] = np.array([1.0])
Gs["G4_gammas"] = np.array([1.0, -1.0])
Gs["cutoff"] = 5.876798323827276

num_images = 100
num_samples = 10
jobs = []
p1 = multiprocessing.Process(target=iterative_sampler, args=(
    images0,
    images_ML_tanh,
    num_images,
    num_samples,
    "MLMD_COCu_pbc_300K_tanh",
    "MD_results/COCu/pbc_300K/tanh/paper/",
    5,
    False,
    300,
    Gs,
    "tanh",
    True,
))
jobs.append(p1)
p1.start()

p2 = multiprocessing.Process(target=iterative_sampler, args=(
    images0,
    images_LJ_tanh,
    num_images,
    num_samples,
    "MLMD_COCu_pbc_300K_tanh",
    "MD_results/COCu/pbc_300K/tanh/paper/",
    5,
    True,
    300,
    Gs,
    "tanh",
    True,
))
jobs.append(p2)
p2.start()
