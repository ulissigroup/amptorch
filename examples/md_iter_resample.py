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
import multiprocessing


def ml_lj(
    IMAGES,
    filename,
    count,
    temp,
    dir="MD_results/",
    const_t=False,
    lj=False,
    fine_tune=None,
    activation_fn="l2amp",
):
    if not os.path.exists(dir):
        os.mkdir(dir)
    elements = ["C", "Cu", "O"]
    G2_etas = np.logspace(np.log10(0.05), np.log10(5.0), num=8)
    G4_etas = np.array([0.005])
    G4_zetas = np.array([1.0, 4.0])
    G4_gammas = np.array([1, -1])
    G = make_symmetry_functions(elements=elements, type="G2", etas=G2_etas)
    G = make_symmetry_functions(
        elements=elements, type="G4", etas=G4_etas, zetas=G4_zetas, gammas=G4_gammas
    )
    # lj optimization
    lj_data = None
    cutoff = 5.876798323827276
    if lj:
        # eV_kcalmol = 0.043372093
        p0 = [
            2.042,
            2.950e-11,
            -1.29265e-2,
            -1.945,
            2.6944e-3,
            6.399e-3,
            2.331,
            5.559e-2,
            2.398e-2,
        ]
        p0 = [
            2.08022879,
            1.89536258e-11,
            -5.47894512e-3,
            -2.10675310,
            1.94943321e-3,
            7.18881277e-3,
            2.29444069,
            1.00651095e-1,
            3.624737e-2,
        ]
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
    # define calculator, model, and descriptor
    # turn force training on by defining a force coefficient>0
    # define the number of cores to parallelize across for fingerprint calculations
    calc = AMP(
        model=AMPModel(
            IMAGES,
            descriptor=Gaussian(Gs=G, cutoff=cutoff),
            cores=1,
            force_coefficient=0.3,
            lj_data=lj_data,
        ),
        label="".join([filename, ".pt"]),
    )

    # calc.model.device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    calc.model.convergence = {"energy": 0.002, "force": 0.02}
    calc.model.lr = 1e-3
    # calc.model.fine_tune = fine_tune
    # calc.model.optimizer = optim.SGD
    if activation_fn == "l2amp":
        calc.model.criterion = CustomLoss
    elif activation_fn == "logcosh":
        calc.model.criterion = LogCoshLoss
    calc.model.val_frac = 0.2
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
    for step in range(count):
        dyn.run(20)
        traj.write(slab)


"""Runs multiple simulations of resampled LJ models and saves corresponding
trajectory files"""


def sampler(
    images,
    sample_images,
    filename,
    dir,
    num_images,
    num_samples,
    i,
    temp,
    lj,
    activation_fn,
    fine_tune=None,
):
    sample_points = random.sample(range(1, num_images), num_samples)
    data = [images[idx] for idx in range(num_images)]
    for idx in sample_points:
        sample_images[idx].set_calculator(EMT())
        data.append(sample_images[idx])
    name = filename + "_%s_iter_%s" % (num_samples, str(i + 1))
    if lj:
        name = filename + "_LJ_%s_iter_%s" % (num_samples, str(i + 1))
    ml_lj(
        data,
        name,
        count=num_images,
        dir=dir,
        temp=temp,
        const_t=True,
        lj=lj,
        fine_tune=fine_tune,
        activation_fn=activation_fn,
    )


def iterative_sampler(
    images, sample_images, sample, filename, dir, iter, lj, activation_fn
):
    resample_images = sample_images
    for i in range(iter):
        sampler(
            images,
            resample_images,
            filename,
            dir=dir,
            num_images=100,
            num_samples=sample,
            i=i,
            temp=300,
            lj=lj,
            activation_fn=activation_fn,
        )
        if lj:
            resample_images = ase.io.read(
                dir + filename + "_LJ_%s_iter_%s.traj" % (sample, str(i + 1)), ":"
            )
        else:
            resample_images = ase.io.read(
                dir + filename + "_%s_iter_%s.traj" % (sample, str(i + 1)), ":"
            )


# define training images
images0 = ase.io.read("../datasets/COCu/COCu_pbc_300K.traj", ":")
images_LJ = ase.io.read(
    "MD_results/COCu/pbc_300K/logcosh/paper/MLMD_COCu_pbc_300K_logcosh_LJ_2.traj", ":"
)
images_ML = ase.io.read(
    "MD_results/COCu/pbc_300K/logcosh/paper/MLMD_COCu_pbc_300K_logcosh_2.traj", ":"
)
images_LJ_l2amp = ase.io.read(
    "MD_results/COCu/pbc_300K/l2amp/paper/MLMD_COCu_pbc_300K_l2amp_LJ_3.traj", ":"
)
images_ML_l2amp = ase.io.read(
    "MD_results/COCu/pbc_300K/l2amp/paper/MLMD_COCu_pbc_300K_l2amp_3.traj", ":"
)

if __name__ == "__main__":
    iterative_sampler(
        images0,
        images_LJ,
        sample=5,
        filename="MLMD_COCu_pbc_300K_logcosh_2",
        dir="MD_results/COCu/pbc_300K/logcosh/paper/",
        iter=5,
        lj=True,
        activation_fn="logcosh",
    )

    iterative_sampler(
        images0,
        images_ML,
        sample=5,
        filename="MLMD_COCu_pbc_300K_logcosh_2",
        dir="MD_results/COCu/pbc_300K/logcosh/paper/",
        iter=5,
        lj=False,
        activation_fn="logcosh",
    )

    iterative_sampler(
        images0,
        images_LJ_l2amp,
        sample=5,
        filename="MLMD_COCu_pbc_300K_l2amp",
        dir="MD_results/COCu/pbc_300K/l2amp/paper/",
        iter=5,
        lj=True,
        activation_fn="l2amp",
    )

    iterative_sampler(
        images0,
        images_ML_l2amp,
        sample=5,
        filename="MLMD_COCu_pbc_300K_l2amp",
        dir="MD_results/COCu/pbc_300K/l2amp/paper/",
        iter=5,
        lj=False,
        activation_fn="l2amp",
    )
