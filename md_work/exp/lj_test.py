import sys
import os
import ase
import torch
from torch.nn import MSELoss
from skorch import NeuralNetRegressor
from skorch.dataset import CVSplit
from skorch.callbacks import Checkpoint, EpochScoring
from skorch.callbacks.lr_scheduler import LRScheduler
from amptorch.gaussian import SNN_Gaussian
from amp.descriptor.gaussian import Gaussian
from amptorch.model import FullNN, CustomLoss, MAELoss
from amptorch.data_preprocess import (
    AtomsDataset,
    factorize_data,
    collate_amp,
    TestDataset,
)
from md_work.md_utils import (
    md_run,
    calculate_energies,
    calculate_forces,
    time_plots,
    kde_plots,
)
from amptorch.skorch_model import AMP
from amptorch.skorch_model.utils import target_extractor, energy_score, forces_score
# from amptorch.lj_model import lj_optim
from amptorch.lj_test import lj_optim
from torch.utils.data import DataLoader
from torch.nn import init
from skorch.utils import to_numpy
import numpy as np
from ase import Atoms, units
from ase.calculators.emt import EMT
import skorch.callbacks.base
import random

# lj optimization
def lj_optimization(images, cutoff, filename, forcesonly):
    a = 12
    p0 = [1.0, 6.3535, 1.0808, 8.5357, 2.1717, 3.7575]
    e0 = [0, 0, 0]
    params_dict = {"C": [], "O": [], "Cu": []}
    lj_model = lj_optim(
        images, p0, e0, params_dict, cutoff, filename, forcesonly=forcesonly,
    )
    fitted_params = lj_model.fit()
    lj_energies, lj_forces, num_atoms = lj_model.lj_pred(
        images, p0, fitted_params, params_dict
    )

def resample_images(base_images, sample_images, num_samples):
    random.seed(3)
    sample_points = random.sample(range(1, len(sample_images)), num_samples)
    print(sample_points)
    file1 = open("resample_log.txt", "a")
    file1.write(str(sample_points)+'\n')
    file1.close()
    images = base_images.copy()
    for idx in sample_points:
        sample_images[idx].set_calculator(EMT())
        images.append(sample_images[idx])
    return images


# Define Training data
# mllj_images = ase.io.read("./lang_results/COCu_lang_2ps_EF_300K-LJ-1.traj", ":")

num_samples = 50
# resampled_lj_images = resample_images(lang_images, mllj_images, num_samples)
# define symmetry functions to be used
cutoff = 5.876798323827276  # EMT asap_cutoff: False

lj_optimization(lang_images, cutoff, 'test_1', forcesonly=True)
# lj_optimization(resampled_lj_images, cutoff, 'test_2', forcesonly=True)
