"""Base AMP calculator class to be utilized similiar to other ASE calculators"""

import copy
import sys
import time
import numpy as np
import os
from torch.utils.data import DataLoader
from amp.utilities import Logger
from amp.descriptor.gaussian import Gaussian
from amp_pytorch.data_preprocess import (
    AtomsDataset,
    TestDataset,
    factorize_data,
    collate_amp,
)
from amp_pytorch.NN_model import FullNN, CustomLoss
from amp_pytorch.trainer import Trainer
from ase.calculators.calculator import Calculator, Parameters
import torch

__author__ = "Muhammed Shuaibi"
__email__ = "mshuaibi@andrew.cmu.edu"


class AMP(Calculator):
    """Atomistics Machine-Learning Potential (AMP) ASE calculator
   Parameters
   ----------
    model : object
        Class representing the regression model. Input arguments include training
        images, descriptor type, and force_coefficient. Model structure and training schemes can be
        modified directly within the class.

    label : str
        Location to save the trained model.

    """

    implemented_properties = ["energy", "forces"]

    def __init__(self, model, label="amptorch.pt"):
        Calculator.__init__(self)

        self.model = model
        self.label = label
        self.fp_scaling = self.model.training_data.fprange
        self.target_sd = self.model.scalings[0]
        self.target_mean = self.model.scalings[1]
        self.parallel = self.model.training_data.parallel
        self.lj = self.model.training_data.lj
        if self.lj:
            self.fitted_params = self.model.lj_data[3]
            self.params_dict = self.model.lj_data[4]
            self.lj_model = self.model.lj_data[5]

    def train(self, overwrite=False):

        self.trained_model = self.model.train()
        if os.path.exists(self.label):
            if overwrite is False:
                print("Could not save! File already exists")
            else:
                torch.save(self.trained_model.state_dict(), self.label)
        else:
            torch.save(self.trained_model.state_dict(), self.label)

    def calculate(self, atoms, properties, system_changes):
        Calculator.calculate(self, atoms, properties, system_changes)
        dataset = TestDataset(
            atoms, self.model.descriptor, self.fp_scaling, self.parallel
        )
        fp_length = dataset.fp_length()
        unique_atoms = dataset.unique()
        architecture = copy.copy(self.model.structure)
        architecture.insert(0, fp_length)
        batch_size = len(dataset)
        dataloader = DataLoader(
            dataset, batch_size, collate_fn=dataset.collate_test, shuffle=False
        )
        model = FullNN(unique_atoms, architecture, "cpu", forcetraining=True)
        model.load_state_dict(torch.load(self.label))

        for batch in dataloader:
            input_data = [batch[0], len(batch[1])]
            for element in unique_atoms:
                input_data[0][element][0] = input_data[0][element][0].requires_grad_(
                    True
                )
            fp_primes = batch[2]
            energy, forces = model(input_data, fp_primes)
        energy = (energy * self.target_sd) + self.target_mean
        energy = np.concatenate(energy.detach().numpy())
        forces = (forces * self.target_sd).detach().numpy()

        if self.lj:
            lj_energy, lj_forces, _ = self.lj_model.lj_pred(
                [atoms], self.fitted_params, self.params_dict
            )
            lj_energy = np.squeeze(lj_energy)
            energy += lj_energy
            forces += lj_forces

        self.results["energy"] = energy
        self.results["forces"] = forces
