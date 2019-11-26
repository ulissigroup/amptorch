"""Base AMP calculator class to be utilized similiar to other ASE calculators"""

import copy
import sys
import numpy as np
import os
from torch.utils.data import DataLoader
from amptorch.utils import Logger, hash_images
from amptorch.gaussian import Gaussian
from amptorch.data_preprocess import AtomsDataset, factorize_data, collate_amp, TestDataset
from amptorch.model import FullNN, CustomLoss
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

    def __init__(self, training_data, model, label):
        Calculator.__init__(self)

        if not os.path.exists("results/trained_models"):
            os.makedirs("results/trained_models")
        self.model = model
        self.testlabel = label
        self.label = "".join(["results/trained_models/", label, ".pt"])
        self.scalings = training_data.scalings
        self.target_sd = self.scalings[0]
        self.target_mean = self.scalings[1]
        self.lj = training_data.lj
        self.Gs = training_data.Gs
        self.fprange = training_data.fprange
        self.descriptor = training_data.base_descriptor
        self.cores = training_data.cores
        self.training_data = training_data
        if self.lj:
            self.fitted_params = self.training_data.lj_data[3]
            self.params_dict = self.training_data.lj_data[4]
            self.lj_model = self.training_data.lj_data[5]

    def train(self, overwrite=True):
        self.model.fit(self.training_data, None)
        if os.path.exists(self.label):
            if overwrite is False:
                print("Could not save! File already exists")
            else:
                self.model.save_params(f_params=self.label)
        else:
            self.model.save_params(f_params=self.label)

    def calculate(self, atoms, properties, system_changes):
        Calculator.calculate(self, atoms, properties, system_changes)
        dataset = TestDataset(
            images=atoms,
            descriptor=self.training_data.base_descriptor,
            Gs=self.Gs,
            fprange=self.fprange,
            label=self.testlabel,
            cores=self.cores,
        )
        unique_atoms = dataset.unique()
        architecture = copy.copy(self.model.module.architecture)
        batch_size = len(dataset)
        dataloader = DataLoader(
            dataset, batch_size, collate_fn=dataset.collate_test, shuffle=False
        )
        model = FullNN(unique_atoms, architecture, "cpu", forcetraining=True)
        model.load_state_dict(torch.load(self.label))
        model.eval()

        for inputs in dataloader:
            for element in unique_atoms:
                inputs[0][element][0] = inputs[0][element][0].requires_grad_(
                    True
                )
            energy, forces = model(inputs)
        energy = (energy * self.target_sd) + self.target_mean
        energy = np.concatenate(energy.detach().numpy())
        forces = (forces * self.target_sd).detach().numpy()

        if self.lj:
            image_hash = hash_images([atoms])
            self.lj_model.neighborlist.calculate_items(image_hash)
            lj_energy, lj_forces, _ = self.lj_model.image_pred(
                atoms, self.fitted_params, self.params_dict
            )
            lj_energy = np.squeeze(lj_energy)
            energy += lj_energy
            forces += lj_forces

        self.results["energy"] = float(energy)
        self.results["forces"] = forces
