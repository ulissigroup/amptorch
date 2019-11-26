"""Base AMP calculator class to be utilized similiar to other ASE calculators"""

import copy
import sys
import numpy as np
import os
from torch.utils.data import DataLoader
from amptorch.utils import Logger, hash_images
from amptorch.gaussian import Gaussian
from amptorch.data_preprocess import (
    AtomsDataset,
    TestDataset,
    factorize_data,
    collate_amp,
)
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

    def __init__(self, model, training_data=None, label=None):
        Calculator.__init__(self)

        if not os.path.exists("results/trained_models"):
            os.mkdir("results/trained_models")
        self.save_logs = model.save_logs
        label = model.label
        self.log = Logger("results/logs/"+label+".txt")
        self.log("Filename: %s" % label)
        self.model = model
        self.label = "".join(["results/trained_models/", label, ".pt"])
        self.fp_scaling = self.model.training_data.fprange
        self.target_sd = self.model.scalings[0]
        self.target_mean = self.model.scalings[1]
        self.lj = self.model.training_data.lj
        self.Gs = self.model.training_data.Gs
        self.log("Symmetry function parameters: %s" % self.Gs)
        if self.lj:
            self.fitted_params = self.model.lj_data[3]
            self.params_dict = self.model.lj_data[4]
            self.lj_model = self.model.lj_data[5]

    def train(self, overwrite=True):
        self.trained_model, value = self.model.train()
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
            images=atoms, descriptor=self.model.descriptor, Gs=self.Gs,
            fprange=self.fp_scaling, label=self.model.label)
        fp_length = dataset.fp_length()
        unique_atoms = dataset.unique()
        architecture = copy.copy(self.model.structure)
        architecture.insert(0, fp_length)
        batch_size = len(dataset)
        dataloader = DataLoader(
            dataset, batch_size, collate_fn=dataset.collate_test, shuffle=False
        )
        if properties == ['energy']:
            model = FullNN(unique_atoms, architecture, "cpu",
                    forcetraining=False)
        elif properties == ['forces']:
            model = FullNN(unique_atoms, architecture, "cpu",
                    forcetraining=True)
        model.load_state_dict(torch.load(self.label))
        model.eval()

        for inputs in dataloader:
            for element in unique_atoms:
                inputs[0][element][0] = inputs[0][element][0].requires_grad_(True)
            energy, forces = model(inputs)
        energy = (energy * self.target_sd) + self.target_mean
        energy = np.concatenate(energy.detach().numpy())
        if properties == ['forces']:
            forces = (forces * self.target_sd).detach().numpy()

        if self.lj:
            image_hash = hash_images([atoms])
            self.lj_model.neighborlist.calculate_items(image_hash)
            lj_energy, lj_forces, _ = self.lj_model.image_pred(
                atoms, self.fitted_params, self.params_dict
            )
            lj_energy = np.squeeze(lj_energy)
            energy += lj_energy
            if properties == ['forces']:
                forces += lj_forces

        self.results["energy"] = float(energy)
        self.results["forces"] = forces
