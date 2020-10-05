"""Base AMP calculator class to be utilized similiar to other ASE calculators"""

import os
import copy
import numpy as np
import time
from torch.utils.data import DataLoader
from amptorch.utils import Logger, hash_images, get_hash
from amptorch.skorch_model.utils import (
    make_force_header,
    make_energy_header,
    make_val_force_header,
    make_val_energy_header,
    log_results,
)
from amptorch.model import BPNN, CustomMSELoss
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

    def __init__(self, training_data, model, label, save_logs=True):
        Calculator.__init__(self)

        os.makedirs("results/", exist_ok=True)
        os.makedirs("results/trained_models", exist_ok=True)
        os.makedirs("results/logs", exist_ok=True)
        self.model = model
        self.testlabel = label
        self.label = "".join(["results/trained_models/", label, ".pt"])
        self.scalings = training_data.scalings
        self.target_ref_per_atom = self.scalings[0]
        self.delta_ref_per_atom = self.scalings[1]
        self.scale = self.scalings[2]
        self.delta = training_data.delta
        self.Gs = training_data.Gs
        self.fprange = training_data.fprange
        self.descriptor = training_data.base_descriptor
        self.cores = training_data.cores
        self.training_data = training_data
        if self.delta:
            self.params = self.training_data.delta_data[3]
            self.delta_model = self.training_data.delta_data[4]

        # TODO make utility logging function
        self.log = Logger("results/logs/{}.txt".format(label))
        if not self.delta:
            self.log(time.asctime())
            self.log("-" * 50)
        self.log("Filename: {}".format(label))
        self.log("Dataset size: {}".format(len(self.training_data)))
        self.log("Target scaling: {}".format(self.scalings))
        self.log("Symmetry function parameters:")
        for i in self.Gs.keys():
            self.log("     {}: {}".format(i, self.Gs[i]))
        self.log("Device: {:s}".format(self.model.device))
        self.log("Model: {}".format(model.module))
        self.log(
            "Architecture:\n   Input Layer - {}\n   # of Hidden Layers - {}\n   Nodes/Layer - {}".format(
                self.model.module.architecture[0],
                self.model.module.architecture[1] - 1,
                self.model.module.architecture[2],
            )
        )
        self.log("Loss Function: {}".format(self.model.criterion))
        self.log(
            "Force coefficient: {}".format(self.model.criterion__force_coefficient)
        )
        self.log("Optimizer: {}".format(self.model.optimizer))
        self.log("Learning Rate: {}".format(self.model.lr))
        self.log("Batch Size: {}".format(self.model.batch_size))
        self.log("Epochs: {}".format(self.model.max_epochs))
        self.log("Shuffle: {}".format(self.model.iterator_train__shuffle))
        if self.model.train_split != 0:
            self.log(
                "Train Split (k-fold if int, fraction if float): {}\n".format(
                    self.model.train_split.cv
                )
            )
        else:
            self.log(
                "Train Split (k-fold if int, fraction if float): {}\n".format(
                    self.model.train_split
                )
            )

    def train(self, overwrite=True):
        # TODO Store training settings (scalings, etc.) alongside model
        # parameters. Necessary for more efficient loading.
        self.model.fit(self.training_data, None)
        log_results(self.model, self.log)
        if os.path.exists(self.label):
            if overwrite is False:
                print("Could not save! File already exists")
            else:
                self.model.save_params(f_params=self.label)
        else:
            self.model.save_params(f_params=self.label)

    def load(self, filename):
        '''
        Loads calculator with previously trained model.
        ------------------

        filename: str.
            Path to trained model'''
        self.model.initialize()
        try:
            self.model.load_params(f_params=filename)
        except:
            raise Exception('File not found or trying to load a model with a different architecture than that defined')

    def calculate(self, atoms, properties, system_changes):
        Calculator.calculate(self, atoms, properties, system_changes)
        dataset = TestDataset(
            images=atoms,
            unique_atoms=self.training_data.elements,
            descriptor=self.training_data.base_descriptor,
            Gs=self.Gs,
            fprange=self.fprange,
            label=self.testlabel,
            cores=self.cores,
        )
        num_atoms = len(atoms)
        unique_atoms = dataset.unique()
        batch_size = len(dataset)
        dataloader = DataLoader(
            dataset, batch_size, collate_fn=dataset.collate_test, shuffle=False
        )
        model = self.model.module
        model.forcetraining = True
        model.load_state_dict(torch.load(self.label))
        model.eval()

        for inputs in dataloader:
            for element in unique_atoms:
                inputs[0][element][0] = inputs[0][element][0].requires_grad_(True)
            energy, forces = model(inputs)
            energy = energy*self.scale.std + self.scale.mean*num_atoms
            forces = self.scale.denorm(forces, energy=False)
        energy = np.concatenate(energy.detach().numpy())
        forces = forces.detach().numpy()

        image_hash = hash_images([atoms])
        if self.delta:
            self.delta_model.neighborlist.calculate_items(image_hash)
            delta_energy, delta_forces, _ = self.delta_model.image_pred(
                atoms, self.params
            )
            delta_energy = np.squeeze(delta_energy)
            energy += delta_energy + num_atoms*(self.target_ref_per_atom - self.delta_ref_per_atom)
            forces += delta_forces

        self.results["energy"] = float(energy)
        self.results["forces"] = forces
