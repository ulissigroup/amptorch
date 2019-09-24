"""Base AMP calculator class to be utilized similiar to other ASE calculators"""

import copy
import numpy as np
import os
import sys
import time
from getpass import getuser
from socket import gethostname
import platform
from torch.utils.data import DataLoader
import amp
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
        self.log = Logger("results/calc-log.txt")
        self.print_header(self.log)
        self.model = model
        self.label = label
        self.fp_scaling = self.model.training_data.fprange
        self.target_sd = self.model.scalings[0]
        self.target_mean = self.model.scalings[1]
        self.parallel = self.model.training_data.parallel

    def train(self, overwrite=False):
        self.trained_model = self.model.train()
        if os.path.exists(self.label):
            self.log("'%s' parameters file already exists!" % self.label)
            if overwrite is False:
                self.log("Could not save parameters, overwrite set to false\n")
                print("Could not save! File already exists")
            else:
                self.log("Overwriting file '%s'..." % self.label)
                torch.save(self.trained_model.state_dict(), self.label)
                self.log("Model parameters successfully saved at '%s'\n" % self.label)
        else:
            self.log("Model parameters successfully saved at '%s'\n" % self.label)
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
        model = FullNN(unique_atoms, architecture, "cpu", forcetraining=True,
                activation_fn=self.model.activation_fn)
        model.load_state_dict(torch.load("amptorch.pt"))

        for batch in dataloader:
            input_data = [batch[0], len(batch[1])]
            for element in unique_atoms:
                input_data[0][element][0] = input_data[0][element][0].requires_grad_(
                    True
                )
            fp_primes = batch[2]
            energy, forces = model(input_data, fp_primes)
        energy = (energy * self.target_sd) + self.target_mean
        forces = forces * self.target_sd

        self.results["energy"] = np.concatenate(energy.detach().numpy())
        self.results["forces"] = forces.detach().numpy()

    def print_header(self, log):
        """Prints header log file"""

        log(logo)
        log("Amp - Atomistic Machine-learning Package")
        log("PyTorch Implementation")
        log(
            "Developed by Muhammed Shuaibi and Zachary Ulissi, Carnegie Mellon University"
        )
        log(
            "Original Amp developed by Andrew Peterson, Alireza Khorshidi, and others, Brown University"
        )
        log("=" * 70)
        log("User: %s" % getuser())
        log("Hostname: %s" % gethostname())
        log("Date: %s" % time.asctime())
        uname = platform.uname()
        log("Architecture: %s" % uname[4])
        log("PyTorch version: %s" % torch.__version__)
        log("Python: v{0}.{1}.{2}".format(*sys.version_info[:3]))
        log("=" * 70)

logo = """
   oo      o       o   oooooo
  o  o     oo     oo   o     o
 o    o    o o   o o   o     o
o      o   o  o o  o   o     o
oooooooo   o   o   o   oooooo
o      o   o       o   o
o      o   o       o   o
o      o   o       o   o
"""
