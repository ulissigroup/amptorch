"""core.py: Defines an AMPtorch instance that allows users to specify the
training images, GPU utilization, and validation data usage. An AMPtorch
instance contains the ability to call forth a 'train' method and subsequent
plotting methods."""

import sys
import time
import os
import tempfile
import shutil
from torch.utils.data import DataLoader
from amp.utilities import Logger
from amp.descriptor.gaussian import Gaussian
from amp_pytorch.data_preprocess import AtomsDataset, factorize_data, collate_amp
from amp_pytorch.NN_model import FullNN, CustomLoss
from amp_pytorch.trainer import train_model, target_scaling, pred_scaling
from ase.calculators.calculator import Calculator, Parameters
import torch

__author__ = "Muhammed Shuaibi"
__email__ = "mshuaibi@andrew.cmu.edu"


class AMPCalc(Calculator):

    implemented_properties = ["energy", "forces"]

    def __init__(self, model, label="amptorch.pt"):
        Calculator.__init__(self)

        self.model = model
        self.label = label

    def train(self, overwrite=False):

        trained_model = self.model.train()
        if os.path.exists(self.label):
            if overwrite is False:
                print('Could not save! File already exists')
            else:
                torch.save(trained_model.state_dict(), self.label)
        else:
            torch.save(trained_model.state_dict(), self.label)

        return trained_model

    def calculate(self)
