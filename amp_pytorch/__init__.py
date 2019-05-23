"""core.py: Defines an AMPtorch instance that allows users to specify the
training images, GPU utilization, and validation data usage. An AMPtorch
instance contains the ability to call forth a 'train' method and subsequent
plotting methods."""

import sys
import time
from torch.utils.data import DataLoader
from amp.utilities import Logger
from amp.descriptor.gaussian import Gaussian
from amp_pytorch.data_preprocess import AtomsDataset, factorize_data, collate_amp
from amp_pytorch.NN_model import FullNN, CustomLoss
from amp_pytorch.trainer import train_model, target_scaling, pred_scaling
from ase.calculators.calculator import Calculator, Parameters
import torch.nn as nn
import torch.optim as optim

__author__ = "Muhammed Shuaibi"
__email__ = "mshuaibi@andrew.cmu.edu"


class AMPCalc(Calculator):

    implemented_properties = ["energy", "forces"]

    def __init__(self, model, label='amptorch'):
        Calculator.__init__(self)

        self.model = model
        self.label = label

    def train(
        self,
        overwrite=False
    ):

        trained_model = self.model.train()
        return trained_model
