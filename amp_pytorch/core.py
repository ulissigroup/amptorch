"""core.py: Defines an AMPtorch instance that allows users to specify the
training images, GPU utilization, and validation data usage. An AMPtorch
instance contains the ability to call forth a 'train' method and subsequent
plotting methods."""

import sys
import time
import os
import torch
import torch.optim as optim
from torch.utils.data import DataLoader
import torch.nn as nn
from amp.utilities import Logger
from amp.descriptor.gaussian import Gaussian
import matplotlib.pyplot as plt
import numpy as np
from amp_pytorch.data_preprocess import AtomsDataset, factorize_data, collate_amp
from amp_pytorch.NN_model import FullNN
from amp_pytorch.trainer import train_model, pred_scaling

__author__ = "Muhammed Shuaibi"
__email__ = "mshuaibi@andrew.cmu.edu"


class AMPtorch:
    def __init__(self, datafile, device="cpu", batch_size=None, val_frac=0):

        if not os.path.exists("results"):
            os.mkdir("results")
        self.log = Logger("results/results-log.txt")
        self.log_epoch = Logger("results/epoch-log.txt")
        self.log(time.asctime())
        self.device = device
        self.filename = datafile
        self.batch_size = batch_size

        self.log("-" * 50)
        self.log("Filename: %s" % self.filename)

        self.training_data = AtomsDataset(self.filename, descriptor=Gaussian())
        test_data = AtomsDataset('../datasets/training.traj', Gaussian())
        self.training_data = [self.training_data[0], test_data[0]]
        self.unique_atoms, _, _, _, _,_ = factorize_data(self.training_data)

        self.dataset_size = len(self.training_data)
        self.validation_frac = val_frac

        if self.batch_size is None:
            self.batch_size = self.dataset_size

        if self.validation_frac != 0:
            samplers = self.training_data.create_splits(
                self.training_data, self.validation_frac
            )
            self.dataset_size = {
                "train": self.data_size - int(self.validation_frac * self.data_size),
                "val": int(self.validation_frac * self.data_size),
            }

            self.log(
                "Training Data = %d Validation Data = %d"
                % (self.dataset_size["train"], self.dataset_size["val"])
            )

            self.atoms_dataloader = {
                x: DataLoader(
                    self.training_data,
                    self.batch_size,
                    collate_fn=collate_amp,
                    sampler=samplers[x],
                )
                for x in ["train", "val"]
            }

        else:
            self.dataset_size = len(self.training_data)
            self.log("Training Data = %d" % self.dataset_size)
            self.atoms_dataloader = DataLoader(
                self.training_data,
                self.batch_size,
                collate_fn=collate_amp,
                shuffle=False,
            )

        self.model = FullNN(self.unique_atoms, self.batch_size, self.device)
        self.model = self.model.to(self.device)

    def train(
        self, criterion=nn.MSELoss(), optimizer_ft=optim.LBFGS, lr=1, rmse_criteria=2e-3
    ):
        """Trains the model under the provided optimizer conditions until
        convergence is reached as specified by the rmse_critieria."""

        criterion = criterion
        self.log("Loss Function: %s" % criterion)
        # Define the optimizer and implement any optimization settings
        optimizer_ft = optimizer_ft(self.model.parameters(), lr)
        self.log("Optimizer Info:\n %s" % optimizer_ft)

        rmse_criteria = rmse_criteria
        self.log("RMSE criteria = {}".format(rmse_criteria))
        self.log("")

        self.model = train_model(
            self.model,
            self.device,
            self.unique_atoms,
            self.dataset_size,
            criterion,
            optimizer_ft,
            self.atoms_dataloader,
            rmse_criteria,
        )
        torch.save(self.model.state_dict(), "results/best_model.pt")
        return self.model

    def parity_plot(self, model):
        """Constructs a parity plot"""

        model.eval()
        predictions = []
        targets = []
        device = self.device
        model = model.to(device)
        with torch.no_grad():
            for sample in self.atoms_dataloader:
                inputs = sample[0]
                for element in self.unique_atoms:
                    inputs[element][0] = inputs[element][0].to(device)
                targets = sample[1]
                targets = targets.to(device)
                predictions = model(inputs)
            scaled_pred = pred_scaling(
                predictions, targets, method="standardize")
            targets = targets.reshape(len(targets), 1)
            data_min = min(targets)
            data_max = max(targets)
            fig = plt.figure(figsize=(7.0, 7.0))
            ax = fig.add_subplot(111)
            targets = targets.detach().cpu().numpy()
            scaled_pred = scaled_pred.detach().cpu().numpy()
            ax.plot(targets, scaled_pred, "bo", markersize=3)
            ax.plot([data_min, data_max], [data_min, data_max], "r-", lw=0.3)
            ax.set_xlabel("ab initio energy, eV")
            ax.set_ylabel("PyTorch energy, eV")
            ax.set_title("Energies")
            fig.savefig("results/parity_plot.pdf")
        plt.show()

    def plot_residuals(self, model):
        """Plots model residuals"""

        model.eval()
        predictions = []
        targets = []
        device = self.device
        model = model.to(device)
        with torch.no_grad():
            for sample in self.atoms_dataloader:
                inputs = sample[0]
                for element in self.unique_atoms:
                    inputs[element][0] = inputs[element][0].to(device)
                targets = sample[1]
                targets = targets.to(device)
                predictions = model(inputs)
        scaled_pred = pred_scaling(predictions, targets, method="standardize")
        targets = targets.reshape(len(targets), 1)
        residuals = targets - scaled_pred
        fig = plt.figure(figsize=(7.0, 7.0))
        ax = fig.add_subplot(111)
        scaled_pred = scaled_pred.detach().cpu().numpy()
        residuals = residuals.detach().cpu().numpy()
        ax.plot(scaled_pred, residuals, "bo", markersize=3)
        ax.set_xlabel("PyTorch energy, eV")
        ax.set_ylabel("residual, eV")
        ax.set_title("Energies")
        fig.savefig("results/residuals_plot.pdf")
        plt.show()
