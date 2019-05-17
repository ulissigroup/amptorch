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
from amp_pytorch.data_preprocess import AtomsDataset, factorize_data, collate_amp
from amp_pytorch.NN_model import FullNN
from amp_pytorch.trainer import train_model, target_scaling, pred_scaling

__author__ = "Muhammed Shuaibi"
__email__ = "mshuaibi@andrew.cmu.edu"


class AMPtorch:
    def __init__(
        self, datafile, device="cpu", batch_size=None, structure=[2, 2], val_frac=0
    ):

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
        # self.training_data = [self.training_data[0], self.training_data[1]]
        # self.unique_atoms = ['O','H']
        # self.fp_length = 20
        self.unique_atoms = self.training_data.unique()
        self.fp_length = self.training_data.fp_length()

        self.dataset_size = len(self.training_data)
        self.validation_frac = val_frac

        if self.batch_size is None:
            self.batch_size = self.dataset_size

        if self.validation_frac != 0:
            samplers = self.training_data.create_splits(
                self.training_data, self.validation_frac
            )
            self.dataset_size = {
                "train": self.dataset_size
                - int(self.validation_frac * self.dataset_size),
                "val": int(self.validation_frac * self.dataset_size),
            }

            self.log(
                "Training Data = %d Validation Data = %d"
                % (self.dataset_size["train"], self.dataset_size["val"])
            )

            if self.batch_size is None:
                self.batch_size = {x: self.dataset_size[x] for x in ["train", "val"]}

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
            self.log("Training Data = %d" % self.dataset_size)
            self.atoms_dataloader = DataLoader(
                self.training_data,
                self.batch_size,
                collate_fn=collate_amp,
                shuffle=False,
            )
        self.architecture = structure
        self.architecture.insert(0, self.fp_length)
        self.model = FullNN(self.unique_atoms, self.architecture, self.device)
        self.model = self.model.to(self.device)

    def train(
        self,
        criterion=nn.MSELoss(),
        optimizer_ft=optim.LBFGS,
        scheduler=None,
        lr=1,
        rmse_criteria=2e-3,
    ):
        """Trains the model under the provided optimizer conditions until
        convergence is reached as specified by the rmse_critieria."""

        criterion = criterion
        self.log("Loss Function: %s" % criterion)
        # Define the optimizer and implement any optimization settings
        optimizer_ft = optimizer_ft(self.model.parameters(), lr)
        self.log("Optimizer Info:\n %s" % optimizer_ft)

        if scheduler:
            scheduler = scheduler(optimizer_ft, step_size=7, gamma=0.1)
        self.log("Scheduler Info: \n %s" % scheduler)

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
            scheduler,
            self.atoms_dataloader,
            rmse_criteria,
        )
        torch.save(self.model.state_dict(), "results/best_model.pt")
        return self.model

    def parity_plot(self, model):
        """Constructs a parity plot"""

        model.eval()
        targets = []
        device = self.device
        model = model.to(device)
        for sample in self.atoms_dataloader:
            inputs = sample[0]
            fp_primes = sample[3]
            for element in self.unique_atoms:
                inputs[element][0] = inputs[element][0].to(device).requires_grad_(True)
            targets = sample[1]
            force_targets = sample[4]
            targets = targets.to(device)
            force_targets = force_targets.to(device)
            inputs = [inputs, len(targets)]
            energy_pred, force_pred = model(inputs, fp_primes)
        scaled_pred = pred_scaling(energy_pred, targets, method="standardize")
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

    def parity_plot_forces(self, model):

        model.eval()
        targets = []
        device = self.device
        model = model.to(device)
        for sample in self.atoms_dataloader:
            inputs = sample[0]
            fp_primes = sample[3]
            for element in self.unique_atoms:
                inputs[element][0] = inputs[element][0].to(device).requires_grad_(True)
            targets = sample[1]
            force_targets = sample[4]
            targets = targets.to(device)
            force_targets = force_targets.to(device)
            inputs = [inputs, len(targets)]
            energy_pred, force_pred = model(inputs, fp_primes)
        _, scaling_factor = target_scaling(targets, method="standardize")
        force_pred = force_pred.reshape(-1, 1)
        force_pred *= scaling_factor
        force_targets = force_targets.reshape(-1, 1)
        data_min = min(force_targets)
        data_max = max(force_targets)
        fig = plt.figure(figsize=(7.0, 7.0))
        ax = fig.add_subplot(111)
        force_targets = force_targets.detach().cpu().numpy()
        force_pred = force_pred.detach().cpu().numpy()
        ax.plot(force_targets, force_pred, "bo", markersize=3)
        ax.plot([data_min, data_max], [data_min, data_max], "r-", lw=0.3)
        ax.set_xlabel("ab initio forces, eV/A")
        ax.set_ylabel("PyTorch forces, eV/A")
        ax.set_title("Forces")
        fig.savefig("results/parity_plot_forces.pdf")
        plt.show()

    def plot_residuals(self, model):
        """Plots model residuals"""

        model.eval()
        targets = []
        device = self.device
        model = model.to(device)
        for sample in self.atoms_dataloader:
            inputs = sample[0]
            fp_primes = sample[3]
            for element in self.unique_atoms:
                inputs[element][0] = inputs[element][0].to(device).requires_grad_(True)
            targets = sample[1]
            force_targets = sample[4]
            targets = targets.to(device)
            force_targets = force_targets.to(device)
            inputs = [inputs, len(targets)]
            energy_pred, force_pred = model(inputs, fp_primes)
        _, scaling_factor = target_scaling(targets, method="standardize")
        force_pred = force_pred.reshape(-1, 1)
        force_pred *= scaling_factor
        force_targets = force_targets.reshape(-1, 1)
        residual = abs(force_targets-force_pred)
        fig = plt.figure(figsize=(7.0, 7.0))
        ax = fig.add_subplot(111)
        force_targets = force_targets.detach().cpu().numpy()
        force_pred = force_pred.detach().cpu().numpy()
        residual = residual.detach().cpu().numpy()
        ax.plot(force_targets, residual, "bo", markersize=3)
        ax.set_xlabel("DFT force, eV/A")
        ax.set_ylabel("|DFT force - PyTorch force|, eV/A")
        ax.set_title("Forces")
        fig.savefig("results/residual_plot_forces.pdf")
        plt.show()
