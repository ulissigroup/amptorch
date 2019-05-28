"""core.py: Defines an AMPtorch instance that allows users to specify the
training images, GPU utilization, and validation data usage. An AMPtorch
instance contains the ability to call forth a 'train' method and subsequent
plotting methods."""

import sys
import time
import os
import copy
import torch
import torch.optim as optim
from torch.utils.data import DataLoader
import torch.nn as nn
from amp.utilities import Logger
from amp.descriptor.gaussian import Gaussian
import matplotlib.pyplot as plt
from amp_pytorch.data_preprocess import AtomsDataset, factorize_data, collate_amp
from amp_pytorch.NN_model import FullNN, CustomLoss
from amp_pytorch.trainer import Trainer

__author__ = "Muhammed Shuaibi"
__email__ = "mshuaibi@andrew.cmu.edu"


class AMPtorch:
    def __init__(
        self,
        datafile,
        device="cpu",
        batch_size=None,
        structure=[3, 5],
        val_frac=0,
        descriptor=Gaussian(),
        criterion=CustomLoss(force_coefficient=0),
        optimizer=optim.LBFGS,
        scheduler=None,
        lr=1,
        criteria={"energy": 0.02, "force": 0.02},
    ):
        if not os.path.exists("results"):
            os.mkdir("results")
        self.log = Logger("results/results-log.txt")
        self.log_epoch = Logger("results/epoch-log.txt")
        self.log(time.asctime())

        self.filename = datafile
        self.device = device
        self.batch_size = batch_size
        self.structure = structure
        self.val_frac = val_frac
        self.descriptor = descriptor
        self.lossfunction = criterion
        self.optimizer = optimizer
        self.scheduler = scheduler
        self.lr = lr
        self.convergence = criteria

        self.training_data = AtomsDataset(self.filename, descriptor=self.descriptor)

        self.log("-" * 50)
        self.log("Filename: %s" % self.filename)

    def train(self):
        """Trains the model under the provided optimizer conditions until
        convergence is reached as specified by the rmse_critieria."""

        training_data = self.training_data
        # print('dataset: %s' %(time.time()-dataset_timer))
        self.unique_atoms = training_data.unique()
        fp_length = training_data.fp_length()

        dataset_size = len(training_data)

        if self.batch_size is None:
            self.batch_size = dataset_size

        if self.val_frac != 0:
            samplers = training_data.create_splits(training_data, self.val_frac)
            dataset_size = {
                "train": dataset_size - int(self.val_frac * dataset_size),
                "val": int(self.val_frac * dataset_size),
            }

            self.log(
                "Training Data = %d Validation Data = %d"
                % (dataset_size["train"], dataset_size["val"])
            )

            if self.batch_size is None:
                self.batch_size = {x: dataset_size[x] for x in ["train", "val"]}

            self.atoms_dataloader = {
                x: DataLoader(
                    training_data,
                    self.batch_size,
                    collate_fn=collate_amp,
                    sampler=samplers[x],
                )
                for x in ["train", "val"]
            }

        else:
            self.log("Training Data = %d" % dataset_size)
            self.atoms_dataloader = DataLoader(
                training_data, self.batch_size, collate_fn=collate_amp, shuffle=False
            )
        architecture = copy.copy(self.structure)
        architecture.insert(0, fp_length)
        model = FullNN(self.unique_atoms, architecture, self.device).to(self.device)

        self.log("Loss Function: %s" % self.lossfunction)
        # Define the optimizer and implement any optimization settings
        optimizer = self.optimizer(model.parameters(), self.lr)
        self.log("Optimizer Info:\n %s" % optimizer)

        if self.scheduler:
            self.scheduler = self.scheduler(optimizer, step_size=7, gamma=0.1)
        self.log("Scheduler Info: \n %s" % self.scheduler)

        self.log("RMSE criteria = {}\n".format(self.convergence))

        self.trainer = Trainer(
            model,
            self.device,
            self.unique_atoms,
            dataset_size,
            self.lossfunction,
            optimizer,
            self.scheduler,
            self.atoms_dataloader,
            self.convergence,
        )

        trained_model = self.trainer.train_model()
        return trained_model

    def parity_plot(self, model):
        """Constructs an energy parity plot"""
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
        """Constructs a forces parity plot"""
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
        """Plots force residuals"""
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
        residual = abs(force_targets - force_pred)
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
