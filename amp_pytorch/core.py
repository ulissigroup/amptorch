"""
core.py: AMPModel defines the core backbone to the training scheme. Model
parameters are defined, a train method carrys out the training, and plotting
methods allow the visualization of the results.
"""

import time
import sys
import os
import copy
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


class AMPModel:
    """Model class used to define the Neural Network architecture and regression
    training scheme

    Parameters
    ----------
    datafile : trajectory file, list, or database
        Training data to be utilized for the regression model.
    device : str
        Hardware to be utilized for training - CPU or GPU.
        default: 'cpu'
    cores : int
        Specify the number of cores to use for parallelization of fingerprint
        calculations. NOTE - To be replaced with PyTorch's parallelization
        scheme upon fingerprint implementation
    envommand: string
        For parallel processing across nodes, a command can be supplied here to
        load the appropriate environment before starting workers.
        default=None. NOTE - See Above.
    batch_size: int
        Number of images in a training batch. None represents no batching,
        treating the entire dataset as one batch. default: None
    structure: list
        Neural Network architecture. First index represents number of layers,
        including the output layer. Second index represents the number of nodes
        in each hidden layer. i.e. [3,5] = 3 layeres (2 hidden layers, 1 output
        layer) and 5 nodes in each hidden layer. default: [3,5]
    val_frac: float
        Proportion of dataset to be used as a validation test set for training
        purposes. default: 0
    descriptor: object
        Descriptor to be utilized to calculate fingerprints and
        fingerprintprimes. default: Gaussian()
    force_coefficient: float
        Define the force coefficient to be utilized in the loss function. A
        coefficient > 0 indicates force training is turned on.
        default: 0
    criterion: object
        Specify the loss function to be optimized.
        default: CustomLoss
    optimizer: object
        Define the training optimizer to be utilized for the regression.
        default: optim.LBFGS
    scheduler: object
        Specify whether a learning rate decay scheme is to be utilized during
        training.
        default: None
    lr: float
        Define the model learning rate. default:1
    criteria: dict
        Define the training convergence criteria.
        default: {'energy':0.02, "force":0.02}

    """

    def __init__(
        self,
        datafile,
        device="cpu",
        cores=1,
        envcommand=None,
        batch_size=None,
        structure=[3, 5],
        val_frac=0,
        descriptor=Gaussian(),
        force_coefficient=0,
        criterion=CustomLoss,
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
        self.batch_size = batch_size
        self.device = device
        self.batch_size = batch_size
        self.structure = structure
        self.val_frac = val_frac
        self.descriptor = descriptor
        self.force_coefficient = force_coefficient
        self.lossfunction = criterion
        self.optimizer = optimizer
        self.scheduler = scheduler
        self.lr = lr
        self.convergence = criteria

        self.forcetraining = False
        if force_coefficient > 0:
            self.forcetraining = True

        self.training_data = AtomsDataset(
            self.filename, descriptor=self.descriptor, cores=cores,
            forcetraining=self.forcetraining, envcommand=envcommand
        )
        self.scalings = self.training_data.scalings()
        self.sd_scaling = self.scalings[0]
        self.mean_scaling = self.scalings[1]

        self.log("-" * 50)
        self.log("Filename: %s" % self.filename)

    def train(self):
        """Trains the model under the provided optimizer conditions until
        convergence is reached as specified by the rmse_critieria."""

        training_data = self.training_data
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
        model = FullNN(
            self.unique_atoms, architecture, self.device, self.forcetraining
        ).to(self.device)

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
            self.lossfunction(force_coefficient=self.force_coefficient),
            optimizer,
            self.scheduler,
            self.atoms_dataloader,
            self.convergence,
            self.scalings,
        )

        self.trained_model = self.trainer.train_model()
        return self.trained_model

    def parity_plot(self, data="energy"):
        """Constructs a parity plot"""
        model = self.trained_model
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
            targets = targets.to(device)
            force_targets = sample[4]
            if self.forcetraining:
                force_targets = force_targets.to(device)
            inputs = [inputs, len(targets)]
            energy_pred, force_pred = model(inputs, fp_primes)
        fig = plt.figure(figsize=(7.0, 7.0))
        ax = fig.add_subplot(111)
        if data == "energy":
            scaled_pred = (energy_pred * self.sd_scaling) + self.mean_scaling
            energy_targets = targets.reshape(-1, 1)
            energy_min = min(targets)
            energy_max = max(targets)
            energy_targets = energy_targets.detach().cpu().numpy()
            scaled_pred = scaled_pred.detach().cpu().numpy()
            ax.plot(energy_targets, scaled_pred, "bo", markersize=3)
            ax.plot([energy_min, energy_max], [energy_min, energy_max], "r-", lw=0.3)
            ax.set_xlabel("ab initio energy, eV")
            ax.set_ylabel("PyTorch energy, eV")
            ax.set_title("Energies")
            fig.savefig("results/parity_plot.pdf")
        if data == "forces":
            force_pred = force_pred.reshape(-1, 1)
            force_pred *= self.sd_scaling
            force_targets = force_targets.reshape(-1, 1)
            force_min = min(force_targets)
            force_max = max(force_targets)
            force_pred = force_pred.detach().cpu().numpy()
            force_targets = force_targets.detach().cpu().numpy()
            ax.plot(force_targets, force_pred, "bo", markersize=3)
            ax.plot([force_min, force_max], [force_min, force_max], "r-", lw=0.3)
            ax.set_xlabel("ab initio forces, eV/A")
            ax.set_ylabel("PyTorch forces, eV/A")
            ax.set_title("Forces")
            fig.savefig("results/parity_plot_forces.pdf")
        plt.show()

    def plot_residuals(self, data="energy"):
        """Constructs a residual plot"""
        model = self.trained_model
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
            if self.forcetraining:
                force_targets = force_targets.to(device)
            inputs = [inputs, len(targets)]
            energy_pred, force_pred = model(inputs, fp_primes)
        fig = plt.figure(figsize=(7.0, 7.0))
        ax = fig.add_subplot(111)
        if data == "energy":
            scaled_pred = (energy_pred * self.sd_scaling) + self.mean_scaling
            energy_targets = targets.reshape(-1, 1)
            energy_targets = energy_targets.detach().cpu().numpy()
            scaled_pred = scaled_pred.detach().cpu().numpy()
            residual = abs(energy_targets - scaled_pred)
            ax.plot(energy_targets, residual, "bo", markersize=3)
            ax.set_xlabel("ab initio energy, eV")
            ax.set_ylabel("|ab initio energy - PyTorch energy|, eV")
            ax.set_title("Energy Residuals")
            fig.savefig("results/residual_plot_e.pdf")
        if data == "forces":
            force_pred = force_pred.reshape(-1, 1)
            force_pred *= self.sd_scaling
            force_targets = force_targets.reshape(-1, 1)
            force_pred = force_pred.detach().cpu().numpy()
            force_targets = force_targets.detach().cpu().numpy()
            residual = abs(force_targets - force_pred)
            ax.plot(force_targets, residual, "bo", markersize=3)
            ax.set_xlabel("ab initio force, eV/A")
            ax.set_ylabel("|ab initio force - PyTorch force|, eV/A")
            ax.set_title("Force Residuals")
            fig.savefig("results/residual_plot_f.pdf")
        plt.show()
