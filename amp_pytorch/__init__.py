import sys
import time
import torch
import torch.nn as nn
import numpy as np
from torch.utils.data import DataLoader
from data_preprocess import AtomsDataset, factorize_data, collate_amp
from amp.utilities import Logger
from amp.descriptor.gaussian import Gaussian
from NN_model import FullNN
from trainer import train_model
import torch.optim as optim
import matplotlib.pyplot as plt


class AMPtorch():

    def __init__(self, datafile='../benchmark_dataset/water.extxyz',
                 device='cpu', val_frac=0):

        self.log = Logger("../benchmark_results/results-log.txt")
        self.log_epoch = Logger("../benchmark_results/epoch-log.txt")
        self.log(time.asctime())
        self.device = device
        self.filename = datafile

        self.log("-" * 50)
        self.log("Filename: %s" % self.filename)

        self.training_data = AtomsDataset(self.filename, descriptor=Gaussian())
        self.unique_atoms, _, _, _ = factorize_data(self.training_data)
        # n_unique_atoms = len(unique_atoms)

        self.batch_size = len(self.training_data)
        self.log("Batch Size = %d" % self.batch_size)
        self.validation_frac = val_frac

        if self.validation_frac != 0:
            samplers = self.training_data.create_splits(
                self.training_data, self.validation_frac)
            dataset_size = {
                "train": (1.0 - self.validation_frac) * len(self.training_data),
                "val": self.validation_frac * len(self.training_data),
            }

            self.log(
                "Training Data = %d Validation Data = %d"
                % (dataset_size["train"], dataset_size["val"])
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
                self.training_data, self.batch_size, collate_fn=collate_amp,
                shuffle=False
            )

        self.model = FullNN(self.unique_atoms, self.batch_size, device='cpu')
        self.model = self.model.to(self.device)

    def train(self, criterion=nn.MSELoss(), optimizer_ft=optim.LBFGS,
              rmse_criteria=2e-3):

        self.criterion = criterion
        self.log("Loss Function: %s" % self.criterion)
        # Define the optimizer and implement any optimization settings
        self.optimizer_ft = optimizer_ft(self.model.parameters(), 1)
        self.log("Optimizer Info:\n %s" % self.optimizer_ft)

        self.rmse_criteria = rmse_criteria
        self.log("RMSE criteria = {}".format(self.rmse_criteria))
        self.log("")

        self.model = train_model(
            self.model,
            self.device,
            self.unique_atoms,
            self.dataset_size,
            self.criterion,
            self.optimizer_ft,
            self.atoms_dataloader,
            rmse_criteria
        )
        torch.save(self.model.state_dict(),
                   "../benchmark_results/benchmark_model.pt")


test = AMPtorch()
test.train()


def parity_plot(training_data):
    loader = DataLoader(training_data, 400,
                        collate_fn=collate_amp, shuffle=False)
    model = FullNN(unique_atoms, 400)
    model.load_state_dict(torch.load(
        "../benchmark_results/benchmark_model.pt"))
    model.eval()
    predictions = []
    targets = []
    # device='cuda:0'
    device = "cpu"
    model = model.to(device)
    with torch.no_grad():
        for sample in loader:
            inputs = sample[0]
            for element in unique_atoms:
                inputs[element][0] = inputs[element][0].to(device)
            targets = sample[1]
            targets = targets.to(device)
            predictions = model(inputs)
        data_max = max(targets)
        data_min = min(targets)
        data_mean = torch.mean(targets)
        data_sd = torch.std(targets, dim=0)
        scale = (predictions * data_sd) + data_mean
        # scale=(predictions*(data_max-data_min))+data_min
        targets = targets.reshape(len(targets), 1)
        # scaled_pred=scaled_pred.reshape(len(targets),1)
        crit = nn.MSELoss()
        loss = crit(scale, targets)
        loss = loss / len(unique_atoms) ** 2
        loss = loss.detach().numpy()
        RMSE = np.sqrt(loss)
        print RMSE
        fig = plt.figure(figsize=(7.0, 7.0))
        ax = fig.add_subplot(111)
        targets = targets.detach().numpy()
        scale = scale.detach().numpy()
        ax.plot(targets, scale, "bo", markersize=3)
        ax.plot([data_min, data_max], [data_min, data_max], "r-", lw=0.3)
        ax.set_xlabel("ab initio energy, eV")
        ax.set_ylabel("PyTorch energy, eV")
        ax.set_title("Energies")
        fig.savefig("../benchmark_results/Plots/PyTorch_Prelims.pdf")
    plt.show()


def plot_hist(training_data):
    loader = DataLoader(
        training_data, 1, collate_fn=collate_amp, shuffle=False)
    model = FullNN(unique_atoms, 1)
    model.load_state_dict(torch.load("benchmark_results/benchmark_model.pt"))
    model.eval()
    predictions = []
    scaled_pred = []
    targets = []
    residuals = []
    # device='cuda:0'
    device = "cpu"
    model = model.to(device)
    for sample in loader:
        inputs = sample[0]
        for element in unique_atoms:
            inputs[element][0] = inputs[element][0].to(device)
        target = sample[1]
        target = target.to(device)
        prediction = model(inputs)
        predictions.append(prediction)
        targets.append(target)
    # data_max = max(targets)
    # data_min = min(targets)
    targets = torch.stack(targets)
    data_mean = torch.mean(targets)
    data_sd = torch.std(targets, dim=0)
    for index, value in enumerate(predictions):
        # scaled_value=(value*(data_max-data_min))+data_min
        scaled_value = (value * data_sd) + data_mean
        scaled_pred.append(scaled_value)
        residual = targets[index] - scaled_value
        residuals.append(residual)
    fig = plt.figure(figsize=(7.0, 7.0))
    ax = fig.add_subplot(111)
    ax.plot(scaled_pred, residuals, "bo", markersize=3)
    ax.set_xlabel("PyTorch energy, eV")
    ax.set_ylabel("residual, eV")
    ax.set_title("Energies")
    fig.savefig("benchmark_results/Plots/PyTorch_Residuals.pdf")
    # plt.show()
