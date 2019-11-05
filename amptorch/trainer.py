"""Trainer class used to train a specified model in accordance with the trainer
arguments"""

from amptorch.NN_model import CustomLoss, weighted_mse_loss
from .utils import (
    Logger,
    make_energy_header,
    make_force_header,
    make_val_energy_header,
    make_val_force_header,
    log_force_results,
    log_energy_results,
)
import copy
import time
import os
import sys
import torch.nn as nn
import numpy as np
import torch
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

"""Trainer class used to train a specified model in accordance with the trainer
arguments"""

__author__ = "Muhammed Shuaibi"
__email__ = "mshuaibi@andrew.cmu.edu"



def plot_grad_flow(named_parameters):
    '''Plots the gradients flowing through different layers in the net during training.
    Can be used for checking for possible gradient vanishing / exploding problems.
    Usage: Plug this function in Trainer class after loss.backwards() as
    "plot_grad_flow(self.model.named_parameters())" to visualize the gradient flow'''
    ave_grads = []
    max_grads = []
    layers = []
    for n, p in named_parameters:
        if(p.requires_grad) and ("bias" not in n):
            layers.append(n)
            ave_grads.append(p.grad.abs().mean())
            max_grads.append(p.grad.abs().max())
    plt.bar(np.arange(len(max_grads)), max_grads,
            alpha=0.1, lw=1, color="c")
    plt.bar(np.arange(len(max_grads)), ave_grads,
            alpha=0.1, lw=1, color="b")
    plt.hlines(0, 0, len(ave_grads)+1, lw=2, color="k")
    plt.xticks(range(0, len(ave_grads), 1), layers, rotation="vertical")
    plt.xlim(left=0, right=len(ave_grads))
    # zoom in on the lower gradient regions
    # plt.ylim(bottom=-0.001, top=0.02)
    plt.xlabel("Layers")
    plt.ylabel("average gradient")
    plt.title("Gradient flow")
    plt.grid(True)
    plt.legend([Line2D([0], [0], color="c", lw=4),
                Line2D([0], [0], color="b", lw=4),
                Line2D([0], [0], color="k", lw=4)], ['max-gradient', 'mean-gradient', 'zero-gradient'])
    plt.show()


class Trainer:

    """Class utilized to train a given model through a training loop method.

    Parameters
    ----------

    model: object
        NeuralNetwork model to be trained
    device: str
        Hardware to be utilized for training - CPU or GPU
    dataset_size: int
        Size of the training dataset
    criterion: object
        Loss function to be optimized.
    Optimizer: object
        Training optimizer to be utilized for the regression
    scheduler: object
        Learning rate decay scheme to be utilized in training
    atoms_dataloader: object
        PyTorch DataLoader object that tells the trainer how to load data for
        training
    convergence_criteria: dict
        Training convergence criteria
    scalings: list
        Scalings applied to the dataset to normalize the energy targets.
        Training scalings will be utilized for calculating predicted
        properities.
    """

    def __init__(
        self,
        model,
        device,
        dataset_size,
        criterion,
        optimizer,
        scheduler,
        atoms_dataloader,
        convergence_criteria,
        scalings,
        label,
        weighted=False,
        save_interval=1000,
        maxtime = None
    ):
        self.model = model
        self.device = device
        self.dataset_size = dataset_size
        self.criterion = criterion
        self.optimizer = optimizer
        self.scheduler = scheduler
        self.atoms_dataloader = atoms_dataloader
        self.convergence_criteria = convergence_criteria
        self.sd_scaling = scalings[0]
        self.mean_scaling = scalings[1]
        self.label = label
        self.weighted = weighted
        self.save_interval = save_interval
        self.maxtime=maxtime
        self.epochs = convergence_criteria["epochs"]
        self.early_stop = convergence_criteria["early_stop"]

    def train_model(self):
        "Training loop"
        forcetraining = False
        if self.criterion.alpha > 0:
            forcetraining = True
        best_train_force_loss = 1e8
        best_train_energy_loss = 1e8
        # dummy variables to track each epochs rmse
        previous_force_rmse = 1e8
        previous_energy_rmse = 1e8
        log = Logger("results/logs/" + self.label + ".txt")
        log_epoch = Logger("results/logs/epochs/" + self.label + "-calc.txt")

        plot_energy_loss = {"train": [], "val": []}
        if forcetraining:
            plot_force_loss = {"train": [], "val": []}

        if isinstance(self.atoms_dataloader, dict):
            validation = True
            best_val_force_loss = 1e8
            best_val_energy_loss = 1e8
            if forcetraining:
                make_val_force_header(log_epoch)
            else:
                make_val_energy_header(log_epoch)
        else:
            validation = False
            if forcetraining:
                make_force_header(log_epoch)
            else:
                make_energy_header(log_epoch)

        since = time.time()
        print("Training Initiated!")
        self.epochs -= 1
        early_stop = self.early_stop
        epoch = 0
        convergence = False
        while not convergence:

            print('epoch: {}'.format(epoch))
            if not os.path.isdir('results/models'):
                os.mkdir('results/models')
            if not (epoch % self.save_interval):
                torch.save(self.model.state_dict(),
                           'results/models/epoch{}.pt'.format(epoch))

            if validation:
                for phase in ["train", "val"]:

                    if phase == "train":
                        if self.scheduler:
                            self.scheduler.step()
                        self.model.train()
                    else:
                        self.model.eval()

                    energy_mse = 0.0
                    force_mse = "N/A"
                    if forcetraining:
                        force_mse = 0.0

                    for data_sample in self.atoms_dataloader[phase]:
                        unique_atoms = data_sample[3]
                        input_data = [data_sample[0], len(data_sample[1]), unique_atoms]
                        target = data_sample[1].requires_grad_(False)
                        batch_size = len(target)
                        target = target.reshape(batch_size, 1).to(self.device)
                        scaled_target = (target - self.mean_scaling) / self.sd_scaling
                        num_of_atoms = (
                            data_sample[2].reshape(batch_size, 1).to(self.device)
                        )
                        if self.weighted:
                            weights = data_sample[-1].to(self.device)
                        for element in unique_atoms:
                            input_data[0][element][0] = (
                                input_data[0][element][0]
                                .to(self.device)
                                .requires_grad_(True)
                            )
                        scaled_target = scaled_target.to(self.device)
                        fp_primes = data_sample[4]

                        if forcetraining:
                            fp_primes = fp_primes.to(self.device)
                            image_forces = data_sample[5].to(self.device)
                            scaled_forces = image_forces / self.sd_scaling

                        def closure():
                            self.optimizer.zero_grad()
                            if forcetraining:
                                energy_pred, force_pred = self.model(
                                    input_data, fp_primes
                                )
                                if self.weighted:
                                    loss = self.criterion(
                                        energy_pred,
                                        scaled_target,
                                        num_of_atoms,
                                        force_pred,
                                        scaled_forces,
                                        self.model,
                                        weights=weights,
                                    )
                                else:
                                    loss = self.criterion(
                                        energy_pred,
                                        scaled_target,
                                        num_of_atoms,
                                        force_pred,
                                        scaled_forces,
                                        self.model,
                                    )
                            else:
                                energy_pred, _ = self.model(input_data)
                                if self.weighted:
                                    loss = self.criterion(
                                        energy_pred, scaled_target, num_of_atoms,
                                        model=self.model,
                                        weights=weights
                                    )
                                else:
                                    loss = self.criterion(
                                        energy_pred, scaled_target, num_of_atoms,
                                        model=self.model
                                    )
                            loss.backward()
                            return loss

                        energy_pred, force_pred = self.model(input_data, fp_primes)
                        mse_loss = nn.MSELoss(reduction="sum")
                        raw_preds = (energy_pred * self.sd_scaling) + self.mean_scaling
                        raw_preds_per_atom = torch.div(raw_preds, num_of_atoms)
                        target_per_atom = torch.div(target, num_of_atoms)
                        energy_loss = mse_loss(raw_preds_per_atom, target_per_atom)
                        w_energy_loss = weighted_mse_loss(raw_preds_per_atom,
                                                          target_per_atom,
                                            torch.FloatTensor([1.] * len(target_per_atom)))
                        energy_mse += torch.tensor(energy_loss.item())

                        if forcetraining:
                            force_pred = force_pred * self.sd_scaling
                            num_atoms_force = torch.cat(
                                [idx.repeat(int(idx)) for idx in num_of_atoms]
                            )
                            num_atoms_force = torch.sqrt(
                                num_atoms_force.reshape(len(num_atoms_force), 1)
                            )
                            force_pred_per_atom = torch.div(force_pred, num_atoms_force)
                            force_targets_per_atom = torch.div(
                                image_forces, num_atoms_force
                            )
                            force_loss = mse_loss(
                                force_pred_per_atom, force_targets_per_atom
                            )
                            if torch.isnan(force_loss):
                                early_stop = True
                            w_f = weighted_mse_loss(force_pred_per_atom, force_targets_per_atom,
                                            torch.FloatTensor([1.] * len(force_pred_per_atom)))
                            force_mse += torch.tensor(force_loss.item())

                        if phase == "train":
                            loss = self.optimizer.step(closure)
                        now = time.asctime()

                    energy_mse /= self.dataset_size[phase]
                    energy_rmse = torch.sqrt(energy_mse)
                    if torch.isnan(energy_rmse):
                        early_stop = True
                    plot_energy_loss[phase].append(energy_rmse)
                    print("%s energy loss: %f" % (phase, energy_rmse))
                    if forcetraining:
                        force_mse /= self.dataset_size[phase]
                        force_rmse = torch.sqrt(force_mse)
                        plot_force_loss[phase].append(force_rmse)
                        print("%s force loss: %f" % (phase, force_rmse))
                        if phase == "train":
                            log_force_results(
                                log_epoch,
                                epoch,
                                now,
                                loss,
                                energy_rmse,
                                force_rmse,
                                phase,
                            )
                        else:
                            log_force_results(
                                log_epoch,
                                epoch,
                                now,
                                "",
                                energy_rmse,
                                force_rmse,
                                phase,
                            )
                        if phase == "train":
                            # early stop when training force error stagnates
                            if abs(force_rmse - previous_force_rmse) <= 1e-5 and \
                                abs(energy_rmse - previous_energy_rmse) <= 1e-5:
                                print('Early Stop Triggered')
                                early_stop = True
                            time_elapsed = time.time() - since
                            if self.maxtime is not None:
                                if time_elapsed >= self.maxtime:
                                    print('Maximum Time Reached') 
                                    early_stop = True

                            previous_force_rmse = force_rmse
                        elif phase == "val":
                            if force_rmse < best_val_force_loss:
                                best_val_energy_loss = energy_rmse
                                best_val_force_loss = force_rmse
                                best_model_wts = copy.deepcopy(self.model.state_dict())
                            energy_convergence = (
                                best_val_force_loss
                                <= self.convergence_criteria["energy"]
                            )
                            force_convergence = (
                                best_val_force_loss
                                <= self.convergence_criteria["force"]
                            )
                            convergence = (
                                (energy_convergence and force_convergence)
                                or (epoch >= self.epochs)
                                or early_stop
                            )

                    else:
                        if phase == "train":
                            log_energy_results(
                                log_epoch, epoch, now, loss, energy_rmse, phase
                            )
                        else:
                            log_energy_results(
                                log_epoch, epoch, now, "", energy_rmse, phase
                            )
                        if phase == "train":
                            # early stop when training energy error stagnates
                            if abs(energy_rmse - previous_energy_rmse) <= 1e-7:
                                print('Maximum Time Reached')
                                early_stop = True
                            previous_energy_rmse = energy_rmse
                        elif phase == "val":
                            if energy_rmse < best_val_energy_loss:
                                best_val_energy_loss = energy_rmse
                                best_model_wts = copy.deepcopy(self.model.state_dict())
                            convergence = (
                                (
                                    best_val_energy_loss
                                    <= self.convergence_criteria["energy"]
                                )
                                or early_stop
                                or (epoch >= self.epochs)
                            )

                print()

            else:
                phase = "train"

                if self.scheduler:
                    self.scheduler.step()
                self.model.train()

                energy_mse = 0.0
                force_mse = "N/A"
                if forcetraining:
                    force_mse = 0.0

                for data_sample in self.atoms_dataloader:
                    unique_atoms = data_sample[3]
                    input_data = [data_sample[0], len(data_sample[1]), unique_atoms]
                    target = data_sample[1].requires_grad_(False)
                    batch_size = len(target)
                    target = target.reshape(batch_size, 1).to(self.device)
                    scaled_target = (target - self.mean_scaling) / self.sd_scaling
                    num_of_atoms = data_sample[2].reshape(batch_size, 1).to(self.device)
                    fp_primes = data_sample[4]
                    if self.weighted:
                        weights = data_sample[-1]

                    if forcetraining:
                        fp_primes = fp_primes.to(self.device)
                        image_forces = data_sample[5].to(self.device)
                        scaled_forces = image_forces / self.sd_scaling
                    for element in unique_atoms:
                        input_data[0][element][0] = (
                            input_data[0][element][0]
                            .to(self.device)
                            .requires_grad_(True)
                        )
                    scaled_target = scaled_target.to(self.device)

                    def closure():
                        self.optimizer.zero_grad()
                        if forcetraining:
                            energy_pred, force_pred = self.model(input_data, fp_primes)
                            if self.weighted:
                                loss = self.criterion(
                                    energy_pred,
                                    scaled_target,
                                    num_of_atoms,
                                    force_pred,
                                    scaled_forces,
                                    self.model,
                                    weights=weights
                                )
                            else:
                                loss = self.criterion(
                                    energy_pred,
                                    scaled_target,
                                    num_of_atoms,
                                    force_pred,
                                    scaled_forces,
                                    self.model,
                                )
                        else:
                            energy_pred, _ = self.model(input_data)
                            if welf.weighted:
                                loss = self.criterion(
                                    energy_pred, scaled_target, num_of_atoms,
                                    model=self.model, weights=weights
                                )
                            else:
                                loss = self.criterion(
                                    energy_pred, scaled_target, num_of_atoms,
                                    model=self.model
                                )
                        loss.backward()
                        return loss

                    mse_loss = nn.MSELoss(reduction="sum")
                    energy_pred, force_pred = self.model(input_data, fp_primes)
                    raw_preds = (energy_pred * self.sd_scaling) + self.mean_scaling
                    raw_preds_per_atom = torch.div(raw_preds, num_of_atoms)
                    target_per_atom = torch.div(target, num_of_atoms)
                    energy_loss = mse_loss(raw_preds_per_atom, target_per_atom)
                    energy_mse += torch.tensor(energy_loss.item())

                    loss = self.optimizer.step(closure)
                    print("loss function value: %f" % loss.item())
                    now = time.asctime()

                    if forcetraining:
                        force_pred = force_pred * self.sd_scaling
                        num_atoms_force = torch.cat(
                            [idx.repeat(int(idx)) for idx in num_of_atoms]
                        )
                        num_atoms_force = torch.sqrt(num_atoms_force).reshape(
                            len(num_atoms_force), 1
                        )
                        force_pred_per_atom = torch.div(force_pred, num_atoms_force)
                        force_targets_per_atom = torch.div(
                            image_forces, num_atoms_force
                        )
                        force_loss = mse_loss(
                            force_pred_per_atom, force_targets_per_atom
                        )
                        # mean over image
                        force_loss /= 3
                        force_mse += torch.tensor(force_loss.item())

                    loss = self.optimizer.step(closure)
                    now = time.asctime()

                energy_mse /= self.dataset_size
                energy_rmse = torch.sqrt(energy_mse)
                if torch.isnan(energy_rmse):
                    early_stop = True
                plot_energy_loss[phase].append(energy_rmse)
                print("energy loss: %f" % energy_rmse)
                if forcetraining:
                    force_mse /= self.dataset_size
                    force_rmse = torch.sqrt(force_mse)
                    if torch.isnan(force_rmse):
                        early_stop = True
                    plot_force_loss[phase].append(force_rmse)
                    print("force loss: %f\n" % force_rmse)
                    log_force_results(
                        log_epoch, epoch, now, loss, energy_rmse, force_rmse, phase
                    )
                    # terminates when error stagnates
                    if abs(force_rmse - previous_force_rmse) <= 1e-5 and \
                            abs(energy_rmse - previous_energy_rmse)<= 1e-5:
                        print('Early Stop Triggered')
                        early_stop = True
                    time_elapsed = time.time() - since
                    if self.maxtime is not None:
                        if time_elapsed >= self.maxtime:
                             early_stop = True
                    if torch.isnan(force_rmse):
                        early_stop = True
                    if force_rmse < best_train_force_loss:
                        best_train_energy_loss = energy_rmse
                        best_train_force_loss = force_rmse
                        best_model_wts = copy.deepcopy(self.model.state_dict())
                    
                    previous_energy_rmse = energy_rmse
                    previous_force_rmse = force_rmse
                    previous_energy_rmse = energy_rmse
                    energy_convergence = (
                        best_train_energy_loss <= self.convergence_criteria["energy"]
                    )
                    force_convergence = (
                        best_train_force_loss <= self.convergence_criteria["force"]
                    )
                    convergence = (
                        (energy_convergence and force_convergence)
                        or early_stop
                        or (epoch >= self.epochs)
                    )
                else:
                    log_energy_results(log_epoch, epoch, now, loss, energy_rmse, phase)
                    # terminates when error stagnates
                    if abs(energy_rmse - previous_energy_rmse) <= 1e-7:
                        print('Early Stop Triggered')
                        convergence = True
                    time_elapsed = time.time() - since
                    if self.maxtime is not None:
                        if time_elapsed >= self.maxtime:
                             early_stop = True
                    elif energy_rmse < best_train_energy_loss:
                        best_train_energy_loss = energy_rmse
                        best_model_wts = copy.deepcopy(self.model.state_dict())
                    previous_energy_rmse = energy_rmse
                    convergence = (
                        (best_train_energy_loss <= self.convergence_criteria["energy"])
                        or early_stop
                        or (epoch >= self.epochs)
                    )

            epoch += 1

        log_epoch("")
        time_elapsed = time.time() - since
        print("Training complete in {} steps".format(epoch))
        print(
            "Training complete in {:.0f}m {:.0f}s".format(
                time_elapsed // 60, time_elapsed % 60
            )
        )

        log("Training complete in {} steps".format(epoch))
        if validation:
            log("Best validation energy loss: {:4f}".format(best_val_energy_loss))
            if forcetraining:
                log("Best validation force loss: {:4f}".format(best_val_force_loss))
        else:
            log("Best training energy loss: {:4f}".format(best_train_energy_loss))
            if forcetraining:
                log("Best training force loss: {:4f}".format(best_train_force_loss))
        log("")
        if not os.path.exists("results/plots/training"):
            os.makedirs("results/plots/training")
        plt.title("RMSE vs. Epoch")
        plt.xlabel("Epoch #")
        plt.ylabel("RMSE")
        plot_epoch_x = list(range(1, epoch + 1))
        plt.plot(plot_epoch_x, plot_energy_loss["train"], label="energy train")
        if validation:
            plt.plot(plot_epoch_x, plot_energy_loss["val"], label="energy val")
        if forcetraining:
            plt.plot(plot_epoch_x, plot_force_loss["train"], label="force train")
            if validation:
                plt.plot(plot_epoch_x, plot_force_loss["val"], label="force val")
        plt.legend()
        plt.savefig("results/plots/training/" + self.label + ".pdf")
        self.model.load_state_dict(best_model_wts)
        # print(force_pred)
        return self.model
