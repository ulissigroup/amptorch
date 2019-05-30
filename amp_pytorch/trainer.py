"""Trainer class used to train a specified model in accordance with the trainer
arguments"""

import sys
import time
import copy
from amp.utilities import Logger
import matplotlib.pyplot as plt
import torch.nn as nn
import torch

__author__ = "Muhammed Shuaibi"
__email__ = "mshuaibi@andrew.cmu.edu"


class Trainer:
    """Class utilized to train a given model through a training loop method.

    Parameters
    ----------

    model: object
        NeuralNetwork model to be trained
    device: str
        Hardware to be utilized for training - CPU or GPU
    unique_atoms: list
        List of unique atoms contained in the training dataset.
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
    rmse_criteria: dict
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
        unique_atoms,
        dataset_size,
        criterion,
        optimizer,
        scheduler,
        atoms_dataloader,
        rmse_criteria,
        scalings,
    ):

        self.model = model
        self.device = device
        self.unique_atoms = unique_atoms
        self.dataset_size = dataset_size
        self.criterion = criterion
        self.optimizer = optimizer
        self.scheduler = scheduler
        self.atoms_dataloader = atoms_dataloader
        self.rmse_criteria = rmse_criteria
        self.sd_scaling = scalings[0]
        self.mean_scaling = scalings[1]

    def train_model(self):
        "Training loop"
        forcetraining = False
        if self.criterion.alpha > 0:
            forcetraining = True
        best_force_loss = 1e8
        best_energy_loss = 1e8
        log = Logger("results/results-log.txt")
        log_epoch = Logger("results/epoch-log.txt")
        log("Model: %s" % self.model)

        since = time.time()
        log_epoch("-" * 50)
        print("Training Initiated!")
        log_epoch("%s Training Initiated!\n" % time.asctime())

        best_model_wts = copy.deepcopy(self.model.state_dict())

        plot_energy_loss = {"train": [], "val": []}
        if forcetraining:
            plot_force_loss = {"train": [], "val": []}

        epoch = 0
        convergence = False
        # while epoch <= 10:
        while not convergence:
            # epoch_timer = time.time()
            log_epoch("{} Epoch {}".format(time.asctime(), epoch + 1))
            log_epoch("-" * 30)

            if isinstance(self.atoms_dataloader, dict):

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
                        input_data = [data_sample[0], len(data_sample[1])]
                        fp_primes = data_sample[3]
                        image_forces = data_sample[4].to(self.device)
                        target = data_sample[1].requires_grad_(False)
                        batch_size = len(target)
                        target = target.reshape(batch_size, 1).to(self.device)
                        scaled_target = (target - self.mean_scaling) / self.sd_scaling
                        scaled_forces = image_forces / self.sd_scaling
                        num_of_atoms = (
                            data_sample[2].reshape(batch_size, 1).to(self.device)
                        )
                        for element in self.unique_atoms:
                            input_data[0][element][0] = (
                                input_data[0][element][0]
                                .to(self.device)
                                .requires_grad_(True)
                            )
                        scaled_target = scaled_target.to(self.device)

                        def closure():
                            self.optimizer.zero_grad()
                            energy_pred, force_pred = self.model(input_data, fp_primes)
                            loss = self.criterion(
                                energy_pred,
                                scaled_target,
                                force_pred,
                                scaled_forces,
                                num_of_atoms,
                            )
                            loss.backward()
                            return loss

                        energy_pred, force_pred = self.model(input_data, fp_primes)
                        if phase == "train":
                            self.optimizer.step(closure)
                        mse_loss = nn.MSELoss(reduction="sum")
                        raw_preds = (energy_pred * self.sd_scaling) + self.mean_scaling
                        raw_preds_per_atom = torch.div(raw_preds, num_of_atoms)
                        target_per_atom = torch.div(target, num_of_atoms)
                        energy_loss = mse_loss(raw_preds_per_atom, target_per_atom)
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
                            force_mse += torch.tensor(force_loss.item())

                    energy_mse /= self.dataset_size[phase]
                    energy_rmse = torch.sqrt(energy_mse)
                    plot_energy_loss[phase].append(torch.log10(energy_rmse))
                    print("%s energy loss: %f" % (phase, energy_rmse))
                    log_epoch("%s energy loss: %f" % (phase, energy_rmse))
                    if forcetraining:
                        force_mse /= self.dataset_size[phase]
                        force_rmse = torch.sqrt(force_mse)
                        plot_force_loss[phase].append(torch.log10(force_rmse))
                        print("%s force loss: %f" % (phase, force_rmse))
                        log_epoch("%s force loss: %f" % (phase, force_rmse))
                        if phase == "val" and (
                            (energy_rmse < best_energy_loss)
                            & (force_rmse < best_force_loss)
                        ):
                            best_energy_loss = energy_rmse
                            best_force_loss = force_rmse
                            best_model_wts = copy.deepcopy(self.model.state_dict())
                        energy_convergence = (
                            best_energy_loss <= self.rmse_criteria["energy"]
                        )
                        force_convergence = (
                            best_force_loss <= self.rmse_criteria["force"]
                        )
                    else:
                        if phase == "val" and energy_rmse < best_energy_loss:
                            best_energy_loss = energy_rmse
                            best_model_wts = copy.deepcopy(self.model.state_dict())
                    convergence = best_energy_loss <= self.rmse_criteria["energy"]
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

                # loading_timer = time.time()
                for data_sample in self.atoms_dataloader:
                    # print('data_loading: %s' % (time.time()-loading_timer))
                    input_data = [data_sample[0], len(data_sample[1])]
                    fp_primes = data_sample[3].to(self.device)
                    image_forces = data_sample[4].to(self.device)
                    target = data_sample[1].requires_grad_(False)
                    batch_size = len(target)
                    target = target.reshape(batch_size, 1).to(self.device)
                    scaled_target = (target - self.mean_scaling) / self.sd_scaling
                    scaled_forces = image_forces / self.sd_scaling
                    num_of_atoms = data_sample[2].reshape(batch_size, 1).to(self.device)
                    for element in self.unique_atoms:
                        input_data[0][element][0] = (
                            input_data[0][element][0]
                            .to(self.device)
                            .requires_grad_(True)
                        )
                    scaled_target = scaled_target.to(self.device)

                    def closure():
                        self.optimizer.zero_grad()
                        energy_pred, force_pred = self.model(input_data, fp_primes)
                        loss = self.criterion(
                            energy_pred,
                            scaled_target,
                            force_pred,
                            scaled_forces,
                            num_of_atoms,
                        )
                        loss.backward()
                        return loss

                    # optimizer_time = time.time()
                    self.optimizer.step(closure)
                    # print('optimizer.step(): %s' %(time.time()-optimizer_time))

                    mse_loss = nn.MSELoss(reduction="sum")
                    energy_pred, force_pred = self.model(input_data, fp_primes)
                    raw_preds = (energy_pred * self.sd_scaling) + self.mean_scaling
                    raw_preds_per_atom = torch.div(raw_preds, num_of_atoms)
                    target_per_atom = torch.div(target, num_of_atoms)
                    energy_loss = mse_loss(raw_preds_per_atom, target_per_atom)
                    energy_mse += torch.tensor(energy_loss.item())

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

                energy_mse /= self.dataset_size
                energy_rmse = torch.sqrt(energy_mse)
                plot_energy_loss[phase].append(torch.log10(energy_rmse))
                print("energy loss: %f" % energy_rmse)
                log_epoch("energy loss: %f" % (energy_rmse))
                if forcetraining:
                    force_mse /= self.dataset_size
                    force_rmse = torch.sqrt(force_mse)
                    plot_force_loss[phase].append(torch.log10(force_rmse))
                    print("force loss: %f\n" % force_rmse)
                    log_epoch("force loss: %f\n" % (force_rmse))
                    if (energy_rmse < best_energy_loss) and (
                        force_rmse < best_force_loss
                    ):
                        best_energy_loss = energy_rmse
                        best_force_loss = force_rmse
                        best_model_wts = copy.deepcopy(self.model.state_dict())
                    energy_convergence = (
                        best_energy_loss <= self.rmse_criteria["energy"]
                    )
                    force_convergence = best_force_loss <= self.rmse_criteria["force"]
                    convergence = energy_convergence and force_convergence
                else:
                    if energy_rmse < best_energy_loss:
                        best_energy_loss = energy_rmse
                        best_model_wts = copy.deepcopy(self.model.state_dict())
                    convergence = best_energy_loss <= self.rmse_criteria["energy"]

            epoch += 1

        time_elapsed = time.time() - since
        print("Training complete in {} steps".format(epoch))
        print(
            "Training complete in {:.0f}m {:.0f}s".format(
                time_elapsed // 60, time_elapsed % 60
            )
        )

        log("Training complete in {} steps".format(epoch))
        log("Best training energy loss: {:4f}\n".format(best_energy_loss))
        if forcetraining:
            log("Best training force loss: {:4f}\n".format(best_force_loss))

        plt.title("RMSE vs. Epoch")
        plt.xlabel("Epoch #")
        plt.ylabel("log(RMSE)")
        plot_epoch_x = list(range(1, epoch + 1))
        plt.plot(plot_epoch_x, plot_energy_loss[phase], label="energy train")
        if forcetraining:
            plt.plot(plot_epoch_x, plot_force_loss[phase], label="force train")
        plt.legend()
        plt.show()
        self.model.load_state_dict(best_model_wts)
        return self.model
