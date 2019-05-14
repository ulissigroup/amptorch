"""GD_trainer.py: Trains a provided model utilizing GD based algorithms
including SGD, Adam, etc. Model convergence is achieved upon reaching the
specified rmse convergence criteria."""

import sys
import time
import copy
from amp.utilities import Logger
import matplotlib.pyplot as plt
import torch.nn as nn
import torch
import numpy as np
from amp_pytorch.NN_model import ForceLossFunction

__author__ = "Muhammed Shuaibi"
__email__ = "mshuaibi@andrew.cmu.edu"


def target_scaling(data, method=None):
    """Scales the target values through a specified method:

        method: 'minmax' = scales values to between[0,1]
                'standardize' = scales values to have a mean of 0 and standard
                deviation of 1
        """
    data_ = data.reshape(-1, 1)
    if method == "minmax":
        data_max = max(data_)
        data_min = min(data_)
        data = (data - data_min) / (data_max - data_min)
    elif method == "standardize":
        data_mean = torch.mean(data_)
        data_sd = torch.std(data_, dim=0)
        data = (data - data_mean) / data_sd
    return data


def pred_scaling(data, target, method=None):
    """Unscales model predictions based off the previous scaling method.
    Unscaling is necessary in order to compute RMSE values off the raw targets
    """

    target = target.reshape(-1, 1)
    if method == "minmax":
        target_max = max(target)
        target_min = min(target)
        data = ((data * (target_max - target_min)) + target_min)
    elif method == "standardize":
        target_mean = torch.mean(target)
        target_sd = torch.std(target, dim=0)
        data = ((data * target_sd) + target_mean)
    return data


def train_model(
    model,
    device,
    unique_atoms,
    dataset_size,
    criterion,
    optimizer,
    scheduler,
    atoms_dataloader,
    rmse_criteria,
):

    log = Logger("results/results-log.txt")
    log_epoch = Logger("results/epoch-log.txt")
    log("Model: %s" % model)

    since = time.time()
    log_epoch("-" * 50)
    print("Training Initiated!")
    log_epoch("%s Training Initiated!" % time.asctime())
    log_epoch("")

    best_model_wts = copy.deepcopy(model.state_dict())
    best_loss = 1e8

    plot_loss_y = {"train": [], "val": []}

    epoch = 0
    while best_loss >= rmse_criteria:
    # while epoch <= 30:
        log_epoch("{} Epoch {}".format(time.asctime(), epoch + 1))
        log_epoch("-" * 30)

        if isinstance(atoms_dataloader, dict):

            for phase in ["train", "val"]:

                if phase == "train":
                    if scheduler:
                        scheduler.step()
                    model.train()
                else:
                    model.eval()

                mse = 0.0

                for data_sample in atoms_dataloader[phase]:
                    input_data = [data_sample[0], len(data_sample[1])]
                    fp_primes = data_sample[3]
                    image_forces = data_sample[4].to(device)
                    target = data_sample[1].requires_grad_(False)
                    batch_size = len(target)
                    target = target.reshape(batch_size, 1).to(device)
                    scaled_target = target_scaling(target, method="standardize")
                    num_of_atoms = data_sample[2].reshape(batch_size, 1).to(device)
                    for element in unique_atoms:
                        input_data[0][element][0] = (
                            input_data[0][element][0].to(device).requires_grad_(True)
                        )
                    scaled_target = scaled_target.to(device)

                    def closure():
                        optimizer.zero_grad()
                        energy_pred, force_pred = model(input_data, fp_primes)
                        CustomLoss = ForceLossFunction()
                        loss = CustomLoss(
                            energy_pred,
                            scaled_target,
                            force_pred,
                            image_forces,
                            num_of_atoms,
                        )
                        loss.backward()
                        return loss

                    energy_pred, force_pred = model(input_data, fp_primes)
                    raw_preds = pred_scaling(energy_pred, target, method="standardize")
                    raw_preds_per_atom = torch.div(raw_preds, num_of_atoms)
                    target_per_atom = torch.div(target, num_of_atoms)
                    MSELoss = nn.MSELoss()
                    actual_loss = MSELoss(raw_preds_per_atom, target_per_atom)

                    if phase == "train":
                        optimizer.step(closure)

                    mse += actual_loss.item()

                mse = mse / dataset_size[phase]
                rmse = np.sqrt(mse)
                epoch_loss = rmse
                plot_loss_y[phase].append(np.log10(rmse))

                print("{} Loss: {:.4f}".format(phase, epoch_loss))

                log_epoch(
                    "{} {} Loss: {:.4f}".format(time.asctime(), phase, epoch_loss)
                )

                if phase == "val" and epoch_loss < best_loss:
                    best_loss = epoch_loss
                    best_model_wts = copy.deepcopy(model.state_dict())

                log_epoch("")

        else:
            phase = "train"

            if scheduler:
                scheduler.step()
            model.train()

            mse = 0.0
            mse_f = 0.0

            for data_sample in atoms_dataloader:
                input_data = [data_sample[0], len(data_sample[1])]
                fp_primes = data_sample[3]
                image_forces = data_sample[4].to(device)
                target = data_sample[1].requires_grad_(False)
                batch_size = len(target)
                target = target.reshape(batch_size, 1).to(device)
                scaled_target = target_scaling(target, method="standardize")
                scaled_forces = target_scaling(image_forces,
                                               method="standardize")
                num_of_atoms = data_sample[2].reshape(batch_size, 1).to(device)
                for element in unique_atoms:
                    input_data[0][element][0] = (
                        input_data[0][element][0].to(device).requires_grad_(True)
                    )
                scaled_target = scaled_target.to(device)

                def closure():
                    optimizer.zero_grad()
                    energy_pred, force_pred = model(input_data, fp_primes)
                    # loss = criterion(energy_pred, scaled_target)
                    # loss = criterion(force_pred, image_forces)
                    CustomLoss = ForceLossFunction()
                    loss = CustomLoss(
                        energy_pred,
                        scaled_target,
                        force_pred,
                        # scaled_forces,
                        image_forces,
                        num_of_atoms,
                    )
                    loss.backward()
                    return loss

                energy_pred, force_pred = model(input_data, fp_primes)
                raw_preds = pred_scaling(energy_pred, target, method="standardize")
                raw_preds_per_atom = torch.div(raw_preds, num_of_atoms)
                target_per_atom = torch.div(target, num_of_atoms)
                actual_loss = criterion(raw_preds_per_atom, target_per_atom)

                # force_pred = pred_scaling(force_pred, image_forces, method="standardize")
                num_atoms_force = torch.cat(
                    [idx.repeat(int(idx)) for idx in num_of_atoms]
                )
                num_atoms_force = torch.sqrt(
                    num_atoms_force.reshape(len(num_atoms_force), 1)
                )
                force_pred_per_atom = torch.div(force_pred, num_atoms_force)
                force_targets_per_atom = torch.div(image_forces, num_atoms_force)
                force_loss = criterion(force_pred_per_atom, force_targets_per_atom)

                optimizer.step(closure)

                mse += actual_loss.item()
                mse_f += force_loss.item()

            mse = mse / dataset_size
            mse_f = mse_f / dataset_size
            rmse = np.sqrt(mse)
            rmse_f = np.sqrt(mse_f)
            epoch_loss = rmse_f
            print("energy: %s" % rmse)
            print("force: %s" % rmse_f)
            print("")
            # print(epoch_loss)
            plot_loss_y[phase].append(np.log10(rmse))

            if epoch_loss < best_loss:
                best_loss = epoch_loss
                best_model_wts = copy.deepcopy(model.state_dict())

            log_epoch("{} {} Loss: {:.4f}".format(time.asctime(), phase, epoch_loss))
            log_epoch("")

        epoch += 1
    # print(target)
    # print(raw_preds)
    # print(image_forces)
    # print(force_pred)

    time_elapsed = time.time() - since
    print("Training complete in {} steps".format(epoch))
    print(
        "Training complete in {:.0f}m {:.0f}s".format(
            time_elapsed // 60, time_elapsed % 60
        )
    )

    log("Training complete in {} steps".format(epoch))
    log("Best training loss: {:4f}".format(best_loss))
    log("")

    plt.title("RMSE vs. Epoch")
    plt.xlabel("Epoch #")
    plt.ylabel("log(RMSE)")
    plot_epoch_x = list(range(1, epoch + 1))
    plt.plot(plot_epoch_x, plot_loss_y[phase], label="train")
    plt.legend()
    plt.show()
    model.load_state_dict(best_model_wts)
    return model
