import sys
import numpy as np
import torch
import torch.nn
import copy
import time
import matplotlib.pyplot as plt
from amp.utilities import Logger


def target_scaling(data, method=None):
    if method == "minmax":
        data_max = max(data)
        data_min = min(data)
        scale = []
        for index, value in enumerate(data):
            scale.append((value - data_min) / (data_max - data_min))
        return torch.stack(scale)
    elif method == "standardize":
        data_mean = torch.mean(data)
        data_sd = torch.std(data, dim=0)
        scale = []
        for index, value in enumerate(data):
            scale.append((value - data_mean) / data_sd)
        return torch.stack(scale)
    else:
        return data


def pred_scaling(data, target, method=None):
    if method == "minmax":
        target_max = max(target)
        target_min = min(target)
        scale = []
        for index, value in enumerate(data):
            scale.append((value * (target_max - target_min)) + target_min)
        return torch.stack(scale)
    elif method == "standardize":
        target_mean = torch.mean(target)
        target_sd = torch.std(target, dim=0)
        scale = []
        for index, value in enumerate(data):
            scale.append((value * target_sd) + target_mean)
        return torch.stack(scale)
    else:
        return data


def train_model(
    model,
    device,
    unique_atoms,
    dataset_size,
    criterion,
    optimizer,
    atoms_dataloader,
    rmse_criteria,
):

    log = Logger("results/results-log.txt")
    log_epoch = Logger("results/epoch-log.txt")
    log("Model: %s" % model)

    model.train()

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
        log_epoch("{} Epoch {}".format(time.asctime(), epoch + 1))
        log_epoch("-" * 30)

        for phase in ["train", "val"]:

            if phase == "train":
                # scheduler.step()
                model.train()
            else:
                model.eval()

            MSE = 0.0

            for data_sample in atoms_dataloader[phase]:
                input_data = data_sample[0]
                target = data_sample[1].requires_grad_(False)
                batch_size = len(target)
                target = target.reshape(batch_size, 1)
                scaled_target = target_scaling(target, method="standardize")
                num_of_atoms = data_sample[2].reshape(batch_size, 1)
                for element in unique_atoms:
                    input_data[element][0] = input_data[element][0].to(device)
                scaled_target = scaled_target.to(device)

                optimizer.zero_grad()

                with torch.set_grad_enabled(phase == "train"):
                    output = model(input_data)
                    loss = criterion(output, scaled_target)
                    if phase == "train":
                        loss.backward()
                        optimizer.step()
                with torch.no_grad():
                    scaled_preds = model(input_data)
                    raw_preds = pred_scaling(scaled_preds, target, method="standardize")
                    raw_preds_per_atom = torch.div(raw_preds, num_of_atoms)
                    target_per_atom = torch.div(target, num_of_atoms)
                    loss = criterion(raw_preds_per_atom, target_per_atom)
                    MSE += loss.item() * batch_size

                # def closure():
                # optimizer.zero_grad()
                # output = model(input_data)
                # loss = criterion(output, scaled_target)
                # loss.backward()
                # return loss

                # optimizer.step(closure)

                # Perform predictions and compute loss
                # with torch.no_grad():
                # scaled_preds = model(input_data)
                # raw_preds = pred_scaling(
                # scaled_preds, target, method='standardize')
                # raw_preds_per_atom = torch.div(raw_preds, num_of_atoms)
                # target_per_atom = torch.div(target, num_of_atoms)
                # loss = criterion(raw_preds_per_atom, target_per_atom)
                # MSE += loss.item()*batch_size

            MSE = MSE / dataset_size[phase]
            RMSE = np.sqrt(MSE)
            epoch_loss = RMSE
            plot_loss_y[phase].append(np.log10(RMSE))

            log_epoch("{} {} Loss: {:.4f}".format(time.asctime(), phase, epoch_loss))

            if phase == "val" and epoch_loss < best_loss:
                best_loss = epoch_loss
                best_model_wts = copy.deepcopy(model.state_dict())
        log_epoch("{} Loss: {:.4f}".format(time.asctime(), epoch_loss))
        epoch += 1

    time_elapsed = time.time() - since
    print("Training complete in {} steps".format(epoch))
    print(
        "Training complete in {:.0f}m {:.0f}s".format(
            time_elapsed // 60, time_elapsed % 60
        )
    )

    # log("Training complete in {} steps".format(epoch))
    # log("Best training loss: {:4f}".format(best_loss))
    # log("")

    # plt.title("RMSE vs. Epoch")
    # plt.xlabel("Epoch #")
    # plt.ylabel("log(RMSE)")
    # plot_epoch_x = list(range(1, epoch + 1))
    # plt.plot(plot_epoch_x, plot_loss_y, label="train")
    # plt.legend()
    # plt.show()
    model.load_state_dict(best_model_wts)
    return model
