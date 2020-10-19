import copy
import os

import matplotlib
import numpy as np
from matplotlib import pylab as plt

matplotlib.use("Agg")


def loss_curve(
    trainer, label, loss_critera=["train_loss", "valid_loss"], extra_label=None
):
    """
    Plot the loss curve for different scoring methods listed in labels.
    Input:
        net: a class instance of skorch.regressor.NeuralNetRegressor
        loss_critera: list of strings in net.history to plot against number of epochs.
    """
    # generate directory to save files if not already
    if not os.path.exists("results/plots/loss_curve"):
        os.makedirs("results/plots/loss_curve")
    hist = net.history
    fig, ax = plt.subplots(figsize=(6.0, 6.0))
    for _ in loss_critera:
        xs = hist[:, "epoch"]
        ys = hist[:, _]
        ax.plot(xs, ys, label=_)
        print((np.amax(ys) / np.amin(ys)))
        if (np.amax(ys) / np.amin(ys)) > 2e2:
            ax.set_yscale("log")

    ax.set_xlabel("epoch")
    ax.set_ylabel("loss")
    ax.set_title("Loss vs. Epochs")
    ax.legend()
    if extra_label is None:
        fig.savefig("results/plots/loss_curve/loss_" + label + ".pdf")
    else:
        fig.savefig(
            "results/plots/loss_curve/loss_" + label + "_" + extra_label + ".pdf"
        )


def train_test_analysis(
    calc, images_train, label, images_test=None, data="energy", legend=None
):
    """
    Train the neural network model with testing sets separated from training dataset.
    Input:
        calc: AMP class
        images_train: List of ase.Atoms objects with potential energy attributes used to train calc.
        images_test: List of ase.Atoms objects with potential energy attributes not evaluated in training of calc.
        data: string. "energy" means to plot energy parity plot
    Output:
        None. Save to results/ directory.
    """
    # deeepcopy train and test images to avoid changing potential energy attributes
    # from single piont calculator
    images_train = copy.copy(images_train)
    images_test = copy.copy(images_test)

    # generate directory to save files if not already
    if not os.path.exists("results/plots/test_split"):
        os.makedirs("results/plots/test_split")
    # sort out customized legend definitions
    if legend == None:
        labels = ["training", "testing"]
    else:
        if type(legend) is list:
            labels = [legend[0], None]
        else:
            labels = legend

    rmse_dict = {}
    if data == "energy":
        # parity plot with training and testing sets
        fig, ax = plt.subplots(figsize=(6.0, 6.0))
        for images_collect in [images_train, images_test]:
            if images_collect is not None:

                targets = []
                preds = []
                for image in images_collect:
                    num_atoms = len(image)  # calc num of atoms in training data
                    targets.append(image.get_potential_energy() / num_atoms)
                    image.set_calculator(calc)
                    preds.append(image.get_potential_energy() / num_atoms)
                targets = np.asarray(targets)
                preds = np.asarray(preds)
                if images_collect == images_train:
                    ax.scatter(
                        targets, preds, marker="v", s=7, alpha=0.6, label=labels[0]
                    )
                    # print set RMSE
                    rmse = np.sqrt(((targets - preds) ** 2).sum() / len(preds))
                    rmse_dict["training RMSE"] = rmse
                    print(
                        "Energy RMSE of training set = {:.4f} eV per atom.".format(rmse)
                    )
                elif images_collect == images_test:
                    ax.scatter(
                        targets, preds, marker="o", s=7, alpha=0.8, label=labels[1]
                    )
                    # print set RMSE
                    rmse = np.sqrt(((targets - preds) ** 2).sum() / len(preds))
                    rmse_dict["testing RMSE"] = rmse
                    print(
                        "Energy RMSE of testing set = {:.4f} eV per atom.".format(rmse)
                    )

        ymin, ymax = ax.get_ylim()
        xmin, xmax = ax.get_xlim()
        _min = np.amin([xmin, ymin])
        _max = np.amax([xmax, ymax])
        ax.plot([_min, _max], [_min, _max], "r-", alpha=0.6, lw=1)
        ax.set_xlabel("simulation energy, eV/atom")
        ax.set_ylabel("AMPTorch energy, eV/atom")
        ax.set_title("Energies")
        ax.legend()
        fig.savefig("results/plots/test_split/parity_E_" + label + ".pdf")

        print("Finished.")

    return rmse_dict


def subsample_analysis(
    calc,
    images_train,
    label,
    images_whole=None,
    images_test=None,
    data="energy",
    legend=None,
):
    """
    Train the neural network model with testing sets separated from training dataset.
    Input:
        calc: AMP class
        images_train: List of ase.Atoms objects with potential energy attributes used to train calc.
        images_test: List of ase.Atoms objects with potential energy attributes not evaluated in training of calc.
        data: string. "energy" means to plot energy parity plot
    Output:
        None. Save to results/ directory.
    """
    # deeepcopy train and test images to avoid changing potential energy attributes
    # from single piont calculator
    images_train = copy.deepcopy(images_train)
    images_test = copy.deepcopy(images_test)
    images_whole = copy.deepcopy(images_whole)

    # generate directory to save files if not already
    if not os.path.exists("results/plots/test_split"):
        os.makedirs("results/plots/test_split")
    # sort out customized legend definitions
    if legend == None:
        labels = ["whole", "train", "test"]
    else:
        labels = legend

    rmse_dict = {}
    if data == "energy":
        # parity plot with training and testing sets
        fig, ax = plt.subplots(figsize=(6.0, 6.0))
        for images_collect in [images_whole, images_train, images_test]:
            if images_collect is not None:
                targets = []
                preds = []
                for image in images_collect:
                    num_atoms = len(image)  # calc num of atoms in training data
                    targets.append(image.get_potential_energy() / num_atoms)
                    image.set_calculator(calc)
                    preds.append(image.get_potential_energy() / num_atoms)
                targets = np.asarray(targets)
                preds = np.asarray(preds)
                if images_collect == images_whole:
                    ax.scatter(
                        targets, preds, marker="v", s=7, alpha=0.4, label=labels[0]
                    )
                    # print set RMSE
                    rmse = np.sqrt(((targets - preds) ** 2).sum() / len(preds))
                    rmse_dict["rmse_whole"] = rmse
                    max_abs_err = np.amax(np.absolute(targets - preds))
                    rmse_dict["maxErr_whole"] = max_abs_err
                    print(
                        "Energy RMSE of the whole dataset = {:.4f} eV per atom.".format(
                            rmse
                        )
                    )
                elif images_collect == images_train:
                    ax.scatter(
                        targets, preds, marker="o", s=7, alpha=0.4, label=labels[1]
                    )
                    # print set RMSE
                    rmse = np.sqrt(((targets - preds) ** 2).sum() / len(preds))
                    rmse_dict["rmse_train"] = rmse
                    max_abs_err = np.amax(np.absolute(targets - preds))
                    rmse_dict["maxErr_train"] = max_abs_err
                    print(
                        "Energy RMSE of subsampled training set = {:.4f} eV per atom.".format(
                            rmse
                        )
                    )
                elif images_collect == images_test:
                    ax.scatter(
                        targets, preds, marker="o", s=7, alpha=0.4, label=labels[2]
                    )
                    # print set RMSE
                    rmse = np.sqrt(((targets - preds) ** 2).sum() / len(preds))
                    rmse_dict["rmse_test"] = rmse
                    max_abs_err = np.amax(np.absolute(targets - preds))
                    rmse_dict["maxErr_test"] = max_abs_err
                    print("Energy RMSE of test set = {:.4f} eV per atom.".format(rmse))
        ymin, ymax = ax.get_ylim()
        xmin, xmax = ax.get_xlim()
        _min = np.amin([xmin, ymin])
        _max = np.amax([xmax, ymax])
        ax.plot([_min, _max], [_min, _max], "r-", alpha=0.6, lw=1)
        ax.set_xlabel("simulation energy, eV/atom")
        ax.set_ylabel("AMPTorch energy, eV/atom")
        ax.set_title("Energies")
        ax.legend()
        fig.savefig("results/plots/test_split/parity_E_" + label + "_subsample.pdf")

        print("Finished.")

    return rmse_dict
