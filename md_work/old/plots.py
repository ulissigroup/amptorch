import sys
import ase
from ase import io
from ase.calculators.emt import EMT
from ase import units
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.stats as stats
from matplotlib import rc
import seaborn as sns
from ase.visualize import view


def kde_plots(emt, forces, labels, filename=None):
    fig, ax = plt.subplots(figsize=(14.15, 10))
    sns.distplot(
        emt,
        hist=False,
        kde=True,
        label="target",
        color="k",
        hist_kws={"edgecolor": "black", "alpha": 0.1},
        kde_kws={"linewidth": 2},
        ax=ax,
    )
    for i, data in enumerate(forces):
        sns.distplot(
            data,
            hist=False,
            label=labels[i],
            kde=True,
            hist_kws={"edgecolor": "black", "alpha": 0.1},
            kde_kws={"linewidth": 2},
        )
    plt.legend(loc="best", fontsize=28)
    plt.xlabel("log max|F|", fontsize=32)
    plt.ylabel("Density Function", fontsize=32)
    plt.title("Force Distributions", fontsize=35)
    plt.xlim(left=-5, right=4)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.show()


def time_rmse_plots(
    target, preds, sampled_points, labels, property="energy", filename=None
):
    fig, ax = plt.subplots(figsize=(14.15, 10))
    time = np.linspace(0, 20 * len(target), 101)
    if property == "force":
        for i, pred in enumerate(preds):
            force_rmse = []
            for k, image in enumerate(target):
                pred_image = pred[k]
                num_atoms = image.shape[1]
                mse = (1 / (3 * num_atoms)) * ((image - pred_image) ** 2).sum()
                force_rmse.append(np.sqrt(mse))
            plt.plot(time, force_rmse, label=labels[i])
            if sampled_points and i == 0:
                sampled_time = [time[t] for t in sampled_points]
                samples = [force_rmse[p] for p in sampled_points]
                plt.plot(sampled_time, samples, "o", label="sampled points")
    elif property == "energy":
        for i, pred in enumerate(preds):
            energy_rmse = []
            for k, image in enumerate(target):
                pred_image = pred[k]
                mse = (1 / (29**2)) * ((image - pred_image) ** 2)
                energy_rmse.append(np.sqrt(mse))
            plt.plot(time, energy_rmse, label=labels[i])
            if sampled_points and i == 0:
                sampled_time = [time[t] for t in sampled_points]
                samples = [energy_rmse[p] for p in sampled_points]
                plt.plot(sampled_time, samples, "o", label="sampled points")
    plt.xlabel("time, (ps)", fontsize=25)
    if property == "energy":
        plt.ylabel("energy error, eV/atom", fontsize=25)
    else:
        plt.ylabel("force error, eV/A", fontsize=25)
    plt.legend(fontsize=25)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylim(ymin=0, ymax=10)
    plt.show()

def time_plots(target, preds, sampled_points, label, property="energy", filename=None):
    fig, ax = plt.subplots(figsize=(14.15, 10))
    time = np.linspace(0, 20 * len(target), 101)
    if sampled_points:
        sampled_time = [time[i] for i in sampled_points]
        samples = [preds[0][i] for i in sampled_points]
        plt.plot(sampled_time, samples, "o", label="sampled points")
    plt.plot(time, target, color="k", label="target")
    plt.plot(time, preds[0], label=label)
    plt.plot(time, preds[1], label=label + "- resample")
    plt.xlabel("time, (ps)", fontsize=25)
    if property == "energy":
        plt.ylabel("energy, eV", fontsize=25)
    else:
        plt.ylabel("logmax|F|", fontsize=25)
    plt.legend(fontsize=25)
    plt.xticks(fontsize=20)
    # plt.ylim(ymin=8.5, ymax=50)
    plt.yticks(fontsize=20)
    plt.show()


def compute_forces(dir, files, type="max"):
    forces = []
    for file in files:
        images = ase.io.read(dir + file, ":")
        data = []
        for image in images:
            image.set_calculator(EMT())
            if type == "max":
                data.append(np.log10(np.amax(np.abs(image.get_forces()))))
            else:
                data.append(image.get_forces())
        if type == "max":
            forces.append((np.array(data)))
        else:
            forces.append(data)
    return forces


def compute_energies(dir, files):
    energies = []
    for file in files:
        images = ase.io.read(dir + file, ":")
        data = []
        for image in images:
            image.set_calculator(EMT())
            data.append(image.get_potential_energy())
        energies.append(data)
    return energies


images_emt = ase.io.read("../../datasets/COCu/COCu_pbc_300K.traj", ":")
images0 = []
forces0 = []
for i in range(101):
    images0.append(images_emt[i])
    forces0.append(np.log10((np.array(np.amax(np.abs(images_emt[i].get_forces()))))))
    # forces0.append(images_emt[i].get_forces())
energies0 = [image.get_potential_energy() for image in images0]

files = [
    "MLMD_COCu_pbc_300K_tanh_LJ_10_iter_1.traj",
    "MLMD_COCu_pbc_300K_tanh_LJ_10_iter_2.traj",
    "MLMD_COCu_pbc_300K_tanh_LJ_10_iter_3.traj",
    "MLMD_COCu_pbc_300K_tanh_LJ_10_iter_4.traj",
    "MLMD_COCu_pbc_300K_tanh_LJ_10_iter_5.traj",
]

forces = compute_forces(dir="MD_results/COCu/pbc_300K/tanh/paper/", files=files, type="max")
# sampled_points = [10, 12, 14, 16, 18, 20, 22, 24, 26, 30]
# sampled_points = [65, 6, 7, 60, 63, 20, 56 ,4, 44, 22]
# sampled_points = None
# energies = compute_energies(dir="MD_results/COCu/pbc_300K/tanh/", files=files)
# time_rmse_plots(
    # forces0, forces, sampled_points, labels=["ML-LJ", "ML-LJ resample"], property="force", filename="test"
# )
# time_rmse_plots(
    # energies0,
    # energies,
    # sampled_points,
    # labels=["ML-LJ", "ML-LJ resample"],
    # property="energy",
    # filename="test",
# )
# rmse_plots(forces0, forces[0], 10, "test")
kde_plots(forces0, forces, ['1',"2", '3', "4", "5"], "test")
