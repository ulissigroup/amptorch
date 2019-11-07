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


def rmse_plots(emt, forces, interval, filename=None):
    fig, ax = plt.subplots(figsize=(14.15, 10))
    rmse = []
    resample_images = [interval * (i + 1) for i in range(len(forces))]
    for file_forces in forces:
        num_atoms = file_forces.shape[1]
        mse = (1 / (3 * num_atoms)) * ((emt - file_forces) ** 2).sum()
        rmse.append(np.sqrt(mse))
    plt.plot(resample_images, rmse)
    plt.xlabel("Images resampled", fontsize=25)
    plt.ylabel("Force RMSE (eV/A)", fontsize=25)
    plt.title("RMSE vs Resampled Images", fontsize=25)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.show()


def time_plots(target, preds, sampled_points, label, property="energy", filename=None):
    fig, ax = plt.subplots(figsize=(14.15, 10))
    time = np.linspace(0, 20 * len(target), 101)
    sampled_time = [time[i] for i in sampled_points]
    samples = [preds[0][i] for i in sampled_points]
    plt.plot(time, target, color="k", label='target')
    plt.plot(time, preds[0], label=label)
    plt.plot(time, preds[1], label=label+" resample")
    plt.plot(sampled_time, samples, 'o', label='sampled points')
    plt.xlabel("time, (ps)", fontsize=25)
    if property == "energy":
        plt.ylabel("energy, eV", fontsize=25)
    else:
        plt.ylabel("logmax|F|", fontsize=25)
    plt.legend(fontsize=25)
    plt.xticks(fontsize=20)
    # plt.ylim(ymin=8.5, ymax=9)
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
    "MLMD_COCu_pbc_300K_tanh-LJ-1.traj",
    # "MLMD_COCu_pbc_300K_l2amp_fixed_10_resample_2.traj",
    # "MLMD_COCu_pbc_300K_l2amp_fixed_LJ_10_resample_2.traj",
    "MLMD_COCu_pbc_300K_tanh2_LJ_10_resample_2.traj",
    # "MLMD_COCu_pbc_300K_tanh_fixed_10_resample_2.traj",
    # "MLMD_COCu_pbc_300K_tanh_fixed_LJ_10_resample_2.traj",
]

forces = compute_forces(
    dir="MD_results/COCu/pbc_300K/tanh/paper/", files=files, type="max"
)
# sampled_points = [10, 12, 14, 16, 18, 20, 22, 24, 26, 30]
# sampled_points = [73, 63, 3, 12, 23, 68, 16, 97, 89, 19]
sampled_points = [79, 77, 35, 17, 70, 41, 9, 66, 28, 75]
energies = compute_energies(dir="MD_results/COCu/pbc_300K/tanh/paper/", files=files)
time_plots(forces0, forces, sampled_points, label='ML-LJ', property="forces", filename="test")
time_plots(energies0, energies, sampled_points, label='ML-LJ', property="energy", filename="test")
# rmse_plots(forces0, forces[0], 10, "test")
# kde_plots(forces0, forces, ['ML', 'ML-LJ'], "test")
