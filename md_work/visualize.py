import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import ase
from ase import units
from ase.md import Langevin
import copy
from ase.calculators.emt import EMT
from ase.md.nvtberendsen import NVTBerendsen
from ase.md.verlet import VelocityVerlet
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
import numpy as np


def calculate_energies(emt_images, ml_images):
    energies_emt = [image.get_potential_energy() for image in emt_images]
    ml_energies_apparent = [image.get_potential_energy() for image in ml_images]
    ml_energies_actual = []
    for image in ml_images:
        image = copy.copy(image)
        image.set_calculator(EMT())
        ml_energies_actual.append(image.get_potential_energy())
    return energies_emt, ml_energies_apparent, ml_energies_actual


def calculate_forces(emt_images, ml_images, type="max"):
    if type == "max":
        forces_emt = [
            np.log10((np.array(np.amax(np.abs(image.get_forces())))))
            for image in emt_images
        ]
        ml_forces_apparent = [
            np.log10((np.array(np.amax(np.abs(image.get_forces())))))
            for image in ml_images
        ]
        ml_forces_actual = []
        for image in ml_images:
            image = copy.copy(image)
            image.set_calculator(EMT())
            ml_forces_actual.append(
                np.log10((np.array(np.amax(np.abs(image.get_forces())))))
            )
    else:
        forces_emt = [image.get_forces() for image in emt_images]
        ml_forces_apparent = [image.get_forces() for image in ml_images]
        ml_forces_actual = []
        for image in ml_images:
            image = copy.copy(image)
            image.set_calculator(EMT())
            ml_forces_actual.append(image.get_forces())
    return forces_emt, ml_forces_apparent, ml_forces_actual


def time_plots(target, preds, label, property="energy", sampled_points=None):
    fig, ax = plt.subplots(figsize=(14.15, 10))
    time = np.linspace(0, len(target), len(target))
    if sampled_points:
        sampled_time = [time[i] for i in sampled_points]
        samples = [preds[1][i] for i in sampled_points]
        plt.plot(sampled_time, samples, "ro", label="sampled points")
    plt.plot(time, target, color="k", label="target")
    colors = ["b", "orange", "g"]
    for idx, data in enumerate(preds):
        if data:
            plt.plot(time, data, color=colors[idx], label=label[idx])
    plt.xlabel("time, (fs)", fontsize=25)
    if property == "energy":
        plt.ylabel("energy, eV", fontsize=25)
    else:
        plt.ylabel("logmax|F|", fontsize=25)
    plt.legend(fontsize=25)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.show()


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
    colors = ["b", "orange", "g"]
    for i, data in enumerate(forces):
        sns.distplot(
            data,
            hist=False,
            label=labels[i],
            kde=True,
            hist_kws={"edgecolor": "black", "alpha": 0.1},
            color=colors[i],
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


def parity_plots(actual_energies, apparent_energies):
    e_min = min(apparent_energies)
    e_max = max(apparent_energies)
    plt.plot([e_min, e_max], [e_min, e_max])
    plt.plot(actual_energies, apparent_energies, "bo")
    plt.show()


emt_images = ase.io.read("../datasets/COCu_lang_1fs_300K.traj", ":")
ml_images = ase.io.read("./lang_results/COCu_lang_2ps_EF_300K-2.traj", ":")
mllj_images = ase.io.read("./lang_results/COCu_lang_2ps_EF_300K-LJ-2.traj", ":")
# mllj_res_images = ase.io.read(
    # "./lang_results/f_lj/COCu_lang_1fs_ampG_res_300K-LJ-2.traj", ":2000"
# )

emt_energies, ml_apparent, ml_actual = calculate_energies(emt_images, ml_images)
_, mllj_apparent, mllj_actual = calculate_energies(emt_images, mllj_images)
# _, mllj_res_apparent, mllj_res_actual = calculate_energies(emt_images, mllj_res_images)

emt_forces, ml_apparent_forces, ml_actual_forces = calculate_forces(
    emt_images, ml_images
)
_, mllj_apparent_forces, mllj_actual_forces = calculate_forces(emt_images, mllj_images)
# _, mllj_res_apparent_forces, mllj_res_actual_forces = calculate_forces(
    # emt_images, mllj_res_images
# )
# sample_points = [488, 1214, 1115, 268, 758, 1876, 1237, 971, 1282, 1190]
time_plots(
    emt_energies,
    [ml_actual, mllj_actual],
    ["ML", "ML-LJ"],
    # [ml_actual, mllj_actual, mllj_res_actual],
    # ["ML", "ML-LJ", "ML-LJ resample"],
    "energy",
    # sampled_points=sample_points,
)

time_plots(
    emt_energies,
    [None, mllj_actual],
    [None, "ML-LJ"],
    # [None, mllj_actual, mllj_res_actual],
    # [None, "ML-LJ", "ML-LJ resample"],
    "energy",
    # sampled_points=sample_points,
)
time_plots(
    emt_forces,
    # [ml_actual_forces, mllj_actual_forces, mllj_res_actual_forces],
    # ["ML", "ML-LJ", "ML-LJ resample"],
    [ml_actual_forces, mllj_actual_forces],
    ["ML", "ML-LJ"],
    "forces",
    # sampled_points=sample_points,
)
kde_plots(
    emt_forces,
    [ml_actual_forces, mllj_actual_forces],
    ["ML", "ML-LJ"],
    # [ml_actual_forces, mllj_actual_forces, mllj_res_actual_forces],
    # ["ML", "ML-LJ", "ML-LJ resample"],
)
