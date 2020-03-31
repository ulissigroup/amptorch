import numpy as np
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase.md.nvtberendsen import NVTBerendsen
from ase.calculators.emt import EMT
import copy
from ase.md import Langevin
from ase import units
import ase
import os
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use("Agg")


def md_run(
    calc, starting_image, temp, dt, count, label, ensemble="NVE",
    printenergy=False, rng=True
):
    """dt (fs)"""
    slab = starting_image.copy()
    slab.set_calculator(calc)
    slab.get_forces(apply_constraint=False)
    if rng:
        np.random.seed(1)
        MaxwellBoltzmannDistribution(slab, temp * units.kB)
    if ensemble == "NVE":
        dyn = VelocityVerlet(slab, dt * units.fs)
    elif ensemble == "nvtberendsen":
        dyn = NVTBerendsen(slab, dt * units.fs, temp, taut=300 * units.fs)
    elif ensemble == "langevin":
        dyn = Langevin(slab, dt * units.fs, temp * units.kB, 0.002)
    traj = ase.io.Trajectory(label + ".traj", "w", slab)
    dyn.attach(traj.write, interval=1)

    def printenergy(a=slab):
        """Function to print( the potential, kinetic, and total energy)"""
        epot = a.get_potential_energy() / len(a)
        ekin = a.get_kinetic_energy() / len(a)
        print(
            "Energy per atom: Epot = %.3feV Ekin = %.3feV (T=%3.0fK) "
            "Etot = %.3feV" % (epot, ekin, ekin / (1.5 * units.kB), epot + ekin)
        )

    if printenergy:
        dyn.attach(printenergy, interval=10)
    dyn.run(count - 1)


def calc_target(emt_images, ftype="max"):
    energies_emt = [image.get_potential_energy(apply_constraint=False) for image in emt_images]
    if ftype == "max":
        forces_emt = [
            np.log10((np.array(np.amax(np.abs(image.get_forces(apply_constraint=False))))))
            for image in emt_images
        ]
    else:
        forces_emt = [image.get_forces(apply_constraint=False) for image in emt_images]
    return energies_emt, forces_emt


def calculate_energies(ml_images):
    ml_energies_apparent = [image.get_potential_energy(apply_constraint=False) for image in ml_images]
    ml_energies_actual = []
    for image in ml_images:
        image = copy.copy(image)
        image.set_calculator(EMT())
        ml_energies_actual.append(image.get_potential_energy(apply_constraint=False))
    return ml_energies_apparent, ml_energies_actual


def calculate_forces(ml_images, type="max"):
    if type == "max":
        ml_forces_apparent = [
            np.log10((np.array(np.amax(np.abs(image.get_forces(apply_constraint=False))))))
            for image in ml_images
        ]
        ml_forces_actual = []
        for image in ml_images:
            image = copy.copy(image)
            image.set_calculator(EMT())
            ml_forces_actual.append(
                np.log10((np.array(np.amax(np.abs(image.get_forces(apply_constraint=False))))))
            )
    else:
        ml_forces_apparent = [image.get_forces(apply_constraint=False) for image in ml_images]
        ml_forces_actual = []
        for image in ml_images:
            image = copy.copy(image)
            image.set_calculator(EMT())
            ml_forces_actual.append(image.get_forces(apply_constraint=False))
    return ml_forces_apparent, ml_forces_actual


def calculate_rmse(target, preds, num_atoms, dtype="energy"):
    if dtype == "energy":
        target = np.array(target)
        preds = np.array(preds)
        energy_rmse = np.sqrt(((preds - target) / num_atoms) ** 2).sum() / len(preds)
        return energy_rmse
    else:
        force_mse = 0
        target_size = len(target)
        for i, j in zip(target, preds):
            i = np.array(i)
            j = np.array(j)
            force_mse += ((i - j) ** 2).sum() / (3 * num_atoms)
        force_rmse = np.sqrt(force_mse / target_size)
        return force_rmse


def time_plots(
    target,
    preds,
    label,
    property="energy",
    sampled_points=None,
    savefig=False,
    savedir="./plots/",
    filename=None,
):
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
    if savefig:
        os.makedirs(savedir, exist_ok=True)
        plt.savefig("".join([savedir, filename, ".png"]))
    plt.show()


def kde_plots(
    emt,
    forces,
    labels,
    dtype="energy",
    savefig=False,
    savedir="./plots/",
    filename=None,
):
    if dtype != "energy":
        emt = np.concatenate(np.array(emt)).flatten()
        forces = [np.concatenate(np.array(image)).flatten() for image in forces]
    fig, ax = plt.subplots(figsize=(14.15, 10))
    sns.distplot(
        emt,
        hist=False,
        kde=True,
        label="target",
        color="k",
        hist_kws={"edgecolor": "black", "alpha": 0.1},
        kde_kws={"linewidth": 2},
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
    if dtype == "energy":
        plt.xlabel("energy, eV", fontsize=32)
        plt.ylabel("Density Function", fontsize=32)
        plt.title("Energy Distributions", fontsize=35)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
    else:
        plt.xlabel("force, eV/A", fontsize=32)
        plt.ylabel("Density Function", fontsize=32)
        plt.title("Force Distributions", fontsize=35)
        plt.xlim(left=-1, right=1)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
    if savefig:
        os.makedirs(savedir, exist_ok=True)
        plt.savefig("".join([savedir, filename, ".png"]))
