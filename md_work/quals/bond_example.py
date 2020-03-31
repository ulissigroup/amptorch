from ase.visualize import view
from amptorch.lj_model import lj_optim
from amptorch.lj_new import lj_optim as lj_optim2
from amptorch.gaussian import Gaussian
from amptorch.core import AMPTorch
from amptorch import AMP
from amptorch.model import CustomLoss
from ase.md import VelocityVerlet, Langevin
from ase.optimize import QuasiNewton
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase import units
from ase.calculators.lj import LennardJones
from ase.calculators.emt import EMT
from ase import Atoms, Atom
from ase.build import molecule
import torch
import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use("Agg")

distances = []
energies = []
forces = []
distances = np.linspace(1, 5, 100)
images = []

offset_atoms = Atoms("CuC", [(0, 0, 0), (0, 0, 1000000)])
offset_atoms.set_calculator(EMT())
e0 = offset_atoms.get_potential_energy()
for displacement in distances:
    atoms = Atoms("CuC", [(0, 0, 0), (0, 0, displacement)])
    atoms.set_cell([10, 10, 10])
    atoms.wrap(pbc=True)
    atoms.set_calculator(EMT())
    energy = atoms.get_potential_energy()
    force = atoms.get_forces()
    energies.append(energy)
    forces.append(force)
    images.append(atoms)

energies = np.array(energies)
forces = np.array([np.amax(np.abs(force)) for force in forces]).reshape(-1, 1)

# training_idx = np.array(random.sample(range(len(k)), 20))
training_idx = np.array([8, 10, 12, 14, 16, 18, 20, 22, 25, 27])
training_images = [images[i] for i in training_idx]
training_energies = energies[training_idx]
training_distances = distances[training_idx]


def train(train_images, test_images, lj=False):
    lj_data = None
    label = "bond_ml.pt"
    lj_test = None
    if lj:
        label = "bond_mllj.pt"
        # p0 = [2.29284676, 0.29639983, 0.08071821, 1.33905162, 0.12290683, 6.41914719]
        p0 = [2.29284676, 0.29639983, 1.33905162, 0.12290683]
        params_dict = {"Cu": [], "C": []}
        element_energies = {}
        for element in ["C", "Cu"]:
            atoms = Atoms(element, cell=[20, 20, 20])
            atoms.set_calculator(EMT())
            element_energies[element] = atoms.get_potential_energy()
        cutoff = 10
        # lj_model = lj_optim(train_images, p0, params_dict, cutoff, "test")
        lj_model = lj_optim2(
            train_images, p0, params_dict, cutoff, "test", element_energies
        )
        fitted_params = lj_model.fit()
        lj_energies, lj_forces, num_atoms = lj_model.lj_pred(
            train_images, fitted_params, params_dict
        )
        # lj_model = lj_optim(test_images, p0, params_dict, cutoff, "test")
        lj_model = lj_optim2(
            test_images, p0, params_dict, cutoff, "test", element_energies
        )
        lj_test, _, _ = lj_model.lj_pred(test_images, fitted_params, params_dict)
        lj_data = [
            lj_energies,
            lj_forces,
            num_atoms,
            fitted_params,
            params_dict,
            lj_model,
        ]

    torch.set_num_threads(1)
    calc = AMP(
        model=AMPTorch(
            train_images,
            descriptor=Gaussian,
            Gs=Gs,
            cores=1,
            force_coefficient=0.3,
            lj_data=lj_data,
            label=label,
            save_logs=False,
        )
    )
    calc.model.lr = 1e-3
    calc.model.convergence = {
        "energy": 0.001,
        "force": 0.001,
        "epochs": 100,
        "early_stop": True,
    }
    calc.model.loader_params = {"shuffle": True, "batch_size": None}
    calc.model.structure = [5, 5]
    calc.model.val_frac = 0
    calc.train(overwrite=True)
    energy_pred = [calc.get_potential_energy(image) for image in test_images]
    return energy_pred, lj_test


Gs = {}
Gs["G2_etas"] = np.logspace(np.log10(0.05), np.log10(5.0), num=4)
Gs["G2_rs_s"] = [0] * 4
Gs["G4_etas"] = [0.005]
Gs["G4_zetas"] = [1.0, 4.0]
Gs["G4_gammas"] = [+1.0, -1.0]
Gs["cutoff"] = 10

ml_energy, _ = train(training_images, images, lj=False)
mllj_energy, lj_energies = train(training_images, images, lj=True)
kb = 8.617e-5
T = 10000
en = np.exp(-1 * energies / (kb * T))
en_sum = en.sum()
boltz = en / en_sum


def plot(name, ml_energy=None, mllj_energy=None, lj_energies=None):
    fig, ax1 = plt.subplots(figsize=(15, 10))
    ax1.set_xlabel(r"Distance, $\AA$", fontsize=32)
    ax1.set_ylabel("Energy, eV", color="k", fontsize=32)
    ax1.plot(distances, energies, color="k", linewidth=3)
    ax1.plot(
        training_distances, training_energies, "ro", label="train data", markersize=8
    )
    ax1.tick_params(axis="y", labelcolor="k")
    # ax2 = ax1.twinx()
    # ax2.set_ylabel("Probability", color="r", fontsize=32)
    # ax2.plot(distances, boltz, "r--")
    # ax2.tick_params(axis="y", labelcolor="r")
    if ml_energy is not None:
        ax1.plot(distances, ml_energy, "m", label="ML", linewidth=3)
    if mllj_energy is not None:
        ax1.plot(distances, mllj_energy, "g", label="ML-LJ", linewidth=3)
    if lj_energies is not None:
        ax1.plot(distances, lj_energies, "b", label="LJ", linewidth=3)
    ax1.set_ylim(top=10)
    ax1.set_ylim(bottom=0)
    ax1.tick_params(axis="both", labelsize=28)
    # ax2.tick_params(axis='both', labelsize=28)
    plt.legend(fontsize=18)
    # ax1.annotate("actual", (4, 0), color="k", fontsize=30)
    # ax1.annotate('Boltz. dist.', (1.9, 2), color='r', fontsize=30)
    # ax1.annotate("train data", (2, -4.7), color="r", fontsize=20)
    # if lj_energy is not None:
    # ax1.annotate("ML-LJ", (4, -1.2), color="g", fontsize=30)
    # ax1.annotate("ML", (4, -2.9), color="m", fontsize=30)
    # elif ml_energy is not None:
    # ax1.annotate("ML", (4, -2.9), color="m", fontsize=30)
    plt.savefig(name, dpi=300)
    plt.show()


plot("test.png", ml_energy, mllj_energy, lj_energies)
