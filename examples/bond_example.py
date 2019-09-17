import sys
import numpy as np
import torch
from ase.build import molecule
from ase import Atoms, Atom
from ase.calculators.emt import EMT
from ase.calculators.lj import LennardJones
from ase import units
import ase.io
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.optimize import QuasiNewton
from ase.md import VelocityVerlet, Langevin
import matplotlib.pyplot as plt
from amp_pytorch.NN_model import CustomLoss
from amp_pytorch import AMP
from amp_pytorch.core import AMPModel
from amp.descriptor.gaussian import Gaussian
from amp_pytorch.lj_model import lj_optim
from ase.visualize import view
import seaborn as sns
import random
import scipy

distances = []
energies = []
forces = []
distances = np.linspace(1, 5, 100)
images = []

offset_atoms = Atoms("CuCu", [(0, 0, 0), (0, 0, 1000000)])
offset_atoms.set_calculator(EMT())
e0 = offset_atoms.get_potential_energy()
for displacement in distances:
    atoms = Atoms("CuCu", [(0, 0, 0), (0, 0, displacement)])
    atoms.set_calculator(EMT())
    atoms.get_potential_energy()
    energy = atoms.get_potential_energy() - e0
    force = atoms.get_forces()
    energies.append(energy)
    forces.append(force)
    images.append(atoms)

energies = np.array(energies)
forces = np.array([np.amax(np.abs(force)) for force in forces]).reshape(-1, 1)

def md_run(images, count, filename):
    traj = ase.io.Trajectory(filename, "w")
    image = images[0]
    image.set_calculator(EMT())
    image.get_potential_energy()
    dyn = QuasiNewton(image, trajectory=(filename[:-5] + "_relax.traj"))
    dyn.run(fmax=0.05)
    traj.write(image)
    MaxwellBoltzmannDistribution(image, 400.0 * units.kB)
    # dyn = Langevin(image, 5*units.fs, 400.0*units.kB, 0.002)
    dyn = VelocityVerlet(image, dt=1.0 * units.fs)
    for step in range(count):
        dyn.run(50)
        traj.write(image)

# training_idx = np.array(random.sample(range(len(k)), 20))
training_idx = np.array([8, 10, 12, 14, 16, 18, 20, 22, 25, 27])
training_images = [images[i] for i in training_idx]
training_energies = energies[training_idx]
training_distances = distances[training_idx]


def train(train_images, test_images, lj=False):
    lj_data = None
    label = "bond_ml.pt"
    if lj:
        label = "bond_mllj.pt"
        eV_kcalmol = 0.043372093
        # p0 = [1.8, 5, 0, 1.8, 5, 0]
        p0 = [1.8, 5, 1]
        p0 = [1.38147567, -0.3407247, -2.15276909]
        # p0 = [1.198875, 29.73048, 0]
        params_dict = {"Cu": [], "C": []}
        cutoff = 10
        lj_model = lj_optim(train_images, p0, params_dict, cutoff)
        # fitted_params = lj_model.fit(method="L-BFGS-B")
        fitted_params = lj_model.fit()
        # fitted_params = p0
        lj_energies, lj_forces, num_atoms = lj_model.lj_pred(
            train_images, fitted_params, params_dict
        )
        # lj_energies, lj_forces, num_atoms = lj_model.lj_pred(
        # test_images, fitted_params, params_dict
        # )
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
        model=AMPModel(
            train_images,
            descriptor=Gaussian(),
            cores=1,
            force_coefficient=0.3,
            lj_data=lj_data,
        ),
        label=label,
    )
    calc.model.convergence = {"energy": 0.001, "force": 0.001}
    calc.model.lr = 0.1
    calc.train(overwrite=True)
    energy_pred = np.concatenate(
        [calc.get_potential_energy(image) - e0 for image in test_images]
    )
    # force_pred = np.array(
    # [np.amax(np.abs(calc.get_forces(image))) for image in test_images]
    # ).reshape(-1, 1)
    force_pred = 0
    return energy_pred, force_pred


# ml_energy, ml_forces = train(training_images, images, lj=False)
mllj_energy, mllj_forces = train(training_images, images, lj=True)
kb = 8.617e-5
T = 10000
en = np.exp(-1 * energies / (kb * T))
en_sum = en.sum()
boltz = en / en_sum

def plot(name, ml_energy=None, lj_energy=None):
    fig, ax1 = plt.subplots(figsize=(15, 10))
    ax1.set_xlabel(r"Distance, $\AA$", fontsize=32)
    ax1.set_ylabel("Energy, eV", color="k", fontsize=32)
    ax1.plot(distances, energies, color="k", linewidth=3)
    ax1.plot(
        training_distances, training_energies, "ro", label="train data",
        markersize=8
    )
    ax1.tick_params(axis="y", labelcolor="k")
    # ax2 = ax1.twinx()
    # ax2.set_ylabel("Probability", color="r", fontsize=32)
    # ax2.plot(distances, boltz, "r--")
    # ax2.tick_params(axis="y", labelcolor="r")
    if ml_energy is not None:
        ax1.plot(distances, ml_energy, "m", label="ML", linewidth=3)
    if lj_energy is not None:
        ax1.plot(distances, mllj_energy, "g", label="ML-LJ", linewidth=3)
    ax1.set_ylim(top=5)
    ax1.set_ylim(bottom=-5)
    ax1.tick_params(axis='both', labelsize=28)
    # ax2.tick_params(axis='both', labelsize=28)
    # plt.legend(fontsize=18)
    ax1.annotate('actual', (4, 0), color='k', fontsize=30)
    # ax1.annotate('Boltz. dist.', (1.9, 2), color='r', fontsize=30)
    ax1.annotate('train data', (2, -4.7), color='r', fontsize=20)
    if lj_energy is not None:
        ax1.annotate('ML-LJ', (4, -1.2), color='g', fontsize=30)
        ax1.annotate('ML', (4, -2.9), color='m', fontsize=30)
    elif ml_energy is not None:
        ax1.annotate('ML', (4, -2.9), color='m', fontsize=30)
    plt.savefig(name, dpi=300)
    plt.show()

# plot('ml_only2.png')
# plot('demo/ml.png', ml_energy,  None)
plot('demo/ml_lj.png', None, mllj_energy)
