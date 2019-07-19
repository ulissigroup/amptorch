import sys
import numpy as np
import torch
from ase.build import molecule
from ase import Atoms
from ase.calculators.emt import EMT
from ase.calculators.lj import LennardJones
from ase import units
import ase.io
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md import VelocityVerlet, Langevin
import matplotlib.pyplot as plt
from amp_pytorch.NN_model import CustomLoss
from amp_pytorch import AMP
from amp_pytorch.core import AMPModel
from amp.descriptor.gaussian import Gaussian
from amp_pytorch.lj_model import lj_optim

distances = []
energies = []
forces = []
distances = np.linspace(2, 5, 100)
images = []
for displacement in distances:
    atoms = Atoms("He2", [(0, 0, 0), (displacement, 0, 0)])
    atoms.set_calculator(LennardJones(sigma=2.2362, epsilon=0.056 * 0.043372093))
    # atoms.set_calculator(EMT())
    atoms.get_potential_energy()
    energy = atoms.get_potential_energy()
    force = atoms.get_forces()
    energies.append(energy)
    forces.append(force)
    images.append(atoms)

forces = np.array([np.amax(np.abs(force)) for force in forces]).reshape(-1, 1)

traj = ase.io.Trajectory("test.traj", "w")
atoms = images[0]
atoms.set_calculator(LennardJones(sigma=2.2362, epsilon=0.056 * 0.043372093))
atoms.get_potential_energy()
MaxwellBoltzmannDistribution(atoms, 300.0 * units.kB)
dyn = VelocityVerlet(atoms, dt=1.0 * units.fs)
for step in range(20):
    dyn.run(50)
    traj.write(atoms)
sys.exit()


def plot(distances, energies, energy_pred, lj_energies):
    lj_energies = lj_energies.reshape(len(energy_pred), 1)
    fig = plt.figure(figsize=(7.0, 7.0))
    ax = fig.add_subplot(111)
    ax.plot(distances, energies, "bo", label="LJ", markersize=2)
    ax.plot(distances, energy_pred, "r-", label="NN_LJ", lw=0.7)
    ax.plot(distances, lj_energies, "g-", label="LJ_optim", lw=0.7)
    ax.legend()
    plt.show()


p0 = [2.362, 0.056 * 0.043372093, 0]
params_dict = {"He": []}
cutoff = 6.5
lj_model = lj_optim(images, p0, params_dict, cutoff)
fitted_params = lj_model.fit()
lj_energies, lj_forces, num_atoms = lj_model.lj_pred(images, fitted_params, params_dict)
lj_data = [lj_energies, lj_forces, num_atoms, fitted_params, params_dict, lj_model]

torch.set_num_threads(1)
calc = AMP(
    model=AMPModel(
        images, descriptor=Gaussian(), cores=1, force_coefficient=0.3, lj_data=lj_data
    ),
    label="amptorch_bond.pt",
)
calc.model.convergence = {"energy": 0.0002, "force": 0.002}
calc.train(overwrite=True)
energy_pred = np.concatenate([calc.get_potential_energy(image) for image in images])

plot(distances, energies, energy_pred, lj_energies)
