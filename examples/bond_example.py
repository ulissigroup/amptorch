import sys
import numpy as np
import torch
from ase.build import molecule
from ase import Atoms
from ase.calculators.emt import EMT
import matplotlib.pyplot as plt
from amp_pytorch.NN_model import CustomLoss
from amp_pytorch import AMP
from amp_pytorch.core import AMPModel
from amp.descriptor.gaussian import Gaussian
from amp_pytorch.lj_model import lj_optim

atoms = molecule('CO')
atoms.set_pbc(False)
displacements = np.linspace(-1, 5, 100)
vec = atoms[1].position - atoms[0].position
images = []
for displacement in displacements:
    atoms = Atoms(atoms)
    atoms[1].position = (atoms[0].position + vec * displacement)
    atoms.set_calculator(EMT())
    atoms.get_potential_energy()
    images.append(atoms)

distances = []
energies = []
for image in images:
    vec = image.positions[1] - image.positions[0]
    distance = np.linalg.norm(vec)
    distances.append(distance)
    energy = image.get_potential_energy()
    energies.append(energy)

def plot(distances, energies, energy_pred, lj_energies):
    fig = plt.figure(figsize=(7., 7.))
    ax = fig.add_subplot(111)
    ax.plot(distances, energies, 'bo', markersize=2)
    ax.plot(distances, energy_pred, 'r-', lw=0.7)
    ax.plot(distances, lj_energies, 'g-', lw=0.7)
    plt.show()

torch.set_num_threads(1)
calc = AMP(model=AMPModel(images, descriptor=Gaussian(), cores=1, force_coefficient=0))
calc.model.convergence = {"energy": 0.02, "force": None}
calc.train(overwrite=True)
energy_pred = np.concatenate([calc.get_potential_energy(image) for image in images])

# p0 = [3.851, 0.105, 1, 3.5, 0.060, 1]
p0 = [.5, .5, .5, .5, .5, .5]
params_dict = {'C': [], 'O': []}
cutoff = 6.5
lj_model = lj_optim(images, p0, params_dict, cutoff)
fitted_params = lj_model.fit(method='L-BFGS-B')
lj_fitted_data = lj_model.lj_pred(fitted_params)
lj_energies = lj_fitted_data[0]

plot(distances, energies, energy_pred, lj_energies)


