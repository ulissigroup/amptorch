import os
import time
import sys
import numpy as np
from ase import Atoms, Atom, units
from ase.build import molecule
from ase.visualize import view
import ase.io
from ase.calculators.emt import EMT
import matplotlib.pyplot as plt

# from asap3 import EMT
from ase.build import fcc110, fcc111, fcc100, add_adsorbate, molecule
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md import VelocityVerlet, Langevin
from ase.optimize import QuasiNewton
from ase.constraints import FixAtoms, Hookean

from amp import Amp
from amp.descriptor.gaussian import Gaussian
from amp.model.neuralnetwork import NeuralNetwork


def generate_data(count, filename, temp, hook, cons_t=False):
    """Generates test or training data with a simple MD simulation."""
    traj = ase.io.Trajectory(filename, "w")
    pair = molecule("CO")
    cons = Hookean(a1=0, a2=1, rt=1.58, k=10.0)
    # pair.set_constraint(cons)
    pair.set_calculator(EMT())
    pair.get_potential_energy()
    dyn = QuasiNewton(pair, trajectory=(filename[:-5] + "_relax.traj"))
    dyn.run(fmax=0.05)
    traj.write(pair)
    MaxwellBoltzmannDistribution(pair, temp * units.kB)
    if cons_t is True:
        dyn = Langevin(pair, 5 * units.fs, temp * units.kB, 0.002)
    else:
        dyn = VelocityVerlet(pair, dt=1.0 * units.fs)
    for step in range(count - 1):
        dyn.run(50)
        traj.write(pair)


generate_data(100, "pairs.traj", temp=300, hook=False, cons_t=False)
k = ase.io.read("pairs.traj", ":")
distances = []
energies = []
for image in k:
    pos1 = image.positions
    vec = pos1[1]-pos1[0]
    d = np.sqrt((vec**2).sum())
    distances.append(d)
    energy = image.get_potential_energy()
    energies.append(energy)
distances.sort()
energies.sort()
plt.plot(distances, energies)
plt.show()
