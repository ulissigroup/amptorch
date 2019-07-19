import os
import time
import sys
import numpy as np
from ase import Atoms, Atom, units
from ase.visualize import view
import ase.io
from ase.calculators.emt import EMT

# from asap3 import EMT
from ase.build import fcc110, fcc111, fcc100, add_adsorbate, molecule
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md import VelocityVerlet, Langevin
from ase.optimize import QuasiNewton
from ase.constraints import FixAtoms, Hookean

from amp import Amp
from amp.descriptor.gaussian import Gaussian
from amp.model.neuralnetwork import NeuralNetwork


def generate_data(count, filename, temp, cons_t=False):
    """Generates test or training data with a simple MD simulation."""
    traj = ase.io.Trajectory(filename, "w")
    slab = fcc100("Cu", size=(3, 3, 3))
    atom = "He"
    atom = Atoms(atom).repeat(10)
    view(atom)
    view(atom)
    sys.exit()
    slab.set_constraint(cons)
    slab.set_calculator(EMT())
    slab.get_potential_energy()
    dyn = QuasiNewton(slab, trajectory=(filename[:-5] + "_relax.traj"))
    dyn.run(fmax=0.05)
    traj.write(slab)
    MaxwellBoltzmannDistribution(slab, temp * units.kB)
    if cons_t is True:
        dyn = Langevin(slab, 5 * units.fs, temp * units.kB, 0.002)
    else:
        dyn = VelocityVerlet(slab, dt=1.0 * units.fs)
    for step in range(count - 1):
        dyn.run(50)
        traj.write(slab)


generate_data(500, "COCu/COCu_conT.traj", temp=300.0, cons_t=True)
generate_data(500, "COCu/COCu.traj", temp=300.0, cons_t=False)
