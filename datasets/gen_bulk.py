import os
import time
import sys
import numpy as np
from ase import Atoms, Atom, units
from ase.build import molecule
from ase.visualize import view
import ase.io
from ase.calculators.emt import EMT

from ase.build import fcc110, fcc111, fcc100, add_adsorbate, molecule, bulk
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md import VelocityVerlet, Langevin, nvtberendsen
from ase.optimize import QuasiNewton
from ase.constraints import FixAtoms, Hookean

from amp import Amp
from amp.descriptor.gaussian import Gaussian
from amp.model.neuralnetwork import NeuralNetwork



def generate_data(count, filename):
    """Generates test or training data with a simple MD simulation."""
    traj = ase.io.Trajectory(filename, "w")
    lc = np.linspace(3, 4.4, count)
    for a in lc:
        bulk_data = bulk('Cu', 'fcc', a)
        bulk_data.set_calculator(EMT())
        traj.write(bulk_data)


generate_data(100, "bulk.traj")


