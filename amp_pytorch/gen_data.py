import os
from ase import Atoms, Atom, units
import ase.io
from ase.calculators.emt import EMT
from ase.build import fcc110, fcc111
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md import VelocityVerlet
from ase.constraints import FixAtoms

from amp import Amp
from amp.descriptor.gaussian import Gaussian
from amp.model.neuralnetwork import NeuralNetwork


def generate_data(count, filename='water_Pt.traj'):
    """Generates test or training data with a simple MD simulation."""
    traj = ase.io.Trajectory(filename, 'w')
    atoms = fcc111('Pt', (2, 2, 2), vacuum=7.)
    atoms.extend(Atoms([Atom('O', atoms[7].position + (0., 0., 2.5))]))
    atoms.set_constraint(FixAtoms(indices=[0, 2]))
    atoms.set_calculator(EMT())
    atoms.get_potential_energy()
    traj.write(atoms)
    MaxwellBoltzmannDistribution(atoms, 300. * units.kB)
    dyn = VelocityVerlet(atoms, dt=1. * units.fs)
    for step in range(count - 1):
        dyn.run(50)
        traj.write(atoms)


generate_data(100)
