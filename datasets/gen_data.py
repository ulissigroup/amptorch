import os
import sys
import numpy as np
from ase import Atoms, Atom, units
from ase.visualize import view
import ase.io
from ase.calculators.emt import EMT
from ase.build import fcc110, fcc111, add_adsorbate, molecule
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md import VelocityVerlet
from ase.constraints import FixAtoms, Hookean

from amp import Amp
from amp.descriptor.gaussian import Gaussian
from amp.model.neuralnetwork import NeuralNetwork


def generate_data(count, filename="H2_Pt.traj"):
    """Generates test or training data with a simple MD simulation."""
    traj = ase.io.Trajectory(filename, "w")
    slab = fcc111("Pt", (3, 3, 3))
    ads = molecule("H2")
    ads.rotate(90, "y")
    add_adsorbate(slab, ads, 4, "ontop")
    add_adsorbate(slab, ads, 4, "ontop", offset=(2, 2))
    slab.center(vacuum=10.0, axis=2)
    cons = FixAtoms(
        indices=[atom.index for atom in slab if (atom.tag == 2 or atom.tag == 3)]
    )
    cons2 = Hookean(a1=27, a2=28, rt=2, k=30)
    cons3 = Hookean(a1=29, a2=30, rt=2, k=30)
    const = [cons, cons2, cons3]
    slab.set_constraint(const)
    slab.set_calculator(EMT())
    slab.get_potential_energy()
    traj.write(slab)
    MaxwellBoltzmannDistribution(slab, 300.0 * units.kB)
    dyn = VelocityVerlet(slab, dt=1.0 * units.fs)
    for step in range(count):
        dyn.run(20)
        traj.write(slab)
        # printenergy(slab)


def printenergy(a):  # store a reference to atoms in the definition.
    """Function to print the potential, kinetic and total energy."""
    epot = a.get_potential_energy() / len(a)
    ekin = a.get_kinetic_energy() / len(a)
    print(
        "Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  "
        "Etot = %.3feV" % (epot, ekin, ekin / (1.5 * units.kB), epot + ekin)
    )


generate_data(1000)
# k = ase.io.read('H2_Pt.traj', ':')
# for image in k:
    # view(image)


