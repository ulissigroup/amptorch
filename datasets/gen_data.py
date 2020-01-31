import os
import time
import sys
import numpy as np
from ase import Atoms, Atom, units
from ase.build import molecule
from ase.visualize import view
import ase.io
from ase.calculators.emt import EMT

from ase.build import fcc110, fcc111, fcc100, add_adsorbate, molecule
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md import VelocityVerlet, Langevin, nvtberendsen
from ase.optimize import QuasiNewton
from ase.constraints import FixAtoms, Hookean

from amp import Amp
from amp.descriptor.gaussian import Gaussian
from amp.model.neuralnetwork import NeuralNetwork


def generate_data(count, filename, temp, hook, dt, ensemble="NVE"):
    """Generates test or training data with a simple MD simulation."""
    slab = fcc100("Cu", size=(3, 3, 3))
    ads = molecule("CO")
    add_adsorbate(slab, ads, 4, offset=(1, 1))
    cons = FixAtoms(
        indices=[atom.index for atom in slab if (atom.tag == 2 or atom.tag == 3)]
    )
    if hook:
        cons2 = Hookean(a1=28, a2=27, rt=1.58, k=10.0)
        slab.set_constraint([cons, cons2])
    else:
        slab.set_constraint(cons)
    slab.center(vacuum=13.0, axis=2)
    slab.set_pbc(True)
    slab.wrap(pbc=[True] * 3)
    slab.set_calculator(EMT())
    slab.get_forces()
    MaxwellBoltzmannDistribution(slab, temp * units.kB)
    if ensemble == "NVE":
        dyn = VelocityVerlet(slab, dt=dt * units.fs)
    elif ensemble == "nvtberendsen":
        dyn = nvtberendsen.NVTBerendsen(slab, dt * units.fs, temp, taut=300 * units.fs)
    elif ensemble == "langevin":
        dyn = Langevin(slab, dt * units.fs, temp * units.kB, 0.002)
    traj = ase.io.Trajectory(filename, "w", slab)
    dyn.attach(traj.write, interval=1)

    def printenergy(a=slab):
        epot = a.get_potential_energy()
        ekin = a.get_kinetic_energy()
        print(
            "Energy per atom: Epot = %.3feV Ekin = %.3feV (T=%3.0fK) "
            "Etot = %.3feV" % (epot, ekin, ekin / (1.5 * units.kB), epot + ekin)
        )

    dyn.attach(printenergy, interval=10)
    dyn.run(count - 1)


generate_data(2000, "COCu_lang_1fs_300K.traj", temp=300.0, hook=False, dt=1, ensemble="langevin")
