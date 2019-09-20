from ase.build import molecule
from ase.calculators.emt import EMT
from ase.optimize import QuasiNewton
from ase.vibrations import Vibrations
from ase.thermochemistry import IdealGasThermo
from ase import io
from ase.units import kB
import sys


# data = io.read("../datasets/COCu/COCu.traj",":")
atoms = molecule('N2')
atoms.set_calculator(EMT())
dyn = QuasiNewton(atoms)
dyn.run(fmax=0.01)
potentialenergy = atoms.get_potential_energy()

vib = Vibrations(atoms)
vib.run()
vib_energies = vib.get_energies()

thermo = IdealGasThermo(vib_energies=vib_energies,
                        potentialenergy=potentialenergy,
                        atoms=atoms,
                        geometry='linear',
                        symmetrynumber=2, spin=0)
G = thermo.get_gibbs_energy(temperature=298.15, pressure=101325.)
