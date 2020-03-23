from ase import Atoms, Atom, units
from ase.build import molecule
from ase.visualize import view
import ase.io
#ASE-EMT Calculator
from ase.calculators.emt import EMT
#ASE-VASP Calculator
from ase.calculators.vasp import Vasp

from ase.build import fcc100, add_adsorbate
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md import VelocityVerlet, Langevin, nvtberendsen
from ase.constraints import FixAtoms

import numpy as np


class MDsimulate:
    def __init__(self, ensemble, dt, temp, count, initial_geometry=None):
        """
        Parameters
        ----------
        ensemble: "NVE", "nvtberendsen", "langevin"
        dt: md time step (fs)
        temp: temperature (K)
        initial_slab: initial geometry to use, if None - will be generated
        """
        self.ensemble = ensemble
        self.dt = dt
        self.temp = temp * units.kB
        self.count = count
        self.starting_geometry = initial_geometry
        if initial_geometry is None:
            self.starting_geometry = self.construct_initial_geometry()

    def construct_initial_geometry(self):
        """Generates initial geometry of system"""
        slab = fcc100("Cu", size=(3, 3, 3))
        ads = molecule("CO")
        add_adsorbate(slab, ads, 2.5, offset=(1, 1))
        cons = FixAtoms(
            indices=[atom.index for atom in slab if (atom.tag == 2 or atom.tag == 3)]
        )
        slab.set_constraint(cons)
        slab.center(vacuum=13.0, axis=2)
        slab.set_pbc(True)
        slab.wrap(pbc=[True] * 3)
        return slab

    def run(self, calc, filename):
        slab = self.starting_geometry.copy()
        slab.set_calculator(calc)
        np.random.seed(1)
        MaxwellBoltzmannDistribution(slab, self.temp * units.kB)
        if self.ensemble == "NVE":
            dyn = VelocityVerlet(slab, self.dt * units.fs)
        elif self.ensemble == "nvtberendsen":
            dyn = nvtberendsen.NVTBerendsen(slab, self.dt * units.fs, self.temp, taut=300 * units.fs)
        elif self.ensemble == "langevin":
            dyn = Langevin(slab, self.dt * units.fs, self.temp * units.kB, 0.002)
        traj = ase.io.Trajectory(filename + ".traj", "w", slab)
        dyn.attach(traj.write, interval=1)

        def printenergy(a=slab):
            """Function to print( the potential, kinetic, and total energy)"""
            epot = a.get_potential_energy() / len(a)
            ekin = a.get_kinetic_energy() / len(a)
            print(
                "Energy per atom: Epot = %.3feV Ekin = %.3feV (T=%3.0fK) "
                "Etot = %.3feV" % (epot, ekin, ekin / (1.5 * units.kB), epot + ekin)
            )

        if printenergy:
            dyn.attach(printenergy, interval=10)
        dyn.run(self.count)

    def get_trajectory(self, filename, start_count, end_count, interval):
        trajectory = ase.io.read(filename+".traj", "{}:{}:{}".format(start_count, end_count, interval))
        return trajectory

def main(calculator='EMT'):
    # Define system to simulate
    slab = fcc100("Cu", size=(3, 3, 3))
    ads = molecule("CO")
    add_adsorbate(slab, ads, 2.5, offset=(1, 1))
    cons = FixAtoms(
        indices=[atom.index for atom in slab if (atom.tag == 2 or atom.tag == 3)]
    )
    slab.set_constraint(cons)
    slab.center(vacuum=13.0, axis=2)
    slab.set_pbc(True)
    slab.wrap(pbc=[True] * 3)

    #Choose calculator
    if calculator == 'EMT':
        #EMT calculator (default)
        calc = EMT()
    elif calculator == 'VASP':
        #VASP calculator
        calc = Vasp(prec='Normal',
            algo='Normal',
            xc='PBE',
            gga='RP',
            encut=400,
            lreal=False,
            ediff=1e-4,
            ispin=1,
            NELM=100,
            lwave=False,
            lcharg=False,
            nsw=0,
            kpoints=(4,4,1))
    else:
        print('Only EMT or VASP calculators allowed.')

    # Define MD settings and run MD
    md_runner = MDsimulate(ensemble="nvtberendsen", dt=1, temp=300,
            count=10000, initial_geometry=slab)
    md_runner.run(calc=calc, filename="COCu_DFT_10ps")

if __name__ == "__main__":
    main(calculator='VASP')
