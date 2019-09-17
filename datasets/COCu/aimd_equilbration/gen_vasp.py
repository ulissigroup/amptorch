from ase.io import read
from ase.units import kJ
from ase.eos import EquationOfState
from ase.calculators.vasp import Vasp
import ase.calculators.vasp as vasp_calculator
from ase.io.trajectory import Trajectory
from ase.io.trajectory import PickleTrajectory
from ase import Atoms 
from ase import Atom
from ase.io import read,write
from ase.constraints import FixAtoms
from ase.optimize import BFGSLineSearch
from ase.optimize import BFGS
import numpy as np
import os

xc = 'PBE' 
### Find Trajectory files or POSCAR
atoms=read('input.POSCAR')




### VASP details
calc = vasp_calculator.Vasp(
			kpar=1,		#2kpoints
			ncore=4,		#2node jobs	
			encut=400,
			xc='PBE',
                        gga='RP',
			kpts=(3,3,1), 
			gamma = False, # Gamma-centered (defaults to Monkhorst-Pack)
			ismear=0,
            		algo = 'fast',
                        nelm=1000,
            		sigma = 0.05,
			ibrion=0,
			lorbit=11,
			potim=0.5,
#			isif=2,
			ediffg=-0.05,  # forces
			ediff=1e-4,  #energy conv. both of these are for the internal relaxation, ie 
			prec='normal',
			lcharg=False,
			lwave=False,
			ispin=1,
			nsw=2000,  # total 1 ps of equilibration
			isym=0,
			smass = 0,
			tebeg = 300,
			teend = 300,
			)


atoms.set_calculator(calc)
energy=atoms.get_potential_energy()
write('final.traj',atoms)
