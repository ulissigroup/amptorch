import ase.db
db = ase.db.connect('argon.db')

known_energies = [row.energy for row in db.select()]

import numpy as np
from scipy.optimize import fmin
from ase.calculators.lj import LennardJones

def my_lj(pars):
    epsilon, sigma = pars
    calc = LennardJones(sigma=sigma, epsilon=epsilon)
    all_atoms = [row.toatoms() for row in db.select()]
    [atoms.set_calculator(calc) for atoms in all_atoms]
    predicted_energies = np.array([atoms.get_potential_energy() for atoms in all_atoms])
    return predicted_energies

def objective(pars):
    known_energies = np.array([row.energy for row in db.select()])
    err = known_energies - my_lj(pars)
    return np.mean(err**2)

LJ_pars = fmin(objective, [0.005, 3.5])
print(LJ_pars)
