import sys
import ase
import numpy as np
from ase import Atoms
from ase.calculators.emt import EMT
from ase.calculators.lj import LennardJones
from ase.neighborlist import NeighborList
from amp_pytorch.neighbors import get_distances
from scipy.optimize import curve_fit, leastsq, fmin, minimize

images = [
    Atoms(
        symbols="PdOPd",
        pbc=np.array([False, False, False], dtype=bool),
        calculator=LennardJones(),
        cell=np.array([[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]]),
        positions=np.array([[0.5, 1.0, 0.5], [1.0, 0.5, 1.0], [1.5, 1.5, 1.5]]),
    ),
    Atoms(
        symbols="PdO",
        pbc=np.array([False, False, False], dtype=bool),
        calculator=LennardJones(),
        cell=np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        positions=np.array([[2.0, 1.0, -1.0], [1.0, 2.0, 1.0]]),
    ),
]

def V_lj(params, data, params_dict):
    idx = 0
    for keys in list(params_dict.keys()):
        params_dict[keys] = params[idx: idx + 2]
        idx += 2
    predicted_energies = []
    for image in data:
        elements = image.get_chemical_symbols()
        pairs, distances = get_distances(image, cutoff)
        V = 0
        for pair, r in zip(pairs, distances[0]):
            atom_i = elements[pair[0]]
            atom_j = elements[pair[1]]
            sig_i = params_dict[atom_i][0]
            sig_j = params_dict[atom_j][0]
            eps_i = params_dict[atom_i][1]
            eps_j = params_dict[atom_j][1]
            # Combining Rules
            sig = (sig_i + sig_j) / 2
            eps = np.sqrt(eps_i * eps_j)
            c6 = (sig/r)**6
            c12 = c6**2
            V_lj = 4 * eps * (c12 - c6)
            V += V_lj
        V /= 2
        predicted_energies.append(V)
    predicted_energies = np.array(predicted_energies).reshape(1, -1)
    return predicted_energies


def objective_fn(params, data, params_dict, target_energies):
    predicted_energies = V_lj(params, data, params_dict)
    print(predicted_energies)
    print(target_energies)
    MSE = np.mean((target_energies - predicted_energies) ** 2)
    return MSE

cutoff = 3
data = [images[1]]
known_energies = np.array([atoms.get_potential_energy() for atoms in data]).reshape(
    1, -1
)

unique = set()
for atoms in data:
    symbols = atoms.get_chemical_symbols()
    unique = unique | set(symbols)
unique_elements = list(unique)
num_lj_params = 2 * len(unique_elements)

# lj parameters initial conditions
params_dict = {"Pd": [], "O": []}
p0 = [1, 1, 1, 1]

assert (len(p0) == num_lj_params), 'Number of initial conditions not equal to \
the number of required LJ parameters'

objective_fn(p0, data, params_dict, known_energies)
# lj = minimize(objective_fn, p0, args=(
    # data, params_dict, known_energies), method="Nelder-Mead")

# idx = 0
# for keys in list(params_dict.keys()):
    # params_dict[keys] = lj.x[idx: idx + 2]
    # idx += 2

# print(params_dict)
