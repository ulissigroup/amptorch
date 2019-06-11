import sys
import ase
import numpy as np
from ase import Atoms
from ase.calculators.emt import EMT
from amp_pytorch.neighbors import get_distances
from scipy.optimize import curve_fit, leastsq, fmin, minimize


images = [
    Atoms(
        symbols="PdOPd2",
        pbc=np.array([False, False, False], dtype=bool),
        calculator=EMT(),
        cell=np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        positions=np.array(
            [[0.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0], [1.0, 0.0, 0.0]]
        ),
    ),
    Atoms(
        symbols="PdOPd2",
        pbc=np.array([False, False, False], dtype=bool),
        calculator=EMT(),
        cell=np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        positions=np.array(
            [[0.0, 1.0, 0.0], [1.0, 2.0, 1.0], [-1.0, 1.0, 2.0], [1.0, 3.0, 2.0]]
        ),
    ),
    Atoms(
        symbols="PdO",
        pbc=np.array([False, False, False], dtype=bool),
        calculator=EMT(),
        cell=np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        positions=np.array([[2.0, 1.0, -1.0], [1.0, 2.0, 1.0]]),
    ),
    Atoms(
        symbols="Pd2O",
        pbc=np.array([False, False, False], dtype=bool),
        calculator=EMT(),
        cell=np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        positions=np.array([[-2.0, -1.0, -1.0], [1.0, 2.0, 1.0], [3.0, 4.0, 4.0]]),
    ),
    Atoms(
        symbols="Cu",
        pbc=np.array([False, False, False], dtype=bool),
        calculator=EMT(),
        cell=np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        positions=np.array([[0.0, 0.0, 0.0]]),
    ),
]

cutoff = 3
data = [images[0], images[1], images[-1]]
energies = [atoms.get_potential_energy() for atoms in data]

unique = set()
for atoms in data:
    symbols = atoms.get_chemical_symbols()
    unique = unique | set(symbols)
unique_elements = list(unique)
num_lj_params = 2 * len(unique_elements)

params_dict = {"Pd": [], "O": [], "Cu": []}
p0 = [1, 1, 1, 1, 1, 1]


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
            V_lj = 4 * eps * ((sig / r) ** 12 - (sig / r) ** 6)
            V += V_lj
        predicted_energies.append(V)
    return predicted_energies 

V_lj(p0, data, params_dict)
sys.exit()


def obj(params, pairs, distances, params_dict):
    known_energies = energies
    err = known_energies - V_lj(params, pairs, distances, params_dict)
    return np.mean(err ** 2)


lj = minimize(obj, p0, args=(pairs, distances, params_dict), method="Nelder-Mead")
print(lj)
