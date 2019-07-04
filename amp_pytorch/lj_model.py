import sys
import time
import numpy as np
from scipy.optimize import minimize
import ase
from ase import Atoms
from ase.build import molecule
from ase.calculators.emt import EMT
from ase.calculators.lj import LennardJones
from ase.neighborlist import NeighborList, NewPrimitiveNeighborList
from amp_pytorch.neighbors import get_distances
from amp.utilities import Logger

log = Logger("lj_log.txt")


def V_lj(params, data, params_dict):
    params_dict = params_to_dict(params, params_dict)
    predicted_energies = []
    predicted_forces = []
    num_atoms = []
    for image in data:
        natoms = len(image)
        num_atoms.append([natoms] * natoms)
        elements = image.get_chemical_symbols()
        pairs, distances = get_distances(image, [(cutoff + 1e-10) / 2] * natoms)
        V = 0
        forces = np.zeros((natoms, 3))
        for pair, d in zip(pairs, distances):
            atom_i = elements[pair[0]]
            atom_j = elements[pair[1]]
            sig_i = params_dict[atom_i][0]
            sig_j = params_dict[atom_j][0]
            eps_i = params_dict[atom_i][1]
            eps_j = params_dict[atom_j][1]
            offset_i = params_dict[atom_i][2]
            offset_j = params_dict[atom_j][2]
            # Combining Rules
            r = np.sqrt((d ** 2).sum())
            sig = (sig_i + sig_j) / 2
            eps = np.sqrt(eps_i * eps_j)
            c6 = (sig / r) ** 6
            c12 = c6 ** 2
            V_lj = 4 * eps * (c12 - c6)
            V_lj += offset_i + offset_j
            f = (24 * eps * (2 * c12 - c6) / r ** 2) * d
            forces[pair[0]] -= f
            forces[pair[1]] += f
            V += V_lj
        V /= 2
        forces /= 2
        predicted_energies.append(V)
        predicted_forces.append(forces)
    predicted_energies = np.array(predicted_energies).reshape(1, -1)
    predicted_forces = np.concatenate(predicted_forces)
    num_atoms = np.concatenate(num_atoms).reshape(-1, 1)
    return predicted_energies, predicted_forces, num_atoms


def objective_fn(params, data, params_dict, target_energies, target_forces):
    predicted_energies, predicted_forces, num_atoms = V_lj(params, data, params_dict)
    MSE_energy = np.mean((target_energies - predicted_energies) ** 2)
    MSE_forces = np.mean((target_forces - predicted_forces) ** 2)
    MSE = MSE_energy + MSE_forces
    return MSE


def lj_param_check(data, p0):
    unique = set()
    for atoms in data:
        symbols = atoms.get_chemical_symbols()
        unique = unique | set(symbols)
    unique_elements = list(unique)
    num_lj_params = 3 * len(unique_elements)
    assert (
        len(p0) == num_lj_params
    ), "Number of initial conditions not equal to \
    the number of required LJ parameters"


def params_to_dict(p0, params_dict):
    idx = 0
    for keys in list(params_dict.keys()):
        params_dict[keys] = p0[idx : idx + 3]
        idx += 3
    return params_dict


def logresults(log, data, cutoff, p0, params_dict, results):
    log("-" * 50)
    log("LJ-Parameter Optimization")
    log("Fits provided data to the Lennard Jones model")
    log("%s" % time.asctime())
    log("Dataset size: %s" % len(data))
    log("cutoff radius: %s" % (cutoff))
    log("inital LJ parameter guess [sig, eps]: %s" % params_to_dict(p0, params_dict))
    log("Optimizer results: \n %s \n" % results)
    # log('Fitted LJ parameters: %s \n' % params_to_dict(results.x, params_dict))


log = Logger("lj_log.txt")
cutoff = 3

images = ase.io.read("../datasets/water.extxyz", ":")
IMAGES = []
for i in range(10):
    IMAGES.append(images[i])

data = IMAGES
known_energies = np.array([(atoms.get_potential_energy()) for atoms in data]).reshape(
    1, -1
)
known_forces = np.concatenate([(atoms.get_forces()) for atoms in data])

# lj parameters initial conditions
params_dict = {"H": [], "O": []}
# p0 structure:  [sigma_element_1, epsilon_element_1, offset_element1, sigma_element_2,...]
p0 = [1, 1, -2700, 1, 1, -2700]
lj_param_check(data, p0)

pred_energy, pred_forces, num_atoms = V_lj(p0, data, params_dict)

lj = minimize(
    objective_fn,
    p0,
    args=(data, params_dict, known_energies, known_forces),
    method="Nelder-Mead",
)

logresults(log, data, cutoff, p0, params_dict, lj)
