import sys
import ase
import numpy as np
from ase import Atoms
from ase.calculators.emt import EMT
from ase.calculators.lj import LennardJones
from ase.neighborlist import NeighborList
from amp_pytorch.neighbors import get_distances
from amp.utilities import Logger
from scipy.optimize import curve_fit, leastsq, fmin, minimize

log = Logger('lj_log.txt')
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
    MSE = np.mean((target_energies - predicted_energies) ** 2)
    return MSE

def lj_param_check(data, p0):
    unique = set()
    for atoms in data:
        symbols = atoms.get_chemical_symbols()
        unique = unique | set(symbols)
    unique_elements = list(unique)
    num_lj_params = 2 * len(unique_elements)
    assert (len(p0) == num_lj_params), 'Number of initial conditions not equal to \
    the number of required LJ parameters'

def params_to_dict(p0, params_dict):
    idx = 0
    for keys in list(params_dict.keys()):
        params_dict[keys] = p0[idx: idx + 2]
        idx +=2
    return params_dict

def logresults(log, data, cutoff, p0, params_dict, results):
    log('LJ-Parameter Optimization')
    log('Fits provided data to the Lennard Jones model')
    log('Dataset size: %s' % len(data))
    log('cutoff radius: %s' % (cutoff))
    log('inital LJ parameter guess [sig, eps]: %s' % params_to_dict(p0, params_dict))
    log('Optimizer results: \n %s \n' % results)
    # log('Fitted LJ parameters: %s \n' % params_to_dict(results.x, params_dict))


log = Logger('lj_log.txt')

images = ase.io.read('../datasets/water.extxyz',':')
IMAGES = []
for i in range(300):
    IMAGES.append(images[i])

cutoff = 5
data = IMAGES
known_energies = np.array([atoms.get_potential_energy() for atoms in data]).reshape(
    1, -1
)

# lj parameters initial conditions
params_dict = {"H": [], "O": []}
p0 = [3.5, 0.005, 3.405, 0.096]
lj_param_check(data, p0)

lj = minimize(objective_fn, p0, args=(
    data, params_dict, known_energies), method="Nelder-Mead")


logresults(log, data, cutoff, p0, params_dict, lj)
