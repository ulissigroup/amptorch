"""An example of how to utilize the package to train on energies and forces"""

import sys
import time
from ase import Atoms
import ase
import ase.db
import numpy as np
from amp_pytorch.lj_model import lj_optim
import ase.io
from ase.calculators.emt import EMT
from ase.calculators.lj import LennardJones


IMAGES = ase.io.read("../datasets/water.extxyz", ":")
energies = []
images = []
cutoff = 10
p0 = [1, 2, 2, 3]

image = IMAGES[0]
image.set_calculator(EMT())
images.append(image)
energies.append(image.get_potential_energy())

HH_distance = 1.6979939322296198
OH_distance = 0.9381181946019991
HO_distance = 0.7598757376276198

HO_sig = 3 / 2.0
HO_eps = np.sqrt(6)


def lj(distance, sig, eps):
    c6 = (sig**2 / distance) ** 3
    c12 = c6 ** 2
    energy = 4 * eps * (c12 - c6)
    return energy

energy = (
    lj(OH_distance, HO_sig, HO_eps)
    + lj(HO_distance, HO_sig, HO_eps)
    + lj(HH_distance, 1, 2)
    + 7.81 + 7.81 + 6.42
)

params_dict = {"H": [], "O": []}
lj_model = lj_optim(images, p0, params_dict, cutoff)
lj_energies, lj_forces, num_atoms = lj_model.lj_pred(images, p0, params_dict)

assert round(energy, 5)  == round(np.float(lj_energies), 5), (
        "Energies do not match!")
print("All tests passed!")
