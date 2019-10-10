import sys
import ase
from ase import io
from ase.calculators.emt import EMT
from ase.calculators.lj import LennardJones
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from collections import defaultdict
from ase.visualize import view
from amp_pytorch.lj_model import lj_optim
from ase.io.trajectory import Trajectory


images_emt = ase.io.read("../datasets/COCu/COCu.traj", ":")
images0 = []
for i in range(101):
    images0.append(images_emt[i])
images1 = ase.io.read("MD_results/MLMD_LJ_COCu.traj", ":")
images2 = ase.io.read("MD_results/MLMD_COCu.traj", ":")

for image in images1:
    image.set_calculator(EMT())
for image in images2:
    image.set_calculator(EMT())

# forces0 = np.array([np.amax(np.abs(image.get_forces())) for image in images0]).reshape(
    # -1, 1
# )
# print(np.mean(forces0))
# forces1 = np.array([np.amax(np.abs(image.get_forces())) for image in images1]).reshape(
    # -1, 1
# )
# print(np.mean(forces1))
# forces2 = np.array([np.amax(np.abs(image.get_forces())) for image in images2]).reshape(
    # -1, 1
# )
# print(np.mean(forces2))
off_atoms = defaultdict(list)
for idx, image in enumerate(images1):
    forces = np.abs(image.get_forces())
    max_force = np.amax(forces, 1).reshape(-1, 1)
    for index, force in enumerate(max_force):
        if force >= 3:
            off_atoms[idx].append(index)

off_image = images1[3]
del off_image[[atom.index for atom in off_image if atom.symbol == "O"]]
del off_image[[atom.index for atom in off_image if ((atom.symbol == "Cu") and
    (atom.index != 22))]]
images = []
energies = []
forces = []
displacements = np.linspace(0.5, 4, 100).reshape(-1, 1)
distances = []
for displacement in displacements:
    new_image = off_image.copy()
    new_image[1].position = new_image[0].position
    new_image[1].position[2] += displacement
    vec = new_image[1].position - new_image[0].position
    d = np.sqrt((vec**2).sum())
    distances.append(d)
    new_image.set_calculator(EMT())
    energy = new_image.get_potential_energy()
    force = np.amax(np.abs(new_image.get_forces()))
    forces.append(force)
    energies.append(energy)
    images.append(new_image)


marks = []
for idx, val in enumerate(forces):
    if abs(val-3.5) <= 0.1:
        marks.append(idx)

fig, ax1 = plt.subplots(figsize=(14.15, 10))
color = 'b'
ax1.set_xlabel(r'Distance, $\AA$', fontsize=28)
ax1.set_ylabel('EMT Energy, eV', color=color, fontsize=28)
ax1.plot(distances, energies, 'b', linewidth=3)
# ax1.set_yscale('log')
ax1.set_ylim(0,50)
ax1.tick_params(axis='y', labelcolor=color, labelsize=28)

ax2 = ax1.twinx()
color = 'g'
ax2.set_ylabel(r'Forces, eV/$\AA$', color=color, fontsize=28)
# ax2.set_yscale('log')
ax2.plot(distances, forces, color=color, linewidth=3)
ax2.set_ylim(0, 150)
# ax2.plot(1.45, 1.1, 'ro', ms=18)
# ax2.plot(2.08, 3.7, 'ro', ms=18, label='_nolegend_')
# ax2.plot(0.67, 120, 'mo', ms=18)
ax2.tick_params(axis='y', labelcolor=color, labelsize=28)
ax1.tick_params(axis='x', labelsize=28)
ax2.tick_params(axis='x', labelsize=28)
fig.legend(['Energy', 'Max |F|', r'5 eV/$\AA$ (ML-LJ)', r'135 eV/$\AA$ (ML)'],
        loc=(.63, .74),
        fontsize=25)
# fig.legend(['Energy', 'Max |F|', r'5 eV/$\AA$ (ML-LJ)'], loc=(.63, .80), fontsize=25)
# fig.legend(['Energy', 'Max |F|'], loc=(.71, .86), fontsize=25)
fig.tight_layout()
plt.show()
# fig.savefig('demo/2atom_ef5135.png', dpi=300)
