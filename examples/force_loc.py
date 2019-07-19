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


def force_dist(forces0, forces1, forces2, data):
    plt.hist(forces0, bins=100, alpha=0.5, label="EMT")
    plt.hist(forces1, bins=100, alpha=0.5, label="ML-LJ")
    if forces2 is not None:
        plt.hist(forces2, bins=100, alpha=0.5, label="ML")
    plt.legend(loc="upper right")
    plt.xlabel("Max Force, eV/A")
    plt.ylabel("Frequency")
    plt.savefig("MD_results/plots/" + data + "_force_dist.png")
    plt.show()


def force_gauss(forces0, forces1, forces2, data):
    forces0 = np.sort(forces0, 0)
    forces1 = np.sort(forces1, 0)
    plt.plot(
        forces0, stats.norm.pdf(forces0, np.mean(forces0), np.std(forces0)), label="EMT"
    )
    plt.plot(
        forces1,
        stats.norm.pdf(forces1, np.mean(forces1), np.std(forces1)),
        label="ML-LJ",
    )
    if forces2 is not None:
        forces2 = np.sort(forces2, 0)
        plt.plot(
            forces2,
            stats.norm.pdf(forces2, np.mean(forces2), np.std(forces2)),
            label="ML",
        )
    plt.legend(loc="upper right")
    plt.xlabel("Max Force, eV/A")
    plt.savefig("MD_results/plots/" + data + "_force_gauss.png")
    plt.show()


images_emt = ase.io.read("../datasets/COCu/COCu.traj", ":")
images0 = []
for i in range(101):
    images0.append(images_emt[i])
images1 = ase.io.read("MD_results/MLMD_LJ_COCu.traj", ":")
# images2 = ase.io.read("MD_results/MLMD_COCu.traj", ":")

for image in images1:
    image.set_calculator(EMT())
# for image in images2:
# image.set_calculator(EMT())

forces0 = np.array([np.amax(np.abs(image.get_forces())) for image in images0]).reshape(
    -1, 1
)
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
    view(new_image)
    sys.exit()
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

fig, ax1 = plt.subplots()
color = 'tab:blue'
ax1.set_xlabel('Distance, A')
ax1.set_ylabel('EMT Energy, eV', color=color)
ax1.plot(distances, energies, 'b--')
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()
color = 'tab:green'
ax2.set_ylabel('Forces, eV/A', color=color)
ax2.plot(distances, forces, color=color)
ax2.axhline(y=5, color='r')
ax2.axhline(y=100, color='purple')
ax2.tick_params(axis='y', labelcolor=color)
fig.legend(['Energy', 'Max |F|', '5 eV/A', '100 eV/A'], loc='upper right')
fig.tight_layout()
plt.show()
# fig.savefig('MD_results/plots/lj_model_error.png')
