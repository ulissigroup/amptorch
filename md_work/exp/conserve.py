import sys
import numpy as np
import matplotlib.pyplot as plt
from ase import io
from ase.calculators.emt import EMT
from amp import Amp

emt_images = io.read("../datasets/COCu_nve_300K.traj", ":400:4")
emt_pot_energies = np.array([image.get_potential_energy(apply_constraint=False)
    for image in emt_images])
emt_ke_energies = np.array([image.get_kinetic_energy() for image in
    emt_images])
emt_tot_energy = (emt_pot_energies + emt_ke_energies).reshape(-1, 1)
print(len(emt_tot_energy))

ml_images = io.read("./nve_results/COCu_nve_ampG_300K-1.traj", ":1000:10")
# ml_images = io.read("../amp_benchmark/amp_lang.traj", ":5000:10")
# calc = Amp.load('../amp_benchmark/calc-untrained-parameters.amp')
# for image in ml_images:
    # image.set_calculator(calc)

ml_pot_energies = np.array([image.get_potential_energy(apply_constraint=False)
    for image in ml_images]).reshape(-1, 1)
ml_ke_energies = np.array([image.get_kinetic_energy() for image in
    ml_images]).reshape(-1, 1)
ml_tot_energy = (ml_pot_energies + ml_ke_energies).reshape(-1, 1)

# mllj_images = io.read("./COCu_LJ_nve_0_1dt_300K.traj", ":5000:5")
# mllj_images = io.read("./COCu_LJ_nve_300K.traj", ":2000:4")
# mllj_pot_energies = np.array([image.get_potential_energy(apply_constraint=False)
    # for image in mllj_images])
# mllj_ke_energies = np.array([image.get_kinetic_energy() for image in
    # mllj_images])
# mllj_tot_energy = (mllj_pot_energies + mllj_ke_energies).reshape(-1, 1)


time = range(100)
plt.plot(time, emt_tot_energy, '-g', label='emt total E')
# plt.plot(time, emt_pot_energies, '-g', label='emt PE')
# plt.plot(time, emt_ke_energies, '--g', label='emt KE', markersize=1)
plt.plot(time, ml_tot_energy, '-b', label='ml total E')
# plt.plot(time, ml_pot_energies, '-b', label='ml PE')
# plt.plot(time, ml_ke_energies, '--b', label='ml KE', markersize=1)
# plt.plot(time, mllj_tot_energy, '-r', label='mllj total E')
# plt.plot(time, mllj_pot_energies, '-r', label='mllj PE')
# plt.plot(time, mllj_ke_energies, '--r', label='mllj KE', markersize=1)

plt.legend()
plt.xlabel('time, fs')
plt.ylabel('energy, eV')
plt.show()
