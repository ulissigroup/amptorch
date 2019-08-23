import sys
import ase
from ase import io
from ase.calculators.emt import EMT
from ase import units
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.stats as stats
from matplotlib import rc
import seaborn as sns
from ase.visualize import view


def kde_plots(data0, data1, data2, filename):
    fig, ax = plt.subplots(figsize=(14.15, 10))
    sns.distplot(
        data0,
        hist=False,
        kde=True,
        # label='EMT',
        hist_kws={"edgecolor": "black", "alpha": 0.1},
        kde_kws={"linewidth": 5},
        ax=ax,
    )
    if data1 is not None:
        sns.distplot(
            data1,
            hist=False,
            # label="ML",
            kde=True,
            # color='darkorange',
            color='g',
            hist_kws={"edgecolor": "black", "alpha": 0.1},
            kde_kws={"linewidth": 5},
        )
    if data2 is not None:
        sns.distplot(
            data2,
            hist=False,
            kde=True,
            color='r',
            # label='ML-LJ',
            hist_kws={"edgecolor": "black", "alpha": 0.1},
            kde_kws={"linewidth": 5, "linestyle": "-"},
        )
    # plt.legend(loc="upper left", fontsize=28)
    # plt.legend(loc="best", fontsize=28)
    # plt.xlabel("log max|F|", fontsize=28)
    plt.xlabel("log max|F|", fontsize=32)
    plt.ylabel("Density Function", fontsize=32)
    # plt.title("Energy Dist.", fontsize=30)
    plt.title("Force Distributions", fontsize=35)
    # plt.title("Force Dist. - EMT evaluated", fontsize=30)
    # plt.title("Force Dist.", fontsize=30)
    plt.xlim(left=-3, right=4)
    # plt.xlim(left=0, right=100)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.ylim(0, 3)
    # plt.annotate('EMT (target)', (-2, .75), fontsize=30, color='steelblue')
    plt.annotate('EMT (target)', (0.3, 1.75), fontsize=30, color='steelblue')
    # plt.annotate('ML-apparent', (-2.1, .75), fontsize=30, color='r')
    # plt.annotate('ML-actual', (2, .65), fontsize=30, color='g')
    # plt.annotate('ML-LJ-apparent', (-2.1, 2.5), fontsize=30, color='r')
    # plt.annotate('ML-LJ', (13, .35), fontsize=30, color='darkorange')
    # plt.annotate('ML-LJ-actual', (.8, 1.), fontsize=30, color='darkorange')
    plt.savefig("MD_results/plots/" + filename + "_kde.png", dpi=300)
    plt.show()


# images_emt = ase.io.read("../datasets/COCu/COCu_pbc.traj", ":")
images_emt = ase.io.read("../datasets/COCu/COCu.traj", ":")
images0 = []
for i in range(101):
    images0.append(images_emt[i])
for image in images0:
    print('PE: %.3f, KE: %.3f, tot: %.3f, T: %3.0fK' %
            (image.get_potential_energy(), image.get_kinetic_energy(),
                (image.get_potential_energy()+image.get_kinetic_energy()),
                image.get_kinetic_energy()/(1.5*units.kB)) )
sys.exit()
# images1 = ase.io.read("MD_results/MLMD_LJ_COCu_2.traj", ":")
images1 = ase.io.read("MD_results/MLMD_COCu.traj", ":")
# images1 = ase.io.read("MD_results/MLMD_LJ_COCu_pbc2.traj", ":")
# images2 = ase.io.read("MD_results/MLMD_COCu_new.traj", ":")
images2 = ase.io.read("MD_results/MLMD_COCu_pbc.traj", ":")
# images2 = ase.io.read("MD_results/MLMD_COCu.traj", ":")
images2_2 = ase.io.read("MD_results/MLMD_COCu.traj", ":")


# j = images1[50]
# j.rotate(-90, 'x', rotate_cell=True)
# j.rotate(30, 'y', rotate_cell=True)
# j.rotate(20, 'x', rotate_cell=True)
# j = j.repeat((3, 3, 1))
# j.wrap()
# ase.io.write('emtmllj.png',j)
# sys.exit()


for image in images1:
    image.set_calculator(EMT())
# for image in images2:
    # image.set_calculator(EMT())

forces0 = np.array([np.amax(np.abs(image.get_forces())) for image in images0]).reshape(
    -1, 1
)
energy0 = np.array([image.get_potential_energy() for image in images0]).reshape(-1, 1)

forces1 = np.array([np.amax(np.abs(image.get_forces())) for image in images1]).reshape(
    -1, 1
)
energy1 = np.array([image.get_potential_energy() for image in images1]).reshape(-1, 1)

forces2 = np.array([np.amax(np.abs(image.get_forces())) for image in images2]).reshape(
    -1, 1
)
energy2 = np.array([image.get_potential_energy() for image in images2]).reshape(-1, 1)
forces2_2 = np.array([np.amax(np.abs(image.get_forces())) for image in
    images2_2]).reshape(
    -1, 1
)

# kde_plots(forces0, forces1, None, "COCu_ML_emt_F")
# kde_plots(np.log10(forces0), np.log10(forces2), np.log10(forces1),
        # "COCu_ALL2_emt_logF")
kde_plots(np.log10(forces0), None, None,
        "COCu_emt")
# kde_plots(forces0, forces1, forces2_2, "COCu_ALL_emt_logF")
# kde_plots(energy0, energy2, energy1, "COCu_ALL2_emt_energy")
# kde_plots(np.log10(energy0), np.log10(energy1), None, "COCu_LJ_emt_E")
# kde_plots(np.log10(energy0), np.log10(energy1), np.log10(energy2), "COCu_ALL_emt_E")
