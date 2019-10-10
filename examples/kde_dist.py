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


def kde_plots(data0, data1, data2, data3=None, data4=None, data5=None,
        filename=None):
    fig, ax = plt.subplots(figsize=(14.15, 10))
    sns.distplot(
        data0,
        hist=False,
        kde=True,
        label='EMT',
        color='b',
        hist_kws={"edgecolor": "black", "alpha": 0.1},
        kde_kws={"linewidth": 2},
        ax=ax,
    )
    if data1 is not None:
        sns.distplot(
            data1,
            hist=False,
            label="iter 1",
            kde=True,
            color='darkorange',
            hist_kws={"edgecolor": "black", "alpha": 0.1},
            kde_kws={"linewidth": 2},
        )
    if data2 is not None:
        sns.distplot(
            data2,
            hist=False,
            kde=True,
            color='r',
            label='iter 2',
            hist_kws={"edgecolor": "black", "alpha": 0.1},
            kde_kws={"linewidth": 2, "linestyle": "-"},
        )
    if data3 is not None:
        sns.distplot(
            data3,
            hist=False,
            kde=True,
            label='iter 3',
            color='g',
            hist_kws={"edgecolor": "black", "alpha": 0.1},
            kde_kws={"linewidth": 2, "linestyle": "-"},
        )
    if data4 is not None:
        sns.distplot(
            data4,
            hist=False,
            kde=True,
            label='iter 4',
            color='y',
            hist_kws={"edgecolor": "black", "alpha": 0.1},
            kde_kws={"linewidth": 2, "linestyle": "-"},
        )
    if data5 is not None:
        sns.distplot(
            data5,
            hist=False,
            kde=True,
            label='iter 5',
            color='k',
            hist_kws={"edgecolor": "black", "alpha": 0.1},
            kde_kws={"linewidth": 2, "linestyle": "-"},
        )
   # plt.legend(loc="upper left", fontsize=28)
    plt.legend(loc="best", fontsize=28)
    # plt.xlabel("log max|F|", fontsize=28)
    plt.xlabel("log max|F|", fontsize=32)
    plt.ylabel("Density Function", fontsize=32)
    # plt.title("Energy Dist.", fontsize=30)
    plt.title("Force Distributions", fontsize=35)
    # plt.title("Force Dist. - EMT evaluated", fontsize=30)
    # plt.title("Force Dist.", fontsize=30)
    plt.xlim(left=-5, right=4)
    # plt.xlim(left=0, right=100)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    # plt.ylim(0, 3)
    # plt.annotate('EMT (target)', (-2, .75), fontsize=30, color='steelblue')
    # plt.annotate('EMT (target)', (0.3, 1.75), fontsize=30, color='steelblue')
    # plt.annotate('ML-apparent', (-2.1, .75), fontsize=30, color='r')
    # plt.annotate('ML-actual', (2, .65), fontsize=30, color='g')
    # plt.annotate('ML-LJ-apparent', (-2.1, 2.5), fontsize=30, color='r')
    # plt.annotate('ML-LJ', (13, .35), fontsize=30, color='darkorange')
    # plt.annotate('ML-LJ-actual', (.8, 1.), fontsize=30, color='darkorange')
    # plt.savefig("MD_results/plots/" + filename + "_kde.png", dpi=300)
    plt.show()

def energy_plots(data0, data1):
    fig, ax = plt.subplots(figsize=(14.15, 10))
    num_images = range(len(data0))
    plt.plot(num_images, data0, label='emt')
    plt.plot(num_images, data1, label='ml')
    plt.legend(fontsize=20)
    plt.xlabel('image', size=20)
    plt.ylabel('Energy, eV', size=20)
    plt.show()

# images_emt = ase.io.read("../datasets/COCu/COCu_pbc.traj", ":")
images_emt = ase.io.read("../datasets/COCu/COCu_pbc_300K.traj", ":")
images0 = []
for i in range(101):
    images0.append(images_emt[i])

images1 = \
ase.io.read("MD_results/COCu/pbc_300K/logcosh/paper/MLMD_COCu_pbc_300K_logcosh_2_5_iter_1.traj", ":")

images2 = \
ase.io.read("MD_results/COCu/pbc_300K/logcosh/paper/MLMD_COCu_pbc_300K_logcosh_2_5_iter_2.traj", ":")

images3 = \
ase.io.read("MD_results/COCu/pbc_300K/logcosh/paper/MLMD_COCu_pbc_300K_logcosh_2_5_iter_3.traj", ":")

images4 = \
ase.io.read("MD_results/COCu/pbc_300K/logcosh/paper/MLMD_COCu_pbc_300K_logcosh_2_5_iter_4.traj", ":")

images5 = \
ase.io.read("MD_results/COCu/pbc_300K/logcosh/paper/MLMD_COCu_pbc_300K_logcosh_2_5_iter_5.traj", ":")


for idx, image in enumerate(images1):
    image.set_calculator(EMT())
    images2[idx].set_calculator(EMT())
    images3[idx].set_calculator(EMT())
    images4[idx].set_calculator(EMT())
    images5[idx].set_calculator(EMT())
# images2 = images2[1:]
# images0 = images0[1:]

forces0 = np.array([np.amax(np.abs(image.get_forces())) for image in images0]).reshape(
    -1, 1
)+1e-10
forces1 = np.array([np.amax(np.abs(image.get_forces())) for image in images1]).reshape(
    -1, 1
)+1e-10
forces2 = np.array([np.amax(np.abs(image.get_forces())) for image in images2]).reshape(
    -1, 1
)+1e-10
forces3 = np.array([np.amax(np.abs(image.get_forces())) for image in images3]).reshape(
    -1, 1
)+1e-10
forces4 = np.array([np.amax(np.abs(image.get_forces())) for image in images4]).reshape(
    -1, 1
)+1e-10
forces5 = np.array([np.amax(np.abs(image.get_forces())) for image in images5]).reshape(
    -1, 1
)+1e-10
# energy0 = np.array([image.get_potential_energy() for image in images0])
# energy1 = np.array([image.get_potential_energy() for image in images1])
# energy2 = np.array([image.get_potential_energy() for image in images2])

# energy_plots(energy0, energy2)
kde_plots(data0=np.log10(forces0), data1=np.log10(forces1), data2=np.log10(forces2),
        data3=np.log10(forces3), data4=np.log10(forces4), data5=np.log10(forces5), filename="test")
# kde_plots(np.log10(forces0), None, None,
        # "COCu_emt")
# kde_plots(forces0, forces1, forces2, "COCu_ALL_emt_logF")
