import sys
import ase
from ase import io
from ase.calculators.emt import EMT
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import seaborn as sns


def kde_plots(data0, data1, data2, filename):
    sns.distplot(
        data0,
        hist=True,
        bins=100,
        kde=True,
        label="EMT",
        hist_kws={"edgecolor": "black", "alpha": 0.1},
        kde_kws={"linewidth": 3},
    )
    sns.distplot(
        data1,
        hist=True,
        bins=100,
        kde=True,
        label="ML-LJ",
        hist_kws={"edgecolor": "black", "alpha": 0.1},
        kde_kws={"linewidth": 3},
    )
    if data2 is not None:
        sns.distplot(
            data2,
            hist=True,
            bins=100,
            kde=True,
            label="ML",
            hist_kws={"edgecolor": "black", "alpha": 0.1},
            kde_kws={"linewidth": 3},
        )
    plt.legend(loc="upper right")
    plt.xlabel("log max|F|")
    plt.ylabel("Density")
    plt.title("Force Dist. - EMT/EMT@ML/EMT@ML-LJ")
    plt.savefig("MD_results/plots/" + filename + "_kde.png", dpi=300)
    plt.show()
    sys.exit()


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

# kde_plots(forces0, forces1, None, "COCu_LJ_forces")
# kde_plots(np.log10(forces0), np.log10(forces1), None, "COCu_LJ_emt_logF")
# kde_plots(np.log10(forces0), np.log10(forces1), np.log10(forces2),
        # "COCu_ALL_emt_logF")
# kde_plots(energy0, energy1, None, "COCu_LJ_energy")
# kde_plots(np.log10(energy0), np.log10(energy1), None, "COCu_LJ_emt_E")
# kde_plots(np.log10(energy0), np.log10(energy1), np.log10(energy2), "COCu_ALL_E")
