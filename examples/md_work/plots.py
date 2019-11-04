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


def kde_plots(emt, forces, filename=None):
    fig, ax = plt.subplots(figsize=(14.15, 10))
    sns.distplot(
        emt,
        hist=False,
        kde=True,
        label="target",
        color="k",
        hist_kws={"edgecolor": "black", "alpha": 0.1},
        kde_kws={"linewidth": 2},
        ax=ax,
    )
    for i, data in enumerate(forces):
        sns.distplot(
            data,
            hist=False,
            label="iter %s" % (i + 1),
            kde=True,
            hist_kws={"edgecolor": "black", "alpha": 0.1},
            kde_kws={"linewidth": 2},
        )
    plt.legend(loc="best", fontsize=28)
    plt.xlabel("log max|F|", fontsize=32)
    plt.ylabel("Density Function", fontsize=32)
    plt.title("Force Distributions", fontsize=35)
    plt.xlim(left=-5, right=4)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.show()

def rmse_plots(emt, forces, interval, filename=None):
    fig, ax = plt.subplots(figsize=(14.15, 10))
    rmse = []
    resample_images = [5*(i+1) for i in range(len(forces))]
    for file_forces in forces:
        num_atoms = file_forces.shape[1]
        mse = (1/(3*num_atoms))*((emt-file_forces)**2).sum()
        rmse.append(np.sqrt(mse))
    plt.plot(resample_images, rmse)
    plt.xlabel("Images resampled", fontsize=25)
    plt.ylabel("Force RMSE (eV/A)", fontsize=25)
    plt.title("RMSE vs Resampled Images", fontsize=25)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.show()

def compute_forces(dir, files, type='max'):
    forces = []
    for file in files:
        images = ase.io.read(dir+file, ":")
        data = []
        for image in images:
            image.set_calculator(EMT())
            if type=='max':
                data.append(np.log10(np.amax(np.abs(image.get_forces()))))
            else:
                data.append(image.get_forces())
        if type=='max':
            forces.append((np.array(data)))
        else:
            forces.append(np.array(data))
    return forces


    forces = []

images_emt = ase.io.read("../../datasets/COCu/COCu_pbc_300K.traj", ":")
images0 = []
forces0 = []
for i in range(101):
    images0.append(images_emt[i])
    # forces0.append(images0[i].get_forces())
    forces0.append(np.log10((np.array(np.amax(np.abs(images_emt[i].get_forces()))))))
forces0 = np.array(forces0)

files = [
    "MLMD_COCu_pbc_300K_l2amp_6SF-1.traj",
    "MLMD_COCu_pbc_300K_l2amp_6SF-LJ-1.traj",
    "MLMD_COCu_pbc_300K_l2amp_6SF-2.traj",
    "MLMD_COCu_pbc_300K_l2amp_6SF-LJ-2.traj"
]

forces = compute_forces(dir="MD_results/",
        files=files, type='max')
# rmse_plots(forces0, forces, 5, "test")
kde_plots(forces0, forces, "test")
