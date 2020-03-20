import multiprocessing as mp
import ase
import time
import torch
from torch.nn import MSELoss
from skorch import NeuralNetRegressor
from skorch.dataset import CVSplit
from skorch.callbacks import Checkpoint, EpochScoring
from skorch.callbacks.lr_scheduler import LRScheduler
from amptorch.gaussian import SNN_Gaussian
from amp.descriptor.gaussian import Gaussian
from amptorch.model import FullNN, CustomLoss, MAELoss
from amptorch.data_preprocess import (
    AtomsDataset,
    factorize_data,
    collate_amp,
    TestDataset,
)
from md_work.md_utils import (
    md_run,
    calculate_energies,
    calculate_forces,
    time_plots,
    kde_plots,
    calculate_rmse,
)
from amptorch.skorch_model import AMP
from amptorch.skorch_model.utils import target_extractor, energy_score, forces_score
from amptorch.lj_model import lj_optim
from amptorch.utils import make_amp_descriptors_simple_nn
from torch.utils.data import DataLoader
from torch.nn import init
from skorch.utils import to_numpy
import numpy as np
import skorch.callbacks.base
import matplotlib.pyplot as plt

def make_figures(calc, images, filename):
    import seaborn as sns

    actual_energies = [
        image.get_potential_energy(apply_constraint=False) for image in images
    ]
    pred_energies = [calc.get_potential_energy(image) for image in images]
    lims = [8, 10]
    energy_fig = sns.jointplot(
        pred_energies, actual_energies, kind="hex", bins="log", extent=lims + lims
    )
    ax1 = energy_fig.ax_joint
    _ = ax1.set_xlim(lims)
    _ = ax1.set_ylim(lims)
    _ = ax1.plot(lims, lims, "--")
    _ = ax1.set_xlabel("ML Energy, eV")
    _ = ax1.set_ylabel("Target Energy, eV")
    energy_fig.savefig("./results/plots/parity/energy_parity_plot_{}".format(filename))

    actual_forces = np.concatenate(
        np.array([image.get_forces(apply_constraint=False) for image in images])
    )
    pred_forces = np.concatenate(np.array([calc.get_forces(image) for image in images]))

    lims = [-5, 5]
    force_fig = sns.jointplot(
        pred_forces, actual_forces, kind="hex", bins="log", extent=lims + lims
    )
    ax2 = force_fig.ax_joint
    _ = ax2.set_xlim(lims)
    _ = ax2.set_ylim(lims)
    _ = ax2.plot(lims, lims, "--")
    _ = ax2.set_xlabel("ML Forces, eV/A")
    _ = ax2.set_ylabel("Target Forces, eV/A")
    force_fig.savefig("./results/plots/parity/forces_parity_plot_{}".format(filename))


def mlmd_run(inputs):
    images, filename, file_dir, count, lj, temp, Gs, forcesonly, scaling, dt, ensemble, save_plots = (
        inputs
    )

    class train_end_load_best_valid_loss(skorch.callbacks.base.Callback):
        def on_train_end(self, net, X, y):
            net.load_params(
                "./results/checkpoints/forces_{}_params.pt".format(filename)
            )

    cutoff = Gs["cutoff"]
    lj_data = None
    if lj:
        a = 12
        p0 = [1.0, 6.3535, 0, 1.0808, 8.5357, 0, 2.1717, 3.7575, 0, a]
        params_dict = {"C": [], "O": [], "Cu": []}
        lj_model = lj_optim(
            images, p0, params_dict, cutoff, filename, forcesonly=forcesonly
        )
        fitted_params = lj_model.fit()
        lj_energies, lj_forces, num_atoms = lj_model.lj_pred(
            images, fitted_params, params_dict
        )
        lj_data = [
            lj_energies,
            lj_forces,
            num_atoms,
            fitted_params,
            params_dict,
            lj_model,
        ]

    forcetraining = True
    training_data = AtomsDataset(
        images,
        SNN_Gaussian,
        Gs,
        forcetraining=forcetraining,
        label=filename,
        cores=1,
        lj_data=lj_data,
        scaling=scaling,
    )
    unique_atoms = training_data.elements
    fp_length = training_data.fp_length
    device = "cpu"

    cp = Checkpoint(
        monitor="forces_score_best",
        fn_prefix="./results/checkpoints/forces_{}_".format(filename),
    )
    load_best_valid_loss = train_end_load_best_valid_loss()
    torch.set_num_threads(1)

    net = NeuralNetRegressor(
        module=FullNN(
            unique_atoms, [fp_length, 3, 20], device, forcetraining=forcetraining
        ),
        criterion=CustomLoss,
        criterion__force_coefficient=0.04,
        optimizer=torch.optim.LBFGS,
        lr=1e-1,
        batch_size=len(training_data),
        max_epochs=200,
        iterator_train__collate_fn=collate_amp,
        iterator_train__shuffle=True,
        iterator_valid__collate_fn=collate_amp,
        iterator_valid__shuffle=False,
        device=device,
        train_split=CVSplit(cv=5, random_state=1),
        callbacks=[
            EpochScoring(
                forces_score,
                on_train=False,
                use_caching=True,
                target_extractor=target_extractor,
            ),
            EpochScoring(
                energy_score,
                on_train=False,
                use_caching=True,
                target_extractor=target_extractor,
            ),
            cp,
            load_best_valid_loss,
        ],
    )
    calc = AMP(training_data, net, label=filename)
    calc.train(overwrite=True)
    if save_plots:
        make_figures(calc, images, filename)
    md_label = file_dir + filename
    md_run(
        calc=calc,
        starting_image=images[0].copy(),
        temp=temp,
        count=count,
        label=md_label,
        dt=dt,
        ensemble=ensemble,
    )
    num_atoms = 29
    md_traj = ase.io.read(md_label + ".traj", ":")
    apparent_energies, actual_energies = calculate_energies(md_traj)
    apparent_forces, actual_forces = calculate_forces(md_traj, type=None)
    e_rmse = calculate_rmse(
        apparent_energies, actual_energies, num_atoms, dtype="energy"
    )
    f_rmse = calculate_rmse(apparent_forces, actual_forces, num_atoms, dtype="forces")
    return [e_rmse, f_rmse]


# Define Training data
ber_images = ase.io.read("../datasets/COCu_ber_100ps_300K.traj", ":5500")

# define symmetry functions to be used
Gs = {}
Gs["G2_etas"] = np.logspace(np.log10(0.05), np.log10(5.0), num=4)
Gs["G2_rs_s"] = [0] * 4
Gs["G4_etas"] = [0.005]
Gs["G4_zetas"] = [1.0, 4.0]
Gs["G4_gammas"] = [+1.0, -1]
Gs["cutoff"] = 5.876798323827276  # EMT asap_cutoff: False

# AtomsDataset(ber_images, SNN_Gaussian, Gs, forcetraining=True,
        # label='learning_curve', cores=10, lj_data=None)

pool = mp.Pool(processes=5)
# data_size = range(100, 5100, 100)
data_size = range(1000, 2000, 100)
filename = "ml_lcurve"
count = 2000
lj = False
temp = 300
forcesonly = False
scaling = "minmax"
dt = 1
ensemble = "nvtberendsen"
save_plot = True

inputs = []
for size in data_size:
    inputs.append(
        (
            ber_images[:size],
            filename + "_{}".format(size),
            "./ber_results/learning_curve/",
            count,
            lj,
            temp,
            Gs,
            forcesonly,
            scaling,
            dt,
            ensemble,
            save_plot,
        )
    )
results = pool.map(mlmd_run, inputs)
e_rmse = [x[0] for x in results]
f_rmse = [x[1] for x in results]
plt.plot(data_size, e_rmse)
plt.xlabel('# of Training Data')
plt.ylabel('Energy RMSE, eV/atom')
plt.savefig('./ber_results/learning_curve/ml_energy_learning_curve.pdf')
plt.plot(data_size, f_rmse)
plt.xlabel('# of Training Data')
plt.ylabel('Force RMSE, eV/Angstrom')
plt.savefig('./ber_results/learning_curve/ml_force_learning_curve.pdf')
