import os
import ase
import torch
import random
from torch.nn import MSELoss
from skorch import NeuralNetRegressor
from skorch.dataset import CVSplit
from skorch.callbacks import Checkpoint, EpochScoring
from skorch.callbacks.lr_scheduler import LRScheduler
from amptorch.gaussian import SNN_Gaussian
from amp.descriptor.gaussian import Gaussian
from amptorch.model import FullNN, CustomLoss, MAELoss
from amptorch.utils import Logger
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
)
from amptorch.skorch_model import AMP
from amptorch.skorch_model.utils import target_extractor, energy_score, forces_score
from amptorch.morse import lj_optim as lj_optim_morse
from torch.utils.data import DataLoader
from torch.nn import init
from skorch.utils import to_numpy
import numpy as np
from ase import Atoms, units
from ase.calculators.emt import EMT
import skorch.callbacks.base


def make_figures(calc, images, filename):
    import seaborn as sns

    os.makedirs("./results/plots/parity", exist_ok=True)

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


def mlmd_run(
    images,
    filename,
    file_dir,
    count,
    temp,
    Gs,
    lj,
    forcesonly,
    scaling,
    dt,
    ensemble,
    rng,
    save_plots,
):
    os.makedirs("./results/checkpoints", exist_ok=True)

    class train_end_load_best_valid_loss(skorch.callbacks.base.Callback):
        def on_train_end(self, net, X, y):
            net.load_params("./results/checkpoints/{}_params.pt".format(filename))

    if not os.path.exists(file_dir):
        os.mkdir(file_dir)
    # lj optimization
    cutoff = Gs["cutoff"]
    lj_data = None
    if lj:
        a = 12
        p0 = [0.972, 6.379, 0.477, 1.09, 8.575, 0.603, 2.168, 3.8386, 1.696]
        params_dict = {"C": [], "O": [], "Cu": []}
        lj_model = lj_optim_morse(
            images, p0, params_dict, cutoff, filename, combo="mean"
        )
        fitted_params = p0
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
        cores=10,
        lj_data=lj_data,
        scaling=scaling,
    )
    unique_atoms = training_data.elements
    fp_length = training_data.fp_length
    device = "cpu"

    LR_schedule = LRScheduler("CosineAnnealingLR", T_max=5)
    cp = Checkpoint(
        monitor="forces_score_best",
        fn_prefix="./results/checkpoints/{}_".format(basename),
    )
    load_best_valid_loss = train_end_load_best_valid_loss()
    torch.set_num_threads(1)

    net = NeuralNetRegressor(
        module=FullNN(
            unique_atoms, [fp_length, 4, 30], device, forcetraining=forcetraining
        ),
        criterion=CustomLoss,
        criterion__force_coefficient=0.04,
        optimizer=torch.optim.LBFGS,
        optimizer__line_search_fn="strong_wolfe",
        lr=1e-1,
        # optimizer=torch.optim.AdamW,
        # lr=1e-3,
        # batch_size=int(len(training_data)/8),
        batch_size=len(training_data),
        max_epochs=300,
        iterator_train__collate_fn=collate_amp,
        iterator_train__shuffle=False,
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
            # LR_schedule,
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
        rng=rng,
    )


def sampler(
    images,
    sample_images,
    count,
    filename,
    file_dir,
    num_samples,
    iteration,
    temp,
    forceonly,
    scaling,
    dt,
    lj,
    Gs,
    ensemble,
    rng,
    save_plots,
):
    random.seed(3)
    sample_points = random.sample(range(1, len(sample_images)), num_samples)
    for idx in sample_points:
        sample_images[idx].set_calculator(EMT())
        images.append(sample_images[idx])
    name = filename + "_{}_iter_{}".format(num_samples, iteration + 1)
    if lj:
        name = filename + "_LJ_{}_iter_{}".format(num_samples, iteration + 1)
    mlmd_run(
        images=images,
        filename=name,
        count=count,
        file_dir=file_dir,
        temp=temp,
        Gs=Gs,
        lj=lj,
        forcesonly=forceonly,
        scaling=scaling,
        dt=dt,
        ensemble=ensemble,
        rng=rng,
        save_plots=save_plots,
    )
    return images


def iterative_sampler(
    images,
    count,
    sample,
    filename,
    file_dir,
    iteration,
    lj,
    temp,
    Gs,
    forceonly,
    scaling,
    dt,
    ensemble,
    rng,
    save_plots,
):
    if lj:
        name = filename + "_LJ"
    else:
        name = filename
    mlmd_run(
        images=images,
        filename=name,
        count=count,
        file_dir=file_dir,
        temp=temp,
        Gs=Gs,
        lj=lj,
        forcesonly=forceonly,
        scaling=scaling,
        dt=dt,
        ensemble=ensemble,
        rng=rng,
        save_plots=save_plots,
    )
    sample_images = ase.io.read(file_dir + name + ".traj", ":")
    for i in range(iteration):
        images = sampler(
            images=images,
            sample_images=sample_images,
            count=count,
            filename=filename,
            file_dir=file_dir,
            num_samples=sample,
            iteration=i,
            temp=temp,
            Gs=Gs,
            lj=lj,
            forceonly=forceonly,
            scaling=scaling,
            dt=dt,
            rng=rng,
            ensemble=ensemble,
            save_plots=save_plots,
        )
        if lj:
            sample_images = ase.io.read(
                file_dir + filename + "_LJ_{}_iter_{}.traj".format(sample, i + 1), ":"
            )
        else:
            sample_images = ase.io.read(
                file_dir + filename + "_{}_iter_{}.traj".format(sample, i + 1), ":"
            )


def generate_data(count, filename, temp, hook, dt, ensemble="NVE"):
    """Generates test or training data with a simple MD simulation."""
    slab = fcc100("Cu", size=(3, 3, 3))
    ads = molecule("CO")
    add_adsorbate(slab, ads, 4, offset=(1, 1))
    cons = FixAtoms(
        indices=[atom.index for atom in slab if (atom.tag == 2 or atom.tag == 3)]
    )
    slab.set_constraint(cons)
    slab.center(vacuum=13.0, axis=2)
    slab.set_pbc(True)
    slab.wrap(pbc=[True] * 3)
    slab.set_calculator(EMT())
    slab.get_forces(apply_constraint=False)
    np.random.seed(1)
    MaxwellBoltzmannDistribution(slab, temp * units.kB)
    if ensemble == "NVE":
        dyn = VelocityVerlet(slab, dt=dt * units.fs)
    elif ensemble == "nvtberendsen":
        dyn = nvtberendsen.NVTBerendsen(slab, dt * units.fs, temp, taut=300 * units.fs)
    elif ensemble == "langevin":
        dyn = Langevin(slab, dt * units.fs, temp * units.kB, 0.002)
    traj = ase.io.Trajectory(filename, "w", slab)
    dyn.attach(traj.write, interval=1)

    def printenergy(a=slab):
        epot = a.get_potential_energy()
        ekin = a.get_kinetic_energy()
        print(
            "Energy per atom: Epot = %.3feV Ekin = %.3feV (T=%3.0fK) "
            "Etot = %.3feV" % (epot, ekin, ekin / (1.5 * units.kB), epot + ekin)
        )

    dyn.attach(printenergy, interval=10)
    dyn.run(count - 1)


def main():
    Gs = {}
    Gs["G2_etas"] = np.logspace(np.log10(0.05), np.log10(5.0), num=4)
    Gs["G2_rs_s"] = [0] * 4
    Gs["G4_etas"] = [0.005]
    Gs["G4_zetas"] = [1.0, 4.0]
    Gs["G4_gammas"] = [+1.0, -1]
    Gs["cutoff"] = 5.876798323827276  # EMT asap_cutoff: False

    num_samples = 100

    basename = "COCu_ml_mdseed_200_5k_2"
    iterative_sampler(
        images=ber_images,
        count=5000,
        sample=num_samples,
        filename=basename,
        file_dir="./morse/",
        iteration=3,
        # lj=True,
        lj=False,
        temp=300,
        Gs=Gs,
        forceonly=True,
        ensemble="nvtberendsen",
        scaling=None,
        # scaling="rel",
        dt=1,
        rng=True,
        save_plots=True,
    )


if __name__ == "__main__":
    main()
