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
from amptorch.model import FullNN, CustomLoss
from amptorch.data_preprocess import AtomsDataset, collate_amp
from md_work.md_utils import md_run
from amptorch.skorch_model import AMP
from amptorch.skorch_model.utils import target_extractor, energy_score, forces_score
from amptorch.morse import morse_potential
from torch.utils.data import DataLoader
from torch.nn import init
import numpy as np
from ase import Atoms, units
from ase.calculators.emt import EMT
from ase.calculators.vasp import Vasp
import skorch.callbacks.base
from md_work.delta_ml.dft.md_utils import MDsimulate


def make_figures(calc, images, filename):
    import seaborn as sns

    os.makedirs("./results/plots/parity", exist_ok=True)

    actual_energies = [
        image.get_potential_energy(apply_constraint=False) for image in images
    ]
    pred_energies = [calc.get_potential_energy(image) for image in images]
    lims = [8, 10]
    energy_fig = sns.jointplot(
        actual_energies, pred_energies, kind="hex", bins="log", extent=lims + lims
    )
    ax1 = energy_fig.ax_joint
    _ = ax1.set_xlim(lims)
    _ = ax1.set_ylim(lims)
    _ = ax1.plot(lims, lims, "--")
    _ = ax1.set_ylabel("ML Energy, eV")
    _ = ax1.set_xlabel("Target Energy, eV")
    energy_fig.savefig("./results/plots/parity/energy_parity_plot_{}".format(filename))

    actual_forces = np.concatenate(
        np.array([image.get_forces(apply_constraint=False) for image in images])
    )
    pred_forces = np.concatenate(np.array([calc.get_forces(image) for image in images]))

    lims = [-5, 5]
    force_fig = sns.jointplot(
        actual_forces, pred_forces, kind="hex", bins="log", extent=lims + lims
    )
    ax2 = force_fig.ax_joint
    _ = ax2.set_xlim(lims)
    _ = ax2.set_ylim(lims)
    _ = ax2.plot(lims, lims, "--")
    _ = ax2.set_ylabel("ML Forces, eV/A")
    _ = ax2.set_xlabel("Target Forces, eV/A")
    force_fig.savefig("./results/plots/parity/forces_parity_plot_{}".format(filename))


def trainer(images, filename, file_dir, Gs, morse, save_plots):
    os.makedirs("./results/checkpoints", exist_ok=True)
    os.makedirs(file_dir, exist_ok=True)

    class train_end_load_best_valid_loss(skorch.callbacks.base.Callback):
        def on_train_end(self, net, X, y):
            net.load_params("./results/checkpoints/{}_params.pt".format(filename))

    cutoff = Gs["cutoff"]
    morse_data = None
    if morse:
        params = {
            "C": {"re": 0.972, "D": 6.379, "sig": 0.477},
            "O": {"re": 1.09, "D": 8.575, "sig": 0.603},
            "Cu": {"re": 2.168, "D": 3.8386, "sig": 1.696},
        }
        morse_model = morse_potential(images, params, cutoff, filename, combo="mean")
        morse_energies, morse_forces, num_atoms = morse_model.morse_pred(images, params)
        morse_data = [morse_energies, morse_forces, num_atoms, params, morse_model]

    forcetraining = True
    training_data = AtomsDataset(
        images,
        SNN_Gaussian,
        Gs,
        forcetraining=forcetraining,
        label=filename,
        cores=10,
        delta_data=morse_data,
    )
    unique_atoms = training_data.elements
    fp_length = training_data.fp_length
    device = "cpu"

    cp = Checkpoint(
        monitor="forces_score_best",
        fn_prefix="./results/checkpoints/{}_".format(filename),
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
        optimizer__line_search_fn="strong_wolfe",
        lr=1e-1,
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
        ],
    )

    calc = AMP(training_data, net, label=filename)
    calc.train(overwrite=True)

    if save_plots:
        make_figures(calc, images, filename)
    return calc


def iterative_sampler(
    images,
    utility_function,
    parent_calc,
    samples_to_retrain,
    filename,
    file_dir,
    iterations,
    morse,
    Gs,
    save_plots,
):
    sample_candidates = None
    for iteration in range(iterations + 1):
        if iteration > 0:
            random.seed(3)
            sample_points = random.sample(
                range(1, len(sample_candidates)), samples_to_retrain
            )
            for idx in sample_points:
                sample_candidates[idx].set_calculator(parent_calc)
                images.append(sample_candidates[idx])
        name = filename + "_{}_iter_{}".format(samples_to_retrain, iteration)
        if morse:
            name = filename + "_morse_{}_iter_{}".format(samples_to_retrain, iteration)
        ml_calc = trainer(
            images=images,
            filename=name,
            file_dir=file_dir,
            Gs=Gs,
            morse=morse,
            save_plots=save_plots,
        )

        md_label = file_dir + name
        utility_function.run(calc=ml_calc, filename=md_label)
        sample_candidates = utility_function.get_trajectory(
            filename=md_label, start_count=0, end_count=-1, interval=1
        )


def main(calculator):
    # define function to be called with ML
    training_data = ase.io.read("../../../datasets/COCu_ber_50ps_300K.traj", ":2000:10")
    md_runner = MDsimulate(
        ensemble="nvtberendsen",
        dt=1,
        temp=300,
        count=5000,
        initial_geometry=training_data[0].copy(),
    )

    if calculator == "EMT":
        parent_calculator = EMT()
    elif calculator == "VASP":
        parent_calculator = Vasp(
            prec="Normal",
            algo="Normal",
            xc="PBE",
            gga="RP",
            encut=400,
            lreal=False,
            ediff=1e-4,
            ispin=1,
            NELM=100,
            lwave=False,
            lcharg=False,
            nsw=0,
            kpoints=(4, 4, 1),
        )

    Gs = {}
    Gs["G2_etas"] = np.logspace(np.log10(0.05), np.log10(5.0), num=4)
    Gs["G2_rs_s"] = [0] * 4
    Gs["G4_etas"] = [0.005]
    Gs["G4_zetas"] = [1.0, 4.0]
    Gs["G4_gammas"] = [+1.0, -1]
    Gs["cutoff"] = 5.876798323827276  # EMT asap_cutoff: False

    num_samples_to_retrain = 100

    basename = "COCu_morse_test"
    iterative_sampler(
        images=training_data,
        utility_function=md_runner,
        parent_calc=parent_calculator,
        samples_to_retrain=num_samples_to_retrain,
        filename=basename,
        file_dir="./morse/",
        iterations=3,
        morse=True,
        Gs=Gs,
        save_plots=True,
    )


if __name__ == "__main__":
    main("EMT")
