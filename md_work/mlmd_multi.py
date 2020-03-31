import os
import ase
import torch
import sklearn
from torch.nn import MSELoss
from skorch import NeuralNetRegressor
from skorch.dataset import CVSplit
from skorch.callbacks import Checkpoint, EpochScoring
from skorch.callbacks.lr_scheduler import LRScheduler
from amptorch.gaussian import SNN_Gaussian
from amp.descriptor.gaussian import Gaussian
from amptorch.model import FullNN, CustomMSELoss, MAELoss
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
from amptorch.lj_model import lj_optim
from amptorch.lj_12_6 import lj_optim as lj_optim12_6
from torch.utils.data import DataLoader
from torch.nn import init
from skorch.utils import to_numpy
import numpy as np
from ase import Atoms, units
from ase.calculators.emt import EMT
import skorch.callbacks.base


class train_end_load_best_valid_loss(skorch.callbacks.base.Callback):
    def on_train_end(self, net, X, y):
        net.load_params("./results/checkpoints/forces_best_params.pt")


def mlmd_run(
    images, filename, dir, count, temp, Gs, lj, forcesonly, scaling, dt, ensemble
):
    if not os.path.exists(dir):
        os.mkdir(dir)

    # lj optimization
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
        cores=10,
        lj_data=lj_data,
        scaling=scaling,
    )
    unique_atoms = training_data.elements
    fp_length = training_data.fp_length
    device = "cpu"

    LR_schedule = LRScheduler("CosineAnnealingLR", T_max=5)
    cp = Checkpoint(
        monitor="forces_score_best", fn_prefix="./results/checkpoints/forces_best_"
    )
    load_best_valid_loss = train_end_load_best_valid_loss()
    torch.set_num_threads(1)

    net = NeuralNetRegressor(
        module=FullNN(
            unique_atoms, [fp_length, 3, 20], device, forcetraining=forcetraining
        ),
        criterion=CustomMSELoss,
        criterion__force_coefficient=0.04,
        optimizer=torch.optim.LBFGS,
        lr=1e-1,
        batch_size=len(training_data),
        # optimizer=torch.optim.Adam,
        # lr=1e-2,
        # batch_size=32,
        # max_epochs=1000,
        iterator_train__collate_fn=collate_amp,
        iterator_train__shuffle=True,
        iterator_valid__collate_fn=collate_amp,
        iterator_valid__shuffle=False,
        device=device,
        train_split=CVSplit(cv=5, random_state=1),
        # train_split=CVSplit(0.2, random_state=1),
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

    print(filename)
    calc = AMP(training_data, net, label=filename)
    calc.train(overwrite=True)

    md_label = dir + filename
    md_run(
        calc=calc,
        starting_image=images[0].copy(),
        temp=temp,
        count=count,
        label=md_label,
        dt=dt,
        ensemble=ensemble,
        rng=False,
    )


def multiple_runs(
    images,
    filename,
    dir,
    ml_run,
    num_iters,
    count,
    temp,
    Gs,
    forcesonly,
    ensemble,
    scaling,
    dt,
):
    for i in range(num_iters):
        if ml_run:
            # ML-MD
            ml_name = filename + "-%s" % str(i + 1)
            mlmd_run(
                images,
                ml_name,
                dir=dir,
                count=count,
                temp=temp,
                Gs=Gs,
                lj=False,
                forcesonly=None,
                scaling=scaling,
                dt=dt,
                ensemble=ensemble,
            )
        # ML-LJ MD
        lj_name = filename + "-LJ-%s" % str(i + 1)
        mlmd_run(
            images,
            lj_name,
            dir=dir,
            count=count,
            temp=temp,
            Gs=Gs,
            lj=True,
            forcesonly=forcesonly,
            scaling=scaling,
            dt=dt,
            ensemble=ensemble,
        )


# Define Training data
# 2ps run of data - 1fs steps
ber_images = ase.io.read("../datasets/COCu_ber_100ps_300K.traj", ":2000:10")

# define symmetry functions to be used
Gs = {}
Gs["G2_etas"] = np.logspace(np.log10(0.05), np.log10(5.0), num=4)
Gs["G2_rs_s"] = [0] * 4
Gs["G4_etas"] = [0.005]
Gs["G4_zetas"] = [1.0, 4.0]
Gs["G4_gammas"] = [+1.0, -1]
Gs["cutoff"] = 5.876798323827276  # EMT asap_cutoff: False

multiple_runs(
    images=ber_images,
    filename="multi_test",
    dir="./ber_results/",
    ml_run=True,
    num_iters=2,
    count=2000,
    temp=300,
    Gs=Gs,
    forcesonly=True,
    ensemble="nvtberendsen",
    scaling="minmax",
    dt=1,
)
