import multiprocessing
import sys
import os
import ase
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
)
from amptorch.skorch_model import AMP
from amptorch.skorch_model.utils import target_extractor, energy_score, forces_score
from amptorch.lj_model import lj_optim
from torch.utils.data import DataLoader
from torch.nn import init
from skorch.utils import to_numpy
import numpy as np
from ase import Atoms, units
from ase.calculators.emt import EMT
import skorch.callbacks.base
import random


class train_end_load_best_valid_loss(skorch.callbacks.base.Callback):
    def on_train_end(self, net, X, y):
        net.load_params(
            "./results/checkpoints/forces_best_{}_params.pt".format(basename)
        )


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
        fitted_params = p0
        # fitted_params = [1.88293222, 1.64105251e-7, -8.27095533e-3,
                # 4.23232128e-1, 7.35126803, -3.01956394e-3, 2.17416181,
                # 1.28998204e1, -8.49985185e-3, 4.11365936e1]
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
    sys.exit()
    unique_atoms = training_data.elements
    fp_length = training_data.fp_length
    device = "cpu"

    LR_schedule = LRScheduler("CosineAnnealingLR", T_max=5)
    cp = Checkpoint(
        monitor="forces_score_best",
        fn_prefix="./results/checkpoints/forces_best_{}_".format(basename),
    )
    load_best_valid_loss = train_end_load_best_valid_loss()
    torch.set_num_threads(1)
    model = FullNN(
            unique_atoms, [fp_length, 3, 5], device, forcetraining=forcetraining
        )
    import torch.nn as nn
    # for name, layer in model.named_modules():
        # if isinstance(layer, nn.Linear):
            # torch.nn.init.constant_(layer.weight, 0.5)
            # torch.nn.init.constant_(layer.bias, 0)

    net = NeuralNetRegressor(
        module=model,
        criterion=CustomLoss,
        criterion__force_coefficient=0.04,
        optimizer=torch.optim.LBFGS,
        # optimizer__line_search_fn="strong_wolfe",
        lr=1,
        batch_size=len(training_data),
        max_epochs=500,
        # optimizer=torch.optim.Adam,
        # lr=1e-2,
        # batch_size=32,
        # max_epochs=1000,
        iterator_train__collate_fn=collate_amp,
        iterator_train__shuffle=False,
        iterator_valid__collate_fn=collate_amp,
        iterator_valid__shuffle=False,
        device=device,
        # train_split=CVSplit(0.1),
        train_split=0,
        callbacks=[
            EpochScoring(
                forces_score,
                # on_train=False,
                on_train=True,
                use_caching=True,
                target_extractor=target_extractor,
            ),
            EpochScoring(
                energy_score,
                # on_train=False,
                on_train=True,
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
    energies = [image.get_potential_energy(apply_constraint=False) for image in
            images]
    pred_energies = [calc.get_potential_energy(image) for image in images]
    print(energies)
    print(pred_energies)
    sys.exit()

    # md_label = dir + filename
    # md_run(
        # calc=calc,
        # starting_image=images[0].copy(),
        # temp=temp,
        # count=count,
        # label=md_label,
        # dt=dt,
        # ensemble=ensemble,
    # )


def multiple_runs(
    images,
    filename,
    dir,
    ml_run,
    lj_run,
    num_iters,
    count,
    temp,
    Gs,
    forcesonly,
    ensemble,
    scaling,
    dt,
):
    file1 = open("resample_log.txt", "a")
    file1.write(filename + "\n")
    file1.close
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
        if lj_run:
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


def resample_images(base_images, sample_images, num_samples):
    random.seed(3)
    sample_points = random.sample(range(1, len(sample_images)), num_samples)
    print(sample_points)
    file1 = open("resample_log.txt", "a")
    file1.write(str(sample_points) + "\n")
    file1.close()
    images = base_images.copy()
    for idx in sample_points:
        sample_images[idx].set_calculator(EMT())
        images.append(sample_images[idx])
    return images


# Define Training data
# lang_images = ase.io.read("../datasets/COCu_lang_2ps_300K.traj", ":2000:20")
# mllj_images = ase.io.read("./lang_results/COCu_lang_2ps_res_EF_300K-LJ-1.traj", ":")
ber_images = ase.io.read("../datasets/COCu_ber_100ps_300K.traj", ":")
# mllj_images = ase.io.read("./ber_results/COCu_ber_2ps_300K-LJ-2.traj", ":")
# mllj_images1 = ase.io.read("./ber_results/COCu_ber_2ps_optim_fixmd_LJ_100_iter_1.traj", ":")
# mllj_images2 = ase.io.read("./ber_results/COCu_ber_2ps_optim_fixmd_LJ_100_iter_2.traj", ":")
# mllj_images3 = ase.io.read("./ber_results/COCu_ber_2ps_optim_fixmd_LJ_100_iter_3.traj", ":")
# mllj_images4 = ase.io.read("./ber_results/COCu_ber_2ps_optim_fixmd_LJ_100_iter_4.traj", ":")

# num_samples = 100
# resampled_lj_images = resample_images(ber_images, mllj_images, num_samples)
# resampled_lj_images1 = resample_images(resampled_lj_images, mllj_images1, num_samples)
# resampled_lj_images2 = resample_images(resampled_lj_images1, mllj_images2, num_samples)
# resampled_lj_images3 = resample_images(resampled_lj_images2, mllj_images3, num_samples)
# resampled_lj_images4 = resample_images(resampled_lj_images3, mllj_images4, num_samples)

# define symmetry functions to be used
Gs = {}
Gs["G2_etas"] = np.logspace(np.log10(0.05), np.log10(5.0), num=4)
Gs["G2_rs_s"] = [0] * 4
Gs["G4_etas"] = [0.005]
Gs["G4_zetas"] = [1.0, 4.0]
Gs["G4_gammas"] = [+1.0, -1]
Gs["cutoff"] = 5.876798323827276  # EMT asap_cutoff: False

basename = "test"
multiple_runs(
        ber_images,
        basename,
        "./",
        False,
        True,
        1,
        2000,
        300,
        Gs,
        False,
        "nvtberendsen",
        None,
        1,
    )
