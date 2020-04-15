import copy
import os
import numpy as np
import torch.multiprocessing as mp
from torch.multiprocessing import Pool
import random

import torch
from torch.nn import init

import ase
from ase.calculators.emt import EMT
from ase.calculators.singlepoint import SinglePointCalculator as sp
from ase.calculators.calculator import Calculator, Parameters
from ase.md.nvtberendsen import NVTBerendsen
from ase.build import fcc100, add_adsorbate, molecule
from ase.constraints import FixAtoms
from ase.optimize import BFGS, QuasiNewton

import skorch
from skorch import NeuralNetRegressor
from skorch.dataset import CVSplit
from skorch.callbacks import Checkpoint, EpochScoring
from skorch.callbacks.lr_scheduler import LRScheduler

from amptorch.gaussian import SNN_Gaussian
from amptorch.skorch_model import AMP
from amptorch.skorch_model.utils import target_extractor, energy_score, forces_score, train_end_load_best_loss
from amptorch.delta_models.morse import morse_potential
from amptorch.utils import make_amp_descriptors_simple_nn
from amptorch.data_preprocess import AtomsDataset, collate_amp
from amptorch.model import FullNN, CustomMSELoss
from amptorch.active_learning.atomistic_methods import MDsimulate, Relaxation
from amptorch.active_learning.ensemble_calc import EnsembleCalc


def bootstrap_ensemble(parent_dataset, resampled_set=None, new_data=None, n_ensembles=1):
    if len(parent_dataset) == 1 and new_data is None :
        ensemble_sets = [parent_dataset.copy() for i in range(n_ensembles)]
        return ensemble_sets, parent_dataset
    ensemble_sets = []
    if new_data is not None and resampled_set is not None:
        n_ensembles = len(resampled_set)
        parent_dataset.append(new_data)
        for i in range(n_ensembles):
            resampled_set[i].append(random.sample(parent_dataset, 1)[0])
            for k in range(len(resampled_set[i])):
                if random.random() < 1 / len(resampled_set[i]):
                    resampled_set[i][k] = new_data
            ensemble_sets.append(resampled_set[i])
    else:
        for i in range(n_ensembles):
            ensemble_sets.append(
                    random.choices(parent_dataset, k=len(parent_dataset))
                    )
    return ensemble_sets, parent_dataset

def calculate_stats(energies, forces):
    energy_median = np.median(energies)
    energy_var = np.var(energies)
    forces_median = np.median(forces, axis=0)
    max_forces_var = np.max(np.var(forces, axis=0))
    return energy_median, forces_median, max_forces_var

def train_calcs(training_data, training_params, ensemble, ncores):
    if ensemble:
        assert(len(training_data) == ensemble), 'Invalid # ensemble training_data'
        trained_calcs = ensemble_trainer(training_data, training_params, ncores)
        return trained_calcs

    calc_parameters = model_trainer(training_data, training_params)
    trained_calc = AMP(*calc_parameters)
    return trained_calc

def ensemble_trainer(ensemble_datasets, training_params, ncores):
    if ncores == 'max':
        ncores = len(ensemble_datasets)
    pool = Pool(ncores)

    input_data = []
    for _, dataset in enumerate(ensemble_datasets):
        parallel_params = training_params.copy()
        parallel_params["filename"] += str(_)
        inputs = (dataset, parallel_params)
        input_data.append(inputs)
    calc_parameters = pool.starmap(model_trainer, input_data)
    trained_calcs = [AMP(*params) for params in calc_parameters]
    ensemble_calc = EnsembleCalc(trained_calcs, training_params)
    return ensemble_calc

def model_trainer(images, training_params):
    Gs = training_params["Gs"]
    morse = training_params["morse"]
    morse_params = training_params["morse_params"]
    forcetraining = training_params["forcetraining"]
    cores = training_params["cores"]
    optimizer = training_params["optimizer"]
    batch_size = training_params["batch_size"]
    criterion = training_params["criterion"]
    num_layers = training_params["num_layers"]
    num_nodes = training_params["num_nodes"]
    force_coefficient = training_params["force_coefficient"]
    learning_rate = training_params["learning_rate"]
    epochs = training_params["epochs"]
    train_split = training_params["test_split"]
    shuffle = training_params["shuffle"]
    filename = training_params["filename"]
    verbose = training_params["verbose"]

    os.makedirs("./results/checkpoints", exist_ok=True)

    cutoff = Gs["cutoff"]
    morse_data = None
    device = "cpu"
    if morse:
        params = morse_params
        morse_model = morse_potential(
                images, params, cutoff, filename, combo="mean"
                )
        morse_energies, morse_forces, num_atoms = morse_model.morse_pred(
                images, params
                )
        morse_data = [morse_energies, morse_forces, num_atoms, params, morse_model]

    forcetraining = forcetraining
    training_data = AtomsDataset(
            images,
            SNN_Gaussian,
            Gs,
            forcetraining=forcetraining,
            label=filename,
            cores=cores,
            delta_data=morse_data,
            )
    unique_atoms = training_data.elements
    fp_length = training_data.fp_length

    load_best_valid_loss = train_end_load_best_loss(filename)
    torch.set_num_threads(1)

    if train_split == 0 or len(images)*train_split < 1:
        train_split = 0
        on_train = True
    else:
        train_split = CVSplit(cv=train_split)
        on_train = False

    if forcetraining:
        force_coefficient = force_coefficient
        callbacks = [
                EpochScoring(
                    forces_score,
                    on_train=on_train,
                    use_caching=True,
                    target_extractor=target_extractor,
                    ),
                EpochScoring(
                    energy_score,
                    on_train=on_train,
                    use_caching=True,
                    target_extractor=target_extractor,
                    ),
                Checkpoint(
                    monitor="forces_score_best",
                    fn_prefix="./results/checkpoints/{}_".format(filename),
                    ),
                load_best_valid_loss,
                ]
    else:
        force_coefficient = 0
        callbacks = [
                EpochScoring(
                    energy_score,
                    on_train=on_train,
                    use_caching=True,
                    target_extractor=target_extractor,
                    ),
                Checkpoint(
                    monitor="energy_score_best",
                    fn_prefix="./results/checkpoints/{}_".format(filename),
                    ),
                load_best_valid_loss,
                ]

    net = NeuralNetRegressor(
            module=FullNN(
                unique_atoms,
                [fp_length, num_layers, num_nodes],
                device,
                forcetraining=forcetraining,
                ),
            criterion=criterion,
            criterion__force_coefficient=force_coefficient,
            optimizer=optimizer,
            lr=learning_rate,
            batch_size=batch_size,
            max_epochs=epochs,
            iterator_train__collate_fn=collate_amp,
            iterator_train__shuffle=shuffle,
            iterator_valid__collate_fn=collate_amp,
            iterator_valid__shuffle=False,
            device=device,
            train_split=train_split,
            callbacks=callbacks,
            verbose=verbose,
            )


    calc = AMP(training_data, net, label=filename)
    calc.train(overwrite=True)

    return [training_data, net, filename]
