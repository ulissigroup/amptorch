import sys
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
from amptorch.active_learning.generator_funcs import MDsimulate, Relaxation


__author__ = "Muhammed Shuaibi"
__email__ = "mshuaibi@andrew.cmu.edu"

class AMPOnlineCalc(Calculator):
    """Atomistics Machine-Learning Potential (AMP) ASE calculator
   Parameters
   ----------
    model : object
        Class representing the regression model. Input arguments include training
        images, descriptor type, and force_coefficient. Model structure and training schemes can be
        modified directly within the class.

    label : str
        Location to save the trained model.

    """

    implemented_properties = ["energy", "forces"]

    def __init__(self,
            parent_dataset,
            parent_calc,
            n_ensembles,
            n_cores,
            training_params):
        Calculator.__init__(self)

        self.n_ensembles = n_ensembles
        self.parent_calc = parent_calc
        self.training_params = training_params
        self.n_cores = n_cores
        self.elements, self.Gs = self.fingerprint_args(parent_dataset)
        self.make_fps(parent_dataset)
        self.ensemble_sets, self.parent_dataset = self.bootstrap_ensemble(
                parent_dataset, n_ensembles=n_ensembles
                )
        self.trained_calcs = self.construct_calc(self.parallel_trainer())
        self.uncertain_tol = training_params["uncertain_tol"]
        self.parent_calls = 0

    def fingerprint_args(self, images):
        elements = np.array([atom.symbol for atoms in images for atom in atoms])
        _, idx = np.unique(elements, return_index=True)
        elements = list(elements[np.sort(idx)])
        Gs = self.training_params["Gs"]
        return elements, Gs

    def make_fps(self, atoms):
        try:
            len(atoms)
        except:
            atoms = [atoms]
            pass
        make_amp_descriptors_simple_nn(atoms, self.Gs, self.elements, cores=1, label="oal")

    def bootstrap_ensemble(self, parent_dataset, resampled_set=None, new_data=None, n_ensembles=1):
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

    def calculate_stats(self, energies, forces):
        energy_median = np.median(energies)
        energy_var = np.var(energies)
        forces_median = np.median(forces, axis=0)
        max_forces_var = np.max(np.var(forces, axis=0))
        return energy_median, forces_median, max_forces_var

    def parallel_trainer(self):
        if self.n_cores == 'max':
            self.n_cores = len(self.ensemble_sets)
        pool = Pool(self.n_cores)

        input_data = []
        for _ in range(len(self.ensemble_sets)):
            parallel_params = self.training_params.copy()
            parallel_params["filename"] += str(_)
            inputs = [self.ensemble_sets[_], parallel_params]
            input_data.append(inputs)
        results = pool.map(train_calc, input_data)
        return results

    def construct_calc(self, calc_parameters):
        calcs = []
        for _ in  range(len(calc_parameters)):
            calc = AMP(
                calc_parameters[_][0], calc_parameters[_][1], calc_parameters[_][2]
            )
            calcs.append(calc)
        return calcs

    def calculate(self, atoms, properties, system_changes):
        Calculator.calculate(self, atoms, properties, system_changes)

        energies = []
        forces = []

        self.make_fps(atoms)
        for calc in self.trained_calcs:
            energies.append(calc.get_potential_energy(atoms))
            forces.append(calc.get_forces(atoms))
        energies = np.array(energies)
        forces = np.array(forces)
        energy_pred, force_pred, uncertainty = self.calculate_stats(energies, forces)

        if uncertainty >= self.uncertain_tol:
            new_data = atoms.copy()
            new_data.set_calculator(self.parent_calc)

            energy_pred = new_data.get_potential_energy(apply_constraint=False)
            force_pred = new_data.get_forces(apply_constraint=False)
            new_data.set_calculator(sp(atoms=new_data, energy=energy_pred,
                forces=force_pred))

            self.ensemble_sets, self.parent_dataset = self.bootstrap_ensemble(
                    self.parent_dataset, self.ensemble_sets, new_data=new_data
                    )

            self.trained_calcs = self.construct_calc(self.parallel_trainer())
            self.parent_calls += 1

        self.results["energy"] = energy_pred
        self.results["forces"] = force_pred

#TODO: construct core trainer
def train_calc(inputs):
    images, training_params = inputs

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
    device = "cpu"

    load_best_valid_loss = train_end_load_best_loss(filename)
    torch.set_num_threads(1)

    if train_split == 0 or len(images)*train_split < 1:
        on_train = True
    else:
        train_split = CVSplit(cv=train_split, random_state=1)
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


if __name__ == "__main__":
    mp.set_start_method("spawn")

    Gs = {}
    Gs["G2_etas"] = np.logspace(np.log10(0.05), np.log10(5.0), num=4)
    Gs["G2_rs_s"] = [0] * 4
    Gs["G4_etas"] = [0.005]
    Gs["G4_zetas"] = [1.0, 4.0]
    Gs["G4_gammas"] = [+1.0, -1]
    Gs["cutoff"] = 5.876798323827276  # EMT asap_cutoff: False

    slab = fcc100("Cu", size=(3, 3, 3))
    adsorbate = molecule("C")
    add_adsorbate(slab, adsorbate, 3, offset=(1,1))
    constraints = FixAtoms(
            indices=[atom.index for atom in slab if (atom.tag == 3)])
    slab.set_constraint(constraints)
    slab.center(vacuum=13.0, axis=2)
    slab.set_pbc(True)
    slab.wrap(pbc=[True] * 3)
    slab.set_calculator(EMT())

    images = [slab]
    morse_params = {
        "C": {"re": 0.972, "D": 6.379, "sig": 0.477},
        "Cu": {"re": 2.168, "D": 3.8386, "sig": 1.696},
    }

    training_params = {
        "uncertain_tol": 0.001,
        "Gs": Gs,
        "morse": True,
        "morse_params": morse_params,
        "forcetraining": True,
        "cores": 10,
        "optimizer": torch.optim.LBFGS,
        "batch_size": 1000,
        "criterion": CustomMSELoss,
        "num_layers": 3,
        "num_nodes": 20,
        "force_coefficient": 0.04,
        "learning_rate": 1e-1,
        "epochs": 20,
        "test_split": 0,
        "shuffle": False,
        "filename": "oal_test",
    }

    structure_optim = Relaxation(slab, BFGS, fmax=0.05, steps=None)
    online_calc = AMPOnlineCalc(images, EMT(), 3, training_params)
    structure_optim.run(online_calc, filename='relax_oal')
