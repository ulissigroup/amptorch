import sys
import os
import numpy as np
import multiprocessing
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
from amptorch.skorch_model.utils import target_extractor, energy_score, forces_score
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

    def __init__(self, parent_dataset, parent_calc, n_ensembles):
        Calculator.__init__(self)

        self.n_ensembles = n_ensembles
        self.parent_calc = parent_calc
        self.ensemble_sets, self.parent_dataset = self.bootstrap_ensemble(
            parent_dataset, n_ensembles=n_ensembles
        )
        self.trained_calcs = self.construct_calc(train(self.ensemble_sets))
        self.uncertain_tol = 0.001
        self.dft_calls = 0

    def bootstrap_ensemble(
        self, parent_dataset, resampled_set=None, new_data=None, n_ensembles=1
    ):
        if len(parent_dataset) == 1 and new_data is None :
            return [parent_dataset.copy()] * n_ensembles, parent_dataset
        ensemble_sets = []
        if new_data is not None and resampled_set is not None:
            n_ensembles = len(resampled_set)
            parent_dataset.append(new_data)
            for i in range(n_ensembles):
                resampled_set[i].append(random.sample(parent_dataset, 1)[0])
                for k in range(0, len(resampled_set[i]) - 1):
                    if random.random() < 1 / len(resampled_set[i]):
                        resampled_set[i][k] = new_data
                ensemble_sets.append(resampled_set[i])
        else:
            for i in range(n_ensembles):
                ensemble_sets.append(
                    random.choices(parent_dataset, k=len(parent_dataset))
                )
        return ensemble_sets, parent_dataset

    def construct_calc(self, calc_parameters):
        calcs = []
        for _ in range(len(calc_parameters)):
            calc = AMP(
                calc_parameters[_][0], calc_parameters[_][1], calc_parameters[_][2]
            )
            calcs.append(calc)
        return calcs

    def calculate_stats(self, energies, forces):
        energy_mean = np.mean(energies)
        energy_var = np.var(energies)
        forces_mean = np.mean(forces, axis=0)
        forces_var = np.var(forces, axis=0)
        return energy_mean, forces_mean, energy_var

    def calculate(self, atoms, properties, system_changes):
        Calculator.calculate(self, atoms, properties, system_changes)

        energies = []
        forces = []

        make_amp_descriptors_simple_nn([atoms], Gs, elements, cores=10, label="test")
        for calc in self.trained_calcs:
            energies.append(calc.get_potential_energy(atoms))
            forces.append(calc.get_forces(atoms))
        energies = np.array(energies)
        forces = np.array(forces)
        energy_pred, force_pred, uncertainty = self.calculate_stats(energies, forces)

        # if uncertainty == self.uncertain_tol:
        if uncertainty >= 0.001:
            new_data = atoms.copy()
            new_data.set_calculator(self.parent_calc)

            sample_energy = new_data.get_potential_energy(apply_constraint=False)
            sample_forces = new_data.get_forces(apply_constraint=False)
            new_data.set_calculator(sp(atoms=new_data, energy=sample_energy,
                forces=sample_forces))

            self.ensemble_sets, self.parent_dataset = self.bootstrap_ensemble(
                self.parent_dataset, self.ensemble_sets, new_data=new_data
            )

            self.trained_calcs = self.construct_calc(train(self.ensemble_sets))
            self.dft_calls += 1

            self.results["energy"] = sample_energy
            self.results["forces"] = sample_forces

        else:
            self.results["energy"] = energy_pred
            self.results["forces"] = force_pred
            print('step')


def train(ensemble_sets):
    pool = multiprocessing.Pool(len(ensemble_sets))

    input_data = []
    for _ in range(len(ensemble_sets)):
        inputs = [ensemble_sets[_], "test" + str(_), "./", Gs, True, True, None]
        input_data.append(inputs)
    results = pool.map(train_calc, input_data)
    return results


def train_calc(inputs):
    images, filename, file_dir, Gs, delta, forcesonly, scaling = inputs
    class train_end_load_best_valid_loss(skorch.callbacks.base.Callback):
            def on_train_end(self, net, X, y):
                net.load_params("./results/checkpoints/{}_params.pt".format(filename))
    cp = Checkpoint(
            monitor="forces_score_best",
            fn_prefix="./results/checkpoints/{}_".format(filename),
        )

    if not os.path.exists(file_dir):
        os.makedirs(file_dir, exist_ok=True)

    forcetraining = True
    training_data = AtomsDataset(
        images,
        SNN_Gaussian,
        Gs,
        forcetraining=forcetraining,
        label=filename,
        cores=1,
        delta_data=None,
    )
    unique_atoms = training_data.elements
    fp_length = training_data.fp_length
    device = "cpu"

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
        max_epochs=100,
        iterator_train__collate_fn=collate_amp,
        iterator_train__shuffle=False,
        iterator_valid__collate_fn=collate_amp,
        iterator_valid__shuffle=False,
        device=device,
        verbose=0,
        # train_split=CVSplit(cv=5, random_state=1),
        train_split=0,
        callbacks=[
            EpochScoring(
                forces_score,
                on_train=True,
                use_caching=True,
                target_extractor=target_extractor,
            ),
            EpochScoring(
                energy_score,
                on_train=True,
                use_caching=True,
                target_extractor=target_extractor,
            ),
        ],
    )
    calc = AMP(training_data, net, label=filename)
    calc.train()
    return [training_data, net, filename]

if __name__ == "__main__":
    multiprocessing.set_start_method("spawn")

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

    elements = np.array([atom.symbol for atoms in images for atom in atoms])
    _, idx = np.unique(elements, return_index=True)
    elements = list(elements[np.sort(idx)])
    make_amp_descriptors_simple_nn(images, Gs, elements, cores=10, label="test")

    structure_optim = Relaxation(slab, BFGS, fmax=0.05, steps=None)
    online_calc = AMPOnlineCalc(images, EMT(), 3)
    structure_optim.run(online_calc, filename='relax_oal')
    print(online_calc.dft_calls)
