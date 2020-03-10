"""Base AMP calculator class to be utilized similiar to other ASE calculators"""

import os
import copy
import numpy as np
import time
from torch.utils.data import DataLoader
from amptorch.utils import Logger, hash_images, get_hash
from amptorch.skorch_model.utils import (
    make_force_header,
    make_energy_header,
    make_val_force_header,
    make_val_energy_header,
    log_results,
)
from amptorch.data_preprocess import (
    AtomsDataset,
    factorize_data,
    collate_amp,
    TestDataset,
)
from amptorch.model import FullNN, CustomLoss
from ase.calculators.calculator import Calculator, Parameters
import torch
import multiprocessing
import random
from amptorch.utils import make_amp_descriptors_simple_nn
import sys
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase.md.nvtberendsen import NVTBerendsen
from ase.calculators.emt import EMT
from ase.calculators.singlepoint import SinglePointCalculator as sp
from ase.md import Langevin
from ase import units
from amptorch.lj_model import lj_optim
import ase
from skorch import NeuralNetRegressor
from skorch.dataset import CVSplit
from skorch.callbacks import Checkpoint, EpochScoring
from skorch.callbacks.lr_scheduler import LRScheduler
from amptorch.gaussian import SNN_Gaussian
from amptorch.skorch_model import AMP
from amptorch.skorch_model.utils import target_extractor, energy_score, forces_score
from torch.nn import init
from skorch.utils import to_numpy


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

    def __init__(self, parent_dataset, n_ensembles):
        Calculator.__init__(self)

        self.n_ensembles = n_ensembles
        self.ensemble_sets, self.parent_dataset = self.bootstrap_ensemble(
            parent_dataset, n_ensembles=n_ensembles
        )
        self.trained_calcs = self.construct_calc(train(self.ensemble_sets))
        self.uncertain_tol = 0.01
        # self.data_to_retrain = []
        # self.train_trigger = 10

    def bootstrap_ensemble(
        self, parent_dataset, resampled_set=None, new_data=None, n_ensembles=1
    ):
        if len(parent_dataset) == 1:
            return [parent_dataset] * n_ensembles
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
        return energy_mean, forces_mean, energy_var

    def calculate(self, atoms, properties, system_changes):
        Calculator.calculate(self, atoms, properties, system_changes)

        energies = []
        forces = []
        for calc in self.trained_calcs:
            energies.append(calc.get_potential_energy(atoms))
            forces.append(calc.get_forces(atoms))
        energies = np.array(energies)
        forces = np.array(forces)
        energy_pred, force_pred, uncertainty = self.calculate_stats(energies, forces)

        # if uncertainty >= self.uncertain_tol:
        if uncertainty != 0:
            new_data = atoms.copy()
            new_data.set_calculator(EMT())
            self.results["energy"] = new_data.get_potential_energy(apply_constraint=False)
            self.results["forces"] = new_data.get_forces(apply_constraint=False)
            self.ensemble_sets, self.parent_dataset = self.bootstrap_ensemble(
                self.parent_dataset, self.ensemble_sets, new_data=new_data
            )
            self.trained_calcs = self.construct_calc(train(self.ensemble_sets))
            sys.exit()
        else:
            self.results["energy"] = float(energy_pred)
            self.results["forces"] = force_pred

def train(ensemble_sets):
    Gs = {}
    Gs["G2_etas"] = np.logspace(np.log10(0.05), np.log10(5.0), num=4)
    Gs["G2_rs_s"] = [0] * 4
    Gs["G4_etas"] = [0.005]
    Gs["G4_zetas"] = [1.0, 4.0]
    Gs["G4_gammas"] = [+1.0, -1]
    Gs["cutoff"] = 5.876798323827276  # EMT asap_cutoff: False

    pool = multiprocessing.Pool(len(ensemble_sets))

    input_data = []
    # results = []
    for _ in range(len(ensemble_sets)):
        inputs = [ensemble_sets[_], "test" + str(_), "./", Gs, True, True, None]
        # results.append(train_calc(inputs))
        # inputs = [ensemble_sets[_], "test" + str(_), "./", Gs, True, True, None]
        input_data.append(inputs)
    results = pool.map(train_calc, input_data)
    return results

def train_calc(inputs):
    images, filename, file_dir, Gs, lj, forcesonly, scaling = inputs
    if not os.path.exists(file_dir):
        os.makedirs(file_dir, exist_ok=True)
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
        # fitted_params = lj_model.fit()
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
        max_epochs=3,
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
        ],
    )
    calc = AMP(training_data, net, label=filename)
    calc.train()
    return [training_data, net, filename]

Gs = {}
Gs["G2_etas"] = np.logspace(np.log10(0.05), np.log10(5.0), num=4)
Gs["G2_rs_s"] = [0] * 4
Gs["G4_etas"] = [0.005]
Gs["G4_zetas"] = [1.0, 4.0]
Gs["G4_gammas"] = [+1.0, -1]
Gs["cutoff"] = 5.876798323827276  # EMT asap_cutoff: False
images = ase.io.read("../../datasets/COCu_ber_100ps_300K.traj", ":10")
elements = np.array([atom.symbol for atoms in images for atom in atoms])
_, idx = np.unique(elements, return_index=True)
elements = list(elements[np.sort(idx)])
make_amp_descriptors_simple_nn(images, Gs, elements, cores=10, label="test")


def md_run2(
    parent_dataset,
    dft_calculator,
    starting_image,
    temp,
    dt,
    count,
    label,
    ensemble="NVE",
    new_data=None,
):
    slab = starting_image.copy()
    slab.set_calculator(AMPOnlineCalc(parent_dataset, 3))
    # MaxwellBoltzmannDistribution(slab, temp * units.kB)
    if ensemble == "NVE":
        dyn = VelocityVerlet(slab, dt * units.fs)
    elif ensemble == "nvtberendsen":
        dyn = NVTBerendsen(slab, dt * units.fs, temp, taut=300 * units.fs)
    elif ensemble == "langevin":
        dyn = Langevin(slab, dt * units.fs, temp * units.kB, 0.002)
    traj = ase.io.Trajectory(label + ".traj", "w", slab)
    dyn.attach(traj.write, interval=1)
    dyn.run(count - 1)


md_run2(images, EMT(), images[0], 300, 1, 5, "test", "nvtberendsen")
