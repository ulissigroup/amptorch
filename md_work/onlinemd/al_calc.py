import os
import numpy as np
import random

import torch

import ase
from ase.calculators.calculator import Calculator
from ase.calculators.emt import EMT

import skorch
from skorch import NeuralNetRegressor
from skorch.dataset import CVSplit
from skorch.callbacks import Checkpoint, EpochScoring

from amptorch.gaussian import SNN_Gaussian
from amptorch.skorch_model import AMP
from amptorch.skorch_model.utils import target_extractor, energy_score, forces_score
from amptorch.data_preprocess import AtomsDataset, collate_amp
from amptorch.model import FullNN, CustomLoss
from amptorch.morse import morse_potential
from md_work.onlinemd.md_utils import MDsimulate


__author__ = "Muhammed Shuaibi"
__email__ = "mshuaibi@andrew.cmu.edu"


class AtomisticActiveLearning(Calculator):
    """Atomistics ASE calculator
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

    def __init__(self, parent_calc, images, filename, file_dir, Gs, morse):
        Calculator.__init__(self)

        self.parent_calc = parent_calc
        self.images = images
        self.filename = filename
        self.file_dir = file_dir
        self.Gs = Gs
        self.morse = morse
        self.ml_calc = None

    def train_calc(self, images):
        os.makedirs("./results/checkpoints", exist_ok=True)
        os.makedirs(self.file_dir, exist_ok=True)

        filename = self.filename

        class train_end_load_best_valid_loss(skorch.callbacks.base.Callback):
            def on_train_end(self, net, X, y):
                net.load_params("./results/checkpoints/{}_params.pt".format(filename))

        cutoff = Gs["cutoff"]
        morse_data = None
        if self.morse:
            params = {
                "C": {"re": 0.972, "D": 6.379, "sig": 0.477},
                "O": {"re": 1.09, "D": 8.575, "sig": 0.603},
                "Cu": {"re": 2.168, "D": 3.8386, "sig": 1.696},
            }
            morse_model = morse_potential(
                images, params, cutoff, filename, combo="mean"
            )
            morse_energies, morse_forces, num_atoms = morse_model.morse_pred(
                images, params
            )
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

        return calc

    def al_random(self, images, sample_candidates, samples_to_retrain):
        random.seed(3)
        sample_points = random.sample(
            range(1, len(sample_candidates)), samples_to_retrain
        )
        for idx in sample_points:
            sample_candidates[idx].set_calculator(self.parent_calc)
            images.append(sample_candidates[idx])
        return images

    def active_learner(self, generating_function, iterations, samples_to_retrain):
        """
        Active learning method that selects data points from a pool of data
        generated by a specified generating function. Currently a random sample
        is used but other, more sophisticated techniques will be implemented.
        Termination is currently achieved based off a fixed number of
        iterations, support for more specific performance metrics will be
        introduced in the near future.

        Parameters
        ----------
        generating_function: Object. A user specified function that the active
        learning framework seeks to improve with a surrogate ML model (i.e. MD,
        NEB, relaxations) with limited data. The function must contain two
        specific methods to be used in this framework:

            (1) def run(calc, filename): A run method that runs the defined
            function and writes the corresponding Atoms object to filename.

            (2) def get_trajectory(filename, start_count, end_count, interval):
                A get_trajectory method that reads the previously written Atoms
                objects with indexing capabilities defined by start_count,
                end_count, and interval. i.e. ase.io.read(filename,
                "{}:{}:{}".format(start_count, end_count, interval))

        iterations: integer. Number of iterations the active learning framework
        is to perform.

        samples_to_retrain. integer. Number of samples to query each iteration
        of the active learning loop."""

        sample_candidates = None
        terminate = False
        iteration = 0

        while not terminate:
            if iteration > 0:
                # active learning random scheme
                self.images = self.al_random(
                    self.images, sample_candidates, samples_to_retrain
                )
            name = self.filename + "_iter_{}".format(iteration)
            # train ml calculator
            print(len(self.images))
            ml_calc = self.train_calc(images=self.images)

            # run generating function using trained ml calculator
            fn_label = self.file_dir + name
            generating_function.run(calc=ml_calc, filename=fn_label)
            # collect resulting trajectory files
            sample_candidates = generating_function.get_trajectory(
                filename=fn_label, start_count=0, end_count=-1, interval=1
            )
            iteration += 1
            # criteria to stop active learning, currently implemented
            # number of iterations
            if iteration > iterations:
                terminate = True
        self.ml_calc = ml_calc

    def calculate(self, atoms, properties, system_changes):
        Calculator.calculate(self, atoms, properties, system_changes)

        energy_pred = self.ml_calc.get_potential_energy(atoms)
        force_pred = self.ml_calc.get_forces(atoms)

        self.results["energy"] = energy_pred
        self.results["forces"] = force_pred


if __name__ == "__main__":
    Gs = {}
    Gs["G2_etas"] = np.logspace(np.log10(0.05), np.log10(5.0), num=4)
    Gs["G2_rs_s"] = [0] * 4
    Gs["G4_etas"] = [0.005]
    Gs["G4_zetas"] = [1.0, 4.0]
    Gs["G4_gammas"] = [+1.0, -1]
    Gs["cutoff"] = 5.876798323827276  # EMT asap_cutoff: False

    images = ase.io.read("../../datasets/COCu_ber_50ps_300K.traj", ":2000:10")
    al_calc = AtomisticActiveLearning(
        parent_calc=EMT(),
        images=images,
        filename="morse_test",
        file_dir="./morse/",
        Gs=Gs,
        morse=True,
    )

    al_calc.active_learner(
        generating_function=MDsimulate(
            ensemble="nvtberendsen",
            dt=1,
            temp=300,
            count=5000,
            initial_geometry=images[0].copy(),
        ),
        iterations=3,
        samples_to_retrain=100,
    )
