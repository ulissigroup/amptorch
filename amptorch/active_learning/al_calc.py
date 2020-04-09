import os
import random
import numpy as np

import torch

from ase.calculators.calculator import Calculator
from ase.calculators.singlepoint import SinglePointCalculator as sp

import skorch
from skorch import NeuralNetRegressor
from skorch.dataset import CVSplit
from skorch.callbacks import Checkpoint, EpochScoring

from amptorch.gaussian import SNN_Gaussian
from amptorch.skorch_model import AMP
from amptorch.skorch_model.utils import target_extractor, energy_score, forces_score, train_end_load_best_loss
from amptorch.data_preprocess import AtomsDataset, collate_amp
from amptorch.model import FullNN, CustomMSELoss
from amptorch.delta_models.morse import morse_potential


__author__ = "Muhammed Shuaibi"
__email__ = "mshuaibi@andrew.cmu.edu"


class AtomisticActiveLearning(Calculator):
    """Active Learning ASE calculator

    Parameters
    ----------
     parent_calc : object
         ASE parent calculator to be called for active learning queries.

     images: list
         Starting set of training images available.

     filename : str
         Label to save model and generated trajectories.

     file_dir: str
         Directory to store results.
     """

    implemented_properties = ["energy", "forces"]

    def __init__(self, parent_calc, images, filename, file_dir):
        Calculator.__init__(self)

        self.parent_calc = parent_calc
        self.images = images
        self.filename = filename
        self.file_dir = file_dir
        self.dft_calls = 0
        self.ml_calc = None

    def train_calc(
        self,
        images,
        training_params,
    ):
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

        os.makedirs("./results/checkpoints", exist_ok=True)
        os.makedirs(self.file_dir, exist_ok=True)

        filename = self.filename

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

        if train_split == 0:
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
        )

        calc = AMP(training_data, net, label=filename)
        calc.train(overwrite=True)

        return calc

    def al_random(self, images, sample_candidates, samples_to_retrain):
        """
        Randomly selects data points from a list of potential candidates to
        query and add to the existing training set.

        Paramaters
        ----------
        images: List. Current training images.

        sample_candidates: List. Potential candidates for query as collected
        from the generating function.

        samples_to_retrain: Integer. Number of samples to be randomly selected
        for query and added to the training set.
        """
        random.seed(3)
        sample_points = random.sample(
            range(1, len(sample_candidates)), samples_to_retrain
        )
        self.dft_calls += len(sample_points)
        for idx in sample_points:
            sample = sample_candidates[idx].copy()
            sample.set_calculator(self.parent_calc)
            sample_energy = sample.get_potential_energy(apply_constraint=False)
            sample_forces = sample.get_forces(apply_constraint=False)
            sample.set_calculator(
                sp(atoms=sample, energy=sample_energy, forces=sample_forces)
            )
            images.append(sample)
        return images

    def active_learner(
        self,
        generating_function,
        training_params,
    ):
        """
        Active learning method that selects data points from a pool of data
        generated by a specified generating function. Currently a random sample
        is used but other, more sophisticated techniques will be implemented.
        Termination is currently achieved based off a fixed number of
        iterations, support for more specific performance metrics will be
        introduced in the near future.

        Parameters
        ----------
        generating_function: Object.
            A user specified function that the active
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

        training_params: Dict.
            Dictionary of training parameters:

            al_convergence: Dict.
                Dictionary of stopping criteria parameters. Currently only # of
                iterations and spot verification is used:

                    Iterations: {"method": "iter", "num_iterations": int}
                        Where "num_iteratiosn" is the number of iterations the AL
                        framework is run before termination.
                    Spot: {"method": "spot", "num2verify": int }
                        Where "num2verify" is the number of points necessary to be randomly sampled
                        and verified with DFT before convergence.
            samples_to_retrain. integer.
                Number of samples to query each iteration
                of the active learning loop.
            Gs: dictionary.
                Dictionary of symmetry functions to be used for
                fingerprinting.
            morse: boolean
                Delta-ML model with morse potential to be used.
            morse_params: dictionary
                If morse True, define morse potential parameters for all elements in
                system.
            forcetraining: boolean
                Use forces in training scheme alongside energies.
            force_coefficient: float
                Define the force coefficient to be utilized in the loss function. A
                coefficient > 0 indicates force training is turned on.
            cores: int
                Number of cores to parallelize across for fingerprint computation
            criterion: object
                Specify the loss function to be optimized.
                default: CustomMSELoss
            num_layers: int
                Number of layers in each atomic neural network.
            num_nodes: int
                Number of nodes in each layer in each atomic neural network.
            learning_rater: float
                Define the model learning rate.
            epochs: int
                Number of training epochs to use.
            train_split: float or int
                Training split to be used for validation. If float, corresponding
                proportion of training set will be used. If int, k-fold cross
                validation will be used."""

        al_convergence = training_params["al_convergence"]
        samples_to_retrain = training_params["samples_to_retrain"]
        test_split = training_params["test_split"]

        sample_candidates = None
        terminate = False
        iteration = 0

        while not terminate:
            #TODO train_split check for initial starting images
            if iteration == 0:
                training_params["test_split"] = 0
            if iteration > 0:
                # active learning random scheme
                self.images = self.al_random(
                    self.images, sample_candidates, samples_to_retrain
                )
                training_params["test_split"] = test_split
            name = self.filename + "_iter_{}".format(iteration)
            # train ml calculator
            self.ml_calc = self.train_calc(
                images=self.images,
                training_params=training_params
            )

            # run generating function using trained ml calculator
            fn_label = self.file_dir + name
            generating_function.run(calc=self.ml_calc, filename=fn_label)
            # collect resulting trajectory files
            sample_candidates = generating_function.get_trajectory(
                filename=fn_label, start_count=0, end_count=-1, interval=1
            )
            iteration += 1
            # criteria to stop active learning
            #TODO Find a better way to structure this.
            method = al_convergence["method"]
            if method == 'iter':
                 kwargs = {"current_i": iteration, "total_i":
                         al_convergence["num_iterations"]}
            elif method == 'spot':
                kwargs = {"images": sample_candidates, "num2verify":
                        al_convergence["num2verify"]}

            terminate = self.termination_criteria(method=method, **kwargs)
        print('Totale # of parent calculator calls: {}'.format(self.dft_calls))

    def termination_criteria(self, method='iter', **kwargs):
        """Criteria for AL termination
       Parameters
       ----------
       method: str
            Method for termination of active learning loop.

            'iter': Terminates after specified number of iterations.

                    args: (current iteration, # of iterations)

            'spot': Terminates after specified number of DFT verifications
                    are consistent with ML predictions.

                    args: (images, # of spot checks)
        """
        terminate = False

        if method == "iter":
            current_i = kwargs["current_i"]
            total_i = kwargs["total_i"]
            if current_i > total_i:
                terminate = True

        if method == "spot":
            images = kwargs["images"]
            num2verify = kwargs["num2verify"]
            points2verify = random.sample(range(len(images)), num2verify)
            for idx in points2verify:
                sample = images[idx].copy()
                sample.set_calculator(self.parent_calc)
                sample_energy = sample.get_potential_energy(apply_constraint=False)
                sample_forces = sample.get_forces(apply_constraint=False)
                ml_energy = self.ml_calc.get_potential_energy(sample)
                ml_forces = self.ml_calc.get_forces(sample)

                energy_consistency = np.abs(sample_energy - ml_energy) <= 0.1
                force_consistency = (np.abs(sample_forces - ml_forces) <=
                        0.1).all()
                consistent = energy_consistency and force_consistency
                if consistent:
                    terminate = True
                    break

                #TODO possibly use these calculations and add them to model as well
        return terminate

    def calculate(self, atoms, properties, system_changes):
        """
        ASE required calculate method to compute specified properties.
        """
        Calculator.calculate(self, atoms, properties, system_changes)

        energy_pred = self.ml_calc.get_potential_energy(atoms)
        force_pred = self.ml_calc.get_forces(atoms)

        self.results["energy"] = energy_pred
        self.results["forces"] = force_pred
