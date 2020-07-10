import copy
import sys
import os
import random
import numpy as np

import torch

import ase.db
from ase.calculators.calculator import Calculator
from ase.calculators.singlepoint import SinglePointCalculator as sp


import skorch
from skorch import NeuralNetRegressor
from skorch.dataset import CVSplit
from skorch.callbacks import Checkpoint, EpochScoring

from amptorch.gaussian import SNN_Gaussian
from amptorch.skorch_model import AMP
from amptorch.skorch_model.utils import (
    target_extractor,
    energy_score,
    forces_score,
    train_end_load_best_loss,
)
from amptorch.data_preprocess import AtomsDataset, collate_amp
from amptorch.model import BPNN, CustomMSELoss
from amptorch.delta_models.morse import morse_potential
from amptorch.active_learning.al_utils import write_to_db
from amptorch.active_learning.trainer import train_calcs
from amptorch.active_learning.bootstrap import bootstrap_ensemble
from amptorch.active_learning.query_methods import termination_criteria


__author__ = "Muhammed Shuaibi"
__email__ = "mshuaibi@andrew.cmu.edu"


class AtomisticActiveLearner:
    """Active Learner

    Parameters
    ----------
     training_data: list
        List of Atoms objects representing the initial dataset.

    training_params: dict
        Dictionary of training parameters and model settings.

    parent_calc: object.
        Calculator to be used for querying calculations.

    ensemble: boolean.
        Whether to train an ensemble of models to make predictions. ensemble
        must be True if uncertainty based query methods are to be used.
     """

    implemented_properties = ["energy", "forces"]

    def __init__(self, training_data, training_params, parent_calc, ensemble=False):
        self.training_data = copy.deepcopy(training_data)
        self.training_params = training_params
        self.parent_calc = parent_calc
        self.ensemble = ensemble
        self.parent_calls = 0
        self.iteration = 0

        if ensemble:
            assert isinstance(ensemble, int) and ensemble > 1, "Invalid ensemble!"
            self.training_data, self.parent_dataset = bootstrap_ensemble(
                self.training_data, n_ensembles=ensemble
            )
        else:
            self.parent_dataset = self.training_data

    def learn(self, atomistic_method, query_strategy):
        db_filepath = self.training_params["db_filepath"]
        db = connect("{}".format(db_filepath))
        def set_database():
            database = []
            for row in db.select():
                a = db.get_atoms(id=row.id)
                database.append(a)
        return database
        al_convergence = self.training_params["al_convergence"]
        samples_to_retrain = self.training_params["samples_to_retrain"]
        filename = self.training_params["filename"]
        file_dir = self.training_params["file_dir"]
        queries_db = ase.db.connect("{}.db".format(filename))
        os.makedirs(file_dir, exist_ok=True)

        terminate = False

        while not terminate:
            fn_label = f"{file_dir}{filename}_iter_{self.iteration}"
            # active learning random scheme
            if self.iteration > 0:
                queried_images = query_strategy(
                    self.parent_dataset,
                    sample_candidates,
                    samples_to_retrain,
                    parent_calc=self.parent_calc,
                )
                write_to_db(queries_db, queried_images)
                self.parent_dataset, self.training_data = self.add_data(queried_images)
                self.parent_calls += len(queried_images)

            # train ml calculator
            if len(set_database()) > 0:
                 trained_calc = train_calcs(
                   training_data=set_database(),
                   training_params=self.training_params,
                   ensemble=self.ensemble,
                   ncores=self.training_params["cores"],
                )
            else:
                 trained_calc = train_calcs(
                     training_data=self.training_data,
                     training_params=self.training_params,
                     ensemble=self.ensemble,
                     ncores=self.training_params["cores"],
                )
            # run atomistic_method using trained ml calculator
            atomistic_method.run(calc=trained_calc, filename=fn_label)
            # collect resulting trajectory files
            sample_candidates = atomistic_method.get_trajectory(filename=fn_label)
            self.iteration += 1
            # criteria to stop active learning
            # TODO Find a better way to structure this.
            method = al_convergence["method"]
            if method == "iter":
                termination_args = {
                    "current_i": self.iteration,
                    "total_i": al_convergence["num_iterations"],
                }
            elif method == "final":
                termination_args = {
                    "images": sample_candidates,
                    "calc": self.parent_calc,
                    "energy_tol": al_convergence["energy_tol"],
                    "force_tol": al_convergence["force_tol"],
                }
                self.parent_calls += 1

            terminate = termination_criteria(
                method=method, termination_args=termination_args
            )

        return trained_calc

    def add_data(self, queried_images):
        if self.ensemble:
            for query in queried_images:
                self.training_data, self.parent_dataset = bootstrap_ensemble(
                    self.parent_dataset,
                    self.training_data,
                    query,
                    n_ensembles=self.ensemble,
                )
        else:
            self.training_data += queried_images
        return self.parent_dataset, self.training_data
