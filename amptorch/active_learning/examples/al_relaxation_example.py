import sys
import numpy as np
import random
import torch

import ase
from ase.calculators.calculator import Calculator
from ase.calculators.emt import EMT
from ase.calculators.singlepoint import SinglePointCalculator as sp
from ase.build import fcc100, add_adsorbate, molecule
from ase.constraints import FixAtoms
from ase.optimize import BFGS, BFGSLineSearch

from amptorch.active_learning.atomistic_methods import MDsimulate, Relaxation
from amptorch.active_learning.learner import AtomisticActiveLearner
from amptorch.active_learning.query_methods import random_query, max_uncertainty
from amptorch.model import CustomMSELoss

import multiprocessing as mp

if __name__ == "__main__":
    random.seed(1)
    mp.set_start_method("spawn")
    # Define initial set of images, can be as few as 1. If 1, make sure to
    slab = fcc100("Cu", size=(3, 3, 3))
    ads = molecule("C")
    add_adsorbate(slab, ads, 3, offset=(1, 1))
    cons = FixAtoms(indices=[atom.index for atom in slab if (atom.tag == 3)])
    slab.set_constraint(cons)
    slab.center(vacuum=13.0, axis=2)
    slab.set_pbc(True)
    slab.wrap(pbc=[True] * 3)
    slab.set_calculator(EMT())
    sample_energy = slab.get_potential_energy(apply_constraint=False)
    sample_forces = slab.get_forces(apply_constraint=False)
    slab.set_calculator(sp(atoms=slab, energy=sample_energy, forces=sample_forces))

    images = [slab]

    # Define symmetry functions
    Gs = {}
    Gs["G2_etas"] = np.logspace(np.log10(0.05), np.log10(5.0), num=4)
    Gs["G2_rs_s"] = [0] * 4
    Gs["G4_etas"] = [0.005]
    Gs["G4_zetas"] = [1.0, 4.0]
    Gs["G4_gammas"] = [+1.0, -1]
    Gs["cutoff"] = 5.876798323827276  # EMT asap_cutoff: False

    training_params = {
        "al_convergence": {"method": "iter", "num_iterations": 3},
        "samples_to_retrain": 5,
        "Gs": Gs,
        "morse": True,
        "forcetraining": True,
        "cores": 10,
        "optimizer": torch.optim.LBFGS,
        "batch_size": 1000,
        "criterion": CustomMSELoss,
        "num_layers": 3,
        "num_nodes": 20,
        "force_coefficient": 0.04,
        "learning_rate": 1e-2,
        "epochs": 100,
        "test_split": 0,
        "shuffle": False,
        "verbose": 1,
        "filename": "relax_example",
        "file_dir": "./",
    }

    # Define AL calculator
    learner = AtomisticActiveLearner(
        training_data=images,
        training_params=training_params,
        parent_calc=EMT(),
        ensemble=False,
    )
    learner.learn(
        atomistic_method=Relaxation(
            initial_geometry=images[0].copy(), optimizer=BFGS, fmax=0.05, steps=50
        ),
        query_strategy=random_query,
    )

    # Calculate true relaxation
    al_iterations = learner.iteration - 1
    file_path = training_params["file_dir"]+training_params["filename"]
    true_relax = Relaxation(slab, BFGS)
    true_relax.run(EMT(), "true_relax")
    parent_calc_traj = true_relax.get_trajectory("true_relax", 0, -1, 1)
    final_ml_traj = ase.io.read("{}_iter_{}.traj".format(file_path, al_iterations), ":")

    # Compute ML predicted energies
    ml_relaxation_energies = [image.get_potential_energy() for image in final_ml_traj]
    # Compute actual (EMT) energies for ML predicted structures
    emt_evaluated_ml_energies = [
        EMT().get_potential_energy(image) for image in final_ml_traj
    ]
    # Compute actual energies for EMT relaxation structures
    emt_relaxation_energies = [
        image.get_potential_energy() for image in parent_calc_traj
    ]
    steps = range(len(final_ml_traj))
    n_samples_iteration = training_params["samples_to_retrain"]
    parent_calls = learner.parent_calls

    def compute_loss(a, b):
        return np.mean(np.sqrt(np.sum((a - b) ** 2, axis=1)))

    initial_structure = images[0].positions
    print(f"Number of AL iterations: {al_iterations}")
    print(f"Number of samples/iteration: {n_samples_iteration}")
    print(f"Total # of queries (parent calls): {parent_calls}\n")

    print(f"Final AL Relaxed Energy: {ml_relaxation_energies[-1]}")
    print(
        f"EMT evaluation at AL structure: {EMT().get_potential_energy(final_ml_traj[-1])}\n"
    )
    al_relaxed_structure = final_ml_traj[-1].positions

    print(f"Total number of EMT steps: {len(emt_relaxation_energies)}")
    print(f"Final EMT Relaxed Energy: {emt_relaxation_energies[-1]}\n")
    emt_relaxed_structure = parent_calc_traj[-1].positions

    initial_structure_error = compute_loss(initial_structure, emt_relaxed_structure)
    relaxed_structure_error = compute_loss(al_relaxed_structure, emt_relaxed_structure)

    print(f"Initial structure error: {initial_structure_error}")
    print(f"AL relaxed structure error: {relaxed_structure_error}")
