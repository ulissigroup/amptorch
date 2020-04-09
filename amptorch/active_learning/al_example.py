import numpy as np
import random
import torch

import ase
from ase.calculators.calculator import Calculator
from ase.calculators.emt import EMT
from ase.calculators.singlepoint import SinglePointCalculator as sp
from ase.build import fcc100, add_adsorbate, molecule
from ase.constraints import FixAtoms
from ase.optimize import BFGS

from amptorch.active_learning.generator_funcs import MDsimulate, Relaxation
from amptorch.active_learning.al_calc import AtomisticActiveLearning
from amptorch.model import CustomMSELoss

if __name__ == "__main__":
	random.seed(0)
	# Define initial set of images, can be as few as 1. If 1, make sure to
	slab = fcc100("Cu", size=(3, 3, 3))
	ads = molecule("C")
	add_adsorbate(slab, ads, 3, offset=(1, 1))
	cons = FixAtoms(
		indices=[atom.index for atom in slab if (atom.tag == 3)]
	)
	slab.set_constraint(cons)
	slab.center(vacuum=13.0, axis=2)
	slab.set_pbc(True)
	slab.wrap(pbc=[True] * 3)
	slab.set_calculator(EMT())

	images = [slab]

	# Define symmetry functions
	Gs = {}
	Gs["G2_etas"] = np.logspace(np.log10(0.05), np.log10(5.0), num=4)
	Gs["G2_rs_s"] = [0] * 4
	Gs["G4_etas"] = [0.005]
	Gs["G4_zetas"] = [1.0, 4.0]
	Gs["G4_gammas"] = [+1.0, -1]
	Gs["cutoff"] = 5.876798323827276  # EMT asap_cutoff: False

	# Define morse parameters if Delta-ML model, o/w morse = False
	morse = True
	morse_params = {
		"C": {"re": 0.972, "D": 6.379, "sig": 0.477},
		"Cu": {"re": 2.168, "D": 3.8386, "sig": 1.696},
	}

	training_params = {
        "al_convergence": {"method": "iter", "num_iterations": 3},
        "samples_to_retrain": 5,
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
        "epochs": 100,
        "test_split": 0.3,
        "shuffle": False
	}

	# Define AL calculator
	al_calc = AtomisticActiveLearning(
		parent_calc=EMT(), images=images, filename="relax_example",
		file_dir="./"
	)
	# Define AL generating function and training scheme.
	al_calc.active_learner(
		generating_function=Relaxation(
			initial_geometry=images[0].copy(),
			optimizer=BFGS,
			fmax=0.05,
			steps=50,
		),
        training_params=training_params
	)


	# Calculate true relaxation
	true_relax = Relaxation(slab, BFGS)
	true_relax.run(EMT(), 'true_relax')
	parent_calc_traj = true_relax.get_trajectory('true_relax', 0, -1, 1)
	final_ml_traj = ase.io.read("./relax_example_iter_3.traj", ":")

	#Compute ML predicted energies
	ml_relaxation_energies = [image.get_potential_energy() for image in final_ml_traj]
	#Compute actual (EMT) energies for ML predicted structures
	emt_evaluated_ml_energies = [EMT().get_potential_energy(image) for image in final_ml_traj]
	#Compute actual energies for EMT relaxation structures
	emt_relaxation_energies = [image.get_potential_energy() for image in parent_calc_traj]
	steps = range(len(final_ml_traj))

	def compute_loss(a, b):
	  return np.mean(np.sqrt(np.sum((a - b)**2, axis=1)))

	initial_structure = images[0].positions
	print('Number of AL iterations: 3\nNumber of samples/iteration: 5\nTotal # of queries (EMT calls): 15 \n')
	print(f'Final AL Relaxed Energy: {ml_relaxation_energies[-1]}')
	print(f'EMT evaluation at AL structure: {EMT().get_potential_energy(final_ml_traj[-1])}\n')
	al_relaxed_structure = final_ml_traj[-1].positions

	print(f'Total number of EMT steps: {len(emt_relaxation_energies)}')
	print(f'Final EMT Relaxed Energy: {emt_relaxation_energies[-1]}\n')
	emt_relaxed_structure = parent_calc_traj[-1].positions


	initial_structure_error = compute_loss(initial_structure, emt_relaxed_structure)
	relaxed_structure_error = compute_loss(al_relaxed_structure, emt_relaxed_structure)

	print(f'Initial structure error: {initial_structure_error}')
	print(f'AL relaxed structure error: {relaxed_structure_error}')
