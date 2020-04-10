import multiprocessing as mp
import numpy as np

import ase
from ase.calculators.emt import EMT
from ase.build import fcc100, add_adsorbate, molecule
from ase.constraints import FixAtoms
from ase.optimize import BFGS, QuasiNewton

import torch

from amptorch.model import CustomMSELoss
from amptorch.active_learning.generator_funcs  import MDsimulate, Relaxation
from amptorch.active_learning.oal_calc  import AMPOnlineCalc



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
        "uncertain_tol": 0.1,
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
        "test_split": 0,
        "shuffle": False,
        "filename": "oal_test",
        "verbose": 0
    }

    structure_optim = Relaxation(slab, BFGS, fmax=0.05, steps=None)
    online_calc = AMPOnlineCalc(parent_dataset=images, parent_calc=EMT(),
            n_ensembles=3, n_cores='max', training_params=training_params)
    structure_optim.run(online_calc, filename='relax_oal')

    # Calculate true relaxation
    true_relax = Relaxation(slab, BFGS)
    true_relax.run(EMT(), 'true_relax')
    parent_calc_traj = true_relax.get_trajectory('true_relax', 0, -1, 1)
    n_parent_calls = online_calc.parent_calls
    final_oal_traj = ase.io.read("./relax_oal.traj", ":")

    #Compute ML predicted energies
    ml_relaxation_energies = [image.get_potential_energy() for image in final_oal_traj]
    #Compute actual (EMT) energies for ML predicted structures
    emt_evaluated_ml_energies = [EMT().get_potential_energy(image) for image in final_oal_traj]
    #Compute actual energies for EMT relaxation structures
    emt_relaxation_energies = [image.get_potential_energy() for image in parent_calc_traj]
    steps = range(len(final_oal_traj))

    def compute_loss(a, b):
      return np.mean(np.sqrt(np.sum((a - b)**2, axis=1)))

    initial_structure = images[0].positions
    print(f'Number of OAL steps: {len(final_oal_traj)}\nTotal # of queries (EMT calls): {n_parent_calls} \n')
    print(f"Final OAL Relaxed Energy: {ml_relaxation_energies[-1]}")
    print(f'EMT evaluation at OAL structure: {EMT().get_potential_energy(final_oal_traj[-1])}\n')
    oal_relaxed_structure = final_oal_traj[-1].positions

    print(f'Total number of EMT steps: {len(emt_relaxation_energies)}')
    print(f'Final EMT Relaxed Energy: {emt_relaxation_energies[-1]}\n')
    emt_relaxed_structure = parent_calc_traj[-1].positions


    initial_structure_error = compute_loss(initial_structure, emt_relaxed_structure)
    relaxed_structure_error = compute_loss(oal_relaxed_structure, emt_relaxed_structure)

    print(f'Initial structure error: {initial_structure_error}')
    print(f'OAL relaxed structure error: {relaxed_structure_error}')
