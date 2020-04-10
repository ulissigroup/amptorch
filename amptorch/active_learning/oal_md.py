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

    md_runner = MDsimulate(ensemble="nvtberendsen", dt=1, temp=300,
            count=2000, initial_geometry=slab)
    online_calc = AMPOnlineCalc(images, EMT(), 5, training_params)
    md_runner.run(online_calc, filename='md_oal')

    # Calculate true relaxation
    md_runner.run(EMT(), 'true_md')
    parent_calc_traj = md_runner.get_trajectory('true_md', 0, -1, 1)
    n_parent_calls = online_calc.parent_calls
    final_oal_traj = ase.io.read("./md_oal.traj", ":")

    #Compute ML predicted energies
    ml_relaxation_energies = [image.get_potential_energy() for image in final_oal_traj]
    #Compute actual (EMT) energies for ML predicted structures
    emt_evaluated_ml_energies = [EMT().get_potential_energy(image) for image in final_oal_traj]
    #Compute actual energies for EMT relaxation structures
    emt_relaxation_energies = [image.get_potential_energy() for image in parent_calc_traj]
