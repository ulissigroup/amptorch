import numpy as np
import random

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
    # change train_split to 0.
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
        al_convergence = {"method": "iter", "num_iterations": 3},
        samples_to_retrain=5,
        Gs=Gs,
        morse=True,
        morse_params=morse_params,
        forcetraining=True,
        cores=10,
        criterion=CustomMSELoss,
        num_layers=3,
        num_nodes=20,
        force_coefficient=0.04,
        learning_rate=1e-1,
        epochs=100,
        train_split=0,
    )


    # Calculate true relaxation
    true_relax = Relaxation(slab, BFGS)
    true_relax.run(EMT(), 'true_relax')
    emt_relax = true_relax.get_trajectory('true_relax', 0, -1, 1)
