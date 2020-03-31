import numpy as np

import ase
from ase.calculators.emt import EMT

from amptorch.model import CustomMSELoss
from amptorch.active_learning.generator_funcs import MDsimulate
from amptorch.active_learning.al_calc import AtomisticActiveLearning


__author__ = "Muhammed Shuaibi"
__email__ = "mshuaibi@andrew.cmu.edu"


if __name__ == "__main__":
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
        "O": {"re": 1.09, "D": 8.575, "sig": 0.603},
        "Cu": {"re": 2.168, "D": 3.8386, "sig": 1.696},
    }

    # Define initial set of images, can be as few as 1. If 1, make sure to
    # change train_split to 0.
    # images = ase.io.read("../../datasets/COCu_ber_50ps_300K.traj", ":2000:10")
    images = ase.io.read("../../datasets/COCu_DFT_10ps.traj", ":100")

    # Define AL calculator
    al_calc = AtomisticActiveLearning(
        parent_calc=EMT(),
        images=images,
        filename="COCu_ml_2ps_300K_noval",
        file_dir="./morse/dft/",
    )

    # Define AL generating function and training scheme.
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
        Gs=Gs,
        morse=False,
        morse_params=morse_params,
        forcetraining=True,
        cores=10,
        criterion=CustomMSELoss,
        num_layers=3,
        num_nodes=20,
        force_coefficient=0.04,
        learning_rate=1e-1,
        epochs=200,
        train_split=0,
    )
