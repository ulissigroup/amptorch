import copy
import numpy as np
import torch
from ase import Atoms
from ase.calculators.emt import EMT

from amptorch.trainer import AtomsTrainer

### Construct test data
distances = np.linspace(2, 5, 100)
images = []
for dist in distances:
    image = Atoms(
        "CuCO",
        [
            (-dist * np.sin(0.65), dist * np.cos(0.65), 0),
            (0, 0, 0),
            (dist * np.sin(0.65), dist * np.cos(0.65), 0),
        ],
    )
    image.set_cell([10, 10, 10])
    image.wrap(pbc=True)
    image.set_calculator(EMT())
    images.append(image)

### Construct parameters
Gs = {
    "default": {
        "G2": {
            "etas": np.logspace(np.log10(0.05), np.log10(5.0), num=4),
            "rs_s": [0],
        },
        "G4": {"etas": [0.005], "zetas": [1.0, 4.0], "gammas": [1.0, -1.0]},
        "cutoff": 6,
    },
}
elements = ["Cu", "C", "O"]

config = {
    "model": {"get_forces": True, "num_layers": 3, "num_nodes": 20},
    "optim": {
        "force_coefficient": 0.04,
        "lr": 1e-2,
        "batch_size": 10,
        "epochs": 300,
        "loss": "mse",
        "metric": "mae",
    },
    "dataset": {
        "raw_data": images,
        "val_split": 0,
        "elements": elements,
        "fp_scheme": "gaussian",
        "fp_params": Gs,
        "scaling": {"type": "standardize"},
    },
    "cmd": {
        "debug": False,
        "seed": 1,
        "identifier": "test",
        "verbose": False,
        "logger": False,
    },
}

true_energies = np.array([image.get_potential_energy() for image in images])
true_forces = np.concatenate(np.array([image.get_forces() for image in images]))


def get_performance_metrics(config):
    trainer = AtomsTrainer(config)
    trainer.train()
    predictions = trainer.predict(images)
    pred_energies = np.array(predictions["energy"])
    pred_forces = np.concatenate(np.array(predictions["forces"]))

    e_mae = np.mean(np.abs(true_energies - pred_energies))
    f_mae = np.mean(np.abs(pred_forces - true_forces))
    assert e_mae < 0.01, "%f !< .01" % e_mae
    assert f_mae < 0.03, "%f !< .03" % f_mae

    return e_mae, f_mae


def test_cutoff_funcs():
    torch.set_num_threads(1)

    config["model"]["get_forces"] = True
    config["optim"]["force_coefficient"] = 0.04
    config["optim"]["loss"] = "mse"

    ### Cosine cutoff function

    cosine_cutoff_params = {
        "cutoff_func": "Cosine",
    }

    config_1 = copy.deepcopy(config)
    config_1["dataset"]["cutoff_params"] = cosine_cutoff_params
    print("training model with Cosine cutoff function")
    print("E_MAE: %f, F_MAE: %f" % get_performance_metrics(config_1))

    ### Polynomial cutoff function (gamma = 2.0)

    polynomial_cutoff_params = {
        "cutoff_func": "Polynomial",
        "gamma": 2.0,
    }

    config_2 = copy.deepcopy(config)
    config["dataset"]["cutoff_params"] = polynomial_cutoff_params
    print("training model with Polynomial cutoff function (gamma = 2.0)")
    print("E_MAE: %f, F_MAE: %f" % get_performance_metrics(config_2))


if __name__ == "__main__":
    print("\n\n--------- Gaussian Cutoff Functions Test ---------\n")
    test_cutoff_funcs()
