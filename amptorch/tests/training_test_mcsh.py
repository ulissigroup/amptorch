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

sigmas = np.logspace(np.log10(0.02), np.log10(1.0), num=4)
MCSHs = {
    "MCSHs": {
        "0": {"groups": [1], "sigmas": sigmas},
        "1": {"groups": [1], "sigmas": sigmas},
        "2": {"groups": [1, 2], "sigmas": sigmas},
        "3": {"groups": [1, 2, 3], "sigmas": sigmas},
    },
    "atom_gaussians": {
        "C": "mcsh_params/C_totaldensity_4.g",
        "O": "mcsh_params/O_totaldensity_5.g",
        "Cu": "mcsh_params/Cu_totaldensity_6.g",
    },
    "cutoff": 10,
}


elements = ["Cu", "C", "O"]
config = {
    "model": {"get_forces": True, "num_layers": 3, "num_nodes": 20},
    "optim": {
        "force_coefficient": 0.04,
        "lr": 1e-3,
        "batch_size": 10,
        "epochs": 300,
        "loss": "mse",
        "metric": "mae",
    },
    "dataset": {
        "raw_data": images,
        "val_split": 0,
        "elements": elements,
        "fp_scheme": "mcsh",
        "fp_params": MCSHs,
        "save_fps": False,
    },
    "cmd": {
        "debug": False,
        "run_dir": "./",
        "seed": 1,
        "identifier": "test",
        "verbose": False,
        "logger": False,
    },
}

true_energies = np.array([image.get_potential_energy() for image in images])
true_forces = np.concatenate(np.array([image.get_forces() for image in images]))


def get_energy_metrics(config):
    trainer = AtomsTrainer(config)
    trainer.train()
    predictions = trainer.predict(images)
    pred_energies = np.array(predictions["energy"])
    mae = np.mean(np.abs(true_energies - pred_energies))
    assert mae < 0.02


def get_force_metrics(config):
    trainer = AtomsTrainer(config)
    trainer.train()
    predictions = trainer.predict(images)
    pred_energies = np.array(predictions["energy"])
    pred_forces = np.concatenate(np.array(predictions["forces"]))

    e_mae = np.mean(np.abs(true_energies - pred_energies))
    f_mae = np.mean(np.abs(pred_forces - true_forces))
    assert e_mae < 0.06
    assert f_mae < 0.06


def test_training():
    torch.set_num_threads(1)

    ### train only
    # energy+forces+mse loss
    config["model"]["get_forces"] = True
    config["optim"]["force_coefficient"] = 0.04
    config["optim"]["loss"] = "mse"
    get_force_metrics(config)
    print("Train energy+forces success!")
    # energy+mae loss
    config["model"]["get_forces"] = False
    config["optim"]["force_coefficient"] = 0
    config["optim"]["loss"] = "mae"
    get_energy_metrics(config)
    print("Train energy only success!")

    ### train+val
    # energy only
    config["model"]["get_forces"] = False
    config["optim"]["force_coefficient"] = 0
    config["optim"]["loss"] = "mae"
    config["dataset"]["val_split"] = 0.1
    get_energy_metrics(config)
    print("Val energy only success!")

    # energy+forces
    config["model"]["get_forces"] = True
    config["optim"]["force_coefficient"] = 0.04
    config["optim"]["loss"] = "mse"
    config["dataset"]["val_split"] = 0.1
    get_force_metrics(config)
    print("Val energy+forces success!")


if __name__ == "__main__":
    test_training()
