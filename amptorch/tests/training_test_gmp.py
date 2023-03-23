import numpy as np
import torch
from ase import Atoms
from ase.calculators.emt import EMT

from amptorch.trainer import AtomsTrainer

### Construct test data
distances = np.linspace(2, 5, 200)
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
nsigmas = 15
max_MCSH_order = 3
max_radial_sigma = 2.0

sigmas = np.linspace(0, max_radial_sigma, nsigmas + 1, endpoint=True)[1:]
GMP = {
    "MCSHs": {"orders": list(range(max_MCSH_order + 1)), "sigmas": sigmas},
    # "atom_gaussians": {
    #     "C": "amptorch/tests/GMP_params/C_pseudodensity_4.g",
    #     "O": "amptorch/tests/GMP_params/O_pseudodensity_4.g",
    #     "Cu": "amptorch/tests/GMP_params/Cu_pseudodensity_4.g",
    # },
    "cutoff": 8,
}

elements = ["Cu", "C", "O"]


def get_config():
    config = {
        "model": {
            "name": "singlenn",
            "get_forces": True,
            "num_layers": 3,
            "num_nodes": 20,
            "batchnorm": True,
            "activation": torch.nn.GELU,
        },
        "optim": {
            "force_coefficient": 0.10,
            "lr": 1e-3,
            "batch_size": 16,
            "epochs": 100,
            "loss": "mse",
            "metric": "mae",
        },
        "dataset": {
            "raw_data": images,
            "fp_scheme": "gmpordernorm",
            "fp_params": GMP,
            "elements": elements,
            "save_fps": False,
            "scaling": {"type": "normalize", "range": (-1, 1)},
            "val_split": 0,
        },
        "cmd": {
            "debug": False,
            "run_dir": "./",
            "seed": 1,
            "identifier": "test",
            "verbose": False,
        },
    }

    return config


true_energies = np.array([image.get_potential_energy() for image in images])
true_forces = np.concatenate(np.array([image.get_forces() for image in images]))


def get_energy_metrics(config):
    trainer = AtomsTrainer(config)
    trainer.train()
    predictions = trainer.predict(images)
    pred_energies = np.array(predictions["energy"])
    mae = np.mean(np.abs(true_energies - pred_energies))
    print("mae={}".format(mae))
    assert mae < 0.10


def get_force_metrics(config):
    trainer = AtomsTrainer(config)
    trainer.train()
    predictions = trainer.predict(images)
    pred_energies = np.array(predictions["energy"])
    pred_forces = np.concatenate(np.array(predictions["forces"]))

    e_mae = np.mean(np.abs(true_energies - pred_energies))
    f_mae = np.mean(np.abs(pred_forces - true_forces))

    print("e_mae=", e_mae)
    print("f_mae=", f_mae)

    assert e_mae < 0.10
    assert f_mae < 0.50


def test_training_gmp():
    torch.set_num_threads(1)

    ### train only
    # energy+forces+mse loss
    config = get_config()
    config["model"]["get_forces"] = True
    config["optim"]["force_coefficient"] = 0.04
    config["optim"]["loss"] = "mse"
    get_force_metrics(config)
    print("Train energy+forces success!")
    # energy+mae loss
    config = get_config()
    config["model"]["get_forces"] = False
    config["optim"]["force_coefficient"] = 0
    config["optim"]["loss"] = "mae"
    get_energy_metrics(config)
    print("Train energy only success!")

    ### train+val
    # energy only
    config = get_config()
    config["model"]["get_forces"] = False
    config["optim"]["force_coefficient"] = 0
    config["optim"]["loss"] = "mae"
    config["dataset"]["val_split"] = 0.1
    get_energy_metrics(config)
    print("Val energy only success!")

    # energy+forces
    config = get_config()
    config["model"]["get_forces"] = True
    config["optim"]["force_coefficient"] = 0.04
    config["optim"]["loss"] = "mse"
    config["dataset"]["val_split"] = 0.1
    get_force_metrics(config)
    print("Val energy+forces success!")


if __name__ == "__main__":
    print("\n\n--------- GMP Training Test ---------\n")
    test_training_gmp()
