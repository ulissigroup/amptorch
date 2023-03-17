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
    "model": {"get_forces": True, "num_layers": 3, "num_nodes": 5},
    "optim": {
        "force_coefficient": 0.04,
        "lr": 1e-3,
        "batch_size": 10,
        "epochs": 5,
        "loss": "mse",
        "metric": "mae",
    },
    "dataset": {
        "raw_data": images,
        "val_split": 0,
        "elements": elements,
        "fp_scheme": "gaussian",
        "fp_params": Gs,
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


def get_metrics(trainer):
    predictions = trainer.predict(images)
    pred_energies = np.array(predictions["energy"])
    pred_forces = np.concatenate(np.array(predictions["forces"]))
    e_mae = np.mean(np.abs(true_energies - pred_energies))
    f_mae = np.mean(np.abs(pred_forces - true_forces))

    return e_mae, f_mae


def test_pretrained():
    torch.set_num_threads(1)
    config_1 = copy.deepcopy(config)
    config_2 = copy.deepcopy(config)

    trainer = AtomsTrainer(config_1)
    trainer.train()
    trained_cpdir = trainer.cp_dir
    e_mae_1, f_mae_1 = get_metrics(trainer)

    config_2["optim"]["epochs"] = 100
    pretrained_trainer = AtomsTrainer(config_2)
    pretrained_trainer.load_pretrained(trained_cpdir)
    e_mae_2, f_mae_2 = get_metrics(pretrained_trainer)

    assert e_mae_1 == e_mae_2, "config - pretrained energy metrics inconsistent!"
    assert f_mae_1 == f_mae_2, "config - pretrained force metrics inconsistent!"

    pretrained_trainer.train()
    e_mae_3, f_mae_3 = get_metrics(pretrained_trainer)
    assert e_mae_3 < e_mae_2, "Retrained metrics are larger!"
    assert f_mae_3 < f_mae_2, "Retrained metrics are larger!"


def test_pretrained_no_config():
    config_1 = copy.deepcopy(config)
    trainer = AtomsTrainer(config_1)
    trainer.train()
    trained_cpdir = trainer.cp_dir
    e_mae_1, f_mae_1 = get_metrics(trainer)

    trainer_2 = AtomsTrainer()
    trainer_2.load_pretrained(trained_cpdir)
    e_mae_2, f_mae_2 = get_metrics(trainer_2)

    assert e_mae_1 == e_mae_2, "configless - pretrained energy metrics inconsistent!"
    assert f_mae_1 == f_mae_2, "configless - pretrained force metrics inconsistent!"


if __name__ == "__main__":
    print("\n\n--------- Pre Trained Test ---------\n")
    test_pretrained()
    test_pretrained_no_config()
