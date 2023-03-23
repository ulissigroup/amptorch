import numpy as np
import torch
from ase import Atoms
from ase.calculators.emt import EMT

from amptorch.ase_utils import AmpTorch
from amptorch.trainer import AtomsTrainer

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
    image.get_potential_energy()
    images.append(image)

Gs = {
    "gaussian": {
        "G2": {
            "etas": np.logspace(np.log10(0.05), np.log10(5.0), num=4),
            "rs_s": [0],
        },
        "G4": {"etas": [0.005], "zetas": [1.0, 4.0], "gammas": [1.0, -1.0]},
        "cutoff": 6,
    },
}

config = {
    "model": {
        "get_forces": True,
        "num_layers": 3,
        "num_nodes": 5,
        "batchnorm": False,
    },
    "optim": {
        "force_coefficient": 0.04,
        "lr": 1e-2,
        "batch_size": 32,
        "epochs": 10,
        "loss": "mse",
        "metric": "mae",
        "gpus": 0,
    },
    "dataset": {
        # "lmdb_path": ["./data.lmdb", "./data2.lmdb"],
        "lmdb_path": ["./data.lmdb"],
        "val_split": 0,
        "cache": "full",
    },
    "cmd": {
        "debug": False,
        "run_dir": "./",
        "dtype": torch.DoubleTensor,
        "seed": 1,
        "identifier": "test",
        "verbose": True,
        # Weights and Biases used for logging - an account(free) is required
        "logger": False,
    },
}

torch.set_num_threads(1)
trainer = AtomsTrainer(config)
trainer.train()

predictions = trainer.predict(images)

true_energies = np.array([image.get_potential_energy() for image in images])
pred_energies = np.array(predictions["energy"])

print("Energy MSE:", np.mean((true_energies - pred_energies) ** 2))
print("Energy MAE:", np.mean(np.abs(true_energies - pred_energies)))

image.set_calculator(AmpTorch(trainer))
image.get_potential_energy()
