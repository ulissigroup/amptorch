import argparse
import numpy as np
import torch
from ase import Atoms
from ase.calculators.emt import EMT
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
    images.append(image)

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


parser = argparse.ArgumentParser()
parser.add_argument("--num_layers", default=3, type=int)
parser.add_argument("--num_nodes", default=5, type=int)
parser.add_argument("--force_coefficient", default=0.04, type=float)
parser.add_argument("--lr", default=1e-2, type=float)
parser.add_argument("--batch_size", default=32, type=int)
parser.add_argument("--epochs", default=10, type=int)
parser.add_argument("--loss", default="mse", type=str)
parser.add_argument("--metric", default="mae", type=str)
parser.add_argument("--gpus", default=0, type=int)
args = parser.parse_args()

config = {
    "model": {
        "get_forces": True,
        "num_layers": args.num_layers,
        "num_nodes": args.num_nodes,
        "batchnorm": False,
    },
    "optim": {
        "force_coefficient": args.force_coefficient,
        "lr": args.lr,
        "batch_size": args.batch_size,
        "epochs": args.epochs,
        "loss": args.loss,
        "metric": args.metric,
        "gpus": args.gpus,
    },
    "dataset": {
        "raw_data": images,
        "val_split": 0.1,
        "fp_params": Gs,
        "save_fps": True,
        # feature scaling to be used - normalize or standardize
        # normalize requires a range to be specified
        "scaling": {"type": "normalize", "range": (0, 1)},
    },
    "cmd": {
        "debug": False,
        "run_dir": "./",
        "seed": 1,
        "identifier": "test",
        "verbose": True,
        # Weights and Biases used for logging - an account(free) is required
        "logger": True,
    },
}

torch.set_num_threads(1)
trainer = AtomsTrainer(config)
trainer.train()
