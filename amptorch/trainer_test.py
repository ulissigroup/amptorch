import ase.io
import numpy as np
from ase import Atoms
from ase.calculators.emt import EMT

from amptorch.trainer import AtomsTrainer

distances = np.linspace(2, 5, 10)
images = []
for l in distances:
    image = Atoms(
        "CuCO",
        [
            (-l * np.sin(0.65), l * np.cos(0.65), 0),
            (0, 0, 0),
            (l * np.sin(0.65), l * np.cos(0.65), 0),
        ],
    )
    image.set_cell([10, 10, 10])
    image.wrap(pbc=True)
    image.set_calculator(EMT())
    images.append(image)

images = ase.io.read(
    "/home/mshuaibi/Documents/amptorch/datasets/COCu_ber_50ps_300K.traj", ":1000"
)

Gs = {
    "default": {
        "G2": {"etas": np.logspace(np.log10(0.05), np.log10(5.0), num=4), "rs_s": [0],},
        "G4": {"etas": [0.005], "zetas": [1.0, 4.0], "gammas": [1.0, -1.0]},
        "cutoff": 6,
    },
}

elements = ["Cu", "C", "O"]
config = {
    "model": {"forcetraining": True, "num_layers": 3, "num_nodes": 5,},
    "optim": {
        "device": "cpu",
        "force_coefficient": 0.04,
        "lr": 1e-2,
        # "batch_size": len(images),
        "batch_size": 32,
        "epochs": 1,
    },
    "dataset": {
        "raw_data": images,
        "val_split": 0,
        "elements": elements,
        "fp_params": Gs,
        "save_fps": True,
    },
    "cmd": {
        "debug": False,
        "run_dir": "./",
        "seed": 1,
        "identifier": "test",
        "verbose": True,
    },
}

trainer = AtomsTrainer(config)
trainer.train()

# predicting
test_dataset = {
        "raw_data": images[:2],
        "val_split": 0,
        "elements": elements,
        "fp_params": Gs,
        "save_fps": True,
    }
predicted = trainer.predict(test_dataset)