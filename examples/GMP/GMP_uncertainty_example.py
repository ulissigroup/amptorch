import numpy as np

from ase import Atoms
from ase.calculators.emt import EMT
from amptorch.trainer import AtomsTrainer
from amptorch.uncertainty import ConformalPredictionLatentSpace
from amptorch.uncertainty.utils import calc_uncertainty_metrics

### Construct/Load system into a list of ASE atoms object
distances = np.linspace(2, 5, 2000)
train_list = []
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
    train_list.append(image)

distances = np.linspace(2.2, 3, 1000)
full_test_list = []
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
    full_test_list.append(image)


# set up the trainer
sigmas = [0.02, 0.2, 0.4, 0.69, 1.1, 1.66, 2.66, 4.4]
GMP = {
    "MCSHs": {
        "0": {"groups": [1], "sigmas": sigmas},
        "1": {"groups": [1], "sigmas": sigmas},
        "2": {"groups": [1, 2], "sigmas": sigmas},
        "3": {"groups": [1, 2, 3], "sigmas": sigmas},
        # "4": {"groups": [1, 2, 3, 4], "sigmas": sigmas},
        # "5": {"groups": [1, 2, 3, 4, 5], "sigmas": sigmas},
        # "6": {"groups": [1, 2, 3, 4, 5, 6, 7], "sigmas": sigmas},
    },
    "atom_gaussians": {
        "C": "./valence_gaussians/C_pseudodensity_4.g",
        "O": "./valence_gaussians/O_pseudodensity_4.g",
        "Cu": "./valence_gaussians/Cu_pseudodensity_4.g",
    },
    "cutoff": 8,
}


elements = ["Cu", "C", "O"]
config = {
    "model": {
        "name": "singlenn",
        "get_forces": False,
        "num_layers": 3,
        "num_nodes": 20,
    },
    "optim": {
        "device": "cpu",
        "force_coefficient": 0.0,
        "lr": 1e-2,
        "batch_size": 10,
        "epochs": 100,
    },
    "dataset": {
        "raw_data": train_list,
        "val_split": 0,
        "elements": elements,
        "fp_scheme": "gmpordernorm",
        "fp_params": GMP,
        "save_fps": False,
    },
    "cmd": {
        "debug": False,
        "run_dir": "./",
        "seed": 1,
        "identifier": "test",
        "verbose": True,
        "logger": False,
    },
}

# train
trainer = AtomsTrainer(config)
trainer.train()

# uncertainty prediction with CP
uncertainty_model = ConformalPredictionLatentSpace()
res = uncertainty_model.fit_predict(trainer, train_list, full_test_list)
print(res["residuals"])
print(res["uncertainty"])

# obtain metrics on uncertainty prediction
calibration, sharpness = calc_uncertainty_metrics(res["residuals"], res["uncertainty"])
