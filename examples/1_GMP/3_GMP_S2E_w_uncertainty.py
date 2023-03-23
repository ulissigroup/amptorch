import numpy as np

from ase import Atoms
from ase.calculators.emt import EMT
from amptorch.trainer import AtomsTrainer
from amptorch.uncertainty import ConformalPredictionLatentSpace
from amptorch.uncertainty.utils import calc_uncertainty_metrics


# # -----------------------------------------------------------------
# This example demonstrates the workflow to use Gaussian Multiple fingerprints
#      with SingleNN atomistic neural network architecture
#      on S2E (structure2energy) task with uncertainty via CP
# # -----------------------------------------------------------------


# # -----------------------------------------------------------------
### Step 1:
###     Construct/Load system into a list of ASE atoms object
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


# # -----------------------------------------------------------------
### Step 2:
### Hyperparameters for fingerprints that needs to be defined

sigmas = [0.02, 0.2, 0.4, 0.69, 1.1, 1.66, 2.66]
max_MCSH_order = 3  # order of angular probes

GMP = {
    "MCSHs": {"orders": list(range(max_MCSH_order + 1)), "sigmas": sigmas},
}


# # -----------------------------------------------------------------
### Step 3:
### Hyperparameters for neural network and optimizers

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

# # -----------------------------------------------------------------
### Step 4:
### Training

trainer = AtomsTrainer(config)
trainer.train()

# # -----------------------------------------------------------------
### Step 5:
### Prediction

# uncertainty prediction with CP
uncertainty_model = ConformalPredictionLatentSpace()
res = uncertainty_model.fit_predict(trainer, train_list, full_test_list)
print(res["residuals"])
print(res["uncertainty"])

# obtain metrics on uncertainty prediction
calibration, sharpness = calc_uncertainty_metrics(res["residuals"], res["uncertainty"])
