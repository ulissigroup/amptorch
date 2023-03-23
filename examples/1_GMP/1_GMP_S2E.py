import numpy as np
from ase import Atoms
from ase.calculators.emt import EMT
import torch
from amptorch.trainer import AtomsTrainer
from amptorch.ase_utils import AmpTorch

# # -----------------------------------------------------------------
# This example demonstrates the workflow to use Gaussian Multiple fingerprints
#      with SingleNN atomistic neural network architecture
#      on S2E (structure2energy) task
# # -----------------------------------------------------------------


# # -----------------------------------------------------------------
### Step 1:
###     Construct/Load system into a list of ASE atoms object

distances = np.linspace(2, 5, 10)
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

###     Alternatively, load the list of ase.Atoms object
###          with ase.io.read

# # -----------------------------------------------------------------
### Step 2:
### Hyperparameters for fingerprints that needs to be defined

nsigmas = 5  # number of radial probes
max_MCSH_order = 3  # order of angular probes
max_radial_sigma = 2.0  # the maximal sigma of gaussian in radial coordiantes

### Construct GMP configuration, no need to change once the hyperparameters are specified.
sigmas = np.linspace(0, max_radial_sigma, nsigmas + 1, endpoint=True)[1:]
GMPs = {
    "MCSHs": {"orders": list(range(max_MCSH_order + 1)), "sigmas": sigmas},
    # "atom_gaussians": {
    #     x: os.path.join(path_to_psp, "{}_pseudodensity.g".format(x)) for x in elements
    # },              # Optional, fitted pseudodensity for valence would be used by default.
    # "cutoff": 10.0, # Optional, a default value will be provided.
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
        # "hidden_layers": [20,20,20], # more flexible way of defining NN, alternative to define both "num_layers" and "num_nodes"
        "activation": torch.nn.GELU,
        "batchnorm": True,
    },
    "optim": {
        "device": "cpu",
        "force_coefficient": 0.0,
        "lr": 1e-2,
        "batch_size": 16,
        "epochs": 100,
    },
    "dataset": {
        "raw_data": images,
        "val_split": 0,
        "fp_scheme": "gmpordernorm",
        "fp_params": GMPs,
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

predictions = trainer.predict(images[:10])

### Error metrics
true_energies = np.array([image.get_potential_energy() for image in images])
pred_energies = np.array(predictions["energy"])

print("Energy MSE:", np.mean((true_energies - pred_energies) ** 2))

# # -----------------------------------------------------------------
### Step 6:
### Use as a calculator for ase
image.set_calculator(AmpTorch(trainer))
image.get_potential_energy()
