import numpy as np
import os
from ase import Atoms
from ase.calculators.emt import EMT
from ase.md.verlet import VelocityVerlet
import torch
from amptorch.trainer import AtomsTrainer

### Construct/Load system into a list of ASE atoms object
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


### Hyperparameters that needs to be defined
elements = ["Cu", "C", "O"]
path_to_psp = "<path>/pseudodensity_psp/"
# path to the GMP pseudopotential (.g)files
# please copy the "pseudodensity_psp" folder to somehere and edit the path to it here

nsigmas = 5
max_MCSH_order = 3
max_radial_sigma = 2.0


### Construct GMP configuration, no need to change once the hyperparameters are specified.
sigmas = np.linspace(0, max_radial_sigma, nsigmas + 1, endpoint=True)[1:]
GMPs = {
    "MCSHs": {"orders": list(range(max_MCSH_order + 1)), "sigmas": sigmas},
    "atom_gaussians": {
        x: os.path.join(path_to_psp, "{}_pseudodensity.g".format(x)) for x in elements
    },
    # "cutoff": 10.0, # don't need to specify, a default value will be provided.
    "square": False,
    "solid_harmonics": True,
}

config = {
    "model": {
        "name": "singlenn",
        "get_forces": False,
        "num_layers": 3,
        "num_nodes": 20,
        # "hidden_layers": [20,20,20], # more flexible way of defining NN
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
        "elements": elements,
        "fp_scheme": "gmpordernorm",
        "fp_params": GMPs,
        "save_fps": True,
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

trainer = AtomsTrainer(config)
trainer.train()

predictions = trainer.predict(images[:10])

true_energies = np.array([image.get_potential_energy() for image in images])
pred_energies = np.array(predictions["energy"])

print("Energy MSE:", np.mean((true_energies - pred_energies) ** 2))
