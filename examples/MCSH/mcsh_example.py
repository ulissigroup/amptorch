import numpy as np
from ase import Atoms
from ase.calculators.emt import EMT

from amptorch.trainer import AtomsTrainer

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


sigmas = np.logspace(np.log10(0.02), np.log10(1.0), num=5)
MCSHs = {
    "MCSHs": {
        "0": {"groups": [1], "sigmas": sigmas},
        "1": {"groups": [1], "sigmas": sigmas},
        "2": {"groups": [1, 2], "sigmas": sigmas},
        "3": {"groups": [1, 2, 3], "sigmas": sigmas},
        "4": {"groups": [1, 2, 3, 4], "sigmas": sigmas},
        "5": {"groups": [1, 2, 3, 4, 5], "sigmas": sigmas},
        "6": {"groups": [1, 2, 3, 4, 5, 6, 7], "sigmas": sigmas},
    },
    "atom_gaussians": {
        "C": "./MCSH_potential/C_totaldensity_4.g",
        "O": "./MCSH_potential/O_totaldensity_5.g",
        "Cu": "./MCSH_potential/Cu_totaldensity_6.g",
    },
    "cutoff": 8,
}


elements = ["Cu", "C", "O"]
config = {
    "model": {"get_forces": False, "num_layers": 3, "num_nodes": 20},
    "optim": {
        "device": "cpu",
        "force_coefficient": 0.0,
        "lr": 1e-2,
        "batch_size": 10,
        "epochs": 100,
    },
    "dataset": {
        "raw_data": images,
        "val_split": 0,
        "elements": elements,
        "fp_scheme": "gmp",
        "fp_params": MCSHs,
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
