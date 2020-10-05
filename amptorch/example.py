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

MCSHs = {   "MCSHs": {   
                        "0": {"groups": [1], "sigmas": [0.1, 0.2, 0.3]},
                        "1": {"groups": [1], "sigmas": [0.1, 0.2, 0.3]},
                        "2": {"groups": [1,2], "sigmas": [0.1, 0.2, 0.3]},
                        "3": {"groups": [1,2,3], "sigmas": [0.1, 0.2, 0.3]},
                        "4": {"groups": [1,2,3,4], "sigmas": [0.1, 0.2, 0.3]},
                        # "5": {"groups": [1,2,3,4,5], "sigmas": [0.1, 0.2, 0.3]},
                        # "6": {"groups": [1,2,3,4,5,6,7], "sigmas": [0.1, 0.2, 0.3]},
                        # "7": {"groups": [1,2,3,4,5,6,7,8], "sigmas": [0.1, 0.2, 0.3]},
                        # "8": {"groups": [1,2,3,4,5,6,7,8,9,10], "sigmas": [0.1, 0.2, 0.3]},
                        # "9": {"groups": [1,2,3,4,5,6,7,8,9,10,11,12], "sigmas": [0.1, 0.2, 0.3]}
                  },
            "atom_gaussians": {
                        "H": "./MCSH_potential/H_pseudodensity_6.g",
                        "O": "./MCSH_potential/O_pseudodensity_6.g",
                        "Fe": "./MCSH_potential/Pt_pseudodensity_8.g"
                  },
            "cutoff": 6
}


elements = ["Cu", "C", "O"]
config = {
    "model": {"get_forces": True, "num_layers": 3, "num_nodes": 5},
    "optim": {
        "device": "cpu",
        "force_coefficient": 0.04,
        "lr": 1e-2,
        "batch_size": 10,
        "epochs": 100,
    },
    "dataset": {
        "raw_data": images,
        "val_split": 0,
        "elements": elements,
        "fp_scheme": "mcsh",
        "fp_params": Gs,
        "save_fps": True,
    },
    "cmd": {
        "debug": False,
        "run_dir": "./",
        "seed": 1,
        "identifier": "test",
        "verbose": True,
        "logger": True,
    },
}

trainer = AtomsTrainer(config)
trainer.train()

predictions = trainer.predict(images[:10])

true_energies = np.array([image.get_potential_energy() for image in images])
pred_energies = np.array(predictions["energy"])

print("Energy MSE:", np.mean((true_energies - pred_energies) ** 2))
