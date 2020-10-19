import numpy as np
from amptorch.ase_utils import AMPtorch
from amptorch.subsample import subsample_traj
from amptorch.trainer import AtomsTrainer
from ase import Atoms
from ase.calculators.emt import EMT

distances = np.linspace(2, 5, 500)
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

elements = ["Cu", "C", "O"]

# apply subsampling
subsampled_images, _ = subsample_traj(images, Gs, elements, cutoff_sig=0.05, rate=0.3)

print("Length of training set before subsample: {}".format(len(images)))
print("Length of training set after subsample: {}".format(len(subsampled_images)))

config = {
    "model": {"get_forces": False, "num_layers": 3, "num_nodes": 5},
    "optim": {
        "device": "cpu",
        "force_coefficient": 0.04,
        "lr": 1e-2,
        "batch_size": 10,
        "epochs": 100,
    },
    "dataset": {
        "raw_data": subsampled_images,
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
        # "logger": True,
    },
}

trainer = AtomsTrainer(config)
trainer.train()

predictions = trainer.predict(subsampled_images)

true_energies = np.array([image.get_potential_energy() for image in subsampled_images])
pred_energies = np.array(predictions["energy"])

print("Energy MSE:", np.mean((true_energies - pred_energies) ** 2))

image.set_calculator(AMPtorch(trainer))
image.get_potential_energy()
