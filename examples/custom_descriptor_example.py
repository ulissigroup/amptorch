import numpy as np
import torch
from amptorch.ase_utils import AMPtorch
from amptorch.descriptor.Gaussian import GaussianDescriptorSet
from amptorch.trainer import AtomsTrainer
from ase import Atoms
from ase.calculators.emt import EMT

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


elements = np.unique([atom.symbol for atom in images[0]])
cutoff = 6.0
cosine_cutoff_params = {"cutoff_func": "cosine"}

Gs = {
    "default": {
        "G2": {
            "etas": np.logspace(np.log10(0.05), np.log10(5.0), num=4),
            "rs_s": [0],
        },
        "G4": {"etas": [0.005], "zetas": [1.0, 4.0], "gammas": [1.0, -1.0]},
        "cutoff": 6.0,
    },
}

gds_combo = GaussianDescriptorSet(elements, cutoff, cosine_cutoff_params)
gds_combo.process_combinatorial_Gs(
    Gs
)  # a GDS can be constructed directly from the `Gs` dict

gds_custom = GaussianDescriptorSet(elements, cutoff, cosine_cutoff_params)
g2_etas = np.logspace(np.log10(0.05), np.log10(5.0), 4)
g2_rs_s = np.zeros(g2_etas.shape)
gds_custom.batch_add_descriptors(
    2, g2_etas, g2_rs_s, []
)  # or the descriptors can be manually added into the GDS
g4_etas = [0.005] * 4
g4_zetas = [1.0, 4.0, 1.0, 4.0]
g4_gammas = [1.0, 1.0, -1.0, -1.0]
gds_custom.batch_add_descriptors(
    4, g4_etas, g4_zetas, g4_gammas
)  # this opens opportunities for finely-tuned, custom descriptors

print("both GaussianDescriptorSets are equivalent!")
print("combo hash", gds_combo.descriptor_setup_hash)
print("custom hash", gds_custom.descriptor_setup_hash)
print("gds_combo == gds_custom:", gds_combo == gds_custom)

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
        "epochs": 100,
        "loss": "mse",
        "metric": "mae",
        "gpus": 0,
    },
    "dataset": {
        "raw_data": images,
        "val_split": 0.1,
        "fp_params": gds_custom,  # either a GDS or the `Gs` dict can be passed here
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

image.set_calculator(AMPtorch(trainer))
image.get_potential_energy()
