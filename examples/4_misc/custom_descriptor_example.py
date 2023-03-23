import numpy as np
import torch
from amptorch.ase_utils import AmpTorch
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

# These are the default Gs (used in train_example.py)
# Gs = {
#     "default": {
#         "G2": {
#             "etas": np.logspace(np.log10(0.05), np.log10(5.0), num=4),
#             "rs_s": [0],
#         },
#         "G4": {"etas": [0.005], "zetas": [1.0, 4.0], "gammas": [1.0, -1.0]},
#         "cutoff": 6.0,
#     },
# }

gds = GaussianDescriptorSet(elements, cutoff, cosine_cutoff_params)
# gds.process_combinatorial_Gs(Gs)
# The GDS can be constructed directly from the `Gs` dict
# Or, the GDS can be manually constructed, as shown below

g2_etas = [0.25, 2.5, 0.25, 2.5]
g2_rs_s = [0.0, 0.0, 3.0, 3.0]
# the "zip" of these lists will be turned into 4 G2 descriptors
gds.batch_add_descriptors(2, g2_etas, g2_rs_s, [])

g4_etas = [0.005, 0.005, 0.01, 0.01]
g4_zetas = [1.0, 4.0, 4.0, 16.0]
g4_gammas = [1.0, 1.0, -1.0, -1.0]
# the "zip" of these lists will be turned into 4 G4 descriptors
gds.batch_add_descriptors(4, g4_etas, g4_zetas, g4_gammas)

# this opens opportunities for creating finely-tuned descriptor sets

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
        "fp_scheme": "gaussian",
        "fp_params": gds,  # either a GDS or the `Gs` dict can be passed here
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

image.set_calculator(AmpTorch(trainer))
image.get_potential_energy()
