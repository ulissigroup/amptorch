import matplotlib.pyplot as plt
import ase
import ase.io
from ase.build import molecule
from ase import Atoms
import numpy as np
import torch
from amptorch.ase_utils import AmpTorch
from amptorch.trainer import AtomsTrainer

#  turn off interactive display
plt.ioff()

# Step 1:
# read all images from the trajectory
images = ase.io.read(
    "./water_2d.traj", index=":"
)  # no train-test/holdout split for demo purpose

# Step 2:
# define hyperparameters for GMP fingerprinting
nsigmas = 10  # number of radial probes
max_MCSH_order = 3  # order of angular probes
max_radial_sigma = 2.0  # the maximal sigma of gaussian in radial coordiantes

sigmas = np.linspace(0, max_radial_sigma, nsigmas + 1, endpoint=True)[1:]
GMPs = {
    "MCSHs": {"orders": list(range(max_MCSH_order + 1)), "sigmas": sigmas},
    # "atom_gaussians": {
    #     x: os.path.join(path_to_psp, "{}_pseudodensity.g".format(x)) for x in elements
    # },              # Optional, fitted pseudodensity for valence would be used by default.
    # "cutoff": 10.0, # Optional, a default value will be provided.
}

# Step 3:
# define the configuration for training with S2EF task template
config = {
    "model": {
        "name": "singlenn",
        "get_forces": True,
        # "num_layers": 3,
        # "num_nodes": 20,
        "hidden_layers": [
            20,
            20,
            20,
        ],  # more flexible way of defining NN, alternative to define both "num_layers" and "num_nodes"
        "activation": torch.nn.Tanh,
        "batchnorm": True,
    },
    "optim": {
        "device": "cpu",
        "force_coefficient": 0.01,
        "lr": 1e-3,
        "batch_size": 16,
        "epochs": 500,
        "loss": "mse",
        "metric": "mae",
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
        "identifier": "2d_water",
        "verbose": True,
        "logger": False,
    },
}

# Step 4:
# train
torch.set_num_threads(1)
trainer = AtomsTrainer(config)
trainer.train()

# Step 5:
# assess training errors
predictions = trainer.predict(images)

true_energies = np.array([image.get_potential_energy() for image in images])
pred_energies = np.array(predictions["energy"])

print("Train Energy MSE:", np.mean((true_energies - pred_energies) ** 2))
print("Train Energy MAE:", np.mean(np.abs(true_energies - pred_energies)))

# to get individual energy through an ase.calculator wrapper class
# images[0].set_calculator(AmpTorch(trainer))
# images[0].get_potential_energy()

# Step 6:
# Use the fitted model to predict on an extended range of the 1D PES for change in O-H bond length
# set up images with one changing bond length

# generate NNFF predictions
distances = np.linspace(0.4, 2.0, 100)
images_1d = []
for dist in distances:
    image = molecule("H2O", vacuum=10.0)
    image.set_cell([10, 10, 10])
    image.set_pbc([1, 1, 1])

    # change bond length
    image.set_distance(0, 2, dist)
    image.set_angle(1, 0, 2, 104.210)
    images_1d.append(image)

predictions = trainer.predict(images_1d)

# get training point

training_angle100 = [
    _ for _ in images if np.isclose(_.get_angle(1, 0, 2), 104.210, atol=1e-3)
]

distances_training = [_.get_distance(0, 2) for _ in training_angle100]
energies_training = [_.get_potential_energy() for _ in training_angle100]

# predict on arbitrary O-H length
plt.scatter(distances, predictions["energy"], label="prediction")
plt.scatter(distances_training, energies_training, label="training")
plt.xlabel("O-H bond length [A]")
plt.ylabel("potential energy [eV]")
plt.legend()

# save figure
plt.savefig("predicted_1D_water_PES.png")
plt.close()
