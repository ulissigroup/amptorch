import numpy as np
from amptorch.descriptor.Gaussian.automated import GaussianDescriptorSet
from amptorch.trainer import AtomsTrainer
import torch
from ase.io import read

distances = np.linspace(2, 5, 100)
images = read("pt3ag4_100_images.traj", index=":")

elements = ["Pt", "Ag"]  # arbitrary

cutoff = 6.0

Gs = {
    "default": {
        "G2": {
            "etas": np.logspace(np.log10(0.05), np.log10(5.0), num=4),
            "rs_s": [0],
        },
        "G4": {"etas": [0.005], "zetas": [1.0, 4.0], "gammas": [1.0, -1.0]},
        "cutoff": cutoff,
    },
}

cosine_cutoff_params = {"cutoff_func": "Cosine"}

config = {
    "model": {"get_forces": True, "num_layers": 3, "num_nodes": 5},
    "optim": {
        "device": "cpu",
        "force_coefficient": 0.04,
        "lr": 1e-2,
        "batch_size": 10,
        "epochs": 250,
        "loss": "mse",
        "metric": "mae",
    },
    "dataset": {
        "raw_data": images,
        "val_split": 0.1,
        "elements": elements,
        "fp_params": Gs,
        "cutoff_params": cosine_cutoff_params,
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

torch.set_num_threads(1)


def test_model(configuration):
    trainer = AtomsTrainer(configuration)
    trainer.train()

    predictions = trainer.predict(images)

    true_energies = np.array([image.get_potential_energy() for image in images])
    pred_energies = np.array(predictions["energy"])

    print("Energy MSE:", np.mean((true_energies - pred_energies) ** 2))
    print("Energy MAE:", np.mean(np.abs(true_energies - pred_energies)))


print("default config")
test_model(config)

print()
gds = GaussianDescriptorSet(
    elements, cutoff=cutoff, cutoff_params={"cutoff_func": "cosine"}
)
Gs_gds = {
    "default": {
        "G4": {"etas": [0.005], "zetas": [1.0, 4.0], "gammas": [1.0, -1.0]},
        "cutoff": cutoff,
    },
}
gds.process_combinatorial_Gs(Gs_gds)
etas = [5.027, 6.143, 3.182, 1.757]
rss = [2.326, 2.628, 2.975, 3.930]
gds.batch_add_descriptors(2, etas, rss, [])
config["dataset"]["fp_params"] = gds
print("automated G2s")
test_model(config)
