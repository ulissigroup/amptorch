import numpy as np
from amptorch.descriptor.Gaussian.automated import GaussianDescriptorSet
from amptorch.trainer import AtomsTrainer
import torch
from ase.io import read


elements = ["Pt", "Ag"]  # arbitrary

cutoff = 6.0

# default descriptors
Gs_default = {
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

# automated descriptors (hardcoded )
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
etas = [23.258, 39.395, 8.552, 2.613]  # [5.027, 6.143, 3.182, 1.757]
rss = [2.390, 2.627, 2.908, 4.009]  # [2.326, 2.628, 2.975, 3.930]
gds.batch_add_descriptors(2, etas, rss, [])

torch.set_num_threads(1)


def get_config(descriptors, num_layers=2, num_nodes=10, epochs=100, verbose=False):
    configs = {
        "model": {"get_forces": True, "num_layers": num_layers, "num_nodes": num_nodes},
        "optim": {
            "device": "cpu",
            "force_coefficient": 0.04,
            "lr": 1e-2,
            "batch_size": 10,
            "epochs": epochs,
            "loss": "mse",
            "metric": "mae",
        },
        "dataset": {
            "raw_data": images,
            "val_split": 0.1,
            "elements": elements,
            "fp_params": descriptors,
            "cutoff_params": cosine_cutoff_params,
            "save_fps": True,
        },
        "cmd": {
            "debug": False,
            "run_dir": "./",
            "seed": 1,
            "identifier": "test",
            "verbose": verbose,
            "logger": False,
        },
    }
    return configs


def test_model(configuration):
    trainer = AtomsTrainer(configuration)
    training_time = trainer.train()

    predictions = trainer.predict(images)

    true_energies = np.array([image.get_potential_energy() for image in images])
    pred_energies = np.array(predictions["energy"])

    mse = np.mean((true_energies - pred_energies) ** 2)
    mae = np.mean(np.abs(true_energies - pred_energies))
    print("Energy MSE:", mse)
    print("Energy MAE:", mae)
    return training_time, mse, mae


images100 = read("pt3ag4_100_images.traj", index=":")
images250 = read("pt3ag4_250_images.traj", index=":")

short_training, long_training = 100, 500

long_arch = 3, 5
wide_arch = 2, 10


trials = [
    (images100, short_training, long_arch),
    (images100, short_training, wide_arch),
    (images100, long_training, wide_arch),
    (images250, short_training, wide_arch),
    (images250, long_training, wide_arch),
]
trials_results = []
for i, (images, epochs, (num_layers, num_nodes)) in enumerate(trials):
    results = []
    print("trial (%d/%d) - %d images, % training epochs, %d layers, %d nodes") % (
        i + 1,
        len(trials),
        len(images),
        epochs,
        num_layers,
        num_nodes,
    )
    print("default config")
    config_default = get_config(Gs_default)
    default_results = list(test_model(config_default))
    results.append(default_results)
    print("automated config")
    config_automated = get_config(gds)
    automated_results = list(test_model(config_automated))
    results.append(automated_results)
    trials_results.append(results)
    print()

trials_results = np.array(trials_results)
np.savetxt("trials_results.csv", trials_results, delimiter=",")
print("benchmarking complete")
