import numpy as np
from amptorch.descriptor.Gaussian.automated import GaussianDescriptorSet
from amptorch.trainer import AtomsTrainer
import torch
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

elements = ["Cu", "C", "O"]  # arbitrary

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
        "epochs": 10,
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
gds.process_combinatorial_Gs(Gs)
print("gds hash:", gds.descriptor_setup_hash)

print()
print("default config --> GDS")
config["dataset"]["fp_params"] = gds
test_model(config)

print()
print("default config --> +1 G2 GDS")
gds.add_g2("Cu", "O", 4.0, 2.5)
print("gds hash:", gds.descriptor_setup_hash)
test_model(config)

print()
print("default config --> +2 G2 GDS")
gds.add_g2("Cu", "O", 6.0, 2.5)
print("gds hash:", gds.descriptor_setup_hash)
test_model(config)
