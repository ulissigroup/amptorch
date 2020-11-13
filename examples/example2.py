import numpy as np
from amptorch.descriptor.Gaussian import GaussianDescriptorSet
from amptorch.trainer import AtomsTrainer

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

elements = ["Cu", "C", "O"]  # arbitrary

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
        "raw_data": [],  # no images required to confirm descriptors + hash match
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

cosine_trainer = AtomsTrainer(config)
cosine_trainer.load_config()
cosine_trainer.load_rng_seed()
cosine_trainer.load_dataset(process=False)
gds = GaussianDescriptorSet(cosine_trainer.elements)
gds.process_combinatorial_Gs(Gs)

gaussian_setup = cosine_trainer.train_dataset.descriptor.descriptor_setup
gaussian_hash = cosine_trainer.train_dataset.descriptor.descriptor_setup_hash

gds_setup = gds.descriptor_setup
gds_hash = gds.descriptor_setup_hash

print("\n#########################")
print("The two descriptor construction methods match?")
for key in gaussian_setup:
    matches = 0
    for g_params, gds_params in zip(gaussian_setup[key], gds_setup[key]):
        matches += np.all(g_params == gds_params)
    print("%s: (%d/%d) matches" % (key, matches, len(gaussian_setup[key])))

print("#########################\n")

print("Gaussian hash")
print(gds_hash)
print("GaussianDescriptorSet hash")
print(gaussian_hash)
print("hashes match:", gaussian_hash == gds_hash)
