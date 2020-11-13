<<<<<<< HEAD
# TODO
#
#
#
# Gs = {
#     "default": {
#         "G2": {
#             "etas": np.logspace(np.log10(0.05), np.log10(5.0), num=4),
#             "rs_s": [0],
#         },
#         "G4": {"etas": [0.005], "zetas": [1.0, 4.0], "gammas": [1.0, -1.0]},
#         "cutoff": 6,
#     },
# }
#
# pass Gs to GaussianDescriptorSet
# pass Gs to Gaussian
#
# show that GDS.to_descriptor_setup() == Gaussian.descriptor_setup
#
#
#
#
#
=======
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


def compare_setups(setup1, setup2):
    print("comparing descriptor setups...")
    assert (
        setup1.keys() == setup2.keys()
    ), "descriptor setups did not have matching elements ({}, {})".format(
        setup1.keys(), setup2.keys()
    )
    for element in setup1:
        assert len(setup1[element]) == len(setup2[element])
        matches = {2: [0, 0], 4: [0, 0], 5: [0, 0]}
        for params1, params2 in zip(setup1[element], setup2[element]):
            matches[params1[0]][1 if np.all(params1 == params2) else 0] += 1
        print("%s descriptors" % element)
        for i in [2, 4, 5]:
            print(
                "\t(%d/%d) G%d descriptors match" % (i, matches[i][0], sum(matches[i]))
            )
        for i in [2, 4, 5]:
            assert matches[i][1] == 0, "only %d/%d %s G%ds matched" % (
                matches[i][0],
                sum(matches[i]),
                element,
                i,
            )


def compare_hashes(hash1, hash2):
    assert hash1 == hash2, "descriptor hashes are not identical ({}, {})".format(
        hash1, hash2
    )


def test_gaussian_descriptor_set():
    cosine_trainer = AtomsTrainer(config)
    cosine_trainer.load_config()
    cosine_trainer.load_rng_seed()
    cosine_trainer.load_dataset(process=False)
    gds = GaussianDescriptorSet(cosine_trainer.elements)
    gds.process_combinatorial_Gs(Gs)

    gaussian_setup = cosine_trainer.train_dataset.descriptor.descriptor_setup
    for element in gaussian_setup.keys():
        gau_descriptors = gaussian_setup[element]
        gds_descriptors = gds_setup[element]
        for gau_d, gds_d, in zip(gau_descriptors, gds_descriptors):
            print('gau', element, gau_d)
            print('gds', element, gds_d)
    gds_setup = gds.descriptor_setup
    compare_setups(gaussian_setup, gds_setup)
    print("Gaussian and GaussianDescriptorSet setups match!")

    gaussian_hash = cosine_trainer.train_dataset.descriptor.descriptor_setup_hash
    gds_hash = gds.descriptor_setup_hash
    compare_hashes(gaussian_hash, gds_hash)
    print("Gaussian and GaussianDescriptorSet hashes match!")


if __name__ == "__main__":
    test_gaussian_descriptor_set()
>>>>>>> fdbcb0c (more debugging)
