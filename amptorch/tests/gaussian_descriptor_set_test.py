import numpy as np
from amptorch.descriptor.Gaussian import GaussianDescriptorSet
from amptorch.dataset import AtomsDataset

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

forcetraining = True
fp_scheme = "gaussian"
fp_params = Gs
save_fps = True
cutoff_params = {"cutoff_func": "Cosine"}

descriptor_setup = (
    fp_scheme,
    fp_params,
    cutoff_params,
    elements,
)

train_dataset = AtomsDataset(
    images=[],
    descriptor_setup=descriptor_setup,
    forcetraining=forcetraining,
    save_fps=True,
    scaling={"type": "normalize", "range": (0, 1)},
    process=False,
)


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
            matches[params1[0]][0 if np.all(params1 == params2) else 1] += 1
        print("%s descriptors" % element)
        for i in [2, 4, 5]:
            print(
                "\t(%d/%d) G%d descriptors match" % (matches[i][0], sum(matches[i]), i)
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
    gds = GaussianDescriptorSet(elements)
    gds.process_combinatorial_Gs(Gs)

    gaussian_setup = train_dataset.descriptor.descriptor_setup
    gds_setup = gds.descriptor_setup

    compare_setups(gaussian_setup, gds_setup)
    print("Gaussian and GaussianDescriptorSet setups match!")

    gaussian_hash = train_dataset.descriptor.descriptor_setup_hash
    gds_hash = gds.descriptor_setup_hash
    compare_hashes(gaussian_hash, gds_hash)
    print("Gaussian and GaussianDescriptorSet hashes match!")


if __name__ == "__main__":
    print("\n\n--------- Gaussian Descriptor Set Test ---------\n")
    test_gaussian_descriptor_set()
