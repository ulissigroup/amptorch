import sys
import ase
from amptorch.utils import make_amp_descriptors_simple_nn
import os
from amp.descriptor.gaussian import Gaussian, make_symmetry_functions
from amp.descriptor.cutoffs import Cosine
from ase.calculators.singlepoint import SinglePointCalculator as sp
from ase.build import molecule
import numpy as np
from pickle import load
from amp.utilities import hash_images as stock_hash
from amptorch.utils import hash_images as new_hash


def test_fp_match():
    for i in range(100):
        # Tests whether the generated fingerprints are consistent with that of AMPs
        atoms = molecule("H2O")
        atoms.set_cell([10, 10, 10])
        atoms.set_pbc = [True] * 3

        atoms.set_calculator(
            sp(atoms=atoms, energy=-1, forces=np.array([[-1, -1, -1], [-1, -1, -1]]))
        )

        Gs = {}
        images = [atoms]
        Gs["G2_etas"] = [0.005] * 2
        Gs["G2_rs_s"] = [0] * 2
        Gs["G4_etas"] = [0.005] * 2
        Gs["G4_zetas"] = [1.0, 4.0]
        Gs["G4_gammas"] = [1.0, -1.0]
        Gs["cutoff"] = 6.5

        elements = ["O", "H"]
        G = make_symmetry_functions(elements=elements, type="G2", etas=Gs["G2_etas"])
        G += make_symmetry_functions(
            elements=elements,
            type="G4",
            etas=Gs["G4_etas"],
            zetas=Gs["G4_zetas"],
            gammas=Gs["G4_gammas"],
        )
        G = {"O": G, "H": G}

        hashes = stock_hash(images)
        amp_hash = list(hashes.keys())[0]

        make_amp_descriptors_simple_nn(images, Gs, elements=elements)
        s_nn_hash = list(new_hash(images, Gs).keys())[0]

        with open("amp-data-fingerprints.ampdb/loose/" + s_nn_hash, "rb") as f:
            simple_nn = load(f)
        os.system("rm amp-data-fingerprints.ampdb/loose/" + s_nn_hash)

        descriptor = Gaussian(elements=elements, Gs=G, cutoff=Cosine(Gs["cutoff"]))
        descriptor.calculate_fingerprints(hashes, calculate_derivatives=True)
        with open("amp-data-fingerprints.ampdb/loose/" + amp_hash, "rb") as f:
            amp = load(f)
        os.system("rm amp-data-fingerprints.ampdb/loose/" + amp_hash)

        for s, am in zip(simple_nn, amp):
            for i, j in zip(s[1], am[1]):
                assert abs(i - j) <= 1e-5, "Fingerprints do not match!"


if __name__ == "__main__":
    test_fp_match()
