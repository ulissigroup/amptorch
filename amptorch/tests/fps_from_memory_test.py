from amptorch.utils import make_amp_descriptors_simple_nn
import os
from amp.descriptor.gaussian import Gaussian, make_symmetry_functions
from amp.descriptor.cutoffs import Cosine
from ase.calculators.singlepoint import SinglePointCalculator as sp
from ase.build import molecule, fcc100, add_adsorbate
from ase.constraints import FixAtoms
import numpy as np
from pickle import load
from amp.utilities import hash_images as stock_hash
from amptorch.utils import hash_images as new_hash
from amptorch.gaussian import SNN_Gaussian
from amptorch.data_preprocess import AtomsDataset, TestDataset
from ase.calculators.emt import EMT


def test_fps_memory():
    slab = fcc100("Cu", size=(3, 3, 3))
    ads = molecule("CO")
    add_adsorbate(slab, ads, 4, offset=(1, 1))
    cons = FixAtoms(
        indices=[atom.index for atom in slab if (atom.tag == 2 or atom.tag == 3)]
    )
    slab.set_constraint(cons)
    slab.center(vacuum=13.0, axis=2)
    slab.set_pbc(True)
    slab.wrap(pbc=[True] * 3)
    slab.set_calculator(EMT())

    images = [slab]

    Gs = {}
    Gs["G2_etas"] = [2]
    Gs["G2_rs_s"] = [0]
    Gs["G4_etas"] = [0.005]
    Gs["G4_zetas"] = [1.0]
    Gs["G4_gammas"] = [1.0]
    Gs["cutoff"] = 6.5

    elements = np.array([atom.symbol for atoms in images for atom in atoms])
    _, idx = np.unique(elements, return_index=True)
    elements = list(elements[np.sort(idx)])

    G = make_symmetry_functions(elements=elements, type="G2", etas=Gs["G2_etas"])
    G += make_symmetry_functions(
        elements=elements,
        type="G4",
        etas=Gs["G4_etas"],
        zetas=Gs["G4_zetas"],
        gammas=Gs["G4_gammas"],
    )
    G = {"O": G, "C": G, "Cu": G}

    snn_hashes = new_hash(images, Gs=Gs)
    base = AtomsDataset(
        images, SNN_Gaussian, Gs, forcetraining=True, label="test", cores=10
    )

    for idx in range(len(images)):
        s_nn_hash = list(snn_hashes.keys())[idx]
        # SimpleNN
        with open("amp-data-fingerprints.ampdb/loose/" + s_nn_hash, "rb") as f:
            simple_nn = load(f)
        os.system("rm amp-data-fingerprints.ampdb/loose/" + s_nn_hash)

        with open("amp-data-fingerprint-primes.ampdb/loose/" + s_nn_hash, "rb") as f:
            simple_nn_prime = load(f)
        os.system("rm amp-data-fingerprint-primes.ampdb/loose/" + s_nn_hash)

        test = TestDataset(
            images[idx],
            base.elements,
            base.base_descriptor,
            Gs,
            base.fprange,
            "test2",
            cores=2,
            save_fps=False,
        )
        test_fp = test.fps
        test_prime = test.fp_primes

        key = simple_nn_prime.keys()

        for s, am in zip(simple_nn, test_fp):
            for i, j in zip(s[1], am[1]):
                assert abs(i - j) <= 1e-4, "Fingerprints do not match! %s, %s" % (i, j)
        for idx in key:
            for s, am in zip(simple_nn_prime[idx], test_prime[idx]):
                assert abs(s - am) <= 1e-4, "Fingerprint primes do not match!"
