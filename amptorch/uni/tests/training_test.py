"""Tests fingerprints generated in training vs calculation, training loop,
convergence, and plotting"""

import os
from pickle import load
from ase import Atoms
from ase.calculators.emt import EMT
import torch
import torch.optim as optim
import numpy as np
from amptorch.uni.NN_uni import CustomLoss
from amptorch import AMP
from amptorch.uni.core_uni import AMPTorch
from amptorch.analysis import parity_plot
from amptorch.gaussian import Gaussian
from amptorch.uni.data_uni import TestDataset


def test_training():
    distances = np.linspace(2, 5, 100)
    label = "example"
    images = []
    energies = []
    forces = []
    for l in distances:
        image = Atoms(
            "CuCO",
            [
                (-l * np.sin(0.65), l * np.cos(0.65), 0),
                (0, 0, 0),
                (l * np.sin(0.65), l * np.cos(0.65), 0),
            ],
        )
        image.set_cell([10, 10, 10])
        image.wrap(pbc=True)
        image.set_calculator(EMT())
        images.append(image)
        energies.append(image.get_potential_energy())
        forces.append(image.get_forces())

    energies = np.array(energies)
    forces = np.concatenate(np.array(forces))
    Gs = {}
    Gs["G2_etas"] = np.logspace(np.log10(0.05), np.log10(5.0), num=2)
    Gs["G2_rs_s"] = [0] * 2
    Gs["G4_etas"] = [0.005]
    Gs["G4_zetas"] = [1.0]
    Gs["G4_gammas"] = [+1.0, -1]
    Gs["cutoff"] = 6.5

    torch.set_num_threads(1)
    calc = AMP(
        model=AMPTorch(
            images,
            descriptor=Gaussian,
            Gs=Gs,
            force_coefficient=0.3,
            label=label,
            save_logs=True,
        )
    )
    calc.model.device = "cpu"
    calc.model.structure = [2, 2]
    calc.model.val_frac = 0
    calc.model.convergence = {
        "energy": 0.005,
        "force": 0.005,
        "early_stop": False,
        "epochs": 1e10,
    }
    calc.model.loader_params = {"batch_size": None, "shuffle": False, "num_workers": 0}
    calc.model.criterion = CustomLoss
    calc.model.optimizer = optim.LBFGS
    calc.model.lr = 1e-2
    calc.model.fine_tune = None

    calc.train(overwrite=True)
    train_hashes = list(calc.model.training_data.hashed_images.keys())
    train_fp_hashes = {}
    for hash in train_hashes:
        with open("amp-data-fingerprints.ampdb/loose/" + hash, "rb") as f:
            fp = load(f)
            train_fp_hashes[hash] = fp
        os.system("rm amp-data-fingerprints.ampdb/loose/" + hash)
    train_prime_hashes = {}
    for hash in train_hashes:
        with open("amp-data-fingerprint-primes.ampdb/loose/" + hash, "rb") as f:
            prime = load(f)
            train_prime_hashes[hash] = prime
        os.system("rm amp-data-fingerprint-primes.ampdb/loose/" + hash)

    dataset = TestDataset(
        images,
        descriptor=calc.model.descriptor,
        Gs=calc.model.training_data.Gs,
        fprange=calc.model.training_data.fprange,
    )

    test_hashes = list(dataset.hashed_images.keys())
    test_fp_hashes = {}
    for hash in test_hashes:
        with open("amp-data-fingerprints.ampdb/loose/" + hash, "rb") as f:
            fp = load(f)
            test_fp_hashes[hash] = fp
        os.system("rm amp-data-fingerprints.ampdb/loose/" + hash)
    test_prime_hashes = {}
    for hash in test_hashes:
        with open("amp-data-fingerprint-primes.ampdb/loose/" + hash, "rb") as f:
            prime = load(f)
            test_prime_hashes[hash] = prime
        os.system("rm amp-data-fingerprint-primes.ampdb/loose/" + hash)

    # test fingerprints are identical
    for train_hash, test_hash in zip(train_hashes, test_hashes):
        for train_fp, test_fp in zip(
            train_fp_hashes[train_hash], test_fp_hashes[test_hash]
        ):
            for i, j in zip(train_fp[1], test_fp[1]):
                assert abs(i - j) <= 1e-5, "Fingerprints do not match!"

    # test fingerprint primes are identical
    for train_hash, test_hash in zip(train_hashes, test_hashes):
        for train_prime, test_prime in zip(
            list(train_prime_hashes[train_hash].values()),
            list(test_prime_hashes[test_hash].values()),
        ):
            for i, j in zip(train_prime, test_prime):
                assert abs(i - j) <= 1e-5, "Fingerprint primes do not match!"

    num_of_atoms = 3
    calculated_energies = np.array(
        [calc.get_potential_energy(image) for image in images]
    )
    energy_rmse = np.sqrt(
        (((calculated_energies - energies) / num_of_atoms) ** 2).sum() / len(images)
    )
    assert (
        energy_rmse <= calc.model.convergence["energy"]
    ), "Energy training convergence not met!"

    calculated_forces = np.concatenate(
        np.array([calc.get_forces(image) for image in images])
    )
    force_rmse = np.sqrt(
        (((calculated_forces - forces)) ** 2).sum() / (3 * num_of_atoms * len(images))
    )
    assert (
        force_rmse <= calc.model.convergence["force"]
    ), "Force training convergence not met!"

    # test plot creation
    parity_plot(calc, images, data="energy", label=label)
    parity_plot(calc, images, data="forces", label=label)
