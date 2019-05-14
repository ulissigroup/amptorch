"""
Exact Gaussian-neural scheme forces and energies of five different non-periodic
configurations and three different periodic configurations have been calculated
in Mathematica, and are given below.  This script checks the values calculated
by the code with and without fortran modules.


FullNN weights must be initialized to 0.5
Output Layer must be Tanh instead of Linear
"""

import sys
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader
import ase
from ase import Atoms
from ase.calculators.emt import EMT
from collections import OrderedDict
from amp import Amp
from amp.descriptor.gaussian import Gaussian
from amp.model.neuralnetwork import NeuralNetwork
from amp.model.neuralnetwork import get_random_scalings 
from amp.utilities import hash_images
from amp.model import calculate_fingerprints_range
from amp_pytorch import core
from amp_pytorch.data_preprocess import AtomsDataset, factorize_data, collate_amp
from amp_pytorch.NN_model import FullNN


def test_non_periodic():
    """Gaussian/Neural non-periodic standard.

    Checks that the answer matches that expected from previous Mathematica
    calculations.
    """

    # Making the list of non-periodic images
    images = [
        Atoms(
            symbols="PdOPd2",
            pbc=np.array([False, False, False], dtype=bool),
            calculator=EMT(),
            cell=np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
            positions=np.array(
                [[0.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0], [1.0, 0.0, 0.0]]
            ),
        ),
        Atoms(
            symbols="PdOPd2",
            pbc=np.array([False, False, False], dtype=bool),
            calculator=EMT(),
            cell=np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
            positions=np.array(
                [[0.0, 1.0, 0.0], [1.0, 2.0, 1.0], [-1.0, 1.0, 2.0], [1.0, 3.0, 2.0]]
            ),
        ),
        Atoms(
            symbols="PdO",
            pbc=np.array([False, False, False], dtype=bool),
            calculator=EMT(),
            cell=np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
            positions=np.array([[2.0, 1.0, -1.0], [1.0, 2.0, 1.0]]),
        ),
        Atoms(
            symbols="Pd2O",
            pbc=np.array([False, False, False], dtype=bool),
            calculator=EMT(),
            cell=np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
            positions=np.array([[-2.0, -1.0, -1.0], [1.0, 2.0, 1.0], [3.0, 4.0, 4.0]]),
        ),
        Atoms(
            symbols="Cu",
            pbc=np.array([False, False, False], dtype=bool),
            calculator=EMT(),
            cell=np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
            positions=np.array([[0.0, 0.0, 0.0]]),
        ),
    ]

    # Parameters
    Gs = {
        "O": [
            {"type": "G2", "element": "Pd", "eta": 0.8},
            {
                "type": "G4",
                "elements": ["Pd", "Pd"],
                "eta": 0.2,
                "gamma": 0.3,
                "zeta": 1,
            },
            {
                "type": "G4",
                "elements": ["O", "Pd"],
                "eta": 0.3,
                "gamma": 0.6,
                "zeta": 0.5,
            },
        ],
        "Pd": [
            {"type": "G2", "element": "Pd", "eta": 0.2},
            {
                "type": "G4",
                "elements": ["Pd", "Pd"],
                "eta": 0.9,
                "gamma": 0.75,
                "zeta": 1.5,
            },
            {
                "type": "G4",
                "elements": ["O", "Pd"],
                "eta": 0.4,
                "gamma": 0.3,
                "zeta": 4,
            },
        ],
        "Cu": [
            {"type": "G2", "element": "Cu", "eta": 0.8},
            {
                "type": "G4",
                "elements": ["Cu", "O"],
                "eta": 0.2,
                "gamma": 0.3,
                "zeta": 1,
            },
            {
                "type": "G4",
                "elements": ["Cu", "Cu"],
                "eta": 0.3,
                "gamma": 0.6,
                "zeta": 0.5,
            },
        ],
    }

    hiddenlayers = {"O": (2,), "Pd": (2,), "Cu": (2,)}

    hashed_images = hash_images(images)
    descriptor = Gaussian(cutoff=6.5, Gs=Gs, fortran=False)
    descriptor.calculate_fingerprints(hashed_images, calculate_derivatives=True)
    fingerprints_range = calculate_fingerprints_range(descriptor, hashed_images)

    weights = OrderedDict(
        [
            (
                "O",
                OrderedDict(
                    [
                        (
                            1,
                            np.matrix([[0.5, 0.5], [0.5, 0.5], [0.5, 0.5], [0.5, 0.5]]),
                        ),
                        (2, np.matrix([[0.5], [0.5], [0.5]])),
                    ]
                ),
            ),
            (
                "Pd",
                OrderedDict(
                    [
                        (
                            1,
                            np.matrix([[0.5, 0.5], [0.5, 0.5], [0.5, 0.5], [0.5, 0.5]]),
                        ),
                        (2, np.matrix([[0.5], [0.5], [0.5]])),
                    ]
                ),
            ),
            (
                "Cu",
                OrderedDict(
                    [
                        (
                            1,
                            np.matrix([[0.5, 0.5], [0.5, 0.5], [0.5, 0.5], [0.5, 0.5]]),
                        ),
                        (2, np.matrix([[0.5], [0.5], [0.5]])),
                    ]
                ),
            ),
        ]
    )

    scalings = OrderedDict(
        [
            ("O", OrderedDict([("intercept", 0), ("slope", 1)])),
            ("Pd", OrderedDict([("intercept", 0), ("slope", 1)])),
            ("Cu", OrderedDict([("intercept", 0), ("slope", 1)])),
        ]
    )

    # Testing pure-python and fortran versions of Gaussian-neural force call
    device = "cpu"
    dataset = AtomsDataset(images, descriptor)
    fp_length = dataset.fp_length()
    unique_atoms = dataset.unique()

    batch_size = len(dataset)
    dataloader = DataLoader(dataset, batch_size, collate_fn=collate_amp, shuffle=False)
    model = FullNN(unique_atoms, [fp_length, 2, 2], device)
    for batch in dataloader:
        input_data = [batch[0], len(batch[1])]
        for element in unique_atoms:
            input_data[0][element][0] = (
                input_data[0][element][0].to(device).requires_grad_(True)
            )
        fp_primes = batch[3]
        energy_pred, force_pred = model(input_data, fp_primes)

    calc = Amp(
        descriptor,
        model=NeuralNetwork(
            hiddenlayers=hiddenlayers,
            weights=weights,
            scalings=scalings,
            activation="tanh",
            fprange=fingerprints_range,
            mode="atom-centered",
            fortran=False,
        ),
    )

    amp_energies = [calc.get_potential_energy(image) for image in images]
    amp_forces = [calc.get_forces(image) for image in images]
    amp_forces = np.concatenate(amp_forces)

    for idx, i in enumerate(amp_energies):
        assert round(i, 4) == round(
            energy_pred.tolist()[idx][0], 4
        ), "The predicted energy of image %i is wrong!" % (idx + 1)
    # for idx, sample in enumerate(amp_forces):
        # for idx_d, value in enumerate(sample):
            # predict = force_pred.tolist()[idx][idx_d]
            # assert abs(value - predict) < 0.00001, (
                # assert round(value, 4)round(force_pred.tolist()[idx][idx_d], 4), (
                # "The predicted force of image % i, direction % i is wrong! Values: %s vs %s"
                # % (idx + 1, idx_d, value, force_pred.tolist()[idx][idx_d])
            # )


def test_periodic():
    """Gaussian/Neural periodic standard.

    Checks that the answer matches that expected from previous Mathematica
    calculations.
    """

    # Making the list of periodic images
    images = [
        Atoms(
            symbols="PdOPd",
            pbc=np.array([True, False, False], dtype=bool),
            calculator=EMT(),
            cell=np.array([[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]]),
            positions=np.array([[0.5, 1.0, 0.5], [1.0, 0.5, 1.0], [1.5, 1.5, 1.5]]),
        ),
        Atoms(
            symbols="PdO",
            pbc=np.array([True, True, False], dtype=bool),
            calculator=EMT(),
            cell=np.array([[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]]),
            positions=np.array([[0.5, 1.0, 0.5], [1.0, 0.5, 1.0]]),
        ),
        Atoms(
            symbols="Cu",
            pbc=np.array([True, True, False], dtype=bool),
            calculator=EMT(),
            cell=np.array([[1.8, 0.0, 0.0], [0.0, 1.8, 0.0], [0.0, 0.0, 1.8]]),
            positions=np.array([[0.0, 0.0, 0.0]]),
        ),
    ]

    # Parameters
    Gs = {
        "O": [
            {"type": "G2", "element": "Pd", "eta": 0.8},
            {
                "type": "G4",
                "elements": ["O", "Pd"],
                "eta": 0.3,
                "gamma": 0.6,
                "zeta": 0.5,
            },
        ],
        "Pd": [
            {"type": "G2", "element": "Pd", "eta": 0.2},
            {
                "type": "G4",
                "elements": ["Pd", "Pd"],
                "eta": 0.9,
                "gamma": 0.75,
                "zeta": 1.5,
            },
        ],
        "Cu": [
            {"type": "G2", "element": "Cu", "eta": 0.8},
            {
                "type": "G4",
                "elements": ["Cu", "Cu"],
                "eta": 0.3,
                "gamma": 0.6,
                "zeta": 0.5,
            },
        ],
    }

    hiddenlayers = {"O": (2,), "Pd": (2,), "Cu": (2,)}

    weights = OrderedDict(
        [
            (
                "O",
                OrderedDict(
                    [
                        (1, np.matrix([[0.5, 0.5], [0.5, 0.5], [0.5, 0.5]])),
                        (2, np.matrix([[0.5], [0.5], [0.5]])),
                    ]
                ),
            ),
            (
                "Pd",
                OrderedDict(
                    [
                        (1, np.matrix([[0.5, 0.5], [0.5, 0.5], [0.5, 0.5]])),
                        (2, np.matrix([[0.5], [0.5], [0.5]])),
                    ]
                ),
            ),
            (
                "Cu",
                OrderedDict(
                    [
                        (1, np.matrix([[0.5, 0.5], [0.5, 0.5], [0.5, 0.5]])),
                        (2, np.matrix([[0.5], [0.5], [0.5]])),
                    ]
                ),
            ),
        ]
    )
    scalings = OrderedDict(
        [
            ("O", OrderedDict([("intercept", 0), ("slope", 1)])),
            ("Pd", OrderedDict([("intercept", 0), ("slope", 1)])),
            ("Cu", OrderedDict([("intercept", 0), ("slope", 1)])),
        ]
    )


    hashed_images = hash_images(images)
    descriptor = Gaussian(cutoff=4, Gs=Gs, fortran=False)
    descriptor.calculate_fingerprints(hashed_images, calculate_derivatives=True)
    fingerprints_range = calculate_fingerprints_range(descriptor, hashed_images)


    # Testing pure-python and fortran versions of Gaussian-neural force call
    device = "cpu"
    dataset = AtomsDataset(images, descriptor)
    fp_length = dataset.fp_length()
    unique_atoms = dataset.unique()
    batch_size = len(dataset)
    dataloader = DataLoader(dataset, batch_size, collate_fn=collate_amp, shuffle=False)
    model = FullNN(unique_atoms, [fp_length, 2, 2], device)
    for batch in dataloader:
        input_data = [batch[0], len(batch[1])]
        for element in unique_atoms:
            input_data[0][element][0] = (
                input_data[0][element][0].to(device).requires_grad_(True)
            )
        fp_primes = batch[3]
        energy_pred, force_pred = model(input_data, fp_primes)

    calc = Amp(
        descriptor,
        model=NeuralNetwork(
            hiddenlayers=hiddenlayers,
            weights=weights,
            scalings=scalings,
            activation="tanh",
            fprange=fingerprints_range,
            mode="atom-centered",
            fortran=False,
        ),
    )

    amp_energies = [calc.get_potential_energy(image) for image in images]
    amp_forces = [calc.get_forces(image) for image in images]
    amp_forces = np.concatenate(amp_forces)

    for idx, i in enumerate(amp_energies):
        assert round(i, 4) == round(
            energy_pred.tolist()[idx][0], 4
        ), "The predicted energy of image %i is wrong!" % (idx + 1)
    for idx, sample in enumerate(amp_forces):
        for idx_d, value in enumerate(sample):
            assert round(value, 4) == round(force_pred.tolist()[idx][idx_d], 4), (
                "The predicted force of image % i, direction % i is wrong!"
                % (idx + 1, idx_d)
            )


if __name__ == "__main__":
    test_non_periodic()
    test_periodic()
