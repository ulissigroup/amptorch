"""
This script tests whether the constructed sparse matrix is accurately created.

"""

import sys
import torch
import torch.nn as nn
from torch.utils.data import DataLoader
import numpy as np
from ase import Atoms
from ase.calculators.emt import EMT
from amp import Amp
from amp.descriptor.gaussian import Gaussian
from amp_pytorch.data_preprocess import AtomsDataset, factorize_data, collate_amp


def test_matrix_construction():

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

    descriptor = Gaussian(cutoff=6.5, Gs=Gs, fortran=False)
    dataset = AtomsDataset(images, descriptor)
    batch_size = len(dataset)
    dataloader = DataLoader(dataset, batch_size, collate_fn=collate_amp, shuffle=False)
    for batch in dataloader:
        sparse_matrix = batch[-2].to_dense()

    idx_shift = 0
    fp_length = 3
    atoms_in_batch = 14
    for idx, sample in enumerate(dataset):
        num_atoms = len(sample[0])
        fingerprint_derivatives = sample[-2]
        for key in list(fingerprint_derivatives.keys()):
            base_atom = key[2] + idx_shift
            wrt_atom = key[0] + idx_shift
            coord = key[4]
            actual_fprime = torch.tensor(fingerprint_derivatives[key])
            constructed_fprime = sparse_matrix[
                base_atom * fp_length:base_atom * fp_length + fp_length,
                wrt_atom + atoms_in_batch * coord,
            ]
            assert (
                actual_fprime.tolist() == constructed_fprime.tolist()
            ), "The derivatives of image %i do not match up!" % (idx + 1)
        idx_shift += num_atoms


if __name__ == "__main__":
    test_matrix_construction()
