import ase
import numpy as np
from ase import Atoms
from ase.calculators.emt import EMT
from amp_pytorch.neighbors import get_distances



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

cutoff = 3
atoms = images[0]
energy = atoms.get_potential_energy()

pairs, distances = get_distances(atoms, cutoff)
