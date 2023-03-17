import numpy as np
import torch
from ase import Atoms
from ase.calculators.emt import EMT
from amptorch import AtomsTrainer, AmpTorch


# adapted from https://gitlab.com/ase/ase/-/blob/master/ase/calculators/test.py#L186-202
def numeric_force(atoms, a, i, d=0.0001):
    """
    Compute numeric force on atom with index a, Cartesian component i,
    with finite step of size d
    """
    p0 = atoms.get_positions()
    p = p0.copy()
    p[a, i] += d
    atoms.set_positions(p, apply_constraint=False)
    eplus = atoms.get_potential_energy()
    p[a, i] -= 2 * d
    atoms.set_positions(p, apply_constraint=False)
    eminus = atoms.get_potential_energy()
    atoms.set_positions(p0, apply_constraint=False)
    return (eminus - eplus) / (2 * d)


# adapted from https://gitlab.com/ase/ase/-/blob/master/ase/calculators/test.py#L186-202
def gradient_test(atoms, indices=None):
    """
    Use numeric_force to compare analytical and numerical forces on atoms

    If indices is None, test is done on all atoms.
    """
    if indices is None:
        indices = range(len(atoms))
    f = atoms.get_forces()[indices]
    print("{0:>16} {1:>20}".format("eps", "max(abs(df))"))
    for eps in np.logspace(-4, -8, 4):
        fn = np.zeros((len(indices), 3))
        for idx, i in enumerate(indices):
            for j in range(3):
                fn[idx, j] = numeric_force(atoms, i, j, eps)
        print("{0:16.12f} {1:20.12f}".format(eps, abs(fn - f).max()))
        assert abs(fn - f).max() <= 2e-5, "Energy/Forces are inconsistent!"
    return f, fn


## Construct test data
def test_energy_force_consistency():
    """Gaussian/Neural non-periodic standard.

    Checks that the answer matches that expected from previous Mathematica
    calculations.
    """

    #: Making the list of non-periodic images
    images = [
        Atoms(
            symbols="PdOPd2",
            pbc=True,
            calculator=EMT(),
            cell=np.array([[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]),
            positions=np.array(
                [[0.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0], [1.0, 0.0, 0.0]]
            ),
        ),
        Atoms(
            symbols="PdOPd2",
            pbc=True,
            calculator=EMT(),
            cell=np.array([[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]),
            positions=np.array(
                [[0.0, 1.0, 0.0], [1.0, 2.0, 1.0], [-1.0, 1.0, 2.0], [1.0, 3.0, 2.0]]
            ),
        ),
        Atoms(
            symbols="PdO",
            pbc=True,
            calculator=EMT(),
            cell=np.array([[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]),
            positions=np.array([[2.0, 1.0, -1.0], [1.0, 2.0, 1.0]]),
        ),
        Atoms(
            symbols="Pd2O",
            pbc=True,
            calculator=EMT(),
            cell=np.array([[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]),
            positions=np.array([[-2.0, -1.0, -1.0], [1.0, 2.0, 1.0], [3.0, 4.0, 4.0]]),
        ),
        Atoms(
            symbols="Cu",
            pbc=True,
            calculator=EMT(),
            cell=np.array([[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]),
            positions=np.array([[0.0, 0.0, 0.0]]),
        ),
    ]

    # Parameters
    Gs = {
        "default": {
            "G2": {
                "etas": [0.2],
                "rs_s": [0],
            },
            "G4": {"etas": [0.4], "zetas": [1], "gammas": [1]},
            "cutoff": 6.5,
        },
    }

    config = {
        "model": {
            "get_forces": True,
            "num_layers": 3,
            "num_nodes": 5,
            "batchnorm": False,
        },
        "optim": {
            "force_coefficient": 0.04,
            "lr": 1e-2,
            "batch_size": 32,
            "epochs": 10,
            "loss": "mse",
            "metric": "mae",
            "gpus": 0,
        },
        "dataset": {
            "raw_data": images,
            "val_split": 0,
            "fp_scheme": "gaussian",
            "fp_params": Gs,
            "save_fps": False,
            # feature scaling to be used - normalize or standardize
            # normalize requires a range to be specified
            "scaling": {"type": "normalize", "range": (-1, 1)},
        },
        "cmd": {
            "debug": False,
            "run_dir": "./",
            "seed": 1,
            "identifier": "test",
            "verbose": True,
            # Weights and Biases used for logging - an account(free) is required
            "logger": False,
            "dtype": torch.DoubleTensor,
        },
    }

    torch.set_num_threads(1)
    trainer = AtomsTrainer(config)
    trainer.load()
    trainer.net.initialize()
    calc = AmpTorch(trainer)
    for image in images:
        image.set_calculator(calc)
        f, fn = gradient_test(image)


if __name__ == "__main__":
    print("\n\n--------- Gaussian Consistency Test ---------\n")
    test_energy_force_consistency()
    print("Success! Energy/Forces are physically consistent: F = -dE/dx")
