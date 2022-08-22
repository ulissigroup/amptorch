import numpy as np

from ase import Atoms
from ase.calculators.emt import EMT
from amptorch.trainer import AtomsTrainer
from amptorch.uncertainty import ConformalPredictionLatentSpace
from amptorch.uncertainty.utils import calc_uncertainty_metrics

### Construct/Load system into a list of ASE atoms object
distances = np.linspace(2, 5, 2000)
train_list = []
for dist in distances:
    image = Atoms(
        "CuCO",
        [
            (-dist * np.sin(0.65), dist * np.cos(0.65), 0),
            (0, 0, 0),
            (dist * np.sin(0.65), dist * np.cos(0.65), 0),
        ],
    )
    image.set_cell([10, 10, 10])
    image.wrap(pbc=True)
    image.set_calculator(EMT())
    train_list.append(image)

distances = np.linspace(2.2, 3, 1000)
full_test_list = []
for dist in distances:
    image = Atoms(
        "CuCO",
        [
            (-dist * np.sin(0.65), dist * np.cos(0.65), 0),
            (0, 0, 0),
            (dist * np.sin(0.65), dist * np.cos(0.65), 0),
        ],
    )
    image.set_cell([10, 10, 10])
    image.wrap(pbc=True)
    image.set_calculator(EMT())
    full_test_list.append(image)

### Construct parameters
nsigmas = 5
max_MCSH_order = 3
max_radial_sigma = 2.0

sigmas = np.linspace(0, max_radial_sigma, nsigmas + 1, endpoint=True)[1:]
GMP = {
    "MCSHs": {"orders": list(range(max_MCSH_order + 1)), "sigmas": sigmas},
    "atom_gaussians": {
        "C": "amptorch/tests/GMP_params/C_pseudodensity_4.g",
        "O": "amptorch/tests/GMP_params/O_pseudodensity_4.g",
        "Cu": "amptorch/tests/GMP_params/Cu_pseudodensity_4.g",
    },
    "cutoff": 8,
}


elements = ["Cu", "C", "O"]

config = {
    "model": {
        "name": "singlenn",
        "get_forces": False,
        "num_layers": 3,
        "num_nodes": 20,
        "batchnorm": True,
    },
    "optim": {
        "device": "cpu",
        "force_coefficient": 0.0,
        "lr": 1e-2,
        "batch_size": 10,
        "epochs": 100,
    },
    "dataset": {
        "raw_data": train_list,
        "val_split": 0,
        "elements": elements,
        "fp_scheme": "gmpordernorm",
        "fp_params": GMP,
        "save_fps": False,
    },
    "cmd": {
        "debug": False,
        "run_dir": "./",
        "seed": 1,
        "identifier": "test",
        "verbose": False,
        "logger": False,
    },
}


def test_cp_uncertainty_calibration():
    # train
    trainer = AtomsTrainer(config)
    trainer.train()
    # uncertainty prediction with CP
    uncertainty_model = ConformalPredictionLatentSpace()
    res = uncertainty_model.fit_predict(trainer, train_list, full_test_list)
    # obtain metrics on uncertainty prediction
    calibration, sharpness = calc_uncertainty_metrics(
        res["residuals"], res["uncertainty"]
    )
    assert np.abs(calibration - res["alpha"]) < 0.05


if __name__ == "__main__":
    print("\n\n--------- GMP Conformal Prediction Calibration Test ---------\n")
    test_cp_uncertainty_calibration()
