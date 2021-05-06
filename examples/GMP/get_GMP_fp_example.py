import numpy as np
from ase import Atoms
from ase.calculators.emt import EMT

from amptorch.dataset import AtomsDataset

distances = np.linspace(2, 5, 10)
images = []
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
    images.append(image)


sigmas = [0.02, 0.2, 0.4, 0.69, 1.1, 1.66, 2.66, 4.4]
MCSHs = {
    "MCSHs": {
        "0": {"groups": [1], "sigmas": sigmas},
        "1": {"groups": [1], "sigmas": sigmas},
        "2": {"groups": [1, 2], "sigmas": sigmas},
        "3": {"groups": [1, 2, 3], "sigmas": sigmas},
        # "4": {"groups": [1, 2, 3, 4], "sigmas": sigmas},
        # "5": {"groups": [1, 2, 3, 4, 5], "sigmas": sigmas},
        # "6": {"groups": [1, 2, 3, 4, 5, 6, 7], "sigmas": sigmas},
    },
    "atom_gaussians": {
        "C": "./valence_gaussians/C_pseudodensity_4.g",
        "O": "./valence_gaussians/O_pseudodensity_4.g",
        "Cu": "./valence_gaussians/Cu_pseudodensity_4.g",
    },
    "cutoff": 8,
}

elements = ["Cu", "C", "O"]
fp_elements = ["Cu", "C", "O"]
ref_elements = ["Cu", "C", "O"]
dataset = AtomsDataset(
    images=images,
    descriptor_setup=("gmp", MCSHs,{}, elements, fp_elements, ref_elements),
    forcetraining=False,
    save_fps=True,
    scaling={"type": "normalize", "range": (0, 1), "elementwise": False}
)

for data in dataset:
    print(data.fingerprint.numpy().shape)
