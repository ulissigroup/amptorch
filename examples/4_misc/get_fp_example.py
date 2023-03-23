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

elements = ["Cu", "C", "O"]
dataset = AtomsDataset(
    images=images,
    descriptor_setup=("gaussian", Gs, {"cutoff_func": "Cosine"}, elements),
    forcetraining=False,
    save_fps=True,
)

for data in dataset:
    print(data.fingerprint.numpy().shape)
