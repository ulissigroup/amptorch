import numpy as np

from amptorch.gaussian import SNN_Gaussian
from amptorch.dataset import (
    AtomsDataset,
    collate_amp
)
from amptorch.model_geometric import BPNN

from ase import Atoms
from ase.calculators.emt import EMT
from ase.io import read

# from torch_geometric.data import DataLoader
from torch.utils.data import DataLoader


# Generate sample dataset
distances = np.linspace(2, 5, 10)
images = []
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

# define symmetry functions to be used
Gs = {}
Gs["G2_etas"] = np.logspace(np.log10(0.05), np.log10(5.0), num=4)
Gs["G2_rs_s"] = [0] * 4
Gs["G4_etas"] = [0.005]
Gs["G4_zetas"] = [1.0, 4.0]
Gs["G4_gammas"] = [+1.0, -1]
Gs["cutoff"] = 6


forcetraining = True
dataset = AtomsDataset(
    images,
    SNN_Gaussian,
    Gs,
    forcetraining=forcetraining,
    label="test",
    cores=4,
)

dataloader = DataLoader(dataset, batch_size=2, shuffle=False, collate_fn=collate_amp)
k = next(iter(dataloader))

unique_atoms = [29, 6, 8]
architecture = [36, 3, 5]
device = "cpu"
forcetraining = True

model = BPNN(unique_atoms, architecture, device, forcetraining)
model(k)
