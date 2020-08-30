import sys
import numpy as np

from amptorch.dataset import AMPTorchDataset, collate_amp
from amptorch.model_geometric import BPNN


from torch.utils.data import DataLoader
from amptorch.descriptor.Gaussian import Gaussian
# from amptorch.descriptor.MCSH import AtomisticMCSH
from amptorch.descriptor.descriptor_calculator import DescriptorCalculator
from ase.io.trajectory import Trajectory
from ase.io import read
from ase import Atoms
from ase.calculators.emt import EMT


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


Gs = {
    "default": {
        "G2": {
            "etas": np.logspace(np.log10(0.05), np.log10(5.0), num=4),
            "rs_s": [0] * 4,
        },
        "G4": {"etas": [0.005], "zetas": [1.0], "gammas": [-1.0, 1.0]},
        "cutoff": 6.5,
    },
}

elements = image.get_chemical_symbols()

descriptor = Gaussian(Gs=Gs, elements=elements)

dataset = AMPTorchDataset(images, descriptor, training_data=True)
sys.exit()


print(len(dataset))
print(dataset[0])
dataloader = DataLoader(dataset, batch_size=2, shuffle=False, collate_fn=collate_amp)
k = next(iter(dataloader))

architecture = [36, 3, 5]
device = "cpu"
forcetraining = True

model = BPNN(unique_atoms, architecture, device, forcetraining)
model(k)
