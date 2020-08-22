import numpy as np

from amptorch.dataset import (
    AMPTorchDataset,
    collate_amp
)
from amptorch.model_geometric import BPNN


# from torch_geometric.data import DataLoader
from torch.utils.data import DataLoader

import numpy as np
from amptorch.amptorch_descriptor.BP_symmetry_function import BPSymmetryFunction
from amptorch.amptorch_descriptor.Atomistic_MCSH import AtomisticMCSH
from amptorch.amptorch_descriptor.descriptor_calculator import DescriptorCalculator
from ase.io.trajectory import Trajectory
from ase.io import read

# large = Trajectory('./large/iron_data.traj')
# trajectories = [large]
# elements = ["H","O","Fe"]
# unique_atoms = [29, 6, 8]

small = read('./small/water.extxyz', index=':')
trajectories = [small]
elements = ["H","O"]
unique_atoms = [6, 8]


Gs = {"G2": {"etas": np.logspace(np.log10(0.05), np.log10(5.0), num=4), "rs_s": [0] * 4},\
      "G4": {"etas": [0.005], "zetas": [1.0], "gammas": [-1.0, 1.0]},\
      "cutoff": 6.5}

descriptor = BPSymmetryFunction(Gs = Gs, elements = elements)


dataset = AMPTorchDataset(trajectories, descriptor, training_data = True)


print(len(dataset))
print(dataset[0])
dataloader = DataLoader(dataset, batch_size=2, shuffle=False, collate_fn=collate_amp)
k = next(iter(dataloader))

architecture = [36, 3, 5]
device = "cpu"
forcetraining = True

model = BPNN(unique_atoms, architecture, device, forcetraining)
model(k)
