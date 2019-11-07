"""An example of how to utilize the package to train on energies and forces"""

import sys
import time
from torch.utils.data import DataLoader
from amptorch.gaussian import Gaussian
from amptorch.data_preprocess import AtomsDataset, collate_amp, factorize_data
# from amptorch.data_backup import AtomsDataset, collate_amp, factorize_data
from ase.io import read

# define training images
label = "example"
# images = read("./trainingset.traj", ":")
images = read("../datasets/COCu/COCu_pbc_300K.traj", ":")

# define symmetry functions to be used
Gs = {}
Gs["G2_etas"] = [0.005, 4.0, 20.0, 80.0]
Gs["G2_rs_s"] = [0] * 4
Gs["G4_etas"] = [0.005]
Gs["G4_zetas"] = [1.0, 4.0]
Gs["G4_gammas"] = [+1.0, -1]
Gs["cutoff"] = 4.0

dataset = AtomsDataset(
    images,
    descriptor=Gaussian,
    Gs=Gs,
    cores=1,
    forcetraining=True,
    lj_data=None,
    label="example",
)

batch = len(dataset)

dataloader = DataLoader(
    dataset, collate_fn=collate_amp, batch_size=batch, shuffle=False, num_workers=0
)

for k in range(100):
    t = time.time()
    for i in dataloader:
        print(time.time() - t)
