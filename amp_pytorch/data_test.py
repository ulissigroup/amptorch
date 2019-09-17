import sys
import ase
import numpy as np
import torch
import time
from torch.utils.data import Dataset, DataLoader
from amp.descriptor.gaussian import Gaussian
from functools import lru_cache
# Pre-process data
from pre_data import AtomsDataset, collate_amp
# Original
# from data_preprocess import AtomsDataset, collate_amp


# images = ase.io.read("../datasets/water/water.extxyz", ":")
images = ase.io.read("../datasets/COCu/COCu_pbc_500.traj", ":")
batch_size = 100
k = time.time()
dataset = AtomsDataset(images[:batch_size], Gaussian(), cores=1, forcetraining=True)
print('Preprocess time: %s' % (time.time()-k))
dataloader = DataLoader(dataset, batch_size=batch_size, collate_fn=collate_amp)

for iter in range(5):
    t = time.time()
    for i, sample in enumerate(dataloader):
        print('Iter %s: time taken = %s' % (iter+1, time.time()-t))

