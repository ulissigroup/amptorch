import sys
import ase
import numpy as np
import torch
import time
from torch.utils.data import Dataset, DataLoader
from amp.descriptor.gaussian import Gaussian
from amp_pytorch.NN_model import FullNN
from functools import lru_cache
# Pre-process data
from pre_data import AtomsDataset, collate_amp
# Original
# from data_preprocess import AtomsDataset, collate_amp


# images = ase.io.read("../datasets/water/water.extxyz", ":")
images = ase.io.read("../datasets/COCu/COCu_pbc_500.traj", ":")
data_size = 100
batches = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
device = 'cpu'

for batch_size in batches:
    k = time.time()
    dataset = AtomsDataset(images[:data_size], Gaussian(), cores=1, forcetraining=True)
    print('Preprocess time: %s ' % (time.time()-k))

    unique_atoms = dataset.unique()
    fp_length = dataset.fp_length
    structure = [fp_length, 5, 5]
    model = FullNN(unique_atoms, structure, device, forcetraining=True).to(device)

    dataloader = DataLoader(dataset, batch_size=batch_size, collate_fn=collate_amp)

    loading_time = 0
    calculation_time = 0
    to = time.time()
    for i, sample in enumerate(dataloader):
        loading_time = time.time() - to
        input_data = [sample[0], len(sample[1])]
        for element in unique_atoms:
            input_data[0][element][0] = (input_data[0][element][0].to(device).requires_grad_(True))
        fp_primes = sample[3]
        f = time.time()
        energy, forces = model(input_data, fp_primes)
        calc_time = time.time() - f
        calculation_time += calc_time
    print('total loading time: %s' % loading_time)
    print('total calculation time: %s' % calculation_time)
