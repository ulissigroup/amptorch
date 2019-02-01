import os
import sys
import numpy as np
from ase import Atoms
import torch

"""Loading Data"""
from torch.utils.data import Dataset, DataLoader
from amp.utilities import hash_images
from amp.utilities import check_images
from amp.descriptor.gaussian import Gaussian

class AtomsDataset(Dataset):
    """Atoms dataset"""

    def __init__(self,filename='sample_training_data.traj',descriptor=Gaussian):
        #raw_data=sys.argv[1]
        raw_data=hash_images(filename)
        check_images(raw_data,forces=False)
        self.raw_data=raw_data
        self.descriptor=Gaussian()

    def __len__(self):
        return len(self.raw_data)

    def __getitem__(self,index):
        hash_name=self.raw_data.keys()[index]
        image=self.raw_data.values()[index]
        self.descriptor.calculate_fingerprints(self.raw_data)
        return


test_data=AtomsDataset()
print(test_data.raw_data.values()[0].get_positions())
print test_data.__getitem__(0).values()[0].get_positions()
