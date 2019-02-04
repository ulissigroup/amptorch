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

    def __init__(self,descriptor,filename='sample_training_data.traj'):
        self.filename=filename
        self.descriptor=descriptor
        self.atom_images=hash_images(filename)
        check_images(self.atom_images,forces=False)
        self.descriptor.calculate_fingerprints(self.atom_images)

    def __len__(self):
        return len(self.atom_images)

    def __getitem__(self,index):
        hash_name=self.atom_images.keys()[index]
        image=self.atom_images.values()[index]
        self.descriptor.calculate_fingerprints(self.atom_images)
        image_fingerprint=self.descriptor.fingerprints[hash_name]
        #total vs potential energy?
        image_energy=self.atom_images[hash_name].get_total_energy()
        return image_fingerprint,image_energy

test_data=AtomsDataset(descriptor=Gaussian())
fp,E=test_data.__getitem__(19)
