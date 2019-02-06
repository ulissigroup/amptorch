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
        image_fingerprint=self.descriptor.fingerprints[hash_name]
        image_potential_energy=self.atom_images[hash_name].get_potential_energy()
        return image_fingerprint,image_potential_energy

training_data=AtomsDataset(descriptor=Gaussian())
sample=training_data[19]
#print(sample)




def identify_unique_atoms(training_data):
    unique_atoms=[]
    for sample in training_data:
        atom_image=sample[0]
        image_potential_energy=sample[1]
        for atom in atom_image:
            element=atom[0]
            if element not in unique_atoms:
                unique_atoms.append(element)
    return unique_atoms

def collate_amp(training_data):
    unique_atoms=identify_unique_atoms(training_data)

unique_atoms=identify_unique_atoms(training_data)
print(unique_atoms)

#collate_amp(training_data)
#atoms_dataloader=DataLoader(test_data,batch_size=1,shuffle=False,sampler=None,batch_sampler=None,num_workers=0,collate_fn=,pin_memory=False,drop_last=False,timeout=0,worker_init_fn=None)
