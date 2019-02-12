import os
import sys
import numpy as np
from ase import Atoms
import torch
import copy
"""Loading Data"""
from torch.utils.data import Dataset, DataLoader
from amp.utilities import hash_images
from amp.utilities import check_images
from amp.descriptor.gaussian import Gaussian

class AtomsDataset(Dataset):
    """
    Atoms dataset
    Parameters: Descriptor type and .traj file name
    Output: Returns, for a given index, the image_fingerprint and image_potential
    energy
    """

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
sample_batch=[training_data[i] for i in range(2)]

def data_factorization(training_data):
    """
    Reads in dataset and factors it into 3 dictionaries:

    unique_atoms = Identifies the unique elements in the dataset

    fingerprint_dict = Extracts the fingerprints for each data sample in the
    dataset

    energy_dict = Extracts the potential energy for a given data sample in the
    dataset
    """
    unique_atoms={}
    fingerprint_dict={}
    energy_dict={}
    #Create empty dictionary to store indices of data
    for index,sample in enumerate(training_data):
        atom_image=sample[0]
        fingerprint_dict[index]=atom_image
        image_potential_energy=sample[1]
        energy_dict[index]=image_potential_energy
        for atom in atom_image:
            element=atom[0]
            if element not in unique_atoms.keys():
                unique_atoms[element]=[]
    return unique_atoms,fingerprint_dict,energy_dict

def collate_amp(training_data):
    """
    Collate function to be utilized by PyTorch's DataLoader .

    Reads in a dataset and outputs a dictionary of dictionaries indexed for
    each data sample, with fingerprints factored for each element.

    e.g.

    {0:{'Cu':[[fingerprints]],'Pt':[[fingerprints]]},1:{'Cu':[[fingerprints]],'Pt':[[fingerprints]]}}

    """


    unique_atoms,fingerprint_dict,energy_dict=data_factorization(training_data)
    # print(fingerprint_dict.values())
    element_specific_fingerprints_idxd={}
    for i in fingerprint_dict.keys():
        element_specific_fingerprints_idxd[i]=copy.deepcopy(unique_atoms)
    for index,fingerprint_sample in enumerate(fingerprint_dict.values()):
        for atom_fingerprint in fingerprint_sample:
            element=atom_fingerprint[0]
            element_specific_fingerprints_idxd[index][element].append(torch.tensor(atom_fingerprint[1]))
    return element_specific_fingerprints_idxd


test=collate_amp(sample_batch)
# dataloader=DataLoader(sample_batch,batch_size=1,collate_fn=collate_amp,shuffle=False)
# for i in dataloader:
    # print i

#atoms_dataloader=DataLoader(test_data,batch_size=1,shuffle=False,sampler=None,batch_sampler=None,num_workers=0,collate_fn=,pin_memory=False,drop_last=False,timeout=0,worker_init_fn=None)


