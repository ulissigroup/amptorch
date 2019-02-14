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
from collections import defaultdict



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
        return {index:(image_fingerprint,image_potential_energy)}

training_data=AtomsDataset(descriptor=Gaussian())
sample_batch=[training_data[1], training_data[0], training_data[3],
        training_data[18] ]
# print(sample_batch)

def data_factorization(training_data):
    """
    Reads in dataset and factors it into 3 dictionaries:

    1. unique_atoms = Identifies the unique elements in the dataset
    2. fingerprint_dict = Extracts the fingerprints for each hashed data sample in the
    dataset
    3. energy_dict = Extracts the potential energy for a given hashed data sample in the
    dataset
    """
    unique_atoms=[]
    fingerprint_set=[]
    sample_indices=[]
    energy_set=[]
    #Create empty dictionary to store indices of data
    for data_sample in training_data:
        sample_index=data_sample.keys()[0]
        sample_indices.append(sample_index)
        atom_image=data_sample[sample_index]
        atom_fingerprint=atom_image[0]
        fingerprint_set.append(atom_fingerprint)
        image_potential_energy=atom_image[1]
        energy_set.append(image_potential_energy)
        for atom in atom_fingerprint:
            element=atom[0]
            if element not in unique_atoms:
                unique_atoms.append(element)
    return unique_atoms,fingerprint_set,energy_set,sample_indices

x,y,z, i=data_factorization(sample_batch)
print x
print ' '
print y
print ' '
print z
print ' '
print i



def collate_amp(training_data):

    unique_atoms,fingerprint_dict,energy_dict=data_factorization(training_data)
    element_specific_fingerprints={}
    for element in unique_atoms:
        element_specific_fingerprints[element]=defaultdict(list)
    for fingerprint_hash in fingerprint_dict.keys():
        for fingerprint in fingerprint_dict[fingerprint_hash]:
            atom_element=fingerprint[0]
            atom_fingerprint=fingerprint[1]
            element_specific_fingerprints[atom_element][fingerprint_hash].append(torch.tensor(atom_fingerprint,dtype=torch.float64))
    for element in unique_atoms:
        element_specific_fingerprints[element]=dict(element_specific_fingerprints[element])
    return element_spe

def data_factorization_old(training_data):
    """
    Reads in dataset and factors it into 3 dictionaries:

    1. unique_atoms = Identifies the unique elements in the dataset
    2. fingerprint_dict = Extracts the fingerprints for each hashed data sample in the
    dataset
    3. energy_dict = Extracts the potential energy for a given hashed data sample in the
    dataset
    """
    unique_atoms={}
    fingerprint_dict={}
    energy_dict={}
    #Create empty dictionary to store indices of data
    for data_sample in training_data:
        atom_hash=data_sample.keys()[0]
        atom_image=data_sample[atom_hash]
        atom_fingerprint=atom_image[0]
        fingerprint_dict[atom_hash]=atom_fingerprint
        image_potential_energy=atom_image[1]
        energy_dict[atom_hash]=image_potential_energy
        for atom in atom_fingerprint:
            element=atom[0]
            if element not in unique_atoms.keys():
                unique_atoms[element]=[]
    return unique_atoms,fingerprint_dict,energy_dict

# x,y,z=data_factorization(sample_batch)

def collate_amp_old(training_data):
    unique_atoms,fingerprint_dict,energy_dict=data_factorization(training_data)
    element_specific_fingerprints={}
    for element in unique_atoms:
        element_specific_fingerprints[element]=defaultdict(list)
    for fingerprint_hash in fingerprint_dict.keys():
        for fingerprint in fingerprint_dict[fingerprint_hash]:
            atom_element=fingerprint[0]
            atom_fingerprint=fingerprint[1]
            element_specific_fingerprints[atom_element][fingerprint_hash].append(torch.tensor(atom_fingerprint,dtype=torch.float64))
    for element in unique_atoms:
        element_specific_fingerprints[element]=dict(element_specific_fingerprints[element])
    return element_specific_fingerprints,energy_dict

# test,energy=collate_amp(sample_batch)
# print(test)
# print(energy)



# dataloader=DataLoader(sample_batch,batch_size=2,collate_fn=collate_amp,shuffle=False)
# for i in dataloader:
    # print i

#atoms_dataloader=DataLoader(test_data,batch_size=1,shuffle=False,sampler=None,batch_sampler=None,num_workers=0,collate_fn=,pin_memory=False,drop_last=False,timeout=0,worker_init_fn=None)


