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
        self.hash_name=hash_name
        self.index=index
        return {hash_name:(image_fingerprint,image_potential_energy)}

training_data=AtomsDataset(descriptor=Gaussian())
sample_batch=[training_data[i] for i in range(2)]

def data_factorization(training_data):
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

x,y,z=data_factorization(sample_batch)
def collate_amp(training_data):
    """
    Collate function to be utilized by PyTorch's DataLoader .

    Reads in a dataset and outputs a dictionary of dictionaries indexed for
    each data sample, with fingerprints factored for each element.

    e.g.

    {0:{'Cu':[[fingerprints]],'Pt':[[fingerprints]]},1:{'Cu':[[fingerprints]],'Pt':[[fingerprints]]}}

    """


    unique_atoms,fingerprint_dict,energy_dict=data_factorization(training_data)
    print(fingerprint_dict.values())
    element_specific_fingerprints_idxd={}
    for i in fingerprint_dict.keys():
        element_specific_fingerprints_idxd[i]=copy.deepcopy(unique_atoms)
    for index,fingerprint_sample in enumerate(fingerprint_dict.values()):
        for atom_fingerprint in fingerprint_sample:
            element=atom_fingerprint[0]
            element_specific_fingerprints_idxd[index][element].append(torch.tensor(atom_fingerprint[1],dtype=torch.float64))
    return element_specific_fingerprints_idxd

def collate_amp_2(training_data):

    unique_atoms,fingerprint_dict,energy_dict=data_factorization(training_data)
    element_specific_fingerprints={}
    for i in unique_atoms:
        element_specific_fingerprints[i]=defaultdict(list)
    for i in fingerprint_dict.keys():
        for fp in fingerprint_dict[i]:
            element_specific_fingerprints[fp[0]][i].append(torch.tensor(fp[1],dtype=torch.float64))
    for i in unique_atoms:
        element_specific_fingerprints[i]=dict(element_specific_fingerprints[i])
    return element_specific_fingerprints,energy_dict

test,energy=collate_amp_2(sample_batch)
print(test)
print(energy)


# test=collate_amp(sample_batch)
# print(test)

# dataloader=DataLoader(sample_batch,batch_size=2,collate_fn=collate_amp,shuffle=False)
# for i in dataloader:
    # print i

#atoms_dataloader=DataLoader(test_data,batch_size=1,shuffle=False,sampler=None,batch_sampler=None,num_workers=0,collate_fn=,pin_memory=False,drop_last=False,timeout=0,worker_init_fn=None)


