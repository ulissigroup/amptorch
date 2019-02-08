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
sample=training_data[0]
sample_batch=[training_data[i] for i in range(2)]


def collate_amp(training_data):
    unique_atoms,fingerprint_dict,energy_dict=identify_unique_atoms(training_data)
    element_specific_fingerprint={}
    for element_type in unique_atoms:
        element_specific_fingerprint[element_type]=[]
    for index,fingerprint_sample in enumerate(fingerprint_dict.values()):
        for atom_fingerprint in fingerprint_sample:
            atom_fingerprint[1].append(index)
            for element in unique_atoms:
                if atom_fingerprint[0]==element:
                    element_specific_fingerprint[element].append(atom_fingerprint[1])
    return element_specific_fingerprint

def identify_unique_atoms(training_data):
    unique_atoms=[]
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
            if element not in unique_atoms:
                unique_atoms.append(element)
    return unique_atoms,fingerprint_dict,energy_dict



unique_atoms, fingerprint_dict, energy =identify_unique_atoms(sample_batch)
#print fingerprint_dict[19]


test=collate_amp(sample_batch)
print(test)
#Cu_test=test['Cu']
# print('')
# Pt_test=test['Pt']
# print(len(Cu_test))
# print(len(Pt_test))
# dataloader=DataLoader(sample_batch,batch_size=1,collate_fn=collate_amp,shuffle=True)
# for i in dataloader:
    # print i
    # sys.exit()
    # print len(i)
    # sys.exit()




#atoms_dataloader=DataLoader(test_data,batch_size=1,shuffle=False,sampler=None,batch_sampler=None,num_workers=0,collate_fn=,pin_memory=False,drop_last=False,timeout=0,worker_init_fn=None)


