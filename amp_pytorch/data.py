import sys
import copy
import os
import numpy as np
import torch
from torch.utils.data import Dataset,SubsetRandomSampler
from amp.utilities import check_images,hash_images
from amp.descriptor.gaussian import Gaussian
from amp.model import calculate_fingerprints_range
import ase
from ase import io

class AtomsDataset(Dataset):
    """
    Atoms dataset
    Parameters: Descriptor type and input images(list,trajectory file, or a
    database)
    Output: Returns, for a given index, the image_fingerprint and image_potential
    energy
    """

    def __init__(self,images,descriptor):
        self.images=images
        self.descriptor=descriptor
        if isinstance(images,str):
            extension=os.path.splitext(images)[1]
            if extension != ('.traj' or '.db'):
                self.atom_images=ase.io.read(images,':')
            else:
                self.atom_images=self.images
        self.hash_images=hash_images(self.atom_images)
        self.descriptor.calculate_fingerprints(self.hash_images)
        #fprange is the min,max across each fingerprint for the entire dataset
        self.fprange=calculate_fingerprints_range(self.descriptor,self.hash_images)

    def __len__(self):
        return len(self.hash_images)

    def __getitem__(self,index):
        hash_name=self.hash_images.keys()[index]
        image_fingerprint=self.descriptor.fingerprints[hash_name]
        fprange=self.fprange
        for i,(symbol,afp) in enumerate(image_fingerprint):
            _afp=copy.copy(afp)
            fprange_atom=fprange[symbol]
            for _ in range(np.shape(_afp)[0]):
                if(fprange_atom[_][1]-fprange_atom[_][0])>(10.**(-8.)):
                    _afp[_]=-1+2.*((_afp[_]-fprange_atom[_][0])/(fprange_atom[_][1]-fprange_atom[_][0]))
            image_fingerprint[i]=(symbol,_afp)
        try:
            image_potential_energy=self.hash_images[hash_name].get_potential_energy(apply_constraint=False)
        except:
            print 'Atoms object has no claculator set! Modify the input images before trying again.'
        return image_fingerprint,image_potential_energy

    def data_split(self,training_data,val_frac):
        dataset_size=len(training_data)
        indices=np.random.permutation(dataset_size)
        split=int(np.floor(val_frac*dataset_size))
        train_idx,val_idx=indices[split:],indices[:split]
        train_sampler=SubsetRandomSampler(train_idx)
        val_sampler=SubsetRandomSampler(val_idx)

        samplers={'train':train_sampler,'val':val_sampler}

        return samplers

def data_factorization(training_data):
    """
    Reads in dataset and factors it into 4 lists:

    1. unique_atoms = Identifies the unique elements in the dataset
    2. fingerprint_dict = Extracts the fingerprints for each hashed data sample in the
    dataset
    3. energy_dict = Extracts the potential energy for a given hashed data sample in the
    dataset
    4. sample_indices = Identifies indices of corresponding fingerprints
    """
    unique_atoms=[]
    fingerprint_dataset=[]
    # sample_indices=[]
    energy_dataset=[]
    #Create empty dictionary to store indices of data
    for data_sample in training_data:
        atom_fingerprint=data_sample[0]
        fingerprint_dataset.append(atom_fingerprint)
        image_potential_energy=data_sample[1]
        energy_dataset.append(image_potential_energy)
        for atom in atom_fingerprint:
            element=atom[0]
            if element not in unique_atoms:
                unique_atoms.append(element)
    return unique_atoms,fingerprint_dataset,energy_dataset

def collate_amp(training_data):
    unique_atoms,fingerprint_dataset,energy_dataset=data_factorization(training_data)
    element_specific_fingerprints={}
    model_input_data=[]
    for element in unique_atoms:
        element_specific_fingerprints[element]=[[],[]]
    for fp_index,sample_fingerprints in enumerate(fingerprint_dataset):
        for fingerprint in sample_fingerprints:
            atom_element=fingerprint[0]
            atom_fingerprint=fingerprint[1]
            element_specific_fingerprints[atom_element][0].append(torch.tensor(atom_fingerprint))
            element_specific_fingerprints[atom_element][1].append(fp_index)
    for element in unique_atoms:
        element_specific_fingerprints[element][0]=torch.stack(element_specific_fingerprints[element][0])
    model_input_data.append(element_specific_fingerprints)
    model_input_data.append(torch.tensor(energy_dataset))
    return model_input_data
