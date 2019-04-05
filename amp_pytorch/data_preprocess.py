import copy
import os
import numpy as np
import torch
from torch.utils.data import Dataset, SubsetRandomSampler
from amp.utilities import hash_images
from amp.model import calculate_fingerprints_range
import ase


class AtomsDataset(Dataset):
    """
    Atoms Dataset

    Parameters:
    images: input training images (list, trajectory file, or database)
    descriptor: descriptor to be utilized for fingerprinting scheme

    Output: Returns, for a given index, the image_fingerprint and
    image_potential energy
    """

    def __init__(self, images, descriptor):
        self.images = images
        self.descriptor = descriptor
        if isinstance(images, str):
            extension = os.path.splitext(images)[1]
            if extension != ('.traj' or '.db'):
                self.atom_images = ase.io.read(images, ':')
            else:
                self.atom_images = self.images
        self.hashed_images = hash_images(self.atom_images)
        self.descriptor.calculate_fingerprints(self.hashed_images)
        self.fprange = calculate_fingerprints_range(
            self.descriptor, self.hashed_images)

    def __len__(self):
        return len(self.hashed_images)

    def __getitem__(self, index):
        hash_name = self.hashed_images.keys()[index]
        image_fingerprint = self.descriptor.fingerprints[hash_name]
        fprange = self.fprange
        # fingerprint scaling to [-1,1]
        for i, (atom, afp) in enumerate(image_fingerprint):
            _afp = copy.copy(afp)
            fprange_atom = fprange[atom]
            for _ in range(np.shape(_afp)[0]):
                if(fprange_atom[_][1]-fprange_atom[_][0]) > (10.**(-8.)):
                    _afp[_] = -1+2.*((_afp[_]-fprange_atom[_][0]) /
                                     (fprange_atom[_][1]-fprange_atom[_][0]))
            image_fingerprint[i] = (atom, _afp)
        try:
            image_potential_energy = self.hashed_images[hash_name].get_potential_energy(
                apply_constraint=False)
        except:
            print 'Atoms object has no claculator set!'
        return image_fingerprint, image_potential_energy

    def create_splits(self, training_data, val_frac):
        dataset_size = len(training_data)
        indices = np.random.permutation(dataset_size)
        split = int(np.floor(val_frac*dataset_size))
        train_idx, val_idx = indices[split:], indices[:split]
        train_sampler = SubsetRandomSampler(train_idx)
        val_sampler = SubsetRandomSampler(val_idx)

        samplers = {'train': train_sampler, 'val': val_sampler}

        return samplers


def factorize_data(training_data):
    """
    Reads in dataset and factors it into 3 lists:

    1. unique_atoms = Identifies the unique elements in the dataset
    2. fingerprint_dataset = Extracts the fingerprints for each hashed data
    sample in the dataset
    3. energy_dataset = Extracts the potential energy for a given hashed data
    sample in the dataset
    """
    unique_atoms = []
    fingerprint_dataset = []
    energy_dataset = []
    num_of_atoms = []
    for image_sample in training_data:
        image_fingerprint = image_sample[0]
        fingerprint_dataset.append(image_fingerprint)
        image_potential_energy = image_sample[1]
        energy_dataset.append(image_potential_energy)
        num_of_atoms.append(float(len(image_fingerprint)))
        for atom in image_fingerprint:
            element = atom[0]
            if element not in unique_atoms:
                unique_atoms.append(element)
    return unique_atoms, fingerprint_dataset, energy_dataset, num_of_atoms


def collate_amp(training_data):
    """
    Reshuffling scheme that reads in raw data and organizes it into element
    specific datasets to be fed into the element specific Neural Nets.
    """

    unique_atoms, fingerprint_dataset, energy_dataset, num_of_atoms = factorize_data(
        training_data)
    element_specific_fingerprints = {}
    model_input_data = []
    for element in unique_atoms:
        element_specific_fingerprints[element] = [[], []]
    for fp_index, sample_fingerprints in enumerate(fingerprint_dataset):
        for fingerprint in sample_fingerprints:
            atom_element = fingerprint[0]
            atom_fingerprint = fingerprint[1]
            element_specific_fingerprints[atom_element][0].append(
                torch.tensor(atom_fingerprint))
            element_specific_fingerprints[atom_element][1].append(fp_index)
    for element in unique_atoms:
        element_specific_fingerprints[element][0] = torch.stack(
            element_specific_fingerprints[element][0])
    model_input_data.append(element_specific_fingerprints)
    model_input_data.append(torch.tensor(energy_dataset))
    model_input_data.append(torch.tensor(num_of_atoms))
    return model_input_data
