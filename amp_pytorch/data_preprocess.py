"""data_preprocess.py: Defines the PyTorch required Dataset class. Takes in raw
training data, computes/scales image fingerprints and potential energies.
Functions are included to factorize the data and organize it accordingly in
'collate_amp' as needed to be fed into the PyTorch required DataLoader class"""

import sys
import copy
import os
import time
import numpy as np
import torch
from torch.utils.data import Dataset, SubsetRandomSampler
from amp.utilities import hash_images
from amp.model import calculate_fingerprints_range
from amp.descriptor.gaussian import Gaussian
import ase

__author__ = "Muhammed Shuaibi"
__email__ = "mshuaibi@andrew.cmu.edu"


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
            if extension != (".traj" or ".db"):
                self.atom_images = ase.io.read(images, ":")
        else:
            self.atom_images = self.images
        self.hashed_images = hash_images(self.atom_images)
        self.descriptor.calculate_fingerprints(
            self.hashed_images, calculate_derivatives=True
        )
        self.fprange = calculate_fingerprints_range(self.descriptor, self.hashed_images)

    def __len__(self):
        return len(self.hashed_images)

    def __getitem__(self, index):
        hash_name = list(self.hashed_images.keys())[index]
        image_fingerprint = self.descriptor.fingerprints[hash_name]
        fprange = self.fprange
        # fingerprint scaling to [-1,1]
        for i, (atom, afp) in enumerate(image_fingerprint):
            _afp = copy.copy(afp)
            fprange_atom = fprange[atom]
            for _ in range(np.shape(_afp)[0]):
                if (fprange_atom[_][1] - fprange_atom[_][0]) > (10.0 ** (-8.0)):
                    _afp[_] = -1 + 2.0 * (
                        (_afp[_] - fprange_atom[_][0])
                        / (fprange_atom[_][1] - fprange_atom[_][0])
                    )
            image_fingerprint[i] = (atom, _afp)
        try:
            image_potential_energy = self.hashed_images[hash_name].get_potential_energy(
                apply_constraint=False
            )
            image_forces = self.hashed_images[hash_name].get_forces(
                apply_constraint=False
            )
            image_primes = self.descriptor.fingerprintprimes[hash_name]
        except NotImplementedError:
            print("Atoms object has no claculator set!")
        return image_fingerprint, image_potential_energy, image_primes, image_forces

    def create_splits(self, training_data, val_frac):
        dataset_size = len(training_data)
        indices = np.random.permutation(dataset_size)
        split = int(val_frac * dataset_size)
        train_idx, val_idx = indices[split:], indices[:split]
        train_sampler = SubsetRandomSampler(train_idx)
        val_sampler = SubsetRandomSampler(val_idx)

        samplers = {"train": train_sampler, "val": val_sampler}

        return samplers


def factorize_data(training_data):
    """
    Reads in dataset and factors it into 3 lists:

    1. unique_atoms = Identifies the unique elements in the dataset

    2. fingerprint_dataset = Extracts the fingerprints for each hashed data
    sample in the dataset

    3. energy_dataset = Extracts the potential energy for a given hashed data
    sample in the dataset

    4. num_of atoms = Number of atoms per image

    5. factored_fingerprints = fingerprint derivatives for each atom across all
    images. Dimensions of said tensor is QxPx3:
        Q - Atoms in batch
        P - Length of fingerprints
        3 - x,y,z components of the fingerprint derivatives. It should be noted
        that for each image containing N atoms, there are 3N directional
        components for each atom. For each atom, N contributions are summed for
        each directional component to arrive at a dimension of 3 for each atom.

    """
    unique_atoms = []
    fingerprint_dataset = []
    energy_dataset = []
    num_of_atoms = []
    fingerprintprimes = torch.tensor([])
    grouped_fp_idxs = torch.tensor([])
    image_forces = []
    fp_length = len(training_data[0][0][0][1])
    for image_sample in training_data:
        image_fingerprint = image_sample[0]
        fingerprint_dataset.append(image_fingerprint)
        image_potential_energy = image_sample[1]
        energy_dataset.append(image_potential_energy)
        num_atom = float(len(image_fingerprint))
        num_of_atoms.append(num_atom)
        image_primes = list(image_sample[2].values())
        fingerprintprimes = torch.cat(
            (fingerprintprimes, torch.tensor(image_primes)), 0
        )
        image_forces.append(torch.from_numpy(image_sample[3]))
        grouped_idx = torch.tensor([num_atom] * int((3 * num_atom)))
        grouped_fp_idxs = torch.cat((grouped_fp_idxs, grouped_idx), 0)
        for atom in image_fingerprint:
            element = atom[0]
            if element not in unique_atoms:
                unique_atoms.append(element)
    image_forces = torch.cat(image_forces)
    atoms_in_batch = int(sum(num_of_atoms))
    fp_groupings = copy.deepcopy(grouped_fp_idxs)
    fp_groupings = [int(i) for i in fp_groupings]
    fp_groupings = [
        idx.repeat(times)
        for idx, times in zip(torch.arange(len(fp_groupings)), fp_groupings)
    ]
    fp_groupings = torch.cat(fp_groupings)
    factored_fingerprints = torch.zeros(len(grouped_fp_idxs), fp_length)
    factored_fingerprints.index_add_(0, fp_groupings, fingerprintprimes)
    factored_fingerprintprimes = factored_fingerprints.reshape(
        atoms_in_batch, 3, fp_length 
    ).transpose(1, 2)
    # factored_fingerprintprimes is of dimensionality QxPx3 with each column
    # representing: dFP/dx dFP/dy dFP/dz
    return unique_atoms, fingerprint_dataset, energy_dataset, num_of_atoms, factored_fingerprintprimes, image_forces


def collate_amp(training_data):
    """
    Reshuffling scheme that reads in raw data and organizes it into element
    specific datasets to be fed into the element specific Neural Nets.
    """

    unique_atoms, fingerprint_dataset, energy_dataset, num_of_atoms, fp_primes, image_forces = factorize_data(
        training_data
    )
    element_specific_fingerprints = {}
    model_input_data = []
    for element in unique_atoms:
        element_specific_fingerprints[element] = [[], []]
    for fp_index, sample_fingerprints in enumerate(fingerprint_dataset):
        for fingerprint in sample_fingerprints:
            atom_element = fingerprint[0]
            atom_fingerprint = fingerprint[1]
            element_specific_fingerprints[atom_element][0].append(
                torch.tensor(atom_fingerprint)
            )
            element_specific_fingerprints[atom_element][1].append(fp_index)
    for element in unique_atoms:
        element_specific_fingerprints[element][0] = torch.stack(
            element_specific_fingerprints[element][0]
        )
    model_input_data.append(element_specific_fingerprints)
    model_input_data.append(torch.tensor(energy_dataset))
    model_input_data.append(torch.tensor(num_of_atoms))
    model_input_data.append(fp_primes)
    model_input_data.append(image_forces)
    return model_input_data
