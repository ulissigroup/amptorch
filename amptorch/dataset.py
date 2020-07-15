"""data_preprocess.py: Defines the PyTorch required Dataset class. Takes in raw
training data, computes/scales image fingerprints and potential energies.
Functions are included to factorize the data and organize it accordingly in
'collate_amp' as needed to be fed into the PyTorch required DataLoader class"""

import copy
import os
import numpy as np
import torch
from torch.utils.data import Dataset, SubsetRandomSampler
from collections import OrderedDict
import ase
from amptorch.gaussian import make_symmetry_functions, SNN_Gaussian
from amptorch.data_utils import Transform
from amptorch.utils import (
    make_amp_descriptors_simple_nn,
    calculate_fingerprints_range,
    hash_images,
    get_hash,
)
from amp.utilities import hash_images as amp_hash
from amp.utilities import get_hash as get_amp_hash

from torch_geometric.data import Data, Batch

__author__ = "Muhammed Shuaibi"
__email__ = "mshuaibi@andrew.cmu.edu"


class AtomsDataset(Dataset):
    def __init__(self, images, descriptor, Gs, forcetraining, label, cores):
        self.images = images
        self.base_descriptor = descriptor
        self.descriptor = descriptor
        self.Gs = Gs
        self.atom_images = self.images
        self.forcetraining = forcetraining
        self.cores = cores

        if isinstance(images, str):
            extension = os.path.splitext(images)[1]
            if extension != (".traj" or ".db"):
                self.atom_images = ase.io.read(images, ":")
        self.elements = self.unique()

        print("Calculating fingerprints...")
        G2_etas = Gs["G2_etas"]
        G2_rs_s = Gs["G2_rs_s"]
        G4_etas = Gs["G4_etas"]
        G4_zetas = Gs["G4_zetas"]
        G4_gammas = Gs["G4_gammas"]
        cutoff = Gs["cutoff"]

        # create simple_nn fingerprints
        self.hashed_images = hash_images(self.atom_images, Gs=Gs)
        make_amp_descriptors_simple_nn(
            self.atom_images, Gs, self.elements, cores=cores, label=label
        )
        G = make_symmetry_functions(
            elements=self.elements, type="G2", etas=G2_etas)
        G += make_symmetry_functions(
            elements=self.elements,
            type="G4",
            etas=G4_etas,
            zetas=G4_zetas,
            gammas=G4_gammas,
        )
        for g in list(G):
            g["Rs"] = G2_rs_s
        self.descriptor = self.descriptor(Gs=G, cutoff=cutoff)
        self.descriptor.calculate_fingerprints(
            self.hashed_images, calculate_derivatives=forcetraining
        )

        print("Fingerprints Calculated!")
        self.fprange = calculate_fingerprints_range(
            self.descriptor, self.hashed_images)

        # perform preprocessing
        self.dataset = self.process()

    def __len__(self):
        return len(self.atom_images)

    def process(self):
        data_list = []
        for atoms in self.images:
            hash_name = get_hash(atoms, self.Gs)

            # scale fingerprints
            image_fingerprint = self.descriptor.fingerprints[hash_name]

            #convert amp structure
            fp_length = len(image_fingerprint[0][1])
            natoms = len(image_fingerprint)
            image_fingerprint = torch.tensor(np.vstack([fp[1] for fp in image_fingerprint]))
            potential_energy = self.hashed_images[hash_name].get_potential_energy(
                apply_constraint=False
            )
            atomic_numbers = torch.tensor(atoms.get_atomic_numbers())

            data = Data(
                fingerprint=image_fingerprint,
                atomic_numbers=atomic_numbers,
                energy=potential_energy,
                natoms=natoms,
            )

            if self.forcetraining:
                image_forces = self.hashed_images[hash_name].get_forces(
                        apply_constraint=False
                    )
                # scale fingerprintprimes
                image_fprimes = self.descriptor.fingerprintprimes[hash_name]
                image_prime_values = list(image_fprimes.values())
                image_prime_keys = list(image_fprimes.keys())
                # convert amp fprime format to matrix
                fingerprintprimes = torch.zeros(
                    fp_length * natoms, 3 * natoms
                )
                for idx, fp_key in enumerate(image_prime_keys):
                    image_prime = torch.tensor(image_prime_values[idx])
                    base_atom = fp_key[2]
                    wrt_atom = fp_key[0]
                    coord = fp_key[4]
                    fingerprintprimes[
                        base_atom * fp_length: base_atom * fp_length + fp_length,
                        wrt_atom * 3 + coord,
                    ] = image_prime

                data.forces = torch.tensor(image_forces)
                data.fprimes = fingerprintprimes

            data_list.append(data)
        return data_list

    def unique(self):
        elements = np.array(
                [atom.symbol for atoms in self.atom_images for atom in atoms])
        _, idx = np.unique(elements, return_index=True)
        elements = list(elements[np.sort(idx)])
        return elements

    def __getitem__(self, index):
        return self.dataset[index]


def collate_amp(data_list):
    batch = Batch.from_data_list(data_list)
    return batch
