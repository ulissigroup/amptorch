"""data_preprocess.py: Defines the PyTorch required Dataset class. Takes in raw
training data, computes/scales image fingerprints and potential energies.
Functions are included to factorize the data and organize it accordingly in
'collate_amp' as needed to be fed into the PyTorch required DataLoader class"""

import sys
import time
import copy
import os
import numpy as np
import torch
from torch.utils.data import Dataset, SubsetRandomSampler
import scipy.sparse as sparse
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
from functools import lru_cache

__author__ = "Muhammed Shuaibi"
__email__ = "mshuaibi@andrew.cmu.edu"


class AtomsDataset(Dataset):
    """
    PyTorch inherited abstract class that specifies how the training data is to be preprocessed and
    passed along to the DataLoader.

    Parameters:
    -----------
    images: list, file, database
        Training data to be utilized for the regression model.

    descriptor: object
        Scheme to be utilized for computing fingerprints.

    Gs: object
        Symmetry function parameters to be used for hashing and fingerprinting.

    forcetraining: float
        Flag to specify whether force training is to be performed - dataset
        will then compute the necessary information - fingerprint derivatives,
        etc.

    cores: int
        Number of cores to parallelize symmetry function
        computation/reorganization.

    delta_data: list
        Energies and forces to be subtracted off from targets, allowing the
        model to learn the difference. default: None

    store_primes: Boolean
        True to save fingerprintprimes matrices for faster preprocessing.
        Default: False


    """

    def __init__(
        self,
        images,
        descriptor,
        Gs,
        cutoff,
        forcetraining,
        label,
        cores,
        delta_data=None,
        store_primes=False,
        logger=None, 
    ):
        self.images = images
        self.base_descriptor = descriptor
        self.descriptor = descriptor
        self.Gs = Gs
        self.cutoff = cutoff
        self.atom_images = self.images
        self.forcetraining = forcetraining
        self.store_primes = store_primes
        self.cores = cores
        self.delta = False
        if delta_data is not None:
            self.delta_data = delta_data
            self.delta_energies = np.array(delta_data[0])
            self.delta_forces = delta_data[1]
            self.num_atoms = np.array(delta_data[2])
            self.delta = True
        if self.store_primes:
            if not os.path.isdir("./stored-primes/"):
                os.mkdir("stored-primes")
        if isinstance(images, str):
            extension = os.path.splitext(images)[1]
            if extension != (".traj" or ".db"):
                self.atom_images = ase.io.read(images, ":")
        self.elements = self.unique()
        # TODO Print log - control verbose
        print("Calculating fingerprints...")
        # G2_etas = Gs["G2_etas"]
        # G2_rs_s = Gs["G2_rs_s"]
        # G4_etas = Gs["G4_etas"]
        # G4_zetas = Gs["G4_zetas"]
        # G4_gammas = Gs["G4_gammas"]
        # cutoff = Gs["cutoff"]
        # create simple_nn fingerprints
        if descriptor == SNN_Gaussian:
<<<<<<< HEAD
            print('SSN Gaussian: %d images' % len(self.atom_images))
            self.hashed_images = hash_images(self.atom_images, Gs=Gs, log=logger or None)
=======
            self.hashed_images = hash_images(self.atom_images, Gs=Gs)
>>>>>>> ulissi/master
            make_amp_descriptors_simple_nn(
                self.atom_images, Gs, self.elements, cores=cores, label=label
            )
            self.isamp_hash = False
        else:  # ignoring non SSN_Gaussian descriptors for now
            print('Amp Gaussian: %d images' % len(self.atom_images))
            self.hashed_images = amp_hash(self.atom_images)
            self.isamp_hash = True
        # G = make_symmetry_functions(elements=self.elements, # type="G2", etas=G2_etas)
        # G += make_symmetry_functions(
        #     elements=self.elements,
        #     type="G4",
        #     etas=G4_etas,
        #     zetas=G4_zetas,
        #     gammas=G4_gammas,
        # )
        # for g in list(G):
        #     g["Rs"] = G2_rs_s
        # self.descriptor = self.descriptor(Gs=G, cutoff=cutoff)
        self.descriptor = self.descriptor(Gs=Gs, cutoff=cutoff)

        self.descriptor.calculate_fingerprints(
            self.hashed_images, calculate_derivatives=forcetraining
        )
        print("Fingerprints Calculated!")
        self.fprange = calculate_fingerprints_range(self.descriptor, self.hashed_images)
        # perform preprocessing
        self.fingerprint_dataset, self.energy_dataset, self.num_of_atoms, self.sparse_fprimes, self.forces_dataset, self.index_hashes, self.scalings, self.rearange_forces = (
            self.preprocess_data()
        )

    def __len__(self):
        return len(self.atom_images)

    def preprocess_data(self):
        # TODO cleanup/optimize
        fingerprint_dataset = []
        fprimes_dataset = []
        energy_dataset = np.array([])
        num_of_atoms = np.array([])
        forces_dataset = []
        index_hashes = []
        self.fp_length = self.fp_length()
        rearange_forces = {}
        n = 0
        for index, atoms_object in enumerate(self.atom_images):
            if self.isamp_hash:
                hash_name = get_amp_hash(atoms_object)
            else:
                hash_name = get_hash(atoms_object, self.Gs)
            index_hashes.append(hash_name)
            image_fingerprint = self.descriptor.fingerprints[hash_name]
            n_atoms = float(len(image_fingerprint))
            num_of_atoms = np.append(num_of_atoms, n_atoms)
            fprange = self.fprange
            atom_order = []
            # fingerprint scaling to [-1,1]
            for i, (atom, afp) in enumerate(image_fingerprint):
                _afp = copy.copy(afp)
                fprange_atom = np.array(fprange[atom])
                for _ in range(np.shape(_afp)[0]):
                    if (fprange_atom[_][1] - fprange_atom[_][0]) > (10.0 ** (-8.0)):
                        _afp[_] = -1 + 2.0 * (
                            (_afp[_] - fprange_atom[_][0])
                            / (fprange_atom[_][1] - fprange_atom[_][0])
                        )
                image_fingerprint[i] = (atom, _afp)
                atom_order.append(atom)
            fingerprint_dataset.append(image_fingerprint)
            image_potential_energy = (
                self.hashed_images[hash_name].get_potential_energy(
                    apply_constraint=False
                )
                / n_atoms
            )
            energy_dataset = np.append(energy_dataset, image_potential_energy)
            if self.forcetraining:
                image_forces = (
                    self.hashed_images[hash_name].get_forces(apply_constraint=False)
                    / n_atoms
                )
                # subtract off delta force contributions
                if self.delta:
                    delta_forces = self.delta_forces[index] / n_atoms
                    image_forces -= delta_forces
                if self.store_primes and os.path.isfile("./stored-primes/" + hash_name):
                    pass
                else:
                    prime_mapping = []
                    for element in self.elements:
                        indices = [i for i, x in enumerate(atom_order) if x == element]
                        prime_mapping += indices
                    new_order = [atom_order[i] for i in prime_mapping]
                    used = set()
                    t = np.array([])
                    for i, x in enumerate(atom_order):
                        for k, l in enumerate(new_order):
                            if (x == l) and (k not in used):
                                used.add(k)
                                t = np.append(t, k)
                                break
                    rearange_forces[index] = t.astype(int)
                    image_primes = self.descriptor.fingerprintprimes[hash_name]
                    # scaling of fingerprint derivatives to be consistent with
                    # fingerprint scaling.
                    _image_primes = copy.copy(image_primes)
                    for _, key in enumerate(list(image_primes.keys())):
                        base_atom = key[3]
                        fprange_atom = np.array(fprange[base_atom])
                        fprange_dif = fprange_atom[:, 1] - fprange_atom[:, 0]
                        fprange_dif[fprange_dif < 10.0 ** (-8.0)] = 2
                        fprime = np.array(image_primes[key])
                        fprime = 2 * fprime / fprange_dif
                        _image_primes[key] = fprime

                    image_prime_values = list(_image_primes.values())
                    image_prime_keys = list(_image_primes.keys())
                    fp_length = len(image_fingerprint[0][1])
                    num_atoms = len(image_fingerprint)
                    fingerprintprimes = torch.zeros(
                        fp_length * num_atoms, 3 * num_atoms
                    )
                    for idx, fp_key in enumerate(image_prime_keys):
                        image_prime = torch.tensor(image_prime_values[idx])
                        base_atom = fp_key[2]
                        wrt_atom = fp_key[0]
                        coord = fp_key[4]
                        fingerprintprimes[
                            base_atom * fp_length : base_atom * fp_length + fp_length,
                            wrt_atom * 3 + coord,
                        ] = image_prime
                    # store primes in a sparse matrix format
                    if self.store_primes:
                        sp_matrix = sparse.coo_matrix(fingerprintprimes)
                        sparse.save_npz(
                            open("./stored-primes/" + hash_name, "wb"), sp_matrix
                        )
                    fprimes_dataset.append(fingerprintprimes)
                forces_dataset.append(torch.from_numpy(image_forces))
        if self.delta:
            self.delta_energies /= num_of_atoms
            target_ref_per_atom = energy_dataset[0]
            delta_ref_per_atom = self.delta_energies[0]
            relative_targets = energy_dataset - target_ref_per_atom
            relative_delta = self.delta_energies - delta_ref_per_atom
            energy_dataset = torch.FloatTensor(relative_targets - relative_delta)
            scalings = [target_ref_per_atom, delta_ref_per_atom]
        else:
            energy_dataset = torch.FloatTensor(energy_dataset)
            scalings = [0, 0]
        scale = Transform(energy_dataset)
        energy_dataset = scale.norm(energy_dataset)
        if self.forcetraining:
            for idx, force in enumerate(forces_dataset):
                forces_dataset[idx] = scale.norm(force, energy=False)
        scalings.append(scale)

        return (
            fingerprint_dataset,
            energy_dataset,
            num_of_atoms,
            fprimes_dataset,
            forces_dataset,
            index_hashes,
            scalings,
            rearange_forces,
        )

    def __getitem__(self, index):
        fingerprint = self.fingerprint_dataset[index]
        energy = self.energy_dataset[index]
        idx_hash = self.index_hashes[index]
        rearange = None
        fprime = None
        forces = None
        if self.forcetraining:
            if self.store_primes:
                fprime = sparse.load_npz(open("./stored-primes/" + idx_hash, "rb"))
                fprime = torch.tensor(fprime.toarray())
            else:
                fprime = self.sparse_fprimes[index]
            forces = self.forces_dataset[index]
            rearange = self.rearange_forces[index]
        return [fingerprint, energy, fprime, forces, self.scalings, rearange]

    def unique(self):
        """Returns the unique elements contained in the training dataset"""
        elements = np.array(
            [atom.symbol for atoms in self.atom_images for atom in atoms]
        )
        _, idx = np.unique(elements, return_index=True)
        elements = list(elements[np.sort(idx)])
        return elements

    def fp_length(self):
        """Computes the fingerprint length of the training images"""
        return len(list(self.fprange.values())[0])

    def create_splits(self, training_data, val_frac, resample=None):
        """Constructs a training and validation sampler to be utilized for
        constructing training and validation sets."""
        dataset_size = len(training_data)
        if resample:
            resample = len(resample)
            dataset_size -= resample
        indices = np.random.permutation(dataset_size)
        split = int(val_frac * dataset_size)
        train_idx, val_idx = indices[split:], indices[:split]
        if resample:
            sample_idx = np.random.permutation(
                np.arange(dataset_size, dataset_size + resample)
            )
            split = int(val_frac * len(sample_idx))
            train_idx = np.concatenate((train_idx, sample_idx[split:]))
            val_idx = np.concatenate((val_idx, sample_idx[:split]))
        train_sampler = SubsetRandomSampler(train_idx)
        val_sampler = SubsetRandomSampler(val_idx)

        samplers = {"train": train_sampler, "val": val_sampler}

        return samplers


def make_sparse(primes):
    primes = primes.to_sparse()
    return primes


def factorize_data(training_data):
    """
    Factorizes the dataset into separate lists.

    1. unique_atoms = Identifies the unique elements in the dataset

    2. fingerprint_dataset = Extracts the fingerprints for each hashed data
    sample in the dataset

    3. energy_dataset = Extracts the ab initio potential energy for each hashed data
    sample in the dataset

    4. num_of atoms = Number of atoms per image

    5. sparse_fprimes = fingerprint derivatives for each atom across all
    images. Dimensions of said tensor is QxPx3:
        Q - Atoms in batch
        P - Length of fingerprints
        3 - x,y,z components of the fingerprint derivatives. It should be noted
        that for each image containing N atoms, there are 3N directional
        components for each atom. For each atom, N contributions are summed for
        each directional component to arrive at a dimension of 3 for each atom.

    6. image_forces = Extracts the ab initio forces for each hashed data sample in the
    dataset.

    """
    forcetraining = False
    if training_data[0][2] is not None:
        forcetraining = True
    # scalings = training_data[0][-1]
    scalings = training_data[0][-2]
    unique_atoms = []
    fingerprint_dataset = []
    energy_dataset = []
    num_of_atoms = []
    if forcetraining:
        total_entries = 0
        previous_entries = 0
    for image in training_data:
        num_atom = float(len(image[0]))
        num_of_atoms.append(num_atom)
        if forcetraining:
            # track the number of entries in the fprimes matrix
            image[2] = make_sparse(image[2])
            total_entries += len(image[2]._values())

    image_forces = torch.tensor([])
    sparse_fprimes = torch.tensor([])
    # Construct a sparse matrix with dimensions PQx3Q, if forcetraining is on.
    if forcetraining:
        image_forces = []
        dim1_start = 0
        dim2_start = 0
        # pre-define matrices filled with zeros
        fprimes_inds = torch.zeros((2, total_entries), dtype=torch.int64)
        fprimes_vals = torch.zeros((total_entries))

    rearange_set = np.zeros((1, int(sum(num_of_atoms))))[0]
    atom_shift = 0
    for idx, image in enumerate(training_data):
        image_fingerprint = image[0]
        fingerprint_dataset.append(image_fingerprint)
        for atom in image_fingerprint:
            element = atom[0]
            if element not in unique_atoms:
                unique_atoms.append(element)
        image_potential_energy = image[1]
        energy_dataset.append(image_potential_energy)
        if forcetraining:
            rearange_set[atom_shift : atom_shift + len(image_fingerprint)] = (
                image[-1] + atom_shift
            )
            atom_shift += len(image_fingerprint)
            fprime = image[2]
            # build the matrix of indices
            # the indices need to be offset by dim1 and dim2
            if not fprime.is_sparse:
                fprime = fprime.to_sparse()
            dim1 = fprime.shape[0]
            dim2 = fprime.shape[1]
            s_fprime_inds = fprime._indices() + torch.LongTensor(
                [[dim1_start], [dim2_start]]
            )
            # build the matrix of values
            s_fprime_vals = fprime._values().type(torch.FloatTensor)
            num_entries = len(s_fprime_vals)
            # fill in the entries
            fprimes_inds[
                :, previous_entries : previous_entries + num_entries
            ] = s_fprime_inds
            fprimes_vals[
                previous_entries : previous_entries + num_entries
            ] = s_fprime_vals
            previous_entries += num_entries
            dim1_start += dim1
            dim2_start += dim2
            image_forces.append((image[3]))
    if forcetraining:
        sparse_fprimes = torch.sparse.FloatTensor(
            fprimes_inds, fprimes_vals, torch.Size([dim1_start, dim2_start])
        )
        image_forces = torch.cat(image_forces).float()

    # unique_atoms = sorted(set(unique_atoms))
    unique_atoms = OrderedDict.fromkeys(unique_atoms, 1)
    rearange_set = torch.LongTensor(rearange_set)
    return (
        unique_atoms,
        fingerprint_dataset,
        energy_dataset,
        num_of_atoms,
        sparse_fprimes,
        image_forces,
        scalings,
        rearange_set,
    )


def collate_amp(training_data):
    """
    Reshuffling scheme that reads in raw data and organizes it into element
    specific datasets to be fed into element specific neural networks.
    """
    (
        unique_atoms,
        fingerprint_dataset,
        energy_dataset,
        num_of_atoms,
        fp_primes,
        image_forces,
        scalings,
        rearange,
    ) = factorize_data(training_data)
    batch_size = len(energy_dataset)
    element_specific_fingerprints = {}
    model_input_data = [[], []]
    for element in unique_atoms:
        element_specific_fingerprints[element] = [[], []]
    for fp_index, sample_fingerprints in enumerate(fingerprint_dataset):
        for fingerprint in sample_fingerprints:
            atom_element = fingerprint[0]
            atom_fingerprint = fingerprint[1]
            element_specific_fingerprints[atom_element][0].append(
                torch.tensor(atom_fingerprint, dtype=torch.get_default_dtype())
            )
            element_specific_fingerprints[atom_element][1].append(fp_index)
    for element in unique_atoms:
        element_specific_fingerprints[element][0] = torch.stack(
            element_specific_fingerprints[element][0]
        )
    model_input_data[0].append(element_specific_fingerprints)
    model_input_data[0].append(batch_size)
    model_input_data[0].append(unique_atoms)
    model_input_data[0].append(fp_primes)
    model_input_data[0].append(rearange)
    model_input_data[1].append(torch.tensor(energy_dataset).reshape(-1, 1))
    model_input_data[1].append(torch.FloatTensor(num_of_atoms).reshape(batch_size, 1))
    model_input_data[1].append(image_forces)
    return model_input_data


class TestDataset(Dataset):
    """
    PyTorch inherited abstract class that specifies how testing data is to be preprocessed and
    passed along to the DataLoader.

    Parameters:
    -----------
    images: list, file, database
        Training data to be utilized for the regression model.

    descriptor: object
        Scheme to be utilized for computing fingerprints

    fprange: Dict
        Fingerprint ranges of the training dataset to be used to scale the test
        dataset fingerprints in the same manner.

    """

    def __init__(
        self, images, unique_atoms, descriptor, Gs, fprange, label="example", cores=1
    ):
        self.images = images
        if type(images) is not list:
            self.images = [images]
        self.descriptor = descriptor
        self.atom_images = self.images
        if isinstance(images, str):
            extension = os.path.splitext(images)[1]
            if extension != (".traj" or ".db"):
                self.atom_images = ase.io.read(images, ":")
        self.fprange = fprange
        self.training_unique_atoms = unique_atoms
        self.hashed_images = amp_hash(self.atom_images)
        if descriptor == SNN_Gaussian:
            self.hashed_images = hash_images(self.atom_images, Gs)
            self.fps, self.fp_primes = make_amp_descriptors_simple_nn(
                self.atom_images,
                Gs,
                self.training_unique_atoms,
                cores=cores,
                label=label,
                save=False,
            )
        self.unique_atoms = self.unique()

    def __len__(self):
        return len(self.hashed_images)

    def __getitem__(self, index):
        image_fingerprint = self.fps
        fprange = self.fprange
        atom_order = []
        # fingerprint scaling to a range of [-1,1].
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
            atom_order.append(atom)
        image_primes = self.fp_primes
        # image_primes = self.descriptor.fingerprintprimes[hash_name]
        # fingerprint derivative scaling to a range of [0,1].
        _image_primes = copy.copy(image_primes)
        for _, key in enumerate(list(image_primes.keys())):
            base_atom = key[3]
            fprange_atom = fprange[base_atom]
            fprime = image_primes[key]
            for i in range(len(fprime)):
                if (fprange_atom[i][1] - fprange_atom[i][0]) > (10.0 ** (-8.0)):
                    fprime[i] = 2.0 * (
                        fprime[i] / (fprange_atom[i][1] - fprange_atom[i][0])
                    )
            _image_primes[key] = fprime
        prime_mapping = []
        for element in self.unique_atoms:
            indices = [i for i, x in enumerate(atom_order) if x == element]
            prime_mapping += indices
        new_order = [atom_order[i] for i in prime_mapping]
        used = set()
        rearange = []
        for i, x in enumerate(atom_order):
            for k, l in enumerate(new_order):
                if (x == l) and (k not in used):
                    used.add(k)
                    rearange.append(k)
                    break

        image_prime_values = list(_image_primes.values())
        image_prime_keys = list(_image_primes.keys())
        fp_length = len(image_fingerprint[0][1])
        num_atoms = len(image_fingerprint)
        fingerprintprimes = torch.zeros(fp_length * num_atoms, 3 * num_atoms)
        for idx, fp_key in enumerate(image_prime_keys):
            image_prime = torch.tensor(image_prime_values[idx])
            base_atom = fp_key[2]
            wrt_atom = fp_key[0]
            coord = fp_key[4]
            fingerprintprimes[
                base_atom * fp_length : base_atom * fp_length + fp_length,
                wrt_atom * 3 + coord,
            ] = image_prime

        return [image_fingerprint, fingerprintprimes, num_atoms, rearange]

    def unique(self):
        elements = np.array(
            [atom.symbol for atoms in self.atom_images for atom in atoms]
        )
        _, idx = np.unique(elements, return_index=True)
        elements = list(elements[np.sort(idx)])
        return elements

    def fp_length(self):
        self.fp_length = len(list(self.fprange.values())[0])
        return len(list(self.fprange.values())[0])

    def collate_test(self, training_data):
        """
        Reshuffling scheme that reads in raw data and organizes it into element
        specific datasets to be fed into the element specific Neural Nets.
        """
        fingerprint_dataset = []
        rearange_set = np.array([])
        num_of_atoms = []
        total_entries = 0
        previous_entries = 0
        atom_shift = 0
        for image in training_data:
            image[1] = image[1].to_sparse()  # presparify the fprimes
            total_entries += len(image[1]._values())
            num_of_atoms.append(image[2])
            rearange_set = np.append(rearange_set, np.array(image[-1]) + atom_shift)
            atom_shift += image[2]
        # Construct a sparse matrix with dimensions PQx3Q
        dim1_start = 0
        dim2_start = 0
        # pre-define matrices filled with zeros
        fprimes_inds = torch.zeros((2, total_entries), dtype=torch.int64)
        fprimes_vals = torch.zeros((total_entries))

        for idx, image in enumerate(training_data):
            fprime = image[1]
            dim1 = fprime.shape[0]
            dim2 = fprime.shape[1]
            if not fprime.is_sparse:
                fprime = fprime.to_sparse()
            s_fprime_inds = fprime._indices() + torch.LongTensor(
                [[dim1_start], [dim2_start]]
            )
            s_fprime_vals = fprime._values().type(torch.FloatTensor)
            num_entries = len(s_fprime_vals)
            # fill in the entries
            fprimes_inds[
                :, previous_entries : previous_entries + num_entries
            ] = s_fprime_inds
            fprimes_vals[
                previous_entries : previous_entries + num_entries
            ] = s_fprime_vals
            previous_entries += num_entries

            dim1_start += dim1
            dim2_start += dim2
            image_fingerprint = image[0]
            fingerprint_dataset.append(image_fingerprint)
        sparse_fprimes = torch.sparse.LongTensor(
            fprimes_inds, fprimes_vals, torch.Size([dim1_start, dim2_start])
        )

        element_specific_fingerprints = {}
        model_input_data = []
        for element in self.unique_atoms:
            element_specific_fingerprints[element] = [[], []]
        for fp_index, sample_fingerprints in enumerate(fingerprint_dataset):
            for fingerprint in sample_fingerprints:
                atom_element = fingerprint[0]
                atom_fingerprint = fingerprint[1]
                element_specific_fingerprints[atom_element][0].append(
                    torch.tensor(atom_fingerprint, dtype=torch.get_default_dtype())
                )
                element_specific_fingerprints[atom_element][1].append(fp_index)
        for element in self.unique_atoms:
            element_specific_fingerprints[element][0] = torch.stack(
                element_specific_fingerprints[element][0]
            )
        batch_size = len(num_of_atoms)
        model_input_data.append(element_specific_fingerprints)
        model_input_data.append(batch_size)
        model_input_data.append(self.unique_atoms)
        model_input_data.append(num_of_atoms)
        model_input_data.append(sparse_fprimes)
        model_input_data.append(torch.LongTensor(rearange_set))

        return model_input_data
