import copy
import os
import numpy as np
import torch
from torch.utils.data import Dataset
from collections import OrderedDict
import itertools
import ase
from amptorch.gaussian import make_symmetry_functions, SNN_Gaussian
from amptorch.data_utils import Transform
from amptorch.utils import (
    make_amp_descriptors_simple_nn,
    calculate_fingerprints_range,
    hash_images,
)
from amp.utilities import hash_images as amp_hash
from amp.utilities import get_hash as get_amp_hash

from torch_geometric.data import Data, Batch
from tqdm import tqdm

# from scipy.sparse import coo_matrix, block_diag

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
        # wrap all fingerprinting work into its own module

        print("Calculating fingerprints...")
        G2_etas = Gs["G2_etas"]
        G2_rs_s = Gs["G2_rs_s"]
        G4_etas = Gs["G4_etas"]
        G4_zetas = Gs["G4_zetas"]
        G4_gammas = Gs["G4_gammas"]
        cutoff = Gs["cutoff"]

        # create fingerprints
        self.hashed_images = amp_hash(self.atom_images)
        G = make_symmetry_functions(elements=self.elements, type="G2", etas=G2_etas)
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
        self.fprange = calculate_fingerprints_range(self.descriptor, self.hashed_images)

        # perform preprocessing
        self.dataset = self.process()

    def __len__(self):
        return len(self.atom_images)

    def process(self):
        data_list = []
        for idx, atoms in tqdm(enumerate(self.images)):
            hash_name = get_amp_hash(atoms)

            # scale fingerprints
            image_fingerprint = self.descriptor.fingerprints[hash_name]
            for i, (atom, afp) in enumerate(image_fingerprint):
                _afp = copy.copy(afp)
                fprange_atom = np.array(self.fprange[atom])
                for _ in range(np.shape(_afp)[0]):
                    if (fprange_atom[_][1] - fprange_atom[_][0]) > (10.0 ** (-8.0)):
                        _afp[_] = -1 + 2.0 * (
                            (_afp[_] - fprange_atom[_][0])
                            / (fprange_atom[_][1] - fprange_atom[_][0])
                        )
                image_fingerprint[i] = (atom, _afp)

            # convert amp structure
            fp_length = len(image_fingerprint[0][1])
            natoms = len(image_fingerprint)
            image_fingerprint = torch.tensor(
                np.vstack([fp[1] for fp in image_fingerprint])
            )
            potential_energy = self.hashed_images[hash_name].get_potential_energy(
                apply_constraint=False
            )
            atomic_numbers = torch.tensor(atoms.get_atomic_numbers())
            image_idx = torch.full((1, natoms), idx, dtype=torch.int64).view(-1)

            data = Data(
                fingerprint=image_fingerprint,
                image_idx=image_idx,
                atomic_numbers=atomic_numbers,
                energy=potential_energy,
                natoms=natoms,
            )

            if self.forcetraining:
                image_forces = self.hashed_images[hash_name].get_forces(
                    apply_constraint=False
                )
                # scale fingerprintprimes
                image_primes = self.descriptor.fingerprintprimes[hash_name]
                _image_primes = copy.copy(image_primes)
                for _, key in enumerate(list(image_primes.keys())):
                    base_atom = key[3]
                    fprange_atom = np.array(self.fprange[base_atom])
                    fprange_dif = fprange_atom[:, 1] - fprange_atom[:, 0]
                    fprange_dif[fprange_dif < 10.0 ** (-8.0)] = 2
                    fprime = np.array(image_primes[key])
                    fprime = 2 * fprime / fprange_dif
                    _image_primes[key] = fprime

                image_prime_values = list(_image_primes.values())
                image_prime_keys = list(_image_primes.keys())
                # convert amp fprime format to matrix
                fingerprintprimes = torch.zeros(fp_length * natoms, 3 * natoms)
                for idx, fp_key in enumerate(image_prime_keys):
                    image_prime = torch.tensor(image_prime_values[idx])
                    base_atom = fp_key[2]
                    wrt_atom = fp_key[0]
                    coord = fp_key[4]
                    fingerprintprimes[
                        base_atom * fp_length : base_atom * fp_length + fp_length,
                        wrt_atom * 3 + coord,
                    ] = image_prime

                data.forces = torch.FloatTensor(image_forces)
                data.fprimes = fingerprintprimes

            data_list.append(data)
            # write dataset as *.pt file for future use
        return data_list

    def unique(self):
        elements = np.array(
            [atom.symbol for atoms in self.atom_images for atom in atoms]
        )
        _, idx = np.unique(elements, return_index=True)
        elements = list(elements[np.sort(idx)])
        return elements

    def __getitem__(self, index):
        return self.dataset[index]


def sparse_block_diag(arrs):
    # Adapted from https://github.com/pytorch/pytorch/issues/31942
    bad_args = [
        k
        for k in range(len(arrs))
        if not (isinstance(arrs[k], torch.Tensor) and arrs[k].ndim == 2)
    ]
    if bad_args:
        raise ValueError(
            "arguments in the following positions must be 2-dimension tensor: %s"
            % bad_args
        )

    shapes = torch.tensor([a.shape for a in arrs])

    i = []
    v = []
    r, c = 0, 0
    for k, (rr, cc) in enumerate(shapes):
        i += [
            torch.LongTensor(
                list(itertools.product(np.arange(r, r+rr), np.arange(c, c+cc)))
            ).t()
        ]
        v += [arrs[k].flatten()]
        r += rr
        c += cc
    out = torch.sparse.DoubleTensor(
        torch.cat(i, dim=1), torch.cat(v), torch.sum(shapes, dim=0).tolist()
    )
    return out


def collate_amp(data_list):
    mtxs = []
    for data in data_list:
        mtxs.append(data.fprimes)
        data.fprimes = None
    batch = Batch.from_data_list(data_list)
    for i, data in enumerate(data_list):
        data.fprimes = mtxs[i]
    block_matrix = sparse_block_diag(mtxs)
    batch.fprimes = block_matrix
    return batch, (batch.energy, batch.forces)
