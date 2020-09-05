import itertools

import numpy as np
import torch
from torch.utils.data import Dataset
from torch_geometric.data import Batch, Data

from .descriptor.descriptor_calculator import DescriptorCalculator


class AtomsDataset(Dataset):
    def __init__(
        self, images, descriptor, forcetraining=True, save_fps=True, cores=1,
    ):

        self.images = images
        self.forcetraining = forcetraining

        self.descriptor_calculator = DescriptorCalculator(
            images=images,
            descriptor=descriptor,
            calc_derviatives=forcetraining,
            save_fps=save_fps,
            cores=cores,
        )

        self.data_list = self.process()

    def process(self):
        descriptor_data = self.descriptor_calculator.prepare_descriptors()
        data_list = []

        for idx, image in enumerate(self.images):
            image_data = descriptor_data[idx]
            energy = image.get_potential_energy(apply_constraint=False)
            forces = image.get_forces(apply_constraint=False)
            atomic_numbers = image.get_atomic_numbers()
            natoms = len(image)
            image_fingerprint = image_data["descriptors"]

            image_fingerprint = torch.tensor(image_fingerprint)
            atomic_numbers = torch.tensor(atomic_numbers)
            image_idx = torch.full((1, natoms), idx, dtype=torch.int64).view(-1)

            data = Data(
                fingerprint=image_fingerprint,
                image_idx=image_idx,
                atomic_numbers=atomic_numbers,
                energy=energy,
                natoms=natoms,
            )
            if self.forcetraining:
                fp_prime_val = image_data["descriptor_primes"]["val"]
                fp_prime_row = image_data["descriptor_primes"]["row"]
                fp_prime_col = image_data["descriptor_primes"]["col"]
                fp_prime_size = image_data["descriptor_primes"]["size"]

                indices = np.vstack((fp_prime_row, fp_prime_col))
                torch_indices = torch.LongTensor(indices)
                torch_values = torch.FloatTensor(fp_prime_val)
                fp_primes = torch.sparse.FloatTensor(
                    torch_indices, torch_values, torch.Size(fp_prime_size)
                )

                data.forces = torch.FloatTensor(forces)
                data.fprimes = fp_primes

            data_list.append(data)

        return data_list

    def __len__(self):
        return len(self.data_list)

    def __getitem__(self, index):
        return self.data[index]


# Adapted from https://github.com/pytorch/pytorch/issues/31942
def sparse_block_diag(arrs):
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
                list(
                    itertools.product(
                        np.arange(r, r + rr.numpy()), np.arange(c, c + cc.numpy())
                    )
                )
            ).t()
        ]
        v += [arrs[k].to_dense().flatten()]
        r += rr.numpy()
        c += cc.numpy()
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
