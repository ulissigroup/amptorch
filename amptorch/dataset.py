import numpy as np
import torch
from torch.utils.data import Dataset
from torch_geometric.data import Batch, Data

from amptorch.data_utils import Normalize, sparse_block_diag
from amptorch.descriptor.descriptor_calculator import DescriptorCalculator


class AtomsDataset(Dataset):
    def __init__(
        self, images, descriptor, forcetraining=True, save_fps=True, cores=1,
    ):
        self.images = images
        self.forcetraining = forcetraining

        self.descriptor_calculator = DescriptorCalculator(
            images=images,
            descriptor=descriptor,
            calc_derivatives=forcetraining,
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

            image_fingerprint = torch.FloatTensor(image_fingerprint)
            atomic_numbers = torch.LongTensor(atomic_numbers)
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

        self.normalizer = Normalize(data_list)
        data_list = self.normalizer.norm(data_list)

        return data_list

    @property
    def input_dim(self):
        return self.data_list[0].fingerprint.shape[1]

    def __len__(self):
        return len(self.data_list)

    def __getitem__(self, index):
        return self.data_list[index]


def data_collater(data_list):
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
