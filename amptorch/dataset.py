from torch.utils.data import Dataset
from .amptorch_descriptor.descriptor_base import AMPTorchDescriptorBase
from .amptorch_descriptor.descriptor_calculator import DescriptorCalculator
from .amptorch_descriptor.constants import ATOM_INDEX_TO_SYMBOL_DICT, ATOM_SYMBOL_TO_INDEX_DICT
from ase import Atoms
from torch_geometric.data import Data, Batch
from scipy.sparse import coo_matrix, vstack
import torch
import numpy as np
import itertools

class AMPTorchDataset(Dataset):
    def __init__(
        self,
        trajs,
        descriptor,
        automatic_calculation = True,
        force_calculation = True,
        # calculate_descriptor_primes = True,
        sparse_prime = True,
        store_descriptors = True,
        PCA_transform = False,
        PCA_ncomponents = 10,
        scaling = True,
        training_data = False,
        result_dir = "./_results_/",
        parallel = False,
        cores = 1,
    ):

        self.force_calculation = force_calculation
        self.PCA_transform = PCA_transform
        self.PCA_ncomponents = PCA_ncomponents
        self.scaling = scaling

        self.descriptor_calculator = DescriptorCalculator(
            trajs = trajs,
            descriptor = descriptor,
            automatic_calculation=False,
            calculate_descriptor_primes=self.force_calculation,
            sparse_prime=sparse_prime,
            store_descriptors=store_descriptors,
            training_data=training_data,
            result_dir=result_dir,
            parallel=parallel,
            cores=cores
        )
        
        if automatic_calculation:
            self.process()

    def process(self):
        self.descriptor_calculator.prepare_descriptors()

        if self.PCA_transform:
            self.descriptor_calculator.calculate_PCA(n_components = self.PCA_ncomponents)
        if self.scaling:
            self.descriptor_calculator.calculate_scaling()

        self._prepare_data()

    def _prepare_data(self, match_max_natoms = True):

        data_list = []
        raw_data = self.descriptor_calculator._get_calculated_descriptors() 
        element_list = self.descriptor_calculator.element_list

        if match_max_natoms:
            max_natoms = self._get_max_natoms(raw_data, element_list)

        for idx, image_data in enumerate(raw_data):
            potential_energy = image_data["energy"]
            image_fp_list = []
            atomic_numbers = []
            natoms = 0

            for element in element_list:
                if element in image_data.keys():
                    atomic_number = ATOM_SYMBOL_TO_INDEX_DICT[element]
                    size_info = image_data[element]["size_info"]
                    num_element = size_info[1]
                    atomic_numbers += [atomic_number] * num_element
                    natoms += num_element
                    image_fp_list.append(image_data[element]["descriptors"])
            
            if match_max_natoms:
                num_dummy_atoms = max_natoms - natoms
                if num_dummy_atoms > 0:
                    num_descriptors = (image_fp_list[0].shape)[1]
                    print("added dummy atoms: {} \t num descriptors: {}".format(num_dummy_atoms, num_descriptors))
                    image_fp_list.append(np.zeros(num_dummy_atoms, num_descriptors))
                    atomic_numbers += [0] * num_dummy_atoms
                    natoms += num_dummy_atoms


            image_fingerprint = torch.tensor(
                np.vstack(image_fp_list)
            )
            atomic_numbers = torch.tensor(atomic_numbers)
            image_idx = torch.full((1, natoms), idx, dtype=torch.int64).view(-1)
            data = Data(
                fingerprint=image_fingerprint,
                image_idx=image_idx,
                atomic_numbers=atomic_numbers,
                energy=potential_energy,
                natoms=natoms,
            )

            if self.force_calculation:
                image_forces_list = []
                image_fp_primes_list = []
                for element in element_list:
                    if element in image_data.keys():
                        forces = image_data[element]["forces"]
                        fp_prime_value = image_data[element]["descriptor_primes"]["value"]
                        fp_prime_row   = image_data[element]["descriptor_primes"]["row"]
                        fp_prime_col   = image_data[element]["descriptor_primes"]["col"]
                        fp_prime_size  = image_data[element]["descriptor_primes"]["size"]
                        #size should be [n_atom_select * n_descriptor, 3 * n_atom_actual]

                        if match_max_natoms:
                            # need to change size to [n_atom_select * n_descriptor, 3 * max_atoms]
                            actual_size = np.array([fp_prime_size[0], max_natoms * 3])
                        else:
                            actual_size = fp_prime_size

                        element_fp_prime_matrix = coo_matrix((fp_prime_value, (fp_prime_row, fp_prime_col)), shape=fp_prime_size)

                        image_forces_list.append(forces)
                        image_fp_primes_list.append(element_fp_prime_matrix)
                
                if match_max_natoms and num_dummy_atoms > 0:
                    image_forces_list.append(np.zeros((num_dummy_atoms, 3)))
                    empty_arr = np.array([])
                    dummy_fp_prime_matrix_size = np.array([num_dummy_atoms, max_natoms * 3])
                    dummy_fp_prime_matrix = coo_matrix((empty_arr, (empty_arr, empty_arr)), shape=fp_prime_size)
                    image_fp_primes_list.append(dummy_fp_prime_matrix)
                    
                    
                image_fp_primes = vstack(image_fp_primes_list)
                # scipy_sparse_fp_prime.data, scipy_sparse_fp_prime.row, scipy_sparse_fp_prime.col, np.array(fp_prime.shape)
                indices = np.vstack((image_fp_primes.row, image_fp_primes.col))
                torch_indices = torch.LongTensor(indices)
                torch_values = torch.FloatTensor(image_fp_primes.data)
                fp_primes_torch_sparse = torch.sparse.FloatTensor(torch_indices, torch_values, torch.Size(image_fp_primes.shape))

                image_forces = np.vstack(image_forces_list)

                data.forces = torch.FloatTensor(image_forces)
                data.fprimes = fp_primes_torch_sparse

                # data.forces = torch.FloatTensor(image_forces)
                # data.fprimes = fingerprintprimes
            
            data_list.append(data)

        self.data_length = len(data_list)
        self.data = data_list

    def _get_max_natoms(self, raw_data, element_list):
        max_natoms = 0
        for idx, image_data in enumerate(raw_data):
            natoms = 0
            for element in element_list:
                if element in image_data.keys():
                    size_info = image_data[element]["size_info"]
                    num_element = size_info[1]
                    natoms += num_element
            if max_natoms < natoms:
                max_natoms = natoms
        return max_natoms

    def __len__(self):
        return self.data_length

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
    print(shapes)

    i = []
    v = []
    r, c = 0, 0
    for k, (rr, cc) in enumerate(shapes):
        # print(r, rr, rr.numpy())
        # print(np.arange(r, r+rr.numpy()))
        # print(np.arange(c, c+cc.numpy()))
        # print(itertools.product(np.arange(r, r+rr.numpy()).astype(int), np.arange(c, c+cc.numpy()).astype(int)))
        # print(list(itertools.product(np.arange(r, r+rr.numpy()), np.arange(c, c+cc.numpy()))))
        i += [
            torch.LongTensor(
                list(itertools.product(np.arange(r, r+rr.numpy()), np.arange(c, c+cc.numpy())))
                # list(itertools.product(np.arange(r, r+rr.numpy()).astype(int), np.arange(c, c+cc.numpy()).astype(int)))
                # list(itertools.product(np.arange(r, r+rr).astype(int), np.arange(c, c+cc).astype(int)))
            ).t()
        ]
        # v += [arrs[k].flatten()]
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