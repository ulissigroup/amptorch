from torch.utils.data import Dataset
from .amptorch_descriptor.descriptor_base import AMPTorchDescriptorBase
from .amptorch_descriptor.descriptor_calculator import DescriptorCalculator
from .amptorch_descriptor.constants import ATOM_INDEX_TO_SYMBOL_DICT, ATOM_SYMBOL_TO_INDEX_DICT
from ase import Atoms
from torch_geometric.data import Data, Batch
from scipy.sparse import coo_matrix, vstack

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

    def _prepare_data(self):
        data_list = []
        raw_data = self.descriptor_calculator._get_calculated_descriptors() 
        element_list = self.descriptor_calculator.element_list
        
        for idx, image_data in enumerate(raw_data):
            potential_energy = image_data["energy"]
            image_fp_list = []
            atomic_numbers = []
            natoms = 0

            for element in element_list:
                atomic_number = ATOM_SYMBOL_TO_INDEX_DICT[element]
                size = image_data[element]["size_info"]
                num_element = size_info[0]
                atomic_numbers += [atomic_number] * num_element
                natoms += num_element
                image_fp_list.append(image_data[element]["descriptors"])

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
                    forces = image_data[element]["forces"]
                    fp_prime_value = image_data[element]["descriptor_primes"]["value"]
                    fp_prime_row   = image_data[element]["descriptor_primes"]["row"]
                    fp_prime_col   = image_data[element]["descriptor_primes"]["col"]
                    fp_prime_size  = image_data[element]["descriptor_primes"]["size"]

                    element_fp_prime_matrix = coo_matrix((fp_prime_value, (fp_prime_row, fp_prime_col)), shape=fp_prime_size)
                    image_fp_primes_list.append(element_fp_prime_matrix)
                    
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

        self.data_length = len(raw_data)
        self.data = data_list

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