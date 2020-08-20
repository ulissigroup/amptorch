from torch.utils.data import Dataset
from .descriptor_base import AMPTorchDescriptorBase
from .descriptor_calculator import DescriptorCalculator
from ase import Atoms

class AMPTorchDataset(Dataset):
    def __init__(
        self,
        trajs,
        descriptor,
        automatic_calculation = True,
        calculate_descriptor_primes = True,
        sparse_prime = True,
        store_descriptors = True,
        training_data = False,
        result_dir = "./_results_/",
        parallel = False,
        cores = 1,
    ):
        self.descriptor_calculator = DescriptorCalculator(
            trajs = trajs,
            descriptor = descriptor,
            automatic_calculation=automatic_calculation,
            calculate_descriptor_primes=calculate_descriptor_primes,
            sparse_prime=sparse_prime,
            store_descriptors=store_descriptors,
            training_data=training_data,
            result_dir=result_dir,
            parallel=parallel,
            cores=cores
        )

    def __len__(self):
        pass

    def __getitem__(self, index):
        pass