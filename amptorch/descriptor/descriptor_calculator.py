from .base_descriptor import BaseDescriptor
from ase import Atoms
from scipy.sparse import coo_matrix, vstack
import numpy as np


class DescriptorCalculator:
    def __init__(
        self,
        images,
        descriptor,
        automatic_calculation=True,
        calculate_descriptor_primes=True,
        sparse_prime=True,
        store_descriptors=True,
        training_data=False,
        parallel=False,
        cores=1,
    ):
        assert isinstance(descriptor, BaseDescriptor)

        self.images = images
        self.descriptor = descriptor
        self.automatic_calculation = automatic_calculation
        self.calculate_descriptor_primes = calculate_descriptor_primes
        self.sparse_prime = sparse_prime
        self.store_descriptors = store_descriptors
        self.training_data = training_data

        self.element_list = self.descriptor._get_element_list()
        self.descriptors_ready = False

        if self.automatic_calculation:
            self.prepare_descriptors()

    def prepare_descriptors(self):
        self.calculated_decsriptor_list = self.descriptor.prepare_fingerprints(
            self.images,
            parallel=None,
            log=None,
            calculate_derivatives=self.calculate_descriptor_primes,
            save=self.store_descriptors,
            get_training_data=self.training_data,
        )

        self.descriptors_ready = True

    def get_descriptors(self, separate_atomtypes=True):
        if not self.descriptors_ready:
            print(
                "ERROR, descriptors not calculated yet, please call prepare_descriptors() function first"
            )
            return None

        if separate_atomtypes:
            result = {}
            for element in self.element_list:
                element_descriptor_list = []
                for calculated_decsriptor in self.calculated_decsriptor_list:
                    if element in calculated_decsriptor.keys():
                        temp = calculated_decsriptor[element]["descriptors"].copy()
                        element_descriptor_list.append(temp)
                result[element] = element_descriptor_list
            return result

        else:
            print(
                "WARNING: atomtype separation turned off, please make sure the dimensions match and you know what you are doing"
            )
            result = []
            for calculated_decsriptor in self.calculated_decsriptor_list:
                descriptors = np.array([])
                for element in self.element_list:
                    if element in calculated_decsriptor.keys():
                        temp = calculated_decsriptor[element]["descriptors"].copy()
                        descriptors = (
                            np.vstack([descriptors, temp]) if descriptors.size else temp
                        )
                        descriptors.append(temp)
                result.append(descriptors)
            return result

    def _get_calculated_descriptors(self):
        return self.calculated_decsriptor_list

    # TODO
    def calculate_PCA(
        self, separate_atomtypes=True, save_models=True, n_components=10, apply_PCA=True
    ):
        raise NotImplementedError

    # TODO
    def calculate_scaling(
        self,
        separate_atomtypes=True,
        save_models=True,
        scaler_min=-1,
        scaler_max=1,
        apply_scaling=True,
    ):
        raise NotImplementedError
