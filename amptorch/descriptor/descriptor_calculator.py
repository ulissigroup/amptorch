import numpy as np

from .base_descriptor import BaseDescriptor


class DescriptorCalculator:
    def __init__(
        self,
        images,
        descriptor,
        calc_derivatives=True,
        save_fps=True,
        verbose=True,
        cores=1,
    ):
        assert isinstance(
            descriptor, BaseDescriptor
        ), "Descriptor must be instance of BaseDescriptor!"

        self.images = images
        self.descriptor = descriptor
        self.calc_derivatives = calc_derivatives
        self.save_fps = save_fps
        self.cores = cores
        self.verbose = verbose

        self.element_list = self.descriptor._get_element_list()
        self.descriptors_ready = False

    def prepare_descriptors(self):
        self.calculated_descriptor_list = self.descriptor.prepare_fingerprints(
            self.images,
            calc_derivatives=self.calc_derivatives,
            save_fps=self.save_fps,
            cores=self.cores,
            verbose=self.verbose,
            log=None,
        )

        self.descriptors_ready = True

        return self.calculated_descriptor_list

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
