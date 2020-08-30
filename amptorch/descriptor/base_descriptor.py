from abc import ABC, abstractmethod
import numpy as np
import h5py

import os
import time
from .util import get_hash, list_symbols_to_indices, list_indices_to_symbols


class BaseDescriptor(ABC):
    def __init__(self):
        super().__init__()
        self.fp_database = "processed/descriptors/"

        # To Be specified/calculated
        self.descriptor_type = "default"
        self.descriptor_setup_hash = "default"

        self.elements = []

    @abstractmethod
    def calculate_fingerprints(self, image, params_set, calculate_derivatives=True):
        # image is a single snapshot
        pass

    @abstractmethod
    def get_descriptor_setup_hash(self):
        # set self.descriptor_setup_hash
        pass

    @abstractmethod
    def save_descriptor_setup(self, filename):
        pass

    @abstractmethod
    def prepare_descriptor_parameters(self):
        # prepare self.params_set
        pass

    def prepare_fingerprints(
        self,
        images,
        parallel=False,
        log=None,
        calculate_derivatives=True,
        save=True,
        get_training_data=False,
    ):
        images_descriptor_list = []

        # if save is true, create directories if not exist
        self._setup_fingerprint_database(save=save)

        for image in images:
            image_hash = get_hash(image)
            image_db_filename = "{}/{}.h5".format(
                self.desc_fp_database_dir,
                image_hash,
            )

            # if save, then read/write from db as needed
            if save:
                temp_descriptor_list = self._compute_fingerprints(
                    image,
                    image_db_filename,
                    parallel=parallel,
                    log=log,
                    calculate_derivatives=calculate_derivatives,
                    save=True,
                    get_training_data=get_training_data,
                )

            # if not save, but db exist, read from db if possible
            elif os.path.exists(self.desc_type_database_dir):
                temp_descriptor_list = self._compute_fingerprints(
                    image,
                    image_db_filename,
                    parallel=parallel,
                    log=log,
                    calculate_derivatives=calculate_derivatives,
                    save=False,
                    get_training_data=get_training_data,
                )

            # if not save, and db not exist, calculat fps on the fly
            else:
                temp_descriptor_list = self._prepare_fingerprints_single_traj_nodb(
                    image,
                    image_db_filename,
                    parallel=parallel,
                    log=log,
                    calculate_derivatives=calculate_derivatives,
                    get_training_data=get_training_data,
                )

            images_descriptor_list += temp_descriptor_list

        return images_descriptor_list

    def _compute_fingerprints(
        self,
        image,
        image_db_filename,
        parallel=False,
        log=None,
        calculate_derivatives=True,
        save=True,
        get_training_data=False,
    ):
        descriptor_list = []

        with h5py.File(image_db_filename, "a") as db:
            image_dict = {}

            symbol_arr = np.array(image.get_chemical_symbols())
            image_dict["atomic_numbers"] = list_symbols_to_indices(symbol_arr)
            num_atoms = len(symbol_arr)
            image_dict["num_atoms"] = num_atoms

            try:
                current_snapshot_grp = db[str(0)]
            except Exception:
                current_snapshot_grp = db.create_group(str(0))

            num_desc_list = []
            index_arr_dict = {}
            num_desc_dict = {}
            fp_dict = {}

            fp_prime_val_dict = {}
            fp_prime_row_dict = {}
            fp_prime_col_dict = {}
            fp_prime_size_dict = {}

            for element in self.elements:
                if element in image.get_chemical_symbols():
                    index_arr = np.arange(num_atoms)[symbol_arr == element]
                    index_arr_dict[element] = index_arr

                    try:
                        current_element_grp = current_snapshot_grp[element]
                    except:
                        current_element_grp = current_snapshot_grp.create_group(element)

                    if calculate_derivatives:
                        try:
                            size_info = np.array(current_element_grp["size_info"])
                            fps = np.array(current_element_grp["fps"])
                            fp_primes_val = np.array(
                                current_element_grp["fp_primes_val"]
                            )
                            fp_primes_row = np.array(
                                current_element_grp["fp_primes_row"]
                            )
                            fp_primes_col = np.array(
                                current_element_grp["fp_primes_col"]
                            )
                            fp_primes_size = np.array(
                                current_element_grp["fp_primes_size"]
                            )
                        except Exception:
                            (
                                size_info,
                                fps,
                                fp_primes_val,
                                fp_primes_row,
                                fp_primes_col,
                                fp_primes_size,
                            ) = self.calculate_fingerprints(
                                image,
                                element,
                                calculate_derivatives=calculate_derivatives,
                            )

                            if save:
                                current_element_grp.create_dataset(
                                    "size_info", data=size_info
                                )
                                current_element_grp.create_dataset(
                                    "fps", data=fps
                                )
                                current_element_grp.create_dataset(
                                    "fp_primes_val", data=fp_primes_val
                                )
                                current_element_grp.create_dataset(
                                    "fp_primes_row", data=fp_primes_row
                                )
                                current_element_grp.create_dataset(
                                    "fp_primes_col", data=fp_primes_col
                                )
                                current_element_grp.create_dataset(
                                    "fp_primes_size", data=fp_primes_size
                                )

                        num_desc_list.append(size_info[2])
                        num_desc_dict[element] = size_info[2]
                        fp_dict[element] = fps

                        fp_prime_val_dict[element] = fp_primes_val
                        fp_prime_row_dict[element] = fp_primes_row
                        fp_prime_col_dict[element] = fp_primes_col
                        fp_prime_size_dict[element] = fp_primes_size

                    else:
                        try:
                            size_info = np.array(current_element_grp["size_info"])
                            fps = np.array(current_element_grp["fps"])
                        except:
                            size_info, fps, _, _, _, _ = self.calculate_fingerprints(
                                image,
                                element,
                                calculate_derivatives=calculate_derivatives,
                            )

                            if save:
                                current_element_grp.create_dataset(
                                    "fps", data=fps
                                )

                        num_desc_list.append(size_info[2])
                        num_desc_dict[element] = size_info[2]
                        fp_dict[element] = fps

                else:
                    print("element not in current image: {}".format(element))

            num_desc_max = np.max(num_desc_list)
            image_fp_array = np.zeros((num_atoms, num_desc_max))
            for element in fp_dict.keys():
                image_fp_array[
                    index_arr_dict[element], : num_desc_dict[element]
                ] = fp_dict[element]

            print("fp size: {}".format(image_fp_array.shape))
            image_dict["descriptors"] = image_fp_array
            image_dict["num_descriptors"] = num_desc_dict

            if calculate_derivatives:
                descriptor_prime_dict = {}
                descriptor_prime_dict["size"] = np.array(
                    [num_desc_max * num_atoms, 3 * num_atoms]
                )
                descriptor_prime_row_list = []
                descriptor_prime_col_list = []
                descriptor_prime_val_list = []
                for element in fp_prime_val_dict.keys():
                    descriptor_prime_row_list.append(
                        self._fp_prime_element_row_index_to_image_row_index(
                            fp_prime_row_dict[element],
                            index_arr_dict[element],
                            num_desc_dict[element],
                            num_desc_max,
                        )
                    )
                    descriptor_prime_col_list.append(
                        fp_prime_col_dict[element]
                    )
                    descriptor_prime_val_list.append(
                        fp_prime_val_dict[element]
                    )
                descriptor_prime_dict["row"] = np.concatenate(
                    descriptor_prime_row_list
                )
                descriptor_prime_dict["col"] = np.concatenate(
                    descriptor_prime_col_list
                )
                descriptor_prime_dict["val"] = np.concatenate(
                    descriptor_prime_val_list
                )
                image_dict["descriptor_primes"] = descriptor_prime_dict

            descriptor_list.append(image_dict)

        return descriptor_list

    def _prepare_fingerprints_single_traj_nodb(
        self,
        traj,
        traj_db_filename,
        parallel=False,
        log=None,
        calculate_derivatives=True,
        get_training_data=False,
    ):
        descriptor_list = []

        Total_Num_Snapshots = len(traj)
        for i, snapshot in enumerate(traj):
            start_time = time.time()
            image_dict = {}
            if get_training_data:
                image_dict["energy"] = snapshot.get_potential_energy()
                image_dict["forces"] = snapshot.get_forces()

            symbol_arr = snapshot.get_chemical_symbols()
            image_dict["atomic_numbers"] = list_symbols_to_indices(symbol_arr)
            num_atoms = len(symbol_arr)
            image_dict["num_atoms"] = num_atoms

            num_desc_list = []
            index_arr_dict = {}
            num_desc_dict = {}
            fp_dict = {}

            fp_prime_val_dict = {}
            fp_prime_row_dict = {}
            fp_prime_col_dict = {}
            fp_prime_size_dict = {}

            for element in self.elements:
                if element in snapshot.get_chemical_symbols():

                    index_arr = np.arange(num_atoms)[symbol_arr == element]
                    index_arr_dict[element] = index_arr

                    if calculate_derivatives:

                        size_info, fps, fp_primes_val, fp_primes_row, fp_primes_col, fp_primes_size = self.calculate_fingerprints(
                            snapshot,
                            element,
                            calculate_derivatives=calculate_derivatives,
                        )

                        num_desc_list.append(size_info[2])
                        num_desc_dict[element] = size_info[2]
                        fp_dict[element] = fps

                        fp_prime_val_dict[element] = fp_primes_val
                        fp_prime_row_dict[element] = fp_primes_row
                        fp_prime_col_dict[element] = fp_primes_col
                        fp_prime_size_dict[element] = fp_primes_size

                    else:

                        size_info, fps, _, _, _, _ = self.calculate_fingerprints(
                            snapshot,
                            element,
                            calculate_derivatives=calculate_derivatives,
                        )

                        num_desc_list.append(size_info[2])
                        num_desc_dict[element] = size_info[2]
                        fp_dict[element] = fps

                else:
                    print("element not in current snapshot: {}".format(element))

            num_desc_max = np.max(num_desc_list)
            image_fp_array = np.zeros((num_atoms, num_desc_max))
            for element in fp_dict.keys():
                image_fp_array[
                    index_arr_dict[element], : num_desc_dict[element]
                ] = fp_dict[element]

            image_dict["descriptors"] = image_fp_array
            image_dict["num_descriptors"] = num_desc_dict

            if calculate_derivatives:
                descriptor_prime_dict = {}
                descriptor_prime_dict["size"] = np.array(
                    [num_desc_max * num_atoms, 3 * num_atoms]
                )
                descriptor_prime_row_list = []
                descriptor_prime_col_list = []
                descriptor_prime_val_list = []
                for element in fp_prime_val_dict.keys():
                    descriptor_prime_row_list.append(
                        self._fp_prime_element_row_index_to_image_row_index(
                            fp_prime_row_dict[element],
                            index_arr_dict[element],
                            num_desc_dict[element],
                            num_desc_max,
                        )
                    )
                    descriptor_prime_col_list.append(fp_prime_col_dict[element])
                    descriptor_prime_val_list.append(fp_prime_val_dict[element])
                descriptor_prime_dict["row"] = np.concatenate(descriptor_prime_row_list)
                descriptor_prime_dict["col"] = np.concatenate(descriptor_prime_col_list)
                descriptor_prime_dict["val"] = np.concatenate(descriptor_prime_val_list)
                image_dict["descriptor_primes"] = descriptor_prime_dict

            descriptor_list.append(image_dict)

            took_time = time.time() - start_time
            print(
                "finished snapshot {}/{}, took time: {}".format(
                    i + 1, Total_Num_Snapshots, took_time
                )
            )

        return descriptor_list

    def _fp_prime_element_row_index_to_image_row_index(
        self, original_rows, index_arr, num_desc, num_desc_max
    ):
        atom_indices_for_specific_element, desc_indices = np.divmod(
            original_rows, num_desc
        )

        atom_indices_in_image = index_arr[atom_indices_for_specific_element]

        new_row = atom_indices_in_image * num_desc_max + desc_indices
        return new_row

    def _setup_fingerprint_database(self, save=False):
        self.get_descriptor_setup_hash()
        self.desc_type_database_dir = "{}/{}".format(
            self.fp_database, self.descriptor_type
        )

        self.desc_fp_database_dir = "{}/{}".format(
            self.desc_type_database_dir, self.descriptor_setup_hash
        )

        if save:
            os.makedirs(self.fp_database, exist_ok=True)
            os.makedirs(self.desc_type_database_dir, exist_ok=True)
            os.makedirs(self.desc_fp_database_dir, exist_ok=True)
            descriptor_setup_filename = "descriptor_log.txt"
            descriptor_setup_path = "{}/{}".format(
                self.desc_fp_database_dir, descriptor_setup_filename
            )
            self.save_descriptor_setup(descriptor_setup_path)

    def _get_element_list(self):
        return self.elements
