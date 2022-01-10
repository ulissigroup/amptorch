import os
from abc import ABC, abstractmethod

import h5py
import numpy as np
from tqdm import tqdm

from .util import get_hash, list_symbols_to_indices, validate_image


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
        self, images, ref_positions_list, calc_derivatives, save_fps, verbose, cores, log
    ):
        images_descriptor_list = []

        # if save is true, create directories if not exist
        self._setup_fingerprint_database(save_fps=save_fps)

        for image, ref_positions in tqdm(
            zip(images, ref_positions_list),
            total=len(images),
            desc="Computing fingerprints",
            disable=not verbose,
        ):
            validate_image(image)
            image_hash = get_hash(image)
            image_db_filename = "{}/{}.h5".format(self.desc_fp_database_dir, image_hash)

            # if save, then read/write from db as needed
            if save_fps:
                try:
                    temp_descriptor_list = self._compute_fingerprints(
                        image,
                        ref_positions,
                        image_db_filename,
                        calc_derivatives=calc_derivatives,
                        save_fps=save_fps,
                        cores=cores,
                        log=log,
                    )
                except Exception:
                    print(
                        "File {} not loaded properly\nProceed to compute in run-time".format(
                            image_db_filename
                        )
                    )
                    temp_descriptor_list = self._compute_fingerprints_nodb(
                        image,
                        ref_positions,
                        calc_derivatives=calc_derivatives,
                        save_fps=save_fps,
                        cores=cores,
                        log=log,
                    )

            # if not save, compute fps on-the-fly
            else:
                temp_descriptor_list = self._compute_fingerprints_nodb(
                    image,
                    ref_positions,
                    calc_derivatives=calc_derivatives,
                    save_fps=save_fps,
                    cores=cores,
                    log=log,
                )

            images_descriptor_list += temp_descriptor_list

        return images_descriptor_list

    def _compute_fingerprints(
        self, image, ref_positions, image_db_filename, calc_derivatives, save_fps, cores, log
    ):
        descriptor_list = []

        with h5py.File(image_db_filename, "a") as db:
            image_dict = {}

            try:
                current_snapshot_grp = db[str(0)]
            except Exception:
                current_snapshot_grp = db.create_group(str(0))


            if calc_derivatives:
                try:
                    size_info = np.array(current_snapshot_grp["size_info"])
                    fps = np.array(current_snapshot_grp["fps"])
                    fp_primes_val = np.array(
                        current_snapshot_grp["fp_primes_val"]
                    )
                    fp_primes_row = np.array(
                        current_snapshot_grp["fp_primes_row"]
                    )
                    fp_primes_col = np.array(
                        current_snapshot_grp["fp_primes_col"]
                    )
                    fp_primes_size = np.array(
                        current_snapshot_grp["fp_primes_size"]
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
                        ref_positions,
                        calc_derivatives=calc_derivatives,
                        log=log,
                    )

                    if save_fps:
                        current_snapshot_grp.create_dataset(
                            "size_info", data=size_info
                        )
                        current_snapshot_grp.create_dataset("fps", data=fps)
                        current_snapshot_grp.create_dataset(
                            "fp_primes_val", data=fp_primes_val
                        )
                        current_snapshot_grp.create_dataset(
                            "fp_primes_row", data=fp_primes_row
                        )
                        current_snapshot_grp.create_dataset(
                            "fp_primes_col", data=fp_primes_col
                        )
                        current_snapshot_grp.create_dataset(
                            "fp_primes_size", data=fp_primes_size
                        )

                image_dict["descriptors"] = fps
                # image_dict["num_descriptors"] = size_info[2]
                fp_prime_dict["row"] = fp_primes_row
                fp_prime_dict["col"] = fp_primes_col
                fp_prime_dict["val"] = fp_primes_val
                image_dict["descriptor_primes"] = fp_prime_dict

            else:
                try:
                    # size_info = np.array(current_element_grp["size_info"])
                    fps = np.array(current_element_grp["fps"])
                except Exception:
                    size_info, fps, _, _, _, _ = self.calculate_fingerprints(
                        image,
                        ref_positions,
                        calc_derivatives=calc_derivatives,
                        log=log,
                    )

                    if save_fps:
                        current_snapshot_grp.create_dataset(
                            "size_info", data=size_info
                        )
                        current_snapshot_grp.create_dataset("fps", data=fps)

                image_dict["descriptors"] = fps

            descriptor_list.append(image_dict)

        return descriptor_list

    def _compute_fingerprints_nodb(
        self, image, ref_positions, calc_derivatives, save_fps, cores, log
    ):
        descriptor_list = []

        image_dict = {}

        if calc_derivatives:

            (
                size_info,
                fps,
                fp_primes_val,
                fp_primes_row,
                fp_primes_col,
                fp_primes_size,
            ) = self.calculate_fingerprints(
                image,
                ref_positions,
                calc_derivatives=calc_derivatives,
                log=log,
            )

            image_dict["descriptors"] = fps
            # image_dict["num_descriptors"] = size_info[2]
            fp_prime_dict["row"] = fp_primes_row
            fp_prime_dict["col"] = fp_primes_col
            fp_prime_dict["val"] = fp_primes_val
            image_dict["descriptor_primes"] = fp_prime_dict

        else:
            
            size_info, fps, _, _, _, _ = self.calculate_fingerprints(
                image,
                ref_positions,
                calc_derivatives=calc_derivatives,
                log=log,
            )

            image_dict["descriptors"] = fps

        descriptor_list.append(image_dict)
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

    def _setup_fingerprint_database(self, save_fps):
        self.get_descriptor_setup_hash()
        self.desc_type_database_dir = "{}/{}".format(
            self.fp_database, self.descriptor_type
        )

        self.desc_fp_database_dir = "{}/{}".format(
            self.desc_type_database_dir, self.descriptor_setup_hash
        )

        if save_fps:
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
