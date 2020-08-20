from abc import ABC, abstractmethod
import numpy as np
import h5py
# import torch
import os
import time
from .util import get_traj_hash
 
class AMPTorchDescriptorBase(ABC):
 
    def __init__(self):
        super().__init__()
        self.fp_database = "./_saved_fingerprints_/"

        # To Be specified/calculated
        self.descriptor_type = "Default"
        self.descriptor_setup_hash = "Default"

        self.elements = []


    
    @abstractmethod
    def calculate_fingerprints(self, image, params_set, calculate_derivatives=True):
        # image is a single snapshot
        pass

    @abstractmethod
    def get_descriptor_setup_hash(self):
        #set self.descriptor_setup_hash
        pass

    @abstractmethod
    def save_descriptor_setup(self, filename):
        pass

    @abstractmethod
    def prepare_descriptor_parameters(self):
        # prepare self.params_set
        pass
    
    def prepare_fingerprints(self, trajs, parallel=False, log=None, calculate_derivatives=True, save=True, get_training_data = False):
        # params_set = self.params_set

        Total_Num_Trajs = len(trajs)

        trajs_descriptor_list = []


        self._setup_fingerprint_database(create_nonexist = save) # if save is true, create directories if not exist

        for traj in trajs:
            traj_start_time = time.time()
            traj_hash = get_traj_hash(traj)
            traj_db_filename = "{}/AmpFP-{}-{}-{}.h5".format(self.desc_fp_database_dir, self.descriptor_type, self.descriptor_setup_hash,traj_hash)

            # if save, then read/write from db as needed
            if save:
                temp_descriptor_list = \
                    self._prepare_fingerprints_single_traj(traj, traj_db_filename, parallel=parallel, log=log, calculate_derivatives=calculate_derivatives, save=True, get_training_data = get_training_data)
            
            # if not save, but db exist, read from db if possible
            elif os.path.exists(self.desc_type_database_dir):
                temp_descriptor_list = \
                    self._prepare_fingerprints_single_traj(traj, traj_db_filename, parallel=parallel, log=log, calculate_derivatives=calculate_derivatives, save=False, get_training_data = get_training_data)

            # if not save, and db not exist, calculat fps on the fly
            else:
                temp_descriptor_list = \
                    self._prepare_fingerprints_single_traj_nodb(traj, traj_db_filename, parallel=parallel, log=log, calculate_derivatives=calculate_derivatives, get_training_data = get_training_data)
            
            trajs_descriptor_list += temp_descriptor_list

            print("finished traj, took {}".format(time.time() - traj_start_time))

        return trajs_descriptor_list


    def _prepare_fingerprints_single_traj(self, traj, traj_db_filename, parallel=False, log=None, calculate_derivatives=True, save=True, get_training_data = False):
        descriptor_list = []
        # fingerprint_prime_list = []

        with h5py.File(traj_db_filename,'a') as db:

            Total_Num_Snapshots = len(traj)
            for i, snapshot in enumerate(traj):
                start_time = time.time()
                image_dict = {}
                if get_training_data:
                    image_dict["energy"] = snapshot.get_potential_energy()

                try:
                    current_snapshot_grp = db[str(i)]
                except:
                    current_snapshot_grp = db.create_group(str(i))
                
                for element in self.elements:
                    if element in snapshot.get_chemical_symbols():
                        # image_dict["descriptors"][element] = {}
                        # image_dict["descriptor_primes"][element] = {}
                        image_dict[element] = {}
                        if get_training_data:
                            image_dict[element]["forces"] = self._get_force_in_element_order(snapshot, element)
                        
                        try:
                            current_element_grp = current_snapshot_grp[element]
                        except:
                            current_element_grp = current_snapshot_grp.create_group(element)

                        if calculate_derivatives:
                            try:
                                size_info = np.array(current_element_grp["size_info"])
                                fps = np.array(current_element_grp["fps"])
                                fp_primes_val = np.array(current_element_grp["fp_primes_val"])
                                fp_primes_row = np.array(current_element_grp["fp_primes_row"])
                                fp_primes_col = np.array(current_element_grp["fp_primes_col"])
                                fp_primes_size = np.array(current_element_grp["fp_primes_size"])
                            except: 
                                size_info, fps, fp_primes_val, fp_primes_row, fp_primes_col, fp_primes_size = \
                                    self.calculate_fingerprints(snapshot, element, calculate_derivatives=calculate_derivatives)

                                if save:
                                    current_element_grp.create_dataset("size_info", data=size_info)
                                    current_element_grp.create_dataset("fps", data=fps)
                                    current_element_grp.create_dataset("fp_primes_val", data=fp_primes_val)
                                    current_element_grp.create_dataset("fp_primes_row", data=fp_primes_row)
                                    current_element_grp.create_dataset("fp_primes_col", data=fp_primes_col)
                                    current_element_grp.create_dataset("fp_primes_size", data=fp_primes_size)

                            # indices = np.vstack((fp_primes_row, fp_primes_col))
                            # torch_indices = torch.LongTensor(indices)
                            # torch_values = torch.FloatTensor(fp_primes_val)
                            # fp_prims_torch_sparse = torch.sparse.FloatTensor(torch_indices, torch_values, torch.Size(fp_primes_size))

                            # trajs_fingerprint_list.append(torch.from_numpy(fps))
                            # trajs_fingerprint_prime_list.append(trajs_fingerprint_list)

                            image_dict[element]["size_info"] = size_info

                            image_dict[element]["descriptors"] = fps

                            image_dict[element]["descriptor_primes"] = {}
                            image_dict[element]["descriptor_primes"]["value"] = fp_primes_val
                            image_dict[element]["descriptor_primes"]["row"]   = fp_primes_row
                            image_dict[element]["descriptor_primes"]["col"]   = fp_primes_col
                            image_dict[element]["descriptor_primes"]["size"]  = fp_primes_size
                        
                        else:
                            try:
                                size_info = np.array(current_element_grp["size_info"])
                                fps = np.array(current_element_grp["fps"])
                            except: 
                                size_info, fps, _, _, _, _ = \
                                    self.calculate_fingerprints(snapshot, element, calculate_derivatives=calculate_derivatives)

                                if save:
                                    current_element_grp.create_dataset("fps", data=fps)

                            # indices = np.vstack((fp_primes_row, fp_primes_col))

                            # trajs_fingerprint_list.append(torch.from_numpy(fps))
                            
                            image_dict[element]["size_info"] = size_info
                            image_dict[element]["descriptors"] = fps
                            
                    else:
                        print("element not in current snapshot: {}".format(element))

                descriptor_list.append(image_dict)
                
                took_time = time.time() - start_time
                print("finished snapshot {}/{}, took time: {}".format(i+1, Total_Num_Snapshots, took_time))
        
        return descriptor_list
    
    def _prepare_fingerprints_single_traj_nodb(self, traj, traj_db_filename, parallel=False, log=None, calculate_derivatives=True, get_training_data = False):
        descriptor_list = []

        Total_Num_Snapshots = len(traj)
        for i, snapshot in enumerate(traj):
            start_time = time.time()
            image_dict = {}
            if get_training_data:
                image_dict["energy"] = snapshot.get_potential_energy()
            
            for element in self.elements:
                if element in snapshot.get_chemical_symbols():
                    image_dict[element] = {}
                    if get_training_data:
                        image_dict[element]["forces"] = self._get_force_in_element_order(snapshot, element)

                    if calculate_derivatives:
                        size_info, fps, fp_primes_val, fp_primes_row, fp_primes_col, fp_primes_size = \
                            self.calculate_fingerprints(snapshot, element, calculate_derivatives=calculate_derivatives)

                        # indices = np.vstack((fp_primes_row, fp_primes_col))
                        # torch_indices = torch.LongTensor(indices)
                        # torch_values = torch.FloatTensor(fp_primes_val)
                        # fp_prims_torch_sparse = torch.sparse.FloatTensor(torch_indices, torch_values, torch.Size(fp_primes_size))

                        # fingerprint_list.append(torch.from_numpy(fps))
                        # fingerprint_prime_list.append(trajs_fingerprint_list)
                        image_dict[element]["size_info"] = size_info
                        image_dict[element]["descriptors"] = fps

                        image_dict[element]["descriptor_primes"] = {}
                        image_dict[element]["descriptor_primes"]["value"] = fp_primes_val
                        image_dict[element]["descriptor_primes"]["row"]   = fp_primes_row
                        image_dict[element]["descriptor_primes"]["col"]   = fp_primes_col
                        image_dict[element]["descriptor_primes"]["size"]  = fp_primes_size
                    
                    else:
                        fps, _, _, _, _ = \
                            self.calculate_fingerprints(snapshot, element, calculate_derivatives=calculate_derivatives)

                        # indices = np.vstack((fp_primes_row, fp_primes_col))

                        # fingerprint_list.append(torch.from_numpy(fps))
                        image_dict[element]["size_info"] = size_info
                        image_dict[element]["descriptors"] = fps
                        
                else:
                    print("element not in current snapshot: {}".format(element))
            
            descriptor_list.append(image_dict)
            took_time = time.time() - start_time
            print("finished snapshot {}/{}, took time: {}".format(i+1, Total_Num_Snapshots, took_time))

        return descriptor_list

    def _get_force_in_element_order(self, image, element):
        sym_arr = np.array(image.get_chemical_symbols())
        index = sym_arr == element
        # energy = image.get_potential_energy()
        forces = image.get_forces()
        return forces[index]

    def _setup_fingerprint_database(self, create_nonexist = True):
        self.get_descriptor_setup_hash()
        self.desc_type_database_dir = "{}/{}".format(self.fp_database, self.descriptor_type)

        self.desc_fp_database_dir = "{}/{}".format(self.desc_type_database_dir, self.descriptor_setup_hash)

        if create_nonexist:
            if not os.path.exists(self.fp_database):
                os.makedirs(self.fp_database)
            
            if not os.path.exists(self.desc_type_database_dir):
                os.makedirs(self.desc_type_database_dir)

            if not os.path.exists(self.desc_fp_database_dir):
                os.makedirs(self.desc_fp_database_dir)

            descriptor_setup_filename = "__descriptor_setup__.txt"
            descriptor_setup_path = "{}/{}".format(self.desc_fp_database_dir, descriptor_setup_filename)
            self.save_descriptor_setup(descriptor_setup_path)

    def _get_element_list(self):
        return self.elements
        
