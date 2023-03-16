import lmdb
import pickle
import numpy as np
import bisect
import torch
from tqdm import tqdm
from torch.utils.data import Dataset
from amptorch.descriptor.Gaussian import Gaussian
from amptorch.descriptor.GMP import GMP
from amptorch.descriptor.GMPOrderNorm import GMPOrderNorm
from torch.utils.data.sampler import Sampler


class AtomsLMDBDataset(Dataset):
    """
    lmdb dataset with no cache

    This is the straight forward yet slow way of using lmdb files
    To access a given image of lmdb files (i.e. the __getitem__ function), it has to go to disk,
    connect to the corresponding lmdb file, and access the desired image.
    Since this is random access for each image, the performance is slow.

    It does support large amount of data, limited only by disk space, and NOT memory (RAM)

    It should be avoid for bad access performance if possible

    Parameter:
    db_paths [str] : a list of strings pointing to the paths of lmdb files.
    """

    def __init__(
        self,
        db_paths,
    ):
        if len(db_paths) < 1:
            raise ValueError("Please provide lmdb file paths")
        self.db_paths = db_paths
        self.envs = []
        self.keys_list = []
        self.length_list = []

        feature_scaler_list = []
        target_scaler_list = []
        descriptor_setup_list = []
        descriptor_list = []
        elements_list = []
        for db_path in self.db_paths:
            temp_env = self.connect_db(db_path)
            self.envs.append(temp_env)
            self.keys_list.append(
                [f"{j}".encode("ascii") for j in range(temp_env.stat()["entries"])]
            )
            with temp_env.begin(write=False) as txn:
                temp_feature_scaler = pickle.loads(
                    txn.get("feature_scaler".encode("ascii"))
                )
                temp_target_scaler = pickle.loads(
                    txn.get("target_scaler".encode("ascii"))
                )
                temp_length = pickle.loads(txn.get("length".encode("ascii")))
                temp_descriptor_setup = pickle.loads(
                    txn.get("descriptor_setup".encode("ascii"))
                )
                temp_descriptor = self.get_descriptor(temp_descriptor_setup)
                temp_elements = pickle.loads(txn.get("elements".encode("ascii")))
                self.length_list.append(temp_length)
                feature_scaler_list.append(temp_feature_scaler)
                target_scaler_list.append(temp_target_scaler)
                descriptor_setup_list.append(temp_descriptor_setup)
                descriptor_list.append(temp_descriptor)
                elements_list.append(temp_elements)

        self._keylen_cumulative = np.cumsum(self.length_list).tolist()
        self.total_length = np.sum(self.length_list)

        # use the scaler/setups from the first lmdb file, but check for consistency across all lmdb files
        self.feature_scaler = feature_scaler_list[0]
        self.target_scaler = target_scaler_list[0]
        self.descriptor_setup = descriptor_setup_list[0]
        self.descriptor = descriptor_list[0]
        self.elements = elements_list[0]
        self.loaded_db_idx = -1

        if len(self.db_paths) > 1:
            if any(
                feature_scaler != self.feature_scaler
                for feature_scaler in feature_scaler_list
            ):
                raise ValueError(
                    "Please make sure all lmdb used the same feature scaler"
                )
            if any(
                target_scaler != self.target_scaler
                for target_scaler in target_scaler_list
            ):
                raise ValueError(
                    "Please make sure all lmdb used the same target scaler"
                )
            if any(descriptor != self.descriptor for descriptor in descriptor_list):
                raise ValueError("Please make sure all lmdb used the same descriptor")
            if any(set(elements) != set(self.elements) for elements in elements_list):
                raise ValueError("Please make sure all lmdb used the same elements")

    def __len__(self):
        return self.total_length

    def __getitem__(self, idx):
        db_idx = bisect.bisect(self._keylen_cumulative, idx)
        if db_idx != 0:
            el_idx = idx - self._keylen_cumulative[db_idx - 1]
        else:
            el_idx = idx

        with self.envs[db_idx].begin(write=False) as txn:
            data = txn.get(self.keys_list[db_idx][el_idx])
            data_object = pickle.loads(data)

        return data_object

    def get_descriptor(self, descriptor_setup):
        fp_scheme, fp_params, cutoff_params, elements = descriptor_setup
        if fp_scheme == "gaussian":
            descriptor = Gaussian(Gs=fp_params, elements=elements, **cutoff_params)
        elif fp_scheme == "gmp":
            descriptor = GMP(MCSHs=fp_params, elements=elements)
        elif fp_scheme == "gmpordernorm":
            descriptor = GMPOrderNorm(MCSHs=fp_params, elements=elements)
        else:
            raise NotImplementedError
        return descriptor

    @property
    def input_dim(self):
        return self[0].fingerprint.shape[1]

    def connect_db(self, lmdb_path):
        env = lmdb.open(
            lmdb_path,
            subdir=False,
            readonly=True,
            lock=False,
            readahead=False,
            meminit=False,
            max_readers=1,
        )
        return env


class AtomsLMDBDatasetPartialCache(Dataset):
    """
    lmdb dataset with partial cache

    This is the optimized way for training with multiple lmdb files
    that CAN NOT be fitted into RAM all at once.
    It assumes the trainer to sequentially look at the images and lmdb files (i.e., first the
    images in lmdb_file1 in order, then images in lmdb_file2 in order, and so on)
    With the above assumption, the dataset load and cache all the images of the current accessing
    lmdb file into RAM, and access the desired image (__getitem__ ) from RAM.
    This is ~ 1-2 orders of magnitudes faster than no cache, because of serial acess of lmdb_files.

    It does support large amount of data, limited only by disk space,
    as long as each lmdb file can be loaded into RAM entirely.

    It has to be used with in-order spliter and randomized dataset.

    Parameter:
    db_paths [str] : a list of strings pointing to the paths of lmdb files.
    """

    def __init__(
        self,
        db_paths,
    ):
        if len(db_paths) < 1:
            raise ValueError("Please provide lmdb file paths")
        self.db_paths = db_paths
        self.envs = []
        self.keys_list = []
        self.length_list = []

        feature_scaler_list = []
        target_scaler_list = []
        descriptor_setup_list = []
        descriptor_list = []
        elements_list = []
        for db_path in self.db_paths:
            temp_env = self.connect_db(db_path)
            self.envs.append(temp_env)
            self.keys_list.append(
                [f"{j}".encode("ascii") for j in range(temp_env.stat()["entries"])]
            )
            with temp_env.begin(write=False) as txn:
                temp_feature_scaler = pickle.loads(
                    txn.get("feature_scaler".encode("ascii"))
                )
                temp_target_scaler = pickle.loads(
                    txn.get("target_scaler".encode("ascii"))
                )
                temp_length = pickle.loads(txn.get("length".encode("ascii")))
                temp_descriptor_setup = pickle.loads(
                    txn.get("descriptor_setup".encode("ascii"))
                )
                temp_descriptor = self.get_descriptor(temp_descriptor_setup)
                temp_elements = pickle.loads(txn.get("elements".encode("ascii")))
                self.length_list.append(temp_length)
                feature_scaler_list.append(temp_feature_scaler)
                target_scaler_list.append(temp_target_scaler)
                descriptor_setup_list.append(temp_descriptor_setup)
                descriptor_list.append(temp_descriptor)
                elements_list.append(temp_elements)

        self._keylen_cumulative = np.cumsum(self.length_list).tolist()
        self.total_length = np.sum(self.length_list)

        # use the scaler/setups from the first lmdb file, but check for consistency across all lmdb files
        self.feature_scaler = feature_scaler_list[0]
        self.target_scaler = target_scaler_list[0]
        self.descriptor_setup = descriptor_setup_list[0]
        self.descriptor = descriptor_list[0]
        self.elements = elements_list[0]
        self.loaded_db_idx = -1

        if len(self.db_paths) > 1:
            if any(
                feature_scaler != self.feature_scaler
                for feature_scaler in feature_scaler_list
            ):
                raise ValueError(
                    "Please make sure all lmdb used the same feature scaler"
                )
            if any(
                target_scaler != self.target_scaler
                for target_scaler in target_scaler_list
            ):
                raise ValueError(
                    "Please make sure all lmdb used the same target scaler"
                )
            if any(descriptor != self.descriptor for descriptor in descriptor_list):
                raise ValueError("Please make sure all lmdb used the same descriptor")
            if any(set(elements) != set(self.elements) for elements in elements_list):
                raise ValueError("Please make sure all lmdb used the same elements")

    def __len__(self):
        return self.total_length

    def __load_dataset__(self, db_idx):
        dataset = []
        with self.envs[db_idx].begin(write=False) as txn:
            for idx in range(self.length_list[db_idx]):
                data = txn.get(self.keys_list[db_idx][idx])
                data_object = pickle.loads(data)
                dataset.append(data_object)

        self.loaded_db_idx = db_idx
        self.loaded_dataset = dataset

    def __getitem__(self, idx):
        db_idx = bisect.bisect(self._keylen_cumulative, idx)
        if db_idx != 0:
            el_idx = idx - self._keylen_cumulative[db_idx - 1]
        else:
            el_idx = idx

        if db_idx != self.loaded_db_idx:
            self.__load_dataset__(db_idx)

        return self.loaded_dataset[el_idx]

    def get_descriptor(self, descriptor_setup):
        fp_scheme, fp_params, cutoff_params, elements = descriptor_setup
        if fp_scheme == "gaussian":
            descriptor = Gaussian(Gs=fp_params, elements=elements, **cutoff_params)
        elif fp_scheme == "gmp":
            descriptor = GMP(MCSHs=fp_params, elements=elements)
        elif fp_scheme == "gmpordernorm":
            descriptor = GMPOrderNorm(MCSHs=fp_params, elements=elements)
        else:
            raise NotImplementedError
        return descriptor

    def get_length_list(self):
        return self.length_list

    @property
    def input_dim(self):
        return self[0].fingerprint.shape[1]

    def connect_db(self, lmdb_path):
        env = lmdb.open(
            lmdb_path,
            subdir=False,
            readonly=True,
            lock=False,
            readahead=False,
            meminit=False,
            max_readers=1,
        )
        return env


class AtomsLMDBDatasetCache(Dataset):
    """
    lmdb dataset with full cache

    This is the fastest way for training with multiple lmdb files
    in case they CAN be fitted into RAM all at once.
    It loads all images into RAM from disk up front.

    It does not large amount of data, as it's limited by RAM size.

    It is the fastest way among the three for trianing, ~3x faster than partial caching.

    Parameter:
    db_paths [str] : a list of strings pointing to the paths of lmdb files.
    """

    def __init__(
        self,
        db_paths,
    ):
        if len(db_paths) < 1:
            raise ValueError("Please provide lmdb file paths")
        self.db_paths = db_paths
        self.envs = []
        self.keys_list = []
        self.length_list = []

        feature_scaler_list = []
        target_scaler_list = []
        descriptor_setup_list = []
        descriptor_list = []
        elements_list = []
        for db_path in self.db_paths:
            temp_env = self.connect_db(db_path)
            self.envs.append(temp_env)
            self.keys_list.append(
                [f"{j}".encode("ascii") for j in range(temp_env.stat()["entries"])]
            )
            with temp_env.begin(write=False) as txn:
                temp_feature_scaler = pickle.loads(
                    txn.get("feature_scaler".encode("ascii"))
                )
                temp_target_scaler = pickle.loads(
                    txn.get("target_scaler".encode("ascii"))
                )
                temp_length = pickle.loads(txn.get("length".encode("ascii")))
                temp_descriptor_setup = pickle.loads(
                    txn.get("descriptor_setup".encode("ascii"))
                )
                temp_descriptor = self.get_descriptor(temp_descriptor_setup)
                temp_elements = pickle.loads(txn.get("elements".encode("ascii")))
                self.length_list.append(temp_length)
                feature_scaler_list.append(temp_feature_scaler)
                target_scaler_list.append(temp_target_scaler)
                descriptor_setup_list.append(temp_descriptor_setup)
                descriptor_list.append(temp_descriptor)
                elements_list.append(temp_elements)

        self._keylen_cumulative = np.cumsum(self.length_list).tolist()
        self.total_length = np.sum(self.length_list)

        # use the scaler/setups from the first lmdb file, but check for consistency across all lmdb files
        self.feature_scaler = feature_scaler_list[0]
        self.target_scaler = target_scaler_list[0]
        self.descriptor_setup = descriptor_setup_list[0]
        self.descriptor = descriptor_list[0]
        self.elements = elements_list[0]

        if len(self.db_paths) > 1:
            if any(
                feature_scaler != self.feature_scaler
                for feature_scaler in feature_scaler_list
            ):
                raise ValueError(
                    "Please make sure all lmdb used the same feature scaler"
                )
            if any(
                target_scaler != self.target_scaler
                for target_scaler in target_scaler_list
            ):
                raise ValueError(
                    "Please make sure all lmdb used the same target scaler"
                )
            if any(descriptor != self.descriptor for descriptor in descriptor_list):
                raise ValueError("Please make sure all lmdb used the same descriptor")
            if any(set(elements) != set(self.elements) for elements in elements_list):
                raise ValueError("Please make sure all lmdb used the same elements")

        self.data_list = []
        for i, env in enumerate(self.envs):
            with self.envs[i].begin(write=False) as txn:
                for idx in tqdm(
                    range(self.length_list[i]),
                    desc="loading from {}".format(self.db_paths[i]),
                    total=self.length_list[i],
                    unit=" images",
                ):
                    data = txn.get(self.keys_list[i][idx])
                    data_object = pickle.loads(data)
                    self.data_list.append(data_object)

    def __len__(self):
        return self.total_length

    def __getitem__(self, idx):
        return self.data_list[idx]

    def get_descriptor(self, descriptor_setup):
        fp_scheme, fp_params, cutoff_params, elements = descriptor_setup
        if fp_scheme == "gaussian":
            descriptor = Gaussian(Gs=fp_params, elements=elements, **cutoff_params)
        elif fp_scheme == "gmp":
            descriptor = GMP(MCSHs=fp_params, elements=elements)
        elif fp_scheme == "gmpordernorm":
            descriptor = GMPOrderNorm(MCSHs=fp_params, elements=elements)
        else:
            raise NotImplementedError
        return descriptor

    @property
    def input_dim(self):
        return self[0].fingerprint.shape[1]

    def connect_db(self, lmdb_path):
        env = lmdb.open(
            lmdb_path,
            subdir=False,
            readonly=True,
            lock=False,
            readahead=False,
            meminit=False,
            max_readers=1,
        )
        return env


class PartialCacheSampler(Sampler):
    """
    Sampling strategy for partial cache scheme.
    """

    def __init__(self, length_list, val_frac):
        len_cumulative = np.cumsum(length_list)
        len_dataset = np.sum(length_list)
        len_val = int(len_dataset * val_frac)
        len_train = len_dataset - len_val
        for i, cum_len in enumerate(len_cumulative):
            if cum_len >= len_train:
                self.length_list = length_list[: i + 1]
                self.length_list[-1] -= cum_len - len_train
                break

        self.num_datasets = len(self.length_list)
        self.start_idx_list = [0] + np.cumsum(self.length_list).tolist()
        self.total_length = np.sum(self.length_list)

    def __iter__(self):
        datapoint_order = []
        dataset_order = torch.randperm(self.num_datasets).tolist()
        for dataset_idx in dataset_order:
            start_idx = self.start_idx_list[dataset_idx]
            datapoint_order += [
                i + start_idx
                for i in torch.randperm(self.length_list[dataset_idx]).tolist()
            ]
        return iter(datapoint_order)


def get_lmdb_dataset(lmdb_paths, cache_type):
    """
    A helper function to assign lmdb dataset types.
    """
    if cache_type == "full":
        return AtomsLMDBDatasetCache(lmdb_paths)
    elif cache_type == "partial":
        return AtomsLMDBDatasetPartialCache(lmdb_paths)
    elif cache_type == "no":
        return AtomsLMDBDataset(lmdb_paths)
    else:
        raise NotImplementedError
