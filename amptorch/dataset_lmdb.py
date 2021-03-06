import lmdb
import pickle
import numpt as np
import bisect
from tqdm import tqdm
from torch.utils.data import Dataset
from amptorch.descriptor.Gaussian import Gaussian
from amptorch.descriptor.MCSH import AtomisticMCSH


class AtomsLMDBDataset(Dataset):
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
        descriptor_list = []
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
                temp_descriptor = self.get_descriptor(
                    pickle.loads(txn.get("descriptor_setup".encode("ascii")))
                )
                self.length_list.append(temp_length)
                feature_scaler_list.append(temp_feature_scaler)
                target_scaler_list.append(temp_target_scaler)
                descriptor_setup.append(temp_descriptor)

        self._keylen_cumulative = np.cumsum(self.length_list).tolist()
        self.total_length = np.sum(self.length_list)
        self.feature_scaler = feature_scaler_list[0]
        self.target_scaler = target_scaler_list[0]
        self.descriptor = descriptor_list[0]

        if len(self.db_paths) > 1:
            if any(
                feature_scalar != self.feature_scaler
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

    def __len__(self):
        return self.total_length

    def __getitem__(self, idx):
        db_idx = bisect.bisect(self._keylen_cumulative, idx)
        if db_idx != 0:
            el_idx = idx - self._keylen_cumulative[db_idx - 1]

        with self.envs[db_idx].begin(write=False) as txn:
            data = txn.get(self.keys_list[el_idx][el_idx])
            data_object = pickle.loads(data)

        return data_object

    def get_descriptor(self, descriptor_setup):
        fp_scheme, fp_params, cutoff_params, elements = descriptor_setup
        if fp_scheme == "gaussian":
            descriptor = Gaussian(Gs=fp_params, elements=elements, **cutoff_params)
        elif fp_scheme == "mcsh":
            descriptor = AtomisticMCSH(MCSHs=fp_params, elements=elements)
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


class AtomsLMDBDatasetCache(Dataset):
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
        descriptor_list = []
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
                temp_descriptor = self.get_descriptor(
                    pickle.loads(txn.get("descriptor_setup".encode("ascii")))
                )
                self.length_list.append(temp_length)
                feature_scaler_list.append(temp_feature_scaler)
                target_scaler_list.append(temp_target_scaler)
                descriptor_setup.append(temp_descriptor)

        self._keylen_cumulative = np.cumsum(self.length_list).tolist()
        self.total_length = np.sum(self.length_list)
        self.feature_scaler = feature_scaler_list[0]
        self.target_scaler = target_scaler_list[0]
        self.descriptor = descriptor_list[0]

        if len(self.db_paths) > 1:
            if any(
                feature_scalar != self.feature_scaler
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
        elif fp_scheme == "mcsh":
            descriptor = AtomisticMCSH(MCSHs=fp_params, elements=elements)
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
