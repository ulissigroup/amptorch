import lmdb
import pickle
from torch.utils.data import Dataset
from amptorch.descriptor.Gaussian import Gaussian
from amptorch.descriptor.MCSH import AtomisticMCSH


class AtomsLMDBDataset(Dataset):
    def __init__(
        self,
        db_path,
    ):
        self.db_path = db_path
        self.env = self.connect_db(self.db_path)
        self.keys = [f"{j}".encode("ascii") for j in range(self.env.stat()["entries"])]
        with self.env.begin(write=False) as txn:
            self.feature_scaler = pickle.loads(
                txn.get("feature_scaler".encode("ascii"))
            )
            self.target_scaler = pickle.loads(txn.get("target_scaler".encode("ascii")))
            self.length = pickle.loads(txn.get("length".encode("ascii")))
            self.elements = pickle.loads(txn.get("elements".encode("ascii")))
            self.descriptor = self.get_descriptor(
                pickle.loads(txn.get("descriptor_setup".encode("ascii")))
            )

    def __len__(self):
        return self.length

    def __getitem__(self, idx):
        with self.env.begin(write=False) as txn:
            data = txn.get(self.keys[idx])
            data_object = pickle.loads(data)

            self.feature_scaler.norm([data_object])
            self.target_scaler.norm([data_object])

        return data_object

    def get_descriptor(self, descriptor_setup):
        fp_scheme, fp_params, elements, cutoff_params = descriptor_setup
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
