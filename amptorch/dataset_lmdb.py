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
        env = self.connect_db(self.db_path)
        self.keys = [f"{j}".encode("ascii") for j in range(env.stat()["entries"])]
        self.feature_scaler = pickle.loads(
            env.begin().get("feature_scaler".encode("ascii"))
        )
        self.target_scaler = pickle.loads(
            env.begin().get("target_scaler".encode("ascii"))
        )
        self.length = pickle.loads(env.begin().get("length".encode("ascii")))
        self.elements = pickle.loads(env.begin().get("elements".encode("ascii")))
        self.descriptor = self.get_descriptor(
            pickle.loads(env.begin().get("descriptor_setup".encode("ascii")))
        )
        env.close()

    def __len__(self):
        return self.length

    def __getitem__(self, idx):
        env = self.connect_db(self.db_path)
        data = env.begin().get(self.keys[idx])
        data_object = pickle.loads(data)
        env.close()

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
            map_size=1099511627776 * 2,
        )
        return env
