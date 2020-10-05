from torch.utils.data import Dataset
from torch_geometric.data import Batch

from amptorch.preprocessing import (AtomsToData, FeatureScaler, TargetScaler,
                                    sparse_block_diag)


class AtomsDataset(Dataset):
    def __init__(
        self,
        images,
        descriptor,
        forcetraining=True,
        save_fps=True,
        cores=1,
    ):
        self.images = images
        self.a2d = AtomsToData(
            descriptor=descriptor,
            r_energy=True,
            r_forces=True,
            save_fps=save_fps,
            fprimes=forcetraining,
            cores=cores,
        )

        self.data_list = self.process()

    def process(self):
        data_list = self.a2d.convert_all(self.images)

        self.feature_scaler = FeatureScaler(data_list)
        self.target_scaler = TargetScaler(data_list)
        self.feature_scaler.norm(data_list)
        self.target_scaler.norm(data_list)

        return data_list

    @property
    def input_dim(self):
        return self.data_list[0].fingerprint.shape[1]

    def __len__(self):
        return len(self.data_list)

    def __getitem__(self, index):
        return self.data_list[index]


def data_collater(data_list, train=True):
    mtxs = []
    for data in data_list:
        mtxs.append(data.fprimes)
        data.fprimes = None
    batch = Batch.from_data_list(data_list)
    for i, data in enumerate(data_list):
        data.fprimes = mtxs[i]
    block_matrix = sparse_block_diag(mtxs)
    batch.fprimes = block_matrix
    if train:
        return batch, (batch.energy, batch.forces)
    else:
        return batch
