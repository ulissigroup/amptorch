import torch


class Normalize:
    """Normalizes an input tensor and later reverts it.
    Adapted from https://github.com/Open-Catalyst-Project/baselines"""

    # TODO: clean normalization
    def __init__(self, data_list):
        if len(data_list) > 1:
            energies = torch.tensor([data.energy for data in data_list])
            fingerprints = torch.cat([data.fingerprint for data in data_list], dim=0)

            self.feature_mean = torch.mean(fingerprints, dim=0)
            self.feature_std = torch.std(fingerprints, dim=0)
            self.feature_std[self.feature_std == 0] = 10e-8

            self.target_mean = torch.mean(energies, dim=0)
            self.target_std = torch.std(energies, dim=0)
            self.target_std[self.target_std == 0] = 10e-8
        else:
            self.feature_mean = torch.mean(fingerprints, dim=0)
            self.feature_std = 1

            self.target_mean = data_list[0].energy
            self.target_std = 1

    def norm(self, data_list, energy=True):
        for data in data_list:
            data.energy = (data.energy - self.target_mean) / self.target_std
            data.forces /= self.target_std
            data.fingerprint = (data.fingerprint - self.feature_mean) / self.feature_std

            scaling_idx = data.fprimes.coalesce().indices()[0] % (
                data.fingerprint.shape[1] - 1
            )
            scalings = self.feature_std[scaling_idx]
            _values = data.fprimes._values() / scalings
            _indices = data.fprimes._indices()
            _size = data.fprimes.size()
            data.fprimes = torch.sparse.FloatTensor(_indices, _values, _size)

        return data_list

    def denorm(self, data_list, energy=True):
        for data in data_list:
            data.energy = (data.energy * self.target_std) + self.target_mean
            data.force = data.force * self.target_std

        return data_list


def sparse_block_diag(arrs):
    # TODO CUDA support
    r = []
    c = []
    v = []
    dim_1, dim_2 = 0, 0
    for k, mtx in enumerate(arrs):
        r += [mtx._indices()[0] + dim_1]
        c += [mtx._indices()[1] + dim_2]
        v += [mtx._values()]
        dim_1 += mtx.shape[0]
        dim_2 += mtx.shape[1]
    r = torch.cat(r, dim=0)
    c = torch.cat(c, dim=0)
    _indices = torch.stack([r, c])
    _values = torch.cat(v)
    _shapes = [dim_1, dim_2]
    out = torch.sparse.DoubleTensor(_indices, _values, _shapes,)

    return out
