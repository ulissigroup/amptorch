import torch


class FeatureScaler:
    """Normalizes an input tensor and later reverts it.
    Adapted from https://github.com/Open-Catalyst-Project/baselines"""

    def __init__(self, data_list, forcetraining):
        self.forcetraining = forcetraining

        if len(data_list) > 1:
            fingerprints = torch.cat([data.fingerprint for data in data_list], dim=0)

            self.feature_max = torch.max(fingerprints, dim=0).values
            self.feature_min = torch.min(fingerprints, dim=0).values

        else:
            self.feature_mean = torch.mean(fingerprints, dim=0)
            self.feature_std = 1

    def norm(self, data_list):
        for data in data_list:
            idx_to_scale = torch.where((self.feature_max - self.feature_min) > 10e-8)[0]
            data.fingerprint[:, idx_to_scale] = -1 + 2 * (
                (data.fingerprint[:, idx_to_scale] - self.feature_min[idx_to_scale])
                / (self.feature_max[idx_to_scale] - self.feature_min[idx_to_scale])
            )

            if self.forcetraining:
                idx_to_scale_prime = data.fprimes._indices()[0] % (
                    data.fingerprint.shape[1] - 1
                )
                nonzero_idx = torch.where(
                    self.feature_max[idx_to_scale_prime]
                    - self.feature_min[idx_to_scale_prime]
                )[0]

                # print("=====================================")
                # print(data.fprimes)
                
                # print(self.feature_max[idx_to_scale_prime][nonzero_idx])
                # print(self.feature_min[idx_to_scale_prime][nonzero_idx])
                data.fprimes._values()[nonzero_idx] *= 2 / (
                    self.feature_max[idx_to_scale_prime][nonzero_idx]
                    - self.feature_min[idx_to_scale_prime][nonzero_idx]
                )
                # print("*************************************")
                # print(data.fprimes)
                _values = data.fprimes._values()
                _indices = data.fprimes._indices()
                _size = data.fprimes.size()
                data.fprimes = torch.sparse.FloatTensor(_indices, _values, _size)

        return data_list


class TargetScaler:
    """Normalizes an input tensor and later reverts it.
    Adapted from https://github.com/Open-Catalyst-Project/baselines"""

    def __init__(self, data_list, forcetraining):
        self.forcetraining = forcetraining

        if len(data_list) > 1:
            energies = torch.tensor([data.energy for data in data_list])

            self.target_mean = torch.mean(energies, dim=0)
            self.target_std = torch.std(energies, dim=0)
        else:
            self.target_mean = data_list[0].energy
            self.target_std = 1

    def norm(self, data_list):
        for data in data_list:
            data.energy = (data.energy - self.target_mean) / self.target_std

            if self.forcetraining:
                data.forces /= self.target_std
        return data_list

    def denorm(self, tensor, pred="energy"):
        if pred == "energy":
            tensor = (tensor * self.target_std) + self.target_mean
        elif pred == "forces":
            tensor = tensor * self.target_std

        return tensor


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
    out = torch.sparse.DoubleTensor(_indices, _values, _shapes)

    return out
