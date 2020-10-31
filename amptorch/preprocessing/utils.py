import torch


class FeatureScaler:
    """
    Normalizes an input tensor and later reverts it.
    Adapted from https://github.com/Open-Catalyst-Project/baselines
    """

    def __init__(self, data_list, forcetraining):
        self.forcetraining = forcetraining
        fingerprints = torch.cat([data.fingerprint for data in data_list], dim=0)
        atomic_numbers = torch.cat([data.atomic_numbers for data in data_list], dim=0)
        self.unique = torch.unique(atomic_numbers).tolist()
        self.scales = {}
        for element in self.unique:
            idx = torch.where(atomic_numbers == element)[0]
            element_fps = fingerprints[idx]
            mean = torch.mean(element_fps, dim=0)
            std = torch.std(element_fps, dim=0)
            std[std == 0] = 1
            self.scales[element] = {"offset": mean, "scale": std}

    def norm(self, data_list, threshold=1e-6):
        for data in data_list:
            fingerprint = data.fingerprint
            atomic_numbers = data.atomic_numbers
            for element in self.unique:
                element_idx = torch.where(atomic_numbers == element)
                element_fp = fingerprint[element_idx]
                element_fp = (
                    element_fp - self.scales[element]["offset"]
                ) / self.scales[element]["scale"]
                fingerprint[element_idx] = element_fp
            if self.forcetraining:
                base_atoms = torch.repeat_interleave(
                    atomic_numbers, data.fingerprint.shape[1]
                )
                fp_idx = data.fprimes._indices()[0]
                fp_idx_to_scale = fp_idx % data.fingerprint.shape[1]
                element_idx = base_atoms[fp_idx].tolist()
                _values = data.fprimes._values()
                for i, element in enumerate(element_idx):
                    _values[i] /= self.scales[element]["scale"][fp_idx_to_scale[i]]

                _indices = data.fprimes._indices()
                _size = data.fprimes.size()
                data.fprimes = torch.sparse.FloatTensor(_indices, _values, _size)
        return data_list


class TargetScaler:
    """
    Normalizes an input tensor and later reverts it.
    Adapted from https://github.com/Open-Catalyst-Project/baselines
    """

    def __init__(self, data_list, forcetraining):
        self.forcetraining = forcetraining

        if len(data_list) > 1:
            energies = torch.tensor([data.energy for data in data_list])

            self.target_mean = torch.mean(energies, dim=0)
            self.target_std = torch.std(energies, dim=0)
        else:
            self.target_mean = 0
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
