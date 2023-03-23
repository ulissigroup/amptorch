import torch
import numpy as np

try:
    shell = get_ipython().__class__.__name__
    if shell == "ZMQInteractiveShell":
        from tqdm.notebook import tqdm
    else:
        from tqdm import tqdm
except NameError:
    from tqdm import tqdm


class FeatureScaler:
    """
    Normalizes an input tensor and later reverts it.
    Adapted from https://github.com/Open-Catalyst-Project/baselines
    """

    def __init__(
        self,
        data_list,
        forcetraining,
        scaling,
    ):
        self.transform = scaling["type"]
        if self.transform not in ["normalize", "standardize"]:
            raise NotImplementedError(f"{self.transform} scaling not supported.")
        if self.transform == "normalize" and "range" not in scaling:
            raise NotImplementedError("Normalization requires desire range.")
        if self.transform == "normalize":
            feature_range = scaling["range"]
        self.forcetraining = forcetraining
        self.elementwise = scaling.get("elementwise", True)
        self.threshold = scaling.get("threshold", 1e-6)
        fingerprints = torch.cat([data.fingerprint for data in data_list], dim=0)
        atomic_numbers = torch.cat([data.atomic_numbers for data in data_list], dim=0)

        if self.elementwise:
            self.unique = torch.unique(atomic_numbers).tolist()
            self.scales = {}
            for element in self.unique:
                idx = torch.where(atomic_numbers == element)[0]
                element_fps = fingerprints[idx]
                if self.transform == "standardize":
                    mean = torch.mean(element_fps, dim=0)
                    std = torch.std(element_fps, dim=0, unbiased=False)
                    std[std < self.threshold] = 1
                    self.scales[element] = {"offset": mean, "scale": std}
                else:
                    fpmin = torch.min(element_fps, dim=0).values
                    fpmax = torch.max(element_fps, dim=0).values
                    data_range = fpmax - fpmin
                    data_range[data_range < self.threshold] = 1
                    scale = (feature_range[1] - feature_range[0]) / (data_range)
                    offset = feature_range[0] - fpmin * scale
                    self.scales[element] = {"offset": offset, "scale": scale}

        else:
            if self.transform == "standardize":
                mean = torch.mean(fingerprints, dim=0)
                std = torch.std(fingerprints, dim=0, unbiased=False)
                std[std < self.threshold] = 1
                self.scale = {"offset": mean, "scale": std}
            else:
                fpmin = torch.min(fingerprints, dim=0).values
                fpmax = torch.max(fingerprints, dim=0).values
                data_range = fpmax - fpmin
                data_range[data_range < self.threshold] = 1
                scale = (feature_range[1] - feature_range[0]) / (data_range)
                offset = feature_range[0] - fpmin * scale
                self.scale = {"offset": offset, "scale": scale}

    def __eq__(self, other):
        """Overrides the default implementation"""
        if isinstance(other, FeatureScaler):
            if (
                self.transform != other.transform
                or self.elementwise != other.elementwise
            ):
                return False
            if self.elementwise:
                for element in self.scales:
                    if element not in other.scales:
                        return False
                    for key in self.scales[element]:
                        if key not in other.scales[element] or not torch.equal(
                            self.scales[element][key], other.scales[element][key]
                        ):
                            return False
            else:
                for key in self.scale:
                    if key not in other.scale or not torch.equal(
                        self.scale[key], other.scale[key]
                    ):
                        return False
            return True
        return NotImplemented

    def norm(self, data_list, disable_tqdm=False):
        if self.elementwise:
            for data in tqdm(
                data_list,
                desc="Scaling Feature data (%s)" % self.transform,
                total=len(data_list),
                unit=" scalings",
                disable=disable_tqdm,
            ):
                fingerprint = data.fingerprint
                atomic_numbers = data.atomic_numbers
                for element in self.unique:
                    element_idx = torch.where(atomic_numbers == element)
                    element_fp = fingerprint[element_idx]
                    if self.transform == "standardize":
                        element_fp = (
                            element_fp - self.scales[element]["offset"]
                        ) / self.scales[element]["scale"]
                    else:
                        element_fp = (
                            element_fp * self.scales[element]["scale"]
                        ) + self.scales[element]["offset"]
                    fingerprint[element_idx] = element_fp
                if self.forcetraining:
                    base_atoms = torch.repeat_interleave(
                        atomic_numbers, data.fingerprint.shape[1]
                    )
                    fp_idx = data.fprimes._indices()[0]
                    fp_idx_to_scale = fp_idx % data.fingerprint.shape[1]
                    element_idx = base_atoms[fp_idx].tolist()
                    _values = data.fprimes._values()

                    dict_elements = {element: [] for element in set(element_idx)}
                    for i, element in enumerate(element_idx):
                        dict_elements[element].append(i)

                    for element, ids in dict_elements.items():
                        scale = self.scales[element]["scale"][fp_idx_to_scale[ids]]
                        if self.transform == "standardize":
                            _values[ids] /= scale
                        else:
                            _values[ids] *= scale

                    _indices = data.fprimes._indices()
                    _size = data.fprimes.size()
                    data.fprimes = torch.sparse.FloatTensor(_indices, _values, _size)

        else:
            for data in tqdm(
                data_list,
                desc="Scaling Feature data (%s)" % self.transform,
                total=len(data_list),
                unit=" scalings",
                disable=disable_tqdm,
            ):
                fingerprint = data.fingerprint
                if self.transform == "standardize":
                    fingerprint = (fingerprint - self.scale["offset"]) / self.scale[
                        "scale"
                    ]
                else:
                    fingerprint = (fingerprint * self.scale["scale"]) + self.scale[
                        "offset"
                    ]
                data.fingerprint = fingerprint

                if self.forcetraining:
                    fp_idx = data.fprimes._indices()[0]
                    fp_idx_to_scale = fp_idx % data.fingerprint.shape[1]
                    _values = data.fprimes._values()

                    scale = self.scale["scale"][fp_idx_to_scale]
                    if self.transform == "standardize":
                        _values /= scale
                    else:
                        _values *= scale
                    _indices = data.fprimes._indices()
                    _size = data.fprimes.size()
                    data.fprimes = torch.sparse.FloatTensor(_indices, _values, _size)

        return data_list


class TargetScaler:
    """
    Normalizes an input tensor and later reverts it (standardize).
    Adapted from https://github.com/Open-Catalyst-Project/baselines
    """

    def __init__(self, data_list, forcetraining):
        self.forcetraining = forcetraining

        energies = torch.tensor([data.energy for data in data_list])
        self.target_mean = torch.mean(energies, dim=0)
        self.target_std = torch.std(energies, dim=0)

        if torch.isnan(self.target_std) or self.target_std == 0:
            self.target_mean = 0
            self.target_std = 1

    def __eq__(self, other):
        """Overrides the default implementation"""
        if isinstance(other, TargetScaler):
            return (
                self.target_mean == other.target_mean
                and self.target_std == other.target_std
            )
        return NotImplemented

    def norm(self, data_list, disable_tqdm=False):
        for data in tqdm(
            data_list,
            desc="Scaling Target data",
            total=len(data_list),
            unit=" scalings",
            disable=disable_tqdm,
        ):
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


class AtomicCorrectionScaler:
    """
    Fit an linear model based on energy and atomic compositions.
    """

    def __init__(self, data_list, load_correction_dictionary=None):
        atom_list = []
        for data in data_list:
            atom_list.extend(data.atomic_numbers.numpy())

        atom_list = np.unique(atom_list)

        atom_dict = {atom: i for i, atom in enumerate(atom_list)}
        atom_dict_rev = {i: atom for i, atom in enumerate(atom_list)}

        print("start preparing linear system")
        num_atom_mat = np.zeros((len(data_list), len(atom_list)))
        energy_vec = np.zeros((len(data_list),))
        for i, data in enumerate(data_list):
            energy_vec[i] = data.energy
            atom_numbers_list = data.atomic_numbers.numpy()
            for atom in atom_numbers_list:
                num_atom_mat[i, atom_dict[atom]] += 1

        # print(num_atom_mat)
        # print(energy_vec)
        lr_result = np.linalg.lstsq(num_atom_mat, energy_vec, rcond=None)
        linear_regression_result = lr_result[0].flatten()
        self.correction_dict = {
            atom_dict_rev[i]: correction
            for i, correction in enumerate(linear_regression_result)
        }
        # print(self.correction_dict)

    def __eq__(self, other):
        """Overrides the default implementation"""
        if isinstance(other, AtomicCorrectionScaler):
            return self.correction_dict == other.correction_dict
        return NotImplemented

    def norm(self, data_list, disable_tqdm=False):
        for data in tqdm(
            data_list,
            desc="Scaling Target data by atomic corrections",
            total=len(data_list),
            unit=" scalings",
            disable=disable_tqdm,
        ):
            atom_numbers_list = data.atomic_numbers.numpy()
            for atom in atom_numbers_list:
                data.energy -= self.correction_dict[atom]
        return data_list

    def denorm(self, tensor, data):
        atom_numbers_list = data.atomic_numbers.numpy()
        for atom in atom_numbers_list:
            tensor += self.correction_dict[atom]
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
