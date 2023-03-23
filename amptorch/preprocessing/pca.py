import torch
from sklearn.decomposition import PCA
from sklearn.preprocessing import normalize
import numpy as np

try:
    shell = get_ipython().__class__.__name__
    if shell == "ZMQInteractiveShell":
        from tqdm.notebook import tqdm
    else:
        from tqdm import tqdm
except NameError:
    from tqdm import tqdm


class PCAReducer:
    def __init__(
        self,
        data_list,
        forcetraining,
        pca_setting,
    ):
        self.forcetraining = forcetraining
        self.elementwise = pca_setting.get("elementwise", False)
        self.num_pc = pca_setting.get("num_pc", 20)
        self.normalize = pca_setting.get("normalize", False)

        fingerprints = torch.cat([data.fingerprint for data in data_list], dim=0)
        # atomic_numbers = torch.cat([data.atomic_numbers for data in data_list], dim=0)

        # fingerprints = normalize(fingerprints, norm='l2', axis=0)

        if self.elementwise:
            raise NotImplementedError

        else:
            if self.normalize:
                mean = torch.mean(fingerprints, dim=0)
                std = torch.std(fingerprints, dim=0, unbiased=False)
                std[std < 1e-8] = 1
                self.scale = {"offset": mean, "scale": std}
                fingerprints_normalized = (
                    fingerprints - self.scale["offset"]
                ) / self.scale["scale"]
                pca_reducer = PCA(n_components=self.num_pc).fit(
                    fingerprints_normalized.numpy()
                )
            else:
                pca_reducer = PCA(n_components=self.num_pc).fit(fingerprints.numpy())
            self.pca_components = torch.tensor(
                np.transpose(pca_reducer.components_), dtype=torch.get_default_dtype()
            )
            self.explained_variance = pca_reducer.explained_variance_
            self.explained_variance_ratio = pca_reducer.explained_variance_ratio_
            print(
                "explained variance ratio: {}\n total: {}".format(
                    self.explained_variance_ratio, np.sum(self.explained_variance_ratio)
                )
            )

    def reduce(self, data_list, disable_tqdm=False):
        if self.elementwise:
            raise NotImplementedError

        else:
            for data in tqdm(
                data_list,
                desc="PCA reducing to: {} components".format(self.num_pc),
                total=len(data_list),
                unit=" images",
                disable=disable_tqdm,
            ):
                fingerprint = data.fingerprint
                # print("size before: {}".format(fingerprint.size()))
                if self.normalize:
                    fingerprint = (fingerprint - self.scale["offset"]) / self.scale[
                        "scale"
                    ]
                fingerprint = torch.matmul(fingerprint, self.pca_components)
                # print("size after: {}".format(fingerprint.size()))
                data.fingerprint = fingerprint

                if self.forcetraining:
                    raise NotImplementedError
                    # fp_idx = data.fprimes._indices()[0]
                    # fp_idx_to_scale = fp_idx % data.fingerprint.shape[1]
                    # _values = data.fprimes._values()

                    # scale = self.scale["scale"][fp_idx_to_scale]
                    # if self.transform == "standardize":
                    #     _values /= scale
                    # else:
                    #     _values *= scale
                    # _indices = data.fprimes._indices()
                    # _size = data.fprimes.size()
                    # data.fprimes = torch.sparse.FloatTensor(_indices, _values, _size)

        return data_list
