import itertools

import torch


class Normalize:
    """Normalizes an input tensor and later reverts it.
    Adapted from https://github.com/Open-Catalyst-Project/baselines"""

    def __init__(self, data_list):
        if len(data_list) > 1:
            energies = torch.tensor([data.energy for data in data_list])
            self.mean = torch.mean(energies, dim=0)
            self.std = torch.std(energies, dim=0)
        else:
            self.mean = data_list[0].energy
            self.std = 1

    def norm(self, data_list, energy=True):
        if energy:
            for data in data_list:
                data.energy = (data.energy - self.mean) / self.std
        else:
            for data in data_list:
                data.force = data.force / self.std

        return data_list

    def denorm(self, data_list, energy=True):
        if energy:
            for data in data_list:
                data.energy = (data.energy * self.std) + self.mean
        else:
            for data in data_list:
                data.force = data.force * self.std

        return data_list


# Adapted from https://github.com/pytorch/pytorch/issues/31942
def sparse_block_diag(arrs):
    bad_args = [
        k
        for k in range(len(arrs))
        if not (isinstance(arrs[k], torch.Tensor) and arrs[k].ndim == 2)
    ]
    if bad_args:
        raise ValueError(
            "arguments in the following positions must be 2-dimension tensor: %s"
            % bad_args
        )

    shapes = torch.tensor([a.shape for a in arrs])

    i = []
    v = []
    r, c = 0, 0
    for k, (rr, cc) in enumerate(shapes):
        i += [
            torch.LongTensor(
                list(
                    itertools.product(torch.arange(r, r + rr), torch.arange(c, c + cc))
                )
            ).t()
        ]
        v += [arrs[k].to_dense().flatten()]
        r += rr
        c += cc
    out = torch.sparse.DoubleTensor(
        torch.cat(i, dim=1), torch.cat(v), torch.sum(shapes, dim=0).tolist()
    )
    return out
