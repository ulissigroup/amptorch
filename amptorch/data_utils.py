import torch


class Normalize:
    """Normalizes an input tensor and later reverts it.
    Adapted from https://github.com/Open-Catalyst-Project/baselines"""
    def __init__(self, tensor):
        self.mean = torch.mean(tensor, dim=0)
        self.std = torch.std(tensor, dim=0)

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
