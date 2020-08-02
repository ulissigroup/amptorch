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
