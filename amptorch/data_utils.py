import torch


class Transform():
    """Normalizes an input tensor and later reverts it.
    Adapted from https://github.com/Open-Catalyst-Project/baselines"""
    def __init__(self, tensor):
        if len(torch.unique(tensor)) > 1:
            self.mean = torch.mean(tensor, dim=0)
            self.std = torch.std(tensor, dim=0)
        else:
            self.mean = tensor[0]
            self.std = 1

    def norm(self, tensor, energy=True):
        return (tensor - self.mean) / self.std if energy else tensor/self.std

    def denorm(self, tensor, energy=True):
        return tensor * self.std + self.mean if energy else tensor * self.std
