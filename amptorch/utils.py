import skorch
import json
import os
from skorch.utils import to_numpy
from torch_geometric.data import Batch
from torch.nn.parallel.scatter_gather import gather
import torch


class InOrderSplit:
    def __init__(self, val_frac):
        self.val_frac = val_frac

    def __call__(self, dataset):
        len_dataset = len(dataset)
        len_val = int(len_dataset * self.val_frac)
        len_train = len_dataset - len_val
        train_idx = list(range(0, len_train))
        val_idx = list(range(len_train, len_dataset))
        train_dataset = torch.utils.data.Subset(dataset, train_idx)
        val_dataset = torch.utils.data.Subset(dataset, val_idx)
        return train_dataset, val_dataset


def target_extractor(y):
    extracted = []
    for batch in y:
        energy_targets = to_numpy(batch[0])
        if len(batch) == 2:
            force_targets = to_numpy(batch[1])
            extracted.append([energy_targets, force_targets])
        elif len(batch) == 1:
            extracted.append([energy_targets, None])
    return extracted


def to_tensor(X, device, accept_sparse=False):
    if isinstance(X[0], Batch):
        return X
    else:
        for i, batch in enumerate(X):
            for j, targets in enumerate(batch):
                X[i][j] = targets.to(device)
        if device != "cpu":
            outputs = gather(X, device)
        else:
            outputs = X[0]
        return outputs


def save_normalizers(normalizers, path):
    tosave = {}
    tosave["feature"] = {
        "type": normalizers["feature"].transform,
        "scales": normalizers["feature"].scales.numpy(),
    }
    tosave["target"] = {
        "mean": normalizers["target"].target_mean.numpy(),
        "stddev": normalizers["target"].target_std.numpy(),
    }
    with open(path, "w", encoding="utf8") as json_file:
        json.dump(tosave, json_file, indent=4)
    return


class train_end_load_best_loss(skorch.callbacks.base.Callback):
    def __init__(self, filename):
        self.filename = filename

    def on_train_end(self, net, X, y):
        net.load_params("./checkpoints/{}/params.pt".format(self.filename))


class check_memory(skorch.callbacks.base.Callback):
    def on_batch_end(self, net, **kwargs):
        print(
            f"Allocated {torch.cuda.memory_allocated() / 1e6} Mb, cached {torch.cuda.memory_cached() / 1e6} Mb"
        )
