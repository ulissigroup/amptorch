import skorch
from torch_geometric.data import Batch
from torch.nn.parallel.scatter_gather import gather


def target_extractor(y):
    return y


def to_tensor(X, device, accept_sparse=False):
    if isinstance(X[0], Batch):
        return X
    else:
        for i, targets in enumerate(X):
            X[i][0] = targets[0].to(device)
            X[i][1] = targets[1].to(device)
        outputs = gather(X, device)
        return (outputs[0], outputs[1])


class train_end_load_best_loss(skorch.callbacks.base.Callback):
    def __init__(self, filename):
        self.filename = filename

    def on_train_end(self, net, X, y):
        net.load_params("./checkpoints/{}/params.pt".format(self.filename))
