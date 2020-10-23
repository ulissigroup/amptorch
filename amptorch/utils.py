import skorch
from skorch.utils import to_numpy


def target_extractor(y):
    # TODO better remove need for numpy call here before GPU support
    if len(y) == 2:
        return (to_numpy(y[0]), to_numpy(y[1]))
    elif len(y) == 1:
        return (to_numpy(y[0]), None)

def to_tensor(X, device, accept_sparse=False):
    return X

class train_end_load_best_loss(skorch.callbacks.base.Callback):
    def __init__(self, filename):
        self.filename = filename

    def on_train_end(self, net, X, y):
        net.load_params("./checkpoints/{}/params.pt".format(self.filename))
