import sys
import torch
from skorch import NeuralNetRegressor
from skorch.dataset import CVSplit
from skorch.callbacks import Checkpoint, EpochScoring
from skorch.callbacks.lr_scheduler import LRScheduler
import skorch.callbacks.base
from amp.descriptor.gaussian import Gaussian
from amptorch.gaussian import SNN_Gaussian
from amptorch.model import FullNN, CustomMSELoss
from amptorch.data_preprocess import AtomsDataset, factorize_data, collate_amp, TestDataset
from amptorch.skorch_model import AMP
from amptorch.skorch_model.utils import target_extractor, energy_score, forces_score
from amptorch.analysis import parity_plot
from torch.utils.data import DataLoader
from torch.nn import init
from skorch.utils import to_numpy
import numpy as np
from ase import Atoms
from ase.calculators.emt import EMT
from ase.io import read

class train_end_load_best_valid_loss(skorch.callbacks.base.Callback):
    def on_train_end(self, net, X, y):
        net.load_params('valid_best_params.pt')

LR_schedule = LRScheduler('CosineAnnealingLR', T_max=5)
# saves best validation loss
cp = Checkpoint(monitor='valid_loss_best', fn_prefix='valid_best_')
# loads best validation loss at the end of training
load_best_valid_loss = train_end_load_best_valid_loss()

distances = np.linspace(2, 5, 100)
label = "example"
images = []
for l in distances:
    image = Atoms(
        "CuCO",
        [
            (-l * np.sin(0.65), l * np.cos(0.65), 0),
            (0, 0, 0),
            (l * np.sin(0.65), l * np.cos(0.65), 0),
        ],
    )
    image.set_cell([10, 10, 10])
    image.wrap(pbc=True)
    image.set_calculator(EMT())
    images.append(image)

# define symmetry functions to be used

Gs = {}
Gs["G2_etas"] = np.logspace(np.log10(0.05), np.log10(5.0), num=4)
Gs["G2_rs_s"] = [0] * 4
Gs["G4_etas"] = [0.005]
Gs["G4_zetas"] = [1.0]
Gs["G4_gammas"] = [+1.0, -1]
Gs["cutoff"] = 6.5

forcetraining = True
training_data = AtomsDataset(images, SNN_Gaussian, Gs, forcetraining=forcetraining,
        label=label, cores=1, delta_data=None)
unique_atoms = training_data.elements
fp_length = training_data.fp_length
device = "cpu"

net = NeuralNetRegressor(
    module=FullNN(unique_atoms, [fp_length, 3, 10], device, forcetraining=forcetraining),
    criterion=CustomLoss,
    criterion__force_coefficient=0.3,
    optimizer=torch.optim.LBFGS,
    optimizer__line_search_fn="strong_wolfe",
    lr=1e-1,
    batch_size=len(images),
    max_epochs=50,
    iterator_train__collate_fn=collate_amp,
    iterator_train__shuffle=False,
    iterator_valid__collate_fn=collate_amp,
    device=device,
    train_split=0,
    callbacks=[
        EpochScoring(
            forces_score,
            on_train=True,
            use_caching=True,
            target_extractor=target_extractor,
        ),
        EpochScoring(
            energy_score,
            on_train=True,
            use_caching=True,
            target_extractor=target_extractor,
        ),
    ],
)
calc = AMP(training_data, net, 'test')
calc.train(overwrite=True)
parity_plot(calc, images, data="energy", label=label)
parity_plot(calc, images, data="forces", label=label)
