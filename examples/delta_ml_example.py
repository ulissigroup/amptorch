"Amptorch example of delta-ml model using a morse potential"

import numpy as np

import skorch.callbacks.base
from skorch import NeuralNetRegressor
from skorch.callbacks import Checkpoint, EpochScoring
from skorch.callbacks.lr_scheduler import LRScheduler

from amptorch.gaussian import SNN_Gaussian
from amptorch.model import FullNN, CustomMSELoss
from amptorch.data_preprocess import (
    AtomsDataset,
    collate_amp,
)
from amptorch.skorch_model import AMP
from amptorch.skorch_model.utils import target_extractor, energy_score, forces_score, train_end_load_best_loss
from amptorch.delta_models.morse import morse_potential

import torch
from torch.nn import init

from ase import Atoms
from ase.calculators.emt import EMT
from ase.io import read


label = "example"

LR_schedule = LRScheduler("CosineAnnealingLR", T_max=5)
# saves best validation loss
cp = Checkpoint(monitor="valid_loss_best", fn_prefix="valid_best_")
# loads best validation loss at the end of training
load_best_valid_loss = train_end_load_best_loss(label)

distances = np.linspace(2, 5, 10)
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


morse_params = {
    "C": {"re": 0.972, "D": 6.379, "sig": 0.477},
    "O": {"re": 1.09, "D": 8.575, "sig": 0.603},
    "Cu": {"re": 2.168, "D": 3.8386, "sig": 1.696},
}
morse_model = morse_potential(images, morse_params, Gs["cutoff"], label, combo="mean")
morse_energies, morse_forces, num_atoms = morse_model.morse_pred(images, morse_params)
morse_data = [morse_energies, morse_forces, num_atoms, morse_params, morse_model]

forcetraining = True
training_data = AtomsDataset(
    images,
    SNN_Gaussian,
    Gs,
    forcetraining=forcetraining,
    label=label,
    cores=1,
    delta_data=morse_data,
)
unique_atoms = training_data.elements
fp_length = training_data.fp_length
device = "cpu"

net = NeuralNetRegressor(
    module=FullNN(
        unique_atoms, [fp_length, 3, 10], device, forcetraining=forcetraining
    ),
    criterion=CustomMSELoss,
    criterion__force_coefficient=0.3,
    optimizer=torch.optim.LBFGS,
    optimizer__line_search_fn="strong_wolfe",
    lr=1e-2,
    batch_size=len(training_data),
    max_epochs=10,
    iterator_train__collate_fn=collate_amp,
    iterator_train__shuffle=True,
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
calc = AMP(training_data, net, "test")
calc.train(overwrite=True)
