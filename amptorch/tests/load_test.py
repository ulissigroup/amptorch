import sys
import torch
import copy
from skorch import NeuralNetRegressor
from skorch.dataset import CVSplit
from skorch.callbacks import Checkpoint, EpochScoring
from skorch.callbacks.lr_scheduler import LRScheduler
import skorch.callbacks.base
from amp.descriptor.gaussian import Gaussian
from amptorch.gaussian import SNN_Gaussian
from amptorch.model import BPNN, CustomMSELoss
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

def test_load():
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
        module=BPNN(unique_atoms, [fp_length, 3, 10], device, forcetraining=forcetraining),
        criterion=CustomMSELoss,
        criterion__force_coefficient=0.3,
        optimizer=torch.optim.LBFGS,
        optimizer__line_search_fn="strong_wolfe",
        lr=1e-1,
        batch_size=len(images),
        max_epochs=5,
        iterator_train__collate_fn=collate_amp,
        iterator_train__shuffle=False,
        iterator_valid__collate_fn=collate_amp,
        device=device,
        train_split=0,
        verbose=0,
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

    net_2 = copy.copy(net)

    calc_1 = AMP(training_data, net, 'test')
    calc_1.train(overwrite=True)

    energy_1 = calc_1.get_potential_energy(images[0])
    forces_1 = calc_1.get_forces(images[0])

    calc_2 = AMP(training_data, net_2, 'test')
    calc_2.load(filename='./results/trained_models/test.pt')
    energy_2 = calc_2.get_potential_energy(images[0])
    forces_2 = calc_2.get_forces(images[0])

    assert energy_1 == energy_2, "Energies do not match!"
    assert (forces_1 == forces_2).all(), "Forces do not match!"
