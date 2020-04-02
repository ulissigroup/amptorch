"""
Force and energy calculation consistency test to be verified against AMP's
calculations
"""
import sys
import os
import numpy as np
from torch.nn import init, Tanh
import torch
from torch.utils.data import DataLoader
from skorch import NeuralNetRegressor
from ase import Atoms
from ase.calculators.emt import EMT
from collections import OrderedDict
from amp import Amp
from amp.descriptor.gaussian import Gaussian, make_symmetry_functions
from amp.model.neuralnetwork import NeuralNetwork
from amptorch.utils import hash_images
from amptorch.gaussian import SNN_Gaussian
from amp.model import calculate_fingerprints_range
from amptorch.data_preprocess import (
    AtomsDataset,
    factorize_data,
    collate_amp,
    TestDataset,
)
from amptorch.model import FullNN, CustomMSELoss, MLP
from amptorch.skorch_model import AMP
from amp.utilities import hash_images as amp_hash
from skorch.utils import to_tensor


def test_calcs():
    """Gaussian/Neural non-periodic standard.

    Checks that the answer matches that expected from previous Mathematica
    calculations.
    """

    #: Making the list of non-periodic images
    images = [
        Atoms(
            symbols="PdOPd2",
            pbc=np.array([False, False, False], dtype=bool),
            calculator=EMT(),
            cell=np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
            positions=np.array(
                [[0.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0], [1.0, 0.0, 0.0]]
            ),
        ),
        Atoms(
            symbols="PdOPd2",
            pbc=np.array([False, False, False], dtype=bool),
            calculator=EMT(),
            cell=np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
            positions=np.array(
                [[0.0, 1.0, 0.0], [1.0, 2.0, 1.0], [-1.0, 1.0, 2.0], [1.0, 3.0, 2.0]]
            ),
        ),
        Atoms(
            symbols="PdO",
            pbc=np.array([False, False, False], dtype=bool),
            calculator=EMT(),
            cell=np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
            positions=np.array([[2.0, 1.0, -1.0], [1.0, 2.0, 1.0]]),
        ),
        Atoms(
            symbols="Pd2O",
            pbc=np.array([False, False, False], dtype=bool),
            calculator=EMT(),
            cell=np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
            positions=np.array([[-2.0, -1.0, -1.0], [1.0, 2.0, 1.0], [3.0, 4.0, 4.0]]),
        ),
        Atoms(
            symbols="Cu",
            pbc=np.array([False, False, False], dtype=bool),
            calculator=EMT(),
            cell=np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
            positions=np.array([[0.0, 0.0, 0.0]]),
        ),
    ]

    # Parameters
    hiddenlayers = {"O": (2,), "Pd": (2,), "Cu": (2,)}

    Gs = {}
    Gs["G2_etas"] = [0.2]
    Gs["G2_rs_s"] = [0]
    Gs["G4_etas"] = [0.4]
    Gs["G4_zetas"] = [1]
    Gs["G4_gammas"] = [1]
    Gs["cutoff"] = 6.5

    elements = ["O", "Pd", "Cu"]

    G = make_symmetry_functions(elements=elements, type="G2", etas=Gs["G2_etas"])
    G += make_symmetry_functions(
        elements=elements,
        type="G4",
        etas=Gs["G4_etas"],
        zetas=Gs["G4_zetas"],
        gammas=Gs["G4_gammas"],
    )
    amp_images = amp_hash(images)
    descriptor = Gaussian(Gs=G, cutoff=Gs["cutoff"])
    descriptor.calculate_fingerprints(amp_images, calculate_derivatives=True)
    fingerprints_range = calculate_fingerprints_range(descriptor, amp_images)
    np.random.seed(1)
    O_weights_1 = np.random.rand(10, 2)
    O_weights_2 = np.random.rand(1, 3).reshape(-1, 1)
    np.random.seed(2)
    Pd_weights_1 = np.random.rand(10, 2)
    Pd_weights_2 = np.random.rand(1, 3).reshape(-1, 1)
    np.random.seed(3)
    Cu_weights_1 = np.random.rand(10, 2)
    Cu_weights_2 = np.random.rand(1, 3).reshape(-1, 1)

    weights = OrderedDict(
        [
            ("O", OrderedDict([(1, O_weights_1), (2, O_weights_2)])),
            ("Pd", OrderedDict([(1, Pd_weights_1), (2, Pd_weights_2)])),
            ("Cu", OrderedDict([(1, Cu_weights_1), (2, Cu_weights_2)])),
        ]
    )

    scalings = OrderedDict(
        [
            ("O", OrderedDict([("intercept", 0), ("slope", 1)])),
            ("Pd", OrderedDict([("intercept", 0), ("slope", 1)])),
            ("Cu", OrderedDict([("intercept", 0), ("slope", 1)])),
        ]
    )

    calc = Amp(
        descriptor,
        model=NeuralNetwork(
            hiddenlayers=hiddenlayers,
            weights=weights,
            scalings=scalings,
            activation="tanh",
            fprange=fingerprints_range,
            mode="atom-centered",
            fortran=False,
        ),
        logging=False,
    )

    amp_energies = [calc.get_potential_energy(image) for image in images]
    amp_forces = [calc.get_forces(image) for image in images]
    amp_forces = np.concatenate(amp_forces)

    torch_O_weights_1 = torch.FloatTensor(O_weights_1[:-1, :]).t()
    torch_O_bias_1 = torch.FloatTensor(O_weights_1[-1, :])
    torch_O_weights_2 = torch.FloatTensor(O_weights_2[:-1, :]).t()
    torch_O_bias_2 = torch.FloatTensor(O_weights_2[-1, :])
    torch_Pd_weights_1 = torch.FloatTensor(Pd_weights_1[:-1, :]).t()
    torch_Pd_bias_1 = torch.FloatTensor(Pd_weights_1[-1, :])
    torch_Pd_weights_2 = torch.FloatTensor(Pd_weights_2[:-1, :]).t()
    torch_Pd_bias_2 = torch.FloatTensor(Pd_weights_2[-1, :])
    torch_Cu_weights_1 = torch.FloatTensor(Cu_weights_1[:-1, :]).t()
    torch_Cu_bias_1 = torch.FloatTensor(Cu_weights_1[-1, :])
    torch_Cu_weights_2 = torch.FloatTensor(Cu_weights_2[:-1, :]).t()
    torch_Cu_bias_2 = torch.FloatTensor(Cu_weights_2[-1, :])

    device = "cpu"
    dataset = AtomsDataset(
        images,
        descriptor=Gaussian,
        cores=1,
        label="consistency",
        Gs=Gs,
        forcetraining=True,
    )

    fp_length = dataset.fp_length
    batch_size = len(dataset)
    dataloader = DataLoader(dataset, batch_size, collate_fn=collate_amp, shuffle=False)
    model = FullNN(elements, [fp_length, 2, 2], device, forcetraining=True)
    model.state_dict()["elementwise_models.O.model_net.0.weight"].copy_(
        torch_O_weights_1
    )
    model.state_dict()["elementwise_models.O.model_net.0.bias"].copy_(torch_O_bias_1)
    model.state_dict()["elementwise_models.O.model_net.2.weight"].copy_(
        torch_O_weights_2
    )
    model.state_dict()["elementwise_models.O.model_net.2.bias"].copy_(torch_O_bias_2)
    model.state_dict()["elementwise_models.Pd.model_net.0.weight"].copy_(
        torch_Pd_weights_1
    )
    model.state_dict()["elementwise_models.Pd.model_net.0.bias"].copy_(torch_Pd_bias_1)
    model.state_dict()["elementwise_models.Pd.model_net.2.weight"].copy_(
        torch_Pd_weights_2
    )
    model.state_dict()["elementwise_models.Pd.model_net.2.bias"].copy_(torch_Pd_bias_2)
    model.state_dict()["elementwise_models.Cu.model_net.0.weight"].copy_(
        torch_Cu_weights_1
    )
    model.state_dict()["elementwise_models.Cu.model_net.0.bias"].copy_(torch_Cu_bias_1)
    model.state_dict()["elementwise_models.Cu.model_net.2.weight"].copy_(
        torch_Cu_weights_2
    )
    model.state_dict()["elementwise_models.Cu.model_net.2.bias"].copy_(torch_Cu_bias_2)
    import torch.nn as nn
    for name, layer in model.named_modules():
        if isinstance(layer, MLP):
            layer.model_net = nn.Sequential(layer.model_net, Tanh())

    for batch in dataloader:
        x = to_tensor(batch[0], device)
        y = to_tensor(batch[1], device)
        energy_pred, force_pred = model(x)
    for idx, i in enumerate(amp_energies):
        assert round(i, 4) == round(
            energy_pred.tolist()[idx][0], 4
        ), "The predicted energy of image %i is wrong!" % (idx + 1)
    print("Energy predictions are correct!")
    for idx, sample in enumerate(amp_forces):
        for idx_d, value in enumerate(sample):
            predict = force_pred.tolist()[idx][idx_d]
            assert abs(value - predict) < 0.0001, (
                "The predicted force of image % i, direction % i is wrong! Values: %s vs %s"
                % (idx + 1, idx_d, value, force_pred.tolist()[idx][idx_d])
            )
    print("Force predictions are correct!")
