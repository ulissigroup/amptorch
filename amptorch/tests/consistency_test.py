import numpy as np
import torch
import torch.nn as nn
from ase import Atoms
from ase.calculators.emt import EMT
from torch.utils.data import DataLoader

from amptorch.dataset import AtomsDataset, DataCollater
from amptorch.model import BPNN, MLP


### Construct test data
def test_energy_force_consistency():
    """Gaussian/Neural non-periodic standard.

    Checks that the answer matches that expected from previous Mathematica
    calculations.
    """

    #: Making the list of non-periodic images
    images = [
        Atoms(
            symbols="PdOPd2",
            pbc=[False, False, False],
            calculator=EMT(),
            cell=np.array([[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]),
            positions=np.array(
                [[0.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0], [1.0, 0.0, 0.0]]
            ),
        ),
        Atoms(
            symbols="PdOPd2",
            pbc=[False, False, False],
            calculator=EMT(),
            cell=np.array([[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]),
            positions=np.array(
                [[0.0, 1.0, 0.0], [1.0, 2.0, 1.0], [-1.0, 1.0, 2.0], [1.0, 3.0, 2.0]]
            ),
        ),
        Atoms(
            symbols="PdO",
            pbc=[False, False, False],
            calculator=EMT(),
            cell=np.array([[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]),
            positions=np.array([[2.0, 1.0, -1.0], [1.0, 2.0, 1.0]]),
        ),
        Atoms(
            symbols="Pd2O",
            pbc=[False, False, False],
            calculator=EMT(),
            cell=np.array([[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]),
            positions=np.array([[-2.0, -1.0, -1.0], [1.0, 2.0, 1.0], [3.0, 4.0, 4.0]]),
        ),
        Atoms(
            symbols="Cu",
            pbc=[False, False, False],
            calculator=EMT(),
            cell=np.array([[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]),
            positions=np.array([[0.0, 0.0, 0.0]]),
        ),
    ]

    # Parameters
    Gs = {
        "default": {
            "G2": {
                "etas": [0.2],
                "rs_s": [0],
            },
            "G4": {"etas": [0.4], "zetas": [1], "gammas": [1]},
            "cutoff": 6.5,
        },
    }

    # "true" energies and forces from a known energy-force consistent MLP
    true_energies = [
        3.156074047088623,
        3.153623580932617,
        1.3491456508636475,
        1.9491342306137085,
        0.6553481817245483,
    ]
    true_forces = [
        [-0.010804446414113045, -0.015696275979280472, -0.018172090873122215],
        [-0.007649174891412258, 0.03893468156456947, -0.011910086497664452],
        [-0.0057038357481360435, -0.007940057665109634, 0.04719368368387222],
        [0.024157455191016197, -0.015298349782824516, -0.017111506313085556],
        [0.005043445620685816, -0.016683951020240784, -0.03857369348406792],
        [0.0007316730916500092, 0.009172249585390091, 0.003559479955583811],
        [-0.03596476465463638, -0.022213507443666458, 0.02469613030552864],
        [0.030189644545316696, 0.029725207015872, 0.010318083688616753],
        [0.017241766676306725, -0.017241766676306725, -0.03448353335261345],
        [-0.017241766676306725, 0.017241766676306725, 0.03448353335261345],
        [-0.01766442134976387, -0.01766442134976387, -0.011776280589401722],
        [-0.005222293548285961, -0.005222293548285961, -0.022553792223334312],
        [0.022886715829372406, 0.022886715829372406, 0.03433007374405861],
        [-0.0, -0.0, -0.0],
    ]

    elements = ["O", "Pd", "Cu"]
    np.random.seed(1)
    O_weights_1 = np.random.rand(10, 2)
    O_weights_2 = np.random.rand(1, 3).reshape(-1, 1)
    np.random.seed(2)
    Pd_weights_1 = np.random.rand(10, 2)
    Pd_weights_2 = np.random.rand(1, 3).reshape(-1, 1)
    np.random.seed(3)
    Cu_weights_1 = np.random.rand(10, 2)
    Cu_weights_2 = np.random.rand(1, 3).reshape(-1, 1)

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

    descriptor_name = "gaussian"
    descriptor_params = Gs
    descriptor_cutoff_params = {"cutoff_func": "cosine"}
    descriptor_elements = elements

    descriptor_setup = (
        descriptor_name,
        descriptor_params,
        descriptor_cutoff_params,
        descriptor_elements,
    )

    dataset = AtomsDataset(
        images,
        descriptor_setup,
        forcetraining=True,
        save_fps=True,
        cores=1,
    )

    collate_fn = DataCollater(train=False, forcetraining=True)
    dataloader = DataLoader(
        dataset, batch_size=len(images), collate_fn=collate_fn, shuffle=False
    )
    batch = next(iter(dataloader))
    input_dim = dataset.input_dim
    elements = [8, 46, 29]
    model = BPNN(elements, input_dim, 2, 2)

    model.state_dict()["elementwise_models.0.model_net.0.weight"].copy_(
        torch_O_weights_1
    )
    model.state_dict()["elementwise_models.0.model_net.0.bias"].copy_(torch_O_bias_1)
    model.state_dict()["elementwise_models.0.model_net.2.weight"].copy_(
        torch_O_weights_2
    )
    model.state_dict()["elementwise_models.0.model_net.2.bias"].copy_(torch_O_bias_2)
    model.state_dict()["elementwise_models.1.model_net.0.weight"].copy_(
        torch_Pd_weights_1
    )
    model.state_dict()["elementwise_models.1.model_net.0.bias"].copy_(torch_Pd_bias_1)
    model.state_dict()["elementwise_models.1.model_net.2.weight"].copy_(
        torch_Pd_weights_2
    )
    model.state_dict()["elementwise_models.1.model_net.2.bias"].copy_(torch_Pd_bias_2)
    model.state_dict()["elementwise_models.2.model_net.0.weight"].copy_(
        torch_Cu_weights_1
    )
    model.state_dict()["elementwise_models.2.model_net.0.bias"].copy_(torch_Cu_bias_1)
    model.state_dict()["elementwise_models.2.model_net.2.weight"].copy_(
        torch_Cu_weights_2
    )
    model.state_dict()["elementwise_models.2.model_net.2.bias"].copy_(torch_Cu_bias_2)

    for name, layer in model.named_modules():
        if isinstance(layer, MLP):
            layer.model_net = nn.Sequential(layer.model_net, nn.Tanh())

    energy_pred, force_pred = model(batch)

    for idx, i in enumerate(true_energies):
        assert round(i, 4) == round(
            energy_pred.tolist()[idx], 4
        ), "The predicted energy of image %i is wrong!" % (idx + 1)
    print("Energy predictions are correct!")
    for idx, sample in enumerate(true_forces):
        for idx_d, value in enumerate(sample):
            predict = force_pred.tolist()[idx][idx_d]
            assert abs(value - predict) < 0.0001, (
                "The predicted force of image % i, direction % i is wrong! Values: %s vs %s"
                % (idx + 1, idx_d, value, force_pred.tolist()[idx][idx_d])
            )
    print("Force predictions are correct!")


if __name__ == "__main__":
    test_energy_force_consistency()
