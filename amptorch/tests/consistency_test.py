import numpy as np
import torch
import torch.nn as nn
from amptorch.dataset import AtomsDataset, DataCollater
from amptorch.model import BPNN, MLP
from ase import Atoms
from ase.calculators.emt import EMT
from torch.utils.data import DataLoader


## Construct test data
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
        -0.6858080625534058,
        -0.3676419258117676,
        -1.020991325378418,
        -1.163346290588379,
        -0.48889073729515076,
    ]
    true_forces = [
        [-0.06926025450229645, -0.15385732054710388, -0.1956554651260376],
        [-0.09186100959777832, 0.40611717104911804, -0.1028066873550415],
        [-0.06731348484754562, -0.06853780150413513, 0.5004026293754578],
        [0.2284347414970398, -0.18372201919555664, -0.20194044709205627],
        [-0.16070456802845, -0.39853358268737793, -0.6236522197723389],
        [0.4219015836715698, 0.1250697821378708, -0.1488376408815384],
        [-0.6119220852851868, -0.3624090552330017, 0.3620723485946655],
        [0.35072505474090576, 0.6358729004859924, 0.4104175567626953],
        [0.002276031533256173, -0.002276031533256173, -0.004552063066512346],
        [-0.002276031533256173, 0.002276031533256173, 0.004552063066512346],
        [-0.001869507716037333, -0.001869507716037333, -0.001246338477358222],
        [0.000417390518123284, 0.000417390518123284, -0.0009318375959992409],
        [0.0014521172270178795, 0.0014521172270178795, 0.002178176073357463],
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
        scaling={"type": "normalize", "range": (-1, 1)},
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
        ), "The predicted energy of image %i is wrong!, %f, %f" % (
            idx + 1,
            round(i, 4),
            round(energy_pred.tolist()[idx], 4),
        )
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
