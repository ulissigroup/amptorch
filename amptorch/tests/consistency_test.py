"""
Force and energy calculation consistency test to be verified against AMP's
calculations
"""

import numpy as np
from torch.nn import init
from torch.utils.data import DataLoader
from ase import Atoms
from ase.calculators.emt import EMT
from collections import OrderedDict
from amp import Amp
from amp.descriptor.gaussian import Gaussian, make_symmetry_functions
from amp.model.neuralnetwork import NeuralNetwork
from amp.utilities import hash_images
from amp.model import calculate_fingerprints_range
from amptorch import core
from amptorch.data_preprocess import AtomsDataset, factorize_data, collate_amp
from amptorch.NN_model import FullNN, Dense
from amptorch import AMP
from amptorch.core import AMPTorch


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
    Gs["G4_zetas"] = [0.6]
    Gs["G4_gammas"] = [0.4]
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
    hashed_images = hash_images(images)
    descriptor = Gaussian(Gs=G, cutoff=Gs["cutoff"])
    descriptor.calculate_fingerprints(hashed_images, calculate_derivatives=True)
    fingerprints_range = calculate_fingerprints_range(descriptor, hashed_images)

    weights = OrderedDict(
        [
            (
                "O",
                OrderedDict(
                    [
                        (
                            1,
                            np.matrix(
                                [
                                    [0.5, 0.5],
                                    [0.5, 0.5],
                                    [0.5, 0.5],
                                    [0.5, 0.5],
                                    [0.5, 0.5],
                                    [0.5, 0.5],
                                    [0.5, 0.5],
                                    [0.5, 0.5],
                                    [0.5, 0.5],
                                    [0.5, 0.5],
                                ]
                            ),
                        ),
                        (2, np.matrix([[0.5], [0.5], [0.5]])),
                    ]
                ),
            ),
            (
                "Pd",
                OrderedDict(
                    [
                        (
                            1,
                            np.matrix(
                                [
                                    [0.5, 0.5],
                                    [0.5, 0.5],
                                    [0.5, 0.5],
                                    [0.5, 0.5],
                                    [0.5, 0.5],
                                    [0.5, 0.5],
                                    [0.5, 0.5],
                                    [0.5, 0.5],
                                    [0.5, 0.5],
                                    [0.5, 0.5],
                                ]
                            ),
                        ),
                        (2, np.matrix([[0.5], [0.5], [0.5]])),
                    ]
                ),
            ),
            (
                "Cu",
                OrderedDict(
                    [
                        (
                            1,
                            np.matrix(
                                [
                                    [0.5, 0.5],
                                    [0.5, 0.5],
                                    [0.5, 0.5],
                                    [0.5, 0.5],
                                    [0.5, 0.5],
                                    [0.5, 0.5],
                                    [0.5, 0.5],
                                    [0.5, 0.5],
                                    [0.5, 0.5],
                                    [0.5, 0.5],
                                ]
                            ),
                        ),
                        (2, np.matrix([[0.5], [0.5], [0.5]])),
                    ]
                ),
            ),
        ]
    )

    scalings = OrderedDict(
        [
            ("O", OrderedDict([("intercept", 0), ("slope", 1)])),
            ("Pd", OrderedDict([("intercept", 0), ("slope", 1)])),
            ("Cu", OrderedDict([("intercept", 0), ("slope", 1)])),
        ]
    )

    # Testing pure-python and fortran versions of Gaussian-neural force call
    device = "cpu"
    dataset = AtomsDataset(
        images, descriptor=Gaussian, Gs=Gs, cores=1, forcetraining=True
    )
    fp_length = dataset.fp_length
    batch_size = len(dataset)
    dataloader = DataLoader(dataset, batch_size, collate_fn=collate_amp, shuffle=False)
    model = FullNN(elements, [fp_length, 2, 2], device, forcetraining=True)
    for name, layer in model.named_modules():
        if isinstance(layer, Dense):
            layer.activation = None
            init.constant_(layer.weight, 0.5)
            init.constant_(layer.bias, 0.5)
    for batch in dataloader:
        input_data = [batch[0], len(batch[1])]
        for element in elements:
            input_data[0][element][0] = (
                input_data[0][element][0].to(device).requires_grad_(True)
            )
        fp_primes = batch[3]
        energy_pred, force_pred = model(input_data, fp_primes)

    calc = Amp(
        descriptor,
        model=NeuralNetwork(
            hiddenlayers=hiddenlayers,
            weights=weights,
            scalings=scalings,
            activation="linear",
            fprange=fingerprints_range,
            mode="atom-centered",
            fortran=False,
        ),
    )

    amp_energies = [calc.get_potential_energy(image) for image in images]
    amp_forces = [calc.get_forces(image) for image in images]
    amp_forces = np.concatenate(amp_forces)

    for idx, i in enumerate(amp_energies):
        assert round(i, 4) == round(
            energy_pred.tolist()[idx][0], 4
        ), "The predicted energy of image %i is wrong!" % (idx + 1)
    print("Energy predictions are correct!")
    for idx, sample in enumerate(amp_forces):
        for idx_d, value in enumerate(sample):
            predict = force_pred.tolist()[idx][idx_d]
            assert abs(value - predict) < 0.00001, (
                "The predicted force of image % i, direction % i is wrong! Values: %s vs %s"
                % (idx + 1, idx_d, value, force_pred.tolist()[idx][idx_d])
            )
    print("Force predictions are correct!")


if __name__ == "__main__":
    test_calcs()
