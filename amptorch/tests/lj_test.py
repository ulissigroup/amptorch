from ase import Atoms
from ase.calculators.lj import LennardJones as LJ
from ase.calculators.emt import EMT
from amptorch.lj_model import lj_optim
from amptorch.gaussian import SNN_Gaussian
from amptorch.data_preprocess import AtomsDataset, collate_amp
from amptorch.core import AMPTorch
from amptorch.model import CustomLoss, FullNN
import numpy as np
import torch
from torch import optim
from skorch import NeuralNetRegressor
from skorch.callbacks import Checkpoint, EpochScoring
from amptorch.skorch_model.utils import forces_score, target_extractor, energy_score


def test_lj():
    atoms = Atoms("CuCuCu", [(0, 0, 1), (0, 0, 1.5), (0, 0, 0.5)])
    parameters = {'epsilon': 1.0, 'sigma': 1.0, 'rc': 10}
    atoms.set_calculator(LJ(**parameters))
    actual_energy = atoms.get_potential_energy()
    actual_forces = atoms.get_forces()

    p0 = [1, 1]
    params_dict = {"Cu": []}

    element_energies = {'Cu': 0}

    atoms = [atoms]
    lj_model = lj_optim(atoms, p0, params_dict, 10, "test", element_energies)
    fitted_params = p0
    lj_energies, lj_forces, num_atoms = lj_model.lj_pred(
        atoms, fitted_params, params_dict
    )

    assert round(actual_energy, 1) == lj_energies[0], "LJ energies don't match!"
    assert actual_forces.all() == lj_forces.all(), "LJ forces don't match!"

    p0 = [1, 1]

    lj_model = lj_optim(atoms, p0, params_dict, 10, "test", element_energies)
    fitted_params = lj_model.fit()
    lj_energies, lj_forces, num_atoms = lj_model.lj_pred(
        atoms, fitted_params, params_dict
    )

    assert round(actual_energy, 1) == round(
        lj_energies[0], 1
    ), "LJ energies don't match!"
    assert actual_forces.all() == lj_forces.all(), "LJ forces don't match!"

def test_ml_lj():
    from amptorch import AMP
    import sys

    distances = np.linspace(2, 5, 10)
    label = "ml_lj"
    images = []
    energies = []
    forces = []
    parameters = {'asap_cutoff': False}
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
        image.set_calculator(EMT(**parameters))
        images.append(image)
        energies.append(image.get_potential_energy())
        forces.append(image.get_forces())
    energies = np.array(energies)
    forces = np.concatenate(np.array(forces))

    element_energies = {}
    for element in ["Cu", "C", "O"]:
        atoms = Atoms(element, cell=[20, 20, 20])
        atoms.set_calculator(EMT())
        element_energies[element] = atoms.get_potential_energy()

    Gs = {}
    Gs["G2_etas"] = np.logspace(np.log10(0.05), np.log10(5.0), num=4)
    Gs["G2_rs_s"] = [0] * 4
    Gs["G4_etas"] = [0.005]
    Gs["G4_zetas"] = [1.0]
    Gs["G4_gammas"] = [+1.0, -1]
    Gs["cutoff"] = 6.5

    p0 = [
        1.33905162,
        0.12290683,
        0.64021468,
        0.08010004,
        2.29284676,
        0.29639983,
        12
    ]
    params_dict = {"C": [], "O": [], "Cu": []}
    lj_model = lj_optim(images, p0, params_dict, Gs["cutoff"], label,
            element_energies)
    fitted_params = lj_model.fit()
    lj_energies, lj_forces, num_atoms = lj_model.lj_pred(
        images, fitted_params, params_dict
    )
    for idx, atoms in enumerate(images):
        model_lj_energy, model_lj_forces, _ = lj_model.lj_pred(
            [atoms], fitted_params, params_dict
        )
        true_lj_energy = lj_energies[idx]
        true_lj_forces = lj_forces[idx]
        assert (
            model_lj_energy == true_lj_energy
        ), "LJ energies calculated by \
        lj_model are inconsistent with predicted lj_pred!"
        assert (
            model_lj_forces.all() == true_lj_forces.all()
        ), "LJ forces calculated by \
        lj_model are inconsistent with lj_pred!"
    lj_model.parity(predicted_energies=lj_energies, predicted_forces=lj_forces)
    lj_data = [lj_energies, lj_forces, num_atoms, fitted_params, params_dict, lj_model]

    torch.set_num_threads(1)

    calc = AMP(
        model=AMPTorch(
            images,
            descriptor=SNN_Gaussian,
            Gs=Gs,
            cores=1,
            force_coefficient=0.3,
            label=label,
            save_logs=True,
            lj_data=lj_data,
        )
    )
    calc.model.device = "cpu"
    calc.model.structure = [2, 2]
    calc.model.val_frac = 0
    calc.model.convergence = {
        "energy": 0.2,
        "force": 0.2,
        "epochs": 1e10,
        "early_stop": False,
    }
    calc.model.loader_params = {"batch_size": None, "shuffle": False, "num_workers": 0}
    calc.model.criterion = CustomLoss
    calc.model.optimizer = optim.LBFGS
    calc.model.lr = 1e-2
    calc.model.fine_tune = None

    calc.train(overwrite=True)
    num_of_atoms = 3
    calculated_energies = np.array(
        [calc.get_potential_energy(image) for idx, image in enumerate(images)]
    )
    energy_rmse = np.sqrt(
        (((calculated_energies - energies) / num_of_atoms) ** 2).sum() / len(images)
    )
    assert (
        energy_rmse <= calc.model.convergence["energy"]
    ), "Energy training convergence not met!"

    calculated_forces = np.concatenate(
        np.array([calc.get_forces(image) for idx, image in enumerate(images)])
    )
    force_rmse = np.sqrt(
        (((calculated_forces - forces)) ** 2).sum() / (3 * num_of_atoms * len(images))
    )
    assert (
        force_rmse <= calc.model.convergence["force"]
    ), "Force training convergence not met!"
test_ml_lj()
def test_skorch_lj():
    from amptorch.skorch_model import AMP

    cp = Checkpoint(monitor="valid_loss_best", fn_prefix="valid_best_")

    distances = np.linspace(2, 5, 10)
    label = "skorch_example"
    images = []
    energies = []
    forces = []
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
        energies.append(image.get_potential_energy())
        forces.append(image.get_forces())

    energies = np.array(energies)
    forces = np.concatenate(np.array(forces))

    element_energies = {}
    for element in ["Cu", "C", "O"]:
        atoms = Atoms(element, cell=[20, 20, 20])
        atoms.set_calculator(EMT())
        element_energies[element] = atoms.get_potential_energy()

    Gs = {}
    Gs["G2_etas"] = np.logspace(np.log10(0.05), np.log10(5.0), num=2)
    Gs["G2_rs_s"] = [0] * 2
    Gs["G4_etas"] = [0.005]
    Gs["G4_zetas"] = [1.0]
    Gs["G4_gammas"] = [+1.0, -1]
    Gs["cutoff"] = 6.5

    p0 = [
        1.33905162,
        0.12290683,
        0.64021468,
        0.08010004,
        2.29284676,
        0.29639983,
    ]
    params_dict = {"C": [], "O": [], "Cu": []}
    lj_model = lj_optim(images, p0, params_dict, Gs["cutoff"], label,
            element_energies)
    fitted_params = lj_model.fit()
    lj_energies, lj_forces, num_atoms = lj_model.lj_pred(
        images, fitted_params, params_dict
    )
    lj_data = [lj_energies, lj_forces, num_atoms, fitted_params, params_dict, lj_model]

    forcetraining = True
    training_data = AtomsDataset(
        images,
        SNN_Gaussian,
        Gs,
        forcetraining=forcetraining,
        label=label,
        cores=1,
        lj_data=lj_data,
    )
    batch_size = len(training_data)
    unique_atoms = training_data.elements
    fp_length = training_data.fp_length
    device = "cpu"

    net = NeuralNetRegressor(
        module=FullNN(
            unique_atoms, [fp_length, 2, 2], device, forcetraining=forcetraining
        ),
        criterion=CustomLoss,
        criterion__force_coefficient=0.3,
        optimizer=torch.optim.LBFGS,
        optimizer__line_search_fn="strong_wolfe",
        lr=1e-2,
        batch_size=batch_size,
        max_epochs=100,
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
    calc = AMP(training_data, net, "test")
    calc.train(overwrite=True)
    num_of_atoms = 3
    calculated_energies = np.array(
        [calc.get_potential_energy(image) for idx, image in enumerate(images)]
    )
    energy_rmse = np.sqrt(
        (((calculated_energies - energies) / num_of_atoms) ** 2).sum() / len(images)
    )
    last_energy_score = net.history[-1]["energy_score"]
    assert round(energy_rmse, 4) == round(
        last_energy_score, 4
    ), "Energy errors incorrect!"
    last_forces_score = net.history[-1]["forces_score"]

    calculated_forces = np.concatenate(
        np.array([calc.get_forces(image) for image in images])
    )
    force_rmse = np.sqrt(
        (((calculated_forces - forces)) ** 2).sum() / (3 * num_of_atoms * len(images))
    )
    assert round(force_rmse, 4) == round(
        last_forces_score, 4
    ), "Force errors incorrect!"
