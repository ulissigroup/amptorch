import torch
from torch import optim
from skorch import NeuralNetRegressor
from skorch.dataset import CVSplit
from skorch.callbacks import Checkpoint, EpochScoring
from amptorch.core import AMPTorch
from amptorch import AMP
from amptorch.gaussian import SNN_Gaussian
from amptorch.skorch_model import AMP as AMP_skorch
from amptorch.skorch_model.utils import target_extractor, energy_score, forces_score
from amptorch.model import FullNN, CustomLoss
from amptorch.data_preprocess import AtomsDataset, collate_amp
import numpy as np
from ase import Atoms
from ase.calculators.emt import EMT
from ase.io import read

def test_skorch_val():
    distances = np.linspace(2, 5, 100)
    label = "example"
    images = []
    energies = []
    forces = []
    for l in distances:
        image = Atoms(
            "CuCCu",
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

    Gs = {}
    Gs["G2_etas"] = np.logspace(np.log10(0.05), np.log10(5.0), num=2)
    Gs["G2_rs_s"] = [0] * 2
    Gs["G4_etas"] = [0.005]
    Gs["G4_zetas"] = [1.0]
    Gs["G4_gammas"] = [+1.0, -1]
    Gs["cutoff"] = 6.5

    forcetraining = True
    training_data = AtomsDataset(
        images,
        SNN_Gaussian,
        Gs,
        forcetraining=forcetraining,
        label=label,
        cores=1,
        delta_data=None,
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
        batch_size=10,
        max_epochs=20,
        iterator_train__collate_fn=collate_amp,
        iterator_train__shuffle=True,
        iterator_valid__collate_fn=collate_amp,
        device=device,
        train_split=CVSplit(0.1, random_state=1),
        callbacks=[
            EpochScoring(
                forces_score,
                on_train=False,
                use_caching=True,
                target_extractor=target_extractor,
            ),
            EpochScoring(
                energy_score,
                on_train=False,
                use_caching=True,
                target_extractor=target_extractor,
            ),
        ],
    )
    val_indices = [80, 84, 33, 81, 93, 17, 36, 82, 69, 65]
    val_images = [images[idx] for idx in val_indices]
    val_energies = energies[val_indices]
    val_forces = np.concatenate(np.array([forces[idx] for idx in val_indices]))
    calc = AMP_skorch(training_data, net, "test")
    calc.train(overwrite=True)
    num_of_atoms = 3

    last_energy_score = net.history[-1]["energy_score"]
    last_forces_score = net.history[-1]["forces_score"]

    calculated_energies = np.array(
        [calc.get_potential_energy(image) for image in val_images]
    )
    energy_rmse = np.sqrt(
        (((calculated_energies - val_energies) / num_of_atoms) ** 2).sum()
        / len(val_images)
    )
    assert round(energy_rmse, 5) == round(
        last_energy_score, 5
    ), "Validation energy errors incorrect!"
    calculated_forces = np.concatenate(
        np.array([calc.get_forces(image) for image in val_images])
    )
    force_rmse = np.sqrt(
        (((calculated_forces - val_forces)) ** 2).sum()
        / (3 * num_of_atoms * len(val_images))
    )
    assert round(force_rmse, 5) == round(
        last_forces_score, 5
    ), "Validation force errors incorrect!"

def test_energy_only_skorch_val():
    distances = np.linspace(2, 5, 100)
    label = "example"
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

    Gs = {}
    Gs["G2_etas"] = np.logspace(np.log10(0.05), np.log10(5.0), num=2)
    Gs["G2_rs_s"] = [0] * 2
    Gs["G4_etas"] = [0.005]
    Gs["G4_zetas"] = [1.0]
    Gs["G4_gammas"] = [+1.0, -1]
    Gs["cutoff"] = 6.5

    forcetraining = False
    training_data = AtomsDataset(
        images,
        SNN_Gaussian,
        Gs,
        forcetraining=forcetraining,
        label=label,
        cores=1,
        delta_data=None,
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
        criterion__force_coefficient=0,
        optimizer=torch.optim.LBFGS,
        optimizer__line_search_fn="strong_wolfe",
        lr=1e-2,
        batch_size=5,
        max_epochs=20,
        iterator_train__collate_fn=collate_amp,
        iterator_train__shuffle=True,
        iterator_valid__collate_fn=collate_amp,
        device=device,
        train_split=CVSplit(cv=0.1, random_state=1),
        callbacks=[
            EpochScoring(
                energy_score,
                on_train=False,
                use_caching=True,
                target_extractor=target_extractor,
            ),
        ],
    )
    val_indices = [80, 84, 33, 81, 93, 17, 36, 82, 69, 65]
    val_images = [images[idx] for idx in val_indices]
    val_energies = energies[val_indices]
    calc = AMP_skorch(training_data, net, "test")
    calc.train(overwrite=True)
    num_of_atoms = 3

    last_energy_score = net.history[-1]["energy_score"]

    calculated_energies = np.array(
        [calc.get_potential_energy(image) for image in val_images]
    )
    energy_rmse = np.sqrt(
        (((calculated_energies - val_energies) / num_of_atoms) ** 2).sum()
        / len(val_images)
    )
    assert round(energy_rmse, 5) == round(
        last_energy_score, 5
    ), "Validation energy errors incorrect!"

def test_val():
    distances = np.linspace(2, 5, 10)
    label = "example"
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
        energies.append(image.get_potential_energy(apply_constraint=False))
        forces.append(image.get_forces(apply_constraint=False))

    energies = np.array(energies)
    forces = np.concatenate(np.array(forces))
    Gs = {}
    Gs["G2_etas"] = np.logspace(np.log10(0.05), np.log10(5.0), num=2)
    Gs["G2_rs_s"] = [0] * 2
    Gs["G4_etas"] = [0.005]
    Gs["G4_zetas"] = [1.0]
    Gs["G4_gammas"] = [+1.0, -1]
    Gs["cutoff"] = 6.5

    torch.set_num_threads(1)
    calc = AMP(
        model=AMPTorch(
            images,
            descriptor=SNN_Gaussian,
            Gs=Gs,
            val_frac=0.2,
            force_coefficient=0.3,
            label=label,
            save_logs=True,
        )
    )
    calc.model.device = "cpu"
    calc.model.structure = [2, 2]
    calc.model.convergence = {
        "energy": 0.3,
        "force": 0.3,
        "early_stop": False,
        "epochs": 20,
    }
    calc.model.loader_params = {"batch_size": None, "shuffle": False, "num_workers": 0}
    calc.model.criterion = CustomLoss
    calc.model.optimizer = optim.LBFGS
    calc.model.lr = 1e-2
    calc.model.fine_tune = None

    calc.train(overwrite=True)
    loaders = calc.model.atoms_dataloader
    val_loader = loaders["val"]
    val_idx = val_loader.sampler.indices
    validation_images = [images[idx] for idx in val_idx]
    val_energies = np.array(
        [
            image.get_potential_energy(apply_constraint=False)
            for image in validation_images
        ]
    )
    val_forces = np.concatenate(
        np.array(
            [image.get_forces(apply_constraint=False) for image in validation_images]
        )
    )

    f_read = open("./results/logs/example.txt", "r")
    lines = f_read.readlines()
    energy_line = lines[-3]
    forces_line = lines[-2]
    reported_energy_loss = float(energy_line.split()[-1])
    reported_forces_loss = float(forces_line.split()[-1])
    f_read.close()

    num_of_atoms = 3
    calculated_energies = np.array(
        [calc.get_potential_energy(image) for image in validation_images]
    )
    energy_rmse = np.sqrt(
        (((calculated_energies - val_energies) / num_of_atoms) ** 2).sum()
        / len(validation_images)
    )
    assert round(energy_rmse, 4) == round(
        reported_energy_loss, 4
    ), "Energy errors are incorrect!"

    calculated_forces = np.concatenate(
        np.array([calc.get_forces(image) for image in validation_images])
    )
    force_rmse = np.sqrt(
        (((calculated_forces - val_forces)) ** 2).sum()
        / (3 * num_of_atoms * len(validation_images))
    )
    assert round(force_rmse, 4) == round(
        reported_forces_loss, 4
    ), "Force errors are incorrect!"

def test_energy_only_val():
    distances = np.linspace(2, 5, 10)
    label = "example"
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
        energies.append(image.get_potential_energy(apply_constraint=False))
        forces.append(image.get_forces(apply_constraint=False))

    energies = np.array(energies)
    forces = np.concatenate(np.array(forces))
    Gs = {}
    Gs["G2_etas"] = np.logspace(np.log10(0.05), np.log10(5.0), num=2)
    Gs["G2_rs_s"] = [0] * 2
    Gs["G4_etas"] = [0.005]
    Gs["G4_zetas"] = [1.0]
    Gs["G4_gammas"] = [+1.0, -1]
    Gs["cutoff"] = 6.5

    torch.set_num_threads(1)
    calc = AMP(
        model=AMPTorch(
            images,
            descriptor=SNN_Gaussian,
            Gs=Gs,
            val_frac=0.2,
            force_coefficient=0,
            label=label,
            save_logs=True,
        )
    )
    calc.model.device = "cpu"
    calc.model.structure = [2, 2]
    calc.model.convergence = {
        "energy": 0.1,
        "force": 0.1,
        "early_stop": False,
        "epochs": 20,
    }
    calc.model.loader_params = {"batch_size": None, "shuffle": False, "num_workers": 0}
    calc.model.criterion = CustomLoss
    calc.model.optimizer = optim.LBFGS
    calc.model.lr = 1e-2
    calc.model.fine_tune = None

    calc.train(overwrite=True)
    loaders = calc.model.atoms_dataloader
    val_loader = loaders["val"]
    val_idx = val_loader.sampler.indices
    validation_images = [images[idx] for idx in val_idx]
    val_energies = np.array(
        [
            image.get_potential_energy(apply_constraint=False)
            for image in validation_images
        ]
    )

    f_read = open("./results/logs/example.txt", "r")
    lines = f_read.readlines()
    energy_line = lines[-2]
    reported_energy_loss = float(energy_line.split()[-1])
    f_read.close()

    num_of_atoms = 3
    calculated_energies = np.array(
        [calc.get_potential_energy(image) for image in validation_images]
    )
    energy_rmse = np.sqrt(
        (((calculated_energies - val_energies) / num_of_atoms) ** 2).sum()
        / len(validation_images)
    )
    assert round(energy_rmse, 4) == round(
        reported_energy_loss, 4
    ), "Energy errors are incorrect!"
