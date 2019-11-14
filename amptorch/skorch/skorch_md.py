import ase
import time
import sys
import torch
from torch.nn import MSELoss
from skorch import NeuralNetRegressor
from skorch.dataset import CVSplit
from skorch.callbacks import Checkpoint, EpochScoring
from amptorch.gaussian import Gaussian
from amptorch.skorch.model_skorch import FullNN, CustomLoss
from amptorch.skorch.skorch_data import AtomsDataset, factorize_data, collate_amp, TestDataset
from amptorch.skorch import AMP
from amptorch.lj_model import lj_optim
from torch.utils.data import DataLoader
from torch.nn import init
from skorch.utils import to_numpy
import numpy as np
from sklearn.pipeline import Pipeline
from ase import Atoms, units
from ase.calculators.emt import EMT
from ase.md import Langevin


def target_extractor(y):
    return (
        (to_numpy(y[0]), to_numpy(y[1]))
        if len(y) == 2
        else (to_numpy(y[0]), to_numpy(y[1]), to_numpy(y[2]))
    )

def energy_score(net, X, y):
    mse_loss = MSELoss(reduction="sum")
    energy_pred = net.infer()[0]
    device = energy_pred.device
    energy_target = torch.tensor(y[0]).to(device)
    num_atoms = torch.tensor(y[1]).to(device)
    dataset_size = len(energy_pred)
    sd_scaling = scalings[0]
    mean_scaling = scalings[1]
    raw_preds = (energy_pred * sd_scaling) + mean_scaling
    raw_preds_per_atom = torch.div(raw_preds, num_atoms)
    raw_targets = (energy_target * sd_scaling) + mean_scaling
    target_per_atom = torch.div(raw_targets, num_atoms)
    energy_loss = mse_loss(raw_preds_per_atom, target_per_atom)
    energy_loss /= dataset_size
    energy_rmse = torch.sqrt(energy_loss)
    return energy_rmse

def forces_score(net, X, y):
    mse_loss = MSELoss(reduction="sum")
    sd_scaling = scalings[0]
    force_pred = net.infer()[1] * sd_scaling
    device = force_pred.device
    num_atoms = torch.tensor(y[1]).to(device)
    force_target = torch.tensor(y[-1], device=device)
    dataset_size = len(num_atoms)
    raw_force_target = force_target * sd_scaling
    num_atoms_force = torch.cat([idx.repeat(int(idx)) for idx in num_atoms])
    num_atoms_force = torch.sqrt(num_atoms_force).reshape(len(num_atoms_force), 1)
    force_pred_per_atom = torch.div(force_pred, num_atoms_force)
    force_targets_per_atom = torch.div(raw_force_target, num_atoms_force)
    force_mse = mse_loss(force_pred_per_atom, force_targets_per_atom)
    force_mse /= 3 * dataset_size
    force_rmse = torch.sqrt(force_mse)
    return force_rmse

# define symmetry functions to be used
Gs = {}
Gs["G2_etas"] = np.logspace(np.log10(0.05), np.log10(5.0), num=4)
Gs["G2_rs_s"] = [0] * 4
Gs["G4_etas"] = [0.005]
Gs["G4_zetas"] = [1.0, 4.0]
Gs["G4_gammas"] = [+1.0, -1]
Gs["cutoff"] = 6.0

label = "example"
images = ase.io.read("../../datasets/COCu/COCu_pbc_300K.traj", ":100")
lj_data = None
cutoff = Gs["cutoff"]
eV_kcalmol = 0.043372093
p0 = [
    1.33905162,
    0.12290683,
    6.41914719,
    0.64021468,
    0.08010004,
    8.26082762,
    2.29284676,
    0.29639983,
    0.08071821,
]
params_dict = {"C": [], "O": [], "Cu": []}
lj_model = lj_optim(label, images, p0, params_dict, cutoff)
fitted_params = lj_model.fit()
# fitted_params = p0
lj_energies, lj_forces, num_atoms = lj_model.lj_pred(
    images, fitted_params, params_dict
)
lj_data = [
    lj_energies,
    lj_forces,
    num_atoms,
    fitted_params,
    params_dict,
    lj_model,
]
 
forcetraining = True
training_data = AtomsDataset(images, Gaussian, Gs, forcetraining=forcetraining,
        label=label, cores=1, lj_data=lj_data)
scalings = training_data.scalings
unique_atoms = training_data.elements
fp_length = training_data.fp_length
device = "cpu"

net = NeuralNetRegressor(
    module=FullNN(unique_atoms, [fp_length, 30, 30], device, forcetraining=forcetraining),
    criterion=CustomLoss,
    criterion__force_coefficient=0.3,
    optimizer=torch.optim.Adam,
    lr=1e-2,
    batch_size=400,
    max_epochs=1000,
    iterator_train__collate_fn=collate_amp,
    iterator_valid__collate_fn=collate_amp,
    device=device,
    train_split=CVSplit(0.1),
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
calc = AMP(training_data, net, 'test2')
calc.train(overwrite=True)

traj = ase.io.Trajectory('test2.traj', "w")
slab = images[0].copy()
slab.set_calculator(calc)
slab.get_forces()
traj.write(slab)
dyn = Langevin(slab, 1.0 * units.fs, 300 * units.kB, 0.002)
time_start = time.time()
for step in range(100):
    dyn.run(20)
    traj.write(slab)
