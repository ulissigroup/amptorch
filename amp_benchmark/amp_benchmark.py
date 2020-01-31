from amp import Amp
from amp.descriptor.gaussian import Gaussian, make_symmetry_functions
from amp.model.neuralnetwork import NeuralNetwork
from amp.regression import Regressor
from amp.model import LossFunction
from amp import analysis
import ase
from ase import io
from ase.md import Langevin
from ase import units
import sys
import numpy as np
from ase import Atoms
from ase.calculators.emt import EMT
from md_work.md_utils import md_run

Gs = {}
Gs["G2_etas"] = np.logspace(np.log10(0.05), np.log10(5.0), num=4)
Gs["G2_rs_s"] = [0] * 4
Gs["G4_etas"] = [0.005]
Gs["G4_zetas"] = [1.0, 4.0]
Gs["G4_gammas"] = np.array([+1.0, -1.0])
Gs["cutoff"] = 5.876798323827276

elements = ["Cu", "C", "O"]

G = make_symmetry_functions(elements=elements, type="G2", etas=Gs["G2_etas"])
G += make_symmetry_functions(
    elements=elements,
    type="G4",
    etas=Gs["G4_etas"],
    zetas=Gs["G4_zetas"],
    gammas=Gs["G4_gammas"],
)
images = ase.io.read("../datasets/COCu_nve_300K.traj", ":2000:20")

calc = Amp(
    descriptor=Gaussian(Gs=G, cutoff=Gs["cutoff"]),
    cores=6,
    model=NeuralNetwork((3, 60), mode="atom-centered"),
    label="calc_3x60",
)
calc.model.lossfunction = LossFunction(
    convergence={"energy_rmse": 0.02, "force_rmse": 0.1}
)
calc.model.regressor = Regressor(optimizer="L-BFGS-B")
# calc.train(images)
calc = Amp.load('./calc_3x60-untrained-parameters.amp')
md_run(calc, images[0], temp=300, dt=0.1, count=10, label="amp_lang_3x60", ensemble="NVE")
