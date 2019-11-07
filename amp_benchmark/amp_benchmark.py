from amp import Amp
from amp.descriptor.gaussian import Gaussian, make_symmetry_functions
from amp.model.neuralnetwork import NeuralNetwork
from amp.model import LossFunction
from amp import analysis
import ase
from ase import io
from ase.md import Langevin
from ase import units
import sys
import numpy as np

Gs = {}
Gs["G2_etas"] = np.logspace(np.log10(0.05), np.log10(5.0), num=4)
Gs["G2_rs_s"] = [0] * 4
Gs["G4_etas"] = [0.005]
Gs["G4_zetas"] = [1.0]
Gs["G4_gammas"] = np.array([+1.0, -1.0])
Gs["cutoff"] = 5.876798323827276

elements = ['C', 'O', 'Cu']

G = make_symmetry_functions(elements=elements, type="G2", etas=Gs["G2_etas"])
G += make_symmetry_functions(
    elements=elements,
    type="G4",
    etas=Gs["G4_etas"],
    zetas=Gs["G4_zetas"],
    gammas=Gs["G4_gammas"],
)

# calc = Amp(descriptor=Gaussian(Gs=G, cutoff=Gs["cutoff"]),
        # model=NeuralNetwork((2, 2), mode='atom-centered'), label="calc")
# calc.model.lossfunction = LossFunction(
    # convergence={"energy_rmse": 0.02, "force_rmse": 0.1}
# )
images = ase.io.read("../datasets/COCu/COCu_pbc_300K.traj", ":")
IMAGES = []
for i in range(100):
    IMAGES.append(images[i])
calc = Amp.load(file='calc.amp')
# calc.load('calc.amp')

"""Generates test or training data with a simple MD simulation."""
traj = ase.io.Trajectory('test.traj', "w")
slab = images[0]
slab.set_calculator(calc)
traj.write(slab)
dyn = Langevin(slab, 1.0 * units.fs, 300 * units.kB, 0.002)
for step in range(100):
    dyn.run(20)
    traj.write(slab)
