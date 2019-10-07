import time
from ase.io import read, write
from amp.model.neuralnetwork import NeuralNetwork
from amp.descriptor.gaussian import Gaussian, make_symmetry_functions
from amp.model import LossFunction
from ase.build import molecule
from ase.calculators.singlepoint import SinglePointCalculator as sp
from amp_simple_nn.convert import make_amp_descriptors_simple_nn
import numpy as np

atoms = molecule('O2')
atoms.set_cell([10,10,10])

atoms.set_calculator(sp(atoms=atoms, energy = -1,
                        forces = np.array([[-1,-1,-1],[-1,-1,-1]])))

#Convergence parameters
energy_rmse = 0.0001
force_rmse = 0.01
energy_maxresid = None
force_maxresid = None
convergence = {'energy_rmse': energy_rmse,
               'force_rmse': force_rmse,
               'energy_maxresid': energy_maxresid,
               'force_maxresid': force_maxresid}

ncores = 1
hiddenlayers = (5,5)
elements = atoms.get_chemical_symbols()

images = [atoms]
g2_etas = [0.005]
g2_rs_s = [0] * 4
g4_etas = [0.005]
g4_zetas = [1., 4.]
g4_gammas = [1., -1.]
cutoff = 4
t=time.time()
make_amp_descriptors_simple_nn(images,g2_etas,g2_rs_s,g4_etas,
                               g4_zetas,g4_gammas,cutoff)
print(time.time()-t)

G = make_symmetry_functions(elements=elements, type='G2',
                             etas=g2_etas)

#Add Rs parameter (0.0 for default) to comply with my version of AMP
#for g in G:
#        g['Rs'] =  0.0


G += make_symmetry_functions(elements=elements, type='G4',
                             etas=g4_etas,
                             zetas=g4_zetas,
                             gammas=g4_gammas)
