from amp import Amp
from ase.io import read, write
from amp.model.neuralnetwork import NeuralNetwork
from amp.descriptor.gaussian import Gaussian, make_symmetry_functions
from amp.model import LossFunction
import os

ncores = 24
hiddenlayers = (10,10)
images = read('images.traj', index = ':')

#Convergence parameters
energy_rmse = 0.0001
force_rmse = 0.01
energy_maxresid = None
force_maxresid = None
convergence = {'energy_rmse': energy_rmse,
               'force_rmse': force_rmse, 
               'energy_maxresid': energy_maxresid, 
               'force_maxresid': force_maxresid}



elements = ['H', 'O', 'Pt']


G = make_symmetry_functions(elements=elements, type='G2',
                             etas=[0.005, 4.0, 20.0, 80.0])

#Add Rs parameter (0.0 for default) to comply with my version of AMP
for g in G:
        g['Rs'] =  0.0


G += make_symmetry_functions(elements=elements, type='G4',
                             etas=[0.005],
                             zetas=[1., 4.],
                             gammas=[+1., -1.])




calc = Amp(descriptor = Gaussian(Gs = G, cutoff = 4.15),
           cores = ncores,
           model = NeuralNetwork(hiddenlayers = hiddenlayers))

calc.model.lossfunction = LossFunction(convergence = convergence, force_coefficient = 0.001)

calc.train(images = images)




