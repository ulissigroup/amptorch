from amp import Amp
from amp.descriptor.gaussian import Gaussian
from amp.model.neuralnetwork import NeuralNetwork
import ase
from ase import io

calc=Amp(descriptor=Gaussian(),model=NeuralNetwork(),label='calc')
images=ase.io.read('water.extxyz',':')
calc.train(images)
