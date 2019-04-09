from amp import Amp
from amp.descriptor.gaussian import Gaussian
from amp.model.neuralnetwork import NeuralNetwork
from amp import analysis
import ase
from ase import io

calc = Amp(descriptor=Gaussian(), model=NeuralNetwork(), label="calc")
images = ase.io.read("defect-trajectory.extxyz", ":")
calc.train(images)
analysis.plot_parity("calc.amp", images)
