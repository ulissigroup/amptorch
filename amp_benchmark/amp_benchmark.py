from amp import Amp
from amp.descriptor.gaussian import Gaussian
from amp.model.neuralnetwork import NeuralNetwork

calc=Amp(descriptor=Gaussian(),model=NeuralNetwork(),label='calc')
calc.train(images='sample_training_data.traj')
