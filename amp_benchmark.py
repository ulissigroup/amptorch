from amp import Amp
from amp.descriptor.gaussian import Gaussian
from amp.model.neuralnetwork import NeuralNetwork


#2 Hidden Layers (5,5), Activation=tanh, Epochs=100
model=NeuralNetwork()
print model.parameters
# calc=Amp(descriptor=Gaussian(),model=NeuralNetwork(),label='calc')
# calc.train(images='sample_training_data.traj')
