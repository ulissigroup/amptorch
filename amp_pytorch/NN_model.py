"""NN_model.py: Constructs a model consisting of element specific Neural
Networks as understood from Behler and Parrinello's works. A model instance is
constructed based off the unique number of atoms in the dataset."""

import sys
import numpy as np
import torch
import torch.nn as nn
from torch.nn import init
from torch.nn import Tanh, Softplus, LeakyReLU
from torch.nn.init import xavier_uniform_, kaiming_uniform_
from amp.utilities import Logger

__author__ = "Muhammed Shuaibi"
__email__ = "mshuaibi@andrew.cmu.edu"


class Dense(nn.Linear):
    """Constructs and applies a dense layer with an activation function (when
    available) y=activation(Ax+b)

    Arguments:
        input_size (int): number of input features
        output_size (int): number of output features
        bias (bool): True to include bias at each neuron. False otherwise
        (Default: True)
        activation (callable): activation function to be utilized
        (Default:None)
    """

    def __init__(self, input_size, output_size, bias=True, activation=None):
        self.activation = activation
        super(Dense, self).__init__(input_size, output_size, bias)

    def reset_parameters(self):
        """Weight initialization scheme"""

        # xavier_uniform_(self.weight,gain=5.0/3)
        kaiming_uniform_(self.weight, nonlinearity="tanh")
        if self.bias is not None:
            fan_in, _ = init._calculate_fan_in_and_fan_out(self.weight)
            bound = 1 / np.sqrt(fan_in)
            init.uniform_(self.bias, -bound, bound)

    def forward(self, inputs):
        neuron_output = super(Dense, self).forward(inputs)
        if self.activation:
            neuron_output = self.activation()(neuron_output)
        return neuron_output


class MLP(nn.Module):
    """Constructs a fully connected neural network model to be utilized for
    each element type'''

    Arguments:
        n_input_nodes: Number of input nodes (Default=20 using BP SF)
        n_output_nodes: Number of output nodes (Default=1)
        n_layers: Total number of layers in the neural network
        n_hidden_size: Number of neurons within each hidden layer
        activation: Activation function to be utilized. (Default=Tanh())
    """

    def __init__(
        self,
        n_input_nodes=20,
        n_output_nodes=1,
        n_layers=3,
        n_hidden_size=5,
        activation=Tanh,
    ):
        super(MLP, self).__init__()
        if isinstance(n_hidden_size, int):
            n_hidden_size = [n_hidden_size] * (n_layers - 1)
        self.n_neurons = [n_input_nodes] + n_hidden_size + [n_output_nodes]
        self.activation = activation
        layers = [
            Dense(self.n_neurons[i], self.n_neurons[i + 1], activation=activation)
            for i in range(n_layers - 1)
        ]
        layers.append(Dense(self.n_neurons[-2], self.n_neurons[-1], activation=None))
        self.model_net = nn.Sequential(*layers)

    def forward(self, inputs):
        """Feeds data forward in the neural network

        Arguments:
            inputs (torch.Tensor): NN inputs
        """

        return self.model_net(inputs)


class FullNN(nn.Module):
    """Combines element specific NNs into a model to predict energy of a given
    structure

    """

    def __init__(self, unique_atoms, device):
        log = Logger("results/results-log.txt")

        super(FullNN, self).__init__()
        self.device = device
        self.unique_atoms = unique_atoms
        n_unique_atoms = len(unique_atoms)
        elementwise_models = nn.ModuleList()
        for n_models in range(n_unique_atoms):
            elementwise_models.append(MLP())
        self.elementwise_models = elementwise_models
        log("Activation Function = %s" % elementwise_models[0].activation)

    def forward(self, input_data):
        batch_size = len(set(list(input_data.values())[0][1]))
        energy_pred = torch.zeros(batch_size, 1)
        energy_pred = energy_pred.to(self.device)
        for index, element in enumerate(self.unique_atoms):
            model_inputs = input_data[element][0]
            contribution_index = input_data[element][1]
            contributions_per_atom = int(len(contribution_index) / batch_size)
            contribution_index = torch.tensor(contribution_index)
            contribution_index = contribution_index.reshape(
                (batch_size, contributions_per_atom)
            )
            atomwise_outputs = self.elementwise_models[index].forward(model_inputs)
            atomwise_outputs = atomwise_outputs.reshape(
                batch_size, contributions_per_atom
            )
            element_pred = torch.sum(atomwise_outputs, 1).reshape(batch_size, 1)
            energy_pred += element_pred
        return energy_pred
