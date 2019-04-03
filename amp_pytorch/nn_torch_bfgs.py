import sys
import numpy as np
import torch
import torch.nn as nn
from torch.nn import init
from torch.nn import Tanh, Softplus, LeakyReLU
from torch.nn.init import xavier_uniform_, kaiming_uniform_
import copy
import time
import matplotlib.pyplot as plt
from amp.utilities import Logger


class Dense(nn.Linear):
    """Constructs and applies a dense layer with an activation function (when
    available) y=activation(Ax+b)

    Arguments:
        input_size (int): number of input features
        output_size (int): number of output features
        bias (bool): True to include bias at each neuron. False otherwise
        (Default: True)
        activation (callable): activation function to be utilized (Default:None)
    """

    def __init__(self, input_size, output_size, bias=True, activation=None):
        self.activation = activation
        super(Dense, self).__init__(input_size, output_size, bias)

    def reset_parameters(self):
        """Weight initialization scheme"""

        # xavier_uniform_(self.weight,gain=5.0/3)
        kaiming_uniform_(self.weight, nonlinearity='tanh')
        if self.bias is not None:
            fan_in, _ = init._calculate_fan_in_and_fan_out(self.weight)
            bound = 1/np.sqrt(fan_in)
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

    def __init__(self, n_input_nodes=20, n_output_nodes=1, n_layers=3, n_hidden_size=5, activation=Tanh):
        super(MLP, self).__init__()
        # if n_hidden_size is None:
        # implement pyramid neuron structure across each layer
        if type(n_hidden_size) is int:
            n_hidden_size = [n_hidden_size] * (n_layers-1)
        self.n_neurons = [n_input_nodes]+n_hidden_size+[n_output_nodes]
        self.activation = activation
        layers = [Dense(self.n_neurons[i], self.n_neurons[i+1],
                        activation=activation) for i in range(n_layers-1)]
        layers.append(
            Dense(self.n_neurons[-2], self.n_neurons[-1], activation=None))
        self.model_net = nn.Sequential(*layers)
        # self.model_net=nn.Sequential(HiddenLayer1,HiddenLayer2,OutputLayer3)

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

    def __init__(self, unique_atoms, batch_size):
        log = Logger('../benchmark_results/results-log.txt')

        super(FullNN, self).__init__()
        # self.device = device
        self.unique_atoms = unique_atoms
        self.batch_size = batch_size
        n_unique_atoms = len(unique_atoms)
        elementwise_models = nn.ModuleList()
        for n_models in range(n_unique_atoms):
            elementwise_models.append(MLP())
        self.elementwise_models = elementwise_models
        log('Activation Function = %s' % elementwise_models[0].activation)

    def forward(self, input_data):
        energy_pred = torch.zeros(self.batch_size, 1)
        # energy_pred = energy_pred.to(self.device)
        for index, element in enumerate(self.unique_atoms):
            model_inputs = input_data[element][0]
            contribution_index = input_data[element][1]
            contributions_per_atom = len(contribution_index)/self.batch_size
            contribution_index = torch.tensor(contribution_index)
            contribution_index = contribution_index.reshape(
                self.batch_size, contributions_per_atom)
            atomwise_outputs = self.elementwise_models[index].forward(
                model_inputs)
            atomwise_outputs = atomwise_outputs.reshape(self.batch_size,
                                                        contributions_per_atom)
            element_pred = torch.sum(
                atomwise_outputs, 1).reshape(self.batch_size, 1)
            energy_pred += element_pred
            # for cindex,atom_output in enumerate(atomwise_outputs):
            # energy_pred[contribution_index[cindex]]+=atom_output
        return energy_pred


def target_scaling(data, method=None):
    if method == 'minmax':
        data_max = max(data)
        data_min = min(data)
        scale = []
        for index, value in enumerate(data):
            scale.append((value-data_min)/(data_max-data_min))
        return torch.stack(scale)
    elif method == 'standardize':
        data_mean = torch.mean(data)
        data_sd = torch.std(data, dim=0)
        scale = []
        for index, value in enumerate(data):
            scale.append((value-data_mean)/data_sd)
        return torch.stack(scale)
    else:
        return data


def pred_scaling(data, target, method=None):
    if method == 'minmax':
        target_max = max(target)
        target_min = min(target)
        scale = []
        for index, value in enumerate(data):
            scale.append((value*(target_max-target_min))+target_min)
        return torch.stack(scale)
    elif method == 'standardize':
        target_mean = torch.mean(target)
        target_sd = torch.std(target, dim=0)
        scale = []
        for index, value in enumerate(data):
            scale.append((value*target_sd)+target_mean)
        return torch.stack(scale)
    else:
        return data


def train_model(model, unique_atoms, dataset_size, criterion, optimizer, atoms_dataloader, num_epochs):
    # device=torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    device = "cpu"

    log = Logger('../benchmark_results/results-log.txt')
    log_epoch = Logger('../benchmark_results/epoch-log.txt')
    log('Model: %s' % model)

    model.train()

    since = time.time()
    log_epoch('-'*50)
    print('Training Initiated!')
    log_epoch('%s Training Initiated!' % time.asctime())
    log_epoch('')

    best_model_wts = copy.deepcopy(model.state_dict())
    best_loss = 1e8

    # plot_epoch_x=list(range(1,num_epochs+1))
    plot_loss_y = []

    epoch = 0
    # for epoch in range(num_epochs):
    while best_loss >= 2e-3:
        log_epoch('{} Epoch {}/{}'.format(time.asctime(), epoch+1, num_epochs))
        log_epoch('-'*30)

        MSE = 0.0

        for data_sample in atoms_dataloader:
            input_data = data_sample[0]
            target = data_sample[1].requires_grad_(False)
            batch_size = len(target)
            target = target.reshape(batch_size, 1)
            scaled_target = target_scaling(target, method='standardize')
            num_of_atoms = data_sample[2].reshape(batch_size, 1)
            # Send inputs and labels to the corresponding device (cpu or gpu)
            for element in unique_atoms:
                input_data[element][0] = input_data[element][0].to(device)
            scaled_target = scaled_target.to(device)

            def closure():
                optimizer.zero_grad()
                output = model(input_data)
                loss = criterion(output, scaled_target)
                loss.backward()
                return loss

            optimizer.step(closure)

            # Perform predictions and compute loss
            with torch.no_grad():
                scaled_preds = model(input_data)
                raw_preds = pred_scaling(
                    scaled_preds, target, method='standardize')
                raw_preds_per_atom = torch.div(raw_preds, num_of_atoms)
                target_per_atom = torch.div(target, num_of_atoms)
                loss = criterion(raw_preds_per_atom, target_per_atom)
                MSE += loss.item()*batch_size

        MSE = MSE/dataset_size
        RMSE = np.sqrt(MSE)
        epoch_loss = RMSE
        print epoch_loss
        plot_loss_y.append(np.log10(RMSE))

        log_epoch('{} Loss: {:.4f}'.format(time.asctime(), epoch_loss))

        if epoch_loss < best_loss:
            best_loss = epoch_loss
            best_model_wts = copy.deepcopy(model.state_dict())
        epoch += 1

    time_elapsed = time.time()-since
    print ('Training complete in {} steps'.format(epoch))
    print('Training complete in {:.0f}m {:.0f}s'.format
          (time_elapsed//60, time_elapsed % 60))

    log('')
    log('Training complete in {:.0f}m {:.0f}s'.format
        (time_elapsed//60, time_elapsed % 60))

    log('Best training loss: {:4f}'.format(best_loss))
    log('')

    plt.title('RMSE vs. Epoch')
    plt.xlabel('Epoch #')
    plt.ylabel('log(RMSE)')
    plot_epoch_x = list(range(1, epoch+1))
    plt.plot(plot_epoch_x, plot_loss_y, label='train')
    plt.legend()
    # plt.show()
    plt.savefig('test.png')
    model.load_state_dict(best_model_wts)
    return model
