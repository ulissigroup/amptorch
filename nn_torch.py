import sys
import numpy as np
import torch
import torch.nn as nn
from torch.nn import init
from torch.nn import Tanh
import copy
from collections import OrderedDict
import time
import matplotlib.pyplot as plt

class Dense(nn.Linear):
    """Constructs and applies a dense layer with an activation function (when
    available) y=activation(Ax+b)

    Arguments:
        input_size (int): number of input features
        output_size (int): number of output features
        bias (bool): True to include bias at each neuron. False otherwise
        (Default: True)
        activation (callable): activation function to be utilized (Default:None)

    (Simplified version of SchNet Pack's implementation')

    """

    def __init__(self,input_size,output_size, bias=True, activation=None):
        self.activation=activation
        super(Dense,self).__init__(input_size,output_size,bias)

    def reset_parameters(self):
        # init.constant_(self.weight,1)
        # init.constant_(self.bias,0)

        super(Dense,self).reset_parameters()

    def forward(self,inputs):
        neuron_output=super(Dense,self).forward(inputs)
        if self.activation:
            neuron_output=self.activation()(neuron_output)
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

    (Simplified version of SchNet Pack's implementation')

    """

    def  __init__(self,n_input_nodes=20,n_output_nodes=1,n_layers=3,n_hidden_size=10,activation=Tanh):
        super(MLP,self).__init__()
        #if n_hidden_size is None:
            #implement pyramid neuron structure across each layer
        if type(n_hidden_size) is int:
            n_hidden_size=[n_hidden_size] * (n_layers-1)
        self.n_neurons=[n_input_nodes]+n_hidden_size+[n_output_nodes]

        layers=[Dense(self.n_neurons[i],self.n_neurons[i+1],bias=True,activation=activation) for i in range(n_layers-1)]
        layers.append(Dense(self.n_neurons[-2],self.n_neurons[-1],activation=None))
        self.model_net=nn.Sequential(*layers)

    def forward(self, inputs):
        """Feeds data forward in the neural network

        Arguments:
            inputs (torch.Tensor): NN inputs
        """

        return self.model_net(inputs)

class FullNN(nn.Module):
    '''Combines element specific NNs into a model to predict energy of a given
    structure

    '''

    def __init__(self,unique_atoms):
        super(FullNN,self).__init__()
        self.unique_atoms=unique_atoms
        n_unique_atoms=len(unique_atoms)
        elementwise_models=nn.ModuleList()
        for n_models in range(n_unique_atoms):
            elementwise_models.append(MLP())
        self.elementwise_models=elementwise_models
        self.model_device=next(elementwise_models.parameters()).device

    def forward(self,data):
        energy_pred=OrderedDict()
        for index,element in enumerate(self.unique_atoms):
            model_inputs=data[element][0]
            model_inputs=model_inputs.to(self.model_device)
            contribution_index=data[element][1]
            targets=data[element][2][0]
            targets=targets.to(self.model_device)
            atomwise_outputs=self.elementwise_models[index].forward(model_inputs)
            for index,atom_output in enumerate(atomwise_outputs):
                if contribution_index[index] not in energy_pred.keys():
                    energy_pred[contribution_index[index]]=atom_output
                else:
                    energy_pred[contribution_index[index]]+=atom_output
        output=energy_pred.values()
        output=torch.stack(output)
        return output,targets

def train_model(model,criterion,optimizer,scheduler,atoms_dataloaders,num_epochs):
    since = time.time()

    best_model_wts=copy.deepcopy(model.state_dict())
    best_loss=100000000

    plot_epoch_x=list(range(1,num_epochs+1))
    plot_loss_y={'train':[],'val':[]}

    for epoch in range(num_epochs):
        print ('Epoch {}/{}'.format(epoch+1,num_epochs))
        print('-'*10)

        for phase in ['train','val']:

            if phase == 'train':
                scheduler.step()
                model.train()
            else:
                model.eval()

            MSE=0.0

            #Iterate over the dataloader
            for data_sample in atoms_dataloaders[phase]:
                #zero the parameter gradients
                optimizer.zero_grad()

                #forward
                with torch.set_grad_enabled(phase == 'train'):
                    outputs,target=model(data_sample)
                    _,preds = torch.max(outputs,1)
                    loss=criterion(outputs,target)

                    #backward + optimize only if in training phase
                    if phase == 'train':
                        loss.backward()
                        optimizer.step()

                MSE += loss.item()

            RMSE=np.sqrt(MSE)
            epoch_loss = RMSE
            plot_loss_y[phase].append(RMSE)

            print('{} Loss: {:.4f}'.format(phase,epoch_loss))

            if phase == 'val' and epoch_loss<best_loss:
                best_loss=epoch_loss
                best_model_wts=copy.deepcopy(model.state_dict())

            # if epoch_loss<best_loss:
                # best_loss=epoch_loss
                # best_model_wts=copy.deepcopy(model.state_dict())

        print('')

    time_elapsed=time.time()-since
    print('Training complete in {:.0f}m {:.0f}s'.format
            (time_elapsed//60,time_elapsed % 60))

    print('Best validation loss: {:4f}'.format(best_loss))

    plt.title('RMSE vs. Epoch')
    plt.xlabel('Epoch #')
    plt.ylabel('RMSE')
    plt.plot(plot_epoch_x,plot_loss_y['train'],label='train')
    plt.plot(plot_epoch_x,plot_loss_y['val'],label='val')
    plt.legend()
    plt.show()

    model.load_state_dict(best_model_wts)
    return model

device=torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

