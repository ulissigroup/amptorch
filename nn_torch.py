import sys
import numpy as np
import torch
import torch.nn as nn
from torch.nn import init
from torch.nn import Tanh, LeakyReLU
import copy
from collections import OrderedDict
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

    (Simplified version of SchNet Pack's implementation')

    """

    def __init__(self,input_size,output_size, bias=True, activation=None):
        self.activation=activation
        super(Dense,self).__init__(input_size,output_size,bias)

    def reset_parameters(self):
        # init.constant_(self.weight,.5)
        # init.constant_(self.bias,1)

        super(Dense,self).reset_parameters()
        # init.constant_(self.bias,1)

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

    def __init__(self,n_input_nodes=20,n_output_nodes=1,n_layers=3,n_hidden_size=5,activation=Tanh):
        super(MLP,self).__init__()
        #if n_hidden_size is None:
            #implement pyramid neuron structure across each layer
        if type(n_hidden_size) is int:
            n_hidden_size=[n_hidden_size] * (n_layers-1)
        self.n_neurons=[n_input_nodes]+n_hidden_size+[n_output_nodes]
        self.activation=activation
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
        log=Logger('benchmark_results/results-log.txt')

        super(FullNN,self).__init__()
        self.unique_atoms=unique_atoms
        n_unique_atoms=len(unique_atoms)
        elementwise_models=nn.ModuleList()
        for n_models in range(n_unique_atoms):
            elementwise_models.append(MLP())
        self.elementwise_models=elementwise_models
        log('Activation Function = %s'%elementwise_models[0].activation)
        # self.model_device=next(elementwise_models.parameters()).device

    def forward(self,data):
        energy_pred=OrderedDict()
        for index,element in enumerate(self.unique_atoms):
            model_inputs=data[element][0]
            contribution_index=data[element][1]
            atomwise_outputs=self.elementwise_models[index].forward(model_inputs)
            for index,atom_output in enumerate(atomwise_outputs):
                if contribution_index[index] not in energy_pred.keys():
                    energy_pred[contribution_index[index]]=atom_output
                else:
                    energy_pred[contribution_index[index]]+=atom_output
        output=energy_pred.values()
        output=torch.stack(output)
        return output

def feature_scaling(data):
    data_max=max(data)
    data_min=min(data)
    scale=[]
    for index,value in enumerate(data):
        scale.append((value-data_min)/(data_max-data_min))
    return torch.stack(scale)

def train_model(model,unique_atoms,dataset_size,criterion,optimizer,scheduler,atoms_dataloader,num_epochs):
    log=Logger('benchmark_results/results-log.txt')
    log_epoch=Logger('benchmark_results/epoch-log.txt')
    log('Model: %s'%model)

    since = time.time()
    log_epoch('-'*50)
    print('Training Initiated!')
    log_epoch('%s Training Initiated!'%time.asctime())
    log_epoch('')

    best_model_wts=copy.deepcopy(model.state_dict())
    best_loss=100000000

    plot_epoch_x=list(range(1,num_epochs+1))
    plot_loss_y={'train':[],'val':[]}

    for epoch in range(num_epochs):
        log_epoch('{} Epoch {}/{}'.format(time.asctime(),epoch+1,num_epochs))
        log_epoch('-'*30)

        if len(atoms_dataloader)==2:

                for phase in ['train','val']:

                    if phase == 'train':
                        scheduler.step()
                        model.train()
                    else:
                        model.eval()

                    MSE=0.0

                    #Iterate over the dataloader
                    for data_sample in atoms_dataloader[phase]:
                        input_data=data_sample[0]
                        target=data_sample[1]

                        #Send inputs and labels to the corresponding device (cpu or gpu)
                        for element in unique_atoms:
                            input_data[element][0]=input_data[element][0].to(device)
                        target=target.to(device)

                        #zero the parameter gradients
                        optimizer.zero_grad()

                        #forward
                        with torch.set_grad_enabled(phase == 'train'):
                            output=model(input_data)
                            # print output.size()
                            _,preds = torch.max(output,1)
                            loss=criterion(output,target)
                            #backward + optimize only if in training phase
                            if phase == 'train':
                                loss.backward()
                                optimizer.step()

                        MSE += loss.item()

                    MSE=MSE/dataset_size[phase]
                    RMSE=np.sqrt(MSE)
                    epoch_loss = RMSE
                    plot_loss_y[phase].append(np.log10(RMSE))

                    log_epoch('{} {} Loss: {:.4f}'.format(time.asctime(),phase,epoch_loss))

                    if phase == 'val' and epoch_loss<best_loss:
                        best_loss=epoch_loss
                        best_model_wts=copy.deepcopy(model.state_dict())
        else:

            phase='train'

            scheduler.step()
            model.train()

            MSE=0.0

            #Iterate over the dataloader
            for data_sample in atoms_dataloader:
                input_data=data_sample[0]
                target=data_sample[1]
                # target=feature_scaling(target)
                # target=target/2000.

                #Send inputs and labels to the corresponding device (cpu or gpu)
                for element in unique_atoms:
                    input_data[element][0]=input_data[element][0].to(device)
                target=target.to(device)
                #zero the parameter gradients
                optimizer.zero_grad()

                #forward
                with torch.set_grad_enabled(True):
                    output=model(input_data)
                    # print output.size()
                    _,preds = torch.max(output,1)
                    loss=criterion(output,target)
                    #backward + optimize only if in training phase
                    loss.backward()
                    optimizer.step()

                MSE += loss.item()

            MSE=MSE/dataset_size
            RMSE=np.sqrt(MSE)
            epoch_loss = RMSE
            print epoch_loss
            plot_loss_y[phase].append(np.log10(RMSE))

            log_epoch('{} {} Loss: {:.4f}'.format(time.asctime(),phase,epoch_loss))

            if epoch_loss<best_loss:
                best_loss=epoch_loss
                best_model_wts=copy.deepcopy(model.state_dict())

        log_epoch('')

    time_elapsed=time.time()-since
    print('Training complete in {:.0f}m {:.0f}s'.format
            (time_elapsed//60,time_elapsed % 60))

    log('')
    log('Training complete in {:.0f}m {:.0f}s'.format
                (time_elapsed//60,time_elapsed % 60))

    if phase=='val':
        log('Best validation loss: {:4f}'.format(best_loss))
    else:
        log('Best training loss: {:4f}'.format(best_loss))

    log('')

    plt.title('RMSE vs. Epoch')
    plt.xlabel('Epoch #')
    plt.ylabel('log(RMSE)')
    plt.plot(plot_epoch_x,plot_loss_y['train'],label='train')
    if phase=='val':
        plt.plot(plot_epoch_x,plot_loss_y['val'],label='val')
    plt.legend()
    plt.show()

    model.load_state_dict(best_model_wts)
    return model

device=torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
