import os
import sys
import numpy as np
from ase import Atoms
import torch
import torch.nn as nn
from torch.nn import init
from torch.nn import Tanh
import copy
import math
from collections import OrderedDict
import time
"""Loading Data"""
from torch.utils.data import Dataset, DataLoader
from amp.utilities import hash_images
from amp.utilities import check_images
from amp.descriptor.gaussian import Gaussian
'''Trainer'''
import torch.optim as optim
from torch.optim import lr_scheduler

class AtomsDataset(Dataset):
    """
    Atoms dataset
    Parameters: Descriptor type and .traj file name
    Output: Returns, for a given index, the image_fingerprint and image_potential
    energy
    """

    def __init__(self,descriptor,filename='sample_training_data.traj'):
        self.filename=filename
        self.descriptor=descriptor
        self.atom_images=hash_images(filename)
        check_images(self.atom_images,forces=False)
        self.descriptor.calculate_fingerprints(self.atom_images)

    def __len__(self):
        return len(self.atom_images)

    def __getitem__(self,index):
        hash_name=self.atom_images.keys()[index]
        image_fingerprint=self.descriptor.fingerprints[hash_name]
        image_potential_energy=self.atom_images[hash_name].get_potential_energy()
        return {index:(image_fingerprint,image_potential_energy)}

def data_factorization(training_data):
    """
    Reads in dataset and factors it into 4 lists:

    1. unique_atoms = Identifies the unique elements in the dataset
    2. fingerprint_dict = Extracts the fingerprints for each hashed data sample in the
    dataset
    3. energy_dict = Extracts the potential energy for a given hashed data sample in the
    dataset
    4. sample_indices = Identifies indices of corresponding fingerprints
    """
    unique_atoms=[]
    fingerprint_dataset=[]
    sample_indices=[]
    energy_dataset=[]
    #Create empty dictionary to store indices of data
    for data_sample in training_data:
        sample_index=data_sample.keys()[0]
        sample_indices.append(sample_index)
        atom_image=data_sample[sample_index]
        atom_fingerprint=atom_image[0]
        fingerprint_dataset.append(atom_fingerprint)
        image_potential_energy=atom_image[1]
        energy_dataset.append(image_potential_energy)
        for atom in atom_fingerprint:
            element=atom[0]
            if element not in unique_atoms:
                unique_atoms.append(element)
    return unique_atoms,fingerprint_dataset,energy_dataset,sample_indices

def collate_amp(training_data):
    unique_atoms,fingerprint_dataset,energy_dataset,sample_indices=data_factorization(training_data)
    # print energy_dataset
    element_specific_fingerprints={}
    for element in unique_atoms:
        element_specific_fingerprints[element]=[[],[],[]]
    for fp_index, sample_fingerprints in enumerate(fingerprint_dataset):
        for fingerprint in sample_fingerprints:
            atom_element=fingerprint[0]
            atom_fingerprint=fingerprint[1]
            element_specific_fingerprints[atom_element][0].append(torch.tensor(atom_fingerprint))
            element_specific_fingerprints[atom_element][1].append(sample_indices[fp_index])
            #INSERT SCALING OF INPUT DATA
    for element in unique_atoms:
        element_specific_fingerprints[element][0]=torch.stack(element_specific_fingerprints[element][0])
        element_specific_fingerprints[element][2].append(torch.tensor(energy_dataset))
    return element_specific_fingerprints

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
        '''For testing purposes let weights initialize to 1 and bias to 0'''

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

    def __init__(self,n_input_nodes=20,n_output_nodes=1,n_layers=2,n_hidden_size=2,activation=Tanh):
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

    def forward(self,data):
        energy_pred=OrderedDict()
        for index,element in enumerate(self.unique_atoms):
            model_inputs=data[element][0]
            model_inputs=model_inputs.to(device)
            contribution_index=data[element][1]
            targets=data[element][2][0]
            targets=targets.to(device)
            atomwise_outputs=self.elementwise_models[index].forward(model_inputs)
            for index,atom_output in enumerate(atomwise_outputs):
                if contribution_index[index] not in energy_pred.keys():
                    energy_pred[contribution_index[index]]=atom_output
                else:
                    energy_pred[contribution_index[index]]+=atom_output
        output=energy_pred.values()
        output=torch.stack(output)
        return output,targets

# device=torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

# training_data=AtomsDataset(descriptor=Gaussian())
# sample_batch=[training_data[1], training_data[0], training_data[3],
        # training_data[18] ]
# unique_atoms,_,_,_=data_factorization(sample_batch)
# n_unique_atoms=len(unique_atoms)

# atoms_dataloader=DataLoader(sample_batch,batch_size=2,collate_fn=collate_amp,drop_last=True,shuffle=False)
# for i in atoms_dataloader:
    # k=i

# print k

# model=FullNN(unique_atoms)
# output,targets=model(k)
# print output
# print targets

# print torch.max(output,1)
# loss=nn.MSELoss()
# test=loss(output,targets)

def train_model(model,criterion,optimizer,scheduler,num_epochs):
    since = time.time()

    best_model_wts=copy.deepcopy(model.state_dict())
    best_loss=100000000

    for epoch in range(num_epochs):
        print ('Epoch {}/{}'.format(epoch+1,num_epochs))
        print('-'*10)

        scheduler.step()
        model.train()

        running_loss=0.0

        #Iterate over the dataloader
        for data_sample in atoms_dataloader:
            #zero the parameter gradients
            optimizer.zero_grad()
            #forward
            with torch.set_grad_enabled(True):
                outputs,target=model(data_sample)
                _,preds = torch.max(outputs,1)
                loss=criterion(outputs,target)

                #backward + optimize only if in training phase
                loss.backward()
                optimizer.step()

            running_loss +=loss.item()

        epoch_loss= running_loss

        print('train Loss: {:.4f}'.format(epoch_loss))

        if epoch_loss<best_loss:
            best_loss=epoch_loss
            best_model_wts=copy.deepcopy(model.state_dict())
        print('')

    time_elapsed=time.time()-since
    print('Training complete in {:.0f}m {:.0f}s'.format
            (time_elapsed//60,time_elapsed % 60))

    print('Best train loss: {:4f}'.format(best_loss))

    model.load_state_dict(best_model_wts)
    return model

device=torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
training_data=AtomsDataset(descriptor=Gaussian())
# sample_batch=[training_data[1], training_data[0], training_data[3],
        # training_data[18] ]
unique_atoms,_,_,_=data_factorization(training_data)
n_unique_atoms=len(unique_atoms)

atoms_dataloader=DataLoader(training_data,batch_size=2,collate_fn=collate_amp,drop_last=True,shuffle=False)


def main():

    model=FullNN(unique_atoms)
    model=model.to(device)
    criterion=nn.MSELoss()

    #Optimize all parameters
    optimizer_ft=optim.SGD(model.parameters(),lr=.001,momentum=0.01)

    #Decay LR over time (factor of .1 every 7 seconds)
    exp_lr_scheduler=lr_scheduler.StepLR(optimizer_ft,step_size=7,gamma=0.1)

    model=train_model(model,criterion,optimizer_ft,exp_lr_scheduler,num_epochs=100)
    torch.save(model.state_dict(),'Atomistic_model.pt')
    # loss=nn.MSELoss()
    # test=loss(output,targets)


if __name__ == '__main__':
    main()
































