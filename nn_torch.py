import os
import sys
import numpy as np
from ase import Atoms
import torch
import torch.nn as nn
from torch.nn import Tanh
import copy
"""Loading Data"""
from torch.utils.data import Dataset, DataLoader
from amp.utilities import hash_images
from amp.utilities import check_images
from amp.descriptor.gaussian import Gaussian
from collections import defaultdict



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

training_data=AtomsDataset(descriptor=Gaussian())
sample_batch=[training_data[1], training_data[0], training_data[3],
        training_data[18] ]
# print(sample_batch)

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
    element_specific_fingerprints={}
    for element in unique_atoms:
        element_specific_fingerprints[element]=[[],[]]
    for fp_index, sample_fingerprints in enumerate(fingerprint_dataset):
        for fingerprint in sample_fingerprints:
            atom_element=fingerprint[0]
            atom_fingerprint=fingerprint[1]
            element_specific_fingerprints[atom_element][0].append(torch.tensor(atom_fingerprint))
            element_specific_fingerprints[atom_element][1].append(torch.tensor(sample_indices[fp_index]))
            #INSERT SCALING OF INPUT DATA
    return element_specific_fingerprints

test=collate_amp(sample_batch)

atoms_dataloader=DataLoader(sample_batch,batch_size=1,collate_fn=collate_amp,shuffle=False)
for i in atoms_dataloader:
    k=i
# """TENSOR SIZING? MAY NEED TO MAKE EACH DATA ENTRY AN ARRAY"""
l=k['Cu'][0][0]
print(l)
# x,_=torch.max(l,0)
# lx=torch.div(l,x)

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
        super(Dense,self).reset_parameters()

    def forward(self,inputs):
        neuron_output=super(Dense,self).forward(inputs)
        if self.activation:
            out=self.activation(out)

        return neuron_output

class AtomisticNN(nn.Module):
    """Constructs a fully connected neural network model

    Arguments:
        n_input_nodes: Number of input nodes (Default=20 using BP SF)
        n_output_nodes: Number of output nodes (Default=1)
        n_layers: Total number of layers in the neural network
        n_hidden_size: Number of neurons within each hidden layer
        activation: Activation function to be utilized. (Default=Tanh())
    """


    def __init__(self,n_input_nodes=20,n_output_nodes=1,n_layers=2,n_hidden_size=1,activation=Tanh()):
        super(AtomisticNN,self).__init__()
        # if n_hidden_size is None:
            #implement pyramid neuron structure across each layer
        if type(n_hidden_size) is int:
            n_hidden_size=[n_hidden_size] * (n_layers-1)
        self.n_neurons=[n_input_nodes]+n_hidden_size+[n_output_nodes]

        layers=[Dense(self.n_neurons[i],self.n_neurons[i+1],bias=True,activation=activation()) for i in range(n_layers-1)]
        layers.append(Dense(self.n_neurons[-2],self.n_neurons[-1],activation=None))

        self.model_net=nn.Sequential(*layers)

    def forward(self, inputs):
        """Feeds data forward in the neural network

        Arguments:
            inputs (torch.Tensor): NN inputs

        """

        return self.model_net(inputs)

AtomisticNN()


#atoms_dataloader=DataLoader(test_data,batch_size=1,shuffle=False,sampler=None,batch_sampler=None,num_workers=0,collate_fn=,pin_memory=False,drop_last=False,timeout=0,worker_init_fn=None)
