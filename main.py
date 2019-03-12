import sys
import time
import torch
import torch.nn as nn
import numpy as np
from torch.utils.data import DataLoader
from data import AtomsDataset,data_factorization,collate_amp
from amp.utilities import Logger
from amp.descriptor.gaussian import Gaussian
from nn_torch import FullNN,train_model
import torch.optim as optim
from torch.optim import lr_scheduler

from ase.build import molecule
from ase import Atoms
from ase.calculators.emt import EMT

log=Logger('benchmark_results/results-log.txt')
log_epoch=Logger('benchmark_results/epoch-log.txt')

log(time.asctime())

device=torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
filename='benchmark_dataset/water.extxyz'

log('-'*50)
log('Filename: %s'%filename)

training_data=AtomsDataset(filename,descriptor=Gaussian())
unique_atoms,_,_,_=data_factorization(training_data)
n_unique_atoms=len(unique_atoms)

validation_frac=.1
samplers=training_data.data_split(training_data,validation_frac)
dataset_sizes={'train':(1.-validation_frac)*len(training_data),'val':validation_frac*len(training_data)}

log('Training Data = %d Validation Data = %d'
        %(dataset_sizes['train'],dataset_sizes['val']))

batch_size=100
atoms_dataloaders={x:DataLoader(training_data,batch_size,collate_fn=collate_amp,sampler=samplers[x],pin_memory=True)
        for x in ['train','val']}
log('Batch Size = %d'%batch_size)

model=FullNN(unique_atoms)
model=model.to(device)
criterion=nn.MSELoss()
log('Loss Function: %s'%criterion)

#Define the optimizer and implement any optimization settings
optimizer_ft=optim.SGD(model.parameters(),lr=.01,momentum=0.9)
# optimizer_ft=optim.Adam(model.parameters(),lr=.1)
# optimizer_ft=optim.SGD(model.parameters(),lr=.001)
log('Optimizer Info:\n %s'%optimizer_ft)

#Define scheduler search strategies
exp_lr_scheduler=lr_scheduler.StepLR(optimizer_ft,step_size=20,gamma=0.1)
log('LR Scheduler Info: \n Step Size = %s \n Gamma = %s'%(exp_lr_scheduler.step_size,exp_lr_scheduler.gamma))

num_epochs=100
log('Number of Epochs = %d'%num_epochs)
log('')
model=train_model(model,unique_atoms,dataset_sizes,criterion,optimizer_ft,exp_lr_scheduler,atoms_dataloaders,num_epochs)
torch.save(model.state_dict(),'benchmark_results/benchmark_model.pt')






def test_model():
    model=FullNN(unique_atoms)
    model=model.to(device)
    model.load_state_dict(torch.load('benchmark_model.pt'))

    for data_sample in atoms_dataloader:
        outputs,target=model(data_sample)
        print outputs
        print target
        print('')
