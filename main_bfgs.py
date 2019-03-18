import sys
import time
import torch
import torch.nn as nn
import numpy as np
from torch.utils.data import DataLoader
from data import AtomsDataset,data_factorization,collate_amp
from amp.utilities import Logger
from amp.descriptor.gaussian import Gaussian
from nn_torch_bfgs import FullNN,train_model
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


batch_size=400
log('Batch Size = %d'%batch_size)
validation_frac=0

if validation_frac!=0:
    samplers=training_data.data_split(training_data,validation_frac)
    dataset_size={'train':(1.-validation_frac)*len(training_data),'val':validation_frac*len(training_data)}
    log('Training Data = %d Validation Data = %d'
            %(dataset_size['train'],dataset_size['val']))
    atoms_dataloader={x:DataLoader(training_data,batch_size,collate_fn=collate_amp,sampler=samplers[x],pin_memory=True) for x in ['train','val']}

else:
    dataset_size=len(training_data)
    log('Training Data = %d'%dataset_size)
    atoms_dataloader=DataLoader(training_data,batch_size,collate_fn=collate_amp,shuffle=True,pin_memory=True)

#Check SD of targets
# for i in atoms_dataloader:
    # print i[1].std(dim=0)

model=FullNN(unique_atoms)
# if torch.cuda.device_count()>1:
    # print('Utilizing',torch.cuda.device_count(),'GPUs!')
    # model=nn.DataParallel(model)
model=model.to(device)
criterion=nn.MSELoss()
log('Loss Function: %s'%criterion)

#Define the optimizer and implement any optimization settings
# optimizer_ft=optim.SGD(model.parameters(),lr=.01,momentum=0.9)
optimizer_ft=optim.LBFGS(model.parameters())

log('Optimizer Info:\n %s'%optimizer_ft)

#Define scheduler search strategies
# exp_lr_scheduler=lr_scheduler.StepLR(optimizer_ft,step_size=20,gamma=0.1)
# log('LR Scheduler Info: \n Step Size = %s \n Gamma = %s'%(exp_lr_scheduler.step_size,exp_lr_scheduler.gamma))

num_epochs=100
log('Number of Epochs = %d'%num_epochs)
log('')
model=train_model(model,unique_atoms,dataset_size,criterion,optimizer_ft,atoms_dataloader,num_epochs)
torch.save(model.state_dict(),'benchmark_results/benchmark_model.pt')

