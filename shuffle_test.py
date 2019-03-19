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
from torch.nn import Tanh
import torch.optim as optim
from torch.optim import lr_scheduler

from ase.build import molecule
from ase import Atoms
from ase.calculators.emt import EMT

filename='benchmark_dataset/water.extxyz'
training_data=AtomsDataset(filename,descriptor=Gaussian())
sample_batch=[training_data[0],training_data[1],training_data[2],training_data[3],training_data[4]]
atoms_loader=DataLoader(sample_batch,batch_size=2,collate_fn=collate_amp,shuffle=False)
unique_atoms,_,_,_=data_factorization(training_data)

'''A test to ensure the shuffling step in the model is working properly.'''
for i in atoms_loader:
    k=i
    ox=k[0]['O'][0]
    hy=k[0]['H'][0]

    model=FullNN(unique_atoms)
    l=list(model.parameters())

    w0=torch.transpose(l[0],0,1)
    b0=l[1]
    ox_con=Tanh()(torch.mm(ox,w0)+b0)
    w1=torch.transpose(l[2],0,1)
    b1=l[3]
    ox_con=Tanh()(torch.mm(ox_con,w1)+b1)
    w2=torch.transpose(l[4],0,1)
    b2=l[5]
    ox_con=torch.mm(ox_con,w2)+b2
    ox_total=ox_con

    hy_total=[]
    for h in hy:
        wh0=torch.transpose(l[6],0,1)
        bh0=l[7]
        hy_con=Tanh()(torch.matmul(h,wh0)+bh0)
        wh1=torch.transpose(l[8],0,1)
        bh1=l[9]
        hy_con=Tanh()(torch.matmul(hy_con,wh1)+bh1)
        wh2=torch.transpose(l[10],0,1)
        bh2=l[11]
        hy_con=torch.matmul(hy_con,wh2)+bh2
        hy_total.append(hy_con)
    hy_indiv_total=torch.stack(hy_total)

    hy_total=[]
    for i in range(len(hy))[::2]:
        hy_total.append(hy_indiv_total[i]+hy_indiv_total[i+1])
    hy_total=torch.stack(hy_total)

    manual_total=(hy_total+ox_total).detach()
    model_total=model(k[0]).detach()
    manual_total=manual_total.tolist()
    model_total=model_total.tolist()
    print manual_total==model_total
    print 'Expected Result:\n %s'%manual_total
    print 'Model Result: \n %s \n'%model_total
