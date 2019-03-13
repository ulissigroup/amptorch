import pytest
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

filename='benchmark_dataset/water.extxyz'
training_data=AtomsDataset(filename,descriptor=Gaussian())
sample_batch=[training_data[0],training_data[1],training_data[2]]
atoms_loader=DataLoader(sample_batch,batch_size=3,collate_fn=collate_amp,shuffle=True)
unique_atoms,_,_,_=data_factorization(training_data)
for i in atoms_loader:
    k=i


# print k

# print k[1]
ox=k[0]['O'][0]
hy=k[0]['H'][0]
model=FullNN(unique_atoms)
l=list(model.parameters())
w0=torch.transpose(l[0],0,1)
b0=l[1]
ox_con=torch.mm(ox,w0)+b0
w1=torch.transpose(l[2],0,1)
b1=l[3]
ox_con=torch.mm(ox_con,w1)+b1
w2=torch.transpose(l[4],0,1)
b2=l[5]
ox_con=torch.mm(ox_con,w2)+b2
ox_total=ox_con

hy_total=[]
for h in hy:
    wh0=torch.transpose(l[6],0,1)
    bh0=l[7]
    hy_con=torch.matmul(h,wh0)+bh0
    wh1=torch.transpose(l[8],0,1)
    bh1=l[9]
    hy_con=torch.matmul(hy_con,wh1)+bh1
    wh2=torch.transpose(l[10],0,1)
    bh2=l[11]
    hy_con=torch.matmul(hy_con,wh2)+bh2
    hy_total.append(hy_con)
hy_indiv=torch.stack(hy_total)
hy_total=[]
for i in range(6)[::2]:
    hy_total.append(hy_indiv[i]+hy_indiv[i+1])
hy_total=torch.stack(hy_total)
# print ox_total
# print hy_total
total=hy_total+ox_total
print total
print model(k[0])
print k[1]



