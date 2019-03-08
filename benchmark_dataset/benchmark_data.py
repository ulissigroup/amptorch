import ase
from ase import io
from amp.utilities import check_images,hash_images
import sys

k=ase.io.read('sample_training_data.traj',':')
print k[0]
k=hash_images(k)
# print k
l=ase.io.read('water.extxyz',':')
print l[0]
l=hash_images(l)
print l[0]





traj=ase.io.Trajectory('test','w')
