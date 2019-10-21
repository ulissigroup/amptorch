import sys
import ase
import ase.io
from ase.visualize import view

k = ase.io.read('50.traj',':')

for i in k:
    view(i)
    sys.exit()
