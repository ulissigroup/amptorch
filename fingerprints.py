import numpy as np
from ase.structure import molecule
from ase import Atoms
atoms=molecule('H2O')
atoms.rotate('y',-np.pi/2.)
atoms.set_pbc(False)
displacements=np.linspace(0.9,8.0,20)
vec=atoms[2].position - atoms[0].position
images=[]
for displacement in displacements:
    atoms=Atoms(atoms)
    atoms[2].position=(atoms[0].position +vec*displacement)
    images.append(atoms)
#print(images)
from amp.descriptor.gaussian import Gaussian
descriptor = Gaussian()
from amp.utilities import hash_images
images = hash_images(images, ordered=True)
descriptor.calculate_fingerprints(images)

from matplotlib import pyplot


