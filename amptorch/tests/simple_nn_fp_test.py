from amptorch.utils import make_amp_descriptors_simple_nn
import os
from amp.descriptor.gaussian import Gaussian, make_symmetry_functions
from amp.descriptor.cutoffs import Cosine
from ase.calculators.singlepoint import SinglePointCalculator as sp
from ase.build import molecule
import numpy as np
from pickle import load
from amp.utilities import hash_images as stock_hash
from amptorch.utils import hash_images as new_hash


atoms = molecule('H2O')
print(atoms)
atoms.set_cell([10,10,10])
atoms.set_pbc = [True] * 3

atoms.set_calculator(sp(atoms=atoms, energy = -1,
                        forces = np.array([[-1,-1,-1],[-1,-1,-1]])))

Gs = {}
images = [atoms]
Gs["G2_etas"] = [0.005]
Gs["G2_rs_s"] = [0]
Gs["G4_etas"] = [0.005]
Gs["G4_zetas"] = [1., 4.]
Gs["G4_gammas"] = [1., -1.]
Gs["cutoff"] = 6.5


elements = ['O', 'H']
G = make_symmetry_functions(elements=elements, type='G2',
                            etas=[0.005])
G += make_symmetry_functions(elements=elements, type='G4',
                             etas=[0.005],
                             zetas=[1., 4.],
                             gammas=[+1., -1.])
G = {'O': G,
     'H':G
     }

hashes = stock_hash(images)

os.system('rm amp-data-fingerprints.ampdb/loose/*')
make_amp_descriptors_simple_nn(images, Gs, elements=elements)


with open('amp-data-fingerprints.ampdb/loose/'+new_hash(atoms, Gs), 'rb') as f:
    simple_nn = load(f)
os.system('rm amp-data-fingerprints.ampdb/loose/'+new_hash(atoms, Gs))

descriptor = Gaussian(elements=['O', 'H'], Gs=G, cutoff=Cosine(Gs["cutoff"]))
descriptor.calculate_fingerprints(hashes, calculate_derivatives=True)
with open('amp-data-fingerprints.ampdb/loose/'+stock_hash(atoms), 'rb') as f:
    amp = load(f)
os.system('rm amp-data-fingerprints.ampdb/loose/'+stock_hash(atoms))

print(simple_nn)
print(amp)

for s,am in zip(simple_nn, amp):
    pass
    #print(np.array(s[1]) - np.array(am[1]))
    #print(s[1])
    #print(am[1])
    #if sum(np.array(s[1]) - np.array(am[1])) != 0.:
    #    raise Exception('simple_nn amp fingerprint comparison failed')


