import ase
from ase import io
from amp.utilities import check_images,hash_images
import sys

l=ase.io.read('water.extxyz',':')
test=l[0]
print test.get_fingerprints()
print l[0]
# print l

