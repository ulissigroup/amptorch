import os
import pickle

import numpy as np
from amp.descriptor.gaussian import Gaussian as Gnew
from amp.utilities import hash_images as amp_hash
from ase import Atoms
from ase.calculators.emt import EMT

from amptorch.dataset_backup import AtomsDataset as backupdataset
from amptorch.descriptor.descriptor_calculator import DescriptorCalculator
from amptorch.descriptor.Gaussian import Gaussian

distances = np.linspace(2, 5, 10)
images = []
for l in distances:
    image = Atoms(
        "CuCO",
        [
            (-l * np.sin(0.65), l * np.cos(0.65), 0),
            (0, 0, 0),
            (l * np.sin(0.65), l * np.cos(0.65), 0),
        ],
    )
    image.set_cell([10, 10, 10])
    image.wrap(pbc=True)
    image.set_calculator(EMT())
    images.append(image)


Gs_new = {
    "default": {
        "G2": {
            "etas": np.logspace(np.log10(0.05), np.log10(5.0), num=4) / 6 ** 2,
            "rs_s": [0],
        },
        "G4": {"etas": [0.005], "zetas": [1.0, 4.0], "gammas": [1.0, -1.0]},
        "cutoff": 6,
    },
}
# Old
Gs = {}
Gs["G2_etas"] = np.logspace(np.log10(0.05), np.log10(5.0), num=4)
Gs["G2_rs_s"] = [0] * 4
Gs["G4_etas"] = [0.005]
Gs["G4_zetas"] = [1.0, 4.0]
Gs["G4_gammas"] = [+1.0, -1]
Gs["cutoff"] = 6

forcetraining = True
dataset = backupdataset(
    images, Gnew, Gs, forcetraining=forcetraining, label="test", cores=4,
)


elements = image.get_chemical_symbols()

descriptor = Gaussian(Gs=Gs_new, elements=elements)
desc_calc = DescriptorCalculator(
    images=images, descriptor=descriptor, calc_derivatives=True, save_fps=True, cores=1,
)
new = desc_calc.prepare_descriptors()

# Old
Gs = {}
Gs["G2_etas"] = np.logspace(np.log10(0.05), np.log10(5.0), num=4)
Gs["G2_rs_s"] = [0] * 4
Gs["G4_etas"] = [0.005]
Gs["G4_zetas"] = [1.0, 4.0]
Gs["G4_gammas"] = [+1.0, -1]
Gs["cutoff"] = 6

image_hashes = list(amp_hash(images).keys())
for i in range(len(images)):
    image_hash = image_hashes[i]
    with open("amp-data-fingerprints.ampdb/loose/" + image_hash, "rb") as f:
        amp_fp = pickle.load(f)
    os.system("rm amp-data-fingerprints.ampdb/loose/" + image_hash)

    with open("amp-data-fingerprint-primes.ampdb/loose/" + image_hash, "rb") as f:
        amp_prime = pickle.load(f)
    os.system("rm amp-data-fingerprint-primes.ampdb/loose/" + image_hash)
