import os
import pickle
import sys
import numpy as np
from amptorch.dataset import AMPTorchDataset, collate_amp
from amptorch.model_geometric import BPNN
from torch.utils.data import DataLoader
from amptorch.descriptor.Gaussian import Gaussian
from amptorch.descriptor.descriptor_calculator import DescriptorCalculator
from amptorch.dataset_backup import AtomsDataset
from amp.descriptor.gaussian import Gaussian as Gnew
from amp.utilities import hash_images as amp_hash
from ase.io.trajectory import Trajectory
from ase.io import read
from ase import Atoms
from ase.calculators.emt import EMT


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


Gs = {
    "default": {
        "G2": {
            "etas": np.logspace(np.log10(0.05), np.log10(5.0), num=4),
            "rs_s": [0] * 4,
        },
        "G4": {"etas": [0.005], "zetas": [1.0, 4.0], "gammas": [1.0, -1.0]},
        "cutoff": 6,
    },
}

elements = image.get_chemical_symbols()

descriptor = Gaussian(Gs=Gs, elements=elements)
desc_calc = DescriptorCalculator(
    images=images,
    descriptor=descriptor,
    automatic_calculation=False,
    calculate_descriptor_primes=True,
    sparse_prime=False,
    store_descriptors=True,
    training_data=False,
    parallel=False,
    cores=1,
)
desc_calc.prepare_descriptors()
new = desc_calc._get_calculated_descriptors()

# Old
Gs = {}
Gs["G2_etas"] = np.logspace(np.log10(0.05), np.log10(5.0), num=4)
Gs["G2_rs_s"] = [0] * 4
Gs["G4_etas"] = [0.005]
Gs["G4_zetas"] = [1.0, 4.0]
Gs["G4_gammas"] = [+1.0, -1]
Gs["cutoff"] = 6


forcetraining = True
dataset = AtomsDataset(
    images,
    Gnew,
    Gs,
    forcetraining=forcetraining,
    label="test",
    cores=4,
)

image_hashes = list(amp_hash(images).keys())
for i in range(len(images)):
    image_hash = image_hashes[i]
    with open("amp-data-fingerprints.ampdb/loose/" + image_hash, "rb") as f:
        amp_fp = pickle.load(f)
    os.system("rm amp-data-fingerprints.ampdb/loose/" + image_hash)

    with open("amp-data-fingerprint-primes.ampdb/loose/" + image_hash, "rb") as f:
        amp_prime = pickle.load(f)
    os.system("rm amp-data-fingerprint-primes.ampdb/loose/" + image_hash)
