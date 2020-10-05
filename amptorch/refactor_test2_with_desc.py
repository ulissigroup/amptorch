import numpy as np
from amp.descriptor.gaussian import Gaussian as Gnew
from amp.utilities import hash_images as amp_hash
from ase import Atoms
from ase.calculators.emt import EMT

from amptorch.dataset import AtomsDataset as NewDataset
from amptorch.dataset import data_collater
from amptorch.dataset_backup import AtomsDataset as OldDataset
from amptorch.dataset_backup import collate_amp as old_collater
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
        "G2": {"etas": np.logspace(np.log10(0.05), np.log10(5.0), num=4), "rs_s": [0],},
        "G4": {"etas": [0.005], "zetas": [1.0, 4.0], "gammas": [1.0, -1.0]},
        "cutoff": 6,
    },
}

elements = ["Cu", "C", "O"]

descriptor = Gaussian(Gs=Gs_new, elements=elements)
newdataset = NewDataset(images, descriptor)
batch_new = data_collater([newdataset[i] for i in range(len(images))])[0]
new = batch_new.fprimes.coalesce().values()

# Old
Gs = {}
Gs["G2_etas"] = np.logspace(np.log10(0.05), np.log10(5.0), num=4)
Gs["G2_rs_s"] = [0] * 4
Gs["G4_etas"] = [0.005]
Gs["G4_zetas"] = [1.0, 4.0]
Gs["G4_gammas"] = [+1.0, -1]
Gs["cutoff"] = 6

forcetraining = True
olddataset = OldDataset(
    images, Gnew, Gs, forcetraining=forcetraining, label="test", cores=4,
)
batch_old = old_collater([olddataset[i] for i in range(len(images))])[0]
old = batch_old.fprimes.coalesce().values()
