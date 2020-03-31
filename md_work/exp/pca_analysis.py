import numpy as np
from amptorch.data_preprocess import AtomsDataset
from amptorch.gaussian import SNN_Gaussian
from ase.io import read
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

label = "pca_analysis"
# images = read("../../datasets/water/water.extxyz", ":")
images = read("../../datasets/COCu_ber_100ps_300K.traj", ":2000")

# define symmetry functions to be used
Gs = {}
Gs["G2_etas"] = np.logspace(np.log10(0.05), np.log10(5.0), num=4)
Gs["G2_rs_s"] = [0] * 4
Gs["G4_etas"] = [0.005]
Gs["G4_zetas"] = [1.0]
Gs["G4_gammas"] = [+1.0, -1]
Gs["cutoff"] = 6.5

forcetraining = True
training_data = AtomsDataset(
    images,
    SNN_Gaussian,
    Gs,
    forcetraining=forcetraining,
    label=label,
    cores=4,
    lj_data=None,
)
fingerprints = training_data.fingerprint_dataset
primes = training_data.sparse_fprimes
energies = np.array(training_data.energy_dataset).reshape(-1,)
forces = training_data.forces_dataset
elements = training_data.elements
fp_len = training_data.fp_length
num_atoms = len(fingerprints[0])

fp_array = np.empty((len(fingerprints), num_atoms*fp_len))
element_fp_array = {}
for idx, image in enumerate(fingerprints):
    i = 0
    for atom in image:
        fp_array[idx, i:fp_len+i] = atom[1]
        i += fp_len

pca = PCA(n_components=2)
projected = pca.fit_transform(fp_array)
print(pca.explained_variance_ratio_)

plt.scatter(projected[:, 0], projected[:, 1], c=energies)
plt.xlabel('PC 1', fontsize=18)
plt.ylabel('PC 2', fontsize=18)
plt.colorbar()
plt.show()
