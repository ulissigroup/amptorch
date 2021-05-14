import numpy as np
from ase import Atoms
from ase.calculators.emt import EMT

from amptorch.dataset import AtomsDataset
from NNSubsampling import subsampling

from pykdtree.kdtree import KDTree
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import Ridge
from tqdm import tqdm

import sys
import os
import pickle

def partition(data, refdata):
    kd_tree = KDTree(refdata,leafsize=6)
    distances, indices = kd_tree.query(data, k=1)
    
    indices, counts = np.unique(indices, return_counts=True)
    count_arr = np.zeros(len(refdata))
    
    for i, index in enumerate(indices):
        count_arr[index] = counts[i]

    max_distance = np.max(distances)
    #print("max distance")
    #print(max_distance)
    
    return count_arr, max_distance

def test_model(count_arr_train, target_train, count_arr_test, target_test, alpha):
    print("=========== start fit alpha {} =============".format(alpha))
    reg = Ridge(alpha=alpha).fit(count_arr_train, target_train)
    print("end fit")
    coef = reg.coef_

            
    prediction_train = reg.predict(count_arr_train)
    error_train = target_train - prediction_train

    print("training MAE: {}".format(np.mean(np.abs(error_train))))
    print("training MSE: {}".format(np.mean(np.square(error_train))))


    prediction_test = reg.predict(count_arr_test)
    error_test = target_test - prediction_test


    print("test MAE: {}".format(np.mean(np.abs(error_test))))
    print("test MSE: {}".format(np.mean(np.square(error_test))))



training_filename = sys.argv[1]
test_filename = sys.argv[2]
cutoff_sig = float(sys.argv[3])

training_images = pickle.load( open( training_filename, "rb" ) )
training_len = len(training_images)

test_images = pickle.load( open( test_filename, "rb" ) )
test_len = len(test_images)

images = training_images + test_images

os.chdir("./partition_models/")

sigmas = [0.02,0.12,0.24,0.36,0.5,0.69,0.92,1.2,1.52,2.0,2.66,3.5,5.0]
MCSHs = {
    "MCSHs": {
        "0": {"groups": [1], "sigmas": sigmas},
        "1": {"groups": [1], "sigmas": sigmas},
        "2": {"groups": [1, 2], "sigmas": sigmas},
        # "3": {"groups": [1, 2, 3], "sigmas": sigmas},
        # "4": {"groups": [1, 2, 3, 4], "sigmas": sigmas},
        # "5": {"groups": [1, 2, 3, 4, 5], "sigmas": sigmas},
        # "6": {"groups": [1, 2, 3, 4, 5, 6, 7], "sigmas": sigmas},
    },
    "atom_gaussians": {
                        "H": "../../valence_gaussians/H_pseudodensity_2.g",
                        "C": "../../valence_gaussians/C_pseudodensity_4.g",
                        "O": "../../valence_gaussians/O_pseudodensity_4.g",
                        "N": "../../valence_gaussians/N_pseudodensity_4.g",
                        "F": "../../valence_gaussians/F_pseudodensity_4.g",
    },
    "cutoff": 8,
}

elements = ["H","C","N","O","F"]
fp_elements = ["H","C","N","O","F"]
ref_elements = ["H","C","N","O","F"]

dataset = AtomsDataset(
    images=images,
    descriptor_setup=("gmp", MCSHs,{}, elements, fp_elements, ref_elements),
    forcetraining=False,
    save_fps=True,
    scaling={"type": "normalize", "range": (0, 1), "elementwise": False}
)

# for data in dataset:
#     print(data.fingerprint.numpy().shape)


overall_fp = np.vstack([data.fingerprint.numpy() for data in dataset[:training_len]])



ref_filename = "Overall_subsampled_{}.pickle".format(cutoff_sig)
try:
    refdata = pickle.load( open( ref_filename, "rb" ))

except:
    refdata = np.array(subsampling(overall_fp,list_desc = [],cutoff_sig=cutoff_sig,rate=0.2, method = "pykdtree",verbose = 2))
    pickle.dump( refdata, open( ref_filename, "wb" ) )
   
num_unique_env = len(refdata)

print("number of refdata: {}".format(len(refdata)))




model_setup = {"refdata_filename": ref_filename, \
                "PCA": False,\
                "PCA_n_components":15}

model = {"model_setup": model_setup}


if model_setup["PCA"]:
    pca = PCA(n_components=model_setup["PCA_n_components"])
    pca.fit(refdata)
    refdata_transformed= pca.transform(refdata)
    model["PCA_model"] = pca
else:
    refdata_transformed = refdata

model["refdata_transformed"] = refdata_transformed

count_arr = np.zeros((training_len, len(refdata)))
target = np.zeros(training_len)
max_distance = 0

#for i, system in enumerate(dataset[:training_len]):
for i, system in tqdm(
    enumerate(dataset[:training_len]),
    desc="partition image",
    total=len(dataset[:training_len]),
    unit=" images",
    ):
    # try:
    #print("start processing system {}".format(i))
    temp_feature_arr = system.fingerprint.numpy()
    if model_setup["PCA"]:
        feature_arr_transformed = model["PCA_model"].transform(temp_feature_arr)
    else:
        feature_arr_transformed = temp_feature_arr
    count_arr[i], temp_max_distance = partition(feature_arr_transformed, refdata_transformed)
    max_distance = max(max_distance, temp_max_distance)
    target[i] = system.energy
    # except:
    #     print("ERROR: system {} not processed".format(system))


count_arr2 = np.zeros((test_len, len(refdata)))
target2 = np.zeros(test_len)

for i, system in tqdm(
    enumerate(dataset[training_len:]),
    desc="partition image",
    total=len(dataset[training_len:]),
    unit=" images",
    ):
    # try:
    #print("start processing system {}".format(i))
    temp_feature_arr = system.fingerprint.numpy()
    if model_setup["PCA"]:
        feature_arr_transformed = model["PCA_model"].transform(temp_feature_arr)
    else:
        feature_arr_transformed = temp_feature_arr
    count_arr2[i], temp_max_distance = partition(feature_arr_transformed, refdata_transformed)
    target2[i] = system.energy
    # except:
    #     print("ERROR: system {} not processed".format(system))


model["max_distance"] = max_distance

for alpha in [1.0, 3.0, 10.0, 30.0, 100.0, 300]:
    test_model(count_arr, target, count_arr2, target2, alpha)



