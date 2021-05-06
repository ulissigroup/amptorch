import numpy as np
from ase import Atoms
from ase.calculators.emt import EMT

from amptorch.dataset import AtomsDataset
from NNSubsampling import subsampling

from pykdtree.kdtree import KDTree
from sklearn.linear_model import LinearRegression

def partition(data, refdata):
    kd_tree = KDTree(refdata,leafsize=6)
    distances, indices = kd_tree.query(data, k=1)
    
    indices, counts = np.unique(indices, return_counts=True)
    count_arr = np.zeros(len(refdata))
    
    for i, index in enumerate(indices):
        count_arr[index] = counts[i]

    max_distance = np.max(distances)
    print("max distance")
    print(max_distance)
    
    return count_arr, max_distance

distances = np.linspace(2, 5, 10)
images = []
for dist in distances:
    image = Atoms(
        "CuCO",
        [
            (-dist * np.sin(0.65), dist * np.cos(0.65), 0),
            (0, 0, 0),
            (dist * np.sin(0.65), dist * np.cos(0.65), 0),
        ],
    )
    image.set_cell([10, 10, 10])
    image.wrap(pbc=True)
    image.set_calculator(EMT())
    images.append(image)


sigmas = [0.02, 0.2, 0.4, 0.69, 1.1, 1.66, 2.66, 4.4]
MCSHs = {
    "MCSHs": {
        "0": {"groups": [1], "sigmas": sigmas},
        "1": {"groups": [1], "sigmas": sigmas},
        "2": {"groups": [1, 2], "sigmas": sigmas},
        "3": {"groups": [1, 2, 3], "sigmas": sigmas},
        # "4": {"groups": [1, 2, 3, 4], "sigmas": sigmas},
        # "5": {"groups": [1, 2, 3, 4, 5], "sigmas": sigmas},
        # "6": {"groups": [1, 2, 3, 4, 5, 6, 7], "sigmas": sigmas},
    },
    "atom_gaussians": {
        "C": "./valence_gaussians/C_pseudodensity_4.g",
        "O": "./valence_gaussians/O_pseudodensity_4.g",
        "Cu": "./valence_gaussians/Cu_pseudodensity_4.g",
    },
    "cutoff": 8,
}

elements = ["Cu", "C", "O"]
fp_elements = ["Cu", "C", "O"]
ref_elements = ["Cu", "C", "O"]
dataset = AtomsDataset(
    images=images,
    descriptor_setup=("gmp", MCSHs,{}, elements, fp_elements, ref_elements),
    forcetraining=False,
    save_fps=True,
    scaling={"type": "normalize", "range": (0, 1), "elementwise": False}
)

# for data in dataset:
#     print(data.fingerprint.numpy().shape)


overall_fp = np.vstack([data.fingerprint.numpy() for data in dataset])


cutoff_sig = 0.1
ref_filename = "Overall_subsampled_{}.pickle".format(cutoff_sig)
ref_data = np.array(subsampling(overall_fp,list_desc = [],cutoff_sig=cutoff_sig,rate=0.2, method = "pykdtree",verbose = 2))

num_unique_env = len(ref_data)

pickle.dump( subsampled, open( ref_filename, "wb" ) )





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

count_arr = np.zeros((len(systems), len(refdata)))
target = np.zeros(len(systems))
max_distance = 0

for i, system in enumerate(dataset):
    # try:
    print("start processing system {}".format(system))
    temp_feature_arr = data.fingerprint.numpy()
    if model_setup["PCA"]:
        feature_arr_transformed = model["PCA_model"].transform(temp_feature_arr)
    else:
        feature_arr_transformed = temp_feature_arr
    count_arr[i], temp_max_distance = partition(feature_arr_transformed, refdata_transformed)
    max_distance = max(max_distance, temp_max_distance)
    target[i] = system.energy
    # except:
    #     print("ERROR: system {} not processed".format(system))

model["max_distance"] = max_distance



reg = LinearRegression().fit(count_arr, target)
coef = reg.coef_

model["regression_model"] = reg

        
prediction = reg.predict(count_arr)
error = target - prediction

print("training MAE: {}".format(np.mean(np.abs(error))))
print("training MSE: {}".format(np.mean(np.square(error))))

pickle.dump( model, open( "cutoff{}_model.pickle".format(cutoff_sig), "wb" ) )
