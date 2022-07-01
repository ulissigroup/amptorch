import numpy as np
from numpy.random import default_rng

from sklearn.neighbors import KDTree
from sklearn.preprocessing import StandardScaler

from scipy.stats import gaussian_kde

def calc_dist(train_X, calib_X, nearest_neighbors=10, metric="minkowski"):
    std_scaler = StandardScaler().fit(train_X)
    train_X_std_scaled = std_scaler.transform(train_X)
    calib_X_std_scaled = std_scaler.transform(calib_X)
    
    kdtree = KDTree(train_X_std_scaled, metric=metric)
    dist, ind = kdtree.query(calib_X_std_scaled, k=nearest_neighbors)
    
    return dist.mean(axis=1)

def predict_data(pred_energies, test_images):    
    true_energies = np.array([image.get_potential_energy() for image in test_images])

    list_of_error_per_atom = []

    for i, image in enumerate(test_images):
        num_atoms = len(image.get_atomic_numbers())
        total_energy_pred = pred_energies[i]
        total_energy_true = true_energies[i]

        error = pred_energies[i] - true_energies[i]
        list_of_error.append(error)

    return list_of_error

class ConformalPrediction():
    def __init__(self, residuals_calib, heurestic_uncertainty_calib, alpha) -> None:
        # score function
        scores = abs(residuals_calib / heurestic_uncertainty_calib)
        scores = np.array(scores)

        n=len(residuals_calib)
        qhat = torch.quantile(torch.from_numpy(scores), np.ceil((n + 1) * (1 - alpha)) / n)
        qhat_value = np.float64(qhat.numpy())
        self.qhat = qhat_value
        pass

    def predict(self, heurestic_uncertainty_test):
        return heurestic_uncertainty_test * self.qhat, self.qhat


def split_test_calib(full_test_X, full_test_y, per_calib, seed=0):
    # sub-select calib and test
    np.random.seed(seed)

    # num_total=5000
    num_total = len(full_test_X)

    num_calib = round(per_calib * num_total)

    rand_idx = np.random.choice(num_total, size=num_total, replace=False)
    calib_idx = rand_idx[:num_calib]
    test_idx = rand_idx[num_calib::]

    calib_X = [full_test_X[_] for _ in calib_idx]
    test_X = [full_test_X[_] for _ in test_idx]

    calib_y = [full_test_y[_] for _ in calib_idx]
    test_y = [full_test_y[_] for _ in test_idx]

    return test_X, test_y, calib_X, calib_y


def prepare_latentNerror_from_trainer(trainer, images, image_type, get_latent=-2):
    predictions = trainer.predict(images, get_latent=get_latent)
    feature = np.array(predictions["latent"])

    # errors for train, test, calib
    if image_type in ["test", "calib", "train"]:
        y_bar = np.array(predictions["energy"])
        list_error = predict_data(y_bar, images)
        list_error = np.array(list_error)
    else:
        raise NotImplementedError  

    return feature, list_error

def calc_uncertainty_metrics(list_epa, uncertainty):
    overconfident_idx = np.argwhere(abs(list_epa) > uncertainty)

    prob_overconfidence = len(overconfident_idx) / len(uncertainty)
    avg_uncertainty = np.mean(uncertainty)
    
    return prob_overconfidence, avg_uncertainty
