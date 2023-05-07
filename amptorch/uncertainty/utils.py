import numpy as np
import torch
from sklearn.preprocessing import StandardScaler
from scipy import stats
from scipy.optimize import minimize


def calc_dist(train_X, calib_X, nearest_neighbors=10):
    """
    Returns the distances of calibration data to training data
    """
    nrows, ncols = train_X.shape
    from sklearn.neighbors import KDTree

    std_scaler = StandardScaler().fit(train_X)
    train_X_std_scaled = std_scaler.transform(train_X)
    calib_X_std_scaled = std_scaler.transform(calib_X)

    kdtree = KDTree(train_X_std_scaled)
    dist, ind = kdtree.query(calib_X_std_scaled, k=nearest_neighbors)

    return dist.mean(axis=1)


class ConformalPrediction:
    """
    Performs quantile regression on score functions to obtain the estimated qhat
        on calibration data and apply to test data during prediction.
    """

    def __init__(self, alpha):
        self.alpha = alpha

    def fit(self, residuals_calib, heurestic_uncertainty_calib) -> None:
        # score function
        scores = abs(residuals_calib / heurestic_uncertainty_calib)
        scores = np.array(scores)

        n = len(residuals_calib)
        qhat = torch.quantile(
            torch.from_numpy(scores), np.ceil((n + 1) * (1 - self.alpha)) / n
        )
        qhat_value = np.float64(qhat.numpy())
        self.qhat = qhat_value
        pass

    def predict(self, heurestic_uncertainty_test):
        cp_uncertainty_test = heurestic_uncertainty_test * self.qhat
        return cp_uncertainty_test, self.qhat


class NegativeLeastLikelihoodEstimator:
    """
    Performs negative least likelihood estimation with assumption that
        errors follow a Gaussian distribution to estimate the standard deviation
        based on latent distances.
    """

    def __init__(self):
        pass

    def fit(self, calib_y, calib_dist, initParams=[5e-8, 4]):
        self.calib_y = calib_y
        self.calib_dist = calib_dist
        results = minimize(self.gaussian, initParams, method="Nelder-Mead")
        self.s1 = results.x[0]
        self.s2 = results.x[1]
        pass

    def predict(self, test_dist):
        test_std = self.s1 + test_dist * self.s2
        return test_std

    def gaussian(self, params, calib_y, calib_dist):
        s1 = params[0]
        s2 = params[1]

        # Calculate negative log likelihood
        nll = -np.sum(
            stats.norm.logpdf(self.calib_y, loc=0, scale=s1 + self.calib_dist * s2)
        )

        return nll


def split_test_calib(full_test_X, full_test_y, per_calib, seed=0):
    """
    Uniformaly sample the test data at random to split a list for test data and a list for calibration.
    """
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


def prepare_latentNerror_from_trainer(trainer, images, get_latent=-2):
    """
    Return the latent representations and respective residuals for data list.
    """
    predictions = trainer.predict(images, get_latent=get_latent)
    feature = np.array(predictions["latent"])

    # errors for train, test, calib
    pred_energies = np.array(predictions["energy"])
    true_energies = np.array([image.get_potential_energy() for image in images])
    list_of_errors = pred_energies - true_energies

    return feature, list_of_errors


def prepare_featureNerror_from_trainer(trainer, images):
    """
    Return the latent representations and respective residuals for data list.
    """
    predictions = trainer.predict(images, get_descriptor=True)
    feature = np.array(predictions["descriptors"])

    # errors for train, test, calib
    pred_energies = np.array(predictions["energy"])
    true_energies = np.array([image.get_potential_energy() for image in images])
    list_of_errors = pred_energies - true_energies

    return feature, list_of_errors


def calc_uncertainty_metrics(list_epa, uncertainty):
    """
    Quantifies how the observed confidence level as sample coverage and the average widths of
        prediction sets.
    """
    overconfident_idx = np.argwhere(abs(list_epa) > uncertainty)

    prob_overconfidence = len(overconfident_idx) / len(uncertainty)
    avg_uncertainty = np.mean(uncertainty)

    return prob_overconfidence, avg_uncertainty
