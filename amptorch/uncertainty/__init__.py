import numpy as np
from .utils import (
    calc_dist,
    ConformalPrediction,
    NegativeLeastLikelihoodEstimator,
    split_test_calib,
    prepare_latentNerror_from_trainer,
    prepare_featureNerror_from_trainer,
    calc_uncertainty_metrics,
)


class EnsembleSDEstimator:
    def __init__(self):
        pass

    def fit_predict(self, list_trainers, train_list, test_list):
        # compute residuals for train and test data
        mat_residuals_test = np.empty((len(test_list), len(list_trainers)))

        for i, trainer in enumerate(list_trainers):
            _, residuals_test = prepare_featureNerror_from_trainer(trainer, test_list)
            mat_residuals_test[:, i] = residuals_test

        test_avg = np.mean(mat_residuals_test, axis=1)
        test_std = np.std(mat_residuals_test, axis=1)
        res = {
            "residuals": test_avg,
            "uncertainty": test_std,
        }
        return res


class ConformalPredictionLatentSpace:
    def __init__(
        self, alpha=0.1, per_calib=0.1, num_nearest_neighbors=10, seed=1, get_latent=-2
    ):
        # set parameters
        self.alpha = alpha
        self.per_calib = per_calib
        self.num_nearest_neighbors = num_nearest_neighbors
        self.seed = seed
        self.get_latent = get_latent

    def fit_predict(self, trainer, train_list, full_test_list):
        # compute the latent representations and residuals for train and test data
        train_X, train_y = prepare_latentNerror_from_trainer(
            trainer, train_list, get_latent=self.get_latent
        )
        full_test_X, full_test_y = prepare_latentNerror_from_trainer(
            trainer, full_test_list, get_latent=self.get_latent
        )

        # predict test and calib
        test_X, test_y, calib_X, calib_y = split_test_calib(
            full_test_X, full_test_y, self.per_calib, seed=self.seed
        )
        train_X = np.array(train_X)
        train_y = np.array(train_y)
        test_X = np.array(test_X)
        test_y = np.array(test_y)
        calib_X = np.array(calib_X)
        calib_y = np.array(calib_y)

        # calculating the distance metric for test+calib
        calib_dist = calc_dist(
            train_X, calib_X, nearest_neighbors=self.num_nearest_neighbors
        )
        test_dist = calc_dist(
            train_X, test_X, nearest_neighbors=self.num_nearest_neighbors
        )

        # conformal prediction for expected confidence level as (1-alpha)100%
        model_cp = ConformalPrediction(alpha=self.alpha)
        model_cp.fit(calib_y, calib_dist)
        test_uncertainty, qhat = model_cp.predict(test_dist)
        res = {
            "residuals": test_y,
            "uncertainty": test_uncertainty,
            "alpha": self.alpha,
        }
        return res


class NegativeLeastLikelihoodLatentSpace:
    def __init__(self, per_calib=0.1, num_nearest_neighbors=10, seed=1, get_latent=-2):
        # set parameters
        self.per_calib = per_calib
        self.num_nearest_neighbors = num_nearest_neighbors
        self.seed = seed
        self.get_latent = get_latent

    def fit_predict(self, trainer, train_list, full_test_list):
        # compute the latent representations and residuals for train and test data
        train_X, train_y = prepare_latentNerror_from_trainer(
            trainer, train_list, get_latent=self.get_latent
        )
        full_test_X, full_test_y = prepare_latentNerror_from_trainer(
            trainer, full_test_list, get_latent=self.get_latent
        )

        # predict test and calib
        test_X, test_y, calib_X, calib_y = split_test_calib(
            full_test_X, full_test_y, self.per_calib, seed=self.seed
        )
        train_X = np.array(train_X)
        train_y = np.array(train_y)
        test_X = np.array(test_X)
        test_y = np.array(test_y)
        calib_X = np.array(calib_X)
        calib_y = np.array(calib_y)

        # calculating the distance metric for test+calib
        calib_dist = calc_dist(
            train_X, calib_X, nearest_neighbors=self.num_nearest_neighbors
        )
        test_dist = calc_dist(
            train_X, test_X, nearest_neighbors=self.num_nearest_neighbors
        )

        # nagetive least likelihood estimation with assumption of Gaussian distribution
        model_nll = NegativeLeastLikelihoodEstimator()
        model_nll.fit(calib_y, calib_dist)
        test_uncertainty = model_nll.predict(test_dist)
        res = {
            "residuals": test_y,
            "uncertainty": test_uncertainty,
        }
        return res


class ConformalPredictionFeatureSpace:
    def __init__(self, alpha=0.1, per_calib=0.1, num_nearest_neighbors=10, seed=1):
        # set parameters
        self.alpha = alpha
        self.per_calib = per_calib

    def fit_predict(self, trainer, train_list, full_test_list):
        # compute the feature representations and residuals for train and test data
        train_X, train_y = prepare_featureNerror_from_trainer(trainer, train_list)
        full_test_X, full_test_y = prepare_featureNerror_from_trainer(
            trainer, full_test_list
        )

        # predict test and calib
        test_X, test_y, calib_X, calib_y = split_test_calib(
            full_test_X, full_test_y, self.per_calib, seed=self.seed
        )
        train_X = np.array(train_X)
        train_y = np.array(train_y)
        test_X = np.array(test_X)
        test_y = np.array(test_y)
        calib_X = np.array(calib_X)
        calib_y = np.array(calib_y)

        # calculating the distance metric for test+calib
        calib_dist = calc_dist(
            train_X, calib_X, nearest_neighbors=self.num_nearest_neighbors
        )
        test_dist = calc_dist(
            train_X, test_X, nearest_neighbors=self.num_nearest_neighbors
        )

        # conformal prediction for expected confidence level as (1-alpha)100%
        model_cp = ConformalPrediction(alpha=self.alpha)
        model_cp.fit(calib_y, calib_dist)
        test_uncertainty, qhat = model_cp.predict(test_dist)
        res = {
            "residuals": test_y,
            "uncertainty": test_uncertainty,
            "alpha": self.alpha,
        }
        return res


class NegativeLeastLikelihoodFeatureSpace:
    def __init__(self, per_calib=0.1, num_nearest_neighbors=10, seed=1, get_latent=-2):
        # set parameters
        self.per_calib = per_calib
        self.num_nearest_neighbors = num_nearest_neighbors
        self.seed = seed
        self.get_latent = get_latent

    def fit_predict(self, trainer, train_list, full_test_list):
        # compute the feature representations and residuals for train and test data
        train_X, train_y = prepare_featureNerror_from_trainer(trainer, train_list)
        full_test_X, full_test_y = prepare_featureNerror_from_trainer(
            trainer, full_test_list
        )

        # predict test and calib
        test_X, test_y, calib_X, calib_y = split_test_calib(
            full_test_X, full_test_y, self.per_calib, seed=self.seed
        )
        train_X = np.array(train_X)
        train_y = np.array(train_y)
        test_X = np.array(test_X)
        test_y = np.array(test_y)
        calib_X = np.array(calib_X)
        calib_y = np.array(calib_y)

        # calculating the distance metric for test+calib
        calib_dist = calc_dist(
            train_X, calib_X, nearest_neighbors=self.num_nearest_neighbors
        )
        test_dist = calc_dist(
            train_X, test_X, nearest_neighbors=self.num_nearest_neighbors
        )

        # nagetive least likelihood estimation with assumption of Gaussian distribution
        model_nll = NegativeLeastLikelihoodEstimator()
        model_nll.fit(calib_y, calib_dist)
        test_uncertainty = model_nll.predict(test_dist)
        res = {
            "residuals": test_y,
            "uncertainty": test_uncertainty,
        }
        return res
