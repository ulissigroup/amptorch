import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator as sp
import copy
import warnings

import os
import copy
import pickle
import timeit
# import random
# numpy.random supports sample with probability
from numpy import random
from sklearn import preprocessing
from sklearn.decomposition import PCA
from ase.calculators.singlepoint import SinglePointCalculator as sp
from amptorch.utils import hash_images, get_hash
from amptorch.data_preprocess import AtomsDataset
from amptorch.gaussian import SNN_Gaussian

"""
All query methods should have the following arguments
    Paramaters
    ----------
    images: List. Current training images.

    sample_candidates: List. Potential candidates for query as collected
    from the generating function.

    samples_to_retrain: Integer. Number of samples to be randomly selected
    for query and added to the training set.

    parent_calc: Object. Parent calculator used to query energy and force
    calculations.
"""


def random_query(_, sample_candidates, samples_to_retrain, parent_calc,
                 params=None, method=None):
    """
    Randomly selects data points from a list of potential candidates to
    query and add to the existing training set.
    """
    if len(sample_candidates) <= samples_to_retrain:
        warnings.warn(
            "# of samples exceeds # of available candidates! Defaulting to all available candidates",
            stacklevel=2,
        )
        samples_to_retrain = len(sample_candidates) - 1
    query_idx = random.sample(range(1, len(sample_candidates)), samples_to_retrain)
    images_to_query = [sample_candidates[idx] for idx in query_idx]
    queried_images = compute_query(images_to_query, parent_calc)
    return queried_images


def max_uncertainty(_, sample_candidates, samples_to_retrain, parent_calc):
    """Selects points with the largest uncertainty"""
    if len(sample_candidates) < samples_to_retrain:
        warnings.warn(
            "# of samples exceeds # of available candidates! Defaulting to all available candidates",
            stacklevel=2,
        )
        samples_to_retrain = len(sample_candidates) - 1
    uncertainty = np.array(
        [atoms.info["uncertainty"][0] for atoms in sample_candidates]
    )
    query_idx = np.argpartition(uncertainty, -1 * samples_to_retrain)[
        -1 * samples_to_retrain :
    ]
    images_to_query = [sample_candidates[idx] for idx in query_idx]
    queried_images = compute_query(images_to_query, parent_calc)
    return queried_images


def compute_query(images_to_calculate, parent_calc):
    queried_images = []
    for image in images_to_calculate:
        image.set_calculator(copy.copy(parent_calc))
        sample_energy = image.get_potential_energy(apply_constraint=False)
        sample_forces = image.get_forces(apply_constraint=False)
        image.set_calculator(
            sp(atoms=image, energy=sample_energy, forces=sample_forces)
        )
        queried_images.append(image)
    return queried_images


def termination_criteria(termination_args, method="iter"):
    """Criteria for AL termination
    Parameters
    ----------
    method: str
        Method for termination of active learning loop.

    'iter': Terminates after specified number of iterations.
        args: (current iteration, # of iterations)
    """
    terminate = False

    if method == "iter":
        current_i = termination_args["current_i"]
        total_i = termination_args["total_i"]
        if current_i > total_i:
            terminate = True
    if method == "final":
        calc = copy.copy(termination_args["calc"])
        final_image = termination_args["images"][-1]
        e_tol = termination_args["energy_tol"]
        f_tol = termination_args["force_tol"]

        ml_energy = final_image.get_potential_energy(apply_constraint=False)
        ml_forces = final_image.get_forces(apply_constraint=False).flatten()

        parent_energy = calc.get_potential_energy(final_image)
        parent_forces = calc.get_forces(final_image).flatten()

        e_terminate = False
        f_terminate = False

        if np.abs(ml_energy-parent_energy)/len(final_image) <= e_tol:
            e_terminate = True
        if np.sum(np.abs(ml_forces-parent_forces))/(3*len(final_image)) <= f_tol:
            f_terminate = True

        terminate = e_terminate and f_terminate
    return terminate


def random_query_with_prob(images, sample_candidates, samples_to_retrain,
                           parent_calc, params, method='distance'):
    """
    Sample images with probabilities being their uncertainties
    """
    sample_candidates = sample_candidates[1:]
    start_time = timeit.default_timer()
    # generate fingerprints for candidates and parent images
    fm_images, fm_candidates = get_fingerprints(images, sample_candidates, params)
    if method == 'distance':
        uncertainties = [compute_uncertainty(fm_images, _) for _ in fm_candidates]
    elif method == 'density':
        uncertainties = get_debu(fm_images, fm_candidates)
    else:
        raise TypeError('Please use the right type of unceratiny method, currently'
            ' only distance or density based methods are supported.')
    # print(min(uncertainties))
    elapsed = timeit.default_timer() - start_time
    print(f"Uncertainties evaluated! Time cost: {elapsed : .2f} seconds")

    prob = np.array(uncertainties)/np.array(uncertainties).sum()
    query_idx = random.choice(
        range(len(sample_candidates)), samples_to_retrain, p=prob
    )
    images_to_query = [sample_candidates[idx] for idx in query_idx]
    queried_images = compute_query(images_to_query, parent_calc)
    return queried_images


def get_fingerprints(images, sample_candidates, params):
    """
    Get the fingerprints of the candidate images    
    """
    AtomsDataset(
                sample_candidates,
                SNN_Gaussian,
                Gs=params["Gs"],
                forcetraining=False,
                label='amp-data',
                cores=params["cores"],
    )
    fp_label = 'amp-data'
    fps_images = []
    for _ in images:
        fp_path = os.path.join(fp_label + "-fingerprints.ampdb", "loose",
                               get_hash(_, Gs=params["Gs"]))
        with open(fp_path, 'rb') as pf:
            fp = pickle.load(pf)
            fps_images.append(np.array([__[-1] for __ in fp ]).reshape(-1))
    fps_candidates = []
    for _ in sample_candidates:
        fp_path = os.path.join(fp_label + "-fingerprints.ampdb", "loose",
                       get_hash(_, Gs=params["Gs"]))
        with open(fp_path, 'rb') as pf:
            fp = pickle.load(pf)
            fps_candidates.append(np.array([__[-1] for __ in fp ]).reshape(-1))
    return fps_images, fps_candidates


def compute_uncertainty(fps_train, fps_target, distance_method='mean',
                   scaler='standard', weight=1., pca=True, 
                   n_components=10, order=2,  scale=0.1):
    """Compute the uncertainty of the target image based off distances
    between a point and a set of points (training set).

    Parameters
    ------------
    fps_train : list or numpy array
    vectors as constructed from the training set.

    fps_target: list or numpy array
    vector(s) as constructed from the target image. Note that at one time
    only the confidence of one image can be evaluated.

    scaler: str (default: "standard")
    algorithm used to scale the raw data

    weight: float, list or array (default: 1.)
    Weights for features.

    pca: boolean (default: True)
    For high dimensional data and highly correlated features, a pca
    transformation is recommended before computing the distances.

    n_components: int (default: 10)
    Number of principal components kept in PCA transformation

    order: int (default: 2)
    order of the norm for calculating distances
    
    scale: float (default: 0.1)
    scale factor for distances evaluation

    distance_method: str (default: 'mean')
    method for using the distances to compute the uncertainty

    Return
    --------
    Uncertainty level of the target image
    """

    fps_train, fps_target = np.array(fps_train), np.array(fps_target)
    size_train, size_image = fps_train.shape, len(fps_target)
    assert size_image == size_train[1]
    v = fps_train
    if fps_train.shape[0] < n_components:
        scaler, pca = None, False
    # [fps_target] convert 1D array into 2D array,
    # Required by the `tranform` method
    if scaler == 'maxabs':
        maxabs_scaler = preprocessing.MaxAbsScaler()
        v_scaled = maxabs_scaler.fit_transform(v)
        fps_target = maxabs_scaler.transform([fps_target])
        scaler = maxabs_scaler
    elif scaler == 'minmax':
        minmax_scaler = preprocessing.MinMaxScaler()
        v_scaled = minmax_scaler.fit_transform(v)
        fps_target = minmax_scaler.transform([fps_target])
        scaler = minmax_scaler
    elif scaler == 'standard':
        standard_scaler = preprocessing.StandardScaler().fit(v)
        v_scaled = standard_scaler.transform(v)
        fps_target = standard_scaler.transform([fps_target])
        scaler = standard_scaler
    else:
        v_scaled, fps_target = fps_train, fps_target
    if pca:
        assert isinstance(n_components, int)
        pca = PCA(n_components=n_components)
        v_scaled = pca.fit_transform(v_scaled)
        fps_target = pca.transform(fps_target)[0]
    v_scaled = np.vstack([v_scaled, fps_target])
    assert isinstance(order, int)
    v_dist = ((weight*(np.abs(v_scaled - v_scaled[-1])**order))
              .sum(axis=1))**(1/order)
    v_dist = v_dist[:-1]
    uncertainty = uncertainty_func(v_dist, scale=scale, method=distance_method)
    return uncertainty


def uncertainty_func(x, scale=1., method='mean'):
    """
    function to calculate uncertainty. It takes as input the (weighted)
    vector distances between a target image and the training images.
    The default is the average of the confidence. The other two methods are
    
    - the max of the confidence

    - the average of the first n maximum confidence, default n is 20% of
    training images
    """
    assert method in ['mean', 'min', 'min_n']
    x = np.exp(-1*scale*x)
    if method == 'mean':
        return np.mean(1/(1+x))
    elif method == 'min':
        return np.max(1/(1+x))
    else:
        n = int(0.2*len(x))
        x = np.sort(x)[:n]
        return np.mean(1/(1+x))


def get_debu(fps_train, fps_target, nbins=50, with_pca=True, n_components=10,
                normalized=True, scale=0.1):
        """histogram density-based uncertainty (DEBU) method

        Parameters
        ------------
        fps_train : list or numpy array
        vectors as constructed from the training set.

        fps_target: list or numpy array
        vector(s) as constructed from the target image. Note that at one time
        only the confidence of one image can be evaluated.

        nbins : int
        number of bins for density estimation

        with_pca: bool
        If the feature vectors are preprocessed with principal component
        analysis

        n_components: int
        Number of component kept for PCA analysis

        scale: float (default: 0.1)
        scale factor for density evaluation

        Return
        --------
        The weighted proportions of features in the interpolation region
        for all target feature vectors.

        TODO: add kernel density estimation
        """
        from sklearn.preprocessing import MinMaxScaler
        nbins = int(nbins)
        fpv = np.array(fps_train)
        fpv_target = np.array(fps_target)
        # apply `MinMaxScaler` on the feature vectors of training images
        scaler = MinMaxScaler()
        fpv = scaler.fit_transform(fpv)
        fpv_target = scaler.transform(fpv_target)
        if fpv.shape[0] < n_components: with_pca, normalized = False, False
        if with_pca:
            pca = PCA(n_components=n_components)
            fpv = pca.fit_transform(fpv)
            fpv_target = pca.transform(fpv_target)
            variances = pca.explained_variance_ratio_
        n_cols = fpv.shape[1]
        # Compute histogram densities and bin_edges of training images
        hist_density, bin_edges = [], []
        for k in range(n_cols):
            hist, edges = np.histogram(fpv[:, k], bins=nbins, density=True)
            hist_density.append(hist)
            bin_edges.append(edges)
        debu = []
        for i in range(fpv_target.shape[0]):
            density_weights = np.zeros(n_cols)
            target = fpv_target[i,:]
            for j in range(n_cols):
                hist_per = hist_density[j]
                edges_per = bin_edges[j][:-1]
                pos = (edges_per - target[j] < 0).sum()
                if with_pca:
                    density_weights[j] = hist_per[pos-1]*variances[j]
                else:
                    density_weights[j] = hist_per[pos-1]
            if normalized:
                hist_density = np.array(hist_density)
                debu.append(sum(density_weights)/
                    sum(hist_density.max(axis=1) * variances))
            else:
                debu.append((density_weights).mean())
        return 1/(1+np.exp(-scale*np.array(debu)))