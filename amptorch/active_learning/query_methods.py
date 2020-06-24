import os
import copy
import pickle
import warnings
import numpy as np   # numpy.random supports sample with probability
from numpy import random
from sklearn import preprocessing
from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.neighbors import KernelDensity
from ase.calculators.singlepoint import SinglePointCalculator as sp
from amptorch.utils import get_hash
from amptorch.gaussian import SNN_Gaussian
from amptorch.data_preprocess import AtomsDataset

"""
All query methods should have the following arguments
    Paramaters
    ----------
    images: List. Current training images. It is also '_' if not explicitly
    indicated.

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
            '# of samples exceeds # of available candidates! '
            'Defaulting to all available candidates', stacklevel=2,
        )
        samples_to_retrain = len(sample_candidates) - 1
    query_idx = random.choice(range(1, len(sample_candidates)),
                              samples_to_retrain)
    images_to_query = [sample_candidates[idx] for idx in query_idx]
    queried_images = compute_query(images_to_query, parent_calc)
    return queried_images


def max_uncertainty(_, sample_candidates, samples_to_retrain, parent_calc):
    """Selects points with the largest uncertainty"""
    if len(sample_candidates) < samples_to_retrain:
        warnings.warn(
            '# of samples exceeds # of available candidates! '
            'Defaulting to all available candidates',
            stacklevel=2,
        )
        samples_to_retrain = len(sample_candidates) - 1
    uncertainty = np.array(
        [atoms.info["uncertainty"][0] for atoms in sample_candidates]
    )
    query_idx = np.argpartition(uncertainty, -1 * samples_to_retrain)[
        -1 * samples_to_retrain:
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

        if np.abs(ml_energy - parent_energy) / len(final_image) <= e_tol:
            e_terminate = True
        if np.sum(np.abs(ml_forces - parent_forces)) / (3 * len(final_image)) <= f_tol:
            f_terminate = True

        terminate = e_terminate and f_terminate
    return terminate


def random_query_with_prob(images, sample_candidates, samples_to_retrain,
                           parent_calc, params, method='distance'):
    """
    Sample images with probabilities being their relative uncertainties.
    """
    sample_candidates = sample_candidates[1:]
    fmaxs = [np.sqrt((np.array(_.get_forces())**2).sum(axis=1)).max()
             for _ in sample_candidates]
    sample_candidates = [candidate for _, candidate in
                         enumerate(sample_candidates)
                         if fmaxs[_] < 5.]  # not trust high fmax candidates
    # generate fingerprints for candidates and parent images
    fm_images, fm_candidates = get_fingerprints(images, sample_candidates,
                                                params)
    if method == 'distance':
        uncertainties = get_dibu(fm_images, fm_candidates)
    elif method == 'density':
        uncertainties = get_debu(fm_images, fm_candidates)
    else:
        raise TypeError('Please use available unceratiny methods, currently'
                        ' only distance/density-based methods are supported.')
    print(uncertainties)
    prob = np.array(uncertainties) / np.array(uncertainties).sum()
    query_idx = random.choice(
        range(len(sample_candidates)), samples_to_retrain, p=prob
    )
    clean_unselected_candidates(sample_candidates, query_idx, params)
    images_to_query = [sample_candidates[idx] for idx in query_idx]
    queried_images = compute_query(images_to_query, parent_calc)
    return queried_images


def get_fingerprints(images, sample_candidates, params):
    """
    Get the fingerprints of the candidate images using simplenn.
    """
    fp_label = 'amp-data'
    AtomsDataset(
        sample_candidates,
        SNN_Gaussian,
        Gs=params["Gs"],
        forcetraining=False,
        label=fp_label,
        cores=params["cores"],
    )
    fps_images = []
    for _ in images:
        fp_path = os.path.join(fp_label + "-fingerprints.ampdb", "loose",
                               get_hash(_, Gs=params["Gs"]))
        with open(fp_path, 'rb') as pf:
            fp = pickle.load(pf)
            fps_images.append(np.array([__[-1] for __ in fp]).reshape(-1))
    fps_candidates = []
    for _ in sample_candidates:
        fp_path = os.path.join(fp_label + "-fingerprints.ampdb", "loose",
                               get_hash(_, Gs=params["Gs"]))
        with open(fp_path, 'rb') as pf:
            fp = pickle.load(pf)
            fps_candidates.append(np.array([__[-1] for __ in fp]).reshape(-1))
    return fps_images, fps_candidates


def clean_unselected_candidates(sample_candidates, query_idx, params,
                                label='amp-data'):
    """
    Clean the fingerprints for unselected candidates. No need to clean
    fingerprint primes because only fingerprints are used for uncertainty
    evaluation.
    """
    unselected_images = [sample_candidates[_]
                         for _ in range(len(sample_candidates))
                         if _ not in query_idx]
    for _ in unselected_images:
        fp_path = os.path.join(label + "-fingerprints.ampdb", "loose",
                               get_hash(_, Gs=params["Gs"]))
        fp_prime_path = os.path.join(label + "-fingerprint-primes.ampdb", "loose",
                                     get_hash(_, Gs=params["Gs"]))
        neighbor_path = os.path.join(label + "-neighborlists.ampdb", "loose",
                                     get_hash(_, Gs=params["Gs"]))
        os.remove(fp_path)
        os.remove(fp_prime_path)
        os.remove(neighbor_path)


def get_dibu(fps_train, fps_target, weighted_method='mean',
             scaler='standard', weight=1., pca=True,
             n_components=10, order=2, scale=1.):
    """
    Compute the uncertainty of the target image based off distances
    between a point and a set of points (training set). DIBU is short
    for distance-based uncertainty metrics.

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

    scale: float (default: 1.)
        scale factor for distances evaluation

    weighted_method: str (default: 'mean')
        how we treat the distances in terms of uncertainty calculation

    Return
    --------
    Uncertainty level of the target image
    """
    fps_train, fps_target = np.array(fps_train), np.array(fps_target)
    assert fps_train.shape[1] == fps_target.shape[1]
    if fps_train.shape[0] < n_components:
        pca = False
    # [fps_target] convert 1D array into 2D array,
    # Required by the `tranform` method
    if scaler == 'standard':
        SCALER = preprocessing.StandardScaler()
    elif scaler == 'minmax':
        SCALER = preprocessing.MinMaxScaler()
    elif scaler == 'maxabs':
        SCALER = preprocessing.MaxAbsScaler()
    fps_train = SCALER.fit_transform(fps_train)
    fps_target = SCALER.transform(fps_target)
    if pca:
        # TruncatedSVD is more efficient than PCA for
        # sparse matrix
        pca = TruncatedSVD(n_components=n_components)
        fps_train = pca.fit_transform(fps_train)
        fps_target = pca.transform(fps_target)
    # distances matrix: each row represents the distance of one
    # target image to all training images
    distances_matrix = np.zeros([fps_target.shape[0], fps_train.shape[0]])
    for _ in range(fps_target.shape[0]):
        fps_stacked = np.vstack([fps_train, fps_target[_]])
        distances = ((weight * (np.abs(fps_stacked - fps_stacked[-1])**order))
                     .sum(axis=1))**(1 / order)
        distances = distances[:-1]
        distances_matrix[_] = distances
    if fps_train.shape[0] <= 2*fps_target.shape[0]:
        distances_matrix /= np.max(distances_matrix, axis=None)
    else:
        distances_matrix /= np.median(distances_matrix, axis=None)
    uncertainties = []
    for dist_array in distances_matrix:
        uncertainty = uncertainty_func(dist_array, scale=scale,
                                       method=weighted_method)
        uncertainties.append(uncertainty)
    return uncertainties


def uncertainty_func(x, scale=1., method='mean'):
    """
    function to calculate uncertainty. It takes as input the (weighted)
    vector distances between a target image and the training images.
    The default is taking average. The other two methods are

    - the maximum uncertainty

    - the average of the first n maximum uncertainty. By default n is 20% of
    training images.
    """
    assert method in ['mean', 'max', 'max_n']
    x = np.tanh(scale * x)
    if method == 'mean':
        return x.mean()
    elif method == 'max':
        return x.max()
    elif method == 'max_n':
        n = int(0.2 * len(x))
        x = np.sort(x)[-n:]
        return x.mean()
    else:
        raise RuntimeError('Methods not available!')


def get_debu(fps_train, fps_target, with_pca=True, n_components=10,
             kernel='gaussian', bandwidth=0.4, normalized=True,
             scale=1.):
    """kernel density-based uncertainty (DEBU) method.

    Parameters
    ------------
    fps_train : list or numpy array
        vectors as constructed from the training set.

    fps_target: list or numpy array
        vector(s) as constructed from the target images.

    with_pca: bool
        whether or not preprocessed with principal component
        analysis

    n_components: int
        number of component kept for PCA analysis

    kernel : str
        kernel to use. Valid kernels are ['gaussian'|'tophat'|
        'epanechnikov'|'linear']

    bandwidth: float or 'auto'
        bandwidth of the kernel. If 'auto' is indicated, it will
        use 1/(n_train_images) as the bandwidth
    
    normalized: bool
        if the uncertainties are normalized based on densities
        under the model

    scale: float (default: 1.)
        scale factor for density evaluation

    Return
    --------
    Density-based uncertainty levels for the target images
    """
    fps_train, fps_target = np.array(fps_train), np.array(fps_target)
    scaler = preprocessing.MinMaxScaler()
    fps_train = scaler.fit_transform(fps_train)
    fps_target = scaler.transform(fps_target)
    if fps_train.shape[0] < n_components:
        with_pca, normalized = False, False
        # use random query when pca is not performed 
        return np.ones(fps_target.shape[0])
    if with_pca:
        pca = TruncatedSVD(n_components=n_components)
        fps_train = pca.fit_transform(fps_train)
        fps_target = pca.transform(fps_target)
        variances_ratio = pca.explained_variance_ratio_
    n_cols = fps_train.shape[1]
    # Compute kernel densities based on training images
    # kdes: the kernel density estimators for each feature (column)
    # scores: the total log probability density under the model
    kdes, scores = [], []
    if bandwidth == 'auto':
        bandwidth = 1/fps_train.shape[0]
    for k in range(n_cols):
        kde = KernelDensity(kernel=kernel, bandwidth=bandwidth)
        X = fps_train[:, k].reshape(-1, 1)
        kde.fit(X)
        kdes.append(kde)
        scores.append(np.exp(kde.score_samples(X)).max())
    densities = np.zeros([fps_target.shape[0], n_cols])
    for j in range(n_cols):
        kde_per = kdes[j]
        target = fps_target[:,j].reshape(-1, 1)
        densities[:,j] = np.exp(kde_per.score_samples(target))
    debus = []
    densities /= np.max(densities, axis=None)
    for density_per in densities:
        if normalized:
            debu = 1 - np.tanh(scale*sum(density_per*variances_ratio) /
                   sum(np.array(scores) * variances_ratio))
            debus.append(debu)
        else:
            debus.append(1 - np.tanh((scale*density_per).mean()))
    return np.array(debus)
