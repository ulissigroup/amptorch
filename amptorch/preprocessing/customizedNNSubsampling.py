# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 14:36:48 2017
@author: Xiangyun Lei
"""

from __future__ import print_function

import random
import time
from collections import Counter
from math import ceil, sqrt

import numpy as np
from sklearn import preprocessing
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

try:
    import cPickle as pickle
except:
    import pickle

try:
    from pyflann import *
except:
    pass
try:
    from pykdtree.kdtree import KDTree
except:
    pass
try:
    from sklearn.neighbors import NearestNeighbors
except:
    pass
try:
    from annoy import AnnoyIndex
except:
    pass
try:
    from scipy.spatial import cKDTree
except:
    pass
try:
    import nmslib
except:
    pass


"""
Helper functions
"""


def get_data_process(li, list_of_index):
    """
    Select features of a list based on index set
    """
    result = []
    for entry in li:
        result.append([entry[i] for i in list_of_index])
    return result


def get_array_based_on_index(li, list_of_index):
    """
    Select entries from a list based on index set
    """
    return np.asarray([li[i] for i in list_of_index])


def remove_list_from_list(a, b):
    """
    Remove entries in list a that's also in list b
    """
    return list(set(a) - set(b))


def chunker(seq, size):
    """
    break a list (seq) into a set of equally sized (size) lists
    """
    return (seq[pos : pos + size] for pos in range(0, len(seq), size))


"""
subsample based on kD-tree
"""


def get_subsampling_index2(
    data_process,
    standard_scale=True,
    cutoff_sig=0.02,
    rate=0.3,
    method="pykdtree",
    verbose=1,
    image_index=[],
):

    """
    Using Nearest-Neighbor search based algorithm, find the list of indices of the subsampled dataset


    Parameters
    -------------
    data_process: List. the list of datapoints, with selected features

    standard_scale [True]: Boolean. Whether to apply standard scaler to the dataset prior to subsampling

    cutoff_sig [0.02]: Float. cutoff significance. the cutoff distance equals to the Euclidean
                       norm of the standard deviations in all dimensions of the data points

    rate [0.3]: Float. possibility of deletion

    method ["pykdtree"]: String. which backend nearest neighbour model to use.
                         possible choices: ["pykdtree", "nmslib", "sklearn", "scipy", "annoy", "flann"]

    verbose [1]: integer. level of verbosity


    Return
    -------------
    overall_keep_list: The list of indices of the final subsampled entries

    """

    if verbose >= 1:
        print("Started NN-subsampling, original length: {}".format(len(data_process)))

    method = method.lower()
    start = time.time()

    if method == "flann":
        if verbose >= 1:
            print("use flann backend")
    elif method == "pykdtree":
        if verbose >= 1:
            print("use pykdtree backend")
    elif method == "sklearn":
        if verbose >= 1:
            print("use slearn nearest neighbors backend")
    elif method == "scipy":
        if verbose >= 1:
            print("use scipy cKDTree backend")
    elif method == "annoy":
        if verbose >= 1:
            print("use annoy backend")
    elif method == "nmslib":
        if verbose >= 1:
            print("use nmslib backend")
    else:
        print("method {} not impletemented".format(method))
        raise NotImplemented

    # apply standard scaling
    if standard_scale:
        if verbose >= 2:
            print("Subample with standard scaled data")
        data_process = StandardScaler().fit_transform(np.asarray(data_process).copy())
    else:
        if verbose >= 2:
            print("Subample with original data")
        data_process = np.asarray(data_process).copy()

    # set cutoff distance
    list_of_descs = zip(*data_process)
    sum_std2 = 0.0
    for descs in list_of_descs:
        temp_std = np.std(descs)
        sum_std2 += temp_std ** 2
    cutoff = cutoff_sig * np.sqrt(sum_std2)

    # initialize the index
    overall_keep_list = np.arange(len(data_process)).tolist()

    keep_going = True
    iter_count = 1
    while keep_going:
        if verbose >= 2:
            print(
                "start iteration {}, total length: {}".format(
                    iter_count, len(overall_keep_list)
                )
            )
        start_cycle = time.time()
        temp_data_process = get_array_based_on_index(
            data_process.copy(), overall_keep_list
        )
        temp_image_index = get_array_based_on_index(image_index, overall_keep_list)

        # build and query nearest neighbour model
        if method == "flann":
            flann = FLANN()
            indices, distances = flann.nn(
                temp_data_process, temp_data_process, 2, algorithm="kmeans"
            )
        elif method == "scipy":
            kd_tree = cKDTree(temp_data_process)
            distances, indices = kd_tree.query(temp_data_process, k=2)
        elif method == "pykdtree":
            kd_tree = KDTree(temp_data_process, leafsize=6)
            distances, indices = kd_tree.query(temp_data_process, k=2)
        elif method == "sklearn":
            nbrs = NearestNeighbors(n_neighbors=2, algorithm="kd_tree", n_jobs=-1).fit(
                temp_data_process
            )
            distances, indices = nbrs.kneighbors(temp_data_process)
        elif method == "annoy":
            annoy = AnnoyIndex(len(temp_data_process[0]), metric="euclidean")
            for i in range(len(temp_data_process)):
                annoy.add_item(i, temp_data_process[i])
            annoy.build(1)
            distances = []
            indices = []
            for i in range(len(temp_data_process)):
                temp_index, temp_dist = annoy.get_nns_by_vector(
                    temp_data_process[i], 2, include_distances=True
                )
                indices.append([i, temp_index[1]])
                distances.append([0.0, temp_dist[1]])
        elif method == "nmslib":
            index = nmslib.init(method="hnsw", space="l2")
            index.addDataPointBatch(temp_data_process)
            index.createIndex(print_progress=False)

            neighbours = index.knnQueryBatch(temp_data_process, k=2)

            distances = []
            indices = []
            for item in neighbours:
                indices.append(item[0])
                distances.append(item[1])

        else:
            raise NotImplemented

        # if distance between each point and its nearest neighbor is below cutoff distance,
        # add the nearest neighbout to the candidate removal list
        remove_index_li = []
        index_li = []

        for index, distance in zip(indices, distances):
            index_li.append(index[0])
            if distance[1] <= cutoff:
                remove_index_li.append(index[1])

        # randomly select datapoints in the candidate removal list (based on rate)
        # and form the final removal list of this iteration
        # stop the cycle if the final removal list is empty
        temp_num = int(ceil(float(len(remove_index_li)) * rate))

        if temp_num == 0:
            keep_going = False
        # remove_index_li = random_subsampling(remove_index_li,temp_num)
        remove_index_li = rank_subsampling(remove_index_li, temp_num, temp_image_index)

        temp_keep_list = remove_list_from_list(index_li, remove_index_li)
        overall_keep_list = [overall_keep_list[i] for i in temp_keep_list]
        try:
            if len(overall_keep_list) == old_overall_keep_list_len:
                keep_going = False
                print("stopped because length didn't change")
        except:
            pass
        if verbose >= 2:
            print(
                "end iteration {}. length: {}\t time:{}".format(
                    iter_count, len(overall_keep_list), time.time() - start_cycle
                )
            )
        iter_count += 1
        old_overall_keep_list_len = len(overall_keep_list)
    if verbose >= 1:
        print(
            "end NN-subsampling. length: {}\t time:{}".format(
                len(overall_keep_list), time.time() - start
            )
        )
    return overall_keep_list


def subsampling(
    data,
    list_desc=[],
    standard_scale=True,
    cutoff_sig=0.05,
    rate=0.3,
    method="pykdtree",
    verbose=1,
    image_index=[],
):

    """
    Run the NN-based subsampling algorithm to a list of data points and
    return the resulting list of subsampled data points

    Parameters
    -------------
    data: List. the original list of data points

    list_desc [[] (empty list)]:  List.
                the indices of descriptors (features) of the datapoints.
                The algorithm would subsample based only on these descriptors
                (although other features will still be kept in the resulting subsampled dataset)
                If the list is empty, then all feature will be taken into account

    standard_scale [True]: Boolean. Whether to apply standard scaler to the dataset prior to subsampling

    cutoff_sig [0.02]: Float. cutoff significance. the cutoff distance equals to the Euclidean
                       norm of the standard deviations in all dimensions of the data points

    rate [0.3]: Float. possibility of deletion

    method ["pykdtree"]: String. which backend nearest neighbour model to use.
                         possible choices: ["pykdtree", "nmslib", "sklearn", "scipy", "annoy", "flann"]

    verbose [1]: integer. level of verbosity


    Return
    -------------
    sampling_result : the result list of subsampled data points


    """

    if len(list_desc) == 0:
        data_process = data
    else:
        data_process = get_data_process(data, list_desc)

    overall_keep_list = get_subsampling_index2(
        data_process,
        standard_scale=standard_scale,
        cutoff_sig=cutoff_sig,
        rate=rate,
        method=method,
        verbose=verbose,
        image_index=image_index,
    )
    sampling_result = [data[i] for i in overall_keep_list]
    image_index_result = [image_index[i] for i in overall_keep_list]
    return sampling_result, image_index_result


def subsampling_with_PCA(
    data,
    list_desc=[],
    standard_scale=True,
    cutoff_sig=0.05,
    rate=0.3,
    start_trial_component=10,
    max_component=30,
    target_variance=0.999999,
    method="pykdtree",
    verbose=1,
):

    """
    Run the NN-based subsampling algorithm to a list of data points and
    return the resulting list of subsampled data points

    The data set will first be transformed by PCA, before running the subsampling algorithm
    The number of PCs kept is the minimal number of PCs that have sum explained variance
    greater than target_variance

    Note that the final resulting list of datapoints (sampling_result) is NOT transformed
    (since we only used the PCA + subsampling alghorithm to find the indices of the datapoints to be kept)

    Parameters
    -------------
    data: List. the original list of data points

    list_desc [[] (empty list)]:  List.
                the indices of descriptors (features) of the datapoints.
                The algorithm would subsample based only on these descriptors
                (although other features will still be kept in the resulting subsampled dataset)
                If the list is empty, then all feature will be taken into account

    standard_scale [True]: Boolean. Whether to apply standard scaler to the dataset prior to subsampling

    cutoff_sig [0.02]: Float. cutoff significance. the cutoff distance equals to the Euclidean
                       norm of the standard deviations in all dimensions of the data points

    rate [0.3]: Float. possibility of deletion

    start_trial_component [10]: Int. minimum number of PCs.
                           if the number of features is below this number, then all features will be kept

    max_component [30]: Int.the maximum number of PCs to be kept,
                        even the target variance has not been reached

    target_variance [0.999999]: Float. the target sum of variance.

    method ["pykdtree"]: String. which backend nearest neighbour model to use.
                         possible choices: ["pykdtree", "nmslib", "sklearn", "scipy", "annoy", "flann"]

    verbose [1]: integer. level of verbosity


    Return
    -------------
    sampling_result : the result list of subsampled data points


    """

    if len(list_desc) == 0:
        data_process = data
    else:
        data_process = get_data_process(data, list_desc)

    if verbose >= 1:
        print("start trial PCA")
    start = time.time()
    pca = PCA(svd_solver="randomized")
    data_pca = pca.fit_transform(data)
    explained_variance_ratio = pca.explained_variance_ratio_

    # determine how many PCs to be kept
    trial_component = start_trial_component - 1
    keep_going = True
    while keep_going:
        trial_component += 1
        sum_explained_variance = sum(explained_variance_ratio[:trial_component])

        if verbose >= 2:
            print(
                "trial components: {} \t explained variance: {}".format(
                    trial_component, sum_explained_variance
                )
            )

        if sum_explained_variance > target_variance:
            keep_going = False

        if trial_component > max_component:
            keep_going = False
            if verbose >= 2:
                print(
                    "stopped PCA at {} components, total explained variance: {}".format(
                        trial_component, sum_explained_variance
                    )
                )

        if trial_component >= len(data_process[0]):
            keep_going = False
            pca_result = data_process

    if verbose >= 1:
        print(
            "end trial PCA, number of PC kept: {} \t took {} s".format(
                trial_component, str(time.time() - start)
            )
        )

    # pass the PCA transformed data to subsampling algorithm, and find the indices of data points to be kept
    overall_keep_list = get_subsampling_index2(
        data_pca[:, :trial_component],
        standard_scale=standard_scale,
        cutoff_sig=cutoff_sig,
        rate=rate,
        method=method,
        verbose=verbose,
    )

    sampling_result = [data[i] for i in overall_keep_list]
    return sampling_result


def batch_subsampling(
    data,
    list_desc=[],
    batch_size=1000000,
    recursive_level=1,
    standard_scale=True,
    cutoff_sig=0.05,
    rate=0.3,
    method="pykdtree",
    verbose=1,
    shuffle=True,
    final_overall_subsample=True,
):

    """
    Subsample with batch
    This is to save the memory if the data set of interest is too large.

    The data set will first be broken down into equally sized batchs (defined by batch size)
    that will be subsampled individually.
    The resulting subsampled datapoints will then be pooled together for a overall subsample

    In case this is not sufficient, multi-level batch subsampling is also allowed.
    So instead of a overall subsample after pooling the resulting subsampled datapoints,
    the pooled data points will again be broken down into batches for a second-level batch subsample.
    This process is repeated multiple times (as defined by recursive_level) before eventually an
    overall subsample is performed


    Parameters
    -------------
    data: List. the original list of data points

    list_desc [[] (empty list)]:  List.
                the indices of descriptors (features) of the datapoints.
                The algorithm would subsample based only on these descriptors
                (although other features will still be kept in the resulting subsampled dataset)
                If the list is empty, then all feature will be taken into account

    batch_size [1000000]: Int. the number of datapoints in each batch

    recursive_level [1]: Int. the number of levels for batch subsampling (as described above)
    standard_scale [True]: Boolean. Whether to apply standard scaler to the dataset prior to subsampling

    cutoff_sig [0.02]: Float. cutoff significance. the cutoff distance equals to the Euclidean
                       norm of the standard deviations in all dimensions of the data points

    rate [0.3]: Float. possibility of deletion

    method ["pykdtree"]: String. which backend nearest neighbour model to use.
                         possible choices: ["pykdtree", "nmslib", "sklearn", "scipy", "annoy", "flann"]

    verbose [1]: integer. level of verbosity

    shuffle [True]: Boolean. whether to shuffle the dataset before breaking down into batchs

    final_overall_subsample [True]: wheter to perform a overall subsample after the levels of batch
                                     subsample are finished

    Return
    -------------
    sampling_result : the result list of subsampled data points

    """

    if verbose >= 1:
        print("\n\nat recursive level {}, length {}".format(recursive_level, len(data)))

    sampling_result = []

    if shuffle:
        random.shuffle(data)

    for data_subgroup in chunker(data, batch_size):
        temp_sampling_result = subsampling(
            data_subgroup,
            list_desc=[],
            standard_scale=standard_scale,
            cutoff_sig=cutoff_sig,
            rate=rate,
            method=method,
            verbose=verbose,
        )
        sampling_result += temp_sampling_result

    if recursive_level == 1:
        if verbose >= 1:
            print(
                "at recursive level 1, length {}, Overall subsample".format(
                    recursive_level, len(sampling_result)
                )
            )
        if final_overall_subsample:
            sampling_result = subsampling(
                sampling_result,
                list_desc=[],
                standard_scale=standard_scale,
                cutoff_sig=cutoff_sig,
                rate=rate,
                method=method,
                verbose=verbose,
            )

    else:
        if verbose >= 1:
            print(
                "end recursive level {}, length {}, Continue".format(
                    recursive_level, len(sampling_result)
                )
            )
        sampling_result = batch_subsampling(
            sampling_result,
            list_desc=[],
            batch_size=batch_size,
            recursive_level=recursive_level - 1,
            standard_scale=standard_scale,
            cutoff_sig=cutoff_sig,
            rate=rate,
            method=method,
            verbose=verbose,
        )
    return sampling_result


def batch_subsampling_with_PCA(
    data,
    list_desc=[],
    batch_size=1000000,
    recursive_level=1,
    start_trial_component=10,
    max_component=30,
    target_variance=0.999999,
    standard_scale=True,
    cutoff_sig=0.05,
    rate=0.3,
    method="pykdtree",
    verbose=1,
    shuffle=True,
    final_overall_subsample=True,
):

    """
    Subsample with batch (with PCA pre-processing)
    This is to save the memory if the data set of interest is too large.

    The data set will first be broken down into equally sized batchs (defined by batch size)
    that will be subsampled individually.
    The resulting subsampled datapoints will then be pooled together for a overall subsample

    In case this is not sufficient, multi-level batch subsampling is also allowed.
    So instead of a overall subsample after pooling the resulting subsampled datapoints,
    the pooled data points will again be broken down into batches for a second-level batch subsample.
    This process is repeated multiple times (as defined by recursive_level) before eventually an
    overall subsample is performed


    Parameters
    -------------
    data: List. the original list of data points

    list_desc [[] (empty list)]:  List.
                the indices of descriptors (features) of the datapoints.
                The algorithm would subsample based only on these descriptors
                (although other features will still be kept in the resulting subsampled dataset)
                If the list is empty, then all feature will be taken into account

    batch_size [1000000]: Int. the number of datapoints in each batch

    recursive_level [1]: Int. the number of levels for batch subsampling (as described above)
    standard_scale [True]: Boolean. Whether to apply standard scaler to the dataset prior to subsampling

    cutoff_sig [0.02]: Float. cutoff significance. the cutoff distance equals to the Euclidean
                       norm of the standard deviations in all dimensions of the data points

    rate [0.3]: Float. possibility of deletion

    start_trial_component [10]: Int. minimum number of PCs.
                           if the number of features is below this number, then all features will be kept

    max_component [30]: Int.the maximum number of PCs to be kept,
                        even the target variance has not been reached

    target_variance [0.999999]: Float. the target sum of variance.

    method ["pykdtree"]: String. which backend nearest neighbour model to use.
                         possible choices: ["pykdtree", "nmslib", "sklearn", "scipy", "annoy", "flann"]

    verbose [1]: integer. level of verbosity

    shuffle [True]: Boolean. whether to shuffle the dataset before breaking down into batchs
    final_overall_subsample [True]: wheter to perform a overall subsample after the levels of batch
                                     subsample are finished

    Return
    -------------
    sampling_result : the result list of subsampled data points

    """
    if verbose >= 1:
        print("at recursive level {}, length {}".format(recursive_level, len(data)))

    sampling_result = []

    if shuffle:
        random.shuffle(data)

    for data_subgroup in chunker(data, batch_size):
        temp_sampling_result = subsampling_with_PCA(
            data_subgroup,
            list_desc=[],
            standard_scale=standard_scale,
            start_trial_component=start_trial_component,
            max_component=max_component,
            target_variance=target_variance,
            cutoff_sig=cutoff_sig,
            rate=rate,
            method=method,
            verbose=verbose,
        )
        sampling_result += temp_sampling_result

    if recursive_level == 1:
        if verbose >= 1:
            print(
                "at recursive level 1, length {}, Overall subsample".format(
                    recursive_level, len(sampling_result)
                )
            )
        if final_overall_subsample:
            sampling_result = subsampling_with_PCA(
                sampling_result,
                list_desc=[],
                standard_scale=standard_scale,
                start_trial_component=start_trial_component,
                max_component=max_component,
                target_variance=target_variance,
                cutoff_sig=cutoff_sig,
                rate=rate,
                method=method,
                verbose=verbose,
            )

    else:
        if verbose >= 1:
            print(
                "end recursive level {}, length {}, Continue".format(
                    recursive_level, len(sampling_result)
                )
            )
        sampling_result = batch_subsampling_with_PCA(
            sampling_result,
            list_desc=[],
            batch_size=batch_size,
            start_trial_component=start_trial_component,
            max_component=max_component,
            target_variance=target_variance,
            recursive_level=recursive_level - 1,
            standard_scale=standard_scale,
            cutoff_sig=cutoff_sig,
            rate=rate,
            method=method,
            verbose=verbose,
        )
    return sampling_result


def random_subsampling(li, num):
    """
    Randomly select certain number (num) of entries from a list
    """
    if num == 0:
        return []
    elif len(li) > num:
        # generate random list of indecies
        index_list = random.sample(range(0, len(li) - 1), num)
        return [li[i] for i in index_list]

    else:
        print("no sampling")
        return li


def rank_subsampling(li, num, image_index):
    """
    Randomly select certain number (num) of entries from a list
    """

    image_index_values, image_index_counts = np.unique(image_index, return_counts=True)

    sorted_index_values = image_index_values[image_index_counts.argsort()]
    sorted_index_counts = image_index_counts[image_index_counts.argsort()]

    keep_image_list = []
    total = 0
    for i, count in enumerate(sorted_index_counts):
        if total + count < num:
            keep_image_list.append(sorted_index_values[i])
            total += count
        elif total + count == num:
            keep_image_list.append(sorted_index_values[i])
            remaining_border_case = 0
            total += count
            break

        else:
            keep_image_list.append(sorted_index_values[i])
            remaining_border_case = total + count - num
            total += count
            break

    result = []
    for i, entry in enumerate(li):
        if image_index[i] in keep_image_list:
            result.append(entry)

        if remaining_border_case > 0 and image_index[i] == keep_image_list[-1]:
            remaining_border_case -= 1
            if remaining_border_case == 0:
                keep_image_list.pop(-1)

    return result
