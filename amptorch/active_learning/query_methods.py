import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator as sp
import random
import copy
import warnings

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


def random_query(_, sample_candidates, samples_to_retrain, parent_calc):
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
