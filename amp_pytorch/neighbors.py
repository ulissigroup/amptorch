import sys
import time
import numpy as np
import torch
from ase.neighborlist import neighbor_list
from ase import Atoms
from amp.utilities import hash_images
from amp.descriptor.gaussian import Gaussian

images = [
    Atoms(
        symbols="PdOPd",
        pbc=np.array([False, False, False], dtype=bool),
        cell=np.array([[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]]),
        positions=np.array([[0.5, 1.0, 0.5], [1.0, 0.5, 1.0], [1.5, 1.5, 1.5]]),
    ),
    Atoms(
        symbols="PdO",
        pbc=np.array([False, False, False], dtype=bool),
        cell=np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        positions=np.array([[2.0, 1.0, -1.0], [1.0, 2.0, 1.0]]),
    ),
]

cutoff = 3
descriptor = Gaussian(cutoff)
atoms = images[0]
hashed_image = hash_images(images)
descriptor.calculate_fingerprints(hashed_image)


def get_distances(atoms, cutoff):
    first_atom_idx, second_atom_idx, shift_vector = neighbor_list(
        "ijS", atoms, cutoff, self_interaction=False
    )
    first_atom_idx = torch.tensor(first_atom_idx)
    second_atom_idx = torch.tensor(second_atom_idx)
    shift_vector = torch.FloatTensor(shift_vector)
    atom_positions = torch.FloatTensor(atoms.positions)
    cell = torch.FloatTensor(atoms.cell)
    first_atom_positions = torch.index_select(atom_positions, 0, first_atom_idx)
    neighbor_positions = torch.index_select(atom_positions, 0, second_atom_idx)
    distance_vec = neighbor_positions - first_atom_positions
    offsets = torch.mm(shift_vector, cell)
    distance_vec += offsets
    distances = torch.norm(distance_vec, dim=1).numpy().reshape(1, -1)

    pairs = torch.cat(
        (first_atom_idx.reshape(-1, 1), second_atom_idx.reshape(-1, 1)), 1
    ).numpy()

    return pairs, distances


pairs, distances = get_distances(atoms, cutoff)


# SchNet implementation of distance calcs
# def get_environment(atoms, cutoff, grid=None):
    # if grid is not None:
        # raise NotImplementedError

    # n_atoms = atoms.get_number_of_atoms()
    # idx_i, idx_j, idx_S = neighbor_list("ijS", atoms, cutoff, self_interaction=False)
    # if idx_i.shape[0] > 0:
        # uidx, n_nbh = np.unique(idx_i, return_counts=True)
        # n_max_nbh = np.max(n_nbh)

        # n_nbh = np.tile(n_nbh[:, np.newaxis], (1, n_max_nbh))
        # nbh_range = np.tile(
            # np.arange(n_max_nbh, dtype=np.int)[np.newaxis], (n_nbh.shape[0], 1)
        # )

        # mask = np.zeros((n_atoms, np.max(n_max_nbh)), dtype=np.bool)
        # mask[uidx, :] = nbh_range < n_nbh
        # neighborhood_idx = -np.ones((n_atoms, np.max(n_max_nbh)), dtype=np.float32)
        # neighborhood_idx[mask] = idx_j

        # offset = np.zeros((n_atoms, np.max(n_max_nbh), 3), dtype=np.float32)
        # offset[mask] = idx_S
    # else:
        # neighborhood_idx = -np.ones((n_atoms, 1), dtype=np.float32)
        # offset = np.zeros((n_atoms, 1, 3), dtype=np.float32)

    # return neighborhood_idx, offset


# neighbors, offset = get_environment(atoms, cutoff)
# positions = torch.FloatTensor(atoms.positions).reshape(1, 3, 3)
# neighbors = torch.LongTensor(neighbors).reshape(1, 3, -1)
# cell = torch.FloatTensor(atoms.cell).reshape(1, 3, 3)
# offset = torch.FloatTensor(offset).reshape(1, 3, -1, 3)
# def get_distance_schnet(positions, neighbors, cell=None, cell_offsets=None):
	# n_batch = positions.size()[0]
	# idx_m = torch.arange(n_batch, dtype=torch.long)[:,
			# None, None]
    # Get atomic positions of all neighboring indices
	# pos_xyz = positions[idx_m, neighbors[:, :, :], :]

    # Subtract positions of central atoms to get distance vectors
	# dist_vec = pos_xyz - positions[:, :, None, :]

    # add cell offset
	# if cell is not None:
		# B, A, N, D = cell_offsets.size()
		# cell_offsets = cell_offsets.view(B, A * N, D)
		# offsets = cell_offsets.bmm(cell)
		# offsets = offsets.view(B, A, N, D)
		# dist_vec += offsets

    # Compute vector lengths
	# distances = torch.norm(dist_vec, 2, 3)
	# return distances

# distance_schnet = get_distance_schnet(positions, neighbors, cell, offset)
