import torch
from ase.neighborlist import neighbor_list

def get_distances(atoms, cutoff):
    first_atom_idx, second_atom_idx, shift_vector = neighbor_list(
        "ijS", atoms, cutoff, self_interaction=False
    )
    _, _, d = neighbor_list(
            "ijd", atoms, cutoff, self_interaction=False
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
    pairs = torch.cat(
        (first_atom_idx.reshape(-1, 1), second_atom_idx.reshape(-1, 1)), 1
    ).numpy()
    distance_vec = distance_vec.numpy()

    return pairs, distance_vec
