import hashlib
import numpy as np
from ase.io.trajectory import Trajectory
from .constants import ATOM_INDEX_TO_SYMBOL_DICT, ATOM_SYMBOL_TO_INDEX_DICT


def _gen_2Darray_for_ffi(arr, ffi, cdata="double"):
    # Function to generate 2D pointer for cffi  
    shape = arr.shape
    arr_p = ffi.new(cdata + " *[%d]" % shape[0])
    for i in range(shape[0]):
        arr_p[i] = ffi.cast(cdata + " *", arr[i].ctypes.data)
    return arr_p

def get_traj_hash(traj):
    
    # assert isinstance(traj, Trajectory)
    string = ""
    atoms = traj[0]
    string += str(atoms.pbc)
    try:
        flattened_cell = atoms.cell.array.flatten()
    except AttributeError:  # older ASE
        flattened_cell = atoms.cell.flatten()
    for number in flattened_cell:
        string += "%.15f" % number
    for number in atoms.get_atomic_numbers():
        string += "%3d" % number
    for number in atoms.get_positions().flatten():
        string += "%.15f" % number
    
    md5 = hashlib.md5(string.encode("utf-8"))
    hash_result = md5.hexdigest()

    return hash_result

def list_symbols_to_indices(list_of_symbols):
    list_indices = []
    for symbol in list_of_symbols:
        list_indices.append(ATOM_SYMBOL_TO_INDEX_DICT[symbol])
    return np.array(list_indices, dtype=np.intc)

def list_indices_to_symbols(list_of_indices):
    list_symbols = []
    for index in list_of_indices:
        list_symbols.append(ATOM_INDEX_TO_SYMBOL_DICT[index])
    return list_symbols