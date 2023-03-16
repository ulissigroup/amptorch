import hashlib

import numpy as np
from scipy import sparse

from ..base_descriptor import BaseDescriptor
from ..constants import ATOM_SYMBOL_TO_INDEX_DICT
from ..util import _gen_2Darray_for_ffi, list_symbols_to_indices
from ._libsymf import ffi, lib
from .descriptor_set import GaussianDescriptorSet


GDS = GaussianDescriptorSet  # so Flake8 allows the commit to proceed, also an easy-access-acronym!


class Gaussian(BaseDescriptor):
    """
    Fingerprinting calculation for Symmetry Functions.

    Args:
    Gs [dict] : a dictionary containing the definition of G2 and G4 parameters.

    elements [dict] : a dictionary of string of chemical elements in the system.

    cutoff_func [str] : defines the form of f_c. Default to "cosine".

    gamma [float] : parameter for polynomial cutoff function. Default to None as the default cutoff function is cosine.
    """

    def __init__(self, Gs, elements, cutoff_func="cosine", gamma=None):
        super().__init__()
        self.descriptor_type = "Gaussian"
        self.Gs = Gs
        self.elements = list(elements)
        self.cutoff_func = cutoff_func.lower()
        if self.cutoff_func not in ["cosine", "polynomial"]:
            raise ValueError('cutoff function must be either "cosine" or "polynomial"')
        if self.cutoff_func == "polynomial":
            if gamma is None:
                raise ValueError(
                    "polynomial cutoff function requires float value > 0. of `gamma`"
                )
            elif gamma <= 0.0:
                raise ValueError("polynomial cutoff function gamma must be > 0.")
        self.gamma = gamma
        self.element_indices = list_symbols_to_indices(elements)

        self.prepare_descriptor_parameters()
        self.get_descriptor_setup_hash()

    def __eq__(self, other):
        """Overrides the default implementation"""
        if isinstance(other, BaseDescriptor):
            if self.descriptor_type != other.descriptor_type:
                return False
            if self.elements != other.elements:
                return False
            for element in self.Gs:
                if element not in other.Gs:
                    return False
                if self.Gs[element]["cutoff"] != other.Gs[element]["cutoff"]:
                    return False
                if "G2" in self.Gs[element]:
                    for key in self.Gs[element]["G2"]:
                        if list(self.Gs[element]["G2"][key]) != list(
                            other.Gs[element]["G2"][key]
                        ):
                            return False
                if "G4" in self.Gs[element]:
                    for key in self.Gs[element]["G4"]:
                        if list(self.Gs[element]["G4"][key]) != list(
                            other.Gs[element]["G4"][key]
                        ):
                            return False
                if "G5" in self.Gs[element]:
                    for key in self.Gs[element]["G5"]:
                        if list(self.Gs[element]["G5"][key]) != list(
                            other.Gs[element]["G5"][key]
                        ):
                            return False
            if self.cutoff_func != other.cutoff_func:
                return False
            return True
        return NotImplemented

    def prepare_descriptor_parameters(self):
        """
        A helper function to prepare the parameters as input to cffi.
        """
        if isinstance(self.Gs, dict):
            self.descriptor_setup = {}
            for element in self.elements:
                if element in self.Gs:
                    self.descriptor_setup[
                        element
                    ] = self._prepare_descriptor_parameters_element(
                        self.Gs[element], self.element_indices
                    )
                elif "default" in self.Gs:
                    self.descriptor_setup[
                        element
                    ] = self._prepare_descriptor_parameters_element(
                        self.Gs["default"], self.element_indices
                    )
                else:
                    raise NotImplementedError(
                        "Symmetry function parameters not defined properly"
                    )
        elif hasattr(self.Gs, "descriptor_setup"):
            self.descriptor_setup = self.Gs.descriptor_setup
        else:
            raise ValueError(
                "Gs must be a dict with descriptor params or a GaussianDescriptorSet object: passed was a (%s)"
                % type(self.Gs)
            )
        self.params_set = dict()
        for element in self.elements:
            element_index = ATOM_SYMBOL_TO_INDEX_DICT[element]
            self.params_set[element_index] = dict()
            params_i = np.asarray(
                self.descriptor_setup[element][:, :3].copy(), dtype=np.intc, order="C"
            )
            params_d = np.asarray(
                self.descriptor_setup[element][:, 3:].copy(),
                dtype=np.float64,
                order="C",
            )
            self.params_set[element_index]["i"] = params_i
            self.params_set[element_index]["d"] = params_d
            self.params_set[element_index]["ip"] = _gen_2Darray_for_ffi(
                self.params_set[element_index]["i"], ffi, "int"
            )
            self.params_set[element_index]["dp"] = _gen_2Darray_for_ffi(
                self.params_set[element_index]["d"], ffi
            )
            self.params_set[element_index]["total"] = np.concatenate(
                (
                    self.params_set[element_index]["i"],
                    self.params_set[element_index]["d"],
                ),
                axis=1,
            )
            self.params_set[element_index]["num"] = len(self.descriptor_setup[element])

        return

    def _prepare_descriptor_parameters_element(self, Gs, element_indices):
        descriptor_setup = {"G2": set(), "G4": set(), "G5": set()}
        cutoff = Gs["cutoff"]
        if "G2" in Gs:
            descriptor_setup["G2"].update(
                [
                    (2, element1, 0, cutoff, eta, rs, 0.0)
                    for eta in np.array(Gs["G2"]["etas"]) / cutoff**2
                    for rs in Gs["G2"]["rs_s"]
                    for element1 in element_indices
                ]
            )

        if "G4" in Gs:
            descriptor_setup["G4"].update(
                [
                    (
                        4,
                        element_indices[i],
                        element_indices[j],
                        cutoff,
                        eta,
                        zeta,
                        gamma,
                    )
                    for eta in (np.array(Gs["G4"]["etas"]) / cutoff**2)
                    for zeta in Gs["G4"]["zetas"]
                    for gamma in Gs["G4"]["gammas"]
                    for i in range(len(element_indices))
                    for j in range(i, len(element_indices))
                ]
            )

        if "G5" in Gs:
            descriptor_setup["G5"].update(
                [
                    (
                        5,
                        element_indices[i],
                        element_indices[j],
                        cutoff,
                        eta,
                        zeta,
                        gamma,
                    )
                    for eta in Gs["G5"]["etas"]
                    for zeta in Gs["G5"]["zetas"]
                    for gamma in Gs["G5"]["gammas"]
                    for i in range(len(element_indices))
                    for j in range(i, len(element_indices))
                ]
            )
        g2s, g4s, g5s = (
            descriptor_setup["G2"],
            descriptor_setup["G4"],
            descriptor_setup["G5"],
        )
        g2s = [list(params) for params in sorted(g2s)]
        g4s = [list(params) for params in sorted(g4s)]
        g5s = [list(params) for params in sorted(g5s)]
        descriptor_setup = np.array(g2s + g4s + g5s)
        return descriptor_setup

    def get_descriptor_setup_hash(self):
        if isinstance(self.Gs, dict):
            string = (
                "cosine"
                if self.cutoff_func == "cosine"
                else "polynomial%.15f" % self.gamma
            )
            for element in self.descriptor_setup.keys():
                string += element
                for desc in self.descriptor_setup[element]:
                    for num in desc:
                        string += "%.15f" % num
            md5 = hashlib.md5(string.encode("utf-8"))
            hash_result = md5.hexdigest()
            self.descriptor_setup_hash = hash_result
        elif hasattr(self.Gs, "descriptor_setup_hash"):
            self.descriptor_setup_hash = self.Gs.descriptor_setup_hash
        else:
            raise ValueError(
                "Gs must be a dict with descriptor params or a GaussianDescriptorSet object: passed was a (%s)"
                % type(self.Gs)
            )

    def save_descriptor_setup(self, filename):
        with open(filename, "w") as out_file:
            for element in self.descriptor_setup.keys():
                out_file.write(
                    "===========\nElement: {} \t num_desc: {}\n".format(
                        element, len(self.descriptor_setup[element])
                    )
                )
                for desc in self.descriptor_setup[element]:
                    out_file.write(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                            int(desc[0]),
                            int(desc[1]),
                            int(desc[2]),
                            desc[3],
                            desc[4],
                            desc[5],
                            desc[6],
                        )
                    )

    def calculate_fingerprints(self, atoms, element, calc_derivatives, log):
        """
        Interfacing with cffi to compute the fingerprints.
        """
        element_index = ATOM_SYMBOL_TO_INDEX_DICT[element]

        symbols = np.array(atoms.get_chemical_symbols())
        atom_num = len(symbols)
        atom_indices = list_symbols_to_indices(symbols)
        unique_atom_indices = np.unique(atom_indices)

        type_num = dict()
        type_idx = dict()

        for atom_index in unique_atom_indices:
            tmp = atom_indices == atom_index
            type_num[atom_index] = np.sum(tmp).astype(np.int64)
            # if atom indexs are sorted by atom type,
            # indexs are sorted in this part.
            # if not, it could generate bug in training process for force training
            type_idx[atom_index] = np.arange(atom_num)[tmp]

        atom_indices_p = ffi.cast("int *", atom_indices.ctypes.data)

        cart = np.copy(atoms.get_positions(wrap=True), order="C")
        scale = np.copy(atoms.get_scaled_positions(wrap=True), order="C")
        cell = np.copy(atoms.cell, order="C")
        pbc = np.copy(atoms.get_pbc()).astype(np.intc)

        cart_p = _gen_2Darray_for_ffi(cart, ffi)
        scale_p = _gen_2Darray_for_ffi(scale, ffi)
        cell_p = _gen_2Darray_for_ffi(cell, ffi)
        pbc_p = ffi.cast("int *", pbc.ctypes.data)

        cal_atoms = np.asarray(type_idx[element_index], dtype=np.intc, order="C")
        cal_num = len(cal_atoms)
        cal_atoms_p = ffi.cast("int *", cal_atoms.ctypes.data)

        size_info = np.array([atom_num, cal_num, self.params_set[element_index]["num"]])

        if calc_derivatives:
            x = np.zeros(
                [cal_num, self.params_set[element_index]["num"]],
                dtype=np.float64,
                order="C",
            )
            dx = np.zeros(
                [cal_num * self.params_set[element_index]["num"], atom_num * 3],
                dtype=np.float64,
                order="C",
            )

            x_p = _gen_2Darray_for_ffi(x, ffi)
            dx_p = _gen_2Darray_for_ffi(dx, ffi)

            errno = (
                lib.calculate_sf_cos(
                    cell_p,
                    cart_p,
                    scale_p,
                    pbc_p,
                    atom_indices_p,
                    atom_num,
                    cal_atoms_p,
                    cal_num,
                    self.params_set[element_index]["ip"],
                    self.params_set[element_index]["dp"],
                    self.params_set[element_index]["num"],
                    x_p,
                    dx_p,
                )
                if self.cutoff_func == "cosine"
                else lib.calculate_sf_poly(
                    cell_p,
                    cart_p,
                    scale_p,
                    pbc_p,
                    atom_indices_p,
                    atom_num,
                    cal_atoms_p,
                    cal_num,
                    self.params_set[element_index]["ip"],
                    self.params_set[element_index]["dp"],
                    self.params_set[element_index]["num"],
                    x_p,
                    dx_p,
                    self.gamma,
                )
            )

            if errno == 1:
                raise NotImplementedError("Descriptor not implemented!")
            fp = np.array(x)
            fp_prime = np.array(dx)
            scipy_sparse_fp_prime = sparse.coo_matrix(fp_prime)

            return (
                size_info,
                fp,
                scipy_sparse_fp_prime.data,
                scipy_sparse_fp_prime.row,
                scipy_sparse_fp_prime.col,
                np.array(fp_prime.shape),
            )

        else:
            x = np.zeros(
                [cal_num, self.params_set[element_index]["num"]],
                dtype=np.float64,
                order="C",
            )
            x_p = _gen_2Darray_for_ffi(x, ffi)

            errno = (
                lib.calculate_sf_cos_noderiv(
                    cell_p,
                    cart_p,
                    scale_p,
                    pbc_p,
                    atom_indices_p,
                    atom_num,
                    cal_atoms_p,
                    cal_num,
                    self.params_set[element_index]["ip"],
                    self.params_set[element_index]["dp"],
                    self.params_set[element_index]["num"],
                    x_p,
                )
                if self.cutoff_func == "cosine"
                else lib.calculate_sf_poly_noderiv(
                    cell_p,
                    cart_p,
                    scale_p,
                    pbc_p,
                    atom_indices_p,
                    atom_num,
                    cal_atoms_p,
                    cal_num,
                    self.params_set[element_index]["ip"],
                    self.params_set[element_index]["dp"],
                    self.params_set[element_index]["num"],
                    x_p,
                    dx_p,
                    self.gamma,
                )
            )

            if errno == 1:
                raise NotImplementedError("Descriptor not implemented!")
            fp = np.array(x)

            return size_info, fp, None, None, None, None
