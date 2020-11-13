import hashlib

import numpy as np
from scipy import sparse

from ..base_descriptor import BaseDescriptor
from ..constants import ATOM_SYMBOL_TO_INDEX_DICT
from ..util import _gen_2Darray_for_ffi, list_symbols_to_indices
from ._libsymf import ffi, lib


class GaussianDescriptorSet:
    # TODO:
    #  - convert default settings to self.full-list descriptor set
    #  - convert element-wise settings to self.full-list descriptor set
    #  - validate full-list descriptor set and save as self.full-list descriptor set
    #  - convert self.full-list descriptor set to various formats required by Gaussian

    def __init__(
        self,
        elements,
        cutoff=6.5,
        cutoff_params={"cutoff_func": "cosine"},
        default_interactions=False,
    ):
        self.elements = elements
        self.element_indices = list_symbols_to_indices(elements)
        self.cutoff = cutoff
        self.cutoff_params = cutoff_params
        self.all_interactions = set()
        self.interactions = {
            element: {"G2": set(), "G4": set(), "G5": set()} for element in elements
        }
        self.descriptor_setup = None
        self.descriptor_setup_hash = None

    def add_g2(self, element_i, element_j, eta=3.0, rs=0.0, cutoff=None, update=True):
        assert element_i in self.elements, f"{element_i} is not in {self.elements}"
        assert element_j in self.elements, f"{element_j} is not in {self.elements}"
        g2_params = (
            2,
            self.element_indices.index(element_j),
            0,
            cutoff or self.cutoff,
            eta,
            rs,
            0.0,
        )
        self.interactions[element_i]["G2"].add(g2_params)
        if update:
            self.update()
        return self.interactions[element_i]["G2"]

    def add_g4(
        self,
        element_i,
        element_j,
        element_k,
        eta=0.005,
        zeta=1.0,
        gamma=1.0,
        cutoff=None,
        update=True,
    ):
        assert element_i in self.elements, f"{element_i} is not in {self.elements}"
        assert element_j in self.elements, f"{element_j} is not in {self.elements}"
        assert element_k in self.elements, f"{element_k} is not in {self.elements}"
        element_j, element_k = sorted([element_j, element_k])
        g4_params = (
            4,
            self.element_indices.index(element_j),
            self.element_indices.index(element_k),
            cutoff or self.cutoff,
            eta / cutoff ** 2,
            zeta,
            gamma,
        )
        self.interactions[element_i]["G4"].add(g4_params)
        if update:
            self.update()
        return self.interactions[element_i]["G4"]

    def add_g5(
        self,
        element_i,
        element_j,
        element_k,
        eta=0.005,
        zeta=1.0,
        gamma=1.0,
        cutoff=None,
        update=True,
    ):
        assert element_i in self.elements, f"{element_i} is not in {self.elements}"
        assert element_j in self.elements, f"{element_j} is not in {self.elements}"
        assert element_k in self.elements, f"{element_k} is not in {self.elements}"
        element_j, element_k = sorted([element_j, element_k])
        g5_params = (
            5,
            self.element_indices.index(element_j),
            self.element_indices.index(element_k),
            cutoff or self.cutoff,
            eta,
            zeta,
            gamma,
        )
        self.interactions[element_i]["G5"].add(g5_params)
        self.descriptor_setup = None
        self.descriptor_hash = None
        if update:
            self.update()
        return self.interactions[element_i]["G5"]

    def update(self):
        self.descriptor_setup = self._get_descriptor_setup()
        self.descriptor_hash = self._get_descriptor_hash()

    def process_combinatorial_Gs(self, Gs):
        for element in self.interactions.keys():
            if element in Gs:
                # self.interactions[element] =
                self._process_element_combinatorial_params(element, Gs[element])
            elif "default" in Gs:
                # self.interactions[element] =
                self._process_element_combinatorial_params(element, Gs["default"])
            else:
                raise NotImplementedError(
                    "Symmetry function parameters not defined properly - element {} passed but not present in {}".format(
                        element, self.elements
                    )
                )
        self.update()

    def _process_element_combinatorial_params(self, element_i, element_Gs):
        cutoff = element_Gs["cutoff"]
        if "G2" in element_Gs:
            for eta in np.array(element_Gs["G2"]["etas"]) / element_Gs["cutoff"] ** 2:
                for rs in element_Gs["G2"]["rs_s"]:
                    for element_j in self.elements:
                        self.add_g2(element_i, element_j, eta, rs, cutoff, update=False)

        if "G4" in element_Gs:
            for eta in np.array(element_Gs["G4"]["etas"]) / element_Gs["cutoff"] ** 2:
                for zeta in element_Gs["G4"]["zetas"]:
                    for gamma in element_Gs["G4"]["gammas"]:
                        for j, element_j in enumerate(self.elements):
                            for element_k in self.elements[j:]:
                                self.add_g4(
                                    element_i,
                                    element_j,
                                    element_k,
                                    eta,
                                    zeta,
                                    gamma,
                                    cutoff,
                                    update=False,
                                )

        if "G5" in element_Gs:
            for eta in element_Gs["G5"]["etas"]:
                for zeta in element_Gs["G5"]["zetas"]:
                    for gamma in element_Gs["G5"]["gammas"]:
                        for j, element_j in enumerate(self.elements):
                            for element_k in self.elements[j:]:
                                self.add_g5(
                                    element_i,
                                    element_j,
                                    element_k,
                                    eta,
                                    zeta,
                                    gamma,
                                    cutoff,
                                    update=False,
                                )

    def _get_descriptor_hash(self):
        string = (
            "cosine"
            if self.cutoff_params["cutoff_func"] == "cosine"
            else "polynomial%.15f" % self.cutoff_params["gamma"]
        )
        descriptor_setup = self.descriptor_setup
        for element in descriptor_setup.keys():
            string += element
            for desc in descriptor_setup[element]:
                for num in desc:
                    string += "%.15f" % num
        md5 = hashlib.md5(string.encode("utf-8"))
        hash_result = md5.hexdigest()
        return hash_result

    def _get_descriptor_setup(self):
        descriptor_setup = {element: None for element in self.elements}
        for element, descriptors in self.interactions.items():
            g2s, g4s, g5s = descriptors["G2"], descriptors["G4"], descriptors["G5"]
            g2s = [list(params) for params in sorted(g2s)]
            g4s = [list(params) for params in sorted(g4s)]
            g5s = [list(params) for params in sorted(g5s)]
            descriptor_setup[element] = np.array(g2s + g4s + g5s)
        return descriptor_setup

    def __eq__(self, other):
        return self.descriptor_setup_hash == other.descriptor_setup_hash

    def __hash__(self):
        return self.descriptor_setup_hash

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        interaction_counts = []
        for element, intxns in self.interactions.items():
            interaction_counts.append(
                "%s: {#G2: %d, #G4: %d, #G5: %d}"
                % (element, len(intxns["G2"]), len(intxns["G4"]), len(intxns["G5"]))
            )
        return "GaussianDescriptorSet(%s)" % "".join(interaction_counts)


class Gaussian(BaseDescriptor):
    def __init__(self, Gs, elements, cutoff_func="cosine", gamma=None):
        super().__init__()
        self.descriptor_type = "Gaussian"
        self.Gs = Gs
        self.elements = elements
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

    def prepare_descriptor_parameters(self):
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
                    for eta in np.array(Gs["G2"]["etas"]) / cutoff ** 2
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
                    for eta in (np.array(Gs["G4"]["etas"]) / cutoff ** 2)
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
        string = (
            "cosine" if self.cutoff_func == "cosine" else "polynomial%.15f" % self.gamma
        )
        for element in self.descriptor_setup.keys():
            string += element
            for desc in self.descriptor_setup[element]:
                for num in desc:
                    string += "%.15f" % num
        md5 = hashlib.md5(string.encode("utf-8"))
        hash_result = md5.hexdigest()
        self.descriptor_setup_hash = hash_result

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
        scale = np.copy(atoms.get_scaled_positions(), order="C")
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
