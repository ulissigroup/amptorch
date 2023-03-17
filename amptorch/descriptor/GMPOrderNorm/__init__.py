import hashlib

import numpy as np
from scipy import sparse

from ..base_descriptor import BaseDescriptor
from ..constants import ATOM_SYMBOL_TO_INDEX_DICT
from ..util import _gen_2Darray_for_ffi, list_symbols_to_indices
from ._libgmpordernorm import ffi, lib


class GMPOrderNorm(BaseDescriptor):
    """
    Fingerprinting calculation for GMP with normalization over orders.

    Args:
        MCSHs [dict] : a dictionary containing the parameters for MCSH definition.

        elements [dict] : a dictionary of string of chemical elements in the system.
    """

    def __init__(
        self,
        MCSHs,
        elements,
        # mode="atom-centered",
    ):
        super().__init__()
        self.descriptor_type = "GMPOrderNorm"
        self.MCSHs = MCSHs
        self.elements = elements
        self.element_indices = list_symbols_to_indices(elements)

        if "cutoff" not in self.MCSHs:
            self.default_cutoff()

        self.prepare_descriptor_parameters()

        self.get_descriptor_setup_hash()

    def __eq__(self, other):
        """Overrides the default implementation"""
        if isinstance(other, BaseDescriptor):
            if self.descriptor_type != other.descriptor_type:
                return False
            if self.descriptor_setup_hash != other.descriptor_setup_hash:
                return False
            return True
        return NotImplemented

    def default_cutoff(self, threshold=1e-8):
        max_sigma = 0.0
        if "MCSHs_detailed_list" in self.MCSHs:
            for detail_setup in self.MCSHs["MCSHs_detailed_list"]:
                sigmas = detail_setup["sigmas"]
                max_sigma = max(max(sigmas), max_sigma)

        else:
            sigmas = self.MCSHs["MCSHs"]["sigmas"]
            max_sigma = max(max(sigmas), max_sigma)

        A = 1.0 / (max_sigma * np.sqrt(2.0 * np.pi))
        alpha = 1.0 / (2.0 * max_sigma * max_sigma)

        cutoff = round(np.sqrt(np.log(threshold / A) / (-alpha)), 2)
        print("default cutoff distance: {} A".format(cutoff))

        self.MCSHs["cutoff"] = cutoff

        return

    def prepare_descriptor_parameters(self):
        descriptor_setup = []
        cutoff = self.MCSHs["cutoff"]
        self.solid_harmonic = self.MCSHs.get("solid_harmonics", True)
        solid_harmonic_i = 1 if self.solid_harmonic else 0
        square = self.MCSHs.get("square", True)
        square_i = 1 if square else 0
        rs_setup = self.MCSHs.get("rs_setup", {"setup": "constant", "rs": 1.0})
        if rs_setup["setup"] == "scale":
            rs_scale = rs_setup["scale_factor"]
        elif rs_setup["setup"] == "constant":
            rs_scale = -1.0 * rs_setup["rs"]
        else:
            rs_scale = -1.0

        if "MCSHs_detailed_list" in self.MCSHs:
            for detail_setup in self.MCSHs["MCSHs_detailed_list"]:
                sigmas = detail_setup["sigmas"]
                descriptor_setup += [
                    [
                        detail_setup["order"],
                        square_i,
                        solid_harmonic_i,
                        sigma,
                        1.0,
                        (1.0 / (sigma * np.sqrt(2.0 * np.pi))) ** 3,
                        1.0 / (2.0 * sigma * sigma),
                        cutoff,
                        (1.0 / (rs_scale * sigma))
                        if rs_scale > 0
                        else (1.0 / (-1.0 * rs_scale)),
                    ]
                    for sigma in sigmas
                ]

        else:
            for i in self.MCSHs["MCSHs"]["orders"]:
                descriptor_setup += [
                    [
                        i,
                        square_i,
                        solid_harmonic_i,
                        sigma,
                        1.0,
                        (1.0 / (sigma * np.sqrt(2.0 * np.pi))) ** 3,
                        1.0 / (2.0 * sigma * sigma),
                        cutoff,
                        (1.0 / (rs_scale * sigma))
                        if rs_scale > 0
                        else (1.0 / (-1.0 * rs_scale)),
                    ]
                    for sigma in self.MCSHs["MCSHs"]["sigmas"]
                ]

        self.descriptor_setup = np.array(descriptor_setup)
        # print(self.descriptor_setup)

        # load pseudo-densities if not provided

        if not self.MCSHs.get("atom_gaussians", None):
            self.load_pseudo_densities(self.elements)

        atomic_gaussian_setup = {}
        for element in self.elements:
            params = list()
            # count = 0
            filename = self.MCSHs["atom_gaussians"][element]
            with open(filename, "r") as fil:
                for line in fil:
                    tmp = line.split()
                    params += [float(tmp[0]), float(tmp[1])]
                    # count += 1
            element_index = ATOM_SYMBOL_TO_INDEX_DICT[element]
            params = np.asarray(params, dtype=np.float64, order="C")
            atomic_gaussian_setup[element_index] = params

        self.atomic_gaussian_setup = atomic_gaussian_setup

        max_gaussian_count = 0
        ngaussian_list = list()
        self.params_set = dict()
        for element_index in self.element_indices:
            self.params_set[element_index] = dict()
            self.params_set[element_index][
                "gaussian_params"
            ] = self.atomic_gaussian_setup[element_index]
            self.params_set[element_index]["gaussian_count"] = int(
                len(self.atomic_gaussian_setup[element_index]) / 2
            )
            ngaussian_list.append(self.params_set[element_index]["gaussian_count"])
            # print("self.params_set[element_index]: {}".format(self.params_set[element_index]))

        # print("ngaussian_list: {}".format(ngaussian_list))
        ngaussian_list = np.asarray(ngaussian_list, dtype=np.intc, order="C")
        max_gaussian_count = np.max(ngaussian_list)
        overall_gaussian_params = list()
        for element_index in self.element_indices:
            temp = np.zeros(max_gaussian_count * 2)
            temp[
                : self.params_set[element_index]["gaussian_count"] * 2
            ] = self.params_set[element_index]["gaussian_params"]
            overall_gaussian_params.append(temp)

        element_index_to_order_list = np.zeros(120, dtype=np.intc)
        for i, element_index in enumerate(self.element_indices):
            # element_index_to_order_list.append(element_index)
            element_index_to_order_list[element_index] = i

        # element_index_to_order_list = np.asarray(element_index_to_order_list, dtype=np.intc, order='C')
        overall_gaussian_params = np.asarray(
            overall_gaussian_params, dtype=np.float64, order="C"
        )
        self.params_set["ngaussians"] = ngaussian_list
        self.params_set["ngaussians_p"] = ffi.cast("int *", ngaussian_list.ctypes.data)
        self.params_set["gaussian_params"] = overall_gaussian_params
        self.params_set["gaussian_params_p"] = _gen_2Darray_for_ffi(
            overall_gaussian_params, ffi
        )
        self.params_set["element_index_to_order"] = element_index_to_order_list
        self.params_set["element_index_to_order_p"] = ffi.cast(
            "int *", element_index_to_order_list.ctypes.data
        )

        # print("ngaussians: {}".format(self.params_set["ngaussians"]))
        # print("gaussian_params: {}".format(self.params_set["gaussian_params"]))

        params_i = np.asarray(
            self.descriptor_setup[:, :3].copy(), dtype=np.intc, order="C"
        )
        params_d = np.asarray(
            self.descriptor_setup[:, 3:].copy(), dtype=np.float64, order="C"
        )
        self.params_set["i"] = params_i
        self.params_set["d"] = params_d

        self.params_set["ip"] = _gen_2Darray_for_ffi(self.params_set["i"], ffi, "int")
        self.params_set["dp"] = _gen_2Darray_for_ffi(self.params_set["d"], ffi)
        self.params_set["total"] = np.concatenate(
            (self.params_set["i"], self.params_set["d"]), axis=1
        )
        self.params_set["num"] = len(self.params_set["total"])

        if "prime_threshold" in self.MCSHs:
            self.params_set["prime_threshold"] = float(self.MCSHs["prime_threshold"])

        self.params_set["log"] = self.MCSHs.get("log", False)

        return

    def load_pseudo_densities(self, elements):
        import os

        # get absolute path to the package
        pkg_path = os.path.dirname(__file__)
        psd_path = os.path.join(pkg_path, "../utils/pseudodensity_psp_v3")

        res = {}
        for element in elements:
            res[element] = os.path.join(psd_path, f"{element}_pseudodensity.g")

        self.MCSHs["atom_gaussians"] = res

    def get_descriptor_setup_hash(self):
        # set self.descriptor_setup_hash
        string = ""
        for desc in self.descriptor_setup:
            for num in desc:
                string += "%.15f" % num
        md5 = hashlib.md5(string.encode("utf-8"))
        hash_result = md5.hexdigest()
        self.descriptor_setup_hash = hash_result

    def save_descriptor_setup(self, filename):
        with open(filename, "w") as out_file:
            for desc in self.descriptor_setup:
                out_file.write(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        int(desc[0]),
                        int(desc[1]),
                        desc[2],
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
        scale = np.copy(atoms.get_scaled_positions(wrap=True), order="C")
        cell = np.copy(atoms.cell, order="C")
        pbc = np.copy(atoms.get_pbc()).astype(np.intc)

        cart_p = _gen_2Darray_for_ffi(cart, ffi)
        scale_p = _gen_2Darray_for_ffi(scale, ffi)
        cell_p = _gen_2Darray_for_ffi(cell, ffi)
        pbc_p = ffi.cast("int *", pbc.ctypes.data)

        cal_atoms = np.asarray(type_idx[element_index], dtype=np.intc, order="C")
        cal_num = len(cal_atoms)
        # print("calculate atom length: {}\ttotal:{}".format(cal_num, atom_num))
        cal_atoms_p = ffi.cast("int *", cal_atoms.ctypes.data)

        # print(self.params_set['i'])
        # print(self.params_set['d'])
        # print(self.params_set['gaussian_params'])
        # print(self.params_set['ngaussians'])
        # print(self.params_set['element_index_to_order'])
        # print(self.params_set['num'])
        # print(atom_indices)

        size_info = np.array([atom_num, cal_num, self.params_set["num"]])

        if calc_derivatives:
            x = np.zeros([cal_num, self.params_set["num"]], dtype=np.float64, order="C")
            dx = np.zeros(
                [cal_num * self.params_set["num"], atom_num * 3],
                dtype=np.float64,
                order="C",
            )

            x_p = _gen_2Darray_for_ffi(x, ffi)
            dx_p = _gen_2Darray_for_ffi(dx, ffi)

            if self.solid_harmonic:
                errno = lib.calculate_solid_gmpordernorm(
                    cell_p,
                    cart_p,
                    scale_p,
                    pbc_p,
                    atom_indices_p,
                    atom_num,
                    cal_atoms_p,
                    cal_num,
                    self.params_set["ip"],
                    self.params_set["dp"],
                    self.params_set["num"],
                    self.params_set["gaussian_params_p"],
                    self.params_set["ngaussians_p"],
                    self.params_set["element_index_to_order_p"],
                    x_p,
                    dx_p,
                )

            else:
                errno = lib.calculate_gmpordernorm(
                    cell_p,
                    cart_p,
                    scale_p,
                    pbc_p,
                    atom_indices_p,
                    atom_num,
                    cal_atoms_p,
                    cal_num,
                    self.params_set["ip"],
                    self.params_set["dp"],
                    self.params_set["num"],
                    self.params_set["gaussian_params_p"],
                    self.params_set["ngaussians_p"],
                    self.params_set["element_index_to_order_p"],
                    x_p,
                    dx_p,
                )

            if errno == 1:
                raise NotImplementedError("Descriptor not implemented!")
            fp = np.array(x, dtype=np.float64)
            fp_prime = np.array(dx, dtype=np.float64)

            fp_prime[np.abs(fp_prime) < 1e-5] = 0.0

            # threshold = 1e-9
            # super_threshold_indices_prime = np.abs(fp_prime) < threshold
            # print("threshold: {} \tnum points set to zero:{} \t outof: {}".format(threshold, np.sum(super_threshold_indices_prime), fp_prime.shape[0] * fp_prime.shape[1]))
            # fp_prime[super_threshold_indices_prime] = 0.0

            scipy_sparse_fp_prime = sparse.coo_matrix(fp_prime)
            # print(fp)
            # print(fp.shape)
            # print(np.sum(super_threshold_indices))
            # print(np.min(np.abs(scipy_sparse_fp_prime.data)))
            # print("density: {}% \n\n----------------------".format(100*len(scipy_sparse_fp_prime.data) / (fp_prime.shape[0] * fp_prime.shape[1])))
            if self.params_set["log"]:
                raise NotImplementedError

            return (
                size_info,
                fp,
                scipy_sparse_fp_prime.data,
                scipy_sparse_fp_prime.row,
                scipy_sparse_fp_prime.col,
                np.array(fp_prime.shape),
            )

        else:
            x = np.zeros([cal_num, self.params_set["num"]], dtype=np.float64, order="C")

            x_p = _gen_2Darray_for_ffi(x, ffi)

            if self.solid_harmonic:
                errno = lib.calculate_solid_gmpordernorm_noderiv(
                    cell_p,
                    cart_p,
                    scale_p,
                    pbc_p,
                    atom_indices_p,
                    atom_num,
                    cal_atoms_p,
                    cal_num,
                    self.params_set["ip"],
                    self.params_set["dp"],
                    self.params_set["num"],
                    self.params_set["gaussian_params_p"],
                    self.params_set["ngaussians_p"],
                    self.params_set["element_index_to_order_p"],
                    x_p,
                )

            else:
                errno = lib.calculate_gmpordernorm_noderiv(
                    cell_p,
                    cart_p,
                    scale_p,
                    pbc_p,
                    atom_indices_p,
                    atom_num,
                    cal_atoms_p,
                    cal_num,
                    self.params_set["ip"],
                    self.params_set["dp"],
                    self.params_set["num"],
                    self.params_set["gaussian_params_p"],
                    self.params_set["ngaussians_p"],
                    self.params_set["element_index_to_order_p"],
                    x_p,
                )

            if errno == 1:
                raise NotImplementedError("Descriptor not implemented!")

            fp = np.array(x, dtype=np.float64)
            if self.params_set["log"]:
                fp = np.abs(fp)
                fp[fp < 1e-8] = 1e-8
                fp = np.log10(fp)

            return size_info, fp, None, None, None, None
