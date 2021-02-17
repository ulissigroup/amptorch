import hashlib

import numpy as np
from scipy import sparse

from ..base_descriptor import BaseDescriptor
from ..constants import ATOM_SYMBOL_TO_INDEX_DICT
from ..util import _gen_2Darray_for_ffi, list_symbols_to_indices
from ._libgmpalign import ffi, lib


class GMPAlign(BaseDescriptor):
    def __init__(
        self,
        MCSHs,
        elements,
        # mode="atom-centered",
    ):
        super().__init__()
        self.descriptor_type = "GMP_Align"
        self.MCSHs = MCSHs
        self.elements = elements
        self.element_indices = list_symbols_to_indices(elements)

        self.prepare_descriptor_parameters()

        self.get_descriptor_setup_hash()

    def prepare_descriptor_parameters(self):
        descriptor_setup = []
        cutoff = self.MCSHs["cutoff"]

        for i in range(10):
            if str(i) in self.MCSHs["MCSHs"].keys():
                descriptor_setup += [
                    [
                        i,
                        group,
                        sigma,
                        1.0,
                        1 / (sigma * np.sqrt(2 * np.pi)),
                        1 / (2 * sigma * sigma),
                        cutoff,
                    ]
                    for group in self.MCSHs["MCSHs"][str(i)]["groups"]
                    for sigma in self.MCSHs["MCSHs"][str(i)]["sigmas"].sort(reverse = True)
                ]

        total_num_desc = 0
        for setup in descriptor_setup:
            total_num_desc += self.group_get_num_desc(setup[0], setup[1])

        self.descriptor_setup = np.array(descriptor_setup)

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
            self.descriptor_setup[:, :2].copy(), dtype=np.intc, order="C"
        )
        params_d = np.asarray(
            self.descriptor_setup[:, 2:].copy(), dtype=np.float64, order="C"
        )
        self.params_set["i"] = params_i
        self.params_set["d"] = params_d

        self.params_set["ip"] = _gen_2Darray_for_ffi(self.params_set["i"], ffi, "int")
        self.params_set["dp"] = _gen_2Darray_for_ffi(self.params_set["d"], ffi)
        self.params_set["total"] = np.concatenate(
            (self.params_set["i"], self.params_set["d"]), axis=1
        )
        self.params_set["num"] = len(self.params_set["total"])

        self.params_set["total_num_desc"] = total_num_desc

        if "prime_threshold" in self.MCSHs:
            self.params_set["prime_threshold"] = float(self.MCSHs["prime_threshold"])

        self.params_set["align"] = self.MCSHs.get("align", True)

        return

    def group_get_num_desc(self, order, group):
        num_desc_dict = {
            0: {1: 1},
            1: {1: 3},
            2: {1: 3, 2: 3},
            3: {1: 3, 2: 6, 3: 1},
            4: {1: 3, 2: 6, 3: 3, 4: 3},

            5: {1: 3, 2: 6, 3: 6, 4: 3, 5: 3},
            6: {1: 3, 2: 6, 3: 6, 4: 3, 5: 3, 6: 6, 7: 1},
            7: {1: 3, 2: 6, 3: 6, 4: 3, 5: 6, 6: 6, 7: 3, 8: 3},
            8: {1: 3, 2: 6, 3: 6, 4: 3, 5: 6, 6: 6, 7: 3, 8: 6, 9: 3, 10: 3},
            9: {1: 3, 2: 6, 3: 6, 4: 3, 5: 6, 6: 6, 7: 6, 8: 6, 9: 3, 10: 3, 11: 6, 12: 1},
        }
        return num_desc_dict[order][group]

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
        scale = np.copy(atoms.get_scaled_positions(), order="C")
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
            raise NotImplementedError

        else:
            x = np.zeros([cal_num, self.params_set["total_num_desc"]], dtype=np.float64, order="C")

            x_p = _gen_2Darray_for_ffi(x, ffi)

            if self.params_set["align"]:
                errno = lib.calculate_gmp_align_noderiv(
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
                    self.params_set["total_num_desc"],
                    self.params_set["gaussian_params_p"],
                    self.params_set["ngaussians_p"],
                    self.params_set["element_index_to_order_p"],
                    x_p,
                )
            else:
                errno = lib.calculate_gmp_no_align_noderiv(
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
                    self.params_set["total_num_desc"],
                    self.params_set["gaussian_params_p"],
                    self.params_set["ngaussians_p"],
                    self.params_set["element_index_to_order_p"],
                    x_p,
                )

            if errno == 1:
                raise NotImplementedError("Descriptor not implemented!")

            fp = np.array(x)

            return size_info, fp, None, None, None, None
