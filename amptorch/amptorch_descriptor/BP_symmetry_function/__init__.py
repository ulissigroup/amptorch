import numpy as np
from scipy import sparse
from ase.calculators.calculator import Parameters
from ._libsymf import lib, ffi
from ..descriptor_base import AMPTorchDescriptorBase
from ..util import _gen_2Darray_for_ffi, list_symbols_to_indices, list_indices_to_symbols
from ..constants import ATOM_INDEX_TO_SYMBOL_DICT, ATOM_SYMBOL_TO_INDEX_DICT


class BPSymmetryFunction(AMPTorchDescriptorBase):

    def __init__(
        self,
        # cutoff=Cosine_cutoff(6.5),
        Gs,
        elements,
        # mode="atom-centered",
    ):
        super().__init__()
        self.descriptor_type = "BP_symmetry_function"
        self.Gs = Gs
        self.elements = elements
        self.element_indices = list_symbols_to_indices(elements)

        self.prepare_descriptor_parameters()

        self.get_descriptor_setup_hash()

    def prepare_descriptor_parameters(self):
        # prepare self.params_set
        descriptor_setup = []
        cutoff = self.Gs["cutoff"]
        # element_indices = list_symbols_to_indices(self.elements)
        if "G2" in self.Gs:
            descriptor_setup += [
                [2, element1, 0, cutoff, eta, rs, 0.0]
                for element1 in self.element_indices
                for eta in self.Gs["G2"]["etas"]
                for rs in self.Gs["G2"]["rs_s"]
            ]
        
        if "G4" in self.Gs:
            descriptor_setup += [
                [4, element1, element2, cutoff, eta, zeta, gamma]
                for element1 in self.element_indices
                for element2 in self.element_indices
                for eta in self.Gs["G4"]["etas"]
                for zeta in self.Gs["G4"]["zetas"]
                for gamma in self.Gs["G4"]["gammas"]
            ]

        if "G5" in self.Gs:
            descriptor_setup += [
                [4, element1, element2, cutoff, eta, zeta, gamma]
                for element1 in self.element_indices
                for element2 in self.element_indices
                for eta in self.Gs["G4"]["etas"]
                for zeta in self.Gs["G4"]["zetas"]
                for gamma in self.Gs["G4"]["gammas"]
            ]

        self.descriptor_setup = np.array(descriptor_setup)

        self.params_set = dict()
        for element_index in self.element_indices:
            self.params_set[element_index] = dict()
            params_i = np.asarray(self.descriptor_setup[:,:3].copy(), dtype=np.intc, order='C')
            params_d = np.asarray(self.descriptor_setup[:,3:].copy(), dtype=np.float64, order='C')
            self.params_set[element_index]["i"] = params_i
            self.params_set[element_index]["d"] = params_d
            self.params_set[element_index]['ip'] = _gen_2Darray_for_ffi(self.params_set[element_index]['i'], ffi, "int")
            self.params_set[element_index]['dp'] = _gen_2Darray_for_ffi(self.params_set[element_index]['d'], ffi)
            self.params_set[element_index]['total'] = np.concatenate((self.params_set[element_index]['i'], self.params_set[element_index]['d']), axis=1)
            self.params_set[element_index]['num'] = len(self.descriptor_setup)

        return
    
    
    def get_descriptor_setup_hash(self):
        #set self.descriptor_setup_hash
        import hashlib
        string = ""
        for desc in self.descriptor_setup:
            for num in desc:
                string += "%.15f" % num
        md5 = hashlib.md5(string.encode("utf-8"))
        hash_result = md5.hexdigest()
        self.descriptor_setup_hash = hash_result

    def save_descriptor_setup(self, filename):
        with open(filename, 'w') as out_file:
            for desc in self.descriptor_setup:
                out_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(int(desc[0]), int(desc[1]), int(desc[2]), desc[3], desc[4], desc[5], desc[6]))


    def calculate_fingerprints(self, atoms, element, log=None, calculate_derivatives=True):
        # params_set = self.params_set
        # print("starting atoms {}".format(atoms) )
        # print("starting element {}".format(element) )
        element_index = ATOM_SYMBOL_TO_INDEX_DICT[element]

        symbols = np.array(atoms.get_chemical_symbols())
        atom_num = len(symbols)
        atom_indices = list_symbols_to_indices(symbols)
        unique_atom_indices = np.unique(atom_indices)

        type_num = dict()
        type_idx = dict()
        
        for atom_index in unique_atom_indices:
            tmp = atom_indices==atom_index
            type_num[atom_index] = np.sum(tmp).astype(np.int64)
            # if atom indexs are sorted by atom type,
            # indexs are sorted in this part.
            # if not, it could generate bug in training process for force training
            type_idx[atom_index] = np.arange(atom_num)[tmp]
        
        # print("type index: {}".format(type_idx))

        atom_indices_p = ffi.cast("int *", atom_indices.ctypes.data)

        cart = np.copy(atoms.get_positions(wrap=True), order='C')
        scale = np.copy(atoms.get_scaled_positions(), order='C')
        cell = np.copy(atoms.cell, order='C')
        pbc = np.copy(atoms.get_pbc()).astype(np.intc)

        cart_p  = _gen_2Darray_for_ffi(cart, ffi)
        scale_p = _gen_2Darray_for_ffi(scale, ffi)
        cell_p  = _gen_2Darray_for_ffi(cell, ffi)
        pbc_p = ffi.cast("int *", pbc.ctypes.data)

        cal_atoms = np.asarray(type_idx[element_index], dtype=np.intc, order='C')
        cal_num = len(cal_atoms)
        cal_atoms_p = ffi.cast("int *", cal_atoms.ctypes.data)

        size_info = np.array([atom_num, cal_num, self.params_set[element_index]['num']])

        if calculate_derivatives:
            x = np.zeros([cal_num, self.params_set[element_index]['num']], dtype=np.float64, order='C')
            dx = np.zeros([cal_num * self.params_set[element_index]['num'], atom_num * 3], dtype=np.float64, order='C')

            x_p = _gen_2Darray_for_ffi(x, ffi)
            dx_p = _gen_2Darray_for_ffi(dx, ffi)

            errno = lib.calculate_sf(cell_p, cart_p, scale_p, pbc_p,\
                            atom_indices_p, atom_num, cal_atoms_p, cal_num, \
                            self.params_set[element_index]['ip'], self.params_set[element_index]['dp'], self.params_set[element_index]['num'], \
                            x_p, dx_p)
            if errno == 1:
                print("ERROR: descriptor not IMPLEMENTED!!")
            fp = np.array(x)
            fp_prime = np.array(dx)
            scipy_sparse_fp_prime = sparse.coo_matrix(fp_prime)
            print("density: {}%".format(100*len(scipy_sparse_fp_prime.data) / (fp_prime.shape[0] * fp_prime.shape[1])))

            return size_info, fp, scipy_sparse_fp_prime.data, scipy_sparse_fp_prime.row, scipy_sparse_fp_prime.col, np.array(fp_prime.shape)
        
        else:
            x = np.zeros([cal_num, self.params_set[element_index]['num']], dtype=np.float64, order='C')
            x_p = _gen_2Darray_for_ffi(x, ffi)

            errno = lib.calculate_sf_no_deriv(cell_p, cart_p, scale_p, pbc_p,\
                            atom_indices_p, atom_num, cal_atoms_p, cal_num, \
                            self.params_set[element_index]['ip'], self.params_set[element_index]['dp'], self.params_set[element_index]['num'], \
                            x_p)
            
            if errno == 1:
                print("ERROR: descriptor not IMPLEMENTED!!")
            fp = np.array(x)

            return size_info, fp, None, None, None, None

