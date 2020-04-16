from __future__ import print_function
from __future__ import division
import os, sys
import tensorflow as tf
import numpy as np
import six
from six.moves import cPickle as pickle
from ase import io
from simple_nn.features.symmetry_function import Symmetry_function
from simple_nn.features.symmetry_function._libsymf import lib, ffi
from simple_nn.utils import _gen_2Darray_for_ffi, compress_outcar, _generate_scale_file, \
    _make_full_featurelist, _make_data_list, _make_str_data_list, pickle_load
from simple_nn.utils import graph as grp
from simple_nn.utils.mpiclass import DummyMPI, MPI4PY
from braceexpand import braceexpand


def _read_params(filename, dir):
    params_i = list()
    params_d = list()
    filename = dir + "/" + filename
    with open(filename, 'r') as fil:
        for line in fil:
            tmp = line.split()
            params_i += [list(map(int, tmp[:3]))]
            params_d += [list(map(float, tmp[3:]))]

    params_i = np.asarray(params_i, dtype=np.intc, order='C')
    params_d = np.asarray(params_d, dtype=np.float64, order='C')

    return params_i, params_d

class Symmetry_function_customized(Symmetry_function):
    def __init__(self, fp_dir, inputs=None, structure_list_path='/str_list'):
        super().__init__(fp_dir, inputs, structure_list_path)
        self.parent = None
        self.dir = fp_dir
        self.key = 'symmetry_function'
        self.default_inputs = {'symmetry_function':
            {
                'params': dict(),
                'refdata_format': 'traj',  # Change default to traj
                'compress_outcar': True,
                'data_per_tfrecord': 150,
                'valid_rate': 0.1,
                'remain_pickle': False,
                'continue': False,
                'add_atom_idx': True,  # For backward compatability
                'num_parallel_calls': 5,
                'atomic_weights': {
                    'type': None,
                    'params': dict(),
                },
                'weight_modifier': {
                    'type': None,
                    'params': dict(),
                },
                'scale_type': 'minmax',
                'scale_scale': 1.0,
            }
        }
        # self.structure_list = './str_list'
        self.structure_list = fp_dir + structure_list_path
        self.pickle_list = fp_dir + '/pickle_list'
        self.train_data_list = fp_dir + '/train_list'
        self.valid_data_list = fp_dir + '/valid_list'
        self.comm = None

    def generate(self, label):

        comm = self.get_comm()

        if comm.rank == 0:
            train_dir = open(self.pickle_list, 'w')

        # Get structure list to calculate
        structures, structure_ind, structure_names, structure_weights = self._parse_strlist()

        # Get parameter list for each atom types
        params_set = dict()
        for item in self.parent.inputs['atom_types']:
            params_set[item] = dict()
            params_set[item]['i'], params_set[item]['d'] = \
                _read_params(self.inputs['params'][item], self.dir)
            params_set[item]['ip'] = _gen_2Darray_for_ffi(params_set[item]['i'], ffi, "int")
            params_set[item]['dp'] = _gen_2Darray_for_ffi(params_set[item]['d'], ffi)
            params_set[item]['total'] = np.concatenate((params_set[item]['i'], params_set[item]['d']), axis=1)
            params_set[item]['num'] = len(params_set[item]['total'])

        data_idx = 1
        for item, ind in zip(structures, structure_ind):
            # FIXME: add another input type

            if len(item) == 1:
                index = 0
                if comm.rank == 0:
                    self.parent.logfile.write('{} 0'.format(item[0]))
            else:
                if ':' in item[1]:
                    index = item[1]
                else:
                    index = int(item[1])

                if comm.rank == 0:
                    self.parent.logfile.write('{} {}'.format(item[0], item[1]))

            if self.inputs['refdata_format'] == 'vasp-out':
                if self.inputs['compress_outcar']:
                    if comm.rank == 0:
                        tmp_name = compress_outcar(item[0])
                    else:
                        tmp_name = None
                    tmp_name = comm.bcast(tmp_name, root=0)
                    comm.barrier()
                    snapshots = io.read(tmp_name, index=index, format=self.inputs['refdata_format'],
                                        force_consistent=True)
                    comm.barrier()
                    if comm.rank == 0:
                        os.remove(tmp_name)
                else:
                    snapshots = io.read(item[0], index=index, format=self.inputs['refdata_format'],
                                        force_consistent=True)
            else:
                snapshots = io.read(item[0], index=index, format=self.inputs['refdata_format'])

            for atoms in snapshots:
                ad_atoms = self.Get_adsorbat_atoms(atoms)
                cart = np.copy(atoms.get_positions(wrap=True), order='C')
                scale = np.copy(atoms.get_scaled_positions(), order='C')
                cell = np.copy(atoms.cell, order='C')

                cart_p = _gen_2Darray_for_ffi(cart, ffi)
                scale_p = _gen_2Darray_for_ffi(scale, ffi)
                cell_p = _gen_2Darray_for_ffi(cell, ffi)

                atom_num = len(atoms.positions)
                symbols = np.array(atoms.get_chemical_symbols())
                atom_i = np.zeros([len(symbols)], dtype=np.intc, order='C')
                type_num = dict()
                type_idx = dict()
                for j, jtem in enumerate(self.parent.inputs['atom_types']):
                    tmp = symbols == jtem
                    atom_i[tmp] = j + 1
                    type_num[jtem] = np.sum(tmp).astype(np.int64)
                    # if atom indexs are sorted by atom type,
                    # indexs are sorted in this part.
                    # if not, it could generate bug in training process for force training
                    type_idx[jtem] = np.arange(atom_num)[tmp]
                atom_i_p = ffi.cast("int *", atom_i.ctypes.data)

                res = dict()
                res['x'] = dict()
                res['dx'] = dict()
                res['params'] = dict()
                res['N'] = type_num
                res['tot_num'] = np.sum(list(type_num.values()))
                res['partition'] = np.ones([res['tot_num']]).astype(np.int32)
                res['E'] = None
                res['F'] = None
                # res['E'] = atoms.get_total_energy()
                # res['F'] = atoms.get_forces()
                res['struct_type'] = structure_names[ind]
                res['struct_weight'] = structure_weights[ind]
                res['atom_idx'] = atom_i

                cal_type_num, cal_type_idx = self.Get_calc_type(ad_atoms, type_idx, self.parent.inputs['atom_types'])

                for j, jtem in enumerate(self.parent.inputs['atom_types']):
                    q = cal_type_num[jtem] // comm.size
                    r = cal_type_num[jtem] % comm.size

                    begin = comm.rank * q + min(comm.rank, r)
                    end = begin + q
                    if r > comm.rank:
                        end += 1

                    cal_atoms = np.asarray(cal_type_idx[jtem][begin:end], dtype=np.intc, order='C')
                    cal_num = len(cal_atoms)
                    cal_atoms_p = ffi.cast("int *", cal_atoms.ctypes.data)

                    x = np.zeros([cal_num, params_set[jtem]['num']], dtype=np.float64, order='C')
                    dx = np.zeros([cal_num, atom_num * params_set[jtem]['num'] * 3], dtype=np.float64, order='C')

                    x_p = _gen_2Darray_for_ffi(x, ffi)
                    dx_p = _gen_2Darray_for_ffi(dx, ffi)

                    errno = lib.calculate_sf(cell_p, cart_p, scale_p, \
                                             atom_i_p, atom_num, cal_atoms_p, cal_num, \
                                             params_set[jtem]['ip'], params_set[jtem]['dp'], params_set[jtem]['num'], \
                                             x_p, dx_p)
                    comm.barrier()
                    if errno == 1:
                        err = "Not implemented symmetry function type."
                        self.parent.logfile.write("Error: {:}\n".format(err))
                        raise NotImplementedError(err)
                    elif errno == 2:
                        err = "Zeta in G4/G5 must be greater or equal to 1.0."
                        self.parent.logfile.write("Error: {:}\n".format(err))
                        raise ValueError(err)
                    else:
                        assert errno == 0

                    if cal_type_num[jtem] != 0:
                        res['x'][jtem] = np.array(comm.gather(x, root=0))
                        res['dx'][jtem] = np.array(comm.gather(dx, root=0))
                        if comm.rank == 0:
                            res['x'][jtem] = np.concatenate(res['x'][jtem], axis=0).reshape(
                                [cal_type_num[jtem], params_set[jtem]['num']])
                            res['dx'][jtem] = np.concatenate(res['dx'][jtem], axis=0). \
                                reshape([cal_type_num[jtem], params_set[jtem]['num'], atom_num, 3])
                            res['partition_' + jtem] = np.ones([cal_type_num[jtem]]).astype(np.int32)
                    else:
                        res['x'][jtem] = np.zeros([0, params_set[jtem]['num']])
                        res['dx'][jtem] = np.zeros([0, params_set[jtem]['num'], atom_num, 3])
                        res['partition_' + jtem] = np.ones([0]).astype(np.int32)
                    res['params'][jtem] = params_set[jtem]['total']

                if comm.rank == 0:
                    data_dir = "./data" + label + "/"
                    if not os.path.exists(data_dir):
                        os.makedirs(data_dir)
                    tmp_filename = os.path.join(data_dir, "data{}.pickle".format(data_idx))
                    # tmp_filename = os.path.join(data_dir, "data{}.tfrecord".format(data_idx))

                    # TODO: add tfrecord writing part
                    # self._write_tfrecords(res, tmp_filename)
                    with open(tmp_filename, "wb") as fil:
                        pickle.dump(res, fil, protocol=2)

                    train_dir.write('{}:{}\n'.format(ind, tmp_filename))
                    tmp_endfile = tmp_filename
                    data_idx += 1

            if comm.rank == 0:
                self.parent.logfile.write(': ~{}\n'.format(tmp_endfile))

        if comm.rank == 0:
            train_dir.close()

    def Get_adsorbat_atoms(self, atoms):
        adsorbate_atoms = []  # the index of adsorbate atoms
        for index, atom in enumerate(atoms):
            if atom.tag == 1:
                adsorbate_atoms.append(index)
        return adsorbate_atoms

    def Get_calc_type(self, ad_atoms, type_idx, atom_types):
        cal_type_num = dict()
        cal_type_idx = dict()
        for j, jitem in enumerate(atom_types):
            num = 0
            idx = []
            for idxi in type_idx[jitem]:
                if idxi in ad_atoms:
                    num = num + 1
                    idx.append(idxi)
            cal_type_num[jitem] = num
            cal_type_idx[jitem] = np.array(idx)
        return cal_type_num, cal_type_idx



