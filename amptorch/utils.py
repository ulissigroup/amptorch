import sys
import copy
import time
import os
import pickle
from pickle import load
from collections import defaultdict, OrderedDict
import shutil
import numpy as np
from ase import io
from ase.db import connect
from simple_nn.features.symmetry_function import Symmetry_function
from .simple_nn_specific_atoms import Symmetry_function_customized

def hash_images(images, Gs=None, log=None, ordered=False):
    """ Converts input images -- which may be a list, a trajectory file, or
    a database -- into a dictionary indexed by their hashes.

    Returns this dictionary. If ordered is True, returns an OrderedDict. When
    duplicate images are encountered (based on encountering an identical hash),
    a warning is written to the logfile. The number of duplicates of each image
    can be accessed by examinging dict_images.metadata['duplicates'], where
    dict_images is the returned dictionary.
    """
    if log is None:
        log = Logger(None)
    if images is None:
        return
    elif hasattr(images, "keys"):
        log(" %i unique images after hashing." % len(images))
        return images  # Apparently already hashed.
    else:
        # Need to be hashed, and possibly read from file.
        if isinstance(images, str):
            log("Attempting to read images from file %s." % images)
            extension = os.path.splitext(images)[1]
            from ase import io

            if extension == ".traj":
                images = io.Trajectory(images, "r")
            elif extension == ".db":
                images = [row.toatoms() for row in connect(images, "db").select(None)]

        # images converted to dictionary form; key is hash of image.
        log("Hashing images...", tic="hash")
        dict_images = MetaDict()
        dict_images.metadata["duplicates"] = {}
        dup = dict_images.metadata["duplicates"]
        if ordered is True:
            dict_images = OrderedDict()
        for image in images:
            hash = get_hash(image, Gs)
            if hash in dict_images.keys():
                log(
                    "Warning: Duplicate image (based on identical hash)."
                    " Was this expected? Hash: %s" % hash
                )
                if hash in dup.keys():
                    dup[hash] += 1
                else:
                    dup[hash] = 2
            dict_images[hash] = image
        log(" %i unique images after hashing." % len(dict_images))
        log("...hashing completed.", toc="hash")
        return dict_images


def calculate_fingerprints_range(fp, images):
    """Calculates the range for the fingerprints corresponding to images,
    stored in fp. fp is a fingerprints object with the fingerprints data
    stored in a dictionary-like object at fp.fingerprints. (Typically this
    is a .utilties.Data structure.) images is a hashed dictionary of atoms
    for which to consider the range.

    In image-centered mode, returns an array of (min, max) values for each
    fingerprint. In atom-centered mode, returns a dictionary of such
    arrays, one per element.
    """
    if fp.parameters.mode == "image-centered":
        raise NotImplementedError()
    elif fp.parameters.mode == "atom-centered":
        fprange = {}
        for hash in images.keys():
            imagefingerprints = fp.fingerprints[hash]
            for element, fingerprint in imagefingerprints:
                if element not in fprange:
                    fprange[element] = [[_, _] for _ in fingerprint]
                else:
                    assert len(fprange[element]) == len(fingerprint)
                    for i, ridge in enumerate(fingerprint):
                        if ridge < fprange[element][i][0]:
                            fprange[element][i][0] = ridge
                        elif ridge > fprange[element][i][1]:
                            fprange[element][i][1] = ridge
    for key, value in fprange.items():
        fprange[key] = value
    return fprange


def make_params_file(
    elements, fp_dir, etas, rs_s, g4_eta=4, cutoff=6.5, g4_zeta=[1.0, 4.0], g4_gamma=[1, -1]
    ):
    """
    makes a params file for simple_NN. This is the file containing
    the descriptors. This function makes g2 descriptos for the eta
    and rs values that are input, and g4 descriptors that are log
    spaced between 10 ** -5 and 10 ** -1. The number of these
    that are made is controlled by the `n_g4_eta` variable
    Parameters:
        elements (list):
            a list of elements for which you'd like to make params
            files for
        etas (list):
            the eta values you'd like to use for the descriptors
        rs_s (list):
            a list corresponding to `etas` that contains the rs
            values for each descriptor
        g4_eta (int or list):
            the number of g4 descriptors you'd like to use. if a
            list is passed in the values of the list will be used
            as eta values
        cutoff (float):
            the distance in angstroms at which you'd like to cut
            off the descriptors
    returns:
        None
    """
    if len(etas) != len(rs_s):
        raise ValueError('the length of the etas list must be equal to the'
                         'length of the rs_s list')
    if type(g4_eta) == int:
        g4_eta = np.logspace(-4, -1, num=g4_eta)
    for element in elements:
        with open(fp_dir+"/params_{}".format(element), "w") as f:
            # G2
            for eta, Rs in zip(etas, rs_s):
                for species in range(1, len(elements) + 1):
                    f.write(
                        "2 {} 0 {} {} {} 0.0\n".format(
                            #species, cutoff, np.round(eta, 6), Rs
                            species, cutoff, eta, Rs
                        )
                    )

            # G4
            for eta in g4_eta:
                for zeta in g4_zeta:
                    for lamda in g4_gamma:
                        for i in range(1, len(elements) + 1):
                            for j in range(i, len(elements) + 1):
                                f.write(
                                    "4 {} {} {} {} {} {}\n".format(
                                        i, j, cutoff, eta, zeta, lamda
                                    )
                                )


def reorganize_simple_nn_derivative(image, specific_atoms, dx_dict):
    """
    reorganizes the fingerprint derivatives from simplen_nn into
    amp format
    Parameters:
        image (ASE atoms object):
            the atoms object used to make the finerprint
        dx_dict (dict):
            a dictionary of the fingerprint derivatives from simple_nn
    """
    # TODO check for bugs
    d = OrderedDict()
    sym_dict = OrderedDict()
    if specific_atoms is True:
        syms = image.get_chemical_symbols()
        for sym in syms:
            sym_dict[sym] = []
        for i, sym in enumerate(syms):
            if image[i].tag == 1:
                sym_dict[sym].append(i)
    else:
        syms = image.get_chemical_symbols()
        for sym in syms:
            sym_dict[sym] = []
        for i, sym in enumerate(syms):
            sym_dict[sym].append(i)
    # the structure is:
    # [elements][atom i][symetry function #][atom j][derivitive in direction]
    for element, full_arr in dx_dict.items():
        for i, arr_t in enumerate(full_arr):
            true_i = sym_dict[element][i]
            for sf in arr_t:
                for j, dir_arr in enumerate(sf):
                    for k, derivative in enumerate(dir_arr):
                        if (j, syms[j], true_i, element, k) not in d:
                            d[(j, syms[j], true_i, element, k)] = []
                        d[(j, syms[j], true_i, element, k)].append(derivative)
    # zero_keys = []
    # for key, derivatives in d.items():
        # zero_check = [a == 0 for a in derivatives]
        # if zero_check == [True] * len(derivatives):
            # zero_keys.append(key)
    # for key in zero_keys:
        # del d[key]
    # d = OrderedDict(d)
    return d


def reorganize_simple_nn_fp(image, specific_atoms, x_dict):
    """
    reorganizes the fingerprints from simplen_nn into
    amp format
    Parameters:
        image (ASE atoms object):
            the atoms object used to make the finerprint
        x_dict (dict):
            a dictionary of the fingerprints from simple_nn
    """
    # TODO check for bugs
    # the structure is:
    # [elements][atom i][symetry function #][fp]
    fp_l = []
    sym_dict = OrderedDict()
    if specific_atoms is True:
        syms = image.get_chemical_symbols()
        for sym in syms:
            sym_dict[sym] = []
        for i, sym in enumerate(syms):
            if image[i].tag == 1:
                sym_dict[sym].append(i)
        for i, sym in enumerate(syms):
            if image[i].tag == 1:
                simple_nn_index = sym_dict[sym].index(i)
                fp = x_dict[sym][simple_nn_index]
                fp_l.append((sym, list(fp)))
    else:
        syms = image.get_chemical_symbols()
        for sym in syms:
            sym_dict[sym] = []
        for i, sym in enumerate(syms):
            sym_dict[sym].append(i)
        for i, sym in enumerate(syms):
            simple_nn_index = sym_dict[sym].index(i)
            fp = x_dict[sym][simple_nn_index]
            fp_l.append((sym, list(fp)))
    return fp_l


def get_hash(atoms, Gs=None):
    import hashlib

    """Creates a unique signature for a particular ASE atoms object.
    This is used to check whether an image has been seen before. This is just
    an md5 hash of a string representation of the atoms object and symmetry
    functions.
    Parameters
    ----------
    atoms : ASE dict
        ASE atoms object.
    Returns
    -------
        Hash string key of 'atoms'.
    """
    string = str(atoms.pbc)
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
    if Gs:
        gs_values = list(Gs.values())
        for number in gs_values[0]:
            string += "%.15f" % number
        for number in gs_values[1]:
            string += "%.15f" % number
        for number in gs_values[2]:
            string += "%.15f" % number
        for number in gs_values[3]:
            string += "%.15f" % number
        for number in gs_values[4]:
            string += "%.15f" % number
        string += "%.15f" % gs_values[5]

    md5 = hashlib.md5(string.encode("utf-8"))
    hash = md5.hexdigest()
    return hash


def factorize_data(traj, Gs):
    new_traj = []
    if os.path.isdir("amp-data-fingerprint-primes.ampdb/"):
        for image in traj:
            hash = get_hash(image, Gs)
            if os.path.isfile(
                "amp-data-fingerprint-primes.ampdb/loose/" + hash
            ) and os.path.isfile("amp-data-fingerprints.ampdb/loose/" + hash):
                pass
            else:
                new_traj.append(image)
    else:
        new_traj = traj
    return new_traj


def convert_simple_nn_fps(traj, Gs, specific_atoms, cores, label, save, delete_old=True):
    from multiprocessing import Pool

    # make the directories
    if not os.path.isdir("./amp-data-fingerprints.ampdb"):
        os.mkdir("./amp-data-fingerprints.ampdb")
    if not os.path.isdir("./amp-data-fingerprints.ampdb/loose"):
        os.mkdir("./amp-data-fingerprints.ampdb/loose")
    if not os.path.isdir("./amp-data-fingerprint-primes.ampdb"):
        os.mkdir("./amp-data-fingerprint-primes.ampdb")
    if not os.path.isdir("./amp-data-fingerprint-primes.ampdb/loose"):
        os.mkdir("amp-data-fingerprint-primes.ampdb/loose")
    # perform the reorganization
    l_trajs = list(enumerate(traj))
    fp_dir = "data"+label
    if len(traj) > 1:
        with Pool(cores) as p:
            l_trajs = [image + (Gs, specific_atoms, label, fp_dir) for image in l_trajs]
            p.map(reorganize, l_trajs)
    else:
        image = (0, traj[0], Gs, specific_atoms, label, fp_dir)
        fps, fp_primes = reorganize(image, save=save)
    if delete_old:
        os.rmdir("./data"+label)
    if save:
        return None, None
    else:
        return fps, fp_primes


def reorganize(inp, delete_old=True, save=True):
    i, image, Gs, specific_atoms, label, fp_dir = inp
    pic = pickle.load(open(fp_dir+"/data{}.pickle".format(i + 1), "rb"))
    im_hash = get_hash(image, Gs)
    x_list = reorganize_simple_nn_fp(image, specific_atoms, pic["x"])
    x_der_dict = reorganize_simple_nn_derivative(image, specific_atoms, pic["dx"])
    if save:
        pickle.dump(x_list, open("./amp-data-fingerprints.ampdb/loose/" + im_hash, "wb"))
        pickle.dump(
            x_der_dict, open("./amp-data-fingerprint-primes.ampdb/loose/" + im_hash, "wb")
        )
    if delete_old:  # in case disk space is an issue
        if os.path.exists(os.path.join("./data" + label,'simple_nn_log')):
            os.remove(os.path.join("./data" + label,'simple_nn_log'))
        os.remove(fp_dir+"/data{}.pickle".format(i + 1))
    return x_list, x_der_dict


class DummySimple_nn(object):
    """
    a dummy class to fool the simple_nn descriptor class into
    thinking it's attached to a simple_nn instance
    """

    def __init__(self, atom_types, dir):
        self.inputs = {
            "generate_features": True,
            "preprocess": False,
            "train_model": True,
            "atom_types": atom_types,
        }
        self.logfile = open(dir+"/simple_nn_log", "w")

def make_simple_nn_fps(traj, Gs, specific_atoms, label, clean_up_directory=True, elements="all"):
    """
    generates descriptors using simple_nn. The files are stored in the
    ./data folder. These descriptors will be in the simple_nn form and
    not immediately useful for other programs
    Parameters:
        traj (list of ASE atoms objects):
            a list of the atoms you'd like to make descriptors for
        descriptors (tuple):
            a tuple containing (g2_etas, g2_rs_s, g4_etas, cutoff, g4_zetas, g4_gammas)
        clean_up_directory (bool):
            if set to True, the input files made by simple_nn will
            be deleted
    returns:
        None
    """
    # handle inputs
    if type(traj) != list:
        traj = [traj]

    G = copy.deepcopy(Gs)
    traj = factorize_data(traj, G)
    calculated = False
    if len(traj) > 0:
        # order descriptors for simple_nn
        cutoff = G["cutoff"]
        G["G2_etas"] = [a / cutoff**2 for a in G["G2_etas"]]
        G["G4_etas"] = [a / cutoff**2 for a in G["G4_etas"]]
        descriptors = (
            G["G2_etas"],
            G["G2_rs_s"],
            G["G4_etas"],
            G["cutoff"],
            G["G4_zetas"],
            G["G4_gammas"],
        )
        fp_dir = "./data"+label+"/"
        if not os.path.isdir(fp_dir):
            os.mkdir(fp_dir)
        # clean up any previous runs
        if os.path.isdir(fp_dir+"/data/"):
            shutil.rmtree(fp_dir+"/data/")

        # set up the input files
        io.write(fp_dir+"/simple_nn_input_traj.traj", traj)
        with open(fp_dir+"/str_list", "w") as f:
            # simple_nn requires this file
            f.write(fp_dir+"/simple_nn_input_traj.traj :")

        if elements == "all":
            atom_types = []
            # TODO rewrite this
            for image in traj:
                atom_types += image.get_chemical_symbols()
                atom_types = list(set(atom_types))
        else:
            atom_types = elements

        make_params_file(atom_types, fp_dir, *descriptors)

        # build the descriptor object
        if specific_atoms is True:
            descriptor = Symmetry_function_customized(fp_dir=fp_dir)
        else:
            descriptor = Symmetry_function(fp_dir=fp_dir)
        params = {a: "params_{}".format(a) for a in atom_types}

        descriptor.inputs = {
            "params": params,
            "refdata_format": "traj",
            "compress_outcar": False,
            "data_per_tfrecord": 150,
            "valid_rate": 0.1,
            "remain_pickle": False,
            "continue": False,
            "add_atom_idx": True,
            "num_parallel_calls": 5,
            "atomic_weights": {"type": None, "params": {}},
            "weight_modifier": {"type": None, "params": {}},
            "scale_type": 'minmax',
            "scale_scale": 1.0,
            "scale_rho": None,
        }
        dummy_class = DummySimple_nn(atom_types=atom_types, dir=fp_dir)
        descriptor.parent = dummy_class

        # generate the descriptors
        descriptor.generate(label)
        if clean_up_directory:
            # clean the folder of all the junk
            files = [
                "simple_nn_input_traj.traj",
                "str_list",
                "pickle_list",
                "simple_nn_log",
            ]
            files += list(params.values())
            dummy_class.logfile.close()
            for file in files:
                os.remove(fp_dir+"/"+file)
        calculated = True
    return traj, calculated

def stored_fps(traj, Gs):
    image_hash = get_hash(traj[0], Gs)
    with open("amp-data-fingerprints.ampdb/loose/"+image_hash, "rb") as f:
        fps = load(f)
    with open("amp-data-fingerprint-primes.ampdb/loose/"+image_hash, "rb") as f:
        fp_primes = load(f)
    return fps, fp_primes

def make_amp_descriptors_simple_nn(atoms, Gs, elements, cores, label, save=True, specific_atoms=False):
    """
    uses simple_nn to make descriptors in the amp format.
    Only creates the same symmetry functions for each element
    for now.
    """
    traj, calculated = make_simple_nn_fps(atoms, Gs, specific_atoms, elements=elements,
            label=label, clean_up_directory=True)
    if calculated:
        fps, fp_primes = convert_simple_nn_fps(traj, Gs, specific_atoms, cores, label, save, delete_old=True)
        return fps, fp_primes
    if save is False and calculated is False:
        fps, fp_primes = stored_fps(atoms, Gs)
        return fps, fp_primes
    else: return None, None


class Logger:
    """Logger that can also deliver timing information.
    Parameters
    ----------
    file : str
        File object or path to the file to write to.  Or set to None for
        a logger that does nothing.
    """

    def __init__(self, file):
        if file is None:
            self.file = None
            return
        if isinstance(file, str):
            self.filename = file
            file = open(file, "a")
        self.file = file
        self.tics = {}

    def tic(self, label=None):
        """Start a timer.

        Parameters
        ----------
        label : str
            Label for managing multiple timers.
        """
        if self.file is None:
            return
        if label:
            self.tics[label] = time.time()
        else:
            self._tic = time.time()

    def __call__(self, message, toc=None, tic=False):
        """Writes message to the log file.

        Parameters
        ---------
        message : str
            Message to be written.
        toc : bool or str
            If toc=True or toc=label, it will append timing information in
            minutes to the timer.
        tic : bool or str
            If tic=True or tic=label, will start the generic timer or a timer
            associated with label. Equivalent to self.tic(label).
        """
        if self.file is None:
            return
        dt = ""
        if toc:
            if toc is True:
                tic = self._tic
            else:
                tic = self.tics[toc]
            dt = (time.time() - tic) / 60.0
            dt = " %.1f min." % dt
        if self.file.closed:
            self.file = open(self.filename, "a")
        self.file.write(message + dt + "\n")
        self.file.flush()
        if tic:
            if tic is True:
                self.tic()
            else:
                self.tic(label=tic)


class MetaDict(dict):
    """Dictionary that can also store metadata. Useful for iamges dictionary
    so that images can still be iterated by keys.
    """

    metadata = {}


def make_force_header(log):
    header = "%5s %24s %12s %12s %12s"
    log(header % ("Epoch", "Time", "Loss", "EnergyRMSE", "ForceRMSE"))
    log(header % ("=" * 5, "=" * 24, "=" * 12, "=" * 12, "=" * 12))


def make_energy_header(log):
    header = "%5s %24s %12s %7s"
    log(header % ("Epoch", "Time", "Loss", "EnergyRMSE"))
    log(header % ("=" * 5, "=" * 24, "=" * 12, "=" * 12))


def make_val_force_header(log):
    header = "%5s %24s %12s %12s %12s %7s"
    log(header % ("Epoch", "Time", "Loss", "EnergyRMSE", "ForceRMSE", "Phase"))
    log(header % ("=" * 5, "=" * 24, "=" * 12, "=" * 12, "=" * 12, "=" * 7))


def make_val_energy_header(log):
    header = "%5s %24s %12s %12s %7s"
    log(header % ("Epoch", "Time", "Loss", "EnergyRMSE", "Phase"))
    log(header % ("=" * 5, "=" * 24, "=" * 12, "=" * 12, "=" * 7))


def log_force_results(log, epoch, now, loss, energy_rmse, force_rmse, phase):
    if type(loss) is str:
        log(
            "%5i %19s %12s %12.4e %12.4e %7s"
            % (epoch, now, loss, energy_rmse, force_rmse, phase)
        )
    else:
        log(
                   "%5i %19s %12.4e %12.4e %12.4e %7s"
            % (epoch, now, loss, energy_rmse, force_rmse, phase)
                )
def log_energy_results(log, epoch, now, loss, energy_rmse, phase):
    if type(loss) is str:
        log("%5i %19s %12s %12.4e %7s" % (epoch, now, loss, energy_rmse, phase))
    else:
        log("%5i %19s %12.4e %12.4e %7s" % (epoch, now, loss, energy_rmse, phase))

def dict2cutoff(dct):
    """This function converts a dictionary (which was created with the
    to_dict method of one of the cutoff classes) into an instantiated
    version of the class. Modeled after ASE's dict2constraint function.
    """
    if len(dct) != 2:
        raise RuntimeError('Cutoff dictionary must have only two values,'
                           ' "name" and "kwargs".')
    return globals()[dct['name']](**dct['kwargs'])

class Cosine(object):
    """Cosine functional form suggested by Behler.

    Parameters
    ---------
    Rc : float
        Radius above which neighbor interactions are ignored.
    """

    def __init__(self, Rc):

        self.Rc = Rc

    def __call__(self, Rij):
        """
        Parameters
        ----------
        Rij : float
            Distance between pair atoms.

        Returns
        -------
        float
            The value of the cutoff function.
        """
        if Rij > self.Rc:
            return 0.
        else:
            return 0.5 * (np.cos(np.pi * Rij / self.Rc) + 1.)

    def prime(self, Rij):
        """Derivative (dfc_dRij) of the Cosine cutoff function with respect to Rij.

        Parameters
        ----------
        Rij : float
            Distance between pair atoms.

        Returns
        -------
        float
            The value of derivative of the cutoff function.
        """
        if Rij > self.Rc:
            return 0.
        else:
            return -0.5 * np.pi / self.Rc * np.sin(np.pi * Rij / self.Rc)

    def todict(self):
        return {'name': 'Cosine',
                'kwargs': {'Rc': self.Rc}}

    def __repr__(self):
        return ('<Cosine cutoff with Rc=%.3f from amp.descriptor.cutoffs>'
                % self.Rc)
