"""Modified utilities and descriptors extracted from AMP's source code"""

import os
import pickle
import tarfile
import time
from .utils import Cosine, dict2cutoff
from ase.calculators.calculator import Parameters
from copy import deepcopy
from .utils import Logger


class SNN_Gaussian(object):
    """Class that calculates Gaussian fingerprints (i.e., Behler-style).

    Parameters
    ----------
    cutoff : object or float
        Cutoff function, typically from amp.descriptor.cutoffs.  Can be also
        fed as a float representing the radius above which neighbor
        interactions are ignored; in this case a cosine cutoff function will be
        employed.  Default is a 6.5-Angstrom cosine cutoff.
    Gs : dict or list
        Dictionary of symbols and lists of dictionaries for making symmetry
        functions. Either auto-genetrated, or given in the following form, for
        example:

               >>> Gs = {"O": [{"type":"G2", "element":"O", "eta":10.},
               ...             {"type":"G4", "elements":["O", "Au"],
               ...              "eta":5., "gamma":1., "zeta":1.0}],
               ...       "Au": [{"type":"G2", "element":"O", "eta":2.},
               ...              {"type":"G4", "elements":["O", "Au"],
               ...               "eta":2., "gamma":1., "zeta":5.0}]}

        You can use amp.model.gaussian.make_symmetry_functions to help create
        these lists of dictionaries.  If you supply a list instead of a
        dictionary, it will assume you want identical symmetry functions for
        each element.
    dblabel : str
        Optional separate prefix/location for database files, including
        fingerprints, fingerprint derivatives, and neighborlists. This file
        location can be shared between calculator instances to avoid
        re-calculating redundant information. If not supplied, just uses the
        value from label.
    elements : list
        List of allowed elements present in the system. If not provided, will
        be found automatically.
    version : str
        Version of fingerprints.
    fortran : bool
        If True, will use fortran modules, if False, will not.
    mode : str
        Can be either 'atom-centered' or 'image-centered'.

    Raises
    ------
        RuntimeError
    """

    def __init__(
        self,
        cutoff=Cosine(6.5),
        Gs=None,
        dblabel=None,
        elements=None,
        version=None,
        fortran=True,
        mode="atom-centered",
    ):

        # Check of the version of descriptor, particularly if restarting.
        compatibleversions = ["2015.12"]
        if (version is not None) and version not in compatibleversions:
            raise RuntimeError(
                "Error: Trying to use Gaussian fingerprints"
                " version %s, but this module only supports"
                " versions %s. You may need an older or "
                " newer version of Amp." % (version, compatibleversions)
            )
        else:
            version = compatibleversions[-1]

        # Check that the mode is atom-centered.
        if mode != "atom-centered":
            raise RuntimeError(
                "Gaussian scheme only works "
                "in atom-centered mode. %s "
                "specified." % mode
            )

        # If the cutoff is provided as a number, Cosine function will be used
        # by default.
        if isinstance(cutoff, int) or isinstance(cutoff, float):
            cutoff = Cosine(cutoff)
        # If the cutoff is provided as a dictionary, assume we need to load it
        # with dict2cutoff.
        if type(cutoff) is dict:
            cutoff = dict2cutoff(cutoff)

        # The parameters dictionary contains the minimum information
        # to produce a compatible descriptor; that is, one that gives
        # an identical fingerprint when fed an ASE image.
        p = self.parameters = Parameters(
            {"importname": ".descriptor.gaussian.Gaussian", "mode": "atom-centered"}
        )
        p.version = version
        p.cutoff = cutoff.todict()
        p.Gs = Gs
        p.elements = elements

        self.dblabel = dblabel
        self.fortran = fortran
        self.parent = None  # Can hold a reference to main Amp instance.

    def tostring(self):
        """Returns an evaluatable representation of the calculator that can
        be used to restart the calculator.
        """
        return self.parameters.tostring()

    def calculate_fingerprints(
        self, images, parallel=None, log=None, calculate_derivatives=False):
        """Calculates the fingerpints of the images, for the ones not already
        done.

        Parameters
        ----------
        images : dict
            Dictionary of images; the key is a unique ID assigned to each
            image and each value is an ASE atoms object. Typically created
            from amp.utilities.hash_images.
        parallel : dict
            Configuration for parallelization. Should be in same form as in
            amp.Amp.
        log : Logger object
            Write function at which to log data. Note this must be a callable
            function.
        calculate_derivatives : bool
            Decides whether or not fingerprintprimes should also be calculated.
        """
        if (self.dblabel is None) and hasattr(self.parent, "dblabel"):
            self.dblabel = self.parent.dblabel
        self.dblabel = "amp-data" if self.dblabel is None else self.dblabel

        p = self.parameters

        if p.elements is None:
            p.elements = set(
                [atom.symbol for atoms in images.values() for atom in atoms]
            )
        p.elements = sorted(p.elements)

        if not hasattr(p.Gs, "keys"):
            p.Gs = {element: deepcopy(p.Gs) for element in p.elements}
        # TODO Ensure this is not needed
        if not hasattr(self, "neighborlist"):
            calc = NeighborlistCalculator(cutoff=p.cutoff["kwargs"]["Rc"])
            self.neighborlist = Data(
                filename="%s-neighborlists" % self.dblabel, calculator=calc
            )
        self.neighborlist.calculate_items(images)

        if not hasattr(self, "fingerprints"):
            self.fingerprints = Data(
                filename="%s-fingerprints" % self.dblabel, calculator=None
            )
        if calculate_derivatives:
            if not hasattr(self, "fingerprintprimes"):
                self.fingerprintprimes = Data(
                    filename="%s-fingerprint-primes" % self.dblabel, calculator=None
                )


# Neighborlist Calculator
class NeighborlistCalculator:
    """For integration with .utilities.Data

    For each image fed to calculate, a list of neighbors with offset distances
    is returned.

    Parameters
    ----------
    cutoff : float
        Radius above which neighbor interactions are ignored.
    """

    def __init__(self, cutoff):
        self.globals = Parameters({"cutoff": cutoff})

    def calculate(self, image, key):
        """For integration with .utilities.Data

        For each image fed to calculate, a list of neighbors with offset
        distances is returned.

        Parameters
        ----------
        image : object
            ASE atoms object.
        key : str
            key of the image after being hashed.
        """
        from ase.neighborlist import NeighborList, NewPrimitiveNeighborList

        cutoff = self.globals.cutoff
        n = NeighborList(cutoffs=[cutoff / 2.0] * len(image),
                self_interaction=False, primitive=NewPrimitiveNeighborList)
        n.update(image)
        return [n.get_neighbors(index) for index in range(len(image))]


class FileDatabase:
    """Using a database file, such as shelve or sqlitedict, that can handle
    multiple processes writing to the file is hard.

    Therefore, we take the stupid approach of having each database entry be
    a separate file. This class behaves essentially like shelve, but saves each
    dictionary entry as a plain pickle file within the directory, with the
    filename corresponding to the dictionary key (which must be a string).

    Like shelve, this also keeps an internal (memory dictionary) representation
    of the variables that have been accessed.

    Also includes an archive feature, where files are instead added to a file
    called 'archive.tar.gz' to save disk space. If an entry exists in both the
    loose and archive formats, the loose is taken to be the new (correct)
    value.
    """

    def __init__(self, filename):
        """Open the filename at specified location. flag is ignored; this
        format is always capable of both reading and writing."""
        if not filename.endswith(os.extsep + "ampdb"):
            filename += os.extsep + "ampdb"
        self.path = filename
        self.loosepath = os.path.join(self.path, "loose")
        self.tarpath = os.path.join(self.path, "archive.tar.gz")
        if not os.path.exists(self.path):
            try:
                os.mkdir(self.path)
            except OSError:
                # Many simultaneous processes might be trying to make the
                # directory at the same time.
                pass
            try:
                os.mkdir(self.loosepath)
            except OSError:
                pass
        self._memdict = {}  # Items already accessed; stored in memory.

    @classmethod
    def open(Cls, filename, flag=None):
        """Open present for compatibility with shelve. flag is ignored; this
        format is always capable of both reading and writing.
        """
        return Cls(filename=filename)

    def close(self):
        """Only present for compatibility with shelve.
        """
        return

    def keys(self):
        """Return list of keys, both of in-memory and out-of-memory
        items.
        """
        keys = os.listdir(self.loosepath)
        if os.path.exists(self.tarpath):
            with tarfile.open(self.tarpath) as tf:
                keys = list(set(keys + tf.getnames()))
        return keys

    def values(self):
        """Return list of values, both of in-memory and out-of-memory
        items. This moves all out-of-memory items into memory.
        """
        keys = self.keys()
        return [self[key] for key in keys]

    def __len__(self):
        return len(self.keys())

    def __setitem__(self, key, value):
        self._memdict[key] = value
        path = os.path.join(self.loosepath, str(key))
        if os.path.exists(path):
            with open(path, "rb") as f:
                contents = self._repeat_read(f)
                if pickle.dumps(contents) == pickle.dumps(value):
                    # Using pickle as a hash...
                    return  # Nothing to update.
        with open(path, "wb") as f:
            pickle.dump(value, f, protocol=0)

    def _repeat_read(self, f, maxtries=5, sleep=0.2):
        """If one process is writing, the other process cannot read without
        errors until it finishes. Reads file-like object f checking for
        errors, and retries up to 'maxtries' times, sleeping 'sleep' sec
        between tries."""
        tries = 0
        while tries < maxtries:
            try:
                contents = pickle.load(f)
            except (UnicodeDecodeError, EOFError, pickle.UnpicklingError):
                time.sleep(0.2)
                tries += 1
            else:
                return contents
        raise IOError("Too many file read attempts.")

    def __getitem__(self, key):
        if key in self._memdict:
            return self._memdict[key]
        keypath = os.path.join(self.loosepath, key)
        if os.path.exists(keypath):
            with open(keypath, "rb") as f:
                return self._repeat_read(f)
        elif os.path.exists(self.tarpath):
            with tarfile.open(self.tarpath) as tf:
                return pickle.load(tf.extractfile(key))
        else:
            raise KeyError(str(key))

    def update(self, newitems):
        for key, value in newitems.items():
            self.__setitem__(key, value)

class Data:
    """Serves as a container (dictionary-like) for (key, value) pairs that
    also serves to calculate them.

    Works by default with python's shelve module, but something that is built
    to share the same commands as shelve will work fine; just specify this in
    dbinstance.

    Designed to hold things like neighborlists, which have a hash, value
    format.

    This will work like a dictionary in that items can be accessed with
    data[key], but other advanced dictionary functions should be accessed with
    through the .d attribute:

    >>> data = Data(...)
    >>> data.open()
    >>> keys = data.d.keys()
    >>> values = data.d.values()
    """

    def __init__(self, filename, db=FileDatabase, calculator=None):
        self.calc = calculator
        self.db = db
        self.filename = filename
        self.d = None

    def calculate_items(self, images):
        d = self.db.open(self.filename, "r")
        calcs_needed = list(set(images.keys()).difference(d.keys()))
        d.close()
        if len(calcs_needed) == 0:
            return
        d = self.db.open(self.filename, "c")
        for key in calcs_needed:
            d[key] = self.calc.calculate(images[key], key)
        d.close()
        self.d = None

    def __getitem__(self, key):
        self.open()
        return self.d[key]

    def close(self):
        """Safely close the database.
        """
        if self.d:
            self.d.close()
        self.d = None

    def open(self, mode="r"):
        """Open the database connection with mode specified.
        """
        if self.d is None:
            self.d = self.db.open(self.filename, mode)

    def __del__(self):
        self.close()


def make_symmetry_functions(elements, type, etas, offsets=None,
                            zetas=None, gammas=None):
    """Helper function to create Gaussian symmetry functions.
    Returns a list of dictionaries with symmetry function parameters
    in the format expected by the Gaussian class.

    Parameters
    ----------
    elements : list of str
        List of element types to be observed in this fingerprint.
    type : str
        Either G2, G4, or G5.
    etas : list of floats
        eta values to use in G2, G4 or G5 fingerprints
    offsets: list of floats
        offset values to use in G2 fingerprints
    zetas : list of floats
        zeta values to use in G4, and G5 fingerprints
    gammas : list of floats
        gamma values to use in G4, and G5 fingerprints

    Returns
    -------
    G : list of dicts
        A list, each item in the list contains a dictionary of fingerprint
        parameters.
    """
    if type == 'G2':
        offsets = [0.] if offsets is None else offsets
        G = [{'type': 'G2', 'element': element, 'eta': eta, 'offset': offset}
             for eta in etas
             for element in elements
             for offset in offsets]
        return G
    elif type == 'G4':
        G = []
        for eta in etas:
            for zeta in zetas:
                for gamma in gammas:
                    for i1, el1 in enumerate(elements):
                        for el2 in elements[i1:]:
                            els = sorted([el1, el2])
                            G.append({'type': 'G4',
                                      'elements': els,
                                      'eta': eta,
                                      'gamma': gamma,
                                      'zeta': zeta})
        return G
    elif type == 'G5':
        G = []
        for eta in etas:
            for zeta in zetas:
                for gamma in gammas:
                    for i1, el1 in enumerate(elements):
                        for el2 in elements[i1:]:
                            els = sorted([el1, el2])
                            G.append({'type': 'G5',
                                      'elements': els,
                                      'eta': eta,
                                      'gamma': gamma,
                                      'zeta': zeta})
        return G
    raise NotImplementedError('Unknown type: {}.'.format(type))
