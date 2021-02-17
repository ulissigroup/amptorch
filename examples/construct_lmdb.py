import pickle
import lmdb
import numpy as np
import ase.io
from tqdm import tqdm
from amptorch.preprocessing import AtomsToData, FeatureScaler, TargetScaler
from amptorch.descriptor.Gaussian import Gaussian


def construct_lmdb(path, lmdb_path="./data.lmdb"):
    """
    path: Path to trajectory/ASE-compatible file
    lmdb_path: Path to store LMDB dataset.
    """
    db = lmdb.open(
        lmdb_path,
        map_size=1099511627776 * 2,
        subdir=False,
        meminit=False,
        map_async=True,
    )

    # Define symmetry functions
    Gs = {
        "default": {
            "G2": {
                "etas": np.logspace(np.log10(0.05), np.log10(5.0), num=4),
                "rs_s": [0],
            },
            "G4": {"etas": [0.005], "zetas": [1.0, 4.0], "gammas": [1.0, -1.0]},
            "cutoff": 6,
        },
    }

    training_atoms = ase.io.read(path, ":")
    elements = np.array([atom.symbol for atoms in training_atoms for atom in atoms])
    elements = np.unique(elements)
    descriptor = Gaussian(Gs=Gs, elements=elements, cutoff_func="Cosine")
    descriptor_setup = ("gaussian", Gs, {"cutoff_func": "Cosine"}, elements)
    forcetraining = True

    a2d = AtomsToData(
        descriptor=descriptor,
        r_energy=True,
        r_forces=True,
        save_fps=False,
        fprimes=forcetraining,
    )

    data_list = []
    idx = 0
    images = ase.io.read(path, ":")
    for image in images:
        do = a2d.convert(image, idx=idx)
        data_list.append(do)
        idx += 1

    scaling = {"type": "normalize", "range": (0, 1)}
    feature_scaler = FeatureScaler(data_list, forcetraining, scaling)
    target_scaler = TargetScaler(data_list, forcetraining)

    feature_scaler.norm(data_list)
    target_scaler.norm(data_list)

    idx = 0
    for do in tqdm(data_list, desc="Writing images to LMDB"):
        txn = db.begin(write=True)
        txn.put(f"{idx}".encode("ascii"), pickle.dumps(do, protocol=-1))
        txn.commit()
        idx += 1

    txn = db.begin(write=True)
    txn.put("feature_scaler".encode("ascii"), pickle.dumps(feature_scaler, protocol=-1))
    txn.commit()

    txn = db.begin(write=True)
    txn.put("target_scaler".encode("ascii"), pickle.dumps(target_scaler, protocol=-1))
    txn.commit()

    txn = db.begin(write=True)
    txn.put("length".encode("ascii"), pickle.dumps(idx, protocol=-1))
    txn.commit()

    txn = db.begin(write=True)
    txn.put("elements".encode("ascii"), pickle.dumps(elements, protocol=-1))
    txn.commit()

    txn = db.begin(write=True)
    txn.put(
        "descriptor_setup".encode("ascii"), pickle.dumps(descriptor_setup, protocol=-1)
    )
    txn.commit()

    db.sync()
    db.close()


if __name__ == "__main__":
    traj_path = "path/to/traj/file.traj"
    construct_lmdb(traj_path, lmdb_path="data.lmdb")
