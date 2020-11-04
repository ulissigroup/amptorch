import pickle
import os
import glob
import lmdb
import numpy as np
import ase.io
from tqdm import tqdm
from amptorch.preprocessing import AtomsToData, FeatureScaler, TargetScaler
from amptorch.descriptor.Gaussian import Gaussian


def construct_lmdb(data_dir, lmdb_path="./data.lmdb"):
    """
    data_dir: Directory containing traj files to construct dataset from
    lmdb_path: Path to store LMDB dataset
    """
    db = lmdb.open(
        lmdb_path,
        map_size=1099511627776 * 2,
        subdir=False,
        meminit=False,
        map_async=True,
    )

    paths = glob.glob(os.path.join(data_dir, "*.traj"))
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

    elements = ["Cu", "C", "O"]
    descriptor = Gaussian(Gs=Gs, elements=elements, cutoff_func="Cosine")
    descriptor_setup = ("gaussian", Gs, elements, {"cutoff_func": "Cosine"})
    scaling = {"type": "normalize", "range": (0, 1)}
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
    for path in tqdm(paths[:2]):
        images = ase.io.read(path, ":")
        for image in images:
            do = a2d.convert(image, idx=idx)
            txn = db.begin(write=True)
            txn.put(f"{idx}".encode("ascii"), pickle.dumps(do, protocol=-1))
            txn.commit()
            data_list.append(do)  # suppress if hitting memory limits

            # get summary statistics for at most 20k random data objects
            # control percentage depending on your dataset size
            # default: sample point with 50% probability

            # unsuppress the following if using a very large dataset
            # if random.randint(0, 100) < 50 and len(data_list) < 20000:
            # data_list.append(do)
            idx += 1

    feature_scaler = FeatureScaler(data_list, forcetraining, scaling)
    txn = db.begin(write=True)
    txn.put("feature_scaler".encode("ascii"), pickle.dumps(feature_scaler, protocol=-1))
    txn.commit()

    target_scaler = TargetScaler(data_list, forcetraining)
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
    data_dir = "/home/jovyan/projects/amptorch/examples/data"  # directory of traj files
    construct_lmdb(data_dir, lmdb_path="data.lmdb")
