import os
import pickle
import lmdb
import numpy as np
import ase.io
import torch
from tqdm import tqdm
from amptorch.preprocessing import AtomsToData, FeatureScaler, TargetScaler
from amptorch.descriptor.GMPOrderNorm import GMPOrderNorm
from ase import Atoms
from ase.calculators.emt import EMT


def construct_lmdb(images, lmdb_path="./data.lmdb", normaliers_path="./normalizers.pt"):
    """
    images: list of ase atoms objects (or trajectory) for fingerprint calculatation
    lmdb_path: Path to store LMDB dataset.
    normaliers_path: path of the scalers, create and store them if not exist
    """

    db = lmdb.open(
        lmdb_path,
        map_size=1099511627776 * 2,
        subdir=False,
        meminit=False,
        map_async=True,
    )

    # Define GMPs
    nsigmas = 5  # number of radial probes
    max_MCSH_order = 3  # order of angular probes
    max_radial_sigma = 2.0  # the maximal sigma of gaussian in radial coordiantes

    ### Construct GMP configuration, no need to change once the hyperparameters are specified.
    sigmas = np.linspace(0, max_radial_sigma, nsigmas + 1, endpoint=True)[1:]
    GMPs = {
        "MCSHs": {"orders": list(range(max_MCSH_order + 1)), "sigmas": sigmas},
    }

    training_atoms = images
    elements = np.array([atom.symbol for atoms in training_atoms for atom in atoms])
    elements = np.unique(elements)
    descriptor = GMPOrderNorm(MCSHs=GMPs, elements=elements)
    descriptor_setup = ("gmpordernorm", GMPs, "NA", elements)
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
    for image in tqdm(
        images,
        desc="calculating fps",
        total=len(images),
        unit=" images",
    ):
        do = a2d.convert(image, idx=idx)
        data_list.append(do)
        idx += 1

    if os.path.isfile(normaliers_path):
        normalizers = torch.load(normaliers_path)
        feature_scaler = normalizers["feature"]
        target_scaler = normalizers["target"]

    else:
        scaling = {"type": "normalize", "range": (-1, 1)}
        feature_scaler = FeatureScaler(data_list, forcetraining, scaling)
        target_scaler = TargetScaler(data_list, forcetraining)
        normalizers = {
            "target": target_scaler,
            "feature": feature_scaler,
        }
        torch.save(normalizers, normaliers_path)

    feature_scaler.norm(data_list)
    target_scaler.norm(data_list)

    print(data_list[0].fingerprint.dtype)

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
    torch.set_default_tensor_type(torch.DoubleTensor)

    images = []

    distances = np.linspace(2, 5, 1000)
    for dist in distances:
        image = Atoms(
            "CuCO",
            [
                (-dist * np.sin(0.65), dist * np.cos(0.65), 0),
                (0, 0, 0),
                (dist * np.sin(0.65), dist * np.cos(0.65), 0),
            ],
        )
        image.set_cell([10, 10, 10])
        image.wrap(pbc=True)
        image.set_calculator(EMT())
        image.get_potential_energy()
        images.append(image)

    construct_lmdb(images, lmdb_path="./data.lmdb", normaliers_path="./normalizers.pt")
