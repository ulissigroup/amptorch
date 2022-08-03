import os
import pickle
import ase.io
import copy
import numpy as np
import torch
import lmdb
from tqdm import tqdm
from ase import Atoms
from ase.calculators.emt import EMT

from amptorch.trainer import AtomsTrainer
from amptorch.preprocessing import AtomsToData, FeatureScaler, TargetScaler
from amptorch.descriptor.Gaussian import Gaussian


### Construct test data
distances = np.linspace(2, 5, 100)
images = []
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
    images.append(image)

### Construct parameters
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

config = {
    "model": {"get_forces": True, "num_layers": 3, "num_nodes": 5},
    "optim": {
        "force_coefficient": 0.04,
        "lr": 1e-3,
        "batch_size": 10,
        "epochs": 5,
        "loss": "mse",
        "metric": "mae",
    },
    "dataset": {
        "lmdb_path": ["./data1.lmdb", "./data2.lmdb"],
        "val_split": 0,
        "cache": "full",
    },
    "cmd": {
        "debug": False,
        "seed": 1,
        "identifier": "test",
        "verbose": False,
        "logger": False,
    },
}

true_energies = np.array([image.get_potential_energy() for image in images])
true_forces = np.concatenate(np.array([image.get_forces() for image in images]))


def construct_lmdb(images, lmdb_path, normaliers_path="./normalizers.pt"):
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

    # Define symmetry functions
    Gs_1 = copy.deepcopy(Gs)

    training_atoms = images
    elements = np.array([atom.symbol for atoms in training_atoms for atom in atoms])
    elements = np.unique(elements)
    descriptor = Gaussian(Gs=Gs_1, elements=elements, cutoff_func="Cosine")
    descriptor_setup = ("gaussian", Gs_1, {"cutoff_func": "Cosine"}, elements)
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
        scaling = {"type": "normalize", "range": (0, 1)}
        feature_scaler = FeatureScaler(data_list, forcetraining, scaling)
        target_scaler = TargetScaler(data_list, forcetraining)
        normalizers = {
            "target": target_scaler,
            "feature": feature_scaler,
        }
        torch.save(normalizers, normaliers_path)

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


def get_metrics(trainer):
    predictions = trainer.predict(images)
    pred_energies = np.array(predictions["energy"])
    pred_forces = np.concatenate(np.array(predictions["forces"]))
    e_mae = np.mean(np.abs(true_energies - pred_energies))
    f_mae = np.mean(np.abs(pred_forces - true_forces))

    return e_mae, f_mae


def test_lmdb_pretrained():
    construct_lmdb(images, "./data1.lmdb")
    construct_lmdb(images, "./data2.lmdb")

    torch.set_num_threads(1)
    config_1 = copy.deepcopy(config)
    config_2 = copy.deepcopy(config)

    trainer = AtomsTrainer(config_1)
    trainer.train()
    trained_cpdir = trainer.cp_dir
    e_mae_1, f_mae_1 = get_metrics(trainer)

    config_2["optim"]["epochs"] = 100
    pretrained_trainer = AtomsTrainer(config_2)
    pretrained_trainer.load_pretrained(trained_cpdir)
    e_mae_2, f_mae_2 = get_metrics(pretrained_trainer)

    assert e_mae_1 == e_mae_2, "config - lmdb pretrained energy metrics inconsistent!"
    assert f_mae_1 == f_mae_2, "config - lmdb pretrained force metrics inconsistent!"

    pretrained_trainer.train()
    e_mae_3, f_mae_3 = get_metrics(pretrained_trainer)
    assert e_mae_3 < e_mae_2, "Retrained metrics are larger!"
    assert f_mae_3 < f_mae_2, "Retrained metrics are larger!"


def test_lmdb_pretrained_no_config():
    construct_lmdb(images, "./data1.lmdb")
    construct_lmdb(images, "./data2.lmdb")

    config_1 = copy.deepcopy(config)
    trainer = AtomsTrainer(config_1)
    trainer.train()
    trained_cpdir = trainer.cp_dir
    e_mae_1, f_mae_1 = get_metrics(trainer)

    trainer_2 = AtomsTrainer()
    trainer_2.load_pretrained(trained_cpdir)
    e_mae_2, f_mae_2 = get_metrics(trainer_2)

    assert (
        e_mae_1 == e_mae_2
    ), "configless - lmdb pretrained energy metrics inconsistent!"
    assert (
        f_mae_1 == f_mae_2
    ), "configless - lmdb pretrained force metrics inconsistent!"


if __name__ == "__main__":
    print("\n\n--------- LMDB Pretrained Test ---------\n")
    test_lmdb_pretrained()
    test_lmdb_pretrained_no_config()
