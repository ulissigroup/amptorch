import time
import datetime
import os
import random
import warnings
import json

import ase.io
import numpy as np
import skorch.net
import torch
from skorch import NeuralNetRegressor
from skorch.callbacks import LRScheduler
from skorch.dataset import CVSplit
from collections import OrderedDict

from amptorch.dataset import AtomsDataset, DataCollater, construct_descriptor
from amptorch.dataset_lmdb import (
    PartialCacheSampler,
    get_lmdb_dataset,
)
from amptorch.descriptor.util import list_symbols_to_indices
from amptorch.metrics import evaluator
from amptorch.model import BPNN, SingleNN, CustomLoss
from amptorch.preprocessing import AtomsToData
from amptorch.utils import (
    to_tensor,
    train_end_load_best_loss,
    save_normalizers,
    check_memory,
    InOrderSplit,
)
from amptorch.data_parallel import DataParallel, ParallelCollater
from amptorch.ase_utils import AmpTorch

try:
    shell = get_ipython().__class__.__name__
    if shell == "ZMQInteractiveShell":
        from tqdm.notebook import tqdm
    else:
        from tqdm import tqdm
except NameError:
    from tqdm import tqdm


class AtomsTrainer:
    """
    Main trainer class to define the atomistic neural network force field for energy (and force prediction).

    Attribute:
    -----------------
    config: [dict]
        A dictionary that defines configuration of the trainer in `model`, `optim`, `dataset` and `cmd` parts. Please refer to `Usage` and `Example` sections in documentation for more information.
    """

    def __init__(self, config=None):
        self.config = config
        self.pretrained = False

    def load(self, load_dataset=True):
        """
        Loading the parameters passed from `config`.
        """
        self.load_config()
        self.load_rng_seed()
        if load_dataset:
            self.load_dataset()
        self.load_model()
        self.load_criterion()
        self.load_optimizer()
        self.load_logger()
        self.load_extras()
        self.load_skorch()

    def load_config(self):
        """
        Set up attributes from input configuration dictionary.
        """
        dtype = self.config["cmd"].get("dtype", torch.DoubleTensor)
        torch.set_default_tensor_type(dtype)
        self.timestamp = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
        self.identifier = self.config["cmd"].get("identifier", False)
        if self.identifier:
            self.identifier = self.timestamp + "-{}".format(self.identifier)
        else:
            self.identifier = self.timestamp

        self.gpus = self.config["optim"].get("gpus", 0)
        if self.gpus > 0:
            self.output_device = 0
            self.device = f"cuda:{self.output_device}"
        else:
            self.device = "cpu"
            self.output_device = -1
        self.debug = self.config["cmd"].get("debug", False)
        run_dir = self.config["cmd"].get("run_dir", "./")
        os.chdir(run_dir)
        if not self.debug:
            self.cp_dir = os.path.join(run_dir, "checkpoints", self.identifier)
            print(f"Results saved to {self.cp_dir}")
            os.makedirs(self.cp_dir, exist_ok=True)

    def load_rng_seed(self):
        """
        Load a fixed random seed for reproducibility.
        """
        seed = self.config["cmd"].get("seed", 0)
        random.seed(seed)
        np.random.seed(seed)
        torch.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False

    def get_unique_elements(self, training_images):
        """
        Get a list of chemical elements in str if not given in config.
        """
        elements = np.array(
            [atom.symbol for atoms in training_images for atom in atoms]
        )
        elements = np.unique(elements)
        return elements

    def load_dataset(self):
        """
        Load dataset either in the format of a list of `ase.Atoms` objects or a path to those objects, or alternatively, the saved lmdb files with generated fingerprints.
        """
        if "lmdb_path" in self.config["dataset"]:
            self.cache = self.config["dataset"].get("cache", "no")
            self.train_dataset = get_lmdb_dataset(
                self.config["dataset"]["lmdb_path"], self.cache
            )
            self.elements = self.train_dataset.elements
            descriptor_setup = self.train_dataset.descriptor_setup
        else:
            training_images = self.config["dataset"]["raw_data"]
            # TODO: Scalability when dataset to large to fit into memory
            if isinstance(training_images, str):
                training_images = ase.io.read(training_images, ":")
            del self.config["dataset"]["raw_data"]

            self.elements = self.config["dataset"].get(
                "elements", self.get_unique_elements(training_images)
            )

            self.forcetraining = self.config["model"].get("get_forces", True)
            self.fp_scheme = (
                self.config["dataset"].get("fp_scheme", "gmpordernorm").lower()
            )
            self.fp_params = self.config["dataset"]["fp_params"]
            self.save_fps = self.config["dataset"].get("save_fps", True)
            self.cutoff_params = self.config["dataset"].get(
                "cutoff_params", {"cutoff_func": "Cosine"}
            )
            descriptor_setup = (
                self.fp_scheme,
                self.fp_params,
                self.cutoff_params,
                self.elements,
            )
            self.train_dataset = AtomsDataset(
                images=training_images,
                descriptor_setup=descriptor_setup,
                forcetraining=self.forcetraining,
                save_fps=self.config["dataset"].get("save_fps", True),
                scaling=self.config["dataset"].get(
                    "scaling",
                    {"type": "normalize", "range": (0, 1), "elementwise": True},
                ),
            )
        self.feature_scaler = self.train_dataset.feature_scaler
        self.target_scaler = self.train_dataset.target_scaler
        self.atomic_correction_scaler = None
        # self.atomic_correction_scaler = self.train_dataset.atomic_correction_scaler
        self.input_dim = self.train_dataset.input_dim
        self.val_split = self.config["dataset"].get("val_split", 0)
        if not self.debug:
            normalizers = {
                "target": self.target_scaler,
                "feature": self.feature_scaler,
            }
            torch.save(normalizers, os.path.join(self.cp_dir, "normalizers.pt"))
            self.config["dataset"]["descriptor"] = descriptor_setup
            self.config["dataset"]["fp_length"] = self.input_dim
            torch.save(self.config, os.path.join(self.cp_dir, "config.pt"))
        print("Loading dataset: {} images".format(len(self.train_dataset)))

    def load_model(self):
        """
        Load the parameters for atomistic neural network models from config.
        """
        elements = list_symbols_to_indices(self.elements)
        model = self.config["model"].get("name", "bpnn").lower()
        if model == "bpnn":
            self.model = BPNN(
                elements=elements, input_dim=self.input_dim, **self.config["model"]
            )
        elif model == "singlenn":
            self.model = SingleNN(
                elements=elements, input_dim=self.input_dim, **self.config["model"]
            )
        else:
            raise NotImplementedError(f"{model} not supported!")
        print("Loading model: {} parameters".format(self.model.num_params))
        self.forcetraining = self.config["model"].get("get_forces", True)
        collate_fn = DataCollater(train=True, forcetraining=self.forcetraining)
        self.parallel_collater = ParallelCollater(self.gpus, collate_fn)
        if self.gpus > 0:
            self.model = DataParallel(
                self.model,
                output_device=self.output_device,
                num_gpus=self.gpus,
            )

    def load_extras(self):
        """
        Load extra commands such as validation splits, CV splits, debug mode and various callbacks.
        """
        callbacks = []
        load_best_loss = train_end_load_best_loss(self.identifier)
        self.val_split = self.config["dataset"].get("val_split", 0)
        self.split_mode = self.config["dataset"].get("val_split_mode", "cv")
        if self.split_mode == "cv":
            self.split = CVSplit(cv=self.val_split) if self.val_split != 0 else 0
        elif self.split_mode == "inorder":
            self.split = InOrderSplit(self.val_split) if self.val_split != 0 else 0
        else:
            raise NotImplementedError

        metrics = evaluator(
            self.val_split,
            self.config["optim"].get("metric", "mae"),
            self.identifier,
            self.forcetraining,
            self.config["optim"].get("cp_metric", "energy"),
        )
        callbacks.extend(metrics)

        if not self.debug:
            callbacks.append(load_best_loss)
        scheduler = self.config["optim"].get("scheduler", None)
        if scheduler:
            scheduler = LRScheduler(scheduler["policy"], **scheduler["params"])
            callbacks.append(scheduler)
        if self.config["cmd"].get("logger", False):
            from skorch.callbacks import WandbLogger

            callbacks.append(
                WandbLogger(
                    self.wandb_run,
                    save_model=False,
                    keys_ignored="dur",
                )
            )

        # early stopping
        if self.config["cmd"].get("early_stopping", False):
            from skorch.callbacks import EarlyStopping

            callbacks.append(
                EarlyStopping(
                    patience=self.config["cmd"].get("early_stoppping_patience", 5)
                )
            )

        self.callbacks = callbacks

    def load_criterion(self):
        """
        Load custom loss function in optimizaiton scheme from config.
        """
        self.criterion = self.config["optim"].get("loss_fn", CustomLoss)

    def load_optimizer(self):
        """
        Load set-up parameters in optimizers from config.
        """
        self.optimizer = {
            "optimizer": self.config["optim"].get("optimizer", torch.optim.Adam)
        }
        optimizer_args = self.config["optim"].get("optimizer_args", False)
        if optimizer_args:
            self.optimizer.update(optimizer_args)

    def load_logger(self):
        """
        Load wandb logger from config.
        """
        if self.config["cmd"].get("logger", False):
            import wandb

            self.wandb_run = wandb.init(
                name=self.identifier,
                config=self.config,
            )

    def load_skorch(self):
        """
        Load the skorch atomistic neural network regression model with parameters from config.
        """
        skorch.net.to_tensor = to_tensor

        if self.config["dataset"].get("cache", None) == "partial":
            self.net = NeuralNetRegressor(
                module=self.model,
                criterion=self.criterion,
                criterion__force_coefficient=self.config["optim"].get(
                    "force_coefficient", 0
                ),
                criterion__loss=self.config["optim"].get("loss", "mse"),
                lr=self.config["optim"].get("lr", 1e-1),
                batch_size=self.config["optim"].get("batch_size", 32),
                max_epochs=self.config["optim"].get("epochs", 100),
                iterator_train__collate_fn=self.parallel_collater,
                iterator_train__sampler=PartialCacheSampler(
                    self.train_dataset.get_length_list(),
                    self.val_split,
                ),
                iterator_train__shuffle=False,
                iterator_train__pin_memory=True,
                iterator_valid__collate_fn=self.parallel_collater,
                iterator_valid__shuffle=False,
                iterator_valid__pin_memory=True,
                device=self.device,
                train_split=self.split,
                callbacks=self.callbacks,
                verbose=self.config["cmd"].get("verbose", True),
                **self.optimizer,
            )
        # turn off pin memory for gaussian symmetry function force training, as torch.geometric outputs error
        elif self.config["dataset"].get(
            "fp_scheme", "gaussian"
        ).lower() == "gaussian" and self.config["model"].get("get_forces", True):
            self.net = NeuralNetRegressor(
                module=self.model,
                criterion=self.criterion,
                criterion__force_coefficient=self.config["optim"].get(
                    "force_coefficient", 0
                ),
                criterion__loss=self.config["optim"].get("loss", "mse"),
                lr=self.config["optim"].get("lr", 1e-1),
                batch_size=self.config["optim"].get("batch_size", 32),
                max_epochs=self.config["optim"].get("epochs", 100),
                iterator_train__collate_fn=self.parallel_collater,
                iterator_train__shuffle=True,
                iterator_train__pin_memory=False,
                iterator_valid__collate_fn=self.parallel_collater,
                iterator_valid__shuffle=False,
                iterator_valid__pin_memory=False,
                device=self.device,
                train_split=self.split,
                callbacks=self.callbacks,
                verbose=self.config["cmd"].get("verbose", True),
                **self.optimizer,
            )
        else:
            self.net = NeuralNetRegressor(
                module=self.model,
                criterion=self.criterion,
                criterion__force_coefficient=self.config["optim"].get(
                    "force_coefficient", 0
                ),
                criterion__loss=self.config["optim"].get("loss", "mse"),
                lr=self.config["optim"].get("lr", 1e-1),
                batch_size=self.config["optim"].get("batch_size", 32),
                max_epochs=self.config["optim"].get("epochs", 100),
                iterator_train__collate_fn=self.parallel_collater,
                iterator_train__shuffle=True,
                iterator_train__pin_memory=True,
                iterator_valid__collate_fn=self.parallel_collater,
                iterator_valid__shuffle=False,
                iterator_valid__pin_memory=True,
                device=self.device,
                train_split=self.split,
                callbacks=self.callbacks,
                verbose=self.config["cmd"].get("verbose", True),
                **self.optimizer,
            )
        print("Loading skorch trainer")

    def train(self, raw_data=None):
        """
        Method used to train the model with defined config by initiating the AtomsTrainer instance with a user-defined config dictionary. Can be fed with a list of `ase.Atoms` objects as training data.
        """
        if raw_data is not None:
            self.config["dataset"]["raw_data"] = raw_data
        if not self.pretrained:
            self.load()

        stime = time.time()
        self.net.fit(self.train_dataset, None)
        elapsed_time = time.time() - stime
        print(f"Training completed in {elapsed_time}s")

    def predict(
        self,
        images,
        disable_tqdm=True,
        get_latent=None,
        get_descriptor=False,
        save_fps=False,
    ):
        """
        Method used to make energy (and force) predictions for input images.

        Attributes:
        --------------------
        images : List[ase.Atoms]
            A list of ase.Atoms objects for the model to make predictions on.

        disabled_tqdm : bool
            Option to disable tqdm for Jupyter notebook compatibility.

        get_latent : int (default to None)
            Option to record the last-layer latent representation for uncertainty quantification purpose. Usually used with -2 to get the last-layer neural network latent representation. If this flag is specified to an integer instead of default None, the output dictionary would contain an extra entry of "latent".

        get_descriptor : bool
            Option to record the feature space representation. If set to True, , the output dictionary would contain an extra entry of "feature".

        save_fps : bool
            Option to save the calculated fingerprints for accelerated computation.

        Output:
        -------------------
        predictions : dict
            A dictionary that contains two basic entries, "energy" and "forces". prediction["energy"] is a list of energies, and prediction["forces"] is a list of force components in Cartesian coordiantes, if only the model is trained with energy and forces.
        """
        if len(images) < 1:
            warnings.warn("No images found!", stacklevel=2)
            return images

        self.descriptor = construct_descriptor(self.config["dataset"]["descriptor"])

        t0 = time.time()
        a2d = AtomsToData(
            descriptor=self.descriptor,
            r_energy=False,
            r_forces=False,
            save_fps=save_fps,
            fprimes=self.forcetraining,
            cores=1,
        )

        data_list = a2d.convert_all(images, disable_tqdm=disable_tqdm)

        self.feature_scaler.norm(data_list, disable_tqdm=disable_tqdm)
        t_fingerPrint = time.time() - t0

        self.net.module.eval()
        collate_fn = DataCollater(train=False, forcetraining=self.forcetraining)

        predictions = {"energy": [], "forces": []}

        t_forwardPass = 0

        # get the latent layer
        if get_latent is not None:
            predictions["latent"] = []

            def hook2get_latent(self, input, output):
                _latent = output.detach().numpy()
                _latent_mean = np.mean(_latent, axis=0)
                predictions["latent"].append(_latent_mean)

            latent_layer = get_latent

            self.net.module.model.model_net[latent_layer].register_forward_hook(
                hook2get_latent
            )

        # get feature descriptor for every image in the trajectory by averaging over atoms.
        if get_descriptor:
            predictions["descriptors"] = []
            for data in data_list:
                _feature = data.fingerprint.cpu().detach().numpy()
                _feature_mean = np.mean(_feature, axis=0)
                predictions["descriptors"].append(_feature_mean)

        # for data in data_list:
        for idx, data in tqdm(
            enumerate(data_list),
            desc="Predicting",
            total=len(data_list),
            unit=" systems",
            disable=disable_tqdm,
        ):
            collated = collate_fn([data]).to(self.device)

            t0 = time.time()
            energy, forces = self.net.module([collated])
            energy = self.target_scaler.denorm(
                energy.detach().cpu(), pred="energy"
            ).tolist()
            # if self.atomic_correction_scaler is not None:
            #     energy = self.atomic_correction_scaler.denorm(energy, data_list[idx])
            forces = self.target_scaler.denorm(
                forces.detach().cpu(), pred="forces"
            ).numpy()
            t_forwardPass += time.time() - t0

            predictions["energy"].extend(energy)
            predictions["forces"].append(forces)

        # time fingerprinting and neural network passing
        predictions["t_fingerPrint"] = t_fingerPrint
        predictions["t_forwardPass"] = t_forwardPass

        return predictions

    def load_pretrained(self, checkpoint_path=None, gpu2cpu=False):
        """
        Load pretrained model with configuration and parameters in the checkpoint.

        Args:
            checkpoint_path: str, Path to checkpoint directory
            gpu2cpu: bool, True if checkpoint was trained with GPUs and you
            wish to load on cpu instead.
        """

        self.pretrained = True
        print(f"Loading checkpoint from {checkpoint_path}")
        assert os.path.isdir(
            checkpoint_path
        ), f"Checkpoint: {checkpoint_path} not found!"
        if not self.config:
            # prediction only
            self.config = torch.load(os.path.join(checkpoint_path, "config.pt"))
            self.config["cmd"]["debug"] = True
            self.elements = self.config["dataset"]["descriptor"][-1]
            self.input_dim = self.config["dataset"]["fp_length"]
            if gpu2cpu:
                self.config["optim"]["gpus"] = 0
            self.load(load_dataset=False)
        else:
            # prediction+retraining
            self.load(load_dataset=True)
        self.net.initialize()

        if gpu2cpu:
            params_path = os.path.join(checkpoint_path, "params_cpu.pt")
            if not os.path.exists(params_path):
                params = torch.load(
                    os.path.join(checkpoint_path, "params.pt"),
                    map_location=torch.device("cpu"),
                )
                new_dict = OrderedDict()
                for k, v in params.items():
                    name = k[7:]
                    new_dict[name] = v
                torch.save(new_dict, params_path)
        else:
            params_path = os.path.join(checkpoint_path, "params.pt")

        try:
            self.net.load_params(
                f_params=params_path,
                f_optimizer=os.path.join(checkpoint_path, "optimizer.pt"),
                f_criterion=os.path.join(checkpoint_path, "criterion.pt"),
                f_history=os.path.join(checkpoint_path, "history.json"),
            )
            normalizers = torch.load(os.path.join(checkpoint_path, "normalizers.pt"))
            self.feature_scaler = normalizers["feature"]
            self.target_scaler = normalizers["target"]
        except NotImplementedError:
            print("Unable to load checkpoint!")

    def get_calc(self):
        """
        Convert the AtomsTrainer class to an `ase.Calculator` class for interfacing with ase.

        Output:
        -----------
        AmpTorch : ase.Calculator class
            After attaching the Calculator to `ase.Atoms` object, the user can use `get_potential_energy()` method to obtain the corresponding energy in ase.
        """
        return AmpTorch(self)
