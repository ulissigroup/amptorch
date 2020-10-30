import datetime
import os
import random
import warnings

import ase.io
import numpy as np
import skorch.net
import torch
from skorch import NeuralNetRegressor
from skorch.callbacks import LRScheduler
from skorch.dataset import CVSplit

from amptorch.dataset import AtomsDataset, DataCollater
from amptorch.descriptor.util import list_symbols_to_indices
from amptorch.metrics import evaluator
from amptorch.model import BPNN, CustomLoss
from amptorch.preprocessing import AtomsToData
from amptorch.utils import to_tensor, train_end_load_best_loss


class AtomsTrainer:
    def __init__(self, config):
        self.config = config
        self.pretrained = False

    def load(self):
        self.load_config()
        self.load_rng_seed()
        self.load_dataset()
        self.load_model()
        self.load_criterion()
        self.load_optimizer()
        self.load_logger()
        self.load_extras()
        self.load_skorch()

    def load_config(self):
        self.timestamp = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
        self.identifier = self.config["cmd"].get("identifier", False)
        if self.identifier:
            self.identifier = self.timestamp + "-{}".format(self.identifier)
        else:
            self.identifier = self.timestamp

        self.device = torch.device(self.config["optim"].get("device", "cpu"))
        self.debug = self.config["cmd"].get("debug", False)
        run_dir = self.config["cmd"].get("run_dir", "./")
        os.chdir(run_dir)
        if not self.debug:
            self.cp_dir = os.path.join(run_dir, "checkpoints", self.identifier)
            print(f"Results saved to {self.cp_dir}")
            os.makedirs(self.cp_dir, exist_ok=True)

    def load_rng_seed(self):
        seed = self.config["cmd"].get("seed", 0)
        random.seed(seed)
        np.random.seed(seed)
        torch.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False

    def get_unique_elements(self, training_images):
        elements = np.array(
            [atom.symbol for atoms in training_images for atom in atoms]
        )
        elements = np.unique(elements)
        return elements

    def load_dataset(self):
        training_images = self.config["dataset"]["raw_data"]
        # TODO: Scalability when dataset to large to fit into memory
        if isinstance(training_images, str):
            training_images = ase.io.read(training_images, ":")
        self.elements = self.config["dataset"].get(
            "elements", self.get_unique_elements(training_images)
        )

        self.forcetraining = self.config["model"].get("get_forces", True)
        self.fp_scheme = self.config["dataset"].get("fp_scheme", "gaussian").lower()
        self.fp_params = self.config["dataset"]["fp_params"]
        self.save_fps = self.config["dataset"].get("save_fps", True)
        self.cutoff_params = self.config["dataset"].get(
            "cutoff_params", {"cutoff_func": "Cosine"}
        )

        self.train_dataset = AtomsDataset(
            images=training_images,
            descriptor_setup=(
                self.fp_scheme,
                self.fp_params,
                self.cutoff_params,
                self.elements,
            ),
            forcetraining=self.forcetraining,
            save_fps=self.save_fps,
        )

        self.target_scaler = self.train_dataset.target_scaler
        if not self.debug:
            normalizers = {"target": self.target_scaler}
            torch.save(normalizers, os.path.join(self.cp_dir, "normalizers.pt"))
        self.input_dim = self.train_dataset.input_dim
        self.val_split = self.config["dataset"].get("val_split", 0)
        print("Loading dataset: {} images".format(len(self.train_dataset)))

    def load_model(self):
        elements = list_symbols_to_indices(self.elements)
        self.model = BPNN(
            elements=elements, input_dim=self.input_dim, **self.config["model"]
        )
        print("Loading model: {} parameters".format(self.model.num_params))

    def load_extras(self):
        callbacks = []
        load_best_loss = train_end_load_best_loss(self.identifier)
        self.split = CVSplit(cv=self.val_split) if self.val_split != 0 else 0

        metrics = evaluator(
            self.val_split,
            self.config["optim"].get("metric", "mae"),
            self.identifier,
            self.forcetraining,
        )
        callbacks.extend(metrics)

        if not self.debug:
            callbacks.append(load_best_loss)
        scheduler = self.config["optim"].get("scheduler", None)
        if scheduler:
            scheduler = LRScheduler(
                scheduler, **self.config["optim"]["scheduler_params"]
            )
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
        self.callbacks = callbacks

    def load_criterion(self):
        self.criterion = self.config["optim"].get("loss_fn", CustomLoss)

    def load_optimizer(self):
        self.optimizer = self.config["optim"].get("optimizer", torch.optim.Adam)

    def load_logger(self):
        if self.config["cmd"].get("logger", False):
            import wandb

            self.wandb_run = wandb.init(
                name=self.identifier,
                config=self.config,
                id=self.timestamp,
            )

    def load_skorch(self):
        skorch.net.to_tensor = to_tensor

        collate_fn = DataCollater(train=True, forcetraining=self.forcetraining)

        self.net = NeuralNetRegressor(
            module=self.model,
            criterion=self.criterion,
            criterion__force_coefficient=self.config["optim"].get(
                "force_coefficient", 0
            ),
            criterion__loss=self.config["optim"].get("loss", "mse"),
            optimizer=self.optimizer,
            lr=self.config["optim"].get("lr", 1e-1),
            batch_size=self.config["optim"].get("batch_size", 32),
            max_epochs=self.config["optim"].get("epochs", 100),
            iterator_train__collate_fn=collate_fn,
            iterator_train__shuffle=True,
            iterator_valid__collate_fn=collate_fn,
            iterator_valid__shuffle=False,
            device=self.device,
            train_split=self.split,
            callbacks=self.callbacks,
            verbose=self.config["cmd"].get("verbose", True),
        )
        print("Loading skorch trainer")

    def train(self, raw_data=None):
        if raw_data is not None:
            self.config["dataset"]["raw_data"] = raw_data
        if not self.pretrained:
            self.load()

        self.net.fit(self.train_dataset, None)

    def predict(self, images, batch_size=32):
        if len(images) < 1:
            warnings.warn("No images found!", stacklevel=2)
            return images

        a2d = AtomsToData(
            descriptor=self.train_dataset.descriptor,
            r_energy=False,
            r_forces=False,
            save_fps=self.save_fps,
            fprimes=self.forcetraining,
            cores=1,
        )

        data_list = a2d.convert_all(images, disable_tqdm=True)

        self.net.module.eval()
        collate_fn = DataCollater(train=False, forcetraining=self.forcetraining)

        predictions = {"energy": [], "forces": []}
        for data in data_list:
            collated = collate_fn([data])
            energy, forces = self.net.module(collated)

            energy = self.target_scaler.denorm(energy, pred="energy").detach().tolist()
            forces = self.target_scaler.denorm(forces, pred="forces").detach().numpy()

            predictions["energy"].extend(energy)
            predictions["forces"].append(forces)

        return predictions

    def load_pretrained(self, checkpoint_path=None):
        print(f"Loading checkpoint from {checkpoint_path}")
        self.load()
        self.net.initialize()
        self.pretrained = True
        try:
            self.net.load_params(
                f_params=os.path.join(checkpoint_path, "params.pt"),
                f_optimizer=os.path.join(checkpoint_path, "optimizer.pt"),
                f_criterion=os.path.join(checkpoint_path, "criterion.pt"),
                f_history=os.path.join(checkpoint_path, "history.json"),
            )
            # TODO(mshuaibi): remove dataset load, use saved normalizers
        except NotImplementedError:
            print("Unable to load checkpoint!")
