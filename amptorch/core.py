"""
core.py: AMPTorch defines the core backbone to the training scheme. Model
parameters are defined, a train method carrys out the training, and plotting
methods allow the visualization of the results.
"""

import time
import os
import copy
import torch
import torch.optim as optim
from torch.utils.data import DataLoader
from amp.descriptor.gaussian import Gaussian
from amptorch.utils import Logger
from .model import BPNN, CustomMSELoss
from .data_preprocess import AtomsDataset, collate_amp
from .trainer import Trainer

__author__ = "Muhammed Shuaibi"
__email__ = "mshuaibi@andrew.cmu.edu"


class AMPTorch:
    """Model class used to define the Neural Network architecture and regression
    training scheme

    Parameters
    ----------
    datafile : trajectory file, list, or database
        Training data to be utilized for the regression model.
    device : str
        Hardware to be utilized for training - CPU or GPU.
        default: 'cpu'
    batch_size: int
        Number of images in a training batch. None represents no batching,
        treating the entire dataset as one batch. default: None
    structure: list
        Neural Network architecture. First index represents number of layers,
        including the output layer. Second index represents the number of nodes
        in each hidden layer. i.e. [3,5] = 3 layeres (2 hidden layers, 1 output
        layer) and 5 nodes in each hidden layer. default: [3,5]
    val_frac: float
        Proportion of dataset to be used as a validation test set for training
        purposes. default: 0
    descriptor: object
        Descriptor to be utilized to calculate fingerprints and
        fingerprintprimes. default: Gaussian
    Gs: object
        Symmetry function parameters to be used.
    cores: int
        Number of cores to parallelize across for Symmetry Function computation
    force_coefficient: float
        Define the force coefficient to be utilized in the loss function. A
        coefficient > 0 indicates force training is turned on.
        default: 0
    criterion: object
        Specify the loss function to be optimized.
        default: CustomMSELoss
    optimizer: object
        Define the training optimizer to be utilized for the regression.
        default: optim.LBFGS
    scheduler: object
        Specify whether a learning rate decay scheme is to be utilized during
        training.
        default: None
    lr: float
        Define the model learning rate. default: 1
    criteria: dict
        Define the training convergence criteria, including - rmse criteria,
        early_stop and number of epochs. Training termination is achieved when
        one of the applied settings are met.
        default: {'energy':0, "force":0, "early_stop": False, "epochs": 1e10 }
    lj_data: list
        Energies and forces to be subtracted off from targets, allowing the
        model to learn the difference. default: None
    fine_tune: str
        model parameters to be loaded and used for transfer learning
    label: str
        Label to be used to save model parameters and logs. default: amptorch
    save_logs: boolean
        True to save logs to a .txt file. False to not save logs. default: True

    """

    def __init__(
        self,
        datafile,
        device="cpu",
        structure=[3, 5],
        val_frac=0,
        descriptor=Gaussian,
        Gs=None,
        cores=1,
        force_coefficient=0,
        criterion=CustomMSELoss,
        optimizer=optim.LBFGS,
        loader_params={"batch_size": None, "shuffle": False, "num_workers": 0},
        resample=None,
        scheduler=None,
        lr=1,
        criteria={"energy": 0, "force": 0,
                  "epochs": 1e10, "early_stop": False},
        lj_data=None,
        fine_tune=None,
        label="amptorch",
        save_logs=True,
    ):
        if not os.path.exists("results/logs/epochs"):
            os.makedirs("results/logs/epochs")
        self.save_logs = save_logs
        self.label = label
        self.log = Logger("results/logs/" + label + ".txt")

        self.filename = datafile
        self.device = device
        self.structure = structure
        self.val_frac = val_frac
        self.loader_params = loader_params
        self.resample = resample
        self.descriptor = descriptor
        self.force_coefficient = force_coefficient
        self.criterion = criterion
        self.optimizer = optimizer
        self.scheduler = scheduler
        self.lr = lr
        self.convergence = criteria
        self.lj_data = lj_data
        self.fine_tune = fine_tune
        self.Gs = Gs

        self.forcetraining = False
        if force_coefficient > 0:
            self.forcetraining = True
        if self.lj_data:
            self.scaling = "rel"
        else:
            self.scaling = None

        self.training_data = AtomsDataset(
            self.filename,
            descriptor=self.descriptor,
            Gs=Gs,
            cores=cores,
            forcetraining=self.forcetraining,
            lj_data=self.lj_data,
            label=label,
            scaling=self.scaling,
        )
        self.scalings = self.training_data.scalings
        self.sd_scaling = self.scalings[0]
        self.mean_scaling = self.scalings[1]

        if not lj_data:
            self.log(time.asctime())
            self.log("-" * 50)
        self.log("LJ Data: %s" % (True if lj_data is not None else None))
        self.log("Force Training: %s - %s" %
                 (self.forcetraining, force_coefficient))

    def train(self):
        """Trains the model under the provided optimizer conditions until
        convergence is reached as specified by the rmse_critieria."""

        training_data = self.training_data
        self.unique_atoms = training_data.unique()
        fp_length = training_data.fp_length
        dataset_size = len(training_data)

        if isinstance(self.val_frac, float):
            samplers = training_data.create_splits(
                training_data, self.val_frac, resample=self.resample)
            dataset_size = {
                "train": dataset_size - int(self.val_frac * dataset_size),
                "val": int(self.val_frac * dataset_size),
            }

            self.log(
                "Training Data = %d Validation Data = %d"
                % (dataset_size["train"], dataset_size["val"])
            )

            loader_dict = {}
            if self.loader_params['batch_size'] is None:
                for x in ['train', 'val']:
                    self.loader_params['batch_size'] = dataset_size[x]
                    loader_dict[x] = self.loader_params.copy()
            else:
                for x in ['train', 'val']:
                    loader_dict[x] = self.loader_params
            self.atoms_dataloader = {
                x: DataLoader(
                    training_data,
                    collate_fn=collate_amp,
                    sampler=samplers[x],
                    **loader_dict[x]
                )
                for x in ["train", "val"]
            }
        else:
            self.log("Training Data = %d" % dataset_size)
            if self.loader_params['batch_size'] is None:
                self.loader_params['batch_size'] = dataset_size
            self.atoms_dataloader = DataLoader(
                training_data, collate_fn=collate_amp, **self.loader_params
            )
            if self.loader_params['batch_size'] is None:
                self.loader_params['batch_size'] = len(training_data)
        self.log("Resampled Points = %s" % self.resample)
        architecture = copy.copy(self.structure)
        architecture.insert(0, fp_length)
        model = BPNN(
            self.unique_atoms, architecture, self.device, self.forcetraining
        ).to(self.device)
        if self.fine_tune is not None:
            model.load_state_dict(torch.load(self.fine_tune))
        self.log("Activation Function: %s" % model.activation_fn)
        self.log("Loss Function: %s" % self.criterion)
        # Define the optimizer and implement any optimization settings
        if self.optimizer == optim.LBFGS:
            optimizer = self.optimizer(
                model.parameters(), self.lr, line_search_fn="strong_wolfe"
            )
        else:
            optimizer = self.optimizer(model.parameters(), self.lr)
        self.log("Optimizer Info:\n %s" % optimizer)

        if self.scheduler:
            self.scheduler = self.scheduler(optimizer, 5)
        self.log("Scheduler Info: %s" % self.scheduler)
        self.log("Convergence criteria = {}".format(self.convergence))
        self.log("Model architecture: %s" % architecture)

        self.trainer = Trainer(
            model,
            self.device,
            dataset_size,
            self.criterion(force_coefficient=self.force_coefficient),
            optimizer,
            self.scheduler,
            self.atoms_dataloader,
            self.convergence,
            self.scalings,
            self.label,
        )
        self.trained_model, value = self.trainer.train_model()
        if not self.save_logs:
            os.remove("results/logs/" + self.label + ".txt")
            os.remove("results/logs/epochs/" + self.label + "-calc.txt")
        return self.trained_model, value
