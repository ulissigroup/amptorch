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
from .utils import Logger
from .NN_model import FullNN, CustomLoss
from .data_preprocess import AtomsDataset, collate_amp
from .weighted_data_preprocess import WeightedAtomsDataset, weighted_collate_amp
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
    cores : int
        Specify the number of cores to use for parallelization of fingerprint
        calculations. NOTE - To be replaced with PyTorch's parallelization
        scheme upon fingerprint implementation
    envommand: string
        For parallel processing across nodes, a command can be supplied here to
        load the appropriate environment before starting workers.
        default=None. NOTE - See Above.
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
    force_coefficient: float
        Define the force coefficient to be utilized in the loss function. A
        coefficient > 0 indicates force training is turned on.
        default: 0
    criterion: object
        Specify the loss function to be optimized.
        default: CustomLoss
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
        Define the training convergence criteria.
        default: {'energy':0, "force":0}
    epochs: int
        If present, train the model for this number of epochs. If epochs or
        criteria are not defined the model will train until error stagnates.
        default: 1e10
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
        cores=1,
        envcommand=None,
        batch_size=None,
        structure=[3, 5],
        val_frac=0,
        descriptor=Gaussian,
        Gs=None,
        force_coefficient=0,
        criterion=CustomLoss,
        optimizer=optim.LBFGS,
        scheduler=None,
        lr=1,
        criteria={"energy": 0, "force": 0},
        epochs=1e10,
        lj_data=None,
        fine_tune=None,
        label='amptorch',
        save_logs=True,
        save_interval=1000,
        store_primes=False,
        weights_dict=None
    ):
        if not os.path.exists("results/logs"):
            os.makedirs("results/logs/epochs")
        self.save_logs = save_logs
        self.save_interval = save_interval
        self.label = label
        self.log = Logger("results/logs/"+label+".txt")

        self.filename = datafile
        self.batch_size = batch_size
        self.device = device
        self.batch_size = batch_size
        self.structure = structure
        self.val_frac = val_frac
        self.descriptor = descriptor
        self.force_coefficient = force_coefficient
        self.criterion = criterion
        self.optimizer = optimizer
        self.scheduler = scheduler
        self.lr = lr
        self.convergence = criteria
        self.epochs = epochs
        self.lj_data = lj_data
        self.fine_tune = fine_tune
        self.Gs = Gs
        self.store_primes = store_primes
        self.weights_dict = weights_dict

        self.forcetraining = False
        if force_coefficient > 0:
            self.forcetraining = True
        self.weighted = (self.weights_dict is not None)

        if not self.weighted:
            self.training_data = AtomsDataset(
                self.filename,
                descriptor=self.descriptor,
                Gs=Gs,
                cores=cores,
                forcetraining=self.forcetraining,
                lj_data=self.lj_data,
                envcommand=envcommand,
                store_primes=store_primes,)
        else:
            self.training_data = WeightedAtomsDataset(
                self.filename,
                descriptor=self.descriptor,
                Gs=Gs,
                cores=cores,
                forcetraining=self.forcetraining,
                lj_data=self.lj_data,
                envcommand=envcommand,
                store_primes=store_primes,
                weights_dict=weights_dict)
        self.scalings = self.training_data.scalings()
        self.sd_scaling = self.scalings[0]
        self.mean_scaling = self.scalings[1]

        if not lj_data:
            self.log(time.asctime())
            self.log("-" * 50)
        self.log("LJ Data: %s" % (True if lj_data is not None else None))
        self.log("Force Training: %s - %s" % (self.forcetraining, force_coefficient))

    def train(self):
        """Trains the model under the provided optimizer conditions until
        convergence is reached as specified by the rmse_critieria."""

        training_data = self.training_data
        self.unique_atoms = training_data.unique()
        fp_length = training_data.fp_length
        dataset_size = len(training_data)
        if self.weighted:
            collate_function = weighted_collate_amp
        else:
            collate_function = collate_amp


        if self.batch_size is None:
            self.batch_size = dataset_size

        if self.val_frac != 0:
            samplers = training_data.create_splits(training_data, self.val_frac)
            dataset_size = {
                "train": dataset_size - int(self.val_frac * dataset_size),
                "val": int(self.val_frac * dataset_size),
            }

            self.log(
                "Training Data = %d Validation Data = %d"
                % (dataset_size["train"], dataset_size["val"])
            )

            if self.batch_size is None:
                self.batch_size = {x: dataset_size[x] for x in ["train", "val"]}

            self.atoms_dataloader = {
                x: DataLoader(
                    training_data,
                    self.batch_size,
                    collate_fn=collate_function,
                    sampler=samplers[x],
                )
                for x in ["train", "val"]
            }

        else:
            self.log("Training Data = %d" % dataset_size)
            self.atoms_dataloader = DataLoader(
                training_data, self.batch_size, collate_fn=collate_function, shuffle=False
            )
        architecture = copy.copy(self.structure)
        architecture.insert(0, fp_length)
        model = FullNN(
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
            self.scheduler = self.scheduler(optimizer, step_size=7, gamma=0.1)
        self.log("Scheduler Info: %s" % self.scheduler)
        self.log("RMSE criteria = {}".format(self.convergence))
        self.log("Model architecture: %s" % architecture)

        self.trainer = Trainer(
            model,
            self.device,
            self.unique_atoms,
            dataset_size,
            self.criterion(force_coefficient=self.force_coefficient),
            optimizer,
            self.scheduler,
            self.atoms_dataloader,
            self.convergence,
            self.epochs,
            self.scalings,
            self.label,
            weighted=self.weighted,
            save_interval=self.save_interval,
        )

        self.trained_model = self.trainer.train_model()
        if not self.save_logs:
            os.remove("results/logs/"+self.label+".txt")
            os.remove("results/logs/epochs/"+self.label+"-calc.txt")
        return self.trained_model
