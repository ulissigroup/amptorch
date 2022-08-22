.. _usage:

==================================
Usage
==================================


Configs
^^^^^^^

To train a model using ``amptorch``, a set of ``configs`` must be
specified to interact with the trainer. An exhaustive list of all
possible flags and their descriptions is provided below:

::

   configs = {
     "model": {
         "num_layers": int,            # No. of hidden layers
         "num_nodes": int,             # No. of nodes per layer
         "get_forces": bool,           # Compute per-atom forces (default: True)
         "batchnorm": bool,            # Enable batch-normalization (default:False)
         "activation": object,         # Activation function (default: nn.Tanh)
         **custom_args                 # Any additional arguments used to customize existing/new models
     },
     "optim": {
         "gpus": int,                  # No. of gpus to use, 0 for cpu (default: 0)
         "force_coefficient": float,   # If force training, coefficient to weight the force component by (default: 0)
         "lr": float,                  # Initial learning rate (default: 1e-1)
         "batch_size": int,            # Batch size (default: 32)
         "epochs": int,                # Max training epochs (default: 100)
         "optimizer": object,          # Training optimizer (default: torch.optim.Adam)
         "loss_fn": object,            # Loss function to optimize (default: CustomLoss)
         "loss": str,                  # Control loss function criterion, "mse" or "mae" (default: "mse")
         "metric": str,                # Metrics to be reported by, "mse" or "mae" (default: "mae")
         "cp_metric": str,             # Property based on which the model is saved. "energy" or "forces" (default: "energy")
         "scheduler": dict,            # Learning rate scheduler to use
                       ##            - {"policy": "StepLR", "params": {"step_size": 10, "gamma": 0.1}}
     },
     "dataset": {
         "raw_data": str or list,      # Path to ASE trajectory or database or list of Atoms objects
         "lmdb_path": str,             # Path to LMDB database file for dataset too large to fit in memory
                           ## Specify either "raw_data" or "lmdb_path"
                       ## LMDB construction can be found in examples/construct_lmdb.py
         "val_split": float,           # Proportion of training set to use for validation
         "elements": list,             # List of unique elements in dataset, optional (default: computes unique elements)
         "fp_scheme": str,             # Fingerprinting scheme to feature dataset, "gaussian" or "gmp" (default: "gaussian")
         "fp_params": dict,            # Fingerprint parameters, see examples for correct layout
         "cutoff_params": dict,        # Cutoff function - polynomial or cosine,
                                       ## Polynomial - {"cutoff_func": "Polynomial", "gamma": 2.0}
                                       ## Cosine     - {"cutoff_func": "Cosine"}
         "save_fps": bool,             # Write calculated fingerprints to disk (default: True)
         "scaling": dict,              # Feature scaling scheme, normalization or standardization
                                       ## normalization (scales features between "range")
                                                     - {"type": "normalize", "range": (0, 1)}
                                       ## standardization (scales data to mean=0, stdev=1)
                                                     - {"type": "standardize"}
     },
     "cmd": {
         "debug": bool,                # Debug mode, does not write/save checkpoints/results (default: False)
         "dtype": object,              # Pytorch level of precision (default: torch.DoubleTensor)
         "run_dir": str,               # Path to run trainer, where logs are to be saved (default: "./")
         "seed": int,                  # Random seed (default: 0)
         "identifier": str,            # Unique identifer to experiment, optional
         "verbose": bool,              # Print training scores (default: True)
         "logger": False,              # Log results to Weights and Biases (https://www.wandb.com/)
                                       ## wandb offers a very clean and flexible interface to monitor results online
                                       ## A free account is necessary to view and log results
     },
   }

Train model
^^^^^^^^^^^

::

   from amptorch import AtomsTrainer

   trainer = AtomsTrainer(configs)
   trainer.train()

Load checkpoints
^^^^^^^^^^^^^^^^

Previously trained models may be loaded as follows:

::

   trainer = AtomsTrainer(configs)
   trainer.load_pretrained(path_to_checkpoint_dir)

Make predictions
^^^^^^^^^^^^^^^^

::

   predictions = trainer.predict(list_of_atoms_objects)
   energies = predictions["energy"]
   forces = predictions["forces"]

Construct AmpTorch-ASE calculator
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To interface with ASE, an ASE calculator may be constructed as follows:

::

   from amptorch import AmpTorch

   calc = AmpTorch(trainer)
   slab.set_calculator(calc)
   energy = slab.get_potential_energy()
   forces = slab.get_forces()
