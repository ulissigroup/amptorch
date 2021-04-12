[![ulissigroup](https://circleci.com/gh/ulissigroup/amptorch.svg?style=svg)](https://app.circleci.com/pipelines/github/ulissigroup/amptorch)
## *AMPtorch*: Atomistic Machine-learning Package - PyTorch

*AMPtorch* is a PyTorch implementation of the [Atomistic Machine-learning Package](https://amp.readthedocs.io/en/latest/) (AMP) code that seeks to provide users with improved performance and flexibility as compared to the original code. The implementation does so by benefiting from state-of-the-art machine learning methods and techniques to be optimized in conjunction with high-throughput supercomputers. *AMPtorch* is built on top of [PyTorch Geometric](https://pytorch-geometric.readthedocs.io/en/latest/) and [Skorch](https://skorch.readthedocs.io/en/stable/).

### Installation

Install dependencies:

1. Ensure conda is up-to-date: ```conda update conda```

2. Create conda environment
- CPU machines: ```conda env create -f env_cpu.yml```
- GPU machines (CUDA 10.2): ```conda env create -f env_gpu.yml```

3. Activate the conda environment `conda activate amptorch` and install the package with `pip install -e .`

4. Install pre-commit hooks: `pre-commit install`

### Usage
#### Configs
To train a model using `amptorch`, a set of `configs` must be specified to interact with the trainer. An exhaustive list of all possible flags and their descriptions is provided below:
```
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
      "dtype": object,              # Pytorch level of precision (default: torch.FloatTensor)
      "run_dir": str,               # Path to run trainer, where logs are to be saved (default: "./")
      "seed": int,                  # Random seed (default: 0)
      "identifier": str,            # Unique identifer to experiment, optional
      "verbose": bool,              # Print training scores (default: True)
      "logger": False,              # Log results to Weights and Biases (https://www.wandb.com/)
                                    ## wandb offers a very clean and flexible interface to monitor results online
                                    ## A free account is necessary to view and log results
  },
}
```
#### Train model
```
from amptorch import AtomsTrainer

trainer = AtomsTrainer(configs)
trainer.train()
```
#### Load checkpoints
Previously trained models may be loaded as follows:
```
trainer = AtomsTrainer(configs)
trainer.load_pretrained(path_to_checkpoint_dir)
```
#### Make predictions
```
predictions = trainer.predict(list_of_atoms_objects)
energies = predictions["energy"]
forces = predictions["forces"]
```
#### Construct AMPtorch-ASE calculator
To interface with ASE, an ASE calculator may be constructed as follows:
```
from amptorch import AMPtorch

calc = AMPtorch(trainer)
slab.set_calculator(calc)
energy = slab.get_potential_energy()
forces = slab.get_forces()
```

### Acknowledgements
- This project is being developed at Carnegie Mellon University in the Department of Chemical Engineering, by Muhammed Shuaibi and Zachary Ulissi, in collaboration with Andrew Peterson, Franklin Goldsmith, Brenda Rubenstein, Andrew Medford, and Adam Willard as part of the Department of Energy's *Bridging the time scale in exascale computing of chemical systems* project. AMPtorch developers include Xiangyun Lei, Ben Comer, Rui Qi Chen, Eric Musa, and Matt Adams.
- Funded by the Department of Energy's Basic Enenergy Science, Computational Chemical Sciences Program Office. Award # DE-SC0019441
- Engineering ideas have been heavily borrowed from our work on the [Open Catalyst Project](https://github.com/Open-Catalyst-Project/baselines)
- Gaussian fingerprints have been adapted from [SIMPLE-NN](https://github.com/MDIL-SNU/SIMPLE-NN)
