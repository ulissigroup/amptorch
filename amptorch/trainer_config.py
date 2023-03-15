from typing import List, Dict, Union, Optional
from dataclasses import dataclass
import os
import ase
import torch


@dataclass
class Config_Model:
    """
    Configuration to define the atomistic neural network model configuration.

    Attributes:
    ------------
        name (str): The model to be used for atomistic neural network force field, "SingleNN" or "BPNN". (default:  "SingleNN")

        num_layers (int): No. of hidden layers

        num_nodes (int): No. of nodes per layer

        get_forces (bool): Whether to compute per-atom forces for force training. Set "force_coefficient" in "Optim" accordingly. (default: True)

        batchnorm (bool): Enable batch-normalization (default:False)

        activation (object): Activation function. Any activation supported by torch. (default: nn.Tanh)

        **custom_args: Any additional arguments used to customize existing/new models
    """

    name: str
    num_layers: int
    num_nodes: int
    get_forces: bool = True
    batchnorm: bool = True
    activation: object = None
    custom_args: Dict[str, Union[int, bool, object]] = None


@dataclass
class Config_Optim:
    """
    Configuration to define the resources, and setting for neural network optimization.

    Attributes:
    ------------
    gpus (int): No. of gpus to use, 0 for cpu (default: 0)

    force_coefficient (float): If force training, coefficient to weight the force component by (default: 0)

    lr (float): Initial learning rate (default: 1e-1)

    batch_size (int): Batch size (default: 32)

    epochs (int): Max training epochs (default: 100)

    optimizer (object): Training optimizer (default: torch.optim.Adam)

    loss_fn (object): Loss function to optimize (default: CustomLoss)

    loss (str): Control loss function criterion, "mse" or "mae" (default: "mse")

    metric (str): Metrics to be reported by, "mse" or "mae" (default: "mae")

    cp_metric (str): Property based on which the model is saved. "energy" or "forces" (default: "energy")

    scheduler (dict): Learning rate scheduler to use
        {"policy": "StepLR", "params": {"step_size": 10, "gamma": 0.1}}
    """

    gpus: int = 0
    force_coefficient: float = 0
    lr: float = 1e-1
    batch_size: int = 32
    epochs: int = 100
    optimizer: object = torch.optim.Adam
    loss: str = "mse"
    metric: str = "mae"
    cp_metric: str = "energy"
    scheduler: object = None


@dataclass
class Config_Dataset:
    """
    Configuration to define the dataset used for training, and featurization scheme to be used. Featurization default to Guassian Multipole Fingerprinting.

    Attributes:
    ------------
    raw_data (str or List[ase.Atoms]]): Path to ASE trajectory or database, or list of Atoms objects.

    lmdb_path (Optional: str): Path to LMDB database file for dataset too large to fit in memory, if raw_data is not provided.

    val_split (float): Proportion of training set to use for validation.

    elements (Optional: List[str]): List of unique elements in dataset, optional. Example: ["H", "O"] for water dataset.

    fp_scheme (str): Fingerprinting scheme to feature dataset, "gmpordernorm" or "gaussian" (default: "gmpordernorm").

    fp_params (dict): Fingerprint parameters, dataclass "GMP_params" or "SF_params".

    cutoff_params (Optional, dict): Cutoff function for Gaussian fingerprinting scheme only - polynomial or cosine,
        polynomial - {"cutoff_func": "Polynomial", "gamma": 2.0},
        cosine - {"cutoff_func": "Cosine"}.

    save_fps (bool): Write calculated fingerprints to disk (default: True).

    scaling (Optional[Dict[str, Any]]): Feature scaling scheme, normalization or standardization,
        normalization (scales features between "range") - {"type": "normalize", "range": (0, 1)},
        standardization (scales data to mean=0, stdev=1) - {"type": "standardize"}.
    """

    raw_data: List
    lmdb_path: str = None
    val_split: float = None
    elements: Optional[List[str]] = None
    fp_scheme: str = "gmpordernorm"
    # fp_params:
    cutoff_params: Optional[object] = None
    save_fps: bool = True
    scaling: Optional[dict]


@dataclass
class GMPs_params:
    """
    Dataclass to represent the Gaussian Momenta Parameters (GMPs) dictionary object as input to Config class fp_params.

    Attributes:
    -----------
    orders : List[int]
        List of integer orders of the spherical harmonics functions. For example, [0, 1, 2, 3] would account up to order 3 from order 0.
    sigmas : List[float]
        List of float values of the widths of the Gaussian basis functions.
    atom_gaussians : Dict[str], optional
        Dictionary mapping element symbols to paths of the corresponding pseudodensity files containing Gaussian basis functions for each element.
    square : bool, optional
        Boolean flag indicating whether to square the solid harmonics functions for smoothness in gradient calculation.
        Default is True.
    solid_harmonics : bool, optional
        Boolean flag indicating whether to use solid harmonics functions instead of real spherical harmonics. Default is True.
    """

    orders: List[int]
    sigmas: List[float]
    atom_gaussians: Dict
    square: bool = True
    solid_harmonics: bool = True

    # def load_atom_gaussians(self, directory: str) -> None:
    #     for element in self.atom_gaussians:
    #         filename = os.path.join(directory, f"{element}_pseudodensity.g")
    #         self.atom_gaussians[element] = filename


@dataclass
class G2Params:
    """
    Dataclass representing parameters for the G2 component of the gaussian symmetry function.

    Attributes:
    -----------
    etas: List[float]
        A list of the exponents used in the G2 function.
    rs_s: List[float]
        A list of radial offsets used in the G2 function.
    """

    etas: List[float]
    rs_s: List[float]


@dataclass
class G4Params:
    """
    Dataclass representing parameters for the G2 component of the gaussian symmetry function.

    Attributes:
    -----------
    etas: List[float]
        A list of exponents used in the G4 function.
    zetas: List[float]
        A list of zetas used in the G4 function.
    gammas: List[float]
        A list of gammas used in the G4 function.
    """

    etas: List[float]
    zetas: List[float]
    gammas: List[float]


@dataclass
class Gaussian_params:
    """
    Dataclass for Gaussian Symmetry Function with G2 and G4 with a cutoff.

    Attributes:
    -----------
    gaussian : Dict
        A dictionary of Gaussian types (G2 and G4) and their parameters.

    cutoff : float
        The cutoff value.
    """

    gaussian: dict
    cutoff: int

    def __init__(self, gaussian, cutoff):
        self.gaussian = gaussian
        self.cutoff = cutoff

    @classmethod
    def from_dict(cls, params_dict):
        gaussian_params = params_dict.get("gaussian", {})
        g2_params_dict = gaussian_params.get("G2", {})
        g4_params_dict = gaussian_params.get("G4", {})
        g2_params = G2Params(
            etas=g2_params_dict.get("etas", []), rs_s=g2_params_dict.get("rs_s", [])
        )
        g4_params = G4Params(
            etas=g4_params_dict.get("etas", []),
            zetas=g4_params_dict.get("zetas", []),
            gammas=g4_params_dict.get("gammas", []),
        )
        return cls(
            gaussian={"G2": g2_params, "G4": g4_params},
            cutoff=params_dict.get("cutoff", 6),
        )


@dataclass
class Config_Cmd:
    """
    Configuration for the extra commands and idenfiers.

    Attributes:
    -----------
        debug : Optional, bool
            If True, enables debug mode, which does not write/save checkpoints/results. Defaults to False.

        dtype : Optional, type
            PyTorch level of precision. Defaults to torch.DoubleTensor.

        run_dir : Optional, str
            Path to the directory where logs are to be saved. Defaults to "./".

        seed : Optional, int
            Random seed to use. If not specified, a random seed of 0 will be used to ensure consistency.

        identifier : Optional, str
            Unique identifier for the experiment. If not specified, the current time will be used.

        verbose : Optional, bool
            If True, print training scores. Defaults to True.

        logger : Optional, bool
            If True, log results to Weights and Biases (https://www.wandb.com/). A free account is necessary to view and log results.
    """

    debug: bool = False
    dtype: Optional(object)
    run_dir: str = "./"
    seed: int = 0
    identifier: str = None
    verbose: bool = True
    logger: bool = False


@dataclass
class Config:
    """
    Input configuration for trainer class to start the training of NNFF.

    Attributes:
    -----------
    model: Model dataclass object or dictionary. Specifying the configuration of neural network force field model.

    optim: Optim dataclass object or dictionary. Specifying the optimizer and resources for optimization.

    dataset: Dataset dataclass object or dictionary. Specifying the dataset used for training and validation, and the fingerprinting scheme.

    cmd: Extra commands for debugging, running directory, identifier.
    """

    model: Config_Model
    optim: Config_Optim
    dataset: Config_Dataset
    cmd: Config_Cmd
