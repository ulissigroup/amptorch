import torch
import torch.nn as nn
from torch.autograd import grad
from torch.nn import Tanh
from torch_scatter import scatter


class MLP(nn.Module):
    """
    Multi-layer perceptron model modified for atomistic input.

    Args:
        n_input_nodes (int): Number of input nodes for the network.
        n_layers (int): Number of hidden layers in the network.
        n_hidden_size (int): Number of hidden units per layer.
        activation (torch.nn.Module): Activation function to use in each layer.
        batchnorm (bool): Whether to use batch normalization after each layer.
        dropout (bool): Whether to use dropout after each layer.
        dropout_rate (float): Dropout rate to use if `dropout` is True.
        hidden_layers (Optional[List[int]]): List of hidden layer sizes. If not None,
            `n_layers` and `n_hidden_size` will be ignored.
        n_output_nodes (int): Number of output nodes for the network.
        initialization (str): Initialization method for the network weights. "xavier" or "zero".
    """

    def __init__(
        self,
        n_input_nodes,
        n_layers,
        n_hidden_size,
        activation,
        batchnorm,
        dropout,
        dropout_rate,
        hidden_layers=None,
        n_output_nodes=1,
        initialization="xavier",
    ):
        super(MLP, self).__init__()
        if hidden_layers is None and isinstance(n_hidden_size, int):
            hidden_layers = [n_hidden_size] * (n_layers)
        else:
            n_layers = len(hidden_layers)
        self.n_neurons = [n_input_nodes] + hidden_layers + [n_output_nodes]
        self.activation = activation
        layers = []
        for _ in range(n_layers):
            layers.append(nn.Linear(self.n_neurons[_], self.n_neurons[_ + 1]))
            if batchnorm:
                layers.append(nn.BatchNorm1d(self.n_neurons[_ + 1]))
            layers.append(activation())
            if dropout:
                layers.append(nn.Dropout(p=dropout_rate))

        layers.append(nn.Linear(self.n_neurons[-2], self.n_neurons[-1]))
        self.model_net = nn.Sequential(*layers)

        print(torch.get_default_dtype())
        if torch.get_default_dtype() == torch.float64:
            self.double()
        # TODO: identify optimal initialization scheme
        self.reset_parameters(initialization.lower())
        # print(self.model_net)

    def reset_parameters(self, initialization):
        if initialization == "xavier":
            print("Use Xavier initialization")
            for m in self.model_net:
                if isinstance(m, torch.nn.Linear):
                    torch.nn.init.xavier_uniform_(m.weight)
                    m.bias.data.fill_(0)
        elif initialization == "zero":
            print("Use constant zero initialization")
            for m in self.model_net:
                if isinstance(m, torch.nn.Linear):
                    torch.nn.init.constant_(m.weight, 0.0)
                    m.bias.data.fill_(0)
        else:
            print("Warning: unrecognized initialization, use default")

    def forward(self, inputs):
        return self.model_net(inputs)


class ElementMask(nn.Module):
    """
    Mask for different chemical element types for BPNN.

    Args:
        elements (List[str]) : a list of strings of unique chemical elements in the system.
    """

    def __init__(self, elements):
        super(ElementMask, self).__init__()
        nelems = len(elements)
        weights = torch.zeros(100, nelems)
        weights[elements, range(nelems)] = 1.0

        self.mask = nn.Embedding(100, nelems)
        self.mask.weight.data = weights

    def forward(self, atomic_numbers):
        return self.mask(atomic_numbers)


class BPNN(nn.Module):
    """
    Atomistic neural network structure described as 2nd generation or Behler-Parrinello neural network for energy (and force) training.

    Args:
    elements : list of str
        List of unique element symbols in the system.
    input_dim : int
        Dimensionality of the input. The dimension depends on the atomistic fingerprinting scheme.
    num_nodes : int, optional (default=20)
        Number of nodes in each hidden layer.
    num_layers : int, optional (default=5)
        Number of hidden layers in the network.
    hidden_layers : list of int, optional (default=None)
        A list of integers, where each element corresponds to the number of nodes in a hidden layer. Overrides num_nodes and num_layers. E.g. [10, 10, 10]
    get_forces : bool, optional (default=True)
        Whether to train with the forces in addition to the energy.
    batchnorm : bool, optional (default=False)
        Whether to se batch normalization in the network.
    dropout : bool, optional (default=False)
        Whether to use to apply dropout in the network.
    dropout_rate : float, optional (default=0.5)
        The dropout probability in [0, 1].
    activation : torch.nn.Module, optional (default=Tanh)
        The activation function to use in the network.
    name : str, optional (default='bpnn')
        Name of the network.
    initialization : str, optional (default='xavier')
        Initialization method to use for weights in the network.

    """

    def __init__(
        self,
        elements,
        input_dim,
        num_nodes=20,
        num_layers=5,
        hidden_layers=None,
        get_forces=True,
        batchnorm=False,
        dropout=False,
        dropout_rate=0.5,
        activation=Tanh,
        name="bpnn",
        initialization="xavier",
    ):
        super(BPNN, self).__init__()
        self.get_forces = get_forces
        self.activation_fn = activation

        n_elements = len(elements)
        self.elementwise_models = nn.ModuleList()
        for element in range(n_elements):
            self.elementwise_models.append(
                MLP(
                    n_input_nodes=input_dim,
                    n_layers=num_layers,
                    n_hidden_size=num_nodes,
                    hidden_layers=hidden_layers,
                    activation=activation,
                    batchnorm=batchnorm,
                    dropout=dropout,
                    dropout_rate=dropout_rate,
                    initialization=initialization,
                )
            )

        self.element_mask = ElementMask(elements)

    def forward(self, batch):
        if isinstance(batch, list):
            batch = batch[0]
        with torch.enable_grad():
            atomic_numbers = batch.atomic_numbers
            fingerprints = batch.fingerprint
            fingerprints.requires_grad = True
            image_idx = batch.batch
            mask = self.element_mask(atomic_numbers)
            o = torch.sum(
                mask
                * torch.cat(
                    [net(fingerprints) for net in self.elementwise_models], dim=1
                ),
                dim=1,
            )
            energy = scatter(o, image_idx, dim=0)

            if self.get_forces:
                gradients = grad(
                    energy,
                    fingerprints,
                    grad_outputs=torch.ones_like(energy),
                    create_graph=True,
                )[0].view(1, -1)

                forces = -1 * torch.sparse.mm(batch.fprimes.t(), gradients.t()).view(
                    -1, 3
                )

            else:
                forces = torch.tensor([], device=energy.device)

            return energy, forces

    @property
    def num_params(self):
        return sum(p.numel() for p in self.parameters())


class SingleNN(nn.Module):
    """
    A modified version of Behler-Parrinello atomistic neural network where all elements shared the same  for energy (and force) training.

    Args:
    elements : list of str
        List of unique element symbols in the system.
    input_dim : int
        Dimensionality of the input. The dimension depends on the atomistic fingerprinting scheme.
    num_nodes : int, optional (default=20)
        Number of nodes in each hidden layer.
    num_layers : int, optional (default=5)
        Number of hidden layers in the network.
    hidden_layers : list of int, optional (default=None)
        A list of integers, where each element corresponds to the number of nodes in a hidden layer. Overrides num_nodes and num_layers. E.g. [10, 10, 10]
    get_forces : bool, optional (default=True)
        Whether to train with the forces in addition to the energy.
    batchnorm : bool, optional (default=False)
        Whether to se batch normalization in the network.
    dropout : bool, optional (default=False)
        Whether to use to apply dropout in the network.
    dropout_rate : float, optional (default=0.5)
        The dropout probability in [0, 1].
    activation : torch.nn.Module, optional (default=Tanh)
        The activation function to use in the network.
    name : str, optional (default='singlenn')
        Name of the network.
    initialization : str, optional (default='xavier')
        Initialization method to use for weights in the network.
    """

    def __init__(
        self,
        elements,
        input_dim,
        num_nodes=20,
        num_layers=5,
        hidden_layers=None,
        get_forces=True,
        batchnorm=False,
        dropout=False,
        dropout_rate=0.5,
        activation=Tanh,
        name="singlenn",
        initialization="xavier",
    ):
        super(SingleNN, self).__init__()
        self.get_forces = get_forces
        self.activation_fn = activation

        self.model = MLP(
            n_input_nodes=input_dim,
            n_layers=num_layers,
            n_hidden_size=num_nodes,
            hidden_layers=hidden_layers,
            activation=activation,
            batchnorm=batchnorm,
            dropout=dropout,
            dropout_rate=dropout_rate,
            initialization=initialization,
        )

    def forward(self, batch):
        if isinstance(batch, list):
            batch = batch[0]
        with torch.enable_grad():
            fingerprints = batch.fingerprint
            fingerprints.requires_grad = True
            image_idx = batch.batch
            sorted_image_idx = torch.unique_consecutive(image_idx)
            o = torch.sum(self.model(fingerprints), dim=1)
            energy = scatter(o, image_idx, dim=0)

            if self.get_forces:
                gradients = grad(
                    energy,
                    fingerprints,
                    grad_outputs=torch.ones_like(energy),
                    create_graph=True,
                )[0].view(1, -1)

                forces = -1 * torch.sparse.mm(batch.fprimes.t(), gradients.t()).view(
                    -1, 3
                )

            else:
                forces = torch.tensor([], device=energy.device)

            return energy, forces

    @property
    def num_params(self):
        return sum(p.numel() for p in self.parameters())


class CustomLoss(nn.Module):
    """
    Customize the loss function based on Parrinello's publication with alpha as the force coefficient.
    """

    def __init__(self, force_coefficient=0, loss="mae"):
        super(CustomLoss, self).__init__()
        self.alpha = force_coefficient
        self.loss = loss

        if self.loss == "mae":
            self.loss = nn.L1Loss()
        elif self.loss == "mse":
            self.loss = nn.MSELoss()
        else:
            raise NotImplementedError(f"{self.loss} loss not available!")

    def forward(self, prediction, target):
        energy_pred = prediction[0]
        energy_target = target[0]
        energy_loss = self.loss(energy_pred, energy_target)
        force_pred = prediction[1]
        if force_pred.nelement() == 0:
            self.alpha = 0

        if self.alpha > 0:
            force_target = target[1]
            force_loss = self.loss(force_pred, force_target)
            loss = 0.5 * (energy_loss + self.alpha * force_loss)
        else:
            loss = 0.5 * energy_loss
        return loss
