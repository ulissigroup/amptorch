import torch
import torch.nn as nn
from torch.nn import Tanh
from torch.autograd import grad
from torch_scatter import scatter

try:
    from torch_sparse import spmm

    spmm_exists = True
except:
    spmm_exists = False
    pass


class MLP(nn.Module):
    def __init__(
        self, n_input_nodes, n_layers, n_hidden_size, activation, n_output_nodes=1
    ):
        super(MLP, self).__init__()
        if isinstance(n_hidden_size, int):
            n_hidden_size = [n_hidden_size] * (n_layers)
        self.n_neurons = [n_input_nodes] + n_hidden_size + [n_output_nodes]
        self.activation = activation
        layers = []
        for _ in range(n_layers - 1):
            layers.append(nn.Linear(self.n_neurons[_], self.n_neurons[_ + 1]))
            layers.append(activation())
        layers.append(nn.Linear(self.n_neurons[-2], self.n_neurons[-1]))
        self.model_net = nn.Sequential(*layers)

    def forward(self, inputs):
        return self.model_net(inputs)


class ElementMask(nn.Module):
    def __init__(self, elements):
        super(ElementMask, self).__init__()
        self.nelems = len(elements)
        self.mask = nn.Embedding(100, self.nelems)

        weights = torch.zeros(100, self.nelems)
        for idx, Z in enumerate(elements):
            weights[Z, idx] = 1.0
        self.mask.weight.data = weights

    def forward(self, atomic_numbers):
        return self.mask(atomic_numbers)


class BPNN(nn.Module):
    def __init__(
        self,
        unique_atoms,
        architecture,
        device,
        forcetraining,
        activation=Tanh,
        require_grd=True,
    ):
        super(BPNN, self).__init__()
        self.device = device
        self.req_grad = require_grd
        self.forcetraining = forcetraining
        self.architecture = architecture
        self.activation_fn = activation

        input_length = architecture[0]
        n_layers = architecture[1]
        n_hidden_size = architecture[2]
        n_elements = len(unique_atoms)
        self.elementwise_models = nn.ModuleList()
        for element in range(n_elements):
            self.elementwise_models.append(
                MLP(
                    n_input_nodes=input_length,
                    n_layers=n_layers,
                    n_hidden_size=n_hidden_size,
                    activation=activation,
                )
            )

        self.element_mask = ElementMask(unique_atoms)

    def forward(self, batch):
        atomic_numbers = batch.atomic_numbers
        fingerprints = batch.fingerprint.float()
        fingerprints.requires_grad = True
        image_idx = batch.image_idx
        mask = self.element_mask(atomic_numbers)
        o = torch.sum(
            mask * torch.cat([net(fingerprints) for net in self.elementwise_models], 1), dim=1
        )
        energy = scatter(o, image_idx)

        if self.forcetraining:
            gradients = grad(
                energy,
                fingerprints,
                grad_outputs=torch.ones_like(energy),
                create_graph=True,
            )[0].view(1, -1)

            forces = -1 * torch.matmul(gradients, batch.fprimes).view(-1, 3)
        return energy, forces


class CustomMSELoss(nn.Module):
    """Custom loss function to be optimized by the regression. Includes aotmic
    energy and force contributions.

    Eq. (26) in A. Khorshidi, A.A. Peterson /
    Computer Physics Communications 207 (2016) 310-324"""

    def __init__(self, force_coefficient=0):
        super(CustomMSELoss, self).__init__()
        self.alpha = force_coefficient

    def forward(self, prediction, target):

        energy_pred = prediction[0]
        energy_targets_per_atom = target[0]
        num_atoms = target[1]
        MSE_loss = nn.MSELoss(reduction="sum")
        energy_pred_per_atom = torch.div(energy_pred, num_atoms)
        energy_loss = MSE_loss(energy_pred_per_atom, energy_targets_per_atom)

        if self.alpha > 0:
            force_pred = prediction[1]
            if force_pred.nelement() == 0:
                raise Exception("Force training disabled. Set force_coefficient to 0")
            force_targets_per_atom = target[-1]
            num_atoms_extended = torch.cat([idx.repeat(int(idx)) for idx in num_atoms])
            num_atoms_extended = torch.sqrt(num_atoms_extended.reshape(-1, 1))
            force_pred_per_atom = torch.div(force_pred, num_atoms_extended)
            force_targets_per_atom = force_targets_per_atom * num_atoms_extended
            force_loss = (self.alpha / 3) * MSE_loss(
                force_pred_per_atom, force_targets_per_atom
            )
            loss = 0.5 * (energy_loss + force_loss)
        else:
            loss = 0.5 * energy_loss
        return loss
