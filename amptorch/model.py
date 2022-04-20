import torch
import torch.nn as nn
from torch.autograd import grad
from torch.nn import Tanh
from torch_scatter import scatter


class MLP(nn.Module):
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
            print("============")
            print(fingerprints)
            weights = batch.weights
            fingerprints.requires_grad = True
            image_idx = batch.batch
            sorted_image_idx = torch.unique_consecutive(image_idx)
            model_out = torch.sum(self.model(fingerprints), dim=1)
            weighted_model_out = torch.mul(model_out, weights)
            energy = scatter(weighted_model_out, image_idx, dim=0, reduce="sum")

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


class VWN(nn.Module):
    def __init__(self):

        super().__init__()
        # initialize parameters with original parameters

        vwn_parameters = torch.tensor(
            [-0.45816529, 0.031091, 0.21370, 7.5957, 3.5876, 1.6382, 0.49294]
        )
        self.vwn_parameters = nn.Parameter(vwn_parameters)

    def forward(self, inputs):
        """
        VWN formalism
        """

        rho = inputs[:, 0]
        print(inputs)
        print(rho)

        C1, gamma, alpha1, beta1, beta2, beta3, beta4 = self.vwn_parameters
        C0I = 0.23873241

        rs = torch.pow(torch.div(C0I, rho), (1.0 / 3.0))
        rtrs = torch.pow(rs, 0.5)

        ex = C1 / rs

        Q0 = -2.0 * gamma * (1.0 + alpha1 * rs)
        Q1 = (
            2.0
            * gamma
            * rtrs
            * (beta1 + rtrs * (beta2 + rtrs * (beta3 + rtrs * beta4)))
        )
        ec = torch.mul(Q0, torch.log(1.0 + torch.div(1.0, Q1)))

        return torch.mul(rho, torch.add(ex, ec))


class SingleNN_deltaLDA(nn.Module):
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
        name="singlenn_deltalda",
        initialization="xavier",
        LDA_type="vwn",
    ):
        super(SingleNN_deltaLDA, self).__init__()
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

        if LDA_type == "vwn":
            self.LDA = VWN()
        else:
            raise NotImplementedError

    def forward(self, batch):
        if isinstance(batch, list):
            batch = batch[0]
        with torch.enable_grad():
            fingerprints = batch.fingerprint
            weights = batch.weights
            fingerprints.requires_grad = True
            image_idx = batch.batch
            sorted_image_idx = torch.unique_consecutive(image_idx)
            
            NN_model_out = torch.sum(self.model(fingerprints), dim=1)

            LDA_out = self.LDA(fingerprints)

            model_out = torch.add(NN_model_out, LDA_out)
            weighted_model_out = torch.mul(model_out, weights)
            energy = scatter(weighted_model_out, image_idx, dim=0, reduce="sum")

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
