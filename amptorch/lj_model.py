import os
import time
from itertools import product
import numpy as np
from scipy.optimize import minimize
from ase.neighborlist import NeighborList, NewPrimitiveNeighborList
from .gaussian import NeighborlistCalculator, Data
from .utils import Logger, hash_images, get_hash
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


class lj_optim:
    def __init__(self, images, params, params_dict, cutoff, filename, fitforces=True):
        if not os.path.exists("results"):
            os.mkdir("results")
        if not os.path.exists("results/logs"):
            os.mkdir("results/logs")
        self.filename = filename
        self.data = images
        self.p0 = params
        self.params_dict = params_dict
        self.fitforces = fitforces
        self.cutoff = cutoff
        self.hashed_images = hash_images(images)
        self.hashed_keys = list(self.hashed_images.keys())
        calc = NeighborlistCalculator(cutoff=cutoff)
        self.neighborlist = Data(filename="amp-data-neighborlists", calculator=calc)
        self.neighborlist.calculate_items(self.hashed_images)

    def fit(self, method="Nelder-Mead"):
        log = Logger("results/logs/" + self.filename + ".txt")
        self.target_energies = np.array(
            [(atoms.get_potential_energy()) for atoms in self.data]
        ).reshape(1, -1)
        self.target_forces = np.concatenate(
            [(atoms.get_forces()) for atoms in self.data]
        )
        print("LJ optimization initiated...")
        s_time = time.time()
        bounds = None
        if method != "Nelder-Mead":
            bounds = [(0, None), (None, None), (None, None)] * len(
                self.params_dict.keys()
            )
        lj_min = minimize(
            self.objective_fn,
            self.p0,
            args=(self.target_energies, self.target_forces),
            tol=1e-2,
            method=method,
            bounds=bounds,
        )
        optim_time = time.time() - s_time
        self.logresults(
            log, self.data, self.cutoff, self.p0, self.params_dict, lj_min, optim_time
        )
        if lj_min.success is True:
            print("Optimizer terminated successfully.")
            return lj_min.x
        else:
            print("Optimizer did not terminate successfully.")
            return lj_min.x

    def image_pred(self, image, p0, params_dict):
        chemical_symbols = np.array(image.get_chemical_symbols())
        unique_symbols = np.unique(chemical_symbols)
        possible_pairs = product(unique_symbols, repeat=2)
        possible_pairs = [list(elem) for elem in list(possible_pairs)]
        e_offset = np.sum(
            [params_dict[symbol][2] for symbol in image.get_chemical_symbols()]
        )
        params = []
        for element in chemical_symbols:
            sig = params_dict[element][0]
            eps = params_dict[element][1]
            params.append(np.array([[sig, eps]]))
        params = np.vstack(np.array(params))

        natoms = len(image)

        image_hash = get_hash(image)
        image_neighbors = self.neighborlist[image_hash]

        positions = image.positions
        cell = image.cell

        energy = 0.0
        forces = np.zeros((natoms, 3))

        for a1 in range(natoms):
            sig_1 = np.abs(params[a1][0])
            eps_1 = params[a1][1]
            neighbors, offsets = image_neighbors[a1]
            cells = np.dot(offsets, cell)
            d = positions[neighbors] + cells - positions[a1]
            sig_n = np.abs(params[neighbors][:, 0])
            eps_n = params[neighbors][:, 1]
            sig = (sig_1 + sig_n) / 2
            eps = np.sqrt(eps_1 * eps_n)
            r2 = (d ** 2).sum(1)
            c6 = (sig ** 2 / r2) ** 3
            c6[r2 > self.cutoff ** 2] = 0.0
            c12 = c6 ** 2
            energy += (4 * eps * (c12 - c6)).sum()
            f = (24 * eps * (2 * c12 - c6) / r2)[:, np.newaxis] * d
            forces[a1] -= f.sum(axis=0)
            for a2, f2 in zip(neighbors, f):
                forces[a2] += f2
        energy += e_offset
        return energy, forces, natoms

    def lj_pred(self, data, p0, params_dict):
        params_dict = self.params_to_dict(p0, params_dict)
        predicted_energies = []
        predicted_forces = []
        num_atoms = []
        for image in data:
            energy, forces, natoms = self.image_pred(image, p0, params_dict)
            predicted_energies.append(energy)
            predicted_forces.append(forces)
            num_atoms.append(natoms)
        predicted_energies = np.concatenate(np.array(predicted_energies).reshape(1, -1))
        predicted_forces = np.concatenate(predicted_forces)
        return predicted_energies, predicted_forces, num_atoms

    def objective_fn(self, params, target_energies, target_forces):
        predicted_energies, predicted_forces, num_atoms = self.lj_pred(
            self.data, params, self.params_dict
        )
        data_size = target_energies.shape[1]
        num_atoms_f = np.array([[i] * i for i in num_atoms]).reshape(-1, 1)
        num_atoms = np.array(num_atoms)
        MSE_energy = (1 / data_size) * (
            ((target_energies - predicted_energies) / num_atoms) ** 2
        ).sum()
        MSE_forces = (1 / data_size) * (
            ((target_forces - predicted_forces) / np.sqrt(3 * num_atoms_f)) ** 2
        ).sum()
        if self.fitforces:
            MSE = MSE_energy + MSE_forces
        else:
            MSE = MSE_energy
        return MSE

    def params_to_dict(self, params, params_dict):
        idx = 0
        for keys in list(params_dict.keys()):
            params_dict[keys] = params[idx : idx + 3]
            idx += 3
        return params_dict

    def logresults(self, log, data, cutoff, p0, params_dict, results, optim_time):
        log("%s" % time.asctime())
        log("-" * 50)
        log("LJ-Parameter Optimization")
        log(
            "inital LJ parameter guess [sig, eps]: %s"
            % self.params_to_dict(p0, params_dict)
        )
        log(
            "Optimizer results: \n fun: %s \n message: %s \n nfev: %s \n nit: %s \n success: %s"
            % (
                results["fun"],
                results["message"],
                results["nfev"],
                results["nit"],
                results["success"],
            )
        )
        log("Fitted LJ parameters: %s \n" % self.params_to_dict(results.x, params_dict))
        log("Optimization time: %s \n" % optim_time)

    def parity(self, predicted_energies, predicted_forces):
        fig = plt.figure(figsize=(7.0, 7.0))
        fig2 = plt.figure(figsize=(7.0, 7.0))
        ax = fig.add_subplot(111)
        ax2 = fig2.add_subplot(111)
        predicted_energies = np.squeeze(predicted_energies)
        predicted_forces = np.squeeze(predicted_forces).reshape(1, -1)
        target_energies = np.squeeze(self.target_energies)
        target_forces = np.squeeze(self.target_forces).reshape(1, -1)
        energy_min = min(target_energies)
        energy_max = max(target_energies)
        force_min = min(target_forces)
        force_max = max(target_forces)
        ax.plot(target_energies, predicted_energies, "bo", markersize=3)
        ax.plot([energy_min, energy_max], [energy_min, energy_max], "r-", lw=0.5)
        ax.set_xlabel("ab initio energy, eV")
        ax.set_ylabel("LJ energy, eV")
        ax.set_title("Energy")
        fig.savefig("results/lj_parity_e.pdf")
        ax2.plot(target_forces, predicted_forces, "bo", markersize=3)
        ax2.plot([force_min, force_max], [force_min, force_max], "r-", lw=0.5)
        ax2.set_xlabel("ab initio force, eV/A")
        ax2.set_ylabel("LJ force, eV/A")
        ax2.set_title("Force")
        fig2.savefig("results/lj_parity_f.pdf")
        plt.show()
