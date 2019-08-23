import sys
import os
import time
import numpy as np
from scipy.optimize import minimize, shgo
import ase
from ase import Atoms
from ase.build import molecule
from ase.calculators.emt import EMT
from ase.calculators.lj import LennardJones
from ase.neighborlist import NeighborList, NewPrimitiveNeighborList
from amp_pytorch.neighbors import get_distances
from amp.utilities import Logger
import matplotlib.pyplot as plt


class lj_optim:
    def __init__(self, data, params, params_dict, cutoff):
        if not os.path.exists("results"):
            os.mkdir("results")
        self.data = data
        self.p0 = params
        self.params_dict = params_dict
        # self.lj_param_check()
        self.cutoff = cutoff

    def fit(self, filename="results/lj_log.txt", method="Nelder-Mead"):
        log = Logger(filename)
        self.target_energies = np.array(
            [(atoms.get_potential_energy()) for atoms in self.data]
        ).reshape(1, -1)
        self.target_forces = np.concatenate(
            [(atoms.get_forces()) for atoms in self.data]
        )
        print("LJ optimization initiated...")
        s_time = time.time()
        bounds = ((0, None), (None, None), (None, None), (0, None), (None,
            None), (None, None), (0, None), (None, None), (None, None))
        # bounds = ((0, None), (None, None), (None, None))
        lj_min = minimize(
            self.objective_fn,
            self.p0,
            args=(self.target_energies, self.target_forces),
            # tol=1e-2,
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

    def lj_pred(self, data, p0, params_dict):
        # p0 = np.concatenate([p0]*len(params_dict.keys()))
        params_dict = self.params_to_dict(p0, params_dict)
        predicted_energies = []
        predicted_forces = []
        num_atoms = []
        for image in data:
            chemical_symbols = image.get_chemical_symbols()
            params = []
            for element in chemical_symbols:
                sig = params_dict[element][0]
                eps = params_dict[element][1]
                offset = params_dict[element][2]
                params.append(np.array([[sig, eps, offset]]))
            params = np.vstack(np.array(params))

            natoms = len(image)
            num_atoms.append(natoms)

            self.n1 = NeighborList([self.cutoff / 2] * natoms, self_interaction=False)
            self.n1.update(image)

            positions = image.positions
            cell = image.cell

            energy = 0.0
            forces = np.zeros((natoms, 3))

            for a1 in range(natoms):
                sig_1 = params[a1][0]
                eps_1 = params[a1][1]
                o1 = params[a1][2]
                neighbors, offsets = self.n1.get_neighbors(a1)
                cells = np.dot(offsets, cell)
                d = positions[neighbors] + cells - positions[a1]

                sig_n = params[neighbors][:, 0]
                eps_n = params[neighbors][:, 1]
                sig = (sig_1 + sig_n) / 2
                eps = np.sqrt(eps_1 * eps_n)
                on = params[neighbors][:, 2]
                r2 = (d ** 2).sum(1)
                c6 = (sig ** 2 / r2) ** 3
                c6[r2 > self.cutoff ** 2] = 0.0
                c12 = c6 ** 2
                energy += (4 * eps * (c12 - c6)).sum() + o1 * len(on) + on.sum()
                # energy += (4 * eps * (c12 - c6)).sum()
                f = (24 * eps * (2 * c12 - c6) / r2)[:, np.newaxis] * d
                forces[a1] -= f.sum(axis=0)
                for a2, f2 in zip(neighbors, f):
                    forces[a2] += f2
            predicted_energies.append(energy)
            predicted_forces.append(forces)
        predicted_energies = np.array(predicted_energies).reshape(1, -1)
        predicted_forces = np.concatenate(predicted_forces)
        return predicted_energies, predicted_forces, num_atoms

    def objective_fn(self, params, target_energies, target_forces):
        predicted_energies, predicted_forces, num_atoms = self.lj_pred(
            self.data, params, self.params_dict
        )
        data_size = target_energies.shape[1]
        num_atoms_f = np.array([[i] * i for i in num_atoms]).reshape(-1, 1)
        num_atoms = np.array(num_atoms)
        MSE_energy = np.mean((target_energies - predicted_energies) ** 2)
        MSE_forces = (1 / data_size) * ((target_forces - predicted_forces) ** 2).sum()
        # MSE_energy = np.mean(((target_energies - predicted_energies) / num_atoms) ** 2)
        # MSE_energy = (1 / data_size) * (
        # ((target_energies - predicted_energies) / num_atoms) ** 2
        # ).sum()
        # MSE_forces = (1 / data_size) * (
        # ((target_forces - predicted_forces) / np.sqrt(3 * num_atoms_f)) ** 2
        # ).sum()
        MSE = MSE_energy + 0.3 * MSE_forces
        return MSE

    def lj_param_check(self):
        unique = set()
        for atoms in self.data:
            symbols = atoms.symbols
            unique = unique | set(symbols)
        unique_elements = list(unique)
        num_lj_params = 3 * len(unique_elements)
        assert (
            len(self.p0) == num_lj_params
        ), "Number of initial conditions not equal to \
        the number of required LJ parameters"

    def params_to_dict(self, params, params_dict):
        idx = 0
        for keys in list(params_dict.keys()):
            params_dict[keys] = params[idx : idx + 3]
            idx += 3
        return params_dict

    def logresults(self, log, data, cutoff, p0, params_dict, results, optim_time):
        log("-" * 50)
        log("LJ-Parameter Optimization")
        log("Fits provided data to the Lennard Jones model")
        log("%s" % time.asctime())
        log("Dataset size: %s" % len(data))
        log("cutoff radius: %s" % (cutoff))
        log(
            "inital LJ parameter guess [sig, eps]: %s"
            % self.params_to_dict(p0, params_dict)
        )
        log("Optimizer results: \n %s \n" % results)
        log("Fitted LJ parameters: %s \n" % self.params_to_dict(results.x, params_dict))
        log("Optimization time: %s" % optim_time)

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
