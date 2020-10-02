import numpy as np

from ase.calculators.calculator import Calculator
from amptorch.utils import make_amp_descriptors_simple_nn


__author__ = "Muhammed Shuaibi"
__email__ = "mshuaibi@andrew.cmu.edu"


class EnsembleCalc(Calculator):
    """Atomistics Machine-Learning Potential (AMP) ASE calculator
   Parameters
   ----------
    model : object
        Class representing the regression model. Input arguments include training
        images, descriptor type, and force_coefficient. Model structure and training schemes can be
        modified directly within the class.

    label : str
        Location to save the trained model.

    """

    implemented_properties = ["energy", "forces", "uncertainty"]

    def __init__(self, trained_calcs, training_params):
        Calculator.__init__(self)
        self.trained_calcs = trained_calcs
        self.training_params = training_params
        
    median_list = [100]
    def calculate_stats(self, energies, forces):
        median_idx = np.argsort(energies)[len(energies) // 2]
        energy_median = energies[median_idx]
        energy_var = np.var(energies)
        forces_median = forces[median_idx]
        max_forces_var = np.max(np.var(forces, axis=0))
        return energy_median, forces_median, max_forces_var

    def fingerprint_args(self, images):
        elements = np.array([atom.symbol for atoms in images for atom in atoms])
        _, idx = np.unique(elements, return_index=True)
        elements = list(elements[np.sort(idx)])
        Gs = self.training_params["Gs"]
        return elements, Gs

    def make_fps(self, atoms):
        if isinstance(atoms, list):
            pass
        else:
            atoms = [atoms]
        elements, Gs = self.fingerprint_args(atoms)
        make_amp_descriptors_simple_nn(atoms, Gs, elements, cores=1, label="oal")

    def calculate(self, atoms, properties, system_changes):
        Calculator.calculate(self, atoms, properties, system_changes)
        energies = []
        forces = []

        self.make_fps(atoms)
        for calc in self.trained_calcs:
            energies.append(calc.get_potential_energy(atoms))
            forces.append(calc.get_forces(atoms))
        energies = np.array(energies)
        forces = np.array(forces)
        energy_pred, force_pred, uncertainty = self.calculate_stats(energies, forces)
        
        self.results["energy"] = energy_pred
        self.results["forces"] = force_pred
        atoms.info["uncertainty"] = np.array([uncertainty])
        
