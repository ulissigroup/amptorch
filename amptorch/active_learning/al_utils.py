import copy
from ase.calculators.calculator import Calculator
from ase.calculators.singlepoint import SinglePointCalculator as sp


class CounterCalc(Calculator):
    implemented_properties = ["energy", "forces", "uncertainty"]
    """Parameters
    --------------
        calc: object. Parent calculator to track force calls."""

    def __init__(self, calc, **kwargs):
        super().__init__()
        self.calc = calc
        self.force_calls = 0

    def calculate(self, atoms, properties, system_changes):
        super().calculate(atoms, properties, system_changes)
        calc = copy.deepcopy(self.calc)
        self.results["energy"] = calc.get_potential_energy(atoms)
        self.results["forces"] = calc.get_forces(atoms)
        self.force_calls += 1


def attach_sp_calc(images):
    "Converts images calculators to single point calculators to avoid double calculations'"
    if isinstance(images, list):
        for image in images:
            construct_sp(image)
    else:
        construct_sp(images)
    return images


def construct_sp(image):
    sample_energy = image.get_potential_energy(apply_constraint=False)
    sample_forces = image.get_forces(apply_constraint=False)
    image.set_calculator(
        sp(atoms=image, energy=sample_energy, forces=sample_forces)
    )
    return image
