from ase.calculators.calculator import Calculator


class AmpTorch(Calculator):
    """
    Create an ase.calculators.calculator.Calculator class to compute the energy (and forces) for the given ase.Atoms object.

    Method:
    ------------
    calculate : Calculates the corresponding energy (and forces) with loaded parameters in the model.
    """

    implemented_properties = ["energy", "forces"]

    def __init__(self, trainer):
        Calculator.__init__(self)

        self.trainer = trainer

    def calculate(self, atoms, properties, system_changes):
        Calculator.calculate(self, atoms, properties, system_changes)

        predictions = self.trainer.predict([atoms])

        self.results["energy"] = predictions["energy"][0]
        self.results["forces"] = predictions["forces"][0]
