from ase import Atoms
from ase.calculators.lj import LennardJones as LJ
from amptorch.lj_model import lj_optim


def test_lj():
    atoms = Atoms("CuCuCu", [(0, 0, 1), (0, 0, 1.5), (0, 0, 0.5)])
    atoms.set_calculator(LJ())
    actual_energy = atoms.get_potential_energy()
    actual_forces = atoms.get_forces()

    p0 = [1, 1, 0]
    params_dict = {"Cu": []}

    atoms = [atoms]
    lj_model = lj_optim("test", atoms, p0, params_dict, 10)
    fitted_params = p0
    lj_energies, lj_forces, num_atoms = lj_model.lj_pred(
        atoms, fitted_params, params_dict
    )

    assert round(actual_energy, 1) == lj_energies[0], "LJ energies don't match!"
    assert actual_forces.all() == lj_forces.all(), "LJ forces don't match!"
