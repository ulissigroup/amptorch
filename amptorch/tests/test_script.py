"""
Test script to be executed before pushing or submitting a PR to master
repository.
"""

import unittest

from .consistency_test import test_energy_force_consistency
from .cutoff_funcs_test import test_cutoff_funcs
from .gaussian_descriptor_set_test import test_gaussian_descriptor_set
from .pretrained_test import test_pretrained, test_pretrained_no_config
from .pretrained_test_lmdb import test_lmdb_pretrained, test_lmdb_pretrained_no_config
from .training_test import test_training
from .training_test_gmp import test_training_gmp
from .cp_uncertainty_calibration_test import test_cp_uncertainty_calibration


class TestMethods(unittest.TestCase):
    def test_energy_force_consistency(self):
        test_energy_force_consistency()

    def test_cosine_and_polynomial_cutoff_funcs(self):
        test_cutoff_funcs()

    def test_gds(self):
        test_gaussian_descriptor_set()

    def test_load_retrain(self):
        test_pretrained()
        test_pretrained_no_config()

    def test_load_retrain_lmdb(self):
        test_lmdb_pretrained()
        test_lmdb_pretrained_no_config()

    def test_training_scenarios(self):
        test_training()

    def test_training_scenarios_gmp(self):
        test_training_gmp()

    def test_uncertainty_cp(self):
        test_cp_uncertainty_calibration()


if __name__ == "__main__":
    unittest.main(warnings="ignore")
