"""
Test script to be executed before pushing or submitting a PR to master
repository.
"""

import unittest

from .consistency_test import test_energy_force_consistency
from .pretrained_test import test_pretrained
from .training_test import test_training


class TestMethods(unittest.TestCase):
    def test_training_scenarios(self):
        test_training()

    def test_load_retrain(self):
        test_pretrained()

    def test_consistency(self):
        test_energy_force_consistency()


if __name__ == "__main__":
    unittest.main(warnings="ignore")
