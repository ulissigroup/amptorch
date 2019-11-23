"""Test script to be executed before pushing or submitting a PR to master
repository."""

import os
import unittest
from consistency_test import test_calcs
from simple_nn_fp_test import test_fp_match
from lj_test import test_lj
from training_test import test_training, test_e_only_training
from skorch_test import test_skorch, test_e_only_skorch
from val_test import (
    test_skorch_val,
    test_val,
    test_energy_only_skorch_val,
    test_energy_only_val,
)


class TestMethods(unittest.TestCase):
    def test_consistency(self):
        test_calcs()
        print("Energy and Force consistency test passed!")

    def test_fps(self):
        test_fp_match()
        print("Fingerprint test passed!")

    def test_lj(self):
        test_lj()
        print("LJ test passed!")

    def test_training(self):
        test_training()
        test_e_only_training()
        print("Custom training test passed!")

    def test_skorch(self):
        test_skorch()
        test_e_only_skorch()
        print("Skorch training test passed!")

    def test_skorch_val(self):
        test_skorch_val()
        test_energy_only_skorch_val()
        print("Skorch validation test passed!")

    def test_val(self):
        test_val()
        test_energy_only_val()
        print("Custom training validation test passed!")


if __name__ == "__main__":
    unittest.main(warnings="ignore")
