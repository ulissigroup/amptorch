"""Test script to be executed before pushing or submitting a PR to master
repository."""

import os
import unittest
from consistency_test import test_calcs
from simple_nn_fp_test import test_fp_match
from lj_test import test_lj
from training_test import test_training


class TestMethods(unittest.TestCase):
    def test_consistency(self):
        test_calcs()
        print("Energy and Force test passed!")

    def test_fps(self):
        test_fp_match()
        print("Fingerprint test passed!")

    def test_lj(self):
        test_lj()
        print("LJ test passed!")

    def test_training(self):
        test_training()
        print("Training test passed!")


if __name__ == "__main__":
    unittest.main(warnings="ignore")
