"""Test script to be executed before pushing or submitting a PR to master
repository."""

import os
import unittest
from consistency_test import test_calcs
from simple_nn_fp_test import test_fp_match


class TestMethods(unittest.TestCase):
    def test_consistency(self):
        test_calcs()
        print("Energy and Force test passed!")

    def test_fps(self):
        test_fp_match()
        print("Fingerprint test passed!")


if __name__ == "__main__":
    unittest.main(warnings="ignore")
