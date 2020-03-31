"""Test script to be executed before pushing or submitting a PR to master
repository."""

import os
import unittest
from consistency_test import test_calcs
from simple_nn_fp_test import test_fp_match
from delta_test import test_skorch_delta
from training_test import test_training, test_e_only_training
from skorch_test import test_skorch, test_e_only_skorch
from fps_from_memory_test import test_fps_memory
from load_test import test_load
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

    def test_delta(self):
        test_skorch_delta()
        print("Delta test passed!")

    def test_skorch(self):
        test_skorch()
        test_e_only_skorch()
        print("Skorch training tests passed!")

    def test_load_fps_from_memory(self):
        test_fps_memory()
        print("Loading fps from memory passed!")

    def test_model_load(self):
        test_load()
        print("Loading trained model test passed!")

    def test_skorch_val(self):
        test_skorch_val()
        test_energy_only_skorch_val()
        print("Skorch validation tests passed!")


if __name__ == "__main__":
    unittest.main(warnings="ignore")
