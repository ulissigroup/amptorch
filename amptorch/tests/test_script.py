"""Test script to be executed before pushing or submitting a PR to master
repository."""

import os
import unittest

from bootstrap_test import test_bootstrap
from consistency_test import test_calcs
from delta_test import test_skorch_delta
from fps_from_memory_test import test_fps_memory
from load_test import test_load
from simple_nn_fp_test import test_fp_match
from skorch_test import test_e_only_skorch, test_skorch
from training_test import test_e_only_training, test_training
from val_test import (test_energy_only_skorch_val, test_energy_only_val,
                      test_skorch_val, test_val)


class TestMethods(unittest.TestCase):
    def test_consistency(self):
        # test_calcs()
        print("Energy and Force consistency test passed!")

    # def test_fps(self):
    # test_fp_match()
    # print("Fingerprint test passed!")

    # def test_delta(self):
    # test_skorch_delta()
    # print("Delta test passed!")

    # def test_skorch(self):
    # test_skorch()
    # test_e_only_skorch()
    # print("Skorch training tests passed!")

    # def test_load_fps_from_memory(self):
    # test_fps_memory()
    # print("Loading fps from memory passed!")

    # def test_model_load(self):
    # test_load()
    # print("Loading trained model test passed!")

    # def test_skorch_val(self):
    # test_skorch_val()
    # test_energy_only_skorch_val()
    # print("Skorch validation tests passed!")

    # def test_active_learning(self):
    # test_bootstrap()
    # TODO write remaining


if __name__ == "__main__":
    unittest.main(warnings="ignore")
