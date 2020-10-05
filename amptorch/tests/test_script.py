"""Test script to be executed before pushing or submitting a PR to master
repository."""

import unittest

from debug_travis import test_travis


class TestMethods(unittest.TestCase):
    def debug_travis(self):
        test_travis()


if __name__ == "__main__":
    unittest.main(warnings="ignore")
