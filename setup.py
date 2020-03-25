#!/usr/bin/env python

from distutils.core import setup
from setuptools import find_packages

setup(
    name="amptorch",
    version="0.1",
    description="A Pytorch Implementation of AMP",
    author="Muhammed Shuaibi",
    author_email="mshuaibi@andrew.cmu.edu",
    url="https://github.com/ulissigroup/amptorch",
    packages=find_packages(),
    install_requires=["spglib", "torch", "scikit-learn==0.21.3", "skorch==0.6.0", "ase", "scipy"],
)
