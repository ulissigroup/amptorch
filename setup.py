#!/usr/bin/env python
from setuptools import setup, find_packages, Extension
# To use a consistent encoding
from codecs import open
from os import path, listdir
from pkg_resources import DistributionNotFound, get_distribution
from subprocess import check_output

install_requires = [
    "torch",
    "spglib",
    # "scikit-learn==0.21.3",
    "skorch==0.6.0",
    "scipy",
    "pandas",
    'numpy',
    'ase>=3.10.0',
    'cffi>=1.0.0',
    'h5py',
]

setup_requires = [
    'cffi>=1.0.0',
]


setup(
    name="amptorch",
    version="0.1",
    description="A Pytorch Implementation of AMP",
    author="Muhammed Shuaibi, Xiangyun Lei",
    author_email="mshuaibi@andrew.cmu.edu, xlei38@gatech.edu",
    url="https://github.com/ulissigroup/amptorch",
    classifiers=[   # Optional
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # other arguments are listed here.

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',

    ],
    packages=find_packages(exclude=['contrib', 'docs', 'tests']),
    package_data={'':['*.cpp', '*.h']},
    python_requires='>=2.7, <4',
    install_requires=install_requires,
    setup_requires=setup_requires,
    cffi_modules=[
        "amptorch/descriptor/Gaussian/libsymf_builder.py:ffibuilder",
        "amptorch/descriptor/MCSH/libmcsh_builder.py:ffibuilder",
    ],
)
