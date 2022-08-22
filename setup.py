#!/usr/bin/env python
from setuptools import find_packages, setup

setup_requires = [
    "cffi>=1.0.0",
]

setup(
    name="amptorch",
    version="0.1",
    description="Atomistic Machine-learning Package - PyTorch",
    author="Muhammed Shuaibi, Yuge Hu, Xiangyun Lei",
    author_email="mshuaibi@andrew.cmu.edu, yhu459@gatech.edu, xlei38@gatech.edu",
    url="https://github.com/ulissigroup/amptorch",
    packages=find_packages(exclude=["contrib", "docs", "tests"]),
    package_data={"": ["*.cpp", "*.h"]},
    python_requires=">=3.6, <4",
    setup_requires=setup_requires,
    cffi_modules=[
        "amptorch/descriptor/Gaussian/libsymf_builder.py:ffibuilder",
        "amptorch/descriptor/GMP/libgmp_builder.py:ffibuilder",
        "amptorch/descriptor/GMPOrderNorm/libgmpordernorm_builder.py:ffibuilder",
    ],
)
