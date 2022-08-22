.. _install:

==================================
Installation
==================================

Install dependencies:

1. Ensure conda is up-to-date: ``conda update conda``

2. Create conda environment

-  CPU machines: ``conda env create -f env_cpu.yml``
-  GPU machines (CUDA 10.2): ``conda env create -f env_gpu.yml``

3. Activate the conda environment ``conda activate amptorch`` and
   install the package with ``pip install -e .``

4. Install pre-commit hooks: ``pre-commit install``
