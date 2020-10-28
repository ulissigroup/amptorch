[![Build Status](https://travis-ci.org/ulissigroup/amptorch.svg?branch=master)](https://travis-ci.org/ulissigroup/amptorch)
[![Coverage Status](https://coveralls.io/repos/github/ulissigroup/amptorch/badge.svg?branch=master)](https://coveralls.io/github/ulissigroup/amptorch?branch=master)
## *AMPtorch*: Atomistic Machine-learning Package - PyTorch

*AMPtorch* is a PyTorch implementation of the [Atomistic Machine-learning Package](https://amp.readthedocs.io/en/latest/) (AMP) code that seeks to provide users with improved performance and flexibility as compared to the original code. The implementation does so by benefiting from state-of-the-art machine learning methods and techniques to be optimized in conjunction with high-throughput supercomputers. *AMPtorch* is built on top of [PyTorch Geometric](https://pytorch-geometric.readthedocs.io/en/latest/) and [Skorch](https://skorch.readthedocs.io/en/stable/).

This project is being developed at Carnegie Mellon University in the Department of Chemical Engineering, by Muhammed Shuaibi and Zachary Ulissi, in collaboration with Brown University's Andrew Peterson as part of the Department of Energy's *Bridging the time scale in exascale computing of chemical systems* project.

### Installation

Install dependencies:

1. Ensure conda is up-to-date: ```conda update conda```

2. Create conda environment
- CPU machines: ```conda env create -f env_cpu.yml```
- GPU machines (CUDA 10.2): ```conda env create -f env_gpu.yml```

3. Activate the conda environment `conda activate amptorch` and install the package with `pip install -e .`

4. Install pre-commit hooks: `pre-commit install`

### Use
#TODO

### Acknowledgements
- Funded by the Department of Energy's Basic Enenergy Science, Computational Chemical Sciences Program Office. Award # DE-SC0019441
- Engineering ideas have been heavily borrowed from our work on the [Open Catalyst Project](https://github.com/Open-Catalyst-Project/baselines)
