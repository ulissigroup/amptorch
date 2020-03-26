[![Build Status](https://travis-ci.org/ulissigroup/amptorch.svg?branch=master)](https://travis-ci.org/ulissigroup/amptorch)
[![Coverage Status](https://coveralls.io/repos/github/ulissigroup/amptorch/badge.svg?branch=master)](https://coveralls.io/github/ulissigroup/amptorch?branch=master)
## *AMPTorch*: Atomistic Machine-learning Package - PyTorch

PyTorch is an open-source machine learning (ML) module for Python that provides flexibility and ease in developing, training, and optimizing various ML models (https://pytorch.org). *AMPTorch* is a PyTorch implementation of the *Amp* code that seeks to provide enhanced performance, speed, and flexibility as compared to the original code. The implementation will do so by benefiting from state-of-the-art machine learning methods and capabilities to be optimized in conjunction with high-throughput supercomputers.

This project is being developed at Carnegie Mellon University in the Department of Chemical Engineering, by Muhammed Shuaibi and Zachary Ulissi, in collaboration with Brown University's Andrew Peterson as part of the Department of Energy's *Bridging the time scale in exascale computing of chemical systems* project. 

### Use

AMPTorch has been made to be similar to *Amp* as to require minimal learning for its currents users. An `example.py` file is provided to help demonstrate the steps required to develop and train a model. A calculator instance `AMP(model, label)` is first constructed, with the model containing several arguments that allow for the customizability of the model/training scheme. For a more detailed look into these options, the `model` class can be found in `core.py`. Once the calculator has been fully defined, training and predicting is performed identical to *AMP* with a `.train()` and `.get_potential_energy()`/`get_forces()` pass.

Fingerprints are not handled directly within the current implementation of Amptorch, and one can rely on either the standard *Amp* code or *SIMPLE-NN* for fingerprint / fingerprint-derivative calculations.  Such calculations will be directly implemented in the near future as to enhance the model's computational performance by taking advantage of PyTorch's ```autograd``` feature. Skorch has been integrated to allow for even greater customizability and flexibility by interfacing with sklearn. An `skorch_example.py` file is provided to demonstrate the steps required to develop and train a model using skorch.

Installation Instructions:

The simplest way to install Amptorch is with pip; you can follow the instructions below to use this method. The SIMPLE-NN line is an optional dependency; note that this will also result in a tensorflow installation.

- ```pip install git+https://github.com/ulissigroup/amptorch.git```
- ```pip install git+https://github.com/mshuaibii/SIMPLE-NN.git```
- ```pip install amp-atomistics```

For more customized installations, for example if you want to be able to easily play with the code of Amptorch without reinstalling, you can do the following:

- Install and make sure you can import (into python3) the following dependencies (names for pip installations are in square brackets):
    - ase [ase]
    - amp [amp-atomistics]
    - scipy [scipy]
    - torch [torch]
    - skorch [skorch]
    - spglib [spglib]

- Clone the repository and modify your environment such that python knows it exists. E.g., the second line might be included in your .bashrc file or in a custom load function.

    - ```git clone https://github.com/ulissigroup/amptorch.git```
    - ```export PYTHONPATH=/path/to/amptorch/:$PYTHONPATH```

### Suggestions and further information

As AMPTorch is a work in progress, please feel free to report any issues/suggestions/feedback/questions and we'll try to address them or get back to you as soon as possible.

*Amp* documentation can be found at: https://amp.readthedocs.io/.

### Acknowledgements 
- Funded by the Department of Energy's Basic Enenergy Science, Computational Chemical Sciences Program Office. Award # DE-SC0019441
- *Amp*: Atomistic Machine-Learning Package, https://bitbucket.org/andrewpeterson/amp
- Khorshidi & Peterson, “Amp: A modular approach to machine learning in atomistic simulations”, Computer Physics Communications 207:310-324, 2016. DOI:10.1016/j.cpc.2016.05.010
- PyTorch, https://pytorch.org/
- Skorch, https://skorch.readthedocs.io/en/stable/
