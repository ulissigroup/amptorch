[![Build Status](https://travis-ci.org/ulissigroup/amptorch.svg?branch=master)](https://travis-ci.org/ulissigroup/amptorch)
[![Coverage Status](https://coveralls.io/repos/github/ulissigroup/amptorch/badge.svg)](https://coveralls.io/github/ulissigroup/amptorch)
## *AMPTorch*: Atomistic Machine-learning Package - PyTorch

PyTorch is an open-source machine learning (ML) module for Python that provides flexibility and ease in developing, training, and optimizing various ML models (https://pytorch.org). *AMPTorch* is a PyTorch implementation of the *Amp* code that seeks to provide enhanced performance, speed, and flexibility as compared to the original code. The implementation will do so by benefiting from state-of-the-art machine learning methods and capabilities to be optimized in conjunction with high-throughput supercomputers.

This project is being developed at Carnegie Mellon University in the Department of Chemical Engineering, by Muhammed Shuaibi and Zachary Ulissi, in collaboration with Brown University's Andrew Peterson as part of the Department of Energy's *Bridging the time scale in exascale computing of chemical systems* project. 

### Use

AMPTorch has been made to be similar to *Amp* as to require minimal learning for its currents users. An `example.py` file is provided to help demonstrate the steps required to develop and train a model. A calculator instance `AMP(model, label)` is first constructed, with the model containing several arguments that allow for the customizability of the model/training scheme. For a more detailed look into these options, the `model` class can be found in `core.py`. Once the calculator has been fully defined, training and predicting is performed identical to *AMP* with a `.train()` and `.get_potential_energy()`/`get_forces()` pass.

The implementation currently relies on *SIMPLE-NN* to carry out fingerprint/fingerprint derivative calculations. Such calculations will be directly implemented in the near future as to enhance the model's computational performance by taking advantage of PyTorch's ```autograd``` feature. Skorch has been integrated to allow for even greater customizability and flexibility by interfacing with sklearn. An `skorch_example.py` file is provided to demonstrate the steps required to develop and train a model using skorch.

Installation Instructions:
- ```git clone https://github.com/ulissigroup/amptorch.git```
- ```pip install git+https://github.com/mshuaibii/SIMPLE-NN.git```
- ```pip install git+https://bitbucket.org/andrewpeterson/amp.git```
- ```pip install .``` 

As AMPTorch is a work in progress, please feel free to report any issues/suggestions/feedback/questions and we'll try to address them or get back to you as soon as possible.

*Amp* documentation can be found at: https://amp.readthedocs.io/.

### Acknowledgements 
- Funded by the Department of Energy's Basic Enenergy Science, Computational Chemical Sciences Program Office. Award # DE-SC0019441
- *Amp*: Atomistic Machine-Learning Package, https://bitbucket.org/andrewpeterson/amp
- Khorshidi & Peterson, “Amp: A modular approach to machine learning in atomistic simulations”, Computer Physics Communications 207:310-324, 2016. DOI:10.1016/j.cpc.2016.05.010
- PyTorch, https://pytorch.org/
- Skorch, https://skorch.readthedocs.io/en/stable/
