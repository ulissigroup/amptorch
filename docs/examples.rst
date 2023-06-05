.. _examples:

==================================
Examples
==================================

--------------------------------------------
Gaussian Multiple (GMP) Descriptors
--------------------------------------------

In addition to conventional Atom-centered Symmetry Functions as
fingerprinting scheme, AmpTorch also support GMP descriptors that uses
multipole expansions to describe the reconstructed electronic density
around every central atom and its neighbors to encode local
environments. Because the formulation of symmetry functions does not
take into element types into account, the interactions among different
elements are divided into different columns as input. As a result, the
number of feature dimensions undesirably increases with the number of
elements present. A major advantage of GMPs is that the input dimensions
remain constant regardless of the number of chemical elements, and
therefore can be adopted for complex datasets. For more technical
details and theorical backgrounds, please refer to *Lei, X., & Medford,
A. J. (2021). A Universal Framework for Featurization of Atomistic
Systems. http://arxiv.org/abs/2102.02390*

For example scripts of using GMP for structure to energy (and forces),
please refer to:

::

   examples/1_GMP

------------------------------------------------------
SingleNN Atomistic Neural Network Structures
------------------------------------------------------


As GMPs encode the information about chemical elements based on
reconstructed electronic environments, GMPs work naturally with the
atomistic Neural Network Structures SingleNN as published by Liu and
Kitchin (*Liu, M., & Kitchin, J. R. (2020). SingleNN: Modified
Behler-Parrinello Neural Network with Shared Weights for Atomistic
Simulations with Transferability. Journal of Physical Chemistry C,
124(32), 17811–17818. https://doi.org/10.1021/acs.jpcc.0c04225*).

To use SingleNN instead of the default Behler-Parrinello
High-dimensional Neural Network scheme, in ``config`` for NN trainer,
define:

::

   config["model"]["name"] == "singlenn"

as shown in ``examples/1_GMP/1_GMP_S2E.py``


----------------------------------------------------------------
lmdb as Database Management Solution for Large Dataset
----------------------------------------------------------------


For AmpTorch to be compatible to train with large datasets such as `Open
Catalyst
Project <https://github.com/Open-Catalyst-Project/baselines>`__, we
leverage ``lmdb``, a Btree-based database management library, to resolve
possible memory issues when it comes to loading and training. It can be
used in either full- or partial-cache fashion depending on whether the
dataset can be fit into RAM altogether.

Examples of using no, full or partial cache can be found in:

::

   examples/3_lmdb

--------------------------------------------------------------------------
Uncertainty Quantification (UQ) via Conformal Prediction (CP)
--------------------------------------------------------------------------


AmpTorch implements UQ as an optional feature during the prediction.
Here we use `conformal prediction
method <https://arxiv.org/abs/2208.08337>`__ with the distances in
neural network’s latent space to output the uncertainty associated with
the predicted energy. CP method ensures *calibration* while showing
advantage of being *sharp* and *scalable* when tested against
benchmarking systems such as MD17, QM9 and OC20 with trained models.

An example python script can be found in:

::

   examples/1_GMP/3_GMP_S2E_w_uncertainty.py
