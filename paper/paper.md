---
title: 'AmpTorch: A Python package for scalable fingerprint-based neural network training on multi-element systems with integrated uncertainty quantification'
tags:
  - Python
  - machine learning interatomic potentials
  - neural networks
  - molecular dynamics
authors:
  - name: Muhammed Shuaibi
    equal-contrib: true
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Yuge Hu
    orcid: 0000-0003-3648-7749
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 2
  - name: Xiangyun Lei
    equal-contrib: false
    affiliation: 2
  - name: Benjamin M. Comer
    equal-contrib: false
    affiliation: 2
  - name: Matt Adams
    equal-contrib: false
    affiliation: 1
  - name: Jacob Paras
    equal_contrib: false
    affiliation: 3
  - name: Rui Qi Chen
    equal-contrib: false
    affiliation: 2
  - name: Eric Musa
    equal-contrib: false
    affiliation: 4
  - name: Joseph Musielewicz
    equal-contrib: false
    affiliation: 1
  - name: Andrew A. Peterson
    orcid: 0000-0003-2855-9482
    corresponding: false
    affiliation: 5
  - name:  Andrew J. Medford
    orcid: 0000-0001-8311-9581
    corresponding: false
    affiliation: 2
  - name: Zachary Ulissi
    orcid: 0000-0002-9401-4918
    corresponding: true # (This is how to denote the corresponding author)
    affiliation: 1
affiliations:
 - name: Department of Chemical Engineering, Carnegie Mellon University, United States
   index: 1
 - name: Department of Chemical and Biomolecular Engineering, Georgia Institute of Technology, United States
   index: 2
 - name: School of Physics and School of Computer Science, Georgia Institute of Technology, United States
   index: 3
 - name: Department of Chemical Engineering, University of Michigan, United States
   index: 4
 - name: School of Engineering, Brown University, United States
   index: 5

date: 7 June 2023
bibliography: paper.bib

doi:
---

# Summary

Machine learning (ML) force fields use ML models to predict the energy (and forces) given a chemical system defined by a set of atomic coordinates. Fingerprint-based and graph-based ML force fields are two major categories of atomistic machine-learning models. Fingerprint-based models convert the atomic coordinates into fixed-length "feature vectors" that are used as inputs to regression models such as neural networks or kernel ridge regression, while graph-based systems convert the coordinates into a molecular graph which is directly input into a deep learning model. Fingerprint-based approaches tend to require less training data and be more interpretable, while graph-based models are typically able to handle more complex chemical systems and are more accurate in the limit of large datasets. Fingerprint-based models are widely used due to their conceptual similarity to traditional force-fields, but most existing packages for fingerprint-based neural network force fields are limited to systems with relatively few ($<5$) elements and relatively small ($<10^6$) datasets.

This work introduces `AmpTorch`, a Python/C++ package that leverages the Gaussian Multipole (GMP) fingerprinting scheme to build atomistic neural network models [@Lei:2022]. It provides an efficient training routine scalable to $\sim10^6$ training points, is compatible with many-element systems ($50+$ elements), and supports statistically rigorous uncertainty quantification (UQ) during the prediction step [@Hu:2022]. The fingerprint-based structure provides relatively simple models compared to graph-based alternatives, and the structure of the fingerprinting scheme allows the package to handle systems with an arbitrary number of elements and arbitrary boundary conditions (isolated, periodic, semi-periodic).


`AmpTorch` is a PyTorch adaptation of the Atomistic Machine-learning Package `AMP` [@Khorshidi:2016]. It maintains the modular structure of `AMP`, with separate modules for generating fingerprints and training neural network models. `AmpTorch` supports standard symmetry function fingerprints and high-dimensional neural network potentials, but is able to scale to larger datasets and more complex chemical systems than `AMP` or other existing codes for feature-based ML potentials (for example, Atom-Centered Symmetry Functions [@Behler:2015] whose fingerprinting dimension scales with the number of elements). `AmpTorch` outsources traditional boilerplate trainer code to [`Skorch`](https://skorch.readthedocs.io/en/stable/) [@Tietz:2017]. `Skorch` serves as a PyTorch wrapper that allows users to easily modify traditional neural network training strategies. `AmpTorch`'s use of `Skorch` allows users to trivially modify training strategies, model architectures, and hyperparameters.

The scalability of `AmpTorch` is a result of its fingerprinting scheme, software design, and database management. The GMP fingerprints have a fixed vector length regardless of the number of elements, and using Gaussian functions allows all necessary integrals to be computed analytically. The analytical solutions are coupled with a C++ implementation, ensuring the efficiency of the fingerprint calculation. The GMP fingerprints are also naturally compatible with the SingleNN [@Liu:2020]  neural network structure, allowing for neural networks that have the same architecture and number of fitted parameters regardless of the number of elements present in the dataset.

In addition, `AmpTorch` leverages the B-tree-based database management library, [`LMDB`](http://www.lmdb.tech/doc/) [@Chu:2010], to resolve possible memory issues when it comes to loading and training on large datasets. This allows `AmpTorch` to be trained on datasets that are too large to fit into temporary memory, although reading and writing data from `LMDB` is slower than using temporary memory.

`AmpTorch` also implements UQ as an optional feature during the prediction. Supported UQ methods include the ensemble method, dropout neural network method, and methods based on latent distances. In particular, `AmpTorch` includes an implementation of a highly scalable new approach that leverages distances in the latent space and the "conformal prediction" statistical technique to provide statistically rigorous error bars on complex pre-trained models.

`AmpTorch` is designed using standards from the `ASE` package for handling atomic structures [@Larsen:2017]. It takes a list of `ase.Atoms` objects with energy (and forces) as input and output `ase.Calculator` object that can be used to compute the energy (and forces) with the trained model. This structure allows for easy integration with active learning packages such as `FINETUNA` [@Musielewicz:2022] and data visualization tools such as `ElectroLens` [@Lei:2019].

# Statement of need

There are numerous software packages for both constructing and applying ML force-fields. These include packages that are based on deep learning and graph convolutional models, kernel-based regression models, and neural networks with fixed feature vectors. Of these packages, only graph convolutional codes have been demonstrated for extremely large and complex datasets with $>$ 10 elements and $>$ 1M training points such as the OC20 dataset in [Open Catalyst Project](https://pubs.acs.org/doi/10.1021/acscatal.0c04525) [@Chanussot:2021]. However, at inference time, adapting these often complex architectures with 1-100M parameter models into atomistic simulation interfaces seeking extremely fast calls is a challenging task. On the other hand, all existing packages for training feature-based models are limited in both the number of elements they can handle (due to the poor scaling of feature vector size) and the number of training points (due to scaling of training kernel-based models and memory management issues in most codes). Thus, `AMPTorch` aims to fill this gap by developing a toolkit to train fingerprint-based neural network force fields on datasets of an arbitrary number of chemical elements. For this reason, we expect that the code will be widely used by researchers seeking to train ML force-field models for complex systems that can be integrated with existing atomistic simulation pipelines.

# Acknowledgements

The  authors  are  grateful  for funding from the U.S. Department of Energy’s Basic Energy Science, Computational Chemical Sciences Program Office, under Award No. DE-SC0019441.
