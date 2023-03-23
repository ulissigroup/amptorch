# Examples

This directory contains example scripts for AmpTorch module.

## 1. GMP

Scripts for running tarining with Gausian Multipole (GMP) fingerprints and SingleNN atomistic neural networks.

- [`1_GMP_S2E.py`](./1_GMP/1_GMP_S2E.py): script for training with GMP+SingleNN on structure to energy (S2E) task.
- [`2_GMP_S2EF.py`](./1_GMP/2_GMP_S2EF.py): script for training with GMP+SingleNN on structure to energy and forces (S2EF) task.
- [`3_GMP_S2E_w_uncertainty.py`](./1_GMP/3_GMP_S2E_w_uncertainty.py): script for training with GMP+SingleNN on structure to energy (S2E) task with uncertainty quantification.

## 2. SF

Scripts for running Symmetry function (SF) and 2nd generation atomistic neural network (aka BPNN) simulations.

- [`1_SF_S2E.py`](./2_SF/1_SF_S2E.py): script for training with SF+BPNN on structure to energy (S2E) task.
- [`2_SF_S2EF.py`](./2_SF/2_SF_S2EF.py): script for training with SF+BPNN on structure to energy and forces (S2EF) task.

## 3. lmdb

Scripts for constructing and training LMDB databases in case of larger datasets for better memory utilization.

- [`1_construct_lmdb.py`](./3_lmdb/1_construct_lmdb.py): script for constructing an LMDB database from input data.
- [`2_train_lmdb_example.py`](./3_lmdb/2_train_lmdb_example.py): script for training a model using an LMDB database.
- [`3_train_lmdb_full_cache_example.py`](./3_lmdb/3_train_lmdb_full_cache_example.py): script for training a model using an LMDB database with a fully cached environment.
- [`4_train_lmdb_partial_cache_example.py`](./3_lmdb/4_train_lmdb_partial_cache_example.py): script for training a model using an LMDB database with a partially cached environment.

## 4. misc

Miscellaneous scripts.

- [`custom_descriptor_example.py`](./4_misc/custom_descriptor_example.py): script for generating custom molecular descriptors for SF+BPNN case.
- [`get_fp_example.py`](./4_misc/get_fp_example.py): script for generating stand-alone fingerprints.
