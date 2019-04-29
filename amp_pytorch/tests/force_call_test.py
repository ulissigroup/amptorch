"""
Exact Gaussian-neural scheme forces and energies of five different non-periodic
configurations and three different periodic configurations have been calculated
in Mathematica, and are given below.  This script checks the values calculated
by the code with and without fortran modules.

"""

import sys
import numpy as np
from ase import Atoms
from ase.calculators.emt import EMT
from collections import OrderedDict
from amp import Amp
from amp.descriptor.gaussian import Gaussian
from amp.model.neuralnetwork import NeuralNetwork
from amp_pytorch import core_test
from amp_pytorch.data_preprocess import AtomsDataset
from amp_pytorch import data_preprocess
import torch
import torch.nn as nn
import torch.optim as optim



def non_periodic_test():
    """Gaussian/Neural non-periodic standard.

    Checks that the answer matches that expected from previous Mathematica
    calculations.
    """

    # Making the list of non-periodic images
    images = [Atoms(symbols='PdOPd2',
                    pbc=np.array([False, False, False], dtype=bool), calculator=EMT(),
                    cell=np.array(
                        [[1.,  0.,  0.],
                         [0.,  1.,  0.],
                            [0.,  0.,  1.]]),
                    positions=np.array(
                        [[0.,  0.,  0.],
                         [0.,  2.,  0.],
                            [0.,  0.,  3.],
                            [1.,  0.,  0.]])),
              Atoms(symbols='PdOPd2',
                    pbc=np.array([False, False, False],
                        dtype=bool), calculator=EMT(),
                    cell=np.array(
                        [[1.,  0.,  0.],
                         [0.,  1.,  0.],
                            [0.,  0.,  1.]]),
                    positions=np.array(
                        [[0.,  1.,  0.],
                         [1.,  2.,  1.],
                            [-1.,  1.,  2.],
                            [1.,  3.,  2.]])),
              Atoms(symbols='PdO',
                    pbc=np.array([False, False, False], dtype=bool), calculator=EMT(),
                    cell=np.array(
                        [[1.,  0.,  0.],
                         [0.,  1.,  0.],
                         [0.,  0.,  1.]]),
                    positions=np.array(
                        [[2.,  1., -1.],
                         [1.,  2.,  1.]])),
              Atoms(symbols='Pd2O',
                    pbc=np.array([False, False, False], dtype=bool),
                    calculator=EMT(),
                    cell=np.array(
                        [[1.,  0.,  0.],
                         [0.,  1.,  0.],
                         [0.,  0.,  1.]]),
                    positions=np.array(
                        [[-2., -1., -1.],
                         [1.,  2.,  1.],
                         [3.,  4.,  4.]])),
              Atoms(symbols='Cu',
                    pbc=np.array([False, False, False], dtype=bool),
                    calculator=EMT(),
                    cell=np.array(
                        [[1.,  0.,  0.],
                         [0.,  1.,  0.],
                         [0.,  0.,  1.]]),
                    positions=np.array(
                        [[0.,  0.,  0.]]))]
    # Correct energies and forces
    for image in images:
        image.get_potential_energy(apply_constraint=False)
        image.get_forces(apply_constraint=False)

    # Parameters
    Gs = {'O': [{'type': 'G2', 'element': 'Pd', 'eta': 0.8},
                {'type': 'G4', 'elements': [
                    'Pd', 'Pd'], 'eta':0.2, 'gamma':0.3, 'zeta':1},
                {'type': 'G4', 'elements': ['O', 'Pd'], 'eta':0.3, 'gamma':0.6,
                 'zeta':0.5}],
          'Pd': [{'type': 'G2', 'element': 'Pd', 'eta': 0.2},
                 {'type': 'G4', 'elements': ['Pd', 'Pd'],
                  'eta':0.9, 'gamma':0.75, 'zeta':1.5},
                 {'type': 'G4', 'elements': ['O', 'Pd'], 'eta':0.4,
                  'gamma':0.3, 'zeta':4}],
          'Cu': [{'type': 'G2', 'element': 'Cu', 'eta': 0.8},
                 {'type': 'G4', 'elements': ['Cu', 'O'],
                  'eta':0.2, 'gamma':0.3, 'zeta':1},
                 {'type': 'G4', 'elements': ['Cu', 'Cu'], 'eta':0.3,
                  'gamma':0.6, 'zeta':0.5}]}

    hiddenlayers = {'O': (2,), 'Pd': (2,), 'Cu': (2,)}

    weights = OrderedDict([('O', OrderedDict([(1, np.matrix([[-2.0, 6.0],
                                                             [3.0, -3.0],
                                                             [1.5, -0.9],
                                                             [-2.5, -1.5]])),
                                              (2, np.matrix([[5.5],
                                                             [3.6],
                                                             [1.4]]))])),
                           ('Pd', OrderedDict([(1, np.matrix([[-1.0, 3.0],
                                                              [2.0, 4.2],
                                                              [1.0, -0.7],
                                                              [-3.0, 2.0]])),
                                               (2, np.matrix([[4.0],
                                                              [0.5],
                                                              [3.0]]))])),
                           ('Cu', OrderedDict([(1, np.matrix([[0.0, 1.0],
                                                              [-1.0, -2.0],
                                                              [2.5, -1.9],
                                                              [-3.5, 0.5]])),
                                               (2, np.matrix([[0.5],
                                                              [1.6],
                                                              [-1.4]]))]))])

    scalings = OrderedDict([('O', OrderedDict([('intercept', -2.3),
                                               ('slope', 4.5)])),
                            ('Pd', OrderedDict([('intercept', 1.6),
                                                ('slope', 2.5)])),
                            ('Cu', OrderedDict([('intercept', -0.3),
                                                ('slope', -0.5)]))])

    fingerprints_range = {"Cu": np.array([[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]),
                          "O": np.array([[0.2139617720858539,
                                          2.258090276328769],
                                         [0.0, 1.085656080548734],
                                         [0.0, 0.0]]),
                          "Pd": np.array([[0.0, 1.4751761770313006],
                                          [0.0, 0.28464992134267897],
                                          [0.0, 0.20167521020630502]])}

    # Testing pure-python and fortran versions of Gaussian-neural force call
    device = "cpu" 
    model = core_test.AMPtorch(images, device,
            descriptor=Gaussian(cutoff=1000, Gs=Gs))

    criterion = nn.MSELoss()
    optimizer = optim.LBFGS
    RMSE=2e-3
    lr=0.8

    train=model.train(criterion,optimizer,lr,RMSE)
    sys.exit()


    """
    for fortran in [False, True]:
        for cores in range(1, 6):
            label = 'call-nonperiodic/%s-%i' % (fortran, cores)
            calc = Amp(descriptor=Gaussian(cutoff=6.5,
                                           Gs=Gs,
                                           fortran=fortran),
                       model=NeuralNetwork(hiddenlayers=hiddenlayers,
                                           weights=weights,
                                           scalings=scalings,
                                           activation='sigmoid',
                                           fprange=fingerprints_range,
                                           mode='atom-centered',
                                           fortran=fortran),
                       label=label,
                       dblabel=label,
                       cores=cores)

            predicted_energies = [calc.get_potential_energy(image) for image in
                                  images]

            for image_no in range(len(predicted_energies)):
                diff = abs(predicted_energies[image_no] -
                           correct_energies[image_no])
                assert (diff < 10.**(-15.)), \
                    'The predicted energy of image %i is wrong!' % (
                        image_no + 1)

            predicted_forces = [calc.get_forces(image) for image in images]

            for image_no in range(len(predicted_forces)):
                for index in range(np.shape(predicted_forces[image_no])[0]):
                    for direction in range(
                            np.shape(predicted_forces[image_no])[1]):
                        diff = abs(predicted_forces[image_no][index][
                            direction] -
                            correct_forces[image_no][index][direction])
                        assert (diff < 10.**(-15.)), \
                            'The predicted %i force of atom %i of image %i ' \
                            'is wrong!' % (direction, index, image_no + 1)

"""
def periodic_test():
    """Gaussian/Neural periodic standard.

    Checks that the answer matches that expected from previous Mathematica
    calculations.
    """

    # Making the list of periodic images
    images = [Atoms(symbols='PdOPd',
                    pbc=np.array([True, False, False], dtype=bool),
                    cell=np.array(
                        [[2.,  0.,  0.],
                         [0.,  2.,  0.],
                         [0.,  0.,  2.]]),
                    positions=np.array(
                        [[0.5,  1., 0.5],
                         [1.,  0.5,  1.],
                         [1.5,  1.5,  1.5]])),
              Atoms(symbols='PdO',
                    pbc=np.array([True, True, False], dtype=bool),
                    cell=np.array(
                        [[2.,  0.,  0.],
                         [0.,  2.,  0.],
                            [0.,  0.,  2.]]),
                    positions=np.array(
                        [[0.5,  1., 0.5],
                         [1.,  0.5,  1.]])),
              Atoms(symbols='Cu',
                    pbc=np.array([True, True, False], dtype=bool),
                    cell=np.array(
                        [[1.8,  0.,  0.],
                         [0.,  1.8,  0.],
                            [0.,  0.,  1.8]]),
                    positions=np.array(
                        [[0.,  0., 0.]]))]

    # Correct energies and forces
    correct_energies = [3.8560954326995978, 1.6120748520627273,
                        0.19433107801410093]
    correct_forces = \
        [[[0.14747720528015523, -3.3010645563584973, 3.3008168318984463],
          [0.03333579762326405, 9.050780376599887, -0.42608278400777605],
            [-0.1808130029034193, -5.7497158202413905, -2.8747340478906698]],
            [[6.5035267996045045 * (10.**(-6.)),
              -6.503526799604495 * (10.**(-6.)),
              0.00010834689201069249],
             [-6.5035267996045045 * (10.**(-6.)),
              6.503526799604495 * (10.**(-6.)),
              -0.00010834689201069249]],
            [[0.0, 0.0, 0.0]]]

    # Parameters
    Gs = {'O': [{'type': 'G2', 'element': 'Pd', 'eta': 0.8},
                {'type': 'G4', 'elements': ['O', 'Pd'], 'eta':0.3, 'gamma':0.6,
                 'zeta':0.5}],
          'Pd': [{'type': 'G2', 'element': 'Pd', 'eta': 0.2},
                 {'type': 'G4', 'elements': ['Pd', 'Pd'],
                  'eta':0.9, 'gamma':0.75, 'zeta':1.5}],
          'Cu': [{'type': 'G2', 'element': 'Cu', 'eta': 0.8},
                 {'type': 'G4', 'elements': ['Cu', 'Cu'], 'eta':0.3,
                          'gamma':0.6, 'zeta':0.5}]}

    hiddenlayers = {'O': (2,), 'Pd': (2,), 'Cu': (2,)}

    weights = OrderedDict([('O', OrderedDict([(1, np.matrix([[-2.0, 6.0],
                                                             [3.0, -3.0],
                                                             [1.5, -0.9]])),
                                              (2, np.matrix([[5.5],
                                                             [3.6],
                                                             [1.4]]))])),
                           ('Pd', OrderedDict([(1, np.matrix([[-1.0, 3.0],
                                                              [2.0, 4.2],
                                                              [1.0, -0.7]])),
                                               (2, np.matrix([[4.0],
                                                              [0.5],
                                                              [3.0]]))])),
                           ('Cu', OrderedDict([(1, np.matrix([[0.0, 1.0],
                                                              [-1.0, -2.0],
                                                              [2.5, -1.9]])),
                                               (2, np.matrix([[0.5],
                                                              [1.6],
                                                              [-1.4]]))]))])

    scalings = OrderedDict([('O', OrderedDict([('intercept', -2.3),
                                               ('slope', 4.5)])),
                            ('Pd', OrderedDict([('intercept', 1.6),
                                                ('slope', 2.5)])),
                            ('Cu', OrderedDict([('intercept', -0.3),
                                                ('slope', -0.5)]))])

    fingerprints_range = {"Cu": np.array([[2.8636310860653253,
                                           2.8636310860653253],
                                          [1.5435994865298275,
                                           1.5435994865298275]]),
                          "O": np.array([[2.9409056366723028,
                                          2.972494902604392],
                                         [1.9522542722823606,
                                          4.0720361595017245]]),
                          "Pd": np.array([[2.4629488092411096,
                                           2.6160138774087125],
                                          [0.27127576524253594,
                                           0.5898312261433813]])}

    # Testing pure-python and fortran versions of Gaussian-neural force call
    for fortran in [False, True]:
        for cores in range(1, 4):
            label = 'call-periodic/%s-%i' % (fortran, cores)
            calc = Amp(descriptor=Gaussian(cutoff=4.,
                                           Gs=Gs,
                                           fortran=fortran),
                       model=NeuralNetwork(hiddenlayers=hiddenlayers,
                                           weights=weights,
                                           scalings=scalings,
                                           activation='tanh',
                                           fprange=fingerprints_range,
                                           mode='atom-centered',
                                           fortran=fortran),
                       label=label,
                       dblabel=label,
                       cores=cores)

            predicted_energies = [calc.get_potential_energy(image) for image in
                                  images]

            for image_no in range(len(predicted_energies)):
                diff = abs(predicted_energies[image_no] -
                           correct_energies[image_no])
                assert (diff < 10.**(-14.)), \
                    'The predicted energy of image %i is wrong!' % (
                        image_no + 1)

            predicted_forces = [calc.get_forces(image) for image in images]

            for image_no in range(len(predicted_forces)):
                for index in range(np.shape(predicted_forces[image_no])[0]):
                    for direction in range(
                            np.shape(predicted_forces[image_no])[1]):
                        diff = abs(predicted_forces[image_no][index][
                            direction] -
                            correct_forces[image_no][index][direction])
                        assert (diff < 10.**(-11.)), \
                            'The predicted %i force of atom %i of image' \
                            ' %i is wrong!' % (direction,
                                               index,
                                               image_no + 1)

if __name__ == '__main__':
    non_periodic_test()
    # periodic_test()

