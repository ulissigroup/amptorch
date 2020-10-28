/*
 Code for calculating symmetry function value
 For that, this code recieve the atomic coordinates and calculate nearest neighbor.
 In LAMMPS code, this part is included in pair_nn.cpp
 */

#include <math.h>
//#include "mpi.h"
#include "symmetry_functions.h"

extern "C" int calculate_sf_cos(double **, double **, double **, int*,
                            int *, int, int*, int,
                            int**, double **, int,
                            double**, double**);

extern "C" int calculate_sf_cos_noderiv(double **, double **, double **, int*,
                            int *, int, int*, int,
                            int**, double **, int,
                            double**);

extern "C" int calculate_sf_poly(double **, double **, double **, int*,
                            int *, int, int*, int,
                            int**, double **, int,
                            double**, double**, double);

extern "C" int calculate_sf_poly_noderiv(double **, double **, double **, int*,
                            int *, int, int*, int,
                            int**, double **, int,
                            double**, double);
