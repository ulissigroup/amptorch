#include <math.h>
//#include "mpi.h"
#include "gmp.h"

extern "C" int calculate_gmp(double **, double **, double **, int*,
                                        int *, int, int*, int,
                                        int**, double **, int, double **, int *, int *,
                                        double**, double**);

extern "C" int calculate_gmp_noderiv(double **, double **, double **, int*,
                                        int *, int, int*, int,
                                        int**, double **, int, double **, int *, int *,
                                        double**);

extern "C" int calculate_gmp_square(double **, double **, double **, int*,
                                        int *, int, int*, int,
                                        int**, double **, int, double **, int *, int *,
                                        double**, double**);

extern "C" int calculate_gmp_square_noderiv(double **, double **, double **, int*,
                                        int *, int, int*, int,
                                        int**, double **, int, double **, int *, int *,
                                        double**);
