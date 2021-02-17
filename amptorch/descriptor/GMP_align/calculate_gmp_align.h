#include <math.h>
//#include "mpi.h"
#include "gmp_align.h"

// extern "C" int calculate_atomistic_mcsh(double **, double **, double **, int*,
//                                         int *, int, int*, int,
//                                         int**, double **, int, double **, int *, int *,
//                                         double**, double**);

extern "C" int calculate_gmp_align_noderiv(double **, double **, double **, int*,
                                        int *, int, int*, int,
                                        int**, double **, int, int, double **, int *, int *,
                                        double**);

extern "C" int calculate_gmp_no_align_noderiv(double **, double **, double **, int*,
                                        int *, int, int*, int,
                                        int**, double **, int, int, double **, int *, int *,
                                        double**);

// extern "C" int calculate_atomistic_mcsh_square(double **, double **, double **, int*,
//                                         int *, int, int*, int,
//                                         int**, double **, int, double **, int *, int *,
//                                         double**, double**);

// extern "C" int calculate_atomistic_mcsh_square_noderiv(double **, double **, double **, int*,
//                                         int *, int, int*, int,
//                                         int**, double **, int, double **, int *, int *,
//                                         double**);
