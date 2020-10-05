#include <math.h>
//#include "mpi.h"
#include "atomistic_mcsh.h"

extern "C" int calculate_atomistic_mcsh(double **, double **, double **, int*,
                                        int *, int, int*, int,
                                        int**, double **, int, double **, int *, int *,
                                        double**, double**);