#include <math.h>
//#include "mpi.h"
#include "gmpordernorm.h"

extern "C" int calculate_gmpordernorm(double **, double **, double **, int*,
                                        int *, int, int*, int,
                                        int**, double **, int, double **, int *, int *,
                                        double**, double**);

extern "C" int calculate_gmpordernorm_noderiv(double **, double **, double **, int*,
                                        int *, int, int*, int,
                                        int**, double **, int, double **, int *, int *,
                                        double**);

extern "C" int calculate_solid_gmpordernorm(double **, double **, double **, int*,
                                        int *, int, int*, int,
                                        int**, double **, int, double **, int *, int *,
                                        double**, double**);

extern "C" int calculate_solid_gmpordernorm_noderiv(double **, double **, double **, int*,
                                        int *, int, int*, int,
                                        int**, double **, int, double **, int *, int *,
                                        double**); 


