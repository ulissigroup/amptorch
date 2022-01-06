#include <math.h>
//#include "mpi.h"
// #include "gmpordernorm.h"
// #include "helper.h"
#include "surface_harmonics.h"
#include "solid_harmonics.h"

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

const int NUM_IMPLEMENTED_TYPE = 73;
const int IMPLEMENTED_MCSH_TYPE[][2] = {
    {0, 0},
    {1, 0},
    {2, 0},
    {3, 0},
    {4, 0},
    {5, 0},
    {6, 0},
    {7, 0},
    {8, 0},
    {9, 0},
    {0, 1},
    {1, 1},
    {2, 1},
    {3, 1},
    {4, 1},
    {5, 1},
    {6, 1},
    {7, 1},
    {8, 1},
    {9, 1},
    {0, 1},
    {1, 1},
    {2, 1},
    {2, 2},
    {3, 1},
    {3, 2},
    {3, 3},
    {4, 1},
    {4, 2},
    {4, 3},
    {4, 4},
    {5, 1},
    {5, 2},
    {5, 3},
    {5, 4},
    {5, 5},
    {6, 1},
    {6, 2},
    {6, 3},
    {6, 4},
    {6, 5},
    {6, 6},
    {6, 7},
    {7, 1},
    {7, 2},
    {7, 3},
    {7, 4},
    {7, 5},
    {7, 6},
    {7, 7},
    {7, 8},
    {8, 1},
    {8, 2},
    {8, 3},
    {8, 4},
    {8, 5},
    {8, 6},
    {8, 7},
    {8, 8},
    {8, 9},
    {8, 10},
    {9, 1},
    {9, 2},
    {9, 3},
    {9, 4},
    {9, 5},
    {9, 6},
    {9, 7},
    {9, 8},
    {9, 9},
    {9, 10},
    {9, 11},
    {9, 12}
}; // Change this when you implement new type!


