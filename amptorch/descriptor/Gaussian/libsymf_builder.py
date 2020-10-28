import cffi

ffibuilder = cffi.FFI()
ffibuilder.cdef(
    """int calculate_sf_cos(double **, double **, double **, int*,
                            int *, int, int*, int,
                            int**, double **, int,
                            double**, double**);

       int calculate_sf_cos_noderiv(double **, double **, double **, int*,
                            int *, int, int*, int,
                            int**, double **, int,
                            double**);

       int calculate_sf_poly(double **, double **, double **, int*,
                            int *, int, int*, int,
                            int**, double **, int,
                            double**, double**, double);

       int calculate_sf_poly_noderiv(double **, double **, double **, int*,
                            int *, int, int*, int,
                            int**, double **, int,
                            double**, double);"""
)
ffibuilder.set_source(
    "amptorch.descriptor.Gaussian._libsymf",
    '#include "calculate_sf.h"',
    sources=[
        "amptorch/descriptor/Gaussian/calculate_sf.cpp",
        "amptorch/descriptor/Gaussian/symmetry_functions.cpp",
    ],
    source_extension=".cpp",
    include_dirs=["amptorch/descriptor/Gaussian/"],
)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
