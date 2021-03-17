import cffi

ffibuilder = cffi.FFI()
ffibuilder.cdef(
    """int calculate_gmp(double **, double **, double **, int*,
                        int *, int, int*, int,
                        int**, double **, int, double **, int*, int*,
                        double**, double**);

        int calculate_gmp_square(double **, double **, double **, int*,
                                int *, int, int*, int,
                                int**, double **, int, double **, int*, int*,
                                double**, double**);

        int calculate_gmp_noderiv(double **, double **, double **, int*,
                                    int *, int, int*, int,
                                    int**, double **, int, double **, int*, int*,
                                    double**);

        int calculate_gmp_square_noderiv(double **, double **, double **, int*,
                                    int *, int, int*, int,
                                    int**, double **, int, double **, int*, int*,
                                    double**);
    """
)
ffibuilder.set_source(
    "amptorch.descriptor.GMP._libgmp",
    '#include "calculate_gmp.h"',
    sources=[
        "amptorch/descriptor/GMP/calculate_gmp.cpp",
        "amptorch/descriptor/GMP/gmp.cpp",
    ],
    source_extension=".cpp",
    include_dirs=["amptorch/descriptor/GMP/"],
)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
