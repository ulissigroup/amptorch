import cffi

ffibuilder = cffi.FFI()
ffibuilder.cdef(
    """
        int calculate_gmp_align_noderiv(double **, double **, double **, int*,
                                    int *, int, int*, int,
                                    int**, double **, int, int, double **, int*, int*,
                                    double**);

        int calculate_gmp_no_align_noderiv(double **, double **, double **, int*,
                                    int *, int, int*, int,
                                    int**, double **, int, int, double **, int*, int*,
                                    double**);

    """
)
ffibuilder.set_source(
    "amptorch.descriptor.GMP_align._libgmpalign",
    '#include "calculate_gmp_align.h"',
    sources=[
        "amptorch/descriptor/GMP_align/calculate_gmp_align.cpp",
        "amptorch/descriptor/GMP_align/gmp_align.cpp",
    ],
    source_extension=".cpp",
    include_dirs=["amptorch/descriptor/GMP_align/"],
)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
