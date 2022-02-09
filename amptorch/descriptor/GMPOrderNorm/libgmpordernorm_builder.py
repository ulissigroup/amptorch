import cffi

ffibuilder = cffi.FFI()
ffibuilder.cdef(
    """int calculate_gmpordernorm(double **, double **, double **, int*,
                        int *, int, int*, int,
                        int**, double **, int, double **, int*, int*,
                        double**, double**);

        int calculate_gmpordernorm_noderiv(double **, double **, double **, int*,
                                    int *, int, int*, int,
                                    int**, double **, int, double **, int*, int*,
                                    double**);

        int calculate_solid_gmpordernorm(double **, double **, double **, int*,
                                    int *, int, int*, int,
                                    int**, double **, int, double **, int*, int*,
                                    double**, double**);

        int calculate_solid_gmpordernorm_noderiv(double **, double **, double **, int*,
                                    int *, int, int*, int,
                                    int**, double **, int, double **, int*, int*,
                                    double**);


    """
)
ffibuilder.set_source(
    "amptorch.descriptor.GMPOrderNorm._libgmpordernorm",
    '#include "calculate_gmpordernorm.h"',
    sources=[
        "amptorch/descriptor/GMPOrderNorm/calculate_gmpordernorm.cpp",
        # "amptorch/descriptor/GMPOrderNorm/gmpordernorm.cpp",
        "amptorch/descriptor/GMPOrderNorm/helper.cpp",
        "amptorch/descriptor/GMPOrderNorm/surface_harmonics.cpp",
        "amptorch/descriptor/GMPOrderNorm/solid_harmonics.cpp",
    ],
    source_extension=".cpp",
    include_dirs=["amptorch/descriptor/GMPOrderNorm/"],
)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
