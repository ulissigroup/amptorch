/*
 Code for calculate symmetry function.
 This code is used for both Python and LAMMPS code
 */
const int IMPLEMENTED_TYPE[] = {2, 4, 5}; // Change this when you implement new symfunc type!

static inline double pow_int(const double &x, const double n) {
    double res,tmp;

    if (x == 0.0) return 0.0; // FIXME: abs(x) < epsilon
    int nn = (n > 0) ? n : -n;
    tmp = x;

    for (res = 1.0; nn != 0; nn >>= 1, tmp *= tmp)
        if (nn & 1) res *= tmp;

    return (n > 0) ? res : 1.0/res;
}

double sigm(double, double &);
double cutf(double);
double dcutf(double, double);
double poly_cutf(double, double);
double dpoly_cutf(double, double, double);
double G2(double, double *, double *, double &);
double G4(double, double, double, double, double *, double *, double *);
double G5(double, double, double, double *, double *, double *);

double G2_noderiv(double, double *, double *, double &);
double G4_noderiv(double, double, double, double, double *, double *, double *);
double G5_noderiv(double, double, double, double *, double *, double *);
