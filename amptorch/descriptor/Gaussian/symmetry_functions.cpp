#define _USE_MATH_DEFINES
#include <math.h>
#include "symmetry_functions.h"

double sigm(double x, double &deriv) {
    double expl = 1./(1.+exp(-x));
    deriv = expl*(1-expl);
    return expl;
}

double cutf(double frac) {
    // frac = dist / cutoff_dist
    if (frac >= 1.0) {
        return 0;
    } else {
        return 0.5 * (1 + cos(M_PI*frac));
    }
}

double dcutf(double dist, double cutd) {
    if (dist/cutd >= 1.0) {
        return 0;
    } else {
        return -0.5 * M_PI * sin(M_PI*dist/cutd) / cutd;
    }
}

double poly_cutf(double frac, double gamma) {
    // frac = dist / cutoff_dist
    if (frac >= 1.0) {
        return 0;
    } else {
        return 1.0 + gamma * pow(frac, gamma+1) - (gamma+1) * pow(frac, gamma);
    }
}

double dpoly_cutf(double dist, double cutd, double gamma) {
    if (dist/cutd >= 1.0) {
        return 0;
    } else {
        return gamma * (gamma+1) / cutd * (pow(dist/cutd, gamma) - pow(dist/cutd, gamma-1));
    }
}



double G2(double Rij, double *precal, double *par, double &deriv) {
    // par[0] = cutoff_dist
    // par[1] = eta
    // par[2] = R_s
    double tmp = Rij-par[2];
    double expl = exp(-par[1]*tmp*tmp);
    deriv = expl*(-2*par[1]*tmp*precal[0] + precal[1]);
    return expl*precal[0];
}

double G4(double Rij, double Rik, double Rjk, double powtwo, \
          double *precal, double *par, double *deriv) {
    // cosv: cos(theta)
    // par[0] = cutoff_dist
    // par[1] = eta
    // par[2] = zeta
    // par[3] = lambda
    double expl = exp(-par[1]*precal[6]) * powtwo;
    double cosv = 1 + par[3]*precal[7];
    //double powcos = pow_int(cosv, par[2]-1);
    double powcos = pow(fabs(cosv), fabs(par[2]-1));

    deriv[0] = expl*powcos*precal[2]*precal[4] * \
               ((-2*par[1]*Rij*precal[0] + precal[1])*cosv + \
               par[2]*par[3]*precal[0]*precal[8]); // ij
    deriv[1] = expl*powcos*precal[0]*precal[4] * \
               ((-2*par[1]*Rik*precal[2] + precal[3])*cosv + \
               par[2]*par[3]*precal[2]*precal[9]); // ik
    deriv[2] = expl*powcos*precal[0]*precal[2] * \
               ((-2*par[1]*Rjk*precal[4] + precal[5])*cosv - \
               par[2]*par[3]*precal[4]*precal[10]); // jk

    return powcos*cosv * expl * precal[0] * precal[2] * precal[4];
}

double G5(double Rij, double Rik, double powtwo, \
          double *precal, double *par, double *deriv) {
    // cosv: cos(theta)
    // par[0] = cutoff_dist
    // par[1] = eta
    // par[2] = zeta
    // par[3] = lambda
    double expl = exp(-par[1]*precal[11]) * powtwo;
    double cosv = 1 + par[3]*precal[7];
    //double powcos = pow_int(cosv, par[2]-1);
    double powcos = pow(fabs(cosv), fabs(par[2]-1));

    deriv[0] = expl*powcos*precal[2] * \
               ((-2*par[1]*Rij*precal[0] + precal[1])*cosv + \
               par[2]*par[3]*precal[0]*precal[8]); // ij
    deriv[1] = expl*powcos*precal[0] * \
               ((-2*par[1]*Rik*precal[2] + precal[3])*cosv + \
               par[2]*par[3]*precal[2]*precal[9]); // ik
    deriv[2] = expl*powcos*precal[0]*precal[2] * \
               -par[2]*par[3]*precal[10]; // jk

    return powcos*cosv * expl * precal[0] * precal[2];
}


double G2_noderiv(double Rij, double *precal, double *par, double &deriv) {
    // par[0] = cutoff_dist
    // par[1] = eta
    // par[2] = R_s
    double tmp = Rij-par[2];
    double expl = exp(-par[1]*tmp*tmp);
    return expl*precal[0];
}

double G4_noderiv(double Rij, double Rik, double Rjk, double powtwo, \
          double *precal, double *par, double *deriv) {
    // cosv: cos(theta)
    // par[0] = cutoff_dist
    // par[1] = eta
    // par[2] = zeta
    // par[3] = lambda
    double expl = exp(-par[1]*precal[6]) * powtwo;
    double cosv = 1 + par[3]*precal[7];
    //double powcos = pow_int(cosv, par[2]-1);
    double powcos = pow(fabs(cosv), fabs(par[2]-1));

    return powcos*cosv * expl * precal[0] * precal[2] * precal[4];
}

double G5_noderiv(double Rij, double Rik, double powtwo, \
          double *precal, double *par, double *deriv) {
    // cosv: cos(theta)
    // par[0] = cutoff_dist
    // par[1] = eta
    // par[2] = zeta
    // par[3] = lambda
    double expl = exp(-par[1]*precal[11]) * powtwo;
    double cosv = 1 + par[3]*precal[7];
    //double powcos = pow_int(cosv, par[2]-1);
    double powcos = pow(fabs(cosv), fabs(par[2]-1));

    return powcos*cosv * expl * precal[0] * precal[2];
}
