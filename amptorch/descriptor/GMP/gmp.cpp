#include <math.h>
#include <stdio.h>
#include "gmp.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif


double calc_C1(double A, double B, double alpha, double beta){
    double temp = sqrt(M_PI / (alpha + beta));
    return A * B * temp * temp * temp;
}

double calc_C2(double alpha, double beta){
    return -1.0 * (alpha * beta / (alpha + beta));
}


double calc_lambda(double alpha, double beta){
    // return  alpha / (alpha + beta);
    return beta / (alpha + beta);
}

double calc_gamma(double alpha, double beta){
    return alpha + beta;
}

// double P2(double lambda_x0_2, double gamma){
//     return (0.5 / gamma) + lambda_x0_2;
// }

// double P3(double lambda_x0_3, double lambda_x0, double gamma){
//     return (1.5 * lambda_x0 / gamma) + lambda_x0_4;
// }

// double P4(double lambda_x0_4, double lambda_x0_2, double gamma){
//     return (0.75 / (gamma * gamma)) + (3.0 * lambda_x0_2 / gamma) + lambda_x0_4;
// }

// double P5(double lambda_x0_5, double lambda_x0_3, double lambda_x0, double gamma){
//     return ((15.0 * lambda_x0) / (4.0 * gamma * gamma)) + (5.0 * lambda_x0_3 / gamma) + lambda_x0_5;
// }

// double P6(double lambda_x0_6, double lambda_x0_4, double lambda_x0_2, double gamma){
//     return (15.0 / (8.0 * gamma * gamma * gamma)) + ((11.25 * lambda_x0_2) / (gamma * gamma)) + (7.5 * lambda_x0_4 / gamma) + lambda_x0_6;
// }

// double P7(double lambda_x0_7, double lambda_x0_5, double lambda_x0_3, double lambda_x0, double gamma){
//     double term1 = 105.0 * lambda_x0 / (8.0 * gamma * gamma * gamma);
//     double term2 = 105.0 * lambda_x0_3 / (4.0 * gamma * gamma);
//     double term3 = 10.5 * lambda_x0_5 / gamma;
//     return term1 + term2 + term3 + lambda_x0_7;
// }

// double P8(double lambda_x0_8, double lambda_x0_6, double lambda_x0_4, double lambda_x0_2, double gamma){
//     double term1 = 105.0 / (16.0 * gamma * gamma * gamma * gamma);
//     double term2 = 105.0 * lambda_x0_2 / (2.0 * gamma * gamma * gamma) ;
//     double term3 = 105.0 * lambda_x0_4 / (2.0 * gamma * gamma);
//     double term4 = 14.0 * lambda_x0_6 / gamma;
//     return term1 + term2 + term3 + term4 + lambda_x0_8;
// }

// double P9(double lambda_x0_7, double lambda_x0_5, double lambda_x0_3, double lambda_x0, double gamma){
//     double term1 = 945.0 * lambda_x0 / (16.0 * gamma * gamma * gamma * gamma);
//     double term2 = 315.0 * lambda_x0_3 / (2.0 * gamma * gamma * gamma);
//     double term3 = 189.0 * lambda_x0_5 / (2.0 * gamma * gamma);
//     double term4 = 18.0 * lambda_x0_7 / gamma;
//     return term1 + term2 + term3 + term4 + lambda_x0_9;
// }

void calc_MCSH_0_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double m_0_1 = C1 * exp( C2 * r0_sqr);

    deriv[0] = m_0_1 * (2.0 * C2 * x0);
    deriv[1] = m_0_1 * (2.0 * C2 * y0);
    deriv[2] = m_0_1 * (2.0 * C2 * z0);

    value[0] = m_0_1;
}

void calc_MCSH_1_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double lambda = calc_lambda(alpha, beta);

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = lambda * x0;
    double temp_y = lambda * y0;
    double temp_z = lambda * z0;

    double miu_1_1_1 = temp * temp_x * inv_rs;
    double miu_1_1_2 = temp * temp_y * inv_rs;
    double miu_1_1_3 = temp * temp_z * inv_rs;

    double temp_dx = lambda;
    double temp_dy = lambda;
    double temp_dz = lambda;

    //deriv consts
    double const_2_C2_x = 2.0 * C2 * x0;
    double const_2_C2_y = 2.0 * C2 * y0;
    double const_2_C2_z = 2.0 * C2 * z0;

    // dmiu1 dx/dy/dz
    deriv[0] = temp * inv_rs * (temp_dx + const_2_C2_x * temp_x);
    deriv[1] = miu_1_1_1 * const_2_C2_y;
    deriv[2] = miu_1_1_1 * const_2_C2_z;

    // dmiu2 dx/dy/dz
    deriv[3] = miu_1_1_2 * const_2_C2_x;
    deriv[4] = temp * inv_rs * (temp_dy + const_2_C2_y * temp_y);
    deriv[5] = miu_1_1_2 * const_2_C2_z;

    // dmiu3 dx/dy/dz
    deriv[6] = miu_1_1_3 * const_2_C2_x;
    deriv[7] = miu_1_1_3 * const_2_C2_y;
    deriv[8] = temp * inv_rs * (temp_dz + const_2_C2_z * temp_z);

    value[0] = miu_1_1_1;
    value[1] = miu_1_1_2;
    value[2] = miu_1_1_3;

}

void calc_MCSH_2_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double inv_rs_2 = inv_rs * inv_rs;

    double lambda_x0 = x0 * lambda;
    double lambda_y0 = y0 * lambda;
    double lambda_z0 = z0 * lambda;

    double lambda_x0_sqr = lambda_x0 * lambda_x0;
    double lambda_y0_sqr = lambda_y0 * lambda_y0;
    double lambda_z0_sqr = lambda_z0 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = ((3.0 * inv_rs_2) / (2.0 * gamma)) - 1.0;

    double temp = C1 * exp( C2 * r0_sqr);


    double temp_x = 3.0 * lambda_x0_sqr * inv_rs_2 + C3;
    double temp_y = 3.0 * lambda_y0_sqr * inv_rs_2 + C3;
    double temp_z = 3.0 * lambda_z0_sqr * inv_rs_2 + C3;

    double temp_dx = 6.0 * inv_rs_2 * lambda * lambda_x0; // = lambda * (3 * 2 * lambda_x0)
    double temp_dy = 6.0 * inv_rs_2 * lambda * lambda_y0;
    double temp_dz = 6.0 * inv_rs_2 * lambda * lambda_z0;

    double miu_2_1_1 = temp * temp_x;
    double miu_2_1_2 = temp * temp_y;
    double miu_2_1_3 = temp * temp_z;

    //deriv consts
    double const_2_C2_x = 2.0 * C2 * x0;
    double const_2_C2_y = 2.0 * C2 * y0;
    double const_2_C2_z = 2.0 * C2 * z0;

    // dmiu1 dx/dy/dz
    deriv[0] = temp * (temp_dx + const_2_C2_x * temp_x);
    deriv[1] = miu_2_1_1 * const_2_C2_y;
    deriv[2] = miu_2_1_1 * const_2_C2_z;

    // dmiu2 dx/dy/dz
    deriv[3] = miu_2_1_2 * const_2_C2_x;
    deriv[4] = temp * (temp_dy + const_2_C2_y * temp_y);
    deriv[5] = miu_2_1_2 * const_2_C2_z;

    // dmiu3 dx/dy/dz
    deriv[6] = miu_2_1_3 * const_2_C2_x;
    deriv[7] = miu_2_1_3 * const_2_C2_y;
    deriv[8] = temp * (temp_dz + const_2_C2_z * temp_z);

    value[0] = miu_2_1_1;
    value[1] = miu_2_1_2;
    value[2] = miu_2_1_3;
}

void calc_MCSH_2_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double inv_rs_2 = inv_rs * inv_rs;
    double x0_sqr = x0*x0;
    double y0_sqr = y0*y0;
    double z0_sqr = z0*z0;

    double temp = C1 * exp( C2 * r0_sqr) * lambda * lambda * 3.0 * inv_rs_2;

    double miu_2_2_1 = temp * x0 * y0;
    double miu_2_2_2 = temp * x0 * z0;
    double miu_2_2_3 = temp * y0 * z0;

    //deriv consts
    double const_2_C2_x = 2.0 * C2 * x0; // 2 * C2 * x0
    double const_2_C2_y = 2.0 * C2 * y0;
    double const_2_C2_z = 2.0 * C2 * z0;

    double const_1_p_2_C2_x2 = 1.0 + 2.0 * C2 * x0_sqr; // 1 + 2*C2*x0*x0
    double const_1_p_2_C2_y2 = 1.0 + 2.0 * C2 * y0_sqr;
    double const_1_p_2_C2_z2 = 1.0 + 2.0 * C2 * z0_sqr;

    // dmiu1 dx/dy/dz
    deriv[0] = temp * y0 * const_1_p_2_C2_x2;
    deriv[1] = temp * x0 * const_1_p_2_C2_y2;
    deriv[2] = miu_2_2_1 * const_2_C2_z;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * z0 * const_1_p_2_C2_x2;
    deriv[4] = miu_2_2_2 * const_2_C2_y;
    deriv[5] = temp * x0 * const_1_p_2_C2_z2;

    // dmiu3 dx/dy/dz
    deriv[6] = miu_2_2_3 * const_2_C2_x;
    deriv[7] = temp * z0 * const_1_p_2_C2_y2;
    deriv[8] = temp * y0 * const_1_p_2_C2_z2;

    value[0] = miu_2_2_1;
    value[1] = miu_2_2_2;
    value[2] = miu_2_2_3;
}

void calc_MCSH_3_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double inv_rs_3 = inv_rs * inv_rs * inv_rs;

    double lambda_x0 = x0 * lambda;
    double lambda_y0 = y0 * lambda;
    double lambda_z0 = z0 * lambda;

    double lambda_x0_sqr = lambda_x0 * lambda_x0;
    double lambda_y0_sqr = lambda_y0 * lambda_y0;
    double lambda_z0_sqr = lambda_z0 * lambda_z0;

    double lambda_x0_3 = lambda_x0_sqr * lambda_x0;
    double lambda_y0_3 = lambda_y0_sqr * lambda_y0;
    double lambda_z0_3 = lambda_z0_sqr * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = ((45.0 * inv_rs_3) / (2.0 * gamma)) - (9.0 * inv_rs);

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 15.0 * inv_rs_3 * lambda_x0_3 + C3 * lambda_x0;
    double temp_y = 15.0 * inv_rs_3 * lambda_y0_3 + C3 * lambda_y0;
    double temp_z = 15.0 * inv_rs_3 * lambda_z0_3 + C3 * lambda_z0;

    double temp_dx = lambda * (45.0 * inv_rs_3 * lambda_x0_sqr + C3);
    double temp_dy = lambda * (45.0 * inv_rs_3 * lambda_y0_sqr + C3);
    double temp_dz = lambda * (45.0 * inv_rs_3 * lambda_z0_sqr + C3);

    double miu_3_1_1 = temp * temp_x;
    double miu_3_1_2 = temp * temp_y;
    double miu_3_1_3 = temp * temp_z;

    //deriv consts
    double const_2_C2_x = 2.0 * C2 * x0;
    double const_2_C2_y = 2.0 * C2 * y0;
    double const_2_C2_z = 2.0 * C2 * z0;

    // dmiu1 dx/dy/dz
    deriv[0] = temp * (temp_dx + const_2_C2_x * temp_x);
    deriv[1] = miu_3_1_1 * const_2_C2_y;
    deriv[2] = miu_3_1_1 * const_2_C2_z;

    // dmiu2 dx/dy/dz
    deriv[3] = miu_3_1_2 * const_2_C2_x;
    deriv[4] = temp * (temp_dy + const_2_C2_y * temp_y);
    deriv[5] = miu_3_1_2 * const_2_C2_z;

    // dmiu3 dx/dy/dz
    deriv[6] = miu_3_1_3 * const_2_C2_x;
    deriv[7] = miu_3_1_3 * const_2_C2_y;
    deriv[8] = temp * (temp_dz + const_2_C2_z * temp_z);

    value[0] = miu_3_1_1;
    value[1] = miu_3_1_2;
    value[2] = miu_3_1_3;
}

void calc_MCSH_3_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double inv_rs_3 = inv_rs * inv_rs * inv_rs;
    double x0_sqr = x0*x0;
    double y0_sqr = y0*y0;
    double z0_sqr = z0*z0;

    double lambda_x0 = x0 * lambda;
    double lambda_y0 = y0 * lambda;
    double lambda_z0 = z0 * lambda;

    double lambda_x0_sqr = lambda_x0 * lambda_x0;
    double lambda_y0_sqr = lambda_y0 * lambda_y0;
    double lambda_z0_sqr = lambda_z0 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = ((15.0 * inv_rs_3) / (2.0 * gamma)) - (3.0 * inv_rs);

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 15.0 * inv_rs_3 * lambda_x0_sqr + C3;
    double temp_y = 15.0 * inv_rs_3 * lambda_y0_sqr + C3;
    double temp_z = 15.0 * inv_rs_3 * lambda_z0_sqr + C3;

    double temp_dx = lambda * (30.0 * inv_rs_3 * lambda_x0);
    double temp_dy = lambda * (30.0 * inv_rs_3 * lambda_y0);
    double temp_dz = lambda * (30.0 * inv_rs_3 * lambda_z0);

    double miu_3_2_1 = temp * lambda_y0 * temp_x;
    double miu_3_2_2 = temp * lambda_x0 * temp_y;
    double miu_3_2_3 = temp * lambda_z0 * temp_x;
    double miu_3_2_4 = temp * lambda_x0 * temp_z;
    double miu_3_2_5 = temp * lambda_z0 * temp_y;
    double miu_3_2_6 = temp * lambda_y0 * temp_z;

    //deriv consts
    double const_2_C2_x = 2.0 * C2 * x0; // 2 * C2 * x0
    double const_2_C2_y = 2.0 * C2 * y0;
    double const_2_C2_z = 2.0 * C2 * z0;

    double const_1_p_2_C2_x2 = 1.0 + 2.0 * C2 * x0_sqr; // 1 + 2*C2*x0*x0
    double const_1_p_2_C2_y2 = 1.0 + 2.0 * C2 * y0_sqr;
    double const_1_p_2_C2_z2 = 1.0 + 2.0 * C2 * z0_sqr;

    // dmiu1 dx/dy/dz
    deriv[0] = temp * lambda_y0 * (temp_dx + const_2_C2_x * temp_x);
    deriv[1] = temp * lambda * temp_x * const_1_p_2_C2_y2;
    deriv[2] = miu_3_2_1 * const_2_C2_z;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * lambda * temp_y * const_1_p_2_C2_x2;
    deriv[4] = temp * lambda_x0 * (temp_dy + const_2_C2_y * temp_y);
    deriv[5] = miu_3_2_2 * const_2_C2_z;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * lambda_z0 * (temp_dx + const_2_C2_x * temp_x);
    deriv[7] = miu_3_2_3 * const_2_C2_y;
    deriv[8] = temp * lambda * temp_x * const_1_p_2_C2_z2;

    // dmiu4 dx/dy/dz
    deriv[9] = temp * lambda * temp_z * const_1_p_2_C2_x2;
    deriv[10] = miu_3_2_4 * const_2_C2_y;
    deriv[11] = temp * lambda_x0 * (temp_dz + const_2_C2_z * temp_z);

    // dmiu5 dx/dy/dz
    deriv[12] = miu_3_2_5 * const_2_C2_x;
    deriv[13] = temp * lambda_z0 * (temp_dy + const_2_C2_y * temp_y);
    deriv[14] = temp * lambda * temp_y * const_1_p_2_C2_z2;

    // dmiu6 dx/dy/dz
    deriv[15] = miu_3_2_6 * const_2_C2_x;
    deriv[16] = temp * lambda * temp_z * const_1_p_2_C2_y2;
    deriv[17] = temp * lambda_y0 * (temp_dz + const_2_C2_z * temp_z);

    value[0] = miu_3_2_1;
    value[1] = miu_3_2_2;
    value[2] = miu_3_2_3;
    value[3] = miu_3_2_4;
    value[4] = miu_3_2_5;
    value[5] = miu_3_2_6;
}

void calc_MCSH_3_3(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double lambda = calc_lambda(alpha, beta);
    double inv_rs_3 = inv_rs * inv_rs * inv_rs;

    double temp =  C1 * exp( C2 * r0_sqr) * lambda * lambda * lambda * 15.0 * inv_rs_3;
    double m_3_3 = temp * x0 * y0 * z0;

    deriv[0] = temp * y0 * z0 * (1.0 + 2.0*C2*x0*x0);
    deriv[1] = temp * x0 * z0 * (1.0 + 2.0*C2*y0*y0);
    deriv[2] = temp * x0 * y0 * (1.0 + 2.0*C2*z0*z0);

    value[0] = m_3_3;
}

void calc_MCSH_4_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double inv_rs_2 = inv_rs * inv_rs;
    double inv_rs_4 = inv_rs_2 * inv_rs_2;

    double lambda_x0 = lambda * x0;
    double lambda_y0 = lambda * y0;
    double lambda_z0 = lambda * z0;

    double lambda_x0_sqr = lambda_x0 * lambda_x0;
    double lambda_y0_sqr = lambda_y0 * lambda_y0;
    double lambda_z0_sqr = lambda_z0 * lambda_z0;

    double lambda_x0_3 = lambda_x0_sqr * lambda_x0;
    double lambda_y0_3 = lambda_y0_sqr * lambda_y0;
    double lambda_z0_3 = lambda_z0_sqr * lambda_z0;

    double lambda_x0_4 = lambda_x0_3 * lambda_x0;
    double lambda_y0_4 = lambda_y0_3 * lambda_y0;
    double lambda_z0_4 = lambda_z0_3 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = ((315.0 * inv_rs_4) / gamma) - (90.0 * inv_rs_2);
    double C4 = ((315.0 * inv_rs_4) / (4.0*gamma*gamma)) - (45.0 * inv_rs_2 / gamma) + 9.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 105.0 * inv_rs_4 * lambda_x0_4 + C3 * lambda_x0_sqr + C4;
    double temp_y = 105.0 * inv_rs_4 * lambda_y0_4 + C3 * lambda_y0_sqr + C4;
    double temp_z = 105.0 * inv_rs_4 * lambda_z0_4 + C3 * lambda_z0_sqr + C4;

    double temp_dx = lambda * (420.0 * inv_rs_4 * lambda_x0_3 + 2.0 * C3 * lambda_x0);
    double temp_dy = lambda * (420.0 * inv_rs_4 * lambda_y0_3 + 2.0 * C3 * lambda_y0);
    double temp_dz = lambda * (420.0 * inv_rs_4 * lambda_z0_3 + 2.0 * C3 * lambda_z0);

    double miu_4_1_1 = temp * temp_x;
    double miu_4_1_2 = temp * temp_y;
    double miu_4_1_3 = temp * temp_z;

    //deriv consts
    double const_2_C2_x = 2.0 * C2 * x0;
    double const_2_C2_y = 2.0 * C2 * y0;
    double const_2_C2_z = 2.0 * C2 * z0;

    // dmiu1 dx/dy/dz
    deriv[0] = temp * (temp_dx + const_2_C2_x * temp_x);
    deriv[1] = miu_4_1_1 * const_2_C2_y;
    deriv[2] = miu_4_1_1 * const_2_C2_z;

    // dmiu2 dx/dy/dz
    deriv[3] = miu_4_1_2 * const_2_C2_x;
    deriv[4] = temp * (temp_dy + const_2_C2_y * temp_y);
    deriv[5] = miu_4_1_2 * const_2_C2_z;

    // dmiu3 dx/dy/dz
    deriv[6] = miu_4_1_3 * const_2_C2_x;
    deriv[7] = miu_4_1_3 * const_2_C2_y;
    deriv[8] = temp * (temp_dz + const_2_C2_z * temp_z);

    value[0] = miu_4_1_1;
    value[1] = miu_4_1_2;
    value[2] = miu_4_1_3;
}

void calc_MCSH_4_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double inv_rs_2 = inv_rs * inv_rs;
    double inv_rs_4 = inv_rs_2 * inv_rs_2;
    double x0_sqr = x0*x0;
    double y0_sqr = y0*y0;
    double z0_sqr = z0*z0;

    double lambda_x0 = lambda * x0;
    double lambda_y0 = lambda * y0;
    double lambda_z0 = lambda * z0;

    double lambda_x0_sqr = lambda_x0 * lambda_x0;
    double lambda_y0_sqr = lambda_y0 * lambda_y0;
    double lambda_z0_sqr = lambda_z0 * lambda_z0;

    double lambda_x0_3 = lambda_x0_sqr * lambda_x0;
    double lambda_y0_3 = lambda_y0_sqr * lambda_y0;
    double lambda_z0_3 = lambda_z0_sqr * lambda_z0;

    double gamma = calc_gamma(alpha, beta);

    double C3 = ((315.0 * inv_rs_4) / (2.0 * gamma)) - (45.0 * inv_rs_2);

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 105.0 * inv_rs_4 * lambda_x0_3 + C3 * lambda_x0;
    double temp_y = 105.0 * inv_rs_4 * lambda_y0_3 + C3 * lambda_y0;
    double temp_z = 105.0 * inv_rs_4 * lambda_z0_3 + C3 * lambda_z0;

    double temp_dx = lambda * (315.0 * inv_rs_4 * lambda_x0_sqr + C3);// lambda * (105.0 * 3.0 * lambda_x0_sqr + C3);
    double temp_dy = lambda * (315.0 * inv_rs_4 * lambda_y0_sqr + C3);
    double temp_dz = lambda * (315.0 * inv_rs_4 * lambda_z0_sqr + C3);

    double miu_4_2_1 = temp * lambda_y0 * temp_x;
    double miu_4_2_2 = temp * lambda_x0 * temp_y;
    double miu_4_2_3 = temp * lambda_z0 * temp_x;
    double miu_4_2_4 = temp * lambda_x0 * temp_z;
    double miu_4_2_5 = temp * lambda_z0 * temp_y;
    double miu_4_2_6 = temp * lambda_y0 * temp_z;

    //deriv consts
    double const_2_C2_x = 2.0 * C2 * x0; // 2 * C2 * x0
    double const_2_C2_y = 2.0 * C2 * y0;
    double const_2_C2_z = 2.0 * C2 * z0;

    double const_1_p_2_C2_x2 = 1.0 + 2.0 * C2 * x0_sqr; // 1 + 2*C2*x0*x0
    double const_1_p_2_C2_y2 = 1.0 + 2.0 * C2 * y0_sqr;
    double const_1_p_2_C2_z2 = 1.0 + 2.0 * C2 * z0_sqr;

    // dmiu1 dx/dy/dz
    deriv[0] = temp * lambda_y0 * (temp_dx + const_2_C2_x * temp_x);
    deriv[1] = temp * lambda * temp_x * const_1_p_2_C2_y2;
    deriv[2] = miu_4_2_1 * const_2_C2_z;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * lambda * temp_y * const_1_p_2_C2_x2;
    deriv[4] = temp * lambda_x0 * (temp_dy + const_2_C2_y * temp_y);
    deriv[5] = miu_4_2_2 * const_2_C2_z;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * lambda_z0 * (temp_dx + const_2_C2_x * temp_x);
    deriv[7] = miu_4_2_3 * const_2_C2_y;
    deriv[8] = temp * lambda * temp_x * const_1_p_2_C2_z2;

    // dmiu4 dx/dy/dz
    deriv[9] = temp * lambda * temp_z * const_1_p_2_C2_x2;
    deriv[10] = miu_4_2_4 * const_2_C2_y;
    deriv[11] = temp * lambda_x0 * (temp_dz + const_2_C2_z * temp_z);

    // dmiu5 dx/dy/dz
    deriv[12] = miu_4_2_5 * const_2_C2_x;
    deriv[13] = temp * lambda_z0 * (temp_dy + const_2_C2_y * temp_y);
    deriv[14] = temp * lambda * temp_y * const_1_p_2_C2_z2;

    // dmiu6 dx/dy/dz
    deriv[15] = miu_4_2_6 * const_2_C2_x;
    deriv[16] = temp * lambda * temp_z * const_1_p_2_C2_y2;
    deriv[17] = temp * lambda_y0 * (temp_dz + const_2_C2_z * temp_z);

    value[0] = miu_4_2_1;
    value[1] = miu_4_2_2;
    value[2] = miu_4_2_3;
    value[3] = miu_4_2_4;
    value[4] = miu_4_2_5;
    value[5] = miu_4_2_6;
}


void calc_MCSH_4_3(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double lambda_sqr = lambda * lambda;
    double inv_rs_2 = inv_rs * inv_rs;
    double inv_rs_4 = inv_rs_2 * inv_rs_2;
    double x0_sqr = x0*x0;
    double y0_sqr = y0*y0;
    double z0_sqr = z0*z0;

    double lambda_x0_sqr = lambda_sqr * x0_sqr;
    double lambda_y0_sqr = lambda_sqr * y0_sqr;
    double lambda_z0_sqr = lambda_sqr * z0_sqr;

    double gamma = calc_gamma(alpha, beta);

    double temp_x_2 = lambda_x0_sqr + (1.0/(2.0*gamma));
    double temp_y_2 = lambda_y0_sqr + (1.0/(2.0*gamma));
    double temp_z_2 = lambda_z0_sqr + (1.0/(2.0*gamma));

    double temp_dx_2 = 2.0 * lambda_sqr * x0;
    double temp_dy_2 = 2.0 * lambda_sqr * y0;
    double temp_dz_2 = 2.0 * lambda_sqr * z0;

    double temp_term1_x = 105.0 * inv_rs_4 * temp_x_2 - (15.0 * inv_rs_2);
    double temp_term1_y = 105.0 * inv_rs_4 * temp_y_2 - (15.0 * inv_rs_2);
    // double temp_term1_z = 105.0 * temp_z_2 - 15.0;

    double temp_dterm1_dx = 105.0 * inv_rs_4 * temp_dx_2;
    double temp_dterm1_dy = 105.0 * inv_rs_4 * temp_dy_2;
    // double temp_dterm1_dz = 105.0 * temp_dz_2;

    double temp_term2_x = -15.0 * inv_rs_2 * temp_x_2 + 3.0;
    double temp_term2_y = -15.0 * inv_rs_2 * temp_y_2 + 3.0;
    // double temp_term2_z = -15.0 * temp_z_2 + 3.0;

    double temp_dterm2_dx = -15.0 * inv_rs_2 * temp_dx_2;
    double temp_dterm2_dy = -15.0 * inv_rs_2 * temp_dy_2;
    // double temp_dterm2_dz = -15.0 * temp_dz_2;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_miu1 = temp_y_2 * temp_term1_x + temp_term2_x;
    double temp_miu2 = temp_z_2 * temp_term1_x + temp_term2_x;
    double temp_miu3 = temp_z_2 * temp_term1_y + temp_term2_y;

    double temp_dmiu1_dx = temp_y_2 * temp_dterm1_dx + temp_dterm2_dx;
    double temp_dmiu1_dy = temp_dy_2 * temp_term1_x;

    double temp_dmiu2_dx = temp_z_2 * temp_dterm1_dx + temp_dterm2_dx;
    double temp_dmiu2_dz = temp_dz_2 * temp_term1_x;

    double temp_dmiu3_dy = temp_z_2 * temp_dterm1_dy + temp_dterm2_dy;
    double temp_dmiu3_dz = temp_dz_2 * temp_term1_y;

    double miu_4_3_1 = temp * temp_miu1;
    double miu_4_3_2 = temp * temp_miu2;
    double miu_4_3_3 = temp * temp_miu3;

    //deriv consts
    double const_2_C2_x = 2.0 * C2 * x0; // 2 * C2 * x0
    double const_2_C2_y = 2.0 * C2 * y0;
    double const_2_C2_z = 2.0 * C2 * z0;

    // dmiu1 dx/dy/dz
    deriv[0] = temp * (temp_dmiu1_dx + const_2_C2_x * temp_miu1);
    deriv[1] = temp * (temp_dmiu1_dy + const_2_C2_y * temp_miu1);
    deriv[2] = miu_4_3_1 * const_2_C2_z;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * (temp_dmiu2_dx + const_2_C2_x * temp_miu2);
    deriv[4] = miu_4_3_2 * const_2_C2_y;
    deriv[5] = temp * (temp_dmiu2_dz + const_2_C2_z * temp_miu2);

    // dmiu3 dx/dy/dz
    deriv[6] = miu_4_3_3 * const_2_C2_x;
    deriv[7] = temp * (temp_dmiu3_dy + const_2_C2_y * temp_miu3);
    deriv[8] = temp * (temp_dmiu3_dz + const_2_C2_z * temp_miu3);

    value[0] = miu_4_3_1;
    value[1] = miu_4_3_2;
    value[2] = miu_4_3_3;
}

void calc_MCSH_4_4(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    // double lambda_sqr = lambda * lambda;
    double inv_rs_2 = inv_rs * inv_rs;
    double inv_rs_4 = inv_rs_2 * inv_rs_2;
    double x0_sqr = x0*x0;
    double y0_sqr = y0*y0;
    double z0_sqr = z0*z0;

    double lambda_x0 = lambda * x0;
    double lambda_y0 = lambda * y0;
    double lambda_z0 = lambda * z0;

    double lambda_x0_sqr = lambda_x0 * lambda_x0;
    double lambda_y0_sqr = lambda_y0 * lambda_y0;
    double lambda_z0_sqr = lambda_z0 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = ((105.0 * inv_rs_4) / (2.0 * gamma)) - (15.0 * inv_rs_2);

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 105.0 * inv_rs_4 * lambda_x0_sqr + C3;
    double temp_y = 105.0 * inv_rs_4 * lambda_y0_sqr + C3;
    double temp_z = 105.0 * inv_rs_4 * lambda_z0_sqr + C3;

    double temp_dx = lambda * (210.0 * inv_rs_4 * lambda_x0); // lambda * (105.0 * 2.0 * lambda_x0);
    double temp_dy = lambda * (210.0 * inv_rs_4 * lambda_y0);
    double temp_dz = lambda * (210.0 * inv_rs_4 * lambda_z0);

    double miu_4_4_1 = temp * lambda_y0 * lambda_z0 * temp_x;
    double miu_4_4_2 = temp * lambda_x0 * lambda_z0 * temp_y;
    double miu_4_4_3 = temp * lambda_x0 * lambda_y0 * temp_z;

    //deriv consts
    double const_2_C2_x = 2.0 * C2 * x0; // 2 * C2 * x0
    double const_2_C2_y = 2.0 * C2 * y0;
    double const_2_C2_z = 2.0 * C2 * z0;

    double const_1_p_2_C2_x2 = 1.0 + 2.0 * C2 * x0_sqr; // 1 + 2*C2*x0*x0
    double const_1_p_2_C2_y2 = 1.0 + 2.0 * C2 * y0_sqr;
    double const_1_p_2_C2_z2 = 1.0 + 2.0 * C2 * z0_sqr;

    // dmiu1 dx/dy/dz
    deriv[0] = temp * lambda_y0 * lambda_z0 * (temp_dx + const_2_C2_x * temp_x);
    deriv[1] = temp * lambda_z0 * temp_x * lambda * const_1_p_2_C2_y2;
    deriv[2] = temp * lambda_y0 * temp_x * lambda * const_1_p_2_C2_z2;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * lambda_z0 * temp_y * lambda * const_1_p_2_C2_x2;
    deriv[4] = temp * lambda_x0 * lambda_z0 * (temp_dy + const_2_C2_y * temp_y);
    deriv[5] = temp * lambda_x0 * temp_y * lambda * const_1_p_2_C2_z2;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * lambda_y0 * temp_z * lambda * const_1_p_2_C2_x2;
    deriv[7] = temp * lambda_x0 * temp_z * lambda * const_1_p_2_C2_y2;
    deriv[8] = temp * lambda_x0 * lambda_y0 * (temp_dz + const_2_C2_z * temp_z);

    value[0] = miu_4_4_1;
    value[1] = miu_4_4_2;
    value[2] = miu_4_4_3;

}


void calc_MCSH_5_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double inv_rs_3 = inv_rs * inv_rs * inv_rs;
    double inv_rs_5 = inv_rs_3 * inv_rs * inv_rs;

    double lambda_x0 = lambda * x0;
    double lambda_y0 = lambda * y0;
    double lambda_z0 = lambda * z0;

    double lambda_x0_sqr = lambda_x0 * lambda_x0;
    double lambda_y0_sqr = lambda_y0 * lambda_y0;
    double lambda_z0_sqr = lambda_z0 * lambda_z0;

    double lambda_x0_3 = lambda_x0_sqr * x0 * lambda;
    double lambda_y0_3 = lambda_y0_sqr * y0 * lambda;
    double lambda_z0_3 = lambda_z0_sqr * z0 * lambda;

    double lambda_x0_4 = lambda_x0_3 * x0 * lambda;
    double lambda_y0_4 = lambda_y0_3 * y0 * lambda;
    double lambda_z0_4 = lambda_z0_3 * z0 * lambda;

    double lambda_x0_5 = lambda_x0_4 * x0 * lambda;
    double lambda_y0_5 = lambda_y0_4 * y0 * lambda;
    double lambda_z0_5 = lambda_z0_4 * z0 * lambda;

    double gamma = calc_gamma(alpha, beta);
    double C3 = ((4725.0 * inv_rs_5) / gamma) - (1050.0 * inv_rs_3);
    double C4 = ((14175.0 * inv_rs_5) / (4.0*gamma*gamma)) - ((1575.0 * inv_rs_3) / gamma) + (225.0 * inv_rs);

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 945.0 * inv_rs_5 * lambda_x0_5 + C3 * lambda_x0_3 + C4 * lambda_x0;
    double temp_y = 945.0 * inv_rs_5 * lambda_y0_5 + C3 * lambda_y0_3 + C4 * lambda_y0;
    double temp_z = 945.0 * inv_rs_5 * lambda_z0_5 + C3 * lambda_z0_3 + C4 * lambda_z0;

    double temp_dx = lambda * (4725.0 * inv_rs_5 * lambda_x0_4 + 3.0 * C3 * lambda_x0_sqr + C4);
    double temp_dy = lambda * (4725.0 * inv_rs_5 * lambda_y0_4 + 3.0 * C3 * lambda_y0_sqr + C4);
    double temp_dz = lambda * (4725.0 * inv_rs_5 * lambda_z0_4 + 3.0 * C3 * lambda_z0_sqr + C4);

    double miu_5_1_1 = temp * temp_x;
    double miu_5_1_2 = temp * temp_y;
    double miu_5_1_3 = temp * temp_z;

    //deriv consts
    double const_2_C2_x = 2.0 * C2 * x0;
    double const_2_C2_y = 2.0 * C2 * y0;
    double const_2_C2_z = 2.0 * C2 * z0;

    // dmiu1 dx/dy/dz
    deriv[0] = temp * (temp_dx + const_2_C2_x * temp_x);
    deriv[1] = miu_5_1_1 * const_2_C2_y;
    deriv[2] = miu_5_1_1 * const_2_C2_z;

    // dmiu2 dx/dy/dz
    deriv[3] = miu_5_1_2 * const_2_C2_x;
    deriv[4] = temp * (temp_dy + const_2_C2_y * temp_y);
    deriv[5] = miu_5_1_2 * const_2_C2_z;

    // dmiu3 dx/dy/dz
    deriv[6] = miu_5_1_3 * const_2_C2_x;
    deriv[7] = miu_5_1_3 * const_2_C2_y;
    deriv[8] = temp * (temp_dz + const_2_C2_z * temp_z);

    value[0] = miu_5_1_1;
    value[1] = miu_5_1_2;
    value[2] = miu_5_1_3;
}

void calc_MCSH_5_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double inv_rs_3 = inv_rs * inv_rs * inv_rs;
    double inv_rs_5 = inv_rs_3 * inv_rs * inv_rs;

    // double lambda_sqr = lambda * lambda;
    double x0_sqr = x0*x0;
    double y0_sqr = y0*y0;
    double z0_sqr = z0*z0;

    double lambda_x0 = x0 * lambda;
    double lambda_y0 = y0 * lambda;
    double lambda_z0 = z0 * lambda;

    double lambda_x0_sqr = lambda_x0 * lambda_x0;
    double lambda_y0_sqr = lambda_y0 * lambda_y0;
    double lambda_z0_sqr = lambda_z0 * lambda_z0;

    double lambda_x0_3 = lambda_x0_sqr * lambda_x0;
    double lambda_y0_3 = lambda_y0_sqr * lambda_y0;
    double lambda_z0_3 = lambda_z0_sqr * lambda_z0;

    double lambda_x0_4 = lambda_x0_3 * lambda_x0;
    double lambda_y0_4 = lambda_y0_3 * lambda_y0;
    double lambda_z0_4 = lambda_z0_3 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = ((2835.0 * inv_rs_5) / gamma) - (630.0 * inv_rs_3);
    double C4 = ((2835.0 * inv_rs_5) / (4.0 * gamma * gamma)) - ((315.0 * inv_rs_3) / gamma) + (45.0 * inv_rs);

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 945.0 * inv_rs_5 * lambda_x0_4 + C3 * lambda_x0_sqr + C4;
    double temp_y = 945.0 * inv_rs_5 * lambda_y0_4 + C3 * lambda_y0_sqr + C4;
    double temp_z = 945.0 * inv_rs_5 * lambda_z0_4 + C3 * lambda_z0_sqr + C4;

    double temp_dx = lambda * (945.0 * 4.0 * inv_rs_5 * lambda_x0_3 + 2.0 * C3 * lambda_x0);
    double temp_dy = lambda * (945.0 * 4.0 * inv_rs_5 * lambda_y0_3 + 2.0 * C3 * lambda_y0);
    double temp_dz = lambda * (945.0 * 4.0 * inv_rs_5 * lambda_z0_3 + 2.0 * C3 * lambda_z0);

    double miu_5_2_1 = temp * lambda_y0 * temp_x;
    double miu_5_2_2 = temp * lambda_x0 * temp_y;
    double miu_5_2_3 = temp * lambda_z0 * temp_x;
    double miu_5_2_4 = temp * lambda_x0 * temp_z;
    double miu_5_2_5 = temp * lambda_z0 * temp_y;
    double miu_5_2_6 = temp * lambda_y0 * temp_z;


    //deriv consts
    double const_2_C2_x = 2.0 * C2 * x0; // 2 * C2 * x0
    double const_2_C2_y = 2.0 * C2 * y0;
    double const_2_C2_z = 2.0 * C2 * z0;

    double const_1_p_2_C2_x2 = 1.0 + 2.0 * C2 * x0_sqr; // 1 + 2*C2*x0*x0
    double const_1_p_2_C2_y2 = 1.0 + 2.0 * C2 * y0_sqr;
    double const_1_p_2_C2_z2 = 1.0 + 2.0 * C2 * z0_sqr;

    // dmiu1 dx/dy/dz
    deriv[0] = temp * lambda_y0 * (temp_dx + const_2_C2_x * temp_x);
    deriv[1] = temp * lambda * temp_x * const_1_p_2_C2_y2;
    deriv[2] = miu_5_2_1 * const_2_C2_z;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * lambda * temp_y * const_1_p_2_C2_x2;
    deriv[4] = temp * lambda_x0 * (temp_dy + const_2_C2_y * temp_y);
    deriv[5] = miu_5_2_2 * const_2_C2_z;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * lambda_z0 * (temp_dx + const_2_C2_x * temp_x);
    deriv[7] = miu_5_2_3 * const_2_C2_y;
    deriv[8] = temp * lambda * temp_x * const_1_p_2_C2_z2;

    // dmiu4 dx/dy/dz
    deriv[9] = temp * lambda * temp_z * const_1_p_2_C2_x2;
    deriv[10] = miu_5_2_4 * const_2_C2_y;
    deriv[11] = temp * lambda_x0 * (temp_dz + const_2_C2_z * temp_z);

    // dmiu5 dx/dy/dz
    deriv[12] = miu_5_2_5 * const_2_C2_x;
    deriv[13] = temp * lambda_z0 * (temp_dy + const_2_C2_y * temp_y);
    deriv[14] = temp * lambda * temp_y * const_1_p_2_C2_z2;

    // dmiu6 dx/dy/dz
    deriv[15] = miu_5_2_6 * const_2_C2_x;
    deriv[16] = temp * lambda * temp_z * const_1_p_2_C2_y2;
    deriv[17] = temp * lambda_y0 * (temp_dz + const_2_C2_z * temp_z);

    value[0] = miu_5_2_1;
    value[1] = miu_5_2_2;
    value[2] = miu_5_2_3;
    value[3] = miu_5_2_4;
    value[4] = miu_5_2_5;
    value[5] = miu_5_2_6;
}

void calc_MCSH_5_3(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double lambda_sqr = lambda * lambda;
    double inv_rs_3 = inv_rs * inv_rs * inv_rs;
    double inv_rs_5 = inv_rs_3 * inv_rs * inv_rs;
    // double x0_sqr = x0*x0;
    // double y0_sqr = y0*y0;
    // double z0_sqr = z0*z0;

    double lambda_x0 = x0 * lambda;
    double lambda_y0 = y0 * lambda;
    double lambda_z0 = z0 * lambda;

    double lambda_x0_sqr = lambda_x0 * lambda_x0;
    double lambda_y0_sqr = lambda_y0 * lambda_y0;
    double lambda_z0_sqr = lambda_z0 * lambda_z0;

    double lambda_x0_3 = lambda_x0_sqr * lambda_x0;
    double lambda_y0_3 = lambda_y0_sqr * lambda_y0;
    double lambda_z0_3 = lambda_z0_sqr * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    // double C3 = (945 / (2 * gamma)) - 105;
    // double C4 = (-315 / ( 2 *gamma)) + 45;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x_2 = lambda_x0_sqr + (1.0/(2.0*gamma));
    double temp_y_2 = lambda_y0_sqr + (1.0/(2.0*gamma));
    double temp_z_2 = lambda_z0_sqr + (1.0/(2.0*gamma));

    double temp_dx_2 = 2.0 * lambda_sqr * x0;
    double temp_dy_2 = 2.0 * lambda_sqr * y0;
    double temp_dz_2 = 2.0 * lambda_sqr * z0;

    double C3_1 = 3.0 / (2.0 * gamma);
    double temp_x_3 = lambda_x0_3 + C3_1 * lambda_x0;
    double temp_y_3 = lambda_y0_3 + C3_1 * lambda_y0;
    double temp_z_3 = lambda_z0_3 + C3_1 * lambda_z0;

    double temp_dx_3 = lambda * (3.0 * lambda_x0_sqr + C3_1);
    double temp_dy_3 = lambda * (3.0 * lambda_y0_sqr + C3_1);
    double temp_dz_3 = lambda * (3.0 * lambda_z0_sqr + C3_1);

    double temp_term1_x = (945.0 * inv_rs_5 * temp_x_3) - (315.0 * inv_rs_3 * lambda_x0);
    double temp_term1_y = (945.0 * inv_rs_5 * temp_y_3) - (315.0 * inv_rs_3 * lambda_y0);
    double temp_term1_z = (945.0 * inv_rs_5 * temp_z_3) - (315.0 * inv_rs_3 * lambda_z0);

    double temp_dterm1_dx = (945.0 * inv_rs_5 * temp_dx_3) - (315.0 * inv_rs_3 * lambda);
    double temp_dterm1_dy = (945.0 * inv_rs_5 * temp_dy_3) - (315.0 * inv_rs_3 * lambda);
    double temp_dterm1_dz = (945.0 * inv_rs_5 * temp_dz_3) - (315.0 * inv_rs_3 * lambda);

    double temp_term2_x = (-105.0 * inv_rs_3 * temp_x_3) + (45.0 * inv_rs * lambda_x0);
    double temp_term2_y = (-105.0 * inv_rs_3 * temp_y_3) + (45.0 * inv_rs * lambda_y0);
    double temp_term2_z = (-105.0 * inv_rs_3 * temp_z_3) + (45.0 * inv_rs * lambda_z0);

    double temp_dterm2_dx = (-105.0 * inv_rs_3 * temp_dx_3) + (45.0 * inv_rs * lambda);
    double temp_dterm2_dy = (-105.0 * inv_rs_3 * temp_dy_3) + (45.0 * inv_rs * lambda);
    double temp_dterm2_dz = (-105.0 * inv_rs_3 * temp_dz_3) + (45.0 * inv_rs * lambda);


    // double temp_term1_x = 945.0 * temp_x_3 - 315.0 * lambda_x0;
    // double temp_term1_y = 945.0 * temp_y_3 - 315.0 * lambda_y0;
    // double temp_term1_z = 945.0 * temp_z_3 - 315.0 * lambda_z0;

    // double temp_dterm1_dx = 945.0 * temp_dx_3 - 315.0 * lambda;
    // double temp_dterm1_dy = 945.0 * temp_dy_3 - 315.0 * lambda;
    // double temp_dterm1_dz = 945.0 * temp_dz_3 - 315.0 * lambda;

    // double temp_term2_x = -105.0 * temp_x_3 + 45.0 * lambda_x0;
    // double temp_term2_y = -105.0 * temp_y_3 + 45.0 * lambda_y0;
    // double temp_term2_z = -105.0 * temp_z_3 + 45.0 * lambda_z0;

    // double temp_dterm2_dx = -105.0 * temp_dx_3 + 45.0 * lambda;
    // double temp_dterm2_dy = -105.0 * temp_dy_3 + 45.0 * lambda;
    // double temp_dterm2_dz = -105.0 * temp_dz_3 + 45.0 * lambda;

    double temp_miu1 = temp_y_2 * temp_term1_x + temp_term2_x;
    double temp_miu2 = temp_x_2 * temp_term1_y + temp_term2_y;
    double temp_miu3 = temp_z_2 * temp_term1_x + temp_term2_x;
    double temp_miu4 = temp_x_2 * temp_term1_z + temp_term2_z;
    double temp_miu5 = temp_z_2 * temp_term1_y + temp_term2_y;
    double temp_miu6 = temp_y_2 * temp_term1_z + temp_term2_z;

    double temp_dmiu1_dx = temp_y_2 * temp_dterm1_dx + temp_dterm2_dx;
    double temp_dmiu1_dy = temp_dy_2 * temp_term1_x;

    double temp_dmiu2_dx = temp_dx_2 * temp_term1_y;
    double temp_dmiu2_dy = temp_x_2 * temp_dterm1_dy + temp_dterm2_dy;

    double temp_dmiu3_dx = temp_z_2 * temp_dterm1_dx + temp_dterm2_dx;
    double temp_dmiu3_dz = temp_dz_2 * temp_term1_x;

    double temp_dmiu4_dx = temp_dx_2 * temp_term1_z;
    double temp_dmiu4_dz = temp_x_2 * temp_dterm1_dz + temp_dterm2_dz;

    double temp_dmiu5_dy = temp_z_2 * temp_dterm1_dy + temp_dterm2_dy;
    double temp_dmiu5_dz = temp_dz_2 * temp_term1_y;

    double temp_dmiu6_dy = temp_dy_2 * temp_term1_z;
    double temp_dmiu6_dz = temp_y_2 * temp_dterm1_dz + temp_dterm2_dz;

    double miu_5_3_1 = temp * temp_miu1;
    double miu_5_3_2 = temp * temp_miu2;
    double miu_5_3_3 = temp * temp_miu3;
    double miu_5_3_4 = temp * temp_miu4;
    double miu_5_3_5 = temp * temp_miu5;
    double miu_5_3_6 = temp * temp_miu6;

    //deriv consts
    double const_2_C2_x = 2.0 * C2 * x0; // 2 * C2 * x0
    double const_2_C2_y = 2.0 * C2 * y0;
    double const_2_C2_z = 2.0 * C2 * z0;

    // dmiu1 dx/dy/dz
    deriv[0] = temp * (temp_dmiu1_dx + const_2_C2_x * temp_miu1);
    deriv[1] = temp * (temp_dmiu1_dy + const_2_C2_y * temp_miu1);
    deriv[2] = miu_5_3_1 * const_2_C2_z;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * (temp_dmiu2_dx + const_2_C2_x * temp_miu2);
    deriv[4] = temp * (temp_dmiu2_dy + const_2_C2_y * temp_miu2);
    deriv[5] = miu_5_3_2 * const_2_C2_z;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * (temp_dmiu3_dx + const_2_C2_x * temp_miu3);
    deriv[7] = miu_5_3_3 * const_2_C2_y;
    deriv[8] = temp * (temp_dmiu3_dz + const_2_C2_z * temp_miu3);

    // dmiu4 dx/dy/dz
    deriv[9] = temp * (temp_dmiu4_dx + const_2_C2_x * temp_miu4);
    deriv[10] = miu_5_3_4 * const_2_C2_y;
    deriv[11] = temp * (temp_dmiu4_dz + const_2_C2_z * temp_miu4);

    // dmiu5 dx/dy/dz
    deriv[12] = miu_5_3_5 * const_2_C2_x;
    deriv[13] = temp * (temp_dmiu5_dy + const_2_C2_y * temp_miu5);
    deriv[14] = temp * (temp_dmiu5_dz + const_2_C2_z * temp_miu5);

    // dmiu6 dx/dy/dz
    deriv[15] = miu_5_3_6 * const_2_C2_x;
    deriv[16] = temp * (temp_dmiu6_dy + const_2_C2_y * temp_miu6);
    deriv[17] = temp * (temp_dmiu6_dz + const_2_C2_z * temp_miu6);

    value[0] = miu_5_3_1;
    value[1] = miu_5_3_2;
    value[2] = miu_5_3_3;
    value[3] = miu_5_3_4;
    value[4] = miu_5_3_5;
    value[5] = miu_5_3_6;
}

void calc_MCSH_5_4(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double lambda_sqr = lambda * lambda;
    double inv_rs_3 = inv_rs * inv_rs * inv_rs;
    double inv_rs_5 = inv_rs_3 * inv_rs * inv_rs;
    double x0_sqr = x0*x0;
    double y0_sqr = y0*y0;
    double z0_sqr = z0*z0;

    double lambda_x0 = x0 * lambda;
    double lambda_y0 = y0 * lambda;
    double lambda_z0 = z0 * lambda;

    double lambda_x0_sqr = lambda_sqr * x0_sqr;
    double lambda_y0_sqr = lambda_sqr * y0_sqr;
    double lambda_z0_sqr = lambda_sqr * z0_sqr;

    double lambda_x0_3 = lambda_x0_sqr * lambda_x0;
    double lambda_y0_3 = lambda_y0_sqr * lambda_y0;
    double lambda_z0_3 = lambda_z0_sqr * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = (2835.0 * inv_rs_5 / (2.0 * gamma)) - (315.0 * inv_rs_3);

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 945.0 * inv_rs_5 * lambda_x0_3 + C3 * lambda_x0;
    double temp_y = 945.0 * inv_rs_5 * lambda_y0_3 + C3 * lambda_y0;
    double temp_z = 945.0 * inv_rs_5 * lambda_z0_3 + C3 * lambda_z0;

    double temp_dx = lambda * (2835.0 * inv_rs_5 * lambda_x0_sqr + C3);
    double temp_dy = lambda * (2835.0 * inv_rs_5 * lambda_y0_sqr + C3);
    double temp_dz = lambda * (2835.0 * inv_rs_5 * lambda_z0_sqr + C3);

    double miu_5_4_1 = temp * lambda_y0 * lambda_z0 * temp_x;
    double miu_5_4_2 = temp * lambda_x0 * lambda_z0 * temp_y;
    double miu_5_4_3 = temp * lambda_x0 * lambda_y0 * temp_z;

    //deriv consts
    double const_2_C2_x = 2.0 * C2 * x0; // 2 * C2 * x0
    double const_2_C2_y = 2.0 * C2 * y0;
    double const_2_C2_z = 2.0 * C2 * z0;

    double const_1_p_2_C2_x2 = 1.0 + 2.0 * C2 * x0_sqr; // 1 + 2*C2*x0*x0
    double const_1_p_2_C2_y2 = 1.0 + 2.0 * C2 * y0_sqr;
    double const_1_p_2_C2_z2 = 1.0 + 2.0 * C2 * z0_sqr;

    // dmiu1 dx/dy/dz
    deriv[0] = temp * lambda_y0 * lambda_z0 * (temp_dx + const_2_C2_x * temp_x);
    deriv[1] = temp * lambda_z0 * temp_x * lambda * const_1_p_2_C2_y2;
    deriv[2] = temp * lambda_y0 * temp_x * lambda * const_1_p_2_C2_z2;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * lambda_z0 * temp_y * lambda * const_1_p_2_C2_x2;
    deriv[4] = temp * lambda_x0 * lambda_z0 * (temp_dy + const_2_C2_y * temp_y);
    deriv[5] = temp * lambda_x0 * temp_y * lambda * const_1_p_2_C2_z2;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * lambda_y0 * temp_z * lambda * const_1_p_2_C2_x2;
    deriv[7] = temp * lambda_x0 * temp_z * lambda * const_1_p_2_C2_y2;
    deriv[8] = temp * lambda_x0 * lambda_y0 * (temp_dz + const_2_C2_z * temp_z);

    value[0] = miu_5_4_1;
    value[1] = miu_5_4_2;
    value[2] = miu_5_4_3;

}

void calc_MCSH_5_5(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double lambda_sqr = lambda * lambda;
    double inv_rs_3 = inv_rs * inv_rs * inv_rs;
    double inv_rs_5 = inv_rs_3 * inv_rs * inv_rs;
    double x0_sqr = x0*x0;
    double y0_sqr = y0*y0;
    double z0_sqr = z0*z0;

    double lambda_x0 = x0 * lambda;
    double lambda_y0 = y0 * lambda;
    double lambda_z0 = z0 * lambda;

    double lambda_x0_sqr = lambda_x0 * lambda_x0;
    double lambda_y0_sqr = lambda_y0 * lambda_y0;
    double lambda_z0_sqr = lambda_z0 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x_2 = lambda_x0_sqr + (1.0/(2.0*gamma));
    double temp_y_2 = lambda_y0_sqr + (1.0/(2.0*gamma));
    double temp_z_2 = lambda_z0_sqr + (1.0/(2.0*gamma));

    double temp_dx_2 = 2.0 * lambda_sqr * x0;
    double temp_dy_2 = 2.0 * lambda_sqr * y0;
    double temp_dz_2 = 2.0 * lambda_sqr * z0;

    double temp_miu1 = (945.0 * inv_rs_5 * temp_x_2 * temp_y_2) - (105.0 * inv_rs_3 * (temp_x_2 + temp_y_2)) + (15.0 * inv_rs);
    double temp_miu2 = (945.0 * inv_rs_5 * temp_x_2 * temp_z_2) - (105.0 * inv_rs_3 * (temp_x_2 + temp_z_2)) + (15.0 * inv_rs);
    double temp_miu3 = (945.0 * inv_rs_5 * temp_y_2 * temp_z_2) - (105.0 * inv_rs_3 * (temp_y_2 + temp_z_2)) + (15.0 * inv_rs);

    double temp_dmiu1_dx = 945.0 * inv_rs_5 * temp_dx_2 * temp_y_2 - 105.0 * inv_rs_3 * temp_dx_2;
    double temp_dmiu1_dy = 945.0 * inv_rs_5 * temp_x_2 * temp_dy_2 - 105.0 * inv_rs_3 * temp_dy_2;

    double temp_dmiu2_dx = 945.0 * inv_rs_5 * temp_dx_2 * temp_z_2 - 105.0 * inv_rs_3 * temp_dx_2;
    double temp_dmiu2_dz = 945.0 * inv_rs_5 * temp_x_2 * temp_dz_2 - 105.0 * inv_rs_3 * temp_dz_2;

    double temp_dmiu3_dy = 945.0 * inv_rs_5 * temp_dy_2 * temp_z_2 - 105.0 * inv_rs_3 * temp_dy_2;
    double temp_dmiu3_dz = 945.0 * inv_rs_5 * temp_y_2 * temp_dz_2 - 105.0 * inv_rs_3 * temp_dz_2;

    double miu_5_5_1 = temp * lambda_z0 * temp_miu1;
    double miu_5_5_2 = temp * lambda_y0 * temp_miu2;
    double miu_5_5_3 = temp * lambda_x0 * temp_miu3;

    //deriv consts
    double const_2_C2_x = 2.0 * C2 * x0; // 2 * C2 * x0
    double const_2_C2_y = 2.0 * C2 * y0;
    double const_2_C2_z = 2.0 * C2 * z0;

    double const_1_p_2_C2_x2 = 1.0 + 2.0 * C2 * x0_sqr; // 1 + 2*C2*x0*x0
    double const_1_p_2_C2_y2 = 1.0 + 2.0 * C2 * y0_sqr;
    double const_1_p_2_C2_z2 = 1.0 + 2.0 * C2 * z0_sqr;

    // dmiu1 dx/dy/dz
    deriv[0] = temp * lambda_z0 * (temp_dmiu1_dx + const_2_C2_x * temp_miu1);
    deriv[1] = temp * lambda_z0 * (temp_dmiu1_dy + const_2_C2_y * temp_miu1);
    deriv[2] = temp * temp_miu1 * lambda * const_1_p_2_C2_z2;

    // dmiu3 dx/dy/dz
    deriv[3] = temp * lambda_y0 * (temp_dmiu2_dx + const_2_C2_x * temp_miu2);
    deriv[4] = temp * temp_miu2 * lambda * const_1_p_2_C2_y2;
    deriv[5] = temp * lambda_y0 * (temp_dmiu2_dz + const_2_C2_z * temp_miu2);

    // dmiu5 dx/dy/dz
    deriv[6] = temp * temp_miu3 * lambda * const_1_p_2_C2_x2;
    deriv[7] = temp * lambda_x0 * (temp_dmiu3_dy + const_2_C2_y * temp_miu3);
    deriv[8] = temp * lambda_x0 * (temp_dmiu3_dz + const_2_C2_z * temp_miu3);

    value[0] = miu_5_5_1;
    value[1] = miu_5_5_2;
    value[2] = miu_5_5_3;
}











void calc_MCSH_0_1_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double m_0_1 = C1 * exp( C2 * r0_sqr);

    value[0] = m_0_1;
}

void calc_MCSH_1_1_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double lambda = calc_lambda(alpha, beta);

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = lambda * x0;
    double temp_y = lambda * y0;
    double temp_z = lambda * z0;

    double miu_1_1_1 = temp * temp_x * inv_rs;
    double miu_1_1_2 = temp * temp_y * inv_rs;
    double miu_1_1_3 = temp * temp_z * inv_rs;

    value[0] = miu_1_1_1;
    value[1] = miu_1_1_2;
    value[2] = miu_1_1_3;

}

void calc_MCSH_2_1_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double inv_rs_2 = inv_rs * inv_rs;

    double lambda_x0 = x0 * lambda;
    double lambda_y0 = y0 * lambda;
    double lambda_z0 = z0 * lambda;

    double lambda_x0_sqr = lambda_x0 * lambda_x0;
    double lambda_y0_sqr = lambda_y0 * lambda_y0;
    double lambda_z0_sqr = lambda_z0 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = ((3.0 * inv_rs_2) / (2.0 * gamma)) - 1.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 3.0 * lambda_x0_sqr * inv_rs_2 + C3;
    double temp_y = 3.0 * lambda_y0_sqr * inv_rs_2 + C3;
    double temp_z = 3.0 * lambda_z0_sqr * inv_rs_2 + C3;

    double miu_2_1_1 = temp * temp_x;
    double miu_2_1_2 = temp * temp_y;
    double miu_2_1_3 = temp * temp_z;

    value[0] = miu_2_1_1;
    value[1] = miu_2_1_2;
    value[2] = miu_2_1_3;
}

void calc_MCSH_2_2_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double inv_rs_2 = inv_rs * inv_rs;

    double temp = C1 * exp( C2 * r0_sqr) * lambda * lambda * 3.0 * inv_rs_2;

    double miu_2_2_1 = temp * x0 * y0;
    double miu_2_2_2 = temp * x0 * z0;
    double miu_2_2_3 = temp * y0 * z0;

    value[0] = miu_2_2_1;
    value[1] = miu_2_2_2;
    value[2] = miu_2_2_3;
}

void calc_MCSH_3_1_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double inv_rs_3 = inv_rs * inv_rs * inv_rs;

    double lambda_x0 = x0 * lambda;
    double lambda_y0 = y0 * lambda;
    double lambda_z0 = z0 * lambda;

    double lambda_x0_sqr = lambda_x0 * lambda_x0;
    double lambda_y0_sqr = lambda_y0 * lambda_y0;
    double lambda_z0_sqr = lambda_z0 * lambda_z0;

    double lambda_x0_3 = lambda_x0_sqr * lambda_x0;
    double lambda_y0_3 = lambda_y0_sqr * lambda_y0;
    double lambda_z0_3 = lambda_z0_sqr * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = ((45.0 * inv_rs_3) / (2.0 * gamma)) - (9.0 * inv_rs);

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 15.0 * inv_rs_3 * lambda_x0_3 + C3 * lambda_x0;
    double temp_y = 15.0 * inv_rs_3  * lambda_y0_3 + C3 * lambda_y0;
    double temp_z = 15.0 * inv_rs_3  * lambda_z0_3 + C3 * lambda_z0;

    double miu_3_1_1 = temp * temp_x;
    double miu_3_1_2 = temp * temp_y;
    double miu_3_1_3 = temp * temp_z;

    value[0] = miu_3_1_1;
    value[1] = miu_3_1_2;
    value[2] = miu_3_1_3;
}

void calc_MCSH_3_2_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double inv_rs_3 = inv_rs * inv_rs * inv_rs;

    double lambda_x0 = x0 * lambda;
    double lambda_y0 = y0 * lambda;
    double lambda_z0 = z0 * lambda;

    double lambda_x0_sqr = lambda_x0 * lambda_x0;
    double lambda_y0_sqr = lambda_y0 * lambda_y0;
    double lambda_z0_sqr = lambda_z0 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = ((15.0 * inv_rs_3) / (2.0 * gamma)) - (3.0 * inv_rs);

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 15.0 * inv_rs_3 * lambda_x0_sqr + C3;
    double temp_y = 15.0 * inv_rs_3 * lambda_y0_sqr + C3;
    double temp_z = 15.0 * inv_rs_3 * lambda_z0_sqr + C3;

    double miu_3_2_1 = temp * lambda_y0 * temp_x;
    double miu_3_2_2 = temp * lambda_x0 * temp_y;
    double miu_3_2_3 = temp * lambda_z0 * temp_x;
    double miu_3_2_4 = temp * lambda_x0 * temp_z;
    double miu_3_2_5 = temp * lambda_z0 * temp_y;
    double miu_3_2_6 = temp * lambda_y0 * temp_z;

    value[0] = miu_3_2_1;
    value[1] = miu_3_2_2;
    value[2] = miu_3_2_3;
    value[3] = miu_3_2_4;
    value[4] = miu_3_2_5;
    value[5] = miu_3_2_6;
}

void calc_MCSH_3_3_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double lambda = calc_lambda(alpha, beta);
    double inv_rs_3 = inv_rs * inv_rs * inv_rs;

    double temp =  C1 * exp( C2 * r0_sqr) * lambda * lambda * lambda * 15.0 * inv_rs_3;
    double m_3_3 = temp * x0 * y0 * z0;

    value[0] = m_3_3;
}

void calc_MCSH_4_1_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double inv_rs_2 = inv_rs * inv_rs;
    double inv_rs_4 = inv_rs_2 * inv_rs_2;

    double lambda_x0 = lambda * x0;
    double lambda_y0 = lambda * y0;
    double lambda_z0 = lambda * z0;

    double lambda_x0_sqr = lambda_x0 * lambda_x0;
    double lambda_y0_sqr = lambda_y0 * lambda_y0;
    double lambda_z0_sqr = lambda_z0 * lambda_z0;

    double lambda_x0_4 = lambda_x0_sqr * lambda_x0_sqr;
    double lambda_y0_4 = lambda_y0_sqr * lambda_y0_sqr;
    double lambda_z0_4 = lambda_z0_sqr * lambda_z0_sqr;

    double gamma = calc_gamma(alpha, beta);
    double C3 = ((315.0 * inv_rs_4) / gamma) - (90.0 * inv_rs_2);
    double C4 = ((315.0 * inv_rs_4) / (4.0*gamma*gamma)) - (45.0 * inv_rs_2 / gamma) + 9.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 105.0 * inv_rs_4 * lambda_x0_4 + C3 * lambda_x0_sqr + C4;
    double temp_y = 105.0 * inv_rs_4 * lambda_y0_4 + C3 * lambda_y0_sqr + C4;
    double temp_z = 105.0 * inv_rs_4 * lambda_z0_4 + C3 * lambda_z0_sqr + C4;

    double miu_4_1_1 = temp * temp_x;
    double miu_4_1_2 = temp * temp_y;
    double miu_4_1_3 = temp * temp_z;

    value[0] = miu_4_1_1;
    value[1] = miu_4_1_2;
    value[2] = miu_4_1_3;
}

void calc_MCSH_4_2_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double inv_rs_2 = inv_rs * inv_rs;
    double inv_rs_4 = inv_rs_2 * inv_rs_2;

    double lambda_x0 = lambda * x0;
    double lambda_y0 = lambda * y0;
    double lambda_z0 = lambda * z0;

    double lambda_x0_sqr = lambda_x0 * lambda_x0;
    double lambda_y0_sqr = lambda_y0 * lambda_y0;
    double lambda_z0_sqr = lambda_z0 * lambda_z0;

    double lambda_x0_3 = lambda_x0_sqr * lambda_x0;
    double lambda_y0_3 = lambda_y0_sqr * lambda_y0;
    double lambda_z0_3 = lambda_z0_sqr * lambda_z0;

    double gamma = calc_gamma(alpha, beta);

    double C3 = ((315.0 * inv_rs_4) / (2.0 * gamma)) - (45.0 * inv_rs_2);

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 105.0 * inv_rs_4 * lambda_x0_3 + C3 * lambda_x0;
    double temp_y = 105.0 * inv_rs_4 * lambda_y0_3 + C3 * lambda_y0;
    double temp_z = 105.0 * inv_rs_4 * lambda_z0_3 + C3 * lambda_z0;

    double miu_4_2_1 = temp * lambda_y0 * temp_x;
    double miu_4_2_2 = temp * lambda_x0 * temp_y;
    double miu_4_2_3 = temp * lambda_z0 * temp_x;
    double miu_4_2_4 = temp * lambda_x0 * temp_z;
    double miu_4_2_5 = temp * lambda_z0 * temp_y;
    double miu_4_2_6 = temp * lambda_y0 * temp_z;

    value[0] = miu_4_2_1;
    value[1] = miu_4_2_2;
    value[2] = miu_4_2_3;
    value[3] = miu_4_2_4;
    value[4] = miu_4_2_5;
    value[5] = miu_4_2_6;
}


void calc_MCSH_4_3_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double lambda_sqr = lambda * lambda;
    double inv_rs_2 = inv_rs * inv_rs;
    double inv_rs_4 = inv_rs_2 * inv_rs_2;
    double x0_sqr = x0*x0;
    double y0_sqr = y0*y0;
    double z0_sqr = z0*z0;

    double lambda_x0_sqr = lambda_sqr * x0_sqr;
    double lambda_y0_sqr = lambda_sqr * y0_sqr;
    double lambda_z0_sqr = lambda_sqr * z0_sqr;

    double gamma = calc_gamma(alpha, beta);

    double temp_x_2 = lambda_x0_sqr + (1.0/(2.0*gamma));
    double temp_y_2 = lambda_y0_sqr + (1.0/(2.0*gamma));
    double temp_z_2 = lambda_z0_sqr + (1.0/(2.0*gamma));

    double temp_term1_x = 105.0 * inv_rs_4 * temp_x_2 - (15.0 * inv_rs_2);
    double temp_term1_y = 105.0 * inv_rs_4 * temp_y_2 - (15.0 * inv_rs_2);
    // double temp_term1_z = 105.0 * temp_z_2 - 15.0;

    double temp_term2_x = -15.0 * inv_rs_2 * temp_x_2 + 3.0;
    double temp_term2_y = -15.0 * inv_rs_2 * temp_y_2 + 3.0;
    // double temp_term2_z = -15.0 * temp_z_2 + 3.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_miu1 = temp_y_2 * temp_term1_x + temp_term2_x;
    double temp_miu2 = temp_z_2 * temp_term1_x + temp_term2_x;
    double temp_miu3 = temp_z_2 * temp_term1_y + temp_term2_y;

    double miu_4_3_1 = temp * temp_miu1;
    double miu_4_3_2 = temp * temp_miu2;
    double miu_4_3_3 = temp * temp_miu3;

    value[0] = miu_4_3_1;
    value[1] = miu_4_3_2;
    value[2] = miu_4_3_3;
}

void calc_MCSH_4_4_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double inv_rs_2 = inv_rs * inv_rs;
    double inv_rs_4 = inv_rs_2 * inv_rs_2;

    double lambda_x0 = lambda * x0;
    double lambda_y0 = lambda * y0;
    double lambda_z0 = lambda * z0;

    double lambda_x0_sqr = lambda_x0 * lambda_x0;
    double lambda_y0_sqr = lambda_y0 * lambda_y0;
    double lambda_z0_sqr = lambda_z0 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = ((105.0 * inv_rs_4) / (2.0 * gamma)) - (15.0 * inv_rs_2);

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 105.0 * inv_rs_4 * lambda_x0_sqr + C3;
    double temp_y = 105.0 * inv_rs_4 * lambda_y0_sqr + C3;
    double temp_z = 105.0 * inv_rs_4 * lambda_z0_sqr + C3;

    double miu_4_4_1 = temp * lambda_y0 * lambda_z0 * temp_x;
    double miu_4_4_2 = temp * lambda_x0 * lambda_z0 * temp_y;
    double miu_4_4_3 = temp * lambda_x0 * lambda_y0 * temp_z;

    value[0] = miu_4_4_1;
    value[1] = miu_4_4_2;
    value[2] = miu_4_4_3;

}


void calc_MCSH_5_1_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double inv_rs_3 = inv_rs * inv_rs * inv_rs;
    double inv_rs_5 = inv_rs_3 * inv_rs * inv_rs;

    double lambda_x0 = lambda * x0;
    double lambda_y0 = lambda * y0;
    double lambda_z0 = lambda * z0;

    double lambda_x0_sqr = lambda_x0 * lambda_x0;
    double lambda_y0_sqr = lambda_y0 * lambda_y0;
    double lambda_z0_sqr = lambda_z0 * lambda_z0;

    double lambda_x0_3 = lambda_x0_sqr * lambda_x0;
    double lambda_y0_3 = lambda_y0_sqr * lambda_y0;
    double lambda_z0_3 = lambda_z0_sqr * lambda_z0;

    double lambda_x0_4 = lambda_x0_3 * lambda_x0;
    double lambda_y0_4 = lambda_y0_3 * lambda_y0;
    double lambda_z0_4 = lambda_z0_3 * lambda_z0;

    double lambda_x0_5 = lambda_x0_4 * lambda_x0;
    double lambda_y0_5 = lambda_y0_4 * lambda_y0;
    double lambda_z0_5 = lambda_z0_4 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = ((4725.0 * inv_rs_5) / gamma) - (1050.0 * inv_rs_3);
    double C4 = ((14175.0 * inv_rs_5) / (4.0*gamma*gamma)) - ((1575.0 * inv_rs_3) / gamma) + (225.0 * inv_rs);

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 945.0 * inv_rs_5 * lambda_x0_5 + C3 * lambda_x0_3 + C4 * lambda_x0;
    double temp_y = 945.0 * inv_rs_5 * lambda_y0_5 + C3 * lambda_y0_3 + C4 * lambda_y0;
    double temp_z = 945.0 * inv_rs_5 * lambda_z0_5 + C3 * lambda_z0_3 + C4 * lambda_z0;

    double miu_5_1_1 = temp * temp_x;
    double miu_5_1_2 = temp * temp_y;
    double miu_5_1_3 = temp * temp_z;

    value[0] = miu_5_1_1;
    value[1] = miu_5_1_2;
    value[2] = miu_5_1_3;
}

void calc_MCSH_5_2_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double inv_rs_3 = inv_rs * inv_rs * inv_rs;
    double inv_rs_5 = inv_rs_3 * inv_rs * inv_rs;

    double lambda_x0 = x0 * lambda;
    double lambda_y0 = y0 * lambda;
    double lambda_z0 = z0 * lambda;

    double lambda_x0_sqr = lambda_x0 * lambda_x0;
    double lambda_y0_sqr = lambda_y0 * lambda_y0;
    double lambda_z0_sqr = lambda_z0 * lambda_z0;

    double lambda_x0_3 = lambda_x0_sqr * lambda_x0;
    double lambda_y0_3 = lambda_y0_sqr * lambda_y0;
    double lambda_z0_3 = lambda_z0_sqr * lambda_z0;

    double lambda_x0_4 = lambda_x0_3 * lambda_x0;
    double lambda_y0_4 = lambda_y0_3 * lambda_y0;
    double lambda_z0_4 = lambda_z0_3 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = ((2835.0 * inv_rs_5) / gamma) - (630.0 * inv_rs_3);
    double C4 = ((2835.0 * inv_rs_5) / (4.0 * gamma * gamma)) - ((315.0 * inv_rs_3) / gamma) + (45.0 * inv_rs);

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 945.0 * inv_rs_5 * lambda_x0_4 + C3 * lambda_x0_sqr + C4;
    double temp_y = 945.0 * inv_rs_5 * lambda_y0_4 + C3 * lambda_y0_sqr + C4;
    double temp_z = 945.0 * inv_rs_5 * lambda_z0_4 + C3 * lambda_z0_sqr + C4;

    double miu_5_2_1 = temp * lambda_y0 * temp_x;
    double miu_5_2_2 = temp * lambda_x0 * temp_y;
    double miu_5_2_3 = temp * lambda_z0 * temp_x;
    double miu_5_2_4 = temp * lambda_x0 * temp_z;
    double miu_5_2_5 = temp * lambda_z0 * temp_y;
    double miu_5_2_6 = temp * lambda_y0 * temp_z;

    value[0] = miu_5_2_1;
    value[1] = miu_5_2_2;
    value[2] = miu_5_2_3;
    value[3] = miu_5_2_4;
    value[4] = miu_5_2_5;
    value[5] = miu_5_2_6;
}

void calc_MCSH_5_3_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double inv_rs_3 = inv_rs * inv_rs * inv_rs;
    double inv_rs_5 = inv_rs_3 * inv_rs * inv_rs;

    double lambda_x0 = x0 * lambda;
    double lambda_y0 = y0 * lambda;
    double lambda_z0 = z0 * lambda;

    double lambda_x0_sqr = lambda_x0 * lambda_x0;
    double lambda_y0_sqr = lambda_y0 * lambda_y0;
    double lambda_z0_sqr = lambda_z0 * lambda_z0;

    double lambda_x0_3 = lambda_x0_sqr * lambda_x0;
    double lambda_y0_3 = lambda_y0_sqr * lambda_y0;
    double lambda_z0_3 = lambda_z0_sqr * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    // double C3 = (945 / (2 * gamma)) - 105;
    // double C4 = (-315 / ( 2 *gamma)) + 45;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x_2 = lambda_x0_sqr + (1.0/(2.0*gamma));
    double temp_y_2 = lambda_y0_sqr + (1.0/(2.0*gamma));
    double temp_z_2 = lambda_z0_sqr + (1.0/(2.0*gamma));

    double C3_1 = 3.0 / (2.0 * gamma);
    double temp_x_3 = lambda_x0_3 + C3_1 * lambda_x0;
    double temp_y_3 = lambda_y0_3 + C3_1 * lambda_y0;
    double temp_z_3 = lambda_z0_3 + C3_1 * lambda_z0;

    double temp_term1_x = (945.0 * inv_rs_5 * temp_x_3) - (315.0 * inv_rs_3 * lambda_x0);
    double temp_term1_y = (945.0 * inv_rs_5 * temp_y_3) - (315.0 * inv_rs_3 * lambda_y0);
    double temp_term1_z = (945.0 * inv_rs_5 * temp_z_3) - (315.0 * inv_rs_3 * lambda_z0);

    double temp_term2_x = (-105.0 * inv_rs_3 * temp_x_3) + (45.0 * inv_rs * lambda_x0);
    double temp_term2_y = (-105.0 * inv_rs_3 * temp_y_3) + (45.0 * inv_rs * lambda_y0);
    double temp_term2_z = (-105.0 * inv_rs_3 * temp_z_3) + (45.0 * inv_rs * lambda_z0);

    double temp_miu1 = temp_y_2 * temp_term1_x + temp_term2_x;
    double temp_miu2 = temp_x_2 * temp_term1_y + temp_term2_y;
    double temp_miu3 = temp_z_2 * temp_term1_x + temp_term2_x;
    double temp_miu4 = temp_x_2 * temp_term1_z + temp_term2_z;
    double temp_miu5 = temp_z_2 * temp_term1_y + temp_term2_y;
    double temp_miu6 = temp_y_2 * temp_term1_z + temp_term2_z;

    double miu_5_3_1 = temp * temp_miu1;
    double miu_5_3_2 = temp * temp_miu2;
    double miu_5_3_3 = temp * temp_miu3;
    double miu_5_3_4 = temp * temp_miu4;
    double miu_5_3_5 = temp * temp_miu5;
    double miu_5_3_6 = temp * temp_miu6;

    value[0] = miu_5_3_1;
    value[1] = miu_5_3_2;
    value[2] = miu_5_3_3;
    value[3] = miu_5_3_4;
    value[4] = miu_5_3_5;
    value[5] = miu_5_3_6;
}

void calc_MCSH_5_4_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double inv_rs_3 = inv_rs * inv_rs * inv_rs;
    double inv_rs_5 = inv_rs_3 * inv_rs * inv_rs;

    double lambda_x0 = x0 * lambda;
    double lambda_y0 = y0 * lambda;
    double lambda_z0 = z0 * lambda;

    double lambda_x0_sqr = lambda_x0 * lambda_x0;
    double lambda_y0_sqr = lambda_y0 * lambda_y0;
    double lambda_z0_sqr = lambda_y0 * lambda_z0;

    double lambda_x0_3 = lambda_x0_sqr * lambda_x0;
    double lambda_y0_3 = lambda_y0_sqr * lambda_y0;
    double lambda_z0_3 = lambda_z0_sqr * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = (2835.0 * inv_rs_5 / (2.0 * gamma)) - (315.0 * inv_rs_3);

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 945.0 * inv_rs_5 * lambda_x0_3 + C3 * lambda_x0;
    double temp_y = 945.0 * inv_rs_5 * lambda_y0_3 + C3 * lambda_y0;
    double temp_z = 945.0 * inv_rs_5 * lambda_z0_3 + C3 * lambda_z0;

    double miu_5_4_1 = temp * lambda_y0 * lambda_z0 * temp_x;
    double miu_5_4_2 = temp * lambda_x0 * lambda_z0 * temp_y;
    double miu_5_4_3 = temp * lambda_x0 * lambda_y0 * temp_z;

    value[0] = miu_5_4_1;
    value[1] = miu_5_4_2;
    value[2] = miu_5_4_3;

}

void calc_MCSH_5_5_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double inv_rs_3 = inv_rs * inv_rs * inv_rs;
    double inv_rs_5 = inv_rs_3 * inv_rs * inv_rs;

    double lambda_x0 = x0 * lambda;
    double lambda_y0 = y0 * lambda;
    double lambda_z0 = z0 * lambda;

    double lambda_x0_sqr = lambda_x0 * lambda_x0;
    double lambda_y0_sqr = lambda_y0 * lambda_y0;
    double lambda_z0_sqr = lambda_z0 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x_2 = lambda_x0_sqr + (1.0/(2.0*gamma));
    double temp_y_2 = lambda_y0_sqr + (1.0/(2.0*gamma));
    double temp_z_2 = lambda_z0_sqr + (1.0/(2.0*gamma));

    double temp_miu1 = (945.0 * inv_rs_5 * temp_x_2 * temp_y_2) - (105.0 * inv_rs_3 * (temp_x_2 + temp_y_2)) + (15.0 * inv_rs);
    double temp_miu2 = (945.0 * inv_rs_5 * temp_x_2 * temp_z_2) - (105.0 * inv_rs_3 * (temp_x_2 + temp_z_2)) + (15.0 * inv_rs);
    double temp_miu3 = (945.0 * inv_rs_5 * temp_y_2 * temp_z_2) - (105.0 * inv_rs_3 * (temp_y_2 + temp_z_2)) + (15.0 * inv_rs);

    double miu_5_5_1 = temp * lambda_z0 * temp_miu1;
    double miu_5_5_2 = temp * lambda_y0 * temp_miu2;
    double miu_5_5_3 = temp * lambda_x0 * temp_miu3;

    value[0] = miu_5_5_1;
    value[1] = miu_5_5_2;
    value[2] = miu_5_5_3;
}




int get_mcsh_type(int mcsh_order, int group_num)
{
    if (mcsh_order == 0) {
        if (group_num == 1) {
            return 1;
        } else {
            return 0;
        }
    } else if (mcsh_order == 1) {
        if (group_num == 1) {
            return 2;
        } else {
            return 0;
        }
    } else if (mcsh_order == 2) {
        if (group_num == 1) {
            return 2;
        } else if (group_num == 2){
            return 2;
        } else {
            return 0;
        }
    } else if (mcsh_order == 3) {
        if (group_num == 1) {
            return 2;
        } else if (group_num == 2){
            return 3;
        } else if (group_num == 3){
            return 1;
        } else {
            return 0;
        }
    } else if (mcsh_order == 4) {
        if (group_num == 1) {
            return 2;
        } else if (group_num == 2){
            return 3;
        } else if (group_num == 3){
            return 2;
        } else if (group_num == 4){
            return 2;
        } else {
            return 0;
        }
    } else if (mcsh_order == 5) {
        if (group_num == 1) {
            return 2;
        } else if (group_num == 2){
            return 3;
        } else if (group_num == 3){
            return 3;
        } else if (group_num == 4){
            return 2;
        } else if (group_num == 5){
            return 2;
        } else {
            return 0;
        }
    // } else if (mcsh_order == 6) {
    //     if (group_num == 1) {
    //         return 2;
    //     } else if (group_num == 2){
    //         return 3;
    //     } else if (group_num == 3){
    //         return 3;
    //     } else if (group_num == 4){
    //         return 2;
    //     } else if (group_num == 5){
    //         return 2;
    //     } else if (group_num == 6){
    //         return 3;
    //     } else if (group_num == 7){
    //         return 1;
    //     } else {
    //         return 0;
    //     }
    // } else if (mcsh_order == 7) {
    //     if (group_num == 1) {
    //         return 2;
    //     } else if (group_num == 2){
    //         return 3;
    //     } else if (group_num == 3){
    //         return 3;
    //     } else if (group_num == 4){
    //         return 2;
    //     } else if (group_num == 5){
    //         return 3;
    //     } else if (group_num == 6){
    //         return 3;
    //     } else if (group_num == 7){
    //         return 2;
    //     } else if (group_num == 8){
    //         return 2;
    //     } else {
    //         return 0;
    //     }
    // } else if (mcsh_order == 8) {
    //     if (group_num == 1) {
    //         return 2;
    //     } else if (group_num == 2){
    //         return 3;
    //     } else if (group_num == 3){
    //         return 3;
    //     } else if (group_num == 4){
    //         return 2;
    //     } else if (group_num == 5){
    //         return 3;
    //     } else if (group_num == 6){
    //         return 3;
    //     } else if (group_num == 7){
    //         return 2;
    //     } else if (group_num == 8){
    //         return 3;
    //     } else if (group_num == 9){
    //         return 2;
    //     } else if (group_num == 10){
    //         return 2;
    //     } else {
    //         return 0;
    //     }
    // } else if (mcsh_order == 9) {
    //     if (group_num == 1) {
    //         return 2;
    //     } else if (group_num == 2){
    //         return 3;
    //     } else if (group_num == 3){
    //         return 3;
    //     } else if (group_num == 4){
    //         return 2;
    //     } else if (group_num == 5){
    //         return 3;
    //     } else if (group_num == 6){
    //         return 3;
    //     } else if (group_num == 7){
    //         return 3;
    //     } else if (group_num == 8){
    //         return 3;
    //     } else if (group_num == 9){
    //         return 2;
    //     } else if (group_num == 10){
    //         return 2;
    //     } else if (group_num == 11){
    //         return 3;
    //     } else if (group_num == 12){
    //         return 1;
    //     } else {
    //         return 0;
    //     }
    } else {
        return 0;
    }
}


GMPFunction get_mcsh_function(int mcsh_order, int group_num)
{
    GMPFunction result;

    if (mcsh_order == 0) {
        if (group_num == 1) {
            result = calc_MCSH_0_1;
        }
    } else if (mcsh_order == 1) {
        if (group_num == 1) {
            result = calc_MCSH_1_1;
        }
    } else if (mcsh_order == 2) {
        if (group_num == 1) {
            result = calc_MCSH_2_1;
        } else if (group_num == 2){
            result = calc_MCSH_2_2;
        }
    } else if (mcsh_order == 3) {
        if (group_num == 1) {
            result = calc_MCSH_3_1;
        } else if (group_num == 2){
            result = calc_MCSH_3_2;
        } else if (group_num == 3){
            result = calc_MCSH_3_3;
        }
    } else if (mcsh_order == 4) {
        if (group_num == 1) {
            result = calc_MCSH_4_1;
        } else if (group_num == 2){
            result = calc_MCSH_4_2;
        } else if (group_num == 3){
            result = calc_MCSH_4_3;
        } else if (group_num == 4){
            result = calc_MCSH_4_4;
        }
    } else if (mcsh_order == 5) {
        if (group_num == 1) {
            result = calc_MCSH_5_1;
        } else if (group_num == 2){
            result = calc_MCSH_5_2;
        } else if (group_num == 3){
            result = calc_MCSH_5_3;
        } else if (group_num == 4){
            result = calc_MCSH_5_4;
        } else if (group_num == 5){
            result = calc_MCSH_5_5;
        }
    } // else if (mcsh_order == 6) {
    //     if (group_num == 1) {
    //         result = calc_MCSH_6_1;
    //     } else if (group_num == 2){
    //         result = calc_MCSH_6_2;
    //     } else if (group_num == 3){
    //         result = calc_MCSH_6_3;
    //     } else if (group_num == 4){
    //         result = calc_MCSH_6_4;
    //     } else if (group_num == 5){
    //         result = calc_MCSH_6_5;
    //     } else if (group_num == 6){
    //         result = calc_MCSH_6_6;
    //     } else if (group_num == 7){
    //         result = calc_MCSH_6_7;
    //     }
    // } else if (mcsh_order == 7) {
    //     if (group_num == 1) {
    //         result = calc_MCSH_7_1;
    //     } else if (group_num == 2){
    //         result = calc_MCSH_7_2;
    //     } else if (group_num == 3){
    //         result = calc_MCSH_7_3;
    //     } else if (group_num == 4){
    //         result = calc_MCSH_7_4;
    //     } else if (group_num == 5){
    //         result = calc_MCSH_7_5;
    //     } else if (group_num == 6){
    //         result = calc_MCSH_7_6;
    //     } else if (group_num == 7){
    //         result = calc_MCSH_7_7;
    //     } else if (group_num == 8){
    //         result = calc_MCSH_7_8;
    //     }
    // } else if (mcsh_order == 8) {
    //     if (group_num == 1) {
    //         result = calc_MCSH_8_1;
    //     } else if (group_num == 2){
    //         result = calc_MCSH_8_2;
    //     } else if (group_num == 3){
    //         result = calc_MCSH_8_3;
    //     } else if (group_num == 4){
    //         result = calc_MCSH_8_4;
    //     } else if (group_num == 5){
    //         result = calc_MCSH_8_5;
    //     } else if (group_num == 6){
    //         result = calc_MCSH_8_6;
    //     } else if (group_num == 7){
    //         result = calc_MCSH_8_7;
    //     } else if (group_num == 8){
    //         result = calc_MCSH_8_8;
    //     } else if (group_num == 9){
    //         result = calc_MCSH_8_9;
    //     } else if (group_num == 10){
    //         result = calc_MCSH_8_10;
    //     }
    // } else if (mcsh_order == 9) {
    //     if (group_num == 1) {
    //         result = calc_MCSH_9_1;
    //     } else if (group_num == 2){
    //         result = calc_MCSH_9_2;
    //     } else if (group_num == 3){
    //         result = calc_MCSH_9_3;
    //     } else if (group_num == 4){
    //         result = calc_MCSH_9_4;
    //     } else if (group_num == 5){
    //         result = calc_MCSH_9_5;
    //     } else if (group_num == 6){
    //         result = calc_MCSH_9_6;
    //     } else if (group_num == 7){
    //         result = calc_MCSH_9_7;
    //     } else if (group_num == 8){
    //         result = calc_MCSH_9_8;
    //     } else if (group_num == 9){
    //         result = calc_MCSH_9_9;
    //     } else if (group_num == 10){
    //         result = calc_MCSH_9_10;
    //     } else if (group_num == 11){
    //         result = calc_MCSH_9_11;
    //     } else if (group_num == 12){
    //         result = calc_MCSH_9_12;
    //     }
    // }

    return result;
}


GMPFunctionNoderiv get_mcsh_function_noderiv(int mcsh_order, int group_num)
{
    GMPFunctionNoderiv result;

    if (mcsh_order == 0) {
        if (group_num == 1) {
            result = calc_MCSH_0_1_noderiv;
        }
    } else if (mcsh_order == 1) {
        if (group_num == 1) {
            result = calc_MCSH_1_1_noderiv;
        }
    } else if (mcsh_order == 2) {
        if (group_num == 1) {
            result = calc_MCSH_2_1_noderiv;
        } else if (group_num == 2){
            result = calc_MCSH_2_2_noderiv;
        }
    } else if (mcsh_order == 3) {
        if (group_num == 1) {
            result = calc_MCSH_3_1_noderiv;
        } else if (group_num == 2){
            result = calc_MCSH_3_2_noderiv;
        } else if (group_num == 3){
            result = calc_MCSH_3_3_noderiv;
        }
    } else if (mcsh_order == 4) {
        if (group_num == 1) {
            result = calc_MCSH_4_1_noderiv;
        } else if (group_num == 2){
            result = calc_MCSH_4_2_noderiv;
        } else if (group_num == 3){
            result = calc_MCSH_4_3_noderiv;
        } else if (group_num == 4){
            result = calc_MCSH_4_4_noderiv;
        }
    } else if (mcsh_order == 5) {
        if (group_num == 1) {
            result = calc_MCSH_5_1_noderiv;
        } else if (group_num == 2){
            result = calc_MCSH_5_2_noderiv;
        } else if (group_num == 3){
            result = calc_MCSH_5_3_noderiv;
        } else if (group_num == 4){
            result = calc_MCSH_5_4_noderiv;
        } else if (group_num == 5){
            result = calc_MCSH_5_5_noderiv;
        }
    } // else if (mcsh_order == 6) {
    //     if (group_num == 1) {
    //         result = calc_MCSH_6_1_noderiv;
    //     } else if (group_num == 2){
    //         result = calc_MCSH_6_2_noderiv;
    //     } else if (group_num == 3){
    //         result = calc_MCSH_6_3_noderiv;
    //     } else if (group_num == 4){
    //         result = calc_MCSH_6_4_noderiv;
    //     } else if (group_num == 5){
    //         result = calc_MCSH_6_5_noderiv;
    //     } else if (group_num == 6){
    //         result = calc_MCSH_6_6_noderiv;
    //     } else if (group_num == 7){
    //         result = calc_MCSH_6_7_noderiv;
    //     }
    // } else if (mcsh_order == 7) {
    //     if (group_num == 1) {
    //         result = calc_MCSH_7_1_noderiv;
    //     } else if (group_num == 2){
    //         result = calc_MCSH_7_2_noderiv;
    //     } else if (group_num == 3){
    //         result = calc_MCSH_7_3_noderiv;
    //     } else if (group_num == 4){
    //         result = calc_MCSH_7_4_noderiv;
    //     } else if (group_num == 5){
    //         result = calc_MCSH_7_5_noderiv;
    //     } else if (group_num == 6){
    //         result = calc_MCSH_7_6_noderiv;
    //     } else if (group_num == 7){
    //         result = calc_MCSH_7_7_noderiv;
    //     } else if (group_num == 8){
    //         result = calc_MCSH_7_8_noderiv;
    //     }
    // } else if (mcsh_order == 8) {
    //     if (group_num == 1) {
    //         result = calc_MCSH_8_1_noderiv;
    //     } else if (group_num == 2){
    //         result = calc_MCSH_8_2_noderiv;
    //     } else if (group_num == 3){
    //         result = calc_MCSH_8_3_noderiv;
    //     } else if (group_num == 4){
    //         result = calc_MCSH_8_4_noderiv;
    //     } else if (group_num == 5){
    //         result = calc_MCSH_8_5_noderiv;
    //     } else if (group_num == 6){
    //         result = calc_MCSH_8_6_noderiv;
    //     } else if (group_num == 7){
    //         result = calc_MCSH_8_7_noderiv;
    //     } else if (group_num == 8){
    //         result = calc_MCSH_8_8_noderiv;
    //     } else if (group_num == 9){
    //         result = calc_MCSH_8_9_noderiv;
    //     } else if (group_num == 10){
    //         result = calc_MCSH_8_10_noderiv;
    //     }
    // } else if (mcsh_order == 9) {
    //     if (group_num == 1) {
    //         result = calc_MCSH_9_1_noderiv;
    //     } else if (group_num == 2){
    //         result = calc_MCSH_9_2_noderiv;
    //     } else if (group_num == 3){
    //         result = calc_MCSH_9_3_noderiv;
    //     } else if (group_num == 4){
    //         result = calc_MCSH_9_4_noderiv;
    //     } else if (group_num == 5){
    //         result = calc_MCSH_9_5_noderiv;
    //     } else if (group_num == 6){
    //         result = calc_MCSH_9_6_noderiv;
    //     } else if (group_num == 7){
    //         result = calc_MCSH_9_7_noderiv;
    //     } else if (group_num == 8){
    //         result = calc_MCSH_9_8_noderiv;
    //     } else if (group_num == 9){
    //         result = calc_MCSH_9_9_noderiv;
    //     } else if (group_num == 10){
    //         result = calc_MCSH_9_10_noderiv;
    //     } else if (group_num == 11){
    //         result = calc_MCSH_9_11_noderiv;
    //     } else if (group_num == 12){
    //         result = calc_MCSH_9_12_noderiv;
    //     }
    // }

    return result;
}