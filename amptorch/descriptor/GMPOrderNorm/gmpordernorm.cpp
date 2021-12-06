#include <math.h>
#include <stdio.h>
#include "gmpordernorm.h"

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


double P1(double lambda, double x0, double gamma){
    return lambda * x0;
}

double P2(double lambda, double x0, double gamma){
    double lambda_x0_2 = lambda * lambda * x0 * x0;
    return (0.5 / gamma) + lambda_x0_2;
}

double P3(double lambda, double x0, double gamma){
    double lambda_x0 = lambda * x0;
    double lambda_x0_3 = lambda_x0 * lambda_x0 * lambda_x0;
    return (1.5 * lambda_x0 / gamma) + lambda_x0_3;
}

double P4(double lambda, double x0, double gamma){
    double lambda_x0_2 = lambda * lambda * x0 * x0;
    double lambda_x0_4 = lambda_x0_2 * lambda_x0_2;
    return (0.75 / (gamma * gamma)) + (3.0 * lambda_x0_2 / gamma) + lambda_x0_4;
}

double P5(double lambda, double x0, double gamma){
    double lambda_x0 = lambda * x0;
    double lambda_x0_2 = lambda_x0 * lambda_x0;
    double lambda_x0_3 = lambda_x0 * lambda_x0_2;
    double lambda_x0_5 = lambda_x0_3 * lambda_x0_2;
    return ((15.0 * lambda_x0) / (4.0 * gamma * gamma)) + (5.0 * lambda_x0_3 / gamma) + lambda_x0_5;
}

double P6(double lambda, double x0, double gamma){
    double lambda_x0_2 = lambda * lambda * x0 * x0;
    double lambda_x0_4 = lambda_x0_2 * lambda_x0_2;
    double lambda_x0_6 = lambda_x0_4 * lambda_x0_2;
    return (15.0 / (8.0 * gamma * gamma * gamma)) + ((11.25 * lambda_x0_2) / (gamma * gamma)) + (7.5 * lambda_x0_4 / gamma) + lambda_x0_6;
}

double P7(double lambda, double x0, double gamma){
    double lambda_x0 = lambda * x0;
    double lambda_x0_2 = lambda_x0 * lambda_x0;
    double lambda_x0_3 = lambda_x0 * lambda_x0_2;
    double lambda_x0_5 = lambda_x0_3 * lambda_x0_2;
    double lambda_x0_7 = lambda_x0_5 * lambda_x0_2;
    double term1 = 105.0 * lambda_x0 / (8.0 * gamma * gamma * gamma);
    double term2 = 105.0 * lambda_x0_3 / (4.0 * gamma * gamma);
    double term3 = 10.5 * lambda_x0_5 / gamma;
    return term1 + term2 + term3 + lambda_x0_7;
}

double P8(double lambda, double x0, double gamma){
    double lambda_x0_2 = lambda * lambda * x0 * x0;
    double lambda_x0_4 = lambda_x0_2 * lambda_x0_2;
    double lambda_x0_6 = lambda_x0_4 * lambda_x0_2;
    double lambda_x0_8 = lambda_x0_6 * lambda_x0_2;
    double term1 = 105.0 / (16.0 * gamma * gamma * gamma * gamma);
    double term2 = 105.0 * lambda_x0_2 / (2.0 * gamma * gamma * gamma) ;
    double term3 = 105.0 * lambda_x0_4 / (2.0 * gamma * gamma);
    double term4 = 14.0 * lambda_x0_6 / gamma;
    return term1 + term2 + term3 + term4 + lambda_x0_8;
}

double P9(double lambda, double x0, double gamma){
    double lambda_x0 = lambda * x0;
    double lambda_x0_2 = lambda_x0 * lambda_x0;
    double lambda_x0_3 = lambda_x0 * lambda_x0_2;
    double lambda_x0_5 = lambda_x0_3 * lambda_x0_2;
    double lambda_x0_7 = lambda_x0_5 * lambda_x0_2;
    double lambda_x0_9 = lambda_x0_7 * lambda_x0_2;
    double term1 = 945.0 * lambda_x0 / (16.0 * gamma * gamma * gamma * gamma);
    double term2 = 315.0 * lambda_x0_3 / (2.0 * gamma * gamma * gamma);
    double term3 = 189.0 * lambda_x0_5 / (2.0 * gamma * gamma);
    double term4 = 18.0 * lambda_x0_7 / gamma;
    return term1 + term2 + term3 + term4 + lambda_x0_9;
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
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);
    double inv_rs_3 = inv_rs * inv_rs * inv_rs;
    double inv_rs_5 = inv_rs_3 * inv_rs * inv_rs;

    double P1x = P1(lambda, x0, gamma);
    double P1y = P1(lambda, y0, gamma);
    double P1z = P1(lambda, z0, gamma);

    double P3x = P3(lambda, x0, gamma);
    double P3y = P3(lambda, y0, gamma);
    double P3z = P3(lambda, z0, gamma);

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double term_x = (945.0 * inv_rs_5 * P5x) - (1050.0 * inv_rs_3 * P3x) + (225.0 * inv_rs * P1x);
    double term_y = (945.0 * inv_rs_5 * P5y) - (1050.0 * inv_rs_3 * P3y) + (225.0 * inv_rs * P1y);
    double term_z = (945.0 * inv_rs_5 * P5z) - (1050.0 * inv_rs_3 * P3z) + (225.0 * inv_rs * P1z);

    double miu_5_1_1 = temp * term_x;
    double miu_5_1_2 = temp * term_y;
    double miu_5_1_3 = temp * term_z;

    value[0] = miu_5_1_1;
    value[1] = miu_5_1_2;
    value[2] = miu_5_1_3;
}

void calc_MCSH_5_2_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);
    double inv_rs_3 = inv_rs * inv_rs * inv_rs;
    double inv_rs_5 = inv_rs_3 * inv_rs * inv_rs;

    double P1x = P1(lambda, x0, gamma);
    double P1y = P1(lambda, y0, gamma);
    double P1z = P1(lambda, z0, gamma);

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double P4x = P4(lambda, x0, gamma);
    double P4y = P4(lambda, y0, gamma);
    double P4z = P4(lambda, z0, gamma);

    double term_x = (945.0 * inv_rs_5 * P4x) - (630.0 * inv_rs_3 * P2x) + (45.0 * inv_rs);
    double term_y = (945.0 * inv_rs_5 * P4y) - (630.0 * inv_rs_3 * P2y) + (45.0 * inv_rs);
    double term_z = (945.0 * inv_rs_5 * P4z) - (630.0 * inv_rs_3 * P2z) + (45.0 * inv_rs);

    double miu_5_2_1 = temp * P1y * term_x;
    double miu_5_2_2 = temp * P1x * term_y;
    double miu_5_2_3 = temp * P1z * term_x;
    double miu_5_2_4 = temp * P1x * term_z;
    double miu_5_2_5 = temp * P1z * term_y;
    double miu_5_2_6 = temp * P1y * term_z;

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
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);
    double inv_rs_3 = inv_rs * inv_rs * inv_rs;
    double inv_rs_5 = inv_rs_3 * inv_rs * inv_rs;

    double P1x = P1(lambda, x0, gamma);
    double P1y = P1(lambda, y0, gamma);
    double P1z = P1(lambda, z0, gamma);

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double P3x = P3(lambda, x0, gamma);
    double P3y = P3(lambda, y0, gamma);
    double P3z = P3(lambda, z0, gamma);

    double temp_miu1 = (945.0 * inv_rs_5 * P3x * P2y) - (105.0 * inv_rs_3 * P3x) - (315.0 * inv_rs_3 * P1x * P2y) + (45.0 * inv_rs * P1x);
    double temp_miu2 = (945.0 * inv_rs_5 * P3y * P2x) - (105.0 * inv_rs_3 * P3y) - (315.0 * inv_rs_3 * P1y * P2x) + (45.0 * inv_rs * P1y);
    double temp_miu3 = (945.0 * inv_rs_5 * P3x * P2z) - (105.0 * inv_rs_3 * P3x) - (315.0 * inv_rs_3 * P1x * P2z) + (45.0 * inv_rs * P1x);
    double temp_miu4 = (945.0 * inv_rs_5 * P3z * P2x) - (105.0 * inv_rs_3 * P3z) - (315.0 * inv_rs_3 * P1z * P2x) + (45.0 * inv_rs * P1z);
    double temp_miu5 = (945.0 * inv_rs_5 * P3y * P2z) - (105.0 * inv_rs_3 * P3y) - (315.0 * inv_rs_3 * P1y * P2z) + (45.0 * inv_rs * P1y);
    double temp_miu6 = (945.0 * inv_rs_5 * P3z * P2y) - (105.0 * inv_rs_3 * P3z) - (315.0 * inv_rs_3 * P1z * P2y) + (45.0 * inv_rs * P1z);

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
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);
    double inv_rs_3 = inv_rs * inv_rs * inv_rs;
    double inv_rs_5 = inv_rs_3 * inv_rs * inv_rs;

    double P1x = P1(lambda, x0, gamma);
    double P1y = P1(lambda, y0, gamma);
    double P1z = P1(lambda, z0, gamma);

    double P3x = P3(lambda, x0, gamma);
    double P3y = P3(lambda, y0, gamma);
    double P3z = P3(lambda, z0, gamma);

    double temp_x = (945.0 * inv_rs_5 * P3x) - (315.0 * inv_rs_3 * P1x);
    double temp_y = (945.0 * inv_rs_5 * P3y) - (315.0 * inv_rs_3 * P1y);
    double temp_z = (945.0 * inv_rs_5 * P3z) - (315.0 * inv_rs_3 * P1z);

    double miu_5_4_1 = temp * P1y * P1z * temp_x;
    double miu_5_4_2 = temp * P1x * P1z * temp_y;
    double miu_5_4_3 = temp * P1x * P1y * temp_z;

    value[0] = miu_5_4_1;
    value[1] = miu_5_4_2;
    value[2] = miu_5_4_3;
}

void calc_MCSH_5_5_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);
    double inv_rs_3 = inv_rs * inv_rs * inv_rs;
    double inv_rs_5 = inv_rs_3 * inv_rs * inv_rs;

    double P1x = P1(lambda, x0, gamma);
    double P1y = P1(lambda, y0, gamma);
    double P1z = P1(lambda, z0, gamma);

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double temp_miu1 = (945.0 * inv_rs_5 * P2x * P2y) - (105.0 * inv_rs_3 * (P2x + P2y)) + (15.0 * inv_rs);
    double temp_miu2 = (945.0 * inv_rs_5 * P2x * P2z) - (105.0 * inv_rs_3 * (P2x + P2z)) + (15.0 * inv_rs);
    double temp_miu3 = (945.0 * inv_rs_5 * P2y * P2z) - (105.0 * inv_rs_3 * (P2y + P2z)) + (15.0 * inv_rs);

    double miu_5_5_1 = temp * P1z * temp_miu1;
    double miu_5_5_2 = temp * P1y * temp_miu2;
    double miu_5_5_3 = temp * P1x * temp_miu3;

    value[0] = miu_5_5_1;
    value[1] = miu_5_5_2;
    value[2] = miu_5_5_3;
}




void calc_MCSH_6_1_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);
    double inv_rs_2 = inv_rs * inv_rs;
    double inv_rs_4 = inv_rs_2 * inv_rs_2;
    double inv_rs_6 = inv_rs_4 * inv_rs_2;

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double P4x = P4(lambda, x0, gamma);
    double P4y = P4(lambda, y0, gamma);
    double P4z = P4(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double term_x = (10395.0 * inv_rs_6 * P6x) - (14175.0 * inv_rs_4 * P4x) + (4725.0 * inv_rs_2 * P2x) - 225.0;
    double term_y = (10395.0 * inv_rs_6 * P6y) - (14175.0 * inv_rs_4 * P4y) + (4725.0 * inv_rs_2 * P2y) - 225.0;
    double term_z = (10395.0 * inv_rs_6 * P6z) - (14175.0 * inv_rs_4 * P4z) + (4725.0 * inv_rs_2 * P2z) - 225.0;

    double miu_1 = temp * term_x;
    double miu_2 = temp * term_y;
    double miu_3 = temp * term_z;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
}

void calc_MCSH_6_2_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);
    double inv_rs_2 = inv_rs * inv_rs;
    double inv_rs_4 = inv_rs_2 * inv_rs_2;
    double inv_rs_6 = inv_rs_4 * inv_rs_2;

    double P1x = P1(lambda, x0, gamma);
    double P1y = P1(lambda, y0, gamma);
    double P1z = P1(lambda, z0, gamma);

    double P3x = P3(lambda, x0, gamma);
    double P3y = P3(lambda, y0, gamma);
    double P3z = P3(lambda, z0, gamma);

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double term_x = (10395.0 * inv_rs_6 * P5x) - (9450.0 * inv_rs_4 * P3x) + (1575.0 * inv_rs_2 * P1x);
    double term_y = (10395.0 * inv_rs_6 * P5y) - (9450.0 * inv_rs_4 * P3y) + (1575.0 * inv_rs_2 * P1y);
    double term_z = (10395.0 * inv_rs_6 * P5z) - (9450.0 * inv_rs_4 * P3z) + (1575.0 * inv_rs_2 * P1z);

    double miu_1 = temp * P1y * term_x;
    double miu_2 = temp * P1x * term_y;
    double miu_3 = temp * P1z * term_x;
    double miu_4 = temp * P1x * term_z;
    double miu_5 = temp * P1z * term_y;
    double miu_6 = temp * P1y * term_z;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;
}


void calc_MCSH_6_3_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);
    double inv_rs_2 = inv_rs * inv_rs;
    double inv_rs_4 = inv_rs_2 * inv_rs_2;
    double inv_rs_6 = inv_rs_4 * inv_rs_2;

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double P4x = P4(lambda, x0, gamma);
    double P4y = P4(lambda, y0, gamma);
    double P4z = P4(lambda, z0, gamma);

    double temp_miu1 = (10395.0 * inv_rs_6 * P4x * P2y) - (945.0 * inv_rs_4 * P4x) - (5670.0 * inv_rs_4 * P2x * P2y) + (630.0 * inv_rs_2 * P2x) + (315.0 * inv_rs_2 * P2y) - 45.0;
    double temp_miu2 = (10395.0 * inv_rs_6 * P4y * P2x) - (945.0 * inv_rs_4 * P4y) - (5670.0 * inv_rs_4 * P2y * P2x) + (630.0 * inv_rs_2 * P2y) + (315.0 * inv_rs_2 * P2x) - 45.0;
    double temp_miu3 = (10395.0 * inv_rs_6 * P4x * P2z) - (945.0 * inv_rs_4 * P4x) - (5670.0 * inv_rs_4 * P2x * P2z) + (630.0 * inv_rs_2 * P2x) + (315.0 * inv_rs_2 * P2z) - 45.0;
    double temp_miu4 = (10395.0 * inv_rs_6 * P4z * P2x) - (945.0 * inv_rs_4 * P4z) - (5670.0 * inv_rs_4 * P2z * P2x) + (630.0 * inv_rs_2 * P2z) + (315.0 * inv_rs_2 * P2x) - 45.0;
    double temp_miu5 = (10395.0 * inv_rs_6 * P4y * P2z) - (945.0 * inv_rs_4 * P4y) - (5670.0 * inv_rs_4 * P2y * P2z) + (630.0 * inv_rs_2 * P2y) + (315.0 * inv_rs_2 * P2z) - 45.0;
    double temp_miu6 = (10395.0 * inv_rs_6 * P4z * P2y) - (945.0 * inv_rs_4 * P4z) - (5670.0 * inv_rs_4 * P2z * P2y) + (630.0 * inv_rs_2 * P2z) + (315.0 * inv_rs_2 * P2y) - 45.0;

    double miu_1 = temp * temp_miu1;
    double miu_2 = temp * temp_miu2;
    double miu_3 = temp * temp_miu3;
    double miu_4 = temp * temp_miu4;
    double miu_5 = temp * temp_miu5;
    double miu_6 = temp * temp_miu6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;
}

void calc_MCSH_6_4_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);
    double inv_rs_2 = inv_rs * inv_rs;
    double inv_rs_4 = inv_rs_2 * inv_rs_2;
    double inv_rs_6 = inv_rs_4 * inv_rs_2;

    double P1x = P1(lambda, x0, gamma);
    double P1y = P1(lambda, y0, gamma);
    double P1z = P1(lambda, z0, gamma);

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double P4x = P4(lambda, x0, gamma);
    double P4y = P4(lambda, y0, gamma);
    double P4z = P4(lambda, z0, gamma);

    double temp_x = (10395.0 * inv_rs_6 * P4x) - (5670.0 * inv_rs_4 * P2x) + (315.0 * inv_rs_2);
    double temp_y = (10395.0 * inv_rs_6 * P4y) - (5670.0 * inv_rs_4 * P2y) + (315.0 * inv_rs_2);
    double temp_z = (10395.0 * inv_rs_6 * P4z) - (5670.0 * inv_rs_4 * P2z) + (315.0 * inv_rs_2);

    double miu_1 = temp * P1y * P1z * temp_x;
    double miu_2 = temp * P1x * P1z * temp_y;
    double miu_3 = temp * P1x * P1y * temp_z;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
}


void calc_MCSH_6_5_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);
    double inv_rs_2 = inv_rs * inv_rs;
    double inv_rs_4 = inv_rs_2 * inv_rs_2;
    double inv_rs_6 = inv_rs_4 * inv_rs_2;

    double P1x = P1(lambda, x0, gamma);
    double P1y = P1(lambda, y0, gamma);
    double P1z = P1(lambda, z0, gamma);

    double P3x = P3(lambda, x0, gamma);
    double P3y = P3(lambda, y0, gamma);
    double P3z = P3(lambda, z0, gamma);

    double temp_1 = (10395.0 * inv_rs_6 * P3x * P3y) - (2835.0 * inv_rs_4 * ((P3x * P1y) + (P1x * P3y))) + (945.0 * inv_rs_2 * P1x * P1y);
    double temp_2 = (10395.0 * inv_rs_6 * P3x * P3z) - (2835.0 * inv_rs_4 * ((P3x * P1z) + (P1x * P3z))) + (945.0 * inv_rs_2 * P1x * P1z);
    double temp_3 = (10395.0 * inv_rs_6 * P3y * P3z) - (2835.0 * inv_rs_4 * ((P3y * P1z) + (P1y * P3z))) + (945.0 * inv_rs_2 * P1y * P1z);

    double miu_1 = temp * temp_1;
    double miu_2 = temp * temp_2;
    double miu_3 = temp * temp_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
}

void calc_MCSH_6_6_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);
    double inv_rs_2 = inv_rs * inv_rs;
    double inv_rs_4 = inv_rs_2 * inv_rs_2;
    double inv_rs_6 = inv_rs_4 * inv_rs_2;

    double P1x = P1(lambda, x0, gamma);
    double P1y = P1(lambda, y0, gamma);
    double P1z = P1(lambda, z0, gamma);

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double P3x = P3(lambda, x0, gamma);
    double P3y = P3(lambda, y0, gamma);
    double P3z = P3(lambda, z0, gamma);

    double term_1 = (10395.0 * inv_rs_6 * P3x * P2y) - (945.0 * inv_rs_4 * P3x) - (2835.0 * inv_rs_4 * P1x * P2y) + (315.0 * inv_rs_2 * P1x);
    double term_2 = (10395.0 * inv_rs_6 * P3y * P2x) - (945.0 * inv_rs_4 * P3y) - (2835.0 * inv_rs_4 * P1y * P2x) + (315.0 * inv_rs_2 * P1y);
    double term_3 = (10395.0 * inv_rs_6 * P3x * P2z) - (945.0 * inv_rs_4 * P3x) - (2835.0 * inv_rs_4 * P1x * P2z) + (315.0 * inv_rs_2 * P1x);
    double term_4 = (10395.0 * inv_rs_6 * P3z * P2x) - (945.0 * inv_rs_4 * P3z) - (2835.0 * inv_rs_4 * P1z * P2x) + (315.0 * inv_rs_2 * P1z);
    double term_5 = (10395.0 * inv_rs_6 * P3y * P2z) - (945.0 * inv_rs_4 * P3y) - (2835.0 * inv_rs_4 * P1y * P2z) + (315.0 * inv_rs_2 * P1y);
    double term_6 = (10395.0 * inv_rs_6 * P3z * P2y) - (945.0 * inv_rs_4 * P3z) - (2835.0 * inv_rs_4 * P1z * P2y) + (315.0 * inv_rs_2 * P1z);

    double miu_1 = temp * P1z * term_1;
    double miu_2 = temp * P1z * term_2;
    double miu_3 = temp * P1y * term_3;
    double miu_4 = temp * P1y * term_4;
    double miu_5 = temp * P1x * term_5;
    double miu_6 = temp * P1x * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;
}

void calc_MCSH_6_7_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);
    double inv_rs_2 = inv_rs * inv_rs;
    double inv_rs_4 = inv_rs_2 * inv_rs_2;
    double inv_rs_6 = inv_rs_4 * inv_rs_2;


    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double temp2 = (10395.0 * inv_rs_6 * P2x * P2y * P2z) - (945.0 * inv_rs_4 * (P2x * P2y + P2x * P2z + P2y * P2z)) + (105.0 * inv_rs_2 * (P2x + P2y + P2z)) - 15.0;

    double m = temp * temp2;

    value[0] = m;
}







void calc_MCSH_7_1_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);
    double inv_rs_2 = inv_rs * inv_rs;
    double inv_rs_3 = inv_rs * inv_rs_2;
    double inv_rs_5 = inv_rs_3 * inv_rs_2;
    double inv_rs_7 = inv_rs_5 * inv_rs_2;

    double P1x = P1(lambda, x0, gamma);
    double P1y = P1(lambda, y0, gamma);
    double P1z = P1(lambda, z0, gamma);

    double P3x = P3(lambda, x0, gamma);
    double P3y = P3(lambda, y0, gamma);
    double P3z = P3(lambda, z0, gamma);

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double term_x = (135135.0 * inv_rs_7 * P7x) - (218295.0 * inv_rs_5 * P5x) + (99225.0 * inv_rs_3 * P3x) - (11025.0 * inv_rs * P1x);
    double term_y = (135135.0 * inv_rs_7 * P7y) - (218295.0 * inv_rs_5 * P5y) + (99225.0 * inv_rs_3 * P3y) - (11025.0 * inv_rs * P1y);
    double term_z = (135135.0 * inv_rs_7 * P7z) - (218295.0 * inv_rs_5 * P5z) + (99225.0 * inv_rs_3 * P3z) - (11025.0 * inv_rs * P1z);

    double miu_1 = temp * term_x;
    double miu_2 = temp * term_y;
    double miu_3 = temp * term_z;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
}

void calc_MCSH_7_2_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);
    double inv_rs_2 = inv_rs * inv_rs;
    double inv_rs_3 = inv_rs * inv_rs_2;
    double inv_rs_5 = inv_rs_3 * inv_rs_2;
    double inv_rs_7 = inv_rs_5 * inv_rs_2;

    double P1x = P1(lambda, x0, gamma);
    double P1y = P1(lambda, y0, gamma);
    double P1z = P1(lambda, z0, gamma);

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double P4x = P4(lambda, x0, gamma);
    double P4y = P4(lambda, y0, gamma);
    double P4z = P4(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double term_x = (135135.0 * inv_rs_7 * P6x) - (155925.0 * inv_rs_5 * P4x) + (42525.0 * inv_rs_3 * P2x) - (1575.0 * inv_rs);
    double term_y = (135135.0 * inv_rs_7 * P6y) - (155925.0 * inv_rs_5 * P4y) + (42525.0 * inv_rs_3 * P2y) - (1575.0 * inv_rs);
    double term_z = (135135.0 * inv_rs_7 * P6z) - (155925.0 * inv_rs_5 * P4z) + (42525.0 * inv_rs_3 * P2z) - (1575.0 * inv_rs);

    double miu_1 = temp * P1y * term_x;
    double miu_2 = temp * P1x * term_y;
    double miu_3 = temp * P1z * term_x;
    double miu_4 = temp * P1x * term_z;
    double miu_5 = temp * P1z * term_y;
    double miu_6 = temp * P1y * term_z;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;
}


void calc_MCSH_7_3_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);
    double inv_rs_2 = inv_rs * inv_rs;
    double inv_rs_3 = inv_rs * inv_rs_2;
    double inv_rs_5 = inv_rs_3 * inv_rs_2;
    double inv_rs_7 = inv_rs_5 * inv_rs_2;

    double P1x = P1(lambda, x0, gamma);
    double P1y = P1(lambda, y0, gamma);
    double P1z = P1(lambda, z0, gamma);

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double P3x = P3(lambda, x0, gamma);
    double P3y = P3(lambda, y0, gamma);
    double P3z = P3(lambda, z0, gamma);

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double temp_miu1 = (135135.0 * inv_rs_7 * P5x * P2y) - (10395.0 * inv_rs_5 * P5x) - (103950.0 * inv_rs_5 * P3x * P2y) + (9450.0 * inv_rs_3 * P3x) + (14175.0 * inv_rs_3 * P1x * P2y) - (1575.0 * inv_rs * P1x);
    double temp_miu2 = (135135.0 * inv_rs_7 * P5y * P2x) - (10395.0 * inv_rs_5 * P5y) - (103950.0 * inv_rs_5 * P3y * P2x) + (9450.0 * inv_rs_3 * P3y) + (14175.0 * inv_rs_3 * P1y * P2x) - (1575.0 * inv_rs * P1y);
    double temp_miu3 = (135135.0 * inv_rs_7 * P5x * P2z) - (10395.0 * inv_rs_5 * P5x) - (103950.0 * inv_rs_5 * P3x * P2z) + (9450.0 * inv_rs_3 * P3x) + (14175.0 * inv_rs_3 * P1x * P2z) - (1575.0 * inv_rs * P1x);
    double temp_miu4 = (135135.0 * inv_rs_7 * P5z * P2x) - (10395.0 * inv_rs_5 * P5z) - (103950.0 * inv_rs_5 * P3z * P2x) + (9450.0 * inv_rs_3 * P3z) + (14175.0 * inv_rs_3 * P1z * P2x) - (1575.0 * inv_rs * P1z);
    double temp_miu5 = (135135.0 * inv_rs_7 * P5y * P2z) - (10395.0 * inv_rs_5 * P5y) - (103950.0 * inv_rs_5 * P3y * P2z) + (9450.0 * inv_rs_3 * P3y) + (14175.0 * inv_rs_3 * P1y * P2z) - (1575.0 * inv_rs * P1y);
    double temp_miu6 = (135135.0 * inv_rs_7 * P5z * P2y) - (10395.0 * inv_rs_5 * P5z) - (103950.0 * inv_rs_5 * P3z * P2y) + (9450.0 * inv_rs_3 * P3z) + (14175.0 * inv_rs_3 * P1z * P2y) - (1575.0 * inv_rs * P1z);

    double miu_1 = temp * temp_miu1;
    double miu_2 = temp * temp_miu2;
    double miu_3 = temp * temp_miu3;
    double miu_4 = temp * temp_miu4;
    double miu_5 = temp * temp_miu5;
    double miu_6 = temp * temp_miu6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;
}

void calc_MCSH_7_4_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);
    double inv_rs_2 = inv_rs * inv_rs;
    double inv_rs_3 = inv_rs * inv_rs_2;
    double inv_rs_5 = inv_rs_3 * inv_rs_2;
    double inv_rs_7 = inv_rs_5 * inv_rs_2;

    double P1x = P1(lambda, x0, gamma);
    double P1y = P1(lambda, y0, gamma);
    double P1z = P1(lambda, z0, gamma);

    double P3x = P3(lambda, x0, gamma);
    double P3y = P3(lambda, y0, gamma);
    double P3z = P3(lambda, z0, gamma);

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double temp_x = (135135.0 * inv_rs_7 * P5x) - (103950.0 * inv_rs_5 * P3x) + (14175.0 * inv_rs_3 * P1x);
    double temp_y = (135135.0 * inv_rs_7 * P5y) - (103950.0 * inv_rs_5 * P3y) + (14175.0 * inv_rs_3 * P1y);
    double temp_z = (135135.0 * inv_rs_7 * P5z) - (103950.0 * inv_rs_5 * P3z) + (14175.0 * inv_rs_3 * P1z);

    double miu_1 = temp * P1y * P1z * temp_x;
    double miu_2 = temp * P1x * P1z * temp_y;
    double miu_3 = temp * P1x * P1y * temp_z;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
}

void calc_MCSH_7_5_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);
    double inv_rs_2 = inv_rs * inv_rs;
    double inv_rs_3 = inv_rs * inv_rs_2;
    double inv_rs_5 = inv_rs_3 * inv_rs_2;
    double inv_rs_7 = inv_rs_5 * inv_rs_2;

    double P1x = P1(lambda, x0, gamma);
    double P1y = P1(lambda, y0, gamma);
    double P1z = P1(lambda, z0, gamma);

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double P3x = P3(lambda, x0, gamma);
    double P3y = P3(lambda, y0, gamma);
    double P3z = P3(lambda, z0, gamma);

    double P4x = P4(lambda, x0, gamma);
    double P4y = P4(lambda, y0, gamma);
    double P4z = P4(lambda, z0, gamma);

    double temp_miu1 = (135135.0 * inv_rs_7 * P4x * P3y) - (31185.0 * inv_rs_5 * P4x * P1y) - (62370.0 * inv_rs_5 * P2x * P3y) + (17010.0 * inv_rs_3 * P2x * P1y) + (2835.0 * inv_rs_3 * P3y) - (945.0 * inv_rs * P1y);
    double temp_miu2 = (135135.0 * inv_rs_7 * P4y * P3x) - (31185.0 * inv_rs_5 * P4y * P1x) - (62370.0 * inv_rs_5 * P2y * P3x) + (17010.0 * inv_rs_3 * P2y * P1x) + (2835.0 * inv_rs_3 * P3x) - (945.0 * inv_rs * P1x);
    double temp_miu3 = (135135.0 * inv_rs_7 * P4x * P3z) - (31185.0 * inv_rs_5 * P4x * P1z) - (62370.0 * inv_rs_5 * P2x * P3z) + (17010.0 * inv_rs_3 * P2x * P1z) + (2835.0 * inv_rs_3 * P3z) - (945.0 * inv_rs * P1z);
    double temp_miu4 = (135135.0 * inv_rs_7 * P4z * P3x) - (31185.0 * inv_rs_5 * P4z * P1x) - (62370.0 * inv_rs_5 * P2z * P3x) + (17010.0 * inv_rs_3 * P2z * P1x) + (2835.0 * inv_rs_3 * P3x) - (945.0 * inv_rs * P1x);
    double temp_miu5 = (135135.0 * inv_rs_7 * P4y * P3z) - (31185.0 * inv_rs_5 * P4y * P1z) - (62370.0 * inv_rs_5 * P2y * P3z) + (17010.0 * inv_rs_3 * P2y * P1z) + (2835.0 * inv_rs_3 * P3z) - (945.0 * inv_rs * P1z);
    double temp_miu6 = (135135.0 * inv_rs_7 * P4z * P3y) - (31185.0 * inv_rs_5 * P4z * P1y) - (62370.0 * inv_rs_5 * P2z * P3y) + (17010.0 * inv_rs_3 * P2z * P1y) + (2835.0 * inv_rs_3 * P3y) - (945.0 * inv_rs * P1y);

    double miu_1 = temp * temp_miu1;
    double miu_2 = temp * temp_miu2;
    double miu_3 = temp * temp_miu3;
    double miu_4 = temp * temp_miu4;
    double miu_5 = temp * temp_miu5;
    double miu_6 = temp * temp_miu6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;
}

void calc_MCSH_7_6_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);
    double inv_rs_2 = inv_rs * inv_rs;
    double inv_rs_3 = inv_rs * inv_rs_2;
    double inv_rs_5 = inv_rs_3 * inv_rs_2;
    double inv_rs_7 = inv_rs_5 * inv_rs_2;

    double P1x = P1(lambda, x0, gamma);
    double P1y = P1(lambda, y0, gamma);
    double P1z = P1(lambda, z0, gamma);

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double P4x = P4(lambda, x0, gamma);
    double P4y = P4(lambda, y0, gamma);
    double P4z = P4(lambda, z0, gamma);

    double term_1 = (135135.0 * inv_rs_7 * P4x * P2y) - (10395.0 * inv_rs_5 * P4x) - (62370.0 * inv_rs_5 * P2x * P2y) + (5670.0 * inv_rs_3 * P2x) + (2835.0 * inv_rs_3 * P2y) - (315.0 * inv_rs);
    double term_2 = (135135.0 * inv_rs_7 * P4y * P2x) - (10395.0 * inv_rs_5 * P4y) - (62370.0 * inv_rs_5 * P2y * P2x) + (5670.0 * inv_rs_3 * P2y) + (2835.0 * inv_rs_3 * P2x) - (315.0 * inv_rs);
    double term_3 = (135135.0 * inv_rs_7 * P4x * P2z) - (10395.0 * inv_rs_5 * P4x) - (62370.0 * inv_rs_5 * P2x * P2z) + (5670.0 * inv_rs_3 * P2x) + (2835.0 * inv_rs_3 * P2z) - (315.0 * inv_rs);
    double term_4 = (135135.0 * inv_rs_7 * P4z * P2x) - (10395.0 * inv_rs_5 * P4z) - (62370.0 * inv_rs_5 * P2z * P2x) + (5670.0 * inv_rs_3 * P2z) + (2835.0 * inv_rs_3 * P2x) - (315.0 * inv_rs);
    double term_5 = (135135.0 * inv_rs_7 * P4y * P2z) - (10395.0 * inv_rs_5 * P4y) - (62370.0 * inv_rs_5 * P2y * P2z) + (5670.0 * inv_rs_3 * P2y) + (2835.0 * inv_rs_3 * P2z) - (315.0 * inv_rs);
    double term_6 = (135135.0 * inv_rs_7 * P4z * P2y) - (10395.0 * inv_rs_5 * P4z) - (62370.0 * inv_rs_5 * P2z * P2y) + (5670.0 * inv_rs_3 * P2z) + (2835.0 * inv_rs_3 * P2y) - (315.0 * inv_rs);

    double miu_1 = temp * P1z * term_1;
    double miu_2 = temp * P1z * term_2;
    double miu_3 = temp * P1y * term_3;
    double miu_4 = temp * P1y * term_4;
    double miu_5 = temp * P1x * term_5;
    double miu_6 = temp * P1x * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;
}

void calc_MCSH_7_7_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);
    double inv_rs_2 = inv_rs * inv_rs;
    double inv_rs_3 = inv_rs * inv_rs_2;
    double inv_rs_5 = inv_rs_3 * inv_rs_2;
    double inv_rs_7 = inv_rs_5 * inv_rs_2;

    double P1x = P1(lambda, x0, gamma);
    double P1y = P1(lambda, y0, gamma);
    double P1z = P1(lambda, z0, gamma);

    double P3x = P3(lambda, x0, gamma);
    double P3y = P3(lambda, y0, gamma);
    double P3z = P3(lambda, z0, gamma);

    double temp_1 = (135135.0 * inv_rs_7 * P3x * P3y) - (31185.0 * inv_rs_5 * ((P3x * P1y) + (P1x * P3y))) + (8505.0 * inv_rs_3 * P1x * P1y);
    double temp_2 = (135135.0 * inv_rs_7 * P3x * P3z) - (31185.0 * inv_rs_5 * ((P3x * P1z) + (P1x * P3z))) + (8505.0 * inv_rs_3 * P1x * P1z);
    double temp_3 = (135135.0 * inv_rs_7 * P3y * P3z) - (31185.0 * inv_rs_5 * ((P3y * P1z) + (P1y * P3z))) + (8505.0 * inv_rs_3 * P1y * P1z);

    double miu_1 = temp * P1z * temp_1;
    double miu_2 = temp * P1y * temp_2;
    double miu_3 = temp * P1x * temp_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
}


void calc_MCSH_7_8_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);
    double inv_rs_2 = inv_rs * inv_rs;
    double inv_rs_3 = inv_rs * inv_rs_2;
    double inv_rs_5 = inv_rs_3 * inv_rs_2;
    double inv_rs_7 = inv_rs_5 * inv_rs_2;

    double P1x = P1(lambda, x0, gamma);
    double P1y = P1(lambda, y0, gamma);
    double P1z = P1(lambda, z0, gamma);

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double P3x = P3(lambda, x0, gamma);
    double P3y = P3(lambda, y0, gamma);
    double P3z = P3(lambda, z0, gamma);

    double temp_1 = (135135.0 * inv_rs_7 * P3x * P2y * P2z) - (10395.0 * inv_rs_5 * P3x * (P2y + P2z)) + (945.0 * inv_rs_3 * P3x) - (31185.0 * inv_rs_5 * P1x * P2y * P2z) + (2835.0 * inv_rs_3 * P1x * (P2y + P2z)) - (315.0 * inv_rs * P1x);
    double temp_2 = (135135.0 * inv_rs_7 * P2x * P3y * P2z) - (10395.0 * inv_rs_5 * P3y * (P2x + P2z)) + (945.0 * inv_rs_3 * P3y) - (31185.0 * inv_rs_5 * P2x * P1y * P2z) + (2835.0 * inv_rs_3 * P1y * (P2x + P2z)) - (315.0 * inv_rs * P1y);
    double temp_3 = (135135.0 * inv_rs_7 * P2x * P2y * P3z) - (10395.0 * inv_rs_5 * P3z * (P2x + P2y)) + (945.0 * inv_rs_3 * P3z) - (31185.0 * inv_rs_5 * P2x * P2y * P1z) + (2835.0 * inv_rs_3 * P1z * (P2x + P2y)) - (315.0 * inv_rs * P1z);

    double miu_1 = temp * temp_1;
    double miu_2 = temp * temp_2;
    double miu_3 = temp * temp_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
}


void calc_MCSH_8_1_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);
    double inv_rs_2 = inv_rs * inv_rs;
    double inv_rs_4 = inv_rs_2 * inv_rs_2;
    double inv_rs_6 = inv_rs_4 * inv_rs_2;
    double inv_rs_8 = inv_rs_6 * inv_rs_2;

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double P4x = P4(lambda, x0, gamma);
    double P4y = P4(lambda, y0, gamma);
    double P4z = P4(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P8x = P8(lambda, x0, gamma);
    double P8y = P8(lambda, y0, gamma);
    double P8z = P8(lambda, z0, gamma);

    double term_x = (2027025.0 * inv_rs_8 * P8x) - (3783780.0 * inv_rs_6 * P6x) + (2182950.0 * inv_rs_4 * P4x) - (396900.0 * inv_rs_2 * P2x) + 11025.0;
    double term_y = (2027025.0 * inv_rs_8 * P8y) - (3783780.0 * inv_rs_6 * P6y) + (2182950.0 * inv_rs_4 * P4y) - (396900.0 * inv_rs_2 * P2y) + 11025.0;
    double term_z = (2027025.0 * inv_rs_8 * P8z) - (3783780.0 * inv_rs_6 * P6z) + (2182950.0 * inv_rs_4 * P4z) - (396900.0 * inv_rs_2 * P2z) + 11025.0;

    double miu_1 = temp * term_x;
    double miu_2 = temp * term_y;
    double miu_3 = temp * term_z;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
}

void calc_MCSH_8_2_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);
    double inv_rs_2 = inv_rs * inv_rs;
    double inv_rs_4 = inv_rs_2 * inv_rs_2;
    double inv_rs_6 = inv_rs_4 * inv_rs_2;
    double inv_rs_8 = inv_rs_6 * inv_rs_2;

    double P1x = P1(lambda, x0, gamma);
    double P1y = P1(lambda, y0, gamma);
    double P1z = P1(lambda, z0, gamma);

    double P3x = P3(lambda, x0, gamma);
    double P3y = P3(lambda, y0, gamma);
    double P3z = P3(lambda, z0, gamma);

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double term_x = (2027025.0 * inv_rs_8 * P7x) - (2837835.0 * inv_rs_6 * P5x) + (1091475.0 * inv_rs_4 * P3x) - (99225.0 * inv_rs_2 * P1x);
    double term_y = (2027025.0 * inv_rs_8 * P7y) - (2837835.0 * inv_rs_6 * P5y) + (1091475.0 * inv_rs_4 * P3y) - (99225.0 * inv_rs_2 * P1y);
    double term_z = (2027025.0 * inv_rs_8 * P7z) - (2837835.0 * inv_rs_6 * P5z) + (1091475.0 * inv_rs_4 * P3z) - (99225.0 * inv_rs_2 * P1z);

    double miu_1 = temp * P1y * term_x;
    double miu_2 = temp * P1x * term_y;
    double miu_3 = temp * P1z * term_x;
    double miu_4 = temp * P1x * term_z;
    double miu_5 = temp * P1z * term_y;
    double miu_6 = temp * P1y * term_z;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;
}


void calc_MCSH_8_3_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);
    double inv_rs_2 = inv_rs * inv_rs;
    double inv_rs_4 = inv_rs_2 * inv_rs_2;
    double inv_rs_6 = inv_rs_4 * inv_rs_2;
    double inv_rs_8 = inv_rs_6 * inv_rs_2;

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double P4x = P4(lambda, x0, gamma);
    double P4y = P4(lambda, y0, gamma);
    double P4z = P4(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double temp_miu1 = (2027025.0 * inv_rs_8 * P6x * P2y) - (135135.0 * inv_rs_6 * P6x) - (2027025.0 * inv_rs_6 * P4x * P2y) + (155925.0 * inv_rs_4 * P4x) + (467775.0 * inv_rs_4 * P2x * P2y) - (42525.0 * inv_rs_2 * P2x) - (14175.0 * inv_rs_2 * P2y) + 1575.0;
    double temp_miu2 = (2027025.0 * inv_rs_8 * P6y * P2x) - (135135.0 * inv_rs_6 * P6y) - (2027025.0 * inv_rs_6 * P4y * P2x) + (155925.0 * inv_rs_4 * P4y) + (467775.0 * inv_rs_4 * P2y * P2x) - (42525.0 * inv_rs_2 * P2y) - (14175.0 * inv_rs_2 * P2x) + 1575.0;
    double temp_miu3 = (2027025.0 * inv_rs_8 * P6x * P2z) - (135135.0 * inv_rs_6 * P6x) - (2027025.0 * inv_rs_6 * P4x * P2z) + (155925.0 * inv_rs_4 * P4x) + (467775.0 * inv_rs_4 * P2x * P2z) - (42525.0 * inv_rs_2 * P2x) - (14175.0 * inv_rs_2 * P2z) + 1575.0;
    double temp_miu4 = (2027025.0 * inv_rs_8 * P6z * P2x) - (135135.0 * inv_rs_6 * P6z) - (2027025.0 * inv_rs_6 * P4z * P2x) + (155925.0 * inv_rs_4 * P4z) + (467775.0 * inv_rs_4 * P2z * P2x) - (42525.0 * inv_rs_2 * P2z) - (14175.0 * inv_rs_2 * P2x) + 1575.0;
    double temp_miu5 = (2027025.0 * inv_rs_8 * P6y * P2z) - (135135.0 * inv_rs_6 * P6y) - (2027025.0 * inv_rs_6 * P4y * P2z) + (155925.0 * inv_rs_4 * P4y) + (467775.0 * inv_rs_4 * P2y * P2z) - (42525.0 * inv_rs_2 * P2y) - (14175.0 * inv_rs_2 * P2z) + 1575.0;
    double temp_miu6 = (2027025.0 * inv_rs_8 * P6z * P2y) - (135135.0 * inv_rs_6 * P6z) - (2027025.0 * inv_rs_6 * P4z * P2y) + (155925.0 * inv_rs_4 * P4z) + (467775.0 * inv_rs_4 * P2z * P2y) - (42525.0 * inv_rs_2 * P2z) - (14175.0 * inv_rs_2 * P2y) + 1575.0;

    double miu_1 = temp * temp_miu1;
    double miu_2 = temp * temp_miu2;
    double miu_3 = temp * temp_miu3;
    double miu_4 = temp * temp_miu4;
    double miu_5 = temp * temp_miu5;
    double miu_6 = temp * temp_miu6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;
}


void calc_MCSH_8_4_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);
    double inv_rs_2 = inv_rs * inv_rs;
    double inv_rs_4 = inv_rs_2 * inv_rs_2;
    double inv_rs_6 = inv_rs_4 * inv_rs_2;
    double inv_rs_8 = inv_rs_6 * inv_rs_2;

    double P1x = P1(lambda, x0, gamma);
    double P1y = P1(lambda, y0, gamma);
    double P1z = P1(lambda, z0, gamma);

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double P4x = P4(lambda, x0, gamma);
    double P4y = P4(lambda, y0, gamma);
    double P4z = P4(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double temp_x = (2027025.0 * inv_rs_8 * P6x) - (2027025.0 * inv_rs_6 * P4x) + (467775.0 * inv_rs_4 * P2x) - (14175.0 * inv_rs_2);
    double temp_y = (2027025.0 * inv_rs_8 * P6y) - (2027025.0 * inv_rs_6 * P4y) + (467775.0 * inv_rs_4 * P2y) - (14175.0 * inv_rs_2);
    double temp_z = (2027025.0 * inv_rs_8 * P6z) - (2027025.0 * inv_rs_6 * P4z) + (467775.0 * inv_rs_4 * P2z) - (14175.0 * inv_rs_2);

    double miu_1 = temp * P1y * P1z * temp_x;
    double miu_2 = temp * P1x * P1z * temp_y;
    double miu_3 = temp * P1x * P1y * temp_z;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
}

void calc_MCSH_8_5_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);
    double inv_rs_2 = inv_rs * inv_rs;
    double inv_rs_4 = inv_rs_2 * inv_rs_2;
    double inv_rs_6 = inv_rs_4 * inv_rs_2;
    double inv_rs_8 = inv_rs_6 * inv_rs_2;

    double P1x = P1(lambda, x0, gamma);
    double P1y = P1(lambda, y0, gamma);
    double P1z = P1(lambda, z0, gamma);

    double P3x = P3(lambda, x0, gamma);
    double P3y = P3(lambda, y0, gamma);
    double P3z = P3(lambda, z0, gamma);

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double temp_miu1 = (2027025.0 * inv_rs_8 * P5x * P3y) - (405405.0 * inv_rs_6 * P5x * P1y) - (1351350.0 * inv_rs_6 * P3x * P3y) + (311850.0 * inv_rs_4 * P3x * P1y) + (155925.0 * inv_rs_4 * P1x * P3y) - (42525.0 * inv_rs_2 * P1x * P1y);
    double temp_miu2 = (2027025.0 * inv_rs_8 * P5y * P3x) - (405405.0 * inv_rs_6 * P5y * P1x) - (1351350.0 * inv_rs_6 * P3y * P3x) + (311850.0 * inv_rs_4 * P3y * P1x) + (155925.0 * inv_rs_4 * P1y * P3x) - (42525.0 * inv_rs_2 * P1y * P1x);
    double temp_miu3 = (2027025.0 * inv_rs_8 * P5x * P3z) - (405405.0 * inv_rs_6 * P5x * P1z) - (1351350.0 * inv_rs_6 * P3x * P3z) + (311850.0 * inv_rs_4 * P3x * P1z) + (155925.0 * inv_rs_4 * P1x * P3z) - (42525.0 * inv_rs_2 * P1x * P1z);
    double temp_miu4 = (2027025.0 * inv_rs_8 * P5z * P3x) - (405405.0 * inv_rs_6 * P5z * P1x) - (1351350.0 * inv_rs_6 * P3z * P3x) + (311850.0 * inv_rs_4 * P3z * P1x) + (155925.0 * inv_rs_4 * P1z * P3x) - (42525.0 * inv_rs_2 * P1z * P1x);
    double temp_miu5 = (2027025.0 * inv_rs_8 * P5y * P3z) - (405405.0 * inv_rs_6 * P5y * P1z) - (1351350.0 * inv_rs_6 * P3y * P3z) + (311850.0 * inv_rs_4 * P3y * P1z) + (155925.0 * inv_rs_4 * P1y * P3z) - (42525.0 * inv_rs_2 * P1y * P1z);
    double temp_miu6 = (2027025.0 * inv_rs_8 * P5z * P3y) - (405405.0 * inv_rs_6 * P5z * P1y) - (1351350.0 * inv_rs_6 * P3z * P3y) + (311850.0 * inv_rs_4 * P3z * P1y) + (155925.0 * inv_rs_4 * P1z * P3y) - (42525.0 * inv_rs_2 * P1z * P1y);

    double miu_1 = temp * temp_miu1;
    double miu_2 = temp * temp_miu2;
    double miu_3 = temp * temp_miu3;
    double miu_4 = temp * temp_miu4;
    double miu_5 = temp * temp_miu5;
    double miu_6 = temp * temp_miu6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;
}

void calc_MCSH_8_6_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);
    double inv_rs_2 = inv_rs * inv_rs;
    double inv_rs_4 = inv_rs_2 * inv_rs_2;
    double inv_rs_6 = inv_rs_4 * inv_rs_2;
    double inv_rs_8 = inv_rs_6 * inv_rs_2;

    double P1x = P1(lambda, x0, gamma);
    double P1y = P1(lambda, y0, gamma);
    double P1z = P1(lambda, z0, gamma);

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double P3x = P3(lambda, x0, gamma);
    double P3y = P3(lambda, y0, gamma);
    double P3z = P3(lambda, z0, gamma);

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double term_1 = (2027025.0 * inv_rs_8 * P5x * P2y) - (135135.0 * inv_rs_6 * P5x) - (1351350.0 * inv_rs_6 * P3x * P2y) + (103950.0 * inv_rs_4 * P3x) + (155925.0 * inv_rs_4 * P1x * P2y) - (14175.0 * inv_rs_2 * P1x);
    double term_2 = (2027025.0 * inv_rs_8 * P5y * P2x) - (135135.0 * inv_rs_6 * P5y) - (1351350.0 * inv_rs_6 * P3y * P2x) + (103950.0 * inv_rs_4 * P3y) + (155925.0 * inv_rs_4 * P1y * P2x) - (14175.0 * inv_rs_2 * P1y);
    double term_3 = (2027025.0 * inv_rs_8 * P5x * P2z) - (135135.0 * inv_rs_6 * P5x) - (1351350.0 * inv_rs_6 * P3x * P2z) + (103950.0 * inv_rs_4 * P3x) + (155925.0 * inv_rs_4 * P1x * P2z) - (14175.0 * inv_rs_2 * P1x);
    double term_4 = (2027025.0 * inv_rs_8 * P5z * P2x) - (135135.0 * inv_rs_6 * P5z) - (1351350.0 * inv_rs_6 * P3z * P2x) + (103950.0 * inv_rs_4 * P3z) + (155925.0 * inv_rs_4 * P1z * P2x) - (14175.0 * inv_rs_2 * P1z);
    double term_5 = (2027025.0 * inv_rs_8 * P5y * P2z) - (135135.0 * inv_rs_6 * P5y) - (1351350.0 * inv_rs_6 * P3y * P2z) + (103950.0 * inv_rs_4 * P3y) + (155925.0 * inv_rs_4 * P1y * P2z) - (14175.0 * inv_rs_2 * P1y);
    double term_6 = (2027025.0 * inv_rs_8 * P5z * P2y) - (135135.0 * inv_rs_6 * P5z) - (1351350.0 * inv_rs_6 * P3z * P2y) + (103950.0 * inv_rs_4 * P3z) + (155925.0 * inv_rs_4 * P1z * P2y) - (14175.0 * inv_rs_2 * P1z);

    double miu_1 = temp * P1z * term_1;
    double miu_2 = temp * P1z * term_2;
    double miu_3 = temp * P1y * term_3;
    double miu_4 = temp * P1y * term_4;
    double miu_5 = temp * P1x * term_5;
    double miu_6 = temp * P1x * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;
}

void calc_MCSH_8_7_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);
    double inv_rs_2 = inv_rs * inv_rs;
    double inv_rs_4 = inv_rs_2 * inv_rs_2;
    double inv_rs_6 = inv_rs_4 * inv_rs_2;
    double inv_rs_8 = inv_rs_6 * inv_rs_2;

    double P1x = P1(lambda, x0, gamma);
    double P1y = P1(lambda, y0, gamma);
    double P1z = P1(lambda, z0, gamma);

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double P4x = P4(lambda, x0, gamma);
    double P4y = P4(lambda, y0, gamma);
    double P4z = P4(lambda, z0, gamma);

    double temp_1 = (2027025.0 * inv_rs_8 * P4x * P4y) - (810810.0 * inv_rs_6 * ((P4x * P2y) + (P2x * P4y))) + (31185.0 * inv_rs_4 * (P4x + P4y)) + (374220.0 * inv_rs_4 * P2x * P2y) - (17010.0 * inv_rs_2 * (P2x + P2y)) + 945.0;
    double temp_2 = (2027025.0 * inv_rs_8 * P4x * P4z) - (810810.0 * inv_rs_6 * ((P4x * P2z) + (P2x * P4z))) + (31185.0 * inv_rs_4 * (P4x + P4z)) + (374220.0 * inv_rs_4 * P2x * P2z) - (17010.0 * inv_rs_2 * (P2x + P2z)) + 945.0;
    double temp_3 = (2027025.0 * inv_rs_8 * P4y * P4z) - (810810.0 * inv_rs_6 * ((P4y * P2z) + (P2y * P4z))) + (31185.0 * inv_rs_4 * (P4y + P4z)) + (374220.0 * inv_rs_4 * P2y * P2z) - (17010.0 * inv_rs_2 * (P2y + P2z)) + 945.0;

    double miu_1 = temp * temp_1;
    double miu_2 = temp * temp_2;
    double miu_3 = temp * temp_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
}

void calc_MCSH_8_8_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);
    double inv_rs_2 = inv_rs * inv_rs;
    double inv_rs_4 = inv_rs_2 * inv_rs_2;
    double inv_rs_6 = inv_rs_4 * inv_rs_2;
    double inv_rs_8 = inv_rs_6 * inv_rs_2;

    double P1x = P1(lambda, x0, gamma);
    double P1y = P1(lambda, y0, gamma);
    double P1z = P1(lambda, z0, gamma);

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double P3x = P3(lambda, x0, gamma);
    double P3y = P3(lambda, y0, gamma);
    double P3z = P3(lambda, z0, gamma);

    double P4x = P4(lambda, x0, gamma);
    double P4y = P4(lambda, y0, gamma);
    double P4z = P4(lambda, z0, gamma);

    double term_1 = (2027025.0 * inv_rs_8 * P4x * P3y) - (405405.0 * inv_rs_6 * P4x * P1y) - (810810.0 * inv_rs_6 * P2x * P3y) + (187110.0 * inv_rs_4 * P2x * P1y) + (31185.0 * inv_rs_4 * P3y) - (8505.0 * inv_rs_2 * P1y);
    double term_2 = (2027025.0 * inv_rs_8 * P4y * P3x) - (405405.0 * inv_rs_6 * P4y * P1x) - (810810.0 * inv_rs_6 * P2y * P3x) + (187110.0 * inv_rs_4 * P2y * P1x) + (31185.0 * inv_rs_4 * P3x) - (8505.0 * inv_rs_2 * P1x);
    double term_3 = (2027025.0 * inv_rs_8 * P4x * P3z) - (405405.0 * inv_rs_6 * P4x * P1z) - (810810.0 * inv_rs_6 * P2x * P3z) + (187110.0 * inv_rs_4 * P2x * P1z) + (31185.0 * inv_rs_4 * P3z) - (8505.0 * inv_rs_2 * P1z);
    double term_4 = (2027025.0 * inv_rs_8 * P4z * P3x) - (405405.0 * inv_rs_6 * P4z * P1x) - (810810.0 * inv_rs_6 * P2z * P3x) + (187110.0 * inv_rs_4 * P2z * P1x) + (31185.0 * inv_rs_4 * P3x) - (8505.0 * inv_rs_2 * P1x);
    double term_5 = (2027025.0 * inv_rs_8 * P4y * P3z) - (405405.0 * inv_rs_6 * P4y * P1z) - (810810.0 * inv_rs_6 * P2y * P3z) + (187110.0 * inv_rs_4 * P2y * P1z) + (31185.0 * inv_rs_4 * P3z) - (8505.0 * inv_rs_2 * P1z);
    double term_6 = (2027025.0 * inv_rs_8 * P4z * P3y) - (405405.0 * inv_rs_6 * P4z * P1y) - (810810.0 * inv_rs_6 * P2z * P3y) + (187110.0 * inv_rs_4 * P2z * P1y) + (31185.0 * inv_rs_4 * P3y) - (8505.0 * inv_rs_2 * P1y);

    double miu_1 = temp * P1z * term_1;
    double miu_2 = temp * P1z * term_2;
    double miu_3 = temp * P1y * term_3;
    double miu_4 = temp * P1y * term_4;
    double miu_5 = temp * P1x * term_5;
    double miu_6 = temp * P1x * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;
}

void calc_MCSH_8_9_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);
    double inv_rs_2 = inv_rs * inv_rs;
    double inv_rs_4 = inv_rs_2 * inv_rs_2;
    double inv_rs_6 = inv_rs_4 * inv_rs_2;
    double inv_rs_8 = inv_rs_6 * inv_rs_2;

    double P1x = P1(lambda, x0, gamma);
    double P1y = P1(lambda, y0, gamma);
    double P1z = P1(lambda, z0, gamma);

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double P4x = P4(lambda, x0, gamma);
    double P4y = P4(lambda, y0, gamma);
    double P4z = P4(lambda, z0, gamma);

    double temp_1 = (2027025.0 * inv_rs_8 * P4x * P2y * P2z) - (135135.0 * inv_rs_6 * P4x * (P2y + P2z)) + (10395.0 * inv_rs_4 * P4x) - (810810.0 * inv_rs_6 * P2x * P2y * P2z) + (62370.0 * inv_rs_4 * P2x * (P2y + P2z)) - (5670.0 * inv_rs_2 * P2x) + (31185.0 * inv_rs_4 * P2y * P2z) - (2835.0 * inv_rs_2 * (P2y + P2z)) + 315.0;
    double temp_2 = (2027025.0 * inv_rs_8 * P2x * P4y * P2z) - (135135.0 * inv_rs_6 * P4y * (P2x + P2z)) + (10395.0 * inv_rs_4 * P4y) - (810810.0 * inv_rs_6 * P2x * P2y * P2z) + (62370.0 * inv_rs_4 * P2y * (P2x + P2z)) - (5670.0 * inv_rs_2 * P2y) + (31185.0 * inv_rs_4 * P2x * P2z) - (2835.0 * inv_rs_2 * (P2x + P2z)) + 315.0;
    double temp_3 = (2027025.0 * inv_rs_8 * P2x * P2y * P4z) - (135135.0 * inv_rs_6 * P4z * (P2x + P2y)) + (10395.0 * inv_rs_4 * P4z) - (810810.0 * inv_rs_6 * P2x * P2y * P2z) + (62370.0 * inv_rs_4 * P2z * (P2x + P2y)) - (5670.0 * inv_rs_2 * P2z) + (31185.0 * inv_rs_4 * P2x * P2y) - (2835.0 * inv_rs_2 * (P2x + P2y)) + 315.0;

    double miu_1 = temp * temp_1;
    double miu_2 = temp * temp_2;
    double miu_3 = temp * temp_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
}

void calc_MCSH_8_10_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);
    double inv_rs_2 = inv_rs * inv_rs;
    double inv_rs_4 = inv_rs_2 * inv_rs_2;
    double inv_rs_6 = inv_rs_4 * inv_rs_2;
    double inv_rs_8 = inv_rs_6 * inv_rs_2;

    double P1x = P1(lambda, x0, gamma);
    double P1y = P1(lambda, y0, gamma);
    double P1z = P1(lambda, z0, gamma);

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double P3x = P3(lambda, x0, gamma);
    double P3y = P3(lambda, y0, gamma);
    double P3z = P3(lambda, z0, gamma);

    double temp_1 = (2027025.0 * inv_rs_8 * P3x * P3y * P2z) - (135135.0 * inv_rs_6 * P3x * P3y) - (405405.0 * inv_rs_6 * P2z * (P3x * P1y + P1x * P3y)) + (31185.0 * inv_rs_4 * (P3x * P1y + P1x * P3y)) + (93555.0 * inv_rs_4 * P1x * P1y * P2z) - (8505.0 * inv_rs_2 * P1x * P1y);
    double temp_2 = (2027025.0 * inv_rs_8 * P3x * P2y * P3z) - (135135.0 * inv_rs_6 * P3x * P3z) - (405405.0 * inv_rs_6 * P2y * (P3x * P1z + P1x * P3z)) + (31185.0 * inv_rs_4 * (P3x * P1z + P1x * P3z)) + (93555.0 * inv_rs_4 * P1x * P2y * P1z) - (8505.0 * inv_rs_2 * P1x * P1z);
    double temp_3 = (2027025.0 * inv_rs_8 * P2x * P3y * P3z) - (135135.0 * inv_rs_6 * P3y * P3z) - (405405.0 * inv_rs_6 * P2x * (P3y * P1z + P1y * P3z)) + (31185.0 * inv_rs_4 * (P3y * P1z + P1y * P3z)) + (93555.0 * inv_rs_4 * P2x * P1y * P1z) - (8505.0 * inv_rs_2 * P1y * P1z);

    double miu_1 = temp * temp_1;
    double miu_2 = temp * temp_2;
    double miu_3 = temp * temp_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
}




// void calc_MCSH_5_1_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
// {
//     // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
//     double C1 = calc_C1(A,B,alpha,beta);
//     double C2 = calc_C2(alpha, beta);

//     double lambda = calc_lambda(alpha, beta);
//     double inv_rs_3 = inv_rs * inv_rs * inv_rs;
//     double inv_rs_5 = inv_rs_3 * inv_rs * inv_rs;

//     double lambda_x0 = lambda * x0;
//     double lambda_y0 = lambda * y0;
//     double lambda_z0 = lambda * z0;

//     double lambda_x0_sqr = lambda_x0 * lambda_x0;
//     double lambda_y0_sqr = lambda_y0 * lambda_y0;
//     double lambda_z0_sqr = lambda_z0 * lambda_z0;

//     double lambda_x0_3 = lambda_x0_sqr * lambda_x0;
//     double lambda_y0_3 = lambda_y0_sqr * lambda_y0;
//     double lambda_z0_3 = lambda_z0_sqr * lambda_z0;

//     double lambda_x0_5 = lambda_x0_3 * lambda_x0_sqr;
//     double lambda_y0_5 = lambda_y0_3 * lambda_y0_sqr;
//     double lambda_z0_5 = lambda_z0_3 * lambda_z0_sqr;

//     double gamma = calc_gamma(alpha, beta);
//     double C3 = ((4725.0 * inv_rs_5) / gamma) - (1050.0 * inv_rs_3);
//     double C4 = ((14175.0 * inv_rs_5) / (4.0*gamma*gamma)) - ((1575.0 * inv_rs_3) / gamma) + (225.0 * inv_rs);

//     double temp = C1 * exp( C2 * r0_sqr);

//     double temp_x = 945.0 * inv_rs_5 * lambda_x0_5 + C3 * lambda_x0_3 + C4 * lambda_x0;
//     double temp_y = 945.0 * inv_rs_5 * lambda_y0_5 + C3 * lambda_y0_3 + C4 * lambda_y0;
//     double temp_z = 945.0 * inv_rs_5 * lambda_z0_5 + C3 * lambda_z0_3 + C4 * lambda_z0;

//     double miu_5_1_1 = temp * temp_x;
//     double miu_5_1_2 = temp * temp_y;
//     double miu_5_1_3 = temp * temp_z;

//     value[0] = miu_5_1_1;
//     value[1] = miu_5_1_2;
//     value[2] = miu_5_1_3;
// }

// void calc_MCSH_5_2_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
// {
//     // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
//     double C1 = calc_C1(A,B,alpha,beta);
//     double C2 = calc_C2(alpha, beta);

//     double lambda = calc_lambda(alpha, beta);
//     double inv_rs_3 = inv_rs * inv_rs * inv_rs;
//     double inv_rs_5 = inv_rs_3 * inv_rs * inv_rs;

//     double lambda_x0 = x0 * lambda;
//     double lambda_y0 = y0 * lambda;
//     double lambda_z0 = z0 * lambda;

//     double lambda_x0_sqr = lambda_x0 * lambda_x0;
//     double lambda_y0_sqr = lambda_y0 * lambda_y0;
//     double lambda_z0_sqr = lambda_z0 * lambda_z0;

//     double lambda_x0_4 = lambda_x0_sqr * lambda_x0_sqr;
//     double lambda_y0_4 = lambda_y0_sqr * lambda_y0_sqr;
//     double lambda_z0_4 = lambda_z0_sqr * lambda_z0_sqr;

//     double gamma = calc_gamma(alpha, beta);
//     double C3 = ((2835.0 * inv_rs_5) / gamma) - (630.0 * inv_rs_3);
//     double C4 = ((2835.0 * inv_rs_5) / (4.0 * gamma * gamma)) - ((315.0 * inv_rs_3) / gamma) + (45.0 * inv_rs);

//     double temp = C1 * exp( C2 * r0_sqr);

//     double temp_x = 945.0 * inv_rs_5 * lambda_x0_4 + C3 * lambda_x0_sqr + C4;
//     double temp_y = 945.0 * inv_rs_5 * lambda_y0_4 + C3 * lambda_y0_sqr + C4;
//     double temp_z = 945.0 * inv_rs_5 * lambda_z0_4 + C3 * lambda_z0_sqr + C4;

//     double miu_5_2_1 = temp * lambda_y0 * temp_x;
//     double miu_5_2_2 = temp * lambda_x0 * temp_y;
//     double miu_5_2_3 = temp * lambda_z0 * temp_x;
//     double miu_5_2_4 = temp * lambda_x0 * temp_z;
//     double miu_5_2_5 = temp * lambda_z0 * temp_y;
//     double miu_5_2_6 = temp * lambda_y0 * temp_z;

//     value[0] = miu_5_2_1;
//     value[1] = miu_5_2_2;
//     value[2] = miu_5_2_3;
//     value[3] = miu_5_2_4;
//     value[4] = miu_5_2_5;
//     value[5] = miu_5_2_6;
// }

// void calc_MCSH_5_3_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
// {
//     // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
//     double C1 = calc_C1(A,B,alpha,beta);
//     double C2 = calc_C2(alpha, beta);

//     double lambda = calc_lambda(alpha, beta);
//     double inv_rs_3 = inv_rs * inv_rs * inv_rs;
//     double inv_rs_5 = inv_rs_3 * inv_rs * inv_rs;

//     double lambda_x0 = x0 * lambda;
//     double lambda_y0 = y0 * lambda;
//     double lambda_z0 = z0 * lambda;

//     double lambda_x0_sqr = lambda_x0 * lambda_x0;
//     double lambda_y0_sqr = lambda_y0 * lambda_y0;
//     double lambda_z0_sqr = lambda_z0 * lambda_z0;

//     double lambda_x0_3 = lambda_x0_sqr * lambda_x0;
//     double lambda_y0_3 = lambda_y0_sqr * lambda_y0;
//     double lambda_z0_3 = lambda_z0_sqr * lambda_z0;

//     double gamma = calc_gamma(alpha, beta);
//     // double C3 = (945 / (2 * gamma)) - 105;
//     // double C4 = (-315 / ( 2 *gamma)) + 45;

//     double temp = C1 * exp( C2 * r0_sqr);

//     double temp_x_2 = lambda_x0_sqr + (1.0/(2.0*gamma));
//     double temp_y_2 = lambda_y0_sqr + (1.0/(2.0*gamma));
//     double temp_z_2 = lambda_z0_sqr + (1.0/(2.0*gamma));

//     double C3_1 = 3.0 / (2.0 * gamma);
//     double temp_x_3 = lambda_x0_3 + C3_1 * lambda_x0;
//     double temp_y_3 = lambda_y0_3 + C3_1 * lambda_y0;
//     double temp_z_3 = lambda_z0_3 + C3_1 * lambda_z0;

//     double temp_term1_x = (945.0 * inv_rs_5 * temp_x_3) - (315.0 * inv_rs_3 * lambda_x0);
//     double temp_term1_y = (945.0 * inv_rs_5 * temp_y_3) - (315.0 * inv_rs_3 * lambda_y0);
//     double temp_term1_z = (945.0 * inv_rs_5 * temp_z_3) - (315.0 * inv_rs_3 * lambda_z0);

//     double temp_term2_x = (-105.0 * inv_rs_3 * temp_x_3) + (45.0 * inv_rs * lambda_x0);
//     double temp_term2_y = (-105.0 * inv_rs_3 * temp_y_3) + (45.0 * inv_rs * lambda_y0);
//     double temp_term2_z = (-105.0 * inv_rs_3 * temp_z_3) + (45.0 * inv_rs * lambda_z0);

//     double temp_miu1 = temp_y_2 * temp_term1_x + temp_term2_x;
//     double temp_miu2 = temp_x_2 * temp_term1_y + temp_term2_y;
//     double temp_miu3 = temp_z_2 * temp_term1_x + temp_term2_x;
//     double temp_miu4 = temp_x_2 * temp_term1_z + temp_term2_z;
//     double temp_miu5 = temp_z_2 * temp_term1_y + temp_term2_y;
//     double temp_miu6 = temp_y_2 * temp_term1_z + temp_term2_z;

//     double miu_5_3_1 = temp * temp_miu1;
//     double miu_5_3_2 = temp * temp_miu2;
//     double miu_5_3_3 = temp * temp_miu3;
//     double miu_5_3_4 = temp * temp_miu4;
//     double miu_5_3_5 = temp * temp_miu5;
//     double miu_5_3_6 = temp * temp_miu6;

//     value[0] = miu_5_3_1;
//     value[1] = miu_5_3_2;
//     value[2] = miu_5_3_3;
//     value[3] = miu_5_3_4;
//     value[4] = miu_5_3_5;
//     value[5] = miu_5_3_6;
// }

// void calc_MCSH_5_4_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
// {
//     // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
//     double C1 = calc_C1(A,B,alpha,beta);
//     double C2 = calc_C2(alpha, beta);

//     double lambda = calc_lambda(alpha, beta);
//     double inv_rs_3 = inv_rs * inv_rs * inv_rs;
//     double inv_rs_5 = inv_rs_3 * inv_rs * inv_rs;

//     double lambda_x0 = x0 * lambda;
//     double lambda_y0 = y0 * lambda;
//     double lambda_z0 = z0 * lambda;

//     double lambda_x0_sqr = lambda_x0 * lambda_x0;
//     double lambda_y0_sqr = lambda_y0 * lambda_y0;
//     double lambda_z0_sqr = lambda_y0 * lambda_z0;

//     double lambda_x0_3 = lambda_x0_sqr * lambda_x0;
//     double lambda_y0_3 = lambda_y0_sqr * lambda_y0;
//     double lambda_z0_3 = lambda_z0_sqr * lambda_z0;

//     double gamma = calc_gamma(alpha, beta);
//     double C3 = (2835.0 * inv_rs_5 / (2.0 * gamma)) - (315.0 * inv_rs_3);

//     double temp = C1 * exp( C2 * r0_sqr);

//     double temp_x = 945.0 * inv_rs_5 * lambda_x0_3 + C3 * lambda_x0;
//     double temp_y = 945.0 * inv_rs_5 * lambda_y0_3 + C3 * lambda_y0;
//     double temp_z = 945.0 * inv_rs_5 * lambda_z0_3 + C3 * lambda_z0;

//     double miu_5_4_1 = temp * lambda_y0 * lambda_z0 * temp_x;
//     double miu_5_4_2 = temp * lambda_x0 * lambda_z0 * temp_y;
//     double miu_5_4_3 = temp * lambda_x0 * lambda_y0 * temp_z;

//     value[0] = miu_5_4_1;
//     value[1] = miu_5_4_2;
//     value[2] = miu_5_4_3;

// }

// void calc_MCSH_5_5_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value)
// {
//     // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
//     double C1 = calc_C1(A,B,alpha,beta);
//     double C2 = calc_C2(alpha, beta);

//     double lambda = calc_lambda(alpha, beta);
//     double inv_rs_3 = inv_rs * inv_rs * inv_rs;
//     double inv_rs_5 = inv_rs_3 * inv_rs * inv_rs;

//     double lambda_x0 = x0 * lambda;
//     double lambda_y0 = y0 * lambda;
//     double lambda_z0 = z0 * lambda;

//     double lambda_x0_sqr = lambda_x0 * lambda_x0;
//     double lambda_y0_sqr = lambda_y0 * lambda_y0;
//     double lambda_z0_sqr = lambda_z0 * lambda_z0;

//     double gamma = calc_gamma(alpha, beta);

//     double temp = C1 * exp( C2 * r0_sqr);

//     double temp_x_2 = lambda_x0_sqr + (1.0/(2.0*gamma));
//     double temp_y_2 = lambda_y0_sqr + (1.0/(2.0*gamma));
//     double temp_z_2 = lambda_z0_sqr + (1.0/(2.0*gamma));

//     double temp_miu1 = (945.0 * inv_rs_5 * temp_x_2 * temp_y_2) - (105.0 * inv_rs_3 * (temp_x_2 + temp_y_2)) + (15.0 * inv_rs);
//     double temp_miu2 = (945.0 * inv_rs_5 * temp_x_2 * temp_z_2) - (105.0 * inv_rs_3 * (temp_x_2 + temp_z_2)) + (15.0 * inv_rs);
//     double temp_miu3 = (945.0 * inv_rs_5 * temp_y_2 * temp_z_2) - (105.0 * inv_rs_3 * (temp_y_2 + temp_z_2)) + (15.0 * inv_rs);

//     double miu_5_5_1 = temp * lambda_z0 * temp_miu1;
//     double miu_5_5_2 = temp * lambda_y0 * temp_miu2;
//     double miu_5_5_3 = temp * lambda_x0 * temp_miu3;

//     value[0] = miu_5_5_1;
//     value[1] = miu_5_5_2;
//     value[2] = miu_5_5_3;
// }


int get_num_groups(int mcsh_order){
    if (mcsh_order == 0) { return 1;}
    else if (mcsh_order == 1) { return 1;}
    else if (mcsh_order == 2) { return 2;}
    else if (mcsh_order == 3) { return 3;}
    else if (mcsh_order == 4) { return 4;}
    else if (mcsh_order == 5) { return 5;}
    else if (mcsh_order == 6) { return 7;}
    else if (mcsh_order == 7) { return 8;}
    else if (mcsh_order == 8) { return 10;}
    else if (mcsh_order == 9) { return 12;}
    else {return 0;}
}

double get_group_coefficients(int mcsh_order, int group_num){
    if (mcsh_order == 0) {
        if (group_num == 1) {
            return 1.0;
        } else {
            return 0.0;
        }
    } else if (mcsh_order == 1) {
        if (group_num == 1) {
            return 1.0;
        } else {
            return 0.0;
        }
    } else if (mcsh_order == 2) {
        if (group_num == 1) {
            return 1.0;
        } else if (group_num == 2){
            return 2.0;
        } else {
            return 0.0;
        }
    } else if (mcsh_order == 3) {
        if (group_num == 1) {
            return 1.0;
        } else if (group_num == 2){
            return 3.0;
        } else if (group_num == 3){
            return 6.0;
        } else {
            return 0.0;
        }
    } else if (mcsh_order == 4) {
        if (group_num == 1) {
            return 1.0; // 4 0 0 -> 24.0 / (24.0 * 1.0 * 1.0)
        } else if (group_num == 2){
            return 4.0; // 3 1 0 -> 24.0 / (6.0 * 1.0 * 1.0)
        } else if (group_num == 3){
            return 6.0; // 2 2 0 -> 24.0 / (2.0 * 2.0 * 1.0)
        } else if (group_num == 4){
            return 12.0; // 2 1 1 -> 24.0 / (2.0 * 1.0 * 1.0)
        } else {
            return 0.0;
        }
    } else if (mcsh_order == 5) {
        if (group_num == 1) {
            return 1.0; // 5 0 0 -> 120.0 / (120.0 * 1.0 * 1.0)
        } else if (group_num == 2){
            return 5.0; // 4 1 0 -> 120.0 / (24.0 * 1.0 * 1.0)
        } else if (group_num == 3){
            return 10.0; // 3 2 0 -> 120.0 / (6.0 * 2.0 * 1.0)
        } else if (group_num == 4){
            return 20.0; // 3 1 1 -> 120.0 / (6.0 * 1.0 * 1.0)
        } else if (group_num == 5){
            return 30.0; // 3 1 1 -> 120.0 / (2.0 * 2.0 * 1.0)
        } else {
            return 0.0;
        }
    } else if (mcsh_order == 6) {
        if (group_num == 1) {
            return 1.0; // 6 0 0 -> 720.0 / (720.0 * 1.0 * 1.0)
        } else if (group_num == 2){
            return 6.0; // 5 1 0 -> 720.0 / (120.0 * 1.0 * 1.0)
        } else if (group_num == 3){
            return 15.0; // 4 2 0 -> 720.0 / (24.0 * 2.0 * 1.0)
        } else if (group_num == 4){
            return 30.0; // 4 1 1 -> 720.0 / (24.0 * 1.0 * 1.0)
        } else if (group_num == 5){
            return 20.0; // 3 3 0 -> 720.0 / (6.0 * 6.0 * 1.0)
        } else if (group_num == 6){
            return 60.0; // 3 2 1 -> 720.0 / (6.0 * 2.0 * 1.0)
        } else if (group_num == 7){
            return 90.0; // 2 2 2 -> 720.0 / (2.0 * 2.0 * 2.0)
        } else {
            return 0.0;
        }
    } else if (mcsh_order == 7) {
        if (group_num == 1) {
            return 1.0; // 7 0 0 -> 5040.0 / (5040.0 * 1.0 * 1.0)
        } else if (group_num == 2){
            return 7.0; // 6 1 0 -> 5040.0 / (720.0 * 1.0 * 1.0)
        } else if (group_num == 3){
            return 21.0; // 5 2 0 -> 5040.0 / (120.0 * 2.0 * 1.0)
        } else if (group_num == 4){
            return 42.0; // 5 1 1 -> 5040.0 / (120.0 * 1.0 * 1.0)
        } else if (group_num == 5){
            return 35.0; // 4 3 0 -> 5040.0 / (24.0 * 6.0 * 1.0)
        } else if (group_num == 6){
            return 105.0; // 4 2 1 -> 5040.0 / (24.0 * 2.0 * 1.0)
        } else if (group_num == 7){
            return 140.0; // 3 3 1 -> 5040.0 / (6.0 * 6.0 * 1.0)
        } else if (group_num == 8){
            return 210.0; // 3 2 2 -> 5040.0 / (6.0 * 2.0 * 2.0)
        } else {
            return 0.0;
        }
    } else if (mcsh_order == 8) {
        if (group_num == 1) {
            return 1.0; // 8 0 0 -> 40320.0 / (40320.0 * 1.0 * 1.0)
        } else if (group_num == 2){
            return 8.0; // 7 1 0 -> 40320.0 / (5040.0 * 1.0 * 1.0)
        } else if (group_num == 3){
            return 28.0; // 6 2 0 -> 40320.0 / (720.0 * 2.0 * 1.0)
        } else if (group_num == 4){
            return 56.0; // 6 1 1 -> 40320.0 / (720.0 * 1.0 * 1.0)
        } else if (group_num == 5){
            return 56.0; // 5 3 0 -> 40320.0 / (120.0 * 6.0 * 1.0)
        } else if (group_num == 6){
            return 168.0; // 5 2 1 -> 40320.0 / (120.0 * 2.0 * 1.0)
        } else if (group_num == 7){
            return 70.0; // 4 4 0 -> 40320.0 / (24.0 * 24.0 * 1.0)
        } else if (group_num == 8){
            return 280.0; // 4 3 1 -> 40320.0 / (24.0 * 6.0 * 1.0)
        } else if (group_num == 9){
            return 420.0; // 4 2 2 -> 40320.0 / (24.0 * 2.0 * 2.0)
        } else if (group_num == 10){
            return 560.0; // 3 3 2 -> 40320.0 / (6.0 * 6.0 * 2.0)
        } else {
            return 0.0;
        }
    } else {
        return 0.0;
    }
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
    } else if (mcsh_order == 6) {
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
        } else if (group_num == 6){
            return 3;
        } else if (group_num == 7){
            return 1;
        } else {
            return 0;
        }
    } else if (mcsh_order == 7) {
        if (group_num == 1) {
            return 2;
        } else if (group_num == 2){
            return 3;
        } else if (group_num == 3){
            return 3;
        } else if (group_num == 4){
            return 2;
        } else if (group_num == 5){
            return 3;
        } else if (group_num == 6){
            return 3;
        } else if (group_num == 7){
            return 2;
        } else if (group_num == 8){
            return 2;
        } else {
            return 0;
        }
    } else if (mcsh_order == 8) {
        if (group_num == 1) {
            return 2;
        } else if (group_num == 2){
            return 3;
        } else if (group_num == 3){
            return 3;
        } else if (group_num == 4){
            return 2;
        } else if (group_num == 5){
            return 3;
        } else if (group_num == 6){
            return 3;
        } else if (group_num == 7){
            return 2;
        } else if (group_num == 8){
            return 3;
        } else if (group_num == 9){
            return 2;
        } else if (group_num == 10){
            return 2;
        } else {
            return 0;
        }
    } else if (mcsh_order == 9) {
        if (group_num == 1) {
            return 2;
        } else if (group_num == 2){
            return 3;
        } else if (group_num == 3){
            return 3;
        } else if (group_num == 4){
            return 2;
        } else if (group_num == 5){
            return 3;
        } else if (group_num == 6){
            return 3;
        } else if (group_num == 7){
            return 3;
        } else if (group_num == 8){
            return 3;
        } else if (group_num == 9){
            return 2;
        } else if (group_num == 10){
            return 2;
        } else if (group_num == 11){
            return 3;
        } else if (group_num == 12){
            return 1;
        } else {
            return 0;
        }
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
    // } else if (mcsh_order == 5) {
    //     if (group_num == 1) {
    //         result = calc_MCSH_5_1;
    //     } else if (group_num == 2){
    //         result = calc_MCSH_5_2;
    //     } else if (group_num == 3){
    //         result = calc_MCSH_5_3;
    //     } else if (group_num == 4){
    //         result = calc_MCSH_5_4;
    //     } else if (group_num == 5){
    //         result = calc_MCSH_5_5;
    //     }
    // } else if (mcsh_order == 6) {
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
    }

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
    } else if (mcsh_order == 6) {
        if (group_num == 1) {
            result = calc_MCSH_6_1_noderiv;
        } else if (group_num == 2){
            result = calc_MCSH_6_2_noderiv;
        } else if (group_num == 3){
            result = calc_MCSH_6_3_noderiv;
        } else if (group_num == 4){
            result = calc_MCSH_6_4_noderiv;
        } else if (group_num == 5){
            result = calc_MCSH_6_5_noderiv;
        } else if (group_num == 6){
            result = calc_MCSH_6_6_noderiv;
        } else if (group_num == 7){
            result = calc_MCSH_6_7_noderiv;
        }
    } else if (mcsh_order == 7) {
        if (group_num == 1) {
            result = calc_MCSH_7_1_noderiv;
        } else if (group_num == 2){
            result = calc_MCSH_7_2_noderiv;
        } else if (group_num == 3){
            result = calc_MCSH_7_3_noderiv;
        } else if (group_num == 4){
            result = calc_MCSH_7_4_noderiv;
        } else if (group_num == 5){
            result = calc_MCSH_7_5_noderiv;
        } else if (group_num == 6){
            result = calc_MCSH_7_6_noderiv;
        } else if (group_num == 7){
            result = calc_MCSH_7_7_noderiv;
        } else if (group_num == 8){
            result = calc_MCSH_7_8_noderiv;
        }
    } else if (mcsh_order == 8) {
        if (group_num == 1) {
            result = calc_MCSH_8_1_noderiv;
        } else if (group_num == 2){
            result = calc_MCSH_8_2_noderiv;
        } else if (group_num == 3){
            result = calc_MCSH_8_3_noderiv;
        } else if (group_num == 4){
            result = calc_MCSH_8_4_noderiv;
        } else if (group_num == 5){
            result = calc_MCSH_8_5_noderiv;
        } else if (group_num == 6){
            result = calc_MCSH_8_6_noderiv;
        } else if (group_num == 7){
            result = calc_MCSH_8_7_noderiv;
        } else if (group_num == 8){
            result = calc_MCSH_8_8_noderiv;
        } else if (group_num == 9){
            result = calc_MCSH_8_9_noderiv;
        } else if (group_num == 10){
            result = calc_MCSH_8_10_noderiv;
        }
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
    }

    return result;
}