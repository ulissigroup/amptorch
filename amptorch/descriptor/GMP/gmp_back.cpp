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
    return  alpha / (alpha + beta);
}

double calc_gamma(double alpha, double beta){
    return alpha + beta;
}


void calc_MCSH_0_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
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

void calc_MCSH_1_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double lambda = calc_lambda(alpha, beta);


    double temp = C1 * exp( C2 * r0_sqr);


    double temp_x = lambda * x0;
    double temp_y = lambda * y0;
    double temp_z = lambda * z0;

    double temp_dx = lambda;
    double temp_dy = lambda;
    double temp_dz = lambda;

    double miu_1_1_1 = temp * temp_x;
    double miu_1_1_2 = temp * temp_y;
    double miu_1_1_3 = temp * temp_z;

    //deriv consts
    double const_2_C2_x = 2.0 * C2 * x0;
    double const_2_C2_y = 2.0 * C2 * y0;
    double const_2_C2_z = 2.0 * C2 * z0;

    // dmiu1 dx/dy/dz
    deriv[0] = temp * (temp_dx + const_2_C2_x * temp_x);
    deriv[1] = miu_1_1_1 * const_2_C2_y;
    deriv[2] = miu_1_1_1 * const_2_C2_z;

    // dmiu2 dx/dy/dz
    deriv[3] = miu_1_1_2 * const_2_C2_x;
    deriv[4] = temp * (temp_dy + const_2_C2_y * temp_y);
    deriv[5] = miu_1_1_2 * const_2_C2_z;

    // dmiu3 dx/dy/dz
    deriv[6] = miu_1_1_3 * const_2_C2_x;
    deriv[7] = miu_1_1_3 * const_2_C2_y;
    deriv[8] = temp * (temp_dz + const_2_C2_z * temp_z);

    value[0] = miu_1_1_1;
    value[1] = miu_1_1_2;
    value[2] = miu_1_1_3;

}

void calc_MCSH_2_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

    double lambda_x0 = x0 * lambda;
    double lambda_y0 = y0 * lambda;
    double lambda_z0 = z0 * lambda;

    double lambda_x0_sqr = lambda_x0 * lambda_x0;
    double lambda_y0_sqr = lambda_y0 * lambda_y0;
    double lambda_z0_sqr = lambda_z0 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = (3.0/(2.0*gamma)) - 1.0;

    double temp = C1 * exp( C2 * r0_sqr);


    double temp_x = 3.0 * lambda_x0_sqr + C3;
    double temp_y = 3.0 * lambda_y0_sqr + C3;
    double temp_z = 3.0 * lambda_z0_sqr + C3;

    double temp_dx = 6.0 * lambda * lambda_x0; // = lambda * (3 * 2 * lambda_x0)
    double temp_dy = 6.0 * lambda * lambda_y0;
    double temp_dz = 6.0 * lambda * lambda_z0;

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

void calc_MCSH_2_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double x0_sqr = x0*x0;
    double y0_sqr = y0*y0;
    double z0_sqr = z0*z0;

    double temp = C1 * exp( C2 * r0_sqr) * lambda * lambda * 3.0;

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

void calc_MCSH_3_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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
    double C3 = (45.0/(2.0*gamma)) - 9.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 15.0 * lambda_x0_3 + C3 * lambda_x0;
    double temp_y = 15.0 * lambda_y0_3 + C3 * lambda_y0;
    double temp_z = 15.0 * lambda_z0_3 + C3 * lambda_z0;

    double temp_dx = lambda * (45.0 * lambda_x0_sqr + C3);
    double temp_dy = lambda * (45.0 * lambda_y0_sqr + C3);
    double temp_dz = lambda * (45.0 * lambda_z0_sqr + C3);

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

void calc_MCSH_3_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
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
    double C3 = (15.0/(2.0*gamma)) - 3.0;

    double temp = C1 * exp( C2 * r0_sqr);// * lambda;

    double temp_x = 15.0 * lambda_x0_sqr + C3;
    double temp_y = 15.0 * lambda_y0_sqr + C3;
    double temp_z = 15.0 * lambda_z0_sqr + C3;

    double temp_dx = lambda * (30.0 * lambda_x0);
    double temp_dy = lambda * (30.0 * lambda_y0);
    double temp_dz = lambda * (30.0 * lambda_z0);

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

void calc_MCSH_3_3(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double lambda = calc_lambda(alpha, beta);

    double temp =  C1 * exp( C2 * r0_sqr) * lambda * lambda * lambda * 15.0;
    double m_3_3 = temp * x0 * y0 * z0;

    deriv[0] = temp * y0 * z0 * (1.0 + 2.0*C2*x0*x0);
    deriv[1] = temp * x0 * z0 * (1.0 + 2.0*C2*y0*y0);
    deriv[2] = temp * x0 * y0 * (1.0 + 2.0*C2*z0*z0);

    value[0] = m_3_3;
}

void calc_MCSH_4_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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
    double C3 = (315.0/gamma) - 90.0;
    double C4 = (315.0/(4.0*gamma*gamma)) - (45.0/gamma) + 9.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 105.0 * lambda_x0_4 + C3 * lambda_x0_sqr + C4;
    double temp_y = 105.0 * lambda_y0_4 + C3 * lambda_y0_sqr + C4;
    double temp_z = 105.0 * lambda_z0_4 + C3 * lambda_z0_sqr + C4;

    double temp_dx = lambda * (420.0 * lambda_x0_3 + 2.0 * C3 * lambda_x0);
    double temp_dy = lambda * (420.0 * lambda_y0_3 + 2.0 * C3 * lambda_y0);
    double temp_dz = lambda * (420.0 * lambda_z0_3 + 2.0 * C3 * lambda_z0);

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

void calc_MCSH_4_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
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

    double C3 = (315.0/(2.0*gamma)) - 45.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 105.0 * lambda_x0_3 + C3 * lambda_x0;
    double temp_y = 105.0 * lambda_y0_3 + C3 * lambda_y0;
    double temp_z = 105.0 * lambda_z0_3 + C3 * lambda_z0;

    double temp_dx = lambda * (315.0 * lambda_x0_sqr + C3);// lambda * (105.0 * 3.0 * lambda_x0_sqr + C3);
    double temp_dy = lambda * (315.0 * lambda_y0_sqr + C3);
    double temp_dz = lambda * (315.0 * lambda_z0_sqr + C3);

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


void calc_MCSH_4_3(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double lambda_sqr = lambda * lambda;
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

    double temp_term1_x = 105.0 * temp_x_2 - 15.0;
    double temp_term1_y = 105.0 * temp_y_2 - 15.0;
    // double temp_term1_z = 105.0 * temp_z_2 - 15.0;

    double temp_dterm1_dx = 105.0 * temp_dx_2;
    double temp_dterm1_dy = 105.0 * temp_dy_2;
    // double temp_dterm1_dz = 105.0 * temp_dz_2;

    double temp_term2_x = -15.0 * temp_x_2 + 3.0;
    double temp_term2_y = -15.0 * temp_y_2 + 3.0;
    // double temp_term2_z = -15.0 * temp_z_2 + 3.0;

    double temp_dterm2_dx = -15.0 * temp_dx_2;
    double temp_dterm2_dy = -15.0 * temp_dy_2;
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

void calc_MCSH_4_4(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    // double lambda_sqr = lambda * lambda;
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
    double C3 = (105.0 / (2.0 * gamma)) - 15.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 105.0 * lambda_x0_sqr + C3;
    double temp_y = 105.0 * lambda_y0_sqr + C3;
    double temp_z = 105.0 * lambda_z0_sqr + C3;

    double temp_dx = lambda * (210.0 * lambda_x0); // lambda * (105.0 * 2.0 * lambda_x0);
    double temp_dy = lambda * (210.0 * lambda_y0);
    double temp_dz = lambda * (210.0 * lambda_z0);

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


void calc_MCSH_5_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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
    double C3 = (4725.0 / gamma) - 1050.0;
    double C4 = (14175.0 / (4.0*gamma*gamma)) - (1575.0 / gamma) + 225.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 945.0 * lambda_x0_5 + C3 * lambda_x0_3 + C4 * lambda_x0;
    double temp_y = 945.0 * lambda_y0_5 + C3 * lambda_y0_3 + C4 * lambda_y0;
    double temp_z = 945.0 * lambda_z0_5 + C3 * lambda_z0_3 + C4 * lambda_z0;

    double temp_dx = lambda * (4725.0 * lambda_x0_4 + 3.0 * C3 * lambda_x0_sqr + C4);
    double temp_dy = lambda * (4725.0 * lambda_y0_4 + 3.0 * C3 * lambda_y0_sqr + C4);
    double temp_dz = lambda * (4725.0 * lambda_z0_4 + 3.0 * C3 * lambda_z0_sqr + C4);

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

void calc_MCSH_5_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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
    double C3 = (945.0 * 3.0 / gamma) - 630.0;
    double C4 = (945.0 * 3.0 / (4.0 * gamma * gamma)) - (630.0 / (2.0 * gamma)) + 45.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 945.0 * lambda_x0_4 + C3 * lambda_x0_sqr + C4;
    double temp_y = 945.0 * lambda_y0_4 + C3 * lambda_y0_sqr + C4;
    double temp_z = 945.0 * lambda_z0_4 + C3 * lambda_z0_sqr + C4;

    double temp_dx = lambda * (945.0 * 4.0 * lambda_x0_3 + 2.0 * C3 * lambda_x0);
    double temp_dy = lambda * (945.0 * 4.0 * lambda_y0_3 + 2.0 * C3 * lambda_y0);
    double temp_dz = lambda * (945.0 * 4.0 * lambda_z0_3 + 2.0 * C3 * lambda_z0);

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

void calc_MCSH_5_3(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double lambda_sqr = lambda * lambda;
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


    double temp_term1_x = 945.0 * temp_x_3 - 315.0 * lambda_x0;
    double temp_term1_y = 945.0 * temp_y_3 - 315.0 * lambda_y0;
    double temp_term1_z = 945.0 * temp_z_3 - 315.0 * lambda_z0;

    double temp_dterm1_dx = 945.0 * temp_dx_3 - 315.0 * lambda;
    double temp_dterm1_dy = 945.0 * temp_dy_3 - 315.0 * lambda;
    double temp_dterm1_dz = 945.0 * temp_dz_3 - 315.0 * lambda;

    double temp_term2_x = -105.0 * temp_x_3 + 45.0 * lambda_x0;
    double temp_term2_y = -105.0 * temp_y_3 + 45.0 * lambda_y0;
    double temp_term2_z = -105.0 * temp_z_3 + 45.0 * lambda_z0;

    double temp_dterm2_dx = -105.0 * temp_dx_3 + 45.0 * lambda;
    double temp_dterm2_dy = -105.0 * temp_dy_3 + 45.0 * lambda;
    double temp_dterm2_dz = -105.0 * temp_dz_3 + 45.0 * lambda;

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

void calc_MCSH_5_4(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double lambda_sqr = lambda * lambda;
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
    double C3 = (2835.0 / (2.0 * gamma)) - 315.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 945.0 * lambda_x0_3 + C3 * lambda_x0;
    double temp_y = 945.0 * lambda_y0_3 + C3 * lambda_y0;
    double temp_z = 945.0 * lambda_z0_3 + C3 * lambda_z0;

    double temp_dx = lambda * (2835.0 * lambda_x0_sqr + C3);
    double temp_dy = lambda * (2835.0 * lambda_y0_sqr + C3);
    double temp_dz = lambda * (2835.0 * lambda_z0_sqr + C3);

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

void calc_MCSH_5_5(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double lambda_sqr = lambda * lambda;
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

    double temp_miu1 = 945.0 * temp_x_2 * temp_y_2 - 105.0 * (temp_x_2 + temp_y_2) + 15.0;
    double temp_miu2 = 945.0 * temp_x_2 * temp_z_2 - 105.0 * (temp_x_2 + temp_z_2) + 15.0;
    double temp_miu3 = 945.0 * temp_y_2 * temp_z_2 - 105.0 * (temp_y_2 + temp_z_2) + 15.0;

    double temp_dmiu1_dx = 945.0 * temp_dx_2 * temp_y_2 - 105.0 * temp_dx_2;
    double temp_dmiu1_dy = 945.0 * temp_x_2 * temp_dy_2 - 105.0 * temp_dy_2;

    double temp_dmiu2_dx = 945.0 * temp_dx_2 * temp_z_2 - 105.0 * temp_dx_2;
    double temp_dmiu2_dz = 945.0 * temp_x_2 * temp_dz_2 - 105.0 * temp_dz_2;

    double temp_dmiu3_dy = 945.0 * temp_dy_2 * temp_z_2 - 105.0 * temp_dy_2;
    double temp_dmiu3_dz = 945.0 * temp_y_2 * temp_dz_2 - 105.0 * temp_dz_2;

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

void calc_MCSH_6_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    // double lambda_sqr = lambda * lambda;
    // double x0_sqr = x0*x0;
    // double y0_sqr = y0*y0;
    // double z0_sqr = z0*z0;

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

    double lambda_x0_6 = lambda_x0_5 * lambda_x0;
    double lambda_y0_6 = lambda_y0_5 * lambda_y0;
    double lambda_z0_6 = lambda_z0_5 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = (10395.0 * 15.0 / (2.0 * gamma)) - 14175.0;
    double C4 = (10395.0 * 45.0 / (4.0 * gamma * gamma)) - (14175.0 * 3.0 / gamma) + 4725.0;
    double C5 = (10395.0 * 15.0 / (8.0 * gamma * gamma * gamma)) - (14175.0 * 3.0 / (4.0 * gamma * gamma)) + (4725.0 / (2.0 * gamma)) - 225.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 10395.0 * lambda_x0_6 + C3 * lambda_x0_4 + C4 * lambda_x0_sqr + C5;
    double temp_y = 10395.0 * lambda_y0_6 + C3 * lambda_y0_4 + C4 * lambda_y0_sqr + C5;
    double temp_z = 10395.0 * lambda_z0_6 + C3 * lambda_z0_4 + C4 * lambda_z0_sqr + C5;

    double temp_dx = lambda * (10395.0 * 6.0 * lambda_x0_5 + 4.0 * C3 * lambda_x0_3 + 2.0 * C4 * lambda_x0);
    double temp_dy = lambda * (10395.0 * 6.0 * lambda_y0_5 + 4.0 * C3 * lambda_y0_3 + 2.0 * C4 * lambda_y0);
    double temp_dz = lambda * (10395.0 * 6.0 * lambda_z0_5 + 4.0 * C3 * lambda_z0_3 + 2.0 * C4 * lambda_z0);

    double miu_6_1_1 = temp * temp_x;
    double miu_6_1_2 = temp * temp_y;
    double miu_6_1_3 = temp * temp_z;

    //deriv consts
    double const_2_C2_x = 2.0 * C2 * x0;
    double const_2_C2_y = 2.0 * C2 * y0;
    double const_2_C2_z = 2.0 * C2 * z0;

    // dmiu1 dx/dy/dz
    deriv[0] = temp * (temp_dx + const_2_C2_x * temp_x);
    deriv[1] = miu_6_1_1 * const_2_C2_y;
    deriv[2] = miu_6_1_1 * const_2_C2_z;

    // dmiu2 dx/dy/dz
    deriv[3] = miu_6_1_2 * const_2_C2_x;
    deriv[4] = temp * (temp_dy + const_2_C2_y * temp_y);
    deriv[5] = miu_6_1_2 * const_2_C2_z;

    // dmiu3 dx/dy/dz
    deriv[6] = miu_6_1_3 * const_2_C2_x;
    deriv[7] = miu_6_1_3 * const_2_C2_y;
    deriv[8] = temp * (temp_dz + const_2_C2_z * temp_z);

    value[0] = miu_6_1_1;
    value[1] = miu_6_1_2;
    value[2] = miu_6_1_3;
}

void calc_MCSH_6_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
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

    double lambda_x0_5 = lambda_x0_4 * lambda_x0;
    double lambda_y0_5 = lambda_y0_4 * lambda_y0;
    double lambda_z0_5 = lambda_z0_4 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = (10395.0 * 5.0 / gamma) - 9450.0;
    double C4 = (10395.0 * 15.0 / ( 4.0 * gamma * gamma)) - (9450.0 * 3.0 / (2.0 * gamma)) + 1575.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 10395.0 * lambda_x0_5 + C3 * lambda_x0_3 + C4 * lambda_x0;
    double temp_y = 10395.0 * lambda_y0_5 + C3 * lambda_y0_3 + C4 * lambda_y0;
    double temp_z = 10395.0 * lambda_z0_5 + C3 * lambda_z0_3 + C4 * lambda_z0;

    double temp_dx = lambda * (10395.0 * 5.0 * lambda_x0_4 + 3.0 * C3 * lambda_x0_sqr + C4);
    double temp_dy = lambda * (10395.0 * 5.0 * lambda_y0_4 + 3.0 * C3 * lambda_y0_sqr + C4);
    double temp_dz = lambda * (10395.0 * 5.0 * lambda_z0_4 + 3.0 * C3 * lambda_z0_sqr + C4);

    double miu_6_2_1 = temp * lambda_y0 * temp_x;
    double miu_6_2_2 = temp * lambda_x0 * temp_y;
    double miu_6_2_3 = temp * lambda_z0 * temp_x;
    double miu_6_2_4 = temp * lambda_x0 * temp_z;
    double miu_6_2_5 = temp * lambda_z0 * temp_y;
    double miu_6_2_6 = temp * lambda_y0 * temp_z;

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
    deriv[2] = miu_6_2_1 * const_2_C2_z;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * lambda * temp_y * const_1_p_2_C2_x2;
    deriv[4] = temp * lambda_x0 * (temp_dy + const_2_C2_y * temp_y);
    deriv[5] = miu_6_2_2 * const_2_C2_z;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * lambda_z0 * (temp_dx + const_2_C2_x * temp_x);
    deriv[7] = miu_6_2_3 * const_2_C2_y;
    deriv[8] = temp * lambda * temp_x * const_1_p_2_C2_z2;

    // dmiu4 dx/dy/dz
    deriv[9] = temp * lambda * temp_z * const_1_p_2_C2_x2;
    deriv[10] = miu_6_2_4 * const_2_C2_y;
    deriv[11] = temp * lambda_x0 * (temp_dz + const_2_C2_z * temp_z);

    // dmiu5 dx/dy/dz
    deriv[12] = miu_6_2_5 * const_2_C2_x;
    deriv[13] = temp * lambda_z0 * (temp_dy + const_2_C2_y * temp_y);
    deriv[14] = temp * lambda * temp_y * const_1_p_2_C2_z2;

    // dmiu6 dx/dy/dz
    deriv[15] = miu_6_2_6 * const_2_C2_x;
    deriv[16] = temp * lambda * temp_z * const_1_p_2_C2_y2;
    deriv[17] = temp * lambda_y0 * (temp_dz + const_2_C2_z * temp_z);

    value[0] = miu_6_2_1;
    value[1] = miu_6_2_2;
    value[2] = miu_6_2_3;
    value[3] = miu_6_2_4;
    value[4] = miu_6_2_5;
    value[5] = miu_6_2_6;
}

void calc_MCSH_6_3(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double lambda_sqr = lambda * lambda;

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

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x_2 = lambda_x0_sqr + (1.0/(2.0*gamma));
    double temp_y_2 = lambda_y0_sqr + (1.0/(2.0*gamma));
    double temp_z_2 = lambda_z0_sqr + (1.0/(2.0*gamma));

    double temp_dx_2 = 2.0 * lambda_sqr * x0;
    double temp_dy_2 = 2.0 * lambda_sqr * y0;
    double temp_dz_2 = 2.0 * lambda_sqr * z0;

    double C4_1 = 3.0 / gamma, C4_2 = 3.0 / (4.0 * gamma * gamma);
    double temp_x_4 = lambda_x0_4 + C4_1 * lambda_x0_sqr + C4_2;
    double temp_y_4 = lambda_y0_4 + C4_1 * lambda_y0_sqr + C4_2;
    double temp_z_4 = lambda_z0_4 + C4_1 * lambda_z0_sqr + C4_2;

    double temp_dx_4 = lambda * (4.0 * lambda_x0_3 + 2.0 * C4_1 * lambda_x0);
    double temp_dy_4 = lambda * (4.0 * lambda_y0_3 + 2.0 * C4_1 * lambda_y0);
    double temp_dz_4 = lambda * (4.0 * lambda_z0_3 + 2.0 * C4_1 * lambda_z0);

    double temp_term1_x = 10395.0 * temp_x_4 - 5670.0 * temp_x_2 + 315.0;
    double temp_term1_y = 10395.0 * temp_y_4 - 5670.0 * temp_y_2 + 315.0;
    double temp_term1_z = 10395.0 * temp_z_4 - 5670.0 * temp_z_2 + 315.0;

    double temp_dterm1_dx = 10395.0 * temp_dx_4 - 5670.0 * temp_dx_2;
    double temp_dterm1_dy = 10395.0 * temp_dy_4 - 5670.0 * temp_dy_2;
    double temp_dterm1_dz = 10395.0 * temp_dz_4 - 5670.0 * temp_dz_2;

    double temp_term2_x = -945.0 * temp_x_4 + 630.0 * temp_x_2 - 45.0;
    double temp_term2_y = -945.0 * temp_y_4 + 630.0 * temp_y_2 - 45.0;
    double temp_term2_z = -945.0 * temp_z_4 + 630.0 * temp_z_2 - 45.0;

    double temp_dterm2_dx = -945.0 * temp_dx_4 + 630.0 * temp_dx_2;
    double temp_dterm2_dy = -945.0 * temp_dy_4 + 630.0 * temp_dy_2;
    double temp_dterm2_dz = -945.0 * temp_dz_4 + 630.0 * temp_dz_2;

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

    double miu_6_3_1 = temp * temp_miu1;
    double miu_6_3_2 = temp * temp_miu2;
    double miu_6_3_3 = temp * temp_miu3;
    double miu_6_3_4 = temp * temp_miu4;
    double miu_6_3_5 = temp * temp_miu5;
    double miu_6_3_6 = temp * temp_miu6;

    //deriv consts
    double const_2_C2_x = 2.0 * C2 * x0; // 2 * C2 * x0
    double const_2_C2_y = 2.0 * C2 * y0;
    double const_2_C2_z = 2.0 * C2 * z0;

    // dmiu1 dx/dy/dz
    deriv[0] = temp * (temp_dmiu1_dx + const_2_C2_x * temp_miu1);
    deriv[1] = temp * (temp_dmiu1_dy + const_2_C2_y * temp_miu1);
    deriv[2] = miu_6_3_1 * const_2_C2_z;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * (temp_dmiu2_dx + const_2_C2_x * temp_miu2);
    deriv[4] = temp * (temp_dmiu2_dy + const_2_C2_y * temp_miu2);
    deriv[5] = miu_6_3_2 * const_2_C2_z;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * (temp_dmiu3_dx + const_2_C2_x * temp_miu3);
    deriv[7] = miu_6_3_3 * const_2_C2_y;
    deriv[8] = temp * (temp_dmiu3_dz + const_2_C2_z * temp_miu3);

    // dmiu4 dx/dy/dz
    deriv[9] = temp * (temp_dmiu4_dx + const_2_C2_x * temp_miu4);
    deriv[10] = miu_6_3_4 * const_2_C2_y;
    deriv[11] = temp * (temp_dmiu4_dz + const_2_C2_z * temp_miu4);

    // dmiu5 dx/dy/dz
    deriv[12] = miu_6_3_5 * const_2_C2_x;
    deriv[13] = temp * (temp_dmiu5_dy + const_2_C2_y * temp_miu5);
    deriv[14] = temp * (temp_dmiu5_dz + const_2_C2_z * temp_miu5);

    // dmiu6 dx/dy/dz
    deriv[15] = miu_6_3_6 * const_2_C2_x;
    deriv[16] = temp * (temp_dmiu6_dy + const_2_C2_y * temp_miu6);
    deriv[17] = temp * (temp_dmiu6_dz + const_2_C2_z * temp_miu6);

    value[0] = miu_6_3_1;
    value[1] = miu_6_3_2;
    value[2] = miu_6_3_3;
    value[3] = miu_6_3_4;
    value[4] = miu_6_3_5;
    value[5] = miu_6_3_6;
}


void calc_MCSH_6_4(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
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
    double C3 = (10395.0 * 3.0 / gamma) - 5670.0;
    double C4 = (10395.0 * 3.0 / (4.0 * gamma * gamma)) - (5670.0 / (2.0 * gamma)) + 315.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 10395.0 * lambda_x0_4 + C3 * lambda_x0_sqr + C4;
    double temp_y = 10395.0 * lambda_y0_4 + C3 * lambda_y0_sqr + C4;
    double temp_z = 10395.0 * lambda_z0_4 + C3 * lambda_z0_sqr + C4;

    double temp_dx = lambda * (10395.0 * 4 * lambda_x0_3 + 2.0 * C3 * lambda_x0);
    double temp_dy = lambda * (10395.0 * 4 * lambda_y0_3 + 2.0 * C3 * lambda_y0);
    double temp_dz = lambda * (10395.0 * 4 * lambda_z0_3 + 2.0 * C3 * lambda_z0);

    double miu_6_4_1 = temp * lambda_y0 * lambda_z0 * temp_x;
    double miu_6_4_2 = temp * lambda_x0 * lambda_z0 * temp_y;
    double miu_6_4_3 = temp * lambda_x0 * lambda_y0 * temp_z;

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

    value[0] = miu_6_4_1;
    value[1] = miu_6_4_2;
    value[2] = miu_6_4_3;
}

void calc_MCSH_6_5(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    // double lambda_sqr = lambda * lambda;

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

    double temp = C1 * exp( C2 * r0_sqr);

    double C3_1 = 3.0 / (2.0 * gamma);
    double temp_x_3 = lambda_x0_3 + C3_1 * lambda_x0;
    double temp_y_3 = lambda_y0_3 + C3_1 * lambda_y0;
    double temp_z_3 = lambda_z0_3 + C3_1 * lambda_z0;

    double temp_dx_3 = lambda * (3.0 * lambda_x0_sqr + C3_1);
    double temp_dy_3 = lambda * (3.0 * lambda_y0_sqr + C3_1);
    double temp_dz_3 = lambda * (3.0 * lambda_z0_sqr + C3_1);

    double temp_miu1 = 10395.0 * temp_x_3 * temp_y_3 - 2835.0 * (temp_x_3 * lambda_y0 + lambda_x0 * temp_y_3) + 945.0 * lambda_x0 * lambda_y0;
    double temp_miu2 = 10395.0 * temp_x_3 * temp_z_3 - 2835.0 * (temp_x_3 * lambda_z0 + lambda_x0 * temp_z_3) + 945.0 * lambda_x0 * lambda_z0;
    double temp_miu3 = 10395.0 * temp_y_3 * temp_z_3 - 2835.0 * (temp_y_3 * lambda_z0 + lambda_y0 * temp_z_3) + 945.0 * lambda_y0 * lambda_z0;

    double temp_dmiu1_dx = 10395.0 * temp_dx_3 * temp_y_3 - 2835.0 * (temp_dx_3 * lambda_y0 + lambda * temp_y_3) + 945.0 * lambda * lambda_y0;
    double temp_dmiu1_dy = 10395.0 * temp_x_3 * temp_dy_3 - 2835.0 * (temp_x_3 * lambda + lambda_x0 * temp_dy_3) + 945.0 * lambda_x0 * lambda;

    double temp_dmiu2_dx = 10395.0 * temp_dx_3 * temp_z_3 - 2835.0 * (temp_dx_3 * lambda_z0 + lambda * temp_z_3) + 945.0 * lambda * lambda_z0;
    double temp_dmiu2_dz = 10395.0 * temp_x_3 * temp_dz_3 - 2835.0 * (temp_x_3 * lambda + lambda_x0 * temp_dz_3) + 945.0 * lambda_x0 * lambda;

    double temp_dmiu3_dy = 10395.0 * temp_dy_3 * temp_z_3 - 2835.0 * (temp_dy_3 * lambda_z0 + lambda * temp_z_3) + 945.0 * lambda * lambda_z0;
    double temp_dmiu3_dz = 10395.0 * temp_y_3 * temp_dz_3 - 2835.0 * (temp_y_3 * lambda + lambda_y0 * temp_dz_3) + 945.0 * lambda_y0 * lambda;

    double miu_6_5_1 = temp * temp_miu1;
    double miu_6_5_2 = temp * temp_miu2;
    double miu_6_5_3 = temp * temp_miu3;

    //deriv consts
    double const_2_C2_x = 2.0 * C2 * x0; // 2 * C2 * x0
    double const_2_C2_y = 2.0 * C2 * y0;
    double const_2_C2_z = 2.0 * C2 * z0;

    // dmiu1 dx/dy/dz
    deriv[0] = temp * (temp_dmiu1_dx + const_2_C2_x * temp_miu1);
    deriv[1] = temp * (temp_dmiu1_dy + const_2_C2_y * temp_miu1);
    deriv[2] = miu_6_5_1 * const_2_C2_z;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * (temp_dmiu2_dx + const_2_C2_x * temp_miu2);
    deriv[4] = miu_6_5_2 * const_2_C2_y;
    deriv[5] = temp * (temp_dmiu2_dz + const_2_C2_z * temp_miu2);

    // dmiu3 dx/dy/dz
    deriv[6] = miu_6_5_3 * const_2_C2_x;
    deriv[7] = temp * (temp_dmiu3_dy + const_2_C2_y * temp_miu3);
    deriv[8] = temp * (temp_dmiu3_dz + const_2_C2_z * temp_miu3);

    value[0] = miu_6_5_1;
    value[1] = miu_6_5_2;
    value[2] = miu_6_5_3;
}

void calc_MCSH_6_6(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double lambda_sqr = lambda * lambda;
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

    double gamma = calc_gamma(alpha, beta);

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

    double temp_term1_x = 10395.0 * temp_x_3 - 2835.0 * lambda_x0;
    double temp_term1_y = 10395.0 * temp_y_3 - 2835.0 * lambda_y0;
    double temp_term1_z = 10395.0 * temp_z_3 - 2835.0 * lambda_z0;

    double temp_dterm1_dx = 10395.0 * temp_dx_3 - 2835.0 * lambda;
    double temp_dterm1_dy = 10395.0 * temp_dy_3 - 2835.0 * lambda;
    double temp_dterm1_dz = 10395.0 * temp_dz_3 - 2835.0 * lambda;

    double temp_term2_x = -945.0 * temp_x_3 + 315.0 * lambda_x0;
    double temp_term2_y = -945.0 * temp_y_3 + 315.0 * lambda_y0;
    double temp_term2_z = -945.0 * temp_z_3 + 315.0 * lambda_z0;

    double temp_dterm2_dx = -945.0 * temp_dx_3 + 315.0 * lambda;
    double temp_dterm2_dy = -945.0 * temp_dy_3 + 315.0 * lambda;
    double temp_dterm2_dz = -945.0 * temp_dz_3 + 315.0 * lambda;

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

    double miu_6_6_1 = temp * lambda_z0 * temp_miu1;
    double miu_6_6_2 = temp * lambda_z0 * temp_miu2;
    double miu_6_6_3 = temp * lambda_y0 * temp_miu3;
    double miu_6_6_4 = temp * lambda_y0 * temp_miu4;
    double miu_6_6_5 = temp * lambda_x0 * temp_miu5;
    double miu_6_6_6 = temp * lambda_x0 * temp_miu6;

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

    // dmiu2 dx/dy/dz
    deriv[3] = temp * lambda_z0 * (temp_dmiu2_dx + const_2_C2_x * temp_miu2);
    deriv[4] = temp * lambda_z0 * (temp_dmiu2_dy + const_2_C2_y * temp_miu2);
    deriv[5] = temp * temp_miu2 * lambda * const_1_p_2_C2_z2;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * lambda_y0 * (temp_dmiu3_dx + const_2_C2_x * temp_miu3);
    deriv[7] = temp * temp_miu3 * lambda * const_1_p_2_C2_y2;
    deriv[8] = temp * lambda_y0 * (temp_dmiu3_dz + const_2_C2_z * temp_miu3);

    // dmiu4 dx/dy/dz
    deriv[9] = temp * lambda_y0 * (temp_dmiu4_dx + const_2_C2_x * temp_miu4);
    deriv[10] = temp * temp_miu4 * lambda * const_1_p_2_C2_y2;
    deriv[11] = temp * lambda_y0 * (temp_dmiu4_dz + const_2_C2_z * temp_miu4);

    // dmiu5 dx/dy/dz
    deriv[12] = temp * temp_miu5 * lambda * const_1_p_2_C2_x2;
    deriv[13] = temp * lambda_x0 * (temp_dmiu5_dy + const_2_C2_y * temp_miu5);
    deriv[14] = temp * lambda_x0 * (temp_dmiu5_dz + const_2_C2_z * temp_miu5);

    // dmiu6 dx/dy/dz
    deriv[15] = temp * temp_miu6 * lambda * const_1_p_2_C2_x2;
    deriv[16] = temp * lambda_x0 * (temp_dmiu6_dy + const_2_C2_y * temp_miu6);
    deriv[17] = temp * lambda_x0 * (temp_dmiu6_dz + const_2_C2_z * temp_miu6);

    value[0] = miu_6_6_1;
    value[1] = miu_6_6_2;
    value[2] = miu_6_6_3;
    value[3] = miu_6_6_4;
    value[4] = miu_6_6_5;
    value[5] = miu_6_6_6;
}


void calc_MCSH_6_7(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double lambda_sqr = lambda * lambda;

    double lambda_x0 = x0 * lambda;
    double lambda_y0 = y0 * lambda;
    double lambda_z0 = z0 * lambda;

    double lambda_x0_sqr = lambda_x0 * lambda_x0;
    double lambda_y0_sqr = lambda_y0 * lambda_y0;
    double lambda_z0_sqr = lambda_z0 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);

    double temp_x = lambda_x0_sqr + (1.0/(2.0*gamma));
    double temp_y = lambda_y0_sqr + (1.0/(2.0*gamma));
    double temp_z = lambda_z0_sqr + (1.0/(2.0*gamma));

    double temp_dx = 2.0 * lambda_sqr * x0;
    double temp_dy = 2.0 * lambda_sqr * y0;
    double temp_dz = 2.0 * lambda_sqr * z0;

    double t1 = 10395.0 * temp_x * temp_y * temp_z;
    double t2 = -945.0 * temp_x * temp_y;
    double t3 = -945.0 * temp_x * temp_z;
    double t4 = -945.0 * temp_y * temp_z;
    double t5 = 105.0 * temp_x;
    double t6 = 105.0 * temp_y;
    double t7 = 105.0 * temp_z;
    double t8 = -15.0;

    double sum_ts = t1 + t2 + t3 + t4 + t5 + t6 + t7 + t8;

    double sum_dtdx = temp_dx * (10395.0 * temp_y * temp_z - 945.0 * (temp_y + temp_z) + 105.0);
    double sum_dtdy = temp_dy * (10395.0 * temp_x * temp_z - 945.0 * (temp_x + temp_z) + 105.0);
    double sum_dtdz = temp_dz * (10395.0 * temp_x * temp_y - 945.0 * (temp_x + temp_y) + 105.0);

    double temp =  C1 * exp( C2 * r0_sqr);
    double m_6_7 = temp * sum_ts;

    deriv[0] = temp * (2.0 * C2 * x0 * sum_ts + sum_dtdx);
    deriv[1] = temp * (2.0 * C2 * y0 * sum_ts + sum_dtdy);
    deriv[2] = temp * (2.0 * C2 * z0 * sum_ts + sum_dtdz);

    value[0] = m_6_7;
}

void calc_MCSH_7_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    // double lambda_sqr = lambda * lambda;


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

    double lambda_x0_6 = lambda_x0_5 * lambda_x0;
    double lambda_y0_6 = lambda_y0_5 * lambda_y0;
    double lambda_z0_6 = lambda_z0_5 * lambda_z0;

    double lambda_x0_7 = lambda_x0_6 * lambda_x0;
    double lambda_y0_7 = lambda_y0_6 * lambda_y0;
    double lambda_z0_7 = lambda_z0_6 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = (135135.0 * 21.0 / (2.0 * gamma)) - 218295.0;
    double C4 = (135135.0 * 105.0 / (4.0 * gamma * gamma)) - (218295.0 * 5.0 / gamma) + 99225.0;
    double C5 = (135135.0 * 105.0 / (8.0 * gamma * gamma * gamma)) - (218295.0 * 15.0 / (4.0 * gamma * gamma)) + (99225.0 * 3.0 / (2.0 * gamma)) - 11025.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 135135.0 * lambda_x0_7 + C3 * lambda_x0_5 + C4 * lambda_x0_3 + C5 * lambda_x0;
    double temp_y = 135135.0 * lambda_y0_7 + C3 * lambda_y0_5 + C4 * lambda_y0_3 + C5 * lambda_y0;
    double temp_z = 135135.0 * lambda_z0_7 + C3 * lambda_z0_5 + C4 * lambda_z0_3 + C5 * lambda_z0;

    double temp_dx = lambda * (135135.0 * 7.0 * lambda_x0_6 + 5.0 * C3 * lambda_x0_4 + 3.0 * C4 * lambda_x0_sqr + C5);
    double temp_dy = lambda * (135135.0 * 7.0 * lambda_y0_6 + 5.0 * C3 * lambda_y0_4 + 3.0 * C4 * lambda_y0_sqr + C5);
    double temp_dz = lambda * (135135.0 * 7.0 * lambda_z0_6 + 5.0 * C3 * lambda_z0_4 + 3.0 * C4 * lambda_z0_sqr + C5);

    double miu_7_1_1 = temp * temp_x;
    double miu_7_1_2 = temp * temp_y;
    double miu_7_1_3 = temp * temp_z;

    //deriv consts
    double const_2_C2_x = 2.0 * C2 * x0;
    double const_2_C2_y = 2.0 * C2 * y0;
    double const_2_C2_z = 2.0 * C2 * z0;

    // dmiu1 dx/dy/dz
    deriv[0] = temp * (temp_dx + const_2_C2_x * temp_x);
    deriv[1] = miu_7_1_1 * const_2_C2_y;
    deriv[2] = miu_7_1_1 * const_2_C2_z;

    // dmiu2 dx/dy/dz
    deriv[3] = miu_7_1_2 * const_2_C2_x;
    deriv[4] = temp * (temp_dy + const_2_C2_y * temp_y);
    deriv[5] = miu_7_1_2 * const_2_C2_z;

    // dmiu3 dx/dy/dz
    deriv[6] = miu_7_1_3 * const_2_C2_x;
    deriv[7] = miu_7_1_3 * const_2_C2_y;
    deriv[8] = temp * (temp_dz + const_2_C2_z * temp_z);

    value[0] = miu_7_1_1;
    value[1] = miu_7_1_2;
    value[2] = miu_7_1_3;
}

void calc_MCSH_7_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
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

    double lambda_x0_5 = lambda_x0_4 * lambda_x0;
    double lambda_y0_5 = lambda_y0_4 * lambda_y0;
    double lambda_z0_5 = lambda_z0_4 * lambda_z0;

    double lambda_x0_6 = lambda_x0_5 * lambda_x0;
    double lambda_y0_6 = lambda_y0_5 * lambda_y0;
    double lambda_z0_6 = lambda_z0_5 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = (135135.0 * 15.0 / (2.0 * gamma)) - 155925.0;
    double C4 = (135135.0 * 45.0 / (4.0 * gamma * gamma)) - (155925.0 * 3.0 / gamma) + 42525.0;
    double C5 = (135135.0 * 15.0 / (8.0 * gamma * gamma * gamma)) - (155925.0 * 3.0 / (4.0 * gamma * gamma)) + (42525.0 / (2.0 * gamma)) - 1575.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 135135.0 * lambda_x0_6 + C3 * lambda_x0_4 + C4 * lambda_x0_sqr + C5;
    double temp_y = 135135.0 * lambda_y0_6 + C3 * lambda_y0_4 + C4 * lambda_y0_sqr + C5;
    double temp_z = 135135.0 * lambda_z0_6 + C3 * lambda_z0_4 + C4 * lambda_z0_sqr + C5;

    double temp_dx = lambda * (135135.0 * 6.0 * lambda_x0_5 + 4.0 * C3 * lambda_x0_3 + 2.0 * C4 * lambda_x0);
    double temp_dy = lambda * (135135.0 * 6.0 * lambda_y0_5 + 4.0 * C3 * lambda_y0_3 + 2.0 * C4 * lambda_y0);
    double temp_dz = lambda * (135135.0 * 6.0 * lambda_z0_5 + 4.0 * C3 * lambda_z0_3 + 2.0 * C4 * lambda_z0);

    double miu_7_2_1 = temp * lambda_y0 * temp_x;
    double miu_7_2_2 = temp * lambda_x0 * temp_y;
    double miu_7_2_3 = temp * lambda_z0 * temp_x;
    double miu_7_2_4 = temp * lambda_x0 * temp_z;
    double miu_7_2_5 = temp * lambda_z0 * temp_y;
    double miu_7_2_6 = temp * lambda_y0 * temp_z;

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
    deriv[2] = miu_7_2_1 * const_2_C2_z;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * lambda * temp_y * const_1_p_2_C2_x2;
    deriv[4] = temp * lambda_x0 * (temp_dy + const_2_C2_y * temp_y);
    deriv[5] = miu_7_2_2 * const_2_C2_z;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * lambda_z0 * (temp_dx + const_2_C2_x * temp_x);
    deriv[7] = miu_7_2_3 * const_2_C2_y;
    deriv[8] = temp * lambda * temp_x * const_1_p_2_C2_z2;

    // dmiu4 dx/dy/dz
    deriv[9] = temp * lambda * temp_z * const_1_p_2_C2_x2;
    deriv[10] = miu_7_2_4 * const_2_C2_y;
    deriv[11] = temp * lambda_x0 * (temp_dz + const_2_C2_z * temp_z);

    // dmiu5 dx/dy/dz
    deriv[12] = miu_7_2_5 * const_2_C2_x;
    deriv[13] = temp * lambda_z0 * (temp_dy + const_2_C2_y * temp_y);
    deriv[14] = temp * lambda * temp_y * const_1_p_2_C2_z2;

    // dmiu6 dx/dy/dz
    deriv[15] = miu_7_2_6 * const_2_C2_x;
    deriv[16] = temp * lambda * temp_z * const_1_p_2_C2_y2;
    deriv[17] = temp * lambda_y0 * (temp_dz + const_2_C2_z * temp_z);

    value[0] = miu_7_2_1;
    value[1] = miu_7_2_2;
    value[2] = miu_7_2_3;
    value[3] = miu_7_2_4;
    value[4] = miu_7_2_5;
    value[5] = miu_7_2_6;
}


void calc_MCSH_7_3(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double lambda_sqr = lambda * lambda;

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

    double lambda_x0_5 = lambda_x0_4 * lambda_x0;
    double lambda_y0_5 = lambda_y0_4 * lambda_y0;
    double lambda_z0_5 = lambda_z0_4 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);

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

    double C5_1 = 5.0 / gamma, C5_2 = 15.0 / (4.0 * gamma * gamma);
    double temp_x_5 = lambda_x0_5 + C5_1 * lambda_x0_3 + C5_2 * lambda_x0;
    double temp_y_5 = lambda_y0_5 + C5_1 * lambda_y0_3 + C5_2 * lambda_y0;
    double temp_z_5 = lambda_z0_5 + C5_1 * lambda_z0_3 + C5_2 * lambda_z0;

    double temp_dx_5 = lambda * (5.0 * lambda_x0_4 + 3.0 * C5_1 * lambda_x0_sqr + C5_2);
    double temp_dy_5 = lambda * (5.0 * lambda_y0_4 + 3.0 * C5_1 * lambda_y0_sqr + C5_2);
    double temp_dz_5 = lambda * (5.0 * lambda_z0_4 + 3.0 * C5_1 * lambda_z0_sqr + C5_2);

    double temp_term1_x = 135135.0 * temp_x_5 - 103950.0 * temp_x_3 + 14175.0 * lambda_x0;
    double temp_term1_y = 135135.0 * temp_y_5 - 103950.0 * temp_y_3 + 14175.0 * lambda_y0;
    double temp_term1_z = 135135.0 * temp_z_5 - 103950.0 * temp_z_3 + 14175.0 * lambda_z0;

    double temp_dterm1_dx = 135135.0 * temp_dx_5 - 103950.0 * temp_dx_3 + 14175.0 * lambda;
    double temp_dterm1_dy = 135135.0 * temp_dy_5 - 103950.0 * temp_dy_3 + 14175.0 * lambda;
    double temp_dterm1_dz = 135135.0 * temp_dz_5 - 103950.0 * temp_dz_3 + 14175.0 * lambda;

    double temp_term2_x = -10395.0 * temp_x_5 + 9450.0 * temp_x_3 - 1575.0 * lambda_x0;
    double temp_term2_y = -10395.0 * temp_y_5 + 9450.0 * temp_y_3 - 1575.0 * lambda_y0;
    double temp_term2_z = -10395.0 * temp_z_5 + 9450.0 * temp_z_3 - 1575.0 * lambda_z0;

    double temp_dterm2_dx = -10395.0 * temp_dx_5 + 9450.0 * temp_dx_3 - 1575.0 * lambda;
    double temp_dterm2_dy = -10395.0 * temp_dy_5 + 9450.0 * temp_dy_3 - 1575.0 * lambda;
    double temp_dterm2_dz = -10395.0 * temp_dz_5 + 9450.0 * temp_dz_3 - 1575.0 * lambda;

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

    double miu_7_3_1 = temp * temp_miu1;
    double miu_7_3_2 = temp * temp_miu2;
    double miu_7_3_3 = temp * temp_miu3;
    double miu_7_3_4 = temp * temp_miu4;
    double miu_7_3_5 = temp * temp_miu5;
    double miu_7_3_6 = temp * temp_miu6;

    //deriv consts
    double const_2_C2_x = 2.0 * C2 * x0; // 2 * C2 * x0
    double const_2_C2_y = 2.0 * C2 * y0;
    double const_2_C2_z = 2.0 * C2 * z0;

    // dmiu1 dx/dy/dz
    deriv[0] = temp * (temp_dmiu1_dx + const_2_C2_x * temp_miu1);
    deriv[1] = temp * (temp_dmiu1_dy + const_2_C2_y * temp_miu1);
    deriv[2] = miu_7_3_1 * const_2_C2_z;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * (temp_dmiu2_dx + const_2_C2_x * temp_miu2);
    deriv[4] = temp * (temp_dmiu2_dy + const_2_C2_y * temp_miu2);
    deriv[5] = miu_7_3_2 * const_2_C2_z;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * (temp_dmiu3_dx + const_2_C2_x * temp_miu3);
    deriv[7] = miu_7_3_3 * const_2_C2_y;
    deriv[8] = temp * (temp_dmiu3_dz + const_2_C2_z * temp_miu3);

    // dmiu4 dx/dy/dz
    deriv[9] = temp * (temp_dmiu4_dx + const_2_C2_x * temp_miu4);
    deriv[10] = miu_7_3_4 * const_2_C2_y;
    deriv[11] = temp * (temp_dmiu4_dz + const_2_C2_z * temp_miu4);

    // dmiu5 dx/dy/dz
    deriv[12] = miu_7_3_5 * const_2_C2_x;
    deriv[13] = temp * (temp_dmiu5_dy + const_2_C2_y * temp_miu5);
    deriv[14] = temp * (temp_dmiu5_dz + const_2_C2_z * temp_miu5);

    // dmiu6 dx/dy/dz
    deriv[15] = miu_7_3_6 * const_2_C2_x;
    deriv[16] = temp * (temp_dmiu6_dy + const_2_C2_y * temp_miu6);
    deriv[17] = temp * (temp_dmiu6_dz + const_2_C2_z * temp_miu6);

    value[0] = miu_7_3_1;
    value[1] = miu_7_3_2;
    value[2] = miu_7_3_3;
    value[3] = miu_7_3_4;
    value[4] = miu_7_3_5;
    value[5] = miu_7_3_6;
}


void calc_MCSH_7_4(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
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

    double lambda_x0_5 = lambda_x0_4 * lambda_x0;
    double lambda_y0_5 = lambda_y0_4 * lambda_y0;
    double lambda_z0_5 = lambda_z0_4 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = (135135.0 * 5.0 / gamma) - 103950.0;
    double C4 = (135135.0 * 15.0 / (4.0 * gamma * gamma)) - (103950.0 * 3.0 / (2.0 * gamma)) + 14175.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 135135.0 * lambda_x0_5 + C3 * lambda_x0_3 + C4 * lambda_x0;
    double temp_y = 135135.0 * lambda_y0_5 + C3 * lambda_y0_3 + C4 * lambda_y0;
    double temp_z = 135135.0 * lambda_z0_5 + C3 * lambda_z0_3 + C4 * lambda_z0;

    double temp_dx = lambda * (135135.0 * 5.0 * lambda_x0_4 + 3.0 * C3 * lambda_x0_sqr + C4);
    double temp_dy = lambda * (135135.0 * 5.0 * lambda_y0_4 + 3.0 * C3 * lambda_y0_sqr + C4);
    double temp_dz = lambda * (135135.0 * 5.0 * lambda_z0_4 + 3.0 * C3 * lambda_z0_sqr + C4);

    double miu_7_4_1 = temp * lambda_y0 * lambda_z0 * temp_x;
    double miu_7_4_2 = temp * lambda_x0 * lambda_z0 * temp_y;
    double miu_7_4_3 = temp * lambda_x0 * lambda_y0 * temp_z;

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

    value[0] = miu_7_4_1;
    value[1] = miu_7_4_2;
    value[2] = miu_7_4_3;
}

void calc_MCSH_7_5(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double lambda_sqr = lambda * lambda;

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

    double C4_1 = 3.0 / gamma, C4_2 = 3.0 / (4.0 * gamma * gamma);
    double temp_x_4 = lambda_x0_4 + C4_1 * lambda_x0_sqr + C4_2;
    double temp_y_4 = lambda_y0_4 + C4_1 * lambda_y0_sqr + C4_2;
    double temp_z_4 = lambda_z0_4 + C4_1 * lambda_z0_sqr + C4_2;

    double temp_dx_4 = lambda * (4.0 * lambda_x0_3 + 2.0 * C4_1 * lambda_x0);
    double temp_dy_4 = lambda * (4.0 * lambda_y0_3 + 2.0 * C4_1 * lambda_y0);
    double temp_dz_4 = lambda * (4.0 * lambda_z0_3 + 2.0 * C4_1 * lambda_z0);


    double temp_term1_x = 135135.0 * temp_x_4 - 62370.0 * temp_x_2 + 2835.0;
    double temp_term1_y = 135135.0 * temp_y_4 - 62370.0 * temp_y_2 + 2835.0;
    double temp_term1_z = 135135.0 * temp_z_4 - 62370.0 * temp_z_2 + 2835.0;

    double temp_dterm1_dx = 135135.0 * temp_dx_4 - 62370.0 * temp_dx_2;
    double temp_dterm1_dy = 135135.0 * temp_dy_4 - 62370.0 * temp_dy_2;
    double temp_dterm1_dz = 135135.0 * temp_dz_4 - 62370.0 * temp_dz_2;

    double temp_term2_x = -31185.0 * temp_x_4 + 17010.0 * temp_x_2 - 945.0;
    double temp_term2_y = -31185.0 * temp_y_4 + 17010.0 * temp_y_2 - 945.0;
    double temp_term2_z = -31185.0 * temp_z_4 + 17010.0 * temp_z_2 - 945.0;

    double temp_dterm2_dx = -31185.0 * temp_dx_4 + 17010.0 * temp_dx_2;
    double temp_dterm2_dy = -31185.0 * temp_dy_4 + 17010.0 * temp_dy_2;
    double temp_dterm2_dz = -31185.0 * temp_dz_4 + 17010.0 * temp_dz_2;

    double temp_miu1 = temp_y_3 * temp_term1_x + lambda_y0 * temp_term2_x;
    double temp_miu2 = temp_x_3 * temp_term1_y + lambda_x0 * temp_term2_y;
    double temp_miu3 = temp_z_3 * temp_term1_x + lambda_z0 * temp_term2_x;
    double temp_miu4 = temp_x_3 * temp_term1_z + lambda_x0 * temp_term2_z;
    double temp_miu5 = temp_z_3 * temp_term1_y + lambda_z0 * temp_term2_y;
    double temp_miu6 = temp_y_3 * temp_term1_z + lambda_y0 * temp_term2_z;

    double temp_dmiu1_dx = temp_y_3 * temp_dterm1_dx + lambda_y0 * temp_dterm2_dx;
    double temp_dmiu1_dy = temp_dy_3 * temp_term1_x + lambda * temp_term2_x;

    double temp_dmiu2_dx = temp_dx_3 * temp_term1_y + lambda * temp_term2_y;
    double temp_dmiu2_dy = temp_x_3 * temp_dterm1_dy + lambda_x0 * temp_dterm2_dy;

    double temp_dmiu3_dx = temp_z_3 * temp_dterm1_dx + lambda_z0 * temp_dterm2_dx;
    double temp_dmiu3_dz = temp_dz_3 * temp_term1_x + lambda * temp_term2_x;

    double temp_dmiu4_dx = temp_dx_3 * temp_term1_z + lambda * temp_term2_z;
    double temp_dmiu4_dz = temp_x_3 * temp_dterm1_dz + lambda_x0 * temp_dterm2_dz;

    double temp_dmiu5_dy = temp_z_3 * temp_dterm1_dy + lambda_z0 * temp_dterm2_dy;
    double temp_dmiu5_dz = temp_dz_3 * temp_term1_y + lambda * temp_term2_y;

    double temp_dmiu6_dy = temp_dy_3 * temp_term1_z + lambda * temp_term2_z;
    double temp_dmiu6_dz = temp_y_3 * temp_dterm1_dz + lambda_y0 * temp_dterm2_dz;

    double miu_7_5_1 = temp * temp_miu1;
    double miu_7_5_2 = temp * temp_miu2;
    double miu_7_5_3 = temp * temp_miu3;
    double miu_7_5_4 = temp * temp_miu4;
    double miu_7_5_5 = temp * temp_miu5;
    double miu_7_5_6 = temp * temp_miu6;

    //deriv consts
    double const_2_C2_x = 2.0 * C2 * x0; // 2 * C2 * x0
    double const_2_C2_y = 2.0 * C2 * y0;
    double const_2_C2_z = 2.0 * C2 * z0;

    // dmiu1 dx/dy/dz
    deriv[0] = temp * (temp_dmiu1_dx + const_2_C2_x * temp_miu1);
    deriv[1] = temp * (temp_dmiu1_dy + const_2_C2_y * temp_miu1);
    deriv[2] = miu_7_5_1 * const_2_C2_z;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * (temp_dmiu2_dx + const_2_C2_x * temp_miu2);
    deriv[4] = temp * (temp_dmiu2_dy + const_2_C2_y * temp_miu2);
    deriv[5] = miu_7_5_2 * const_2_C2_z;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * (temp_dmiu3_dx + const_2_C2_x * temp_miu3);
    deriv[7] = miu_7_5_3 * const_2_C2_y;
    deriv[8] = temp * (temp_dmiu3_dz + const_2_C2_z * temp_miu3);

    // dmiu4 dx/dy/dz
    deriv[9] = temp * (temp_dmiu4_dx + const_2_C2_x * temp_miu4);
    deriv[10] = miu_7_5_4 * const_2_C2_y;
    deriv[11] = temp * (temp_dmiu4_dz + const_2_C2_z * temp_miu4);

    // dmiu5 dx/dy/dz
    deriv[12] = miu_7_5_5 * const_2_C2_x;
    deriv[13] = temp * (temp_dmiu5_dy + const_2_C2_y * temp_miu5);
    deriv[14] = temp * (temp_dmiu5_dz + const_2_C2_z * temp_miu5);

    // dmiu6 dx/dy/dz
    deriv[15] = miu_7_5_6 * const_2_C2_x;
    deriv[16] = temp * (temp_dmiu6_dy + const_2_C2_y * temp_miu6);
    deriv[17] = temp * (temp_dmiu6_dz + const_2_C2_z * temp_miu6);

    value[0] = miu_7_5_1;
    value[1] = miu_7_5_2;
    value[2] = miu_7_5_3;
    value[3] = miu_7_5_4;
    value[4] = miu_7_5_5;
    value[5] = miu_7_5_6;
}

void calc_MCSH_7_6(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double lambda_sqr = lambda * lambda;
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

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x_2 = lambda_x0_sqr + (1.0/(2.0*gamma));
    double temp_y_2 = lambda_y0_sqr + (1.0/(2.0*gamma));
    double temp_z_2 = lambda_z0_sqr + (1.0/(2.0*gamma));

    double temp_dx_2 = 2.0 * lambda_sqr * x0;
    double temp_dy_2 = 2.0 * lambda_sqr * y0;
    double temp_dz_2 = 2.0 * lambda_sqr * z0;

    double C4_1 = 3.0 / gamma, C4_2 = 3.0 / (4.0 * gamma * gamma);
    double temp_x_4 = lambda_x0_4 + C4_1 * lambda_x0_sqr + C4_2;
    double temp_y_4 = lambda_y0_4 + C4_1 * lambda_y0_sqr + C4_2;
    double temp_z_4 = lambda_z0_4 + C4_1 * lambda_z0_sqr + C4_2;

    double temp_dx_4 = lambda * (4.0 * lambda_x0_3 + 2.0 * C4_1 * lambda_x0);
    double temp_dy_4 = lambda * (4.0 * lambda_y0_3 + 2.0 * C4_1 * lambda_y0);
    double temp_dz_4 = lambda * (4.0 * lambda_z0_3 + 2.0 * C4_1 * lambda_z0);

    double temp_term1_x = 135135.0 * temp_x_4 - 62370.0 * temp_x_2 + 2835.0;
    double temp_term1_y = 135135.0 * temp_y_4 - 62370.0 * temp_y_2 + 2835.0;
    double temp_term1_z = 135135.0 * temp_z_4 - 62370.0 * temp_z_2 + 2835.0;

    double temp_dterm1_dx = 135135.0 * temp_dx_4 - 62370.0 * temp_dx_2;
    double temp_dterm1_dy = 135135.0 * temp_dy_4 - 62370.0 * temp_dy_2;
    double temp_dterm1_dz = 135135.0 * temp_dz_4 - 62370.0 * temp_dz_2;

    double temp_term2_x = -10395.0 * temp_x_4 + 5670.0 * temp_x_2 - 315.0;
    double temp_term2_y = -10395.0 * temp_y_4 + 5670.0 * temp_y_2 - 315.0;
    double temp_term2_z = -10395.0 * temp_z_4 + 5670.0 * temp_z_2 - 315.0;

    double temp_dterm2_dx = -10395.0 * temp_dx_4 + 5670.0 * temp_dx_2;
    double temp_dterm2_dy = -10395.0 * temp_dy_4 + 5670.0 * temp_dy_2;
    double temp_dterm2_dz = -10395.0 * temp_dz_4 + 5670.0 * temp_dz_2;

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

    double miu_7_6_1 = temp * lambda_z0 * temp_miu1;
    double miu_7_6_2 = temp * lambda_z0 * temp_miu2;
    double miu_7_6_3 = temp * lambda_y0 * temp_miu3;
    double miu_7_6_4 = temp * lambda_y0 * temp_miu4;
    double miu_7_6_5 = temp * lambda_x0 * temp_miu5;
    double miu_7_6_6 = temp * lambda_x0 * temp_miu6;

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

    // dmiu2 dx/dy/dz
    deriv[3] = temp * lambda_z0 * (temp_dmiu2_dx + const_2_C2_x * temp_miu2);
    deriv[4] = temp * lambda_z0 * (temp_dmiu2_dy + const_2_C2_y * temp_miu2);
    deriv[5] = temp * temp_miu2 * lambda * const_1_p_2_C2_z2;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * lambda_y0 * (temp_dmiu3_dx + const_2_C2_x * temp_miu3);
    deriv[7] = temp * temp_miu3 * lambda * const_1_p_2_C2_y2;
    deriv[8] = temp * lambda_y0 * (temp_dmiu3_dz + const_2_C2_z * temp_miu3);

    // dmiu4 dx/dy/dz
    deriv[9] = temp * lambda_y0 * (temp_dmiu4_dx + const_2_C2_x * temp_miu4);
    deriv[10] = temp * temp_miu4 * lambda * const_1_p_2_C2_y2;
    deriv[11] = temp * lambda_y0 * (temp_dmiu4_dz + const_2_C2_z * temp_miu4);

    // dmiu5 dx/dy/dz
    deriv[12] = temp * temp_miu5 * lambda * const_1_p_2_C2_x2;
    deriv[13] = temp * lambda_x0 * (temp_dmiu5_dy + const_2_C2_y * temp_miu5);
    deriv[14] = temp * lambda_x0 * (temp_dmiu5_dz + const_2_C2_z * temp_miu5);

    // dmiu6 dx/dy/dz
    deriv[15] = temp * temp_miu6 * lambda * const_1_p_2_C2_x2;
    deriv[16] = temp * lambda_x0 * (temp_dmiu6_dy + const_2_C2_y * temp_miu6);
    deriv[17] = temp * lambda_x0 * (temp_dmiu6_dz + const_2_C2_z * temp_miu6);

    value[0] = miu_7_6_1;
    value[1] = miu_7_6_2;
    value[2] = miu_7_6_3;
    value[3] = miu_7_6_4;
    value[4] = miu_7_6_5;
    value[5] = miu_7_6_6;
}


void calc_MCSH_7_7(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    // double lambda_sqr = lambda * lambda;

    double x0_sqr = x0 * x0;
    double y0_sqr = y0 * y0;
    double z0_sqr = z0 * z0;

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

    double temp = C1 * exp( C2 * r0_sqr);

    double C3_1 = 3.0 / (2.0 * gamma);
    double temp_x_3 = lambda_x0_3 + C3_1 * lambda_x0;
    double temp_y_3 = lambda_y0_3 + C3_1 * lambda_y0;
    double temp_z_3 = lambda_z0_3 + C3_1 * lambda_z0;

    double temp_dx_3 = lambda * (3.0 * lambda_x0_sqr + C3_1);
    double temp_dy_3 = lambda * (3.0 * lambda_y0_sqr + C3_1);
    double temp_dz_3 = lambda * (3.0 * lambda_z0_sqr + C3_1);


    double temp_term1_x = 135135.0 * temp_x_3 - 31185.0 * lambda_x0;
    double temp_term1_y = 135135.0 * temp_y_3 - 31185.0 * lambda_y0;

    double temp_dterm1_dx = 135135.0 * temp_dx_3 - 31185.0 * lambda;
    double temp_dterm1_dy = 135135.0 * temp_dy_3 - 31185.0 * lambda;

    double temp_term2_x = -31185.0 * temp_x_3 + 8505.0 * lambda_x0;
    double temp_term2_y = -31185.0 * temp_y_3 + 8505.0 * lambda_y0;

    double temp_dterm2_dx = -31185.0 * temp_dx_3 + 8505.0 * lambda;
    double temp_dterm2_dy = -31185.0 * temp_dy_3 + 8505.0 * lambda;

    double temp_miu1 = temp_y_3 * temp_term1_x + lambda_y0 * temp_term2_x;
    double temp_miu2 = temp_z_3 * temp_term1_x + lambda_z0 * temp_term2_x;
    double temp_miu3 = temp_z_3 * temp_term1_y + lambda_z0 * temp_term2_y;

    double temp_dmiu1_dx = temp_y_3 * temp_dterm1_dx + lambda_y0 * temp_dterm2_dx;
    double temp_dmiu1_dy = temp_dy_3 * temp_term1_x + lambda * temp_term2_x;

    double temp_dmiu2_dx = temp_z_3 * temp_dterm1_dx + lambda_z0 * temp_dterm2_dx;
    double temp_dmiu2_dz = temp_dz_3 * temp_term1_x + lambda * temp_term2_x;

    double temp_dmiu3_dy = temp_z_3 * temp_dterm1_dy + lambda_z0 * temp_dterm2_dy;
    double temp_dmiu3_dz = temp_dz_3 * temp_term1_y + lambda * temp_term2_y;

    double miu_7_7_1 = temp * lambda_z0 * temp_miu1;
    double miu_7_7_2 = temp * lambda_y0 * temp_miu2;
    double miu_7_7_3 = temp * lambda_x0 * temp_miu3;

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

    value[0] = miu_7_7_1;
    value[1] = miu_7_7_2;
    value[2] = miu_7_7_3;
}

void calc_MCSH_7_8(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double lambda_sqr = lambda * lambda;

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

    double temp_term1_x = 135135.0 * temp_x_3 - 31185.0 * lambda_x0;
    double temp_term1_y = 135135.0 * temp_y_3 - 31185.0 * lambda_y0;
    double temp_term1_z = 135135.0 * temp_z_3 - 31185.0 * lambda_z0;

    double temp_dterm1_dx = 135135.0 * temp_dx_3 - 31185.0 * lambda;
    double temp_dterm1_dy = 135135.0 * temp_dy_3 - 31185.0 * lambda;
    double temp_dterm1_dz = 135135.0 * temp_dz_3 - 31185.0 * lambda;

    double temp_term2_x = -10395.0 * temp_x_3 + 2835.0 * lambda_x0;
    double temp_term2_y = -10395.0 * temp_y_3 + 2835.0 * lambda_y0;
    double temp_term2_z = -10395.0 * temp_z_3 + 2835.0 * lambda_z0;

    double temp_dterm2_dx = -10395.0 * temp_dx_3 + 2835.0 * lambda;
    double temp_dterm2_dy = -10395.0 * temp_dy_3 + 2835.0 * lambda;
    double temp_dterm2_dz = -10395.0 * temp_dz_3 + 2835.0 * lambda;

    double temp_term3_x = 945.0 * temp_x_3 - 315.0 * lambda_x0;
    double temp_term3_y = 945.0 * temp_y_3 - 315.0 * lambda_y0;
    double temp_term3_z = 945.0 * temp_z_3 - 315.0 * lambda_z0;

    double temp_dterm3_dx = 945.0 * temp_dx_3 - 315.0 * lambda;
    double temp_dterm3_dy = 945.0 * temp_dy_3 - 315.0 * lambda;
    double temp_dterm3_dz = 945.0 * temp_dz_3 - 315.0 * lambda;

    double temp_miu1 = temp_y_2 * temp_z_2 * temp_term1_x + (temp_y_2 + temp_z_2) * temp_term2_x + temp_term3_x;
    double temp_miu2 = temp_x_2 * temp_z_2 * temp_term1_y + (temp_x_2 + temp_z_2) * temp_term2_y + temp_term3_y;
    double temp_miu3 = temp_x_2 * temp_y_2 * temp_term1_z + (temp_x_2 + temp_y_2) * temp_term2_z + temp_term3_z;

    double temp_dmiu1_dx = temp_y_2 * temp_z_2 * temp_dterm1_dx + (temp_y_2 + temp_z_2) * temp_dterm2_dx + temp_dterm3_dx;
    double temp_dmiu1_dy = temp_dy_2 * temp_z_2 * temp_term1_x + temp_dy_2 * temp_term2_x;
    double temp_dmiu1_dz = temp_y_2 * temp_dz_2 * temp_term1_x + temp_dz_2 * temp_term2_x;

    double temp_dmiu2_dx = temp_dx_2 * temp_z_2 * temp_term1_y + temp_dx_2 * temp_term2_y;
    double temp_dmiu2_dy = temp_x_2 * temp_z_2 * temp_dterm1_dy + (temp_x_2 + temp_z_2) * temp_dterm2_dy + temp_dterm3_dy;
    double temp_dmiu2_dz = temp_x_2 * temp_dz_2 * temp_term1_y + temp_dz_2 * temp_term2_y;

    double temp_dmiu3_dx = temp_dx_2 * temp_y_2 * temp_term1_z + temp_dx_2 * temp_term2_z;
    double temp_dmiu3_dy = temp_x_2 * temp_dy_2 * temp_term1_z + temp_dy_2 * temp_term2_z;
    double temp_dmiu3_dz = temp_x_2 * temp_y_2 * temp_dterm1_dz + (temp_x_2 + temp_y_2) * temp_dterm2_dz + temp_dterm3_dz;


    double miu_7_8_1 = temp * temp_miu1;
    double miu_7_8_2 = temp * temp_miu2;
    double miu_7_8_3 = temp * temp_miu3;

    //deriv consts
    double const_2_C2_x = 2.0 * C2 * x0; // 2 * C2 * x0
    double const_2_C2_y = 2.0 * C2 * y0;
    double const_2_C2_z = 2.0 * C2 * z0;

    // dmiu1 dx/dy/dz
    deriv[0] = temp * (temp_dmiu1_dx + const_2_C2_x * temp_miu1);
    deriv[1] = temp * (temp_dmiu1_dy + const_2_C2_y * temp_miu1);
    deriv[2] = temp * (temp_dmiu1_dz + const_2_C2_z * temp_miu1);

    // dmiu2 dx/dy/dz
    deriv[3] = temp * (temp_dmiu2_dx + const_2_C2_x * temp_miu2);
    deriv[4] = temp * (temp_dmiu2_dy + const_2_C2_y * temp_miu2);
    deriv[5] = temp * (temp_dmiu2_dz + const_2_C2_z * temp_miu2);

    // dmiu3 dx/dy/dz
    deriv[6] = temp * (temp_dmiu3_dx + const_2_C2_x * temp_miu3);
    deriv[7] = temp * (temp_dmiu3_dy + const_2_C2_y * temp_miu3);
    deriv[8] = temp * (temp_dmiu3_dz + const_2_C2_z * temp_miu3);

    value[0] = miu_7_8_1;
    value[1] = miu_7_8_2;
    value[2] = miu_7_8_3;
}

void calc_MCSH_8_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    // double lambda_sqr = lambda * lambda;
    // double x0_sqr = x0*x0;
    // double y0_sqr = y0*y0;
    // double z0_sqr = z0*z0;

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

    double lambda_x0_6 = lambda_x0_5 * lambda_x0;
    double lambda_y0_6 = lambda_y0_5 * lambda_y0;
    double lambda_z0_6 = lambda_z0_5 * lambda_z0;

    double lambda_x0_7 = lambda_x0_6 * lambda_x0;
    double lambda_y0_7 = lambda_y0_6 * lambda_y0;
    double lambda_z0_7 = lambda_z0_6 * lambda_z0;

    double lambda_x0_8 = lambda_x0_7 * lambda_x0;
    double lambda_y0_8 = lambda_y0_7 * lambda_y0;
    double lambda_z0_8 = lambda_z0_7 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = (2027025.0 * 14.0 / gamma) - 3783780.0;
    double C4 = (2027025.0 * 105.0 / (2.0 * gamma * gamma)) - (3783780.0 * 15.0 / (2.0 * gamma)) + 2182950.0;
    double C5 = (2027025.0 * 105.0 / (2.0 * gamma * gamma * gamma)) - (3783780.0 * 45.0 / (4.0 * gamma * gamma)) + (2182950.0 * 3.0 / gamma) - 396900.0;
    double C6 = (2027025.0 * 105.0 / (16.0 * gamma * gamma * gamma * gamma)) - (3783780.0 * 15.0 / (8.0 * gamma * gamma * gamma)) + (2182950.0 * 3.0 / (4.0 * gamma *gamma)) - (396900.0 / (2.0 * gamma)) + 11025.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 2027025.0 * lambda_x0_8 + C3 * lambda_x0_6 + C4 * lambda_x0_4 + C5 * lambda_x0_sqr + C6;
    double temp_y = 2027025.0 * lambda_y0_8 + C3 * lambda_y0_6 + C4 * lambda_y0_4 + C5 * lambda_y0_sqr + C6;
    double temp_z = 2027025.0 * lambda_z0_8 + C3 * lambda_z0_6 + C4 * lambda_z0_4 + C5 * lambda_z0_sqr + C6;

    double temp_dx = lambda * (2027025.0 * 8.0 * lambda_x0_7 + 6.0 * C3 * lambda_x0_5 + 4.0 * C4 * lambda_x0_3 + 2.0 * C5 * lambda_x0);
    double temp_dy = lambda * (2027025.0 * 8.0 * lambda_y0_7 + 6.0 * C3 * lambda_y0_5 + 4.0 * C4 * lambda_y0_3 + 2.0 * C5 * lambda_y0);
    double temp_dz = lambda * (2027025.0 * 8.0 * lambda_z0_7 + 6.0 * C3 * lambda_z0_5 + 4.0 * C4 * lambda_z0_3 + 2.0 * C5 * lambda_z0);

    double miu_8_1_1 = temp * temp_x;
    double miu_8_1_2 = temp * temp_y;
    double miu_8_1_3 = temp * temp_z;

    //deriv consts
    double const_2_C2_x = 2.0 * C2 * x0;
    double const_2_C2_y = 2.0 * C2 * y0;
    double const_2_C2_z = 2.0 * C2 * z0;

    // dmiu1 dx/dy/dz
    deriv[0] = temp * (temp_dx + const_2_C2_x * temp_x);
    deriv[1] = miu_8_1_1 * const_2_C2_y;
    deriv[2] = miu_8_1_1 * const_2_C2_z;

    // dmiu2 dx/dy/dz
    deriv[3] = miu_8_1_2 * const_2_C2_x;
    deriv[4] = temp * (temp_dy + const_2_C2_y * temp_y);
    deriv[5] = miu_8_1_2 * const_2_C2_z;

    // dmiu3 dx/dy/dz
    deriv[6] = miu_8_1_3 * const_2_C2_x;
    deriv[7] = miu_8_1_3 * const_2_C2_y;
    deriv[8] = temp * (temp_dz + const_2_C2_z * temp_z);

    value[0] = miu_8_1_1;
    value[1] = miu_8_1_2;
    value[2] = miu_8_1_3;
}

void calc_MCSH_8_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
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

    double lambda_x0_5 = lambda_x0_4 * lambda_x0;
    double lambda_y0_5 = lambda_y0_4 * lambda_y0;
    double lambda_z0_5 = lambda_z0_4 * lambda_z0;

    double lambda_x0_6 = lambda_x0_5 * lambda_x0;
    double lambda_y0_6 = lambda_y0_5 * lambda_y0;
    double lambda_z0_6 = lambda_z0_5 * lambda_z0;

    double lambda_x0_7 = lambda_x0_6 * lambda_x0;
    double lambda_y0_7 = lambda_y0_6 * lambda_y0;
    double lambda_z0_7 = lambda_z0_6 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = (2027025.0 * 21.0 / (2.0 * gamma)) - 2837835.0;
    double C4 = (2027025.0 * 105.0 / (4.0 * gamma * gamma)) - (2837835.0 * 5.0 / gamma) + 1091475.0;
    double C5 = (2027025.0 * 105.0 / (8.0 * gamma * gamma * gamma)) - (2837835.0 * 15.0 / (4.0 * gamma * gamma)) + (1091475.0 * 3.0 / (2.0 * gamma)) - 99225.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 2027025.0 * lambda_x0_7 + C3 * lambda_x0_5 + C4 * lambda_x0_3 + C5 * lambda_x0;
    double temp_y = 2027025.0 * lambda_y0_7 + C3 * lambda_y0_5 + C4 * lambda_y0_3 + C5 * lambda_y0;
    double temp_z = 2027025.0 * lambda_z0_7 + C3 * lambda_z0_5 + C4 * lambda_z0_3 + C5 * lambda_z0;

    double temp_dx = lambda * (2027025.0 * 7.0 * lambda_x0_6 + 5.0 * C3 * lambda_x0_4 + 3.0 * C4 * lambda_x0_sqr + C5);
    double temp_dy = lambda * (2027025.0 * 7.0 * lambda_y0_6 + 5.0 * C3 * lambda_y0_4 + 3.0 * C4 * lambda_y0_sqr + C5);
    double temp_dz = lambda * (2027025.0 * 7.0 * lambda_z0_6 + 5.0 * C3 * lambda_z0_4 + 3.0 * C4 * lambda_z0_sqr + C5);

    double miu_8_2_1 = temp * lambda_y0 * temp_x;
    double miu_8_2_2 = temp * lambda_x0 * temp_y;
    double miu_8_2_3 = temp * lambda_z0 * temp_x;
    double miu_8_2_4 = temp * lambda_x0 * temp_z;
    double miu_8_2_5 = temp * lambda_z0 * temp_y;
    double miu_8_2_6 = temp * lambda_y0 * temp_z;

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
    deriv[2] = miu_8_2_1 * const_2_C2_z;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * lambda * temp_y * const_1_p_2_C2_x2;
    deriv[4] = temp * lambda_x0 * (temp_dy + const_2_C2_y * temp_y);
    deriv[5] = miu_8_2_2 * const_2_C2_z;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * lambda_z0 * (temp_dx + const_2_C2_x * temp_x);
    deriv[7] = miu_8_2_3 * const_2_C2_y;
    deriv[8] = temp * lambda * temp_x * const_1_p_2_C2_z2;

    // dmiu4 dx/dy/dz
    deriv[9] = temp * lambda * temp_z * const_1_p_2_C2_x2;
    deriv[10] = miu_8_2_4 * const_2_C2_y;
    deriv[11] = temp * lambda_x0 * (temp_dz + const_2_C2_z * temp_z);

    // dmiu5 dx/dy/dz
    deriv[12] = miu_8_2_5 * const_2_C2_x;
    deriv[13] = temp * lambda_z0 * (temp_dy + const_2_C2_y * temp_y);
    deriv[14] = temp * lambda * temp_y * const_1_p_2_C2_z2;

    // dmiu6 dx/dy/dz
    deriv[15] = miu_8_2_6 * const_2_C2_x;
    deriv[16] = temp * lambda * temp_z * const_1_p_2_C2_y2;
    deriv[17] = temp * lambda_y0 * (temp_dz + const_2_C2_z * temp_z);

    value[0] = miu_8_2_1;
    value[1] = miu_8_2_2;
    value[2] = miu_8_2_3;
    value[3] = miu_8_2_4;
    value[4] = miu_8_2_5;
    value[5] = miu_8_2_6;
}


void calc_MCSH_8_3(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double lambda_sqr = lambda * lambda;

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

    double lambda_x0_5 = lambda_x0_4 * lambda_x0;
    double lambda_y0_5 = lambda_y0_4 * lambda_y0;
    double lambda_z0_5 = lambda_z0_4 * lambda_z0;

    double lambda_x0_6 = lambda_x0_5 * lambda_x0;
    double lambda_y0_6 = lambda_y0_5 * lambda_y0;
    double lambda_z0_6 = lambda_z0_5 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x_2 = lambda_x0_sqr + (1.0/(2.0*gamma));
    double temp_y_2 = lambda_y0_sqr + (1.0/(2.0*gamma));
    double temp_z_2 = lambda_z0_sqr + (1.0/(2.0*gamma));

    double temp_dx_2 = 2.0 * lambda_sqr * x0;
    double temp_dy_2 = 2.0 * lambda_sqr * y0;
    double temp_dz_2 = 2.0 * lambda_sqr * z0;

    double C4_1 = 3.0 / gamma, C4_2 = 3.0 / (4.0 * gamma * gamma);
    double temp_x_4 = lambda_x0_4 + C4_1 * lambda_x0_sqr + C4_2;
    double temp_y_4 = lambda_y0_4 + C4_1 * lambda_y0_sqr + C4_2;
    double temp_z_4 = lambda_z0_4 + C4_1 * lambda_z0_sqr + C4_2;

    double temp_dx_4 = lambda * (4.0 * lambda_x0_3 + 2.0 * C4_1 * lambda_x0);
    double temp_dy_4 = lambda * (4.0 * lambda_y0_3 + 2.0 * C4_1 * lambda_y0);
    double temp_dz_4 = lambda * (4.0 * lambda_z0_3 + 2.0 * C4_1 * lambda_z0);

    double C6_1 = 15.0 / (2.0 * gamma), C6_2 = 45.0 / (4.0 * gamma * gamma), C6_3 =  15.0 / (8.0 * gamma * gamma * gamma);
    double temp_x_6 = lambda_x0_6 + C6_1 * lambda_x0_4 + C6_2 * lambda_x0_sqr + C6_3;
    double temp_y_6 = lambda_y0_6 + C6_1 * lambda_y0_4 + C6_2 * lambda_y0_sqr + C6_3;
    double temp_z_6 = lambda_z0_6 + C6_1 * lambda_z0_4 + C6_2 * lambda_z0_sqr + C6_3;

    double temp_dx_6 = lambda * (6.0 * lambda_x0_5 + 4.0 * C6_1 * lambda_x0_3 + 2.0 * C6_2 * lambda_x0);
    double temp_dy_6 = lambda * (6.0 * lambda_y0_5 + 4.0 * C6_1 * lambda_y0_3 + 2.0 * C6_2 * lambda_y0);
    double temp_dz_6 = lambda * (6.0 * lambda_z0_5 + 4.0 * C6_1 * lambda_z0_3 + 2.0 * C6_2 * lambda_z0);

    double temp_term1_x = 2027025.0 * temp_x_6 - 2027025.0 * temp_x_4 + 467775.0 * temp_x_2 - 14175.0;
    double temp_term1_y = 2027025.0 * temp_y_6 - 2027025.0 * temp_y_4 + 467775.0 * temp_y_2 - 14175.0;
    double temp_term1_z = 2027025.0 * temp_z_6 - 2027025.0 * temp_z_4 + 467775.0 * temp_z_2 - 14175.0;

    double temp_dterm1_dx = 2027025.0 * temp_dx_6 - 2027025.0 * temp_dx_4 + 467775.0 * temp_dx_2;
    double temp_dterm1_dy = 2027025.0 * temp_dy_6 - 2027025.0 * temp_dy_4 + 467775.0 * temp_dy_2;
    double temp_dterm1_dz = 2027025.0 * temp_dz_6 - 2027025.0 * temp_dz_4 + 467775.0 * temp_dz_2;

    double temp_term2_x = -135135.0 * temp_x_6 + 155925.0 * temp_x_4 - 42525.0 * temp_x_2 + 1575.0;
    double temp_term2_y = -135135.0 * temp_y_6 + 155925.0 * temp_y_4 - 42525.0 * temp_y_2 + 1575.0;
    double temp_term2_z = -135135.0 * temp_z_6 + 155925.0 * temp_z_4 - 42525.0 * temp_z_2 + 1575.0;

    double temp_dterm2_dx = -135135.0 * temp_dx_6 + 155925.0 * temp_dx_4 - 42525.0 * temp_dx_2;
    double temp_dterm2_dy = -135135.0 * temp_dy_6 + 155925.0 * temp_dy_4 - 42525.0 * temp_dy_2;
    double temp_dterm2_dz = -135135.0 * temp_dz_6 + 155925.0 * temp_dz_4 - 42525.0 * temp_dz_2;

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

    double miu_8_3_1 = temp * temp_miu1;
    double miu_8_3_2 = temp * temp_miu2;
    double miu_8_3_3 = temp * temp_miu3;
    double miu_8_3_4 = temp * temp_miu4;
    double miu_8_3_5 = temp * temp_miu5;
    double miu_8_3_6 = temp * temp_miu6;

    //deriv consts
    double const_2_C2_x = 2.0 * C2 * x0; // 2 * C2 * x0
    double const_2_C2_y = 2.0 * C2 * y0;
    double const_2_C2_z = 2.0 * C2 * z0;

    // dmiu1 dx/dy/dz
    deriv[0] = temp * (temp_dmiu1_dx + const_2_C2_x * temp_miu1);
    deriv[1] = temp * (temp_dmiu1_dy + const_2_C2_y * temp_miu1);
    deriv[2] = miu_8_3_1 * const_2_C2_z;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * (temp_dmiu2_dx + const_2_C2_x * temp_miu2);
    deriv[4] = temp * (temp_dmiu2_dy + const_2_C2_y * temp_miu2);
    deriv[5] = miu_8_3_2 * const_2_C2_z;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * (temp_dmiu3_dx + const_2_C2_x * temp_miu3);
    deriv[7] = miu_8_3_3 * const_2_C2_y;
    deriv[8] = temp * (temp_dmiu3_dz + const_2_C2_z * temp_miu3);

    // dmiu4 dx/dy/dz
    deriv[9] = temp * (temp_dmiu4_dx + const_2_C2_x * temp_miu4);
    deriv[10] = miu_8_3_4 * const_2_C2_y;
    deriv[11] = temp * (temp_dmiu4_dz + const_2_C2_z * temp_miu4);

    // dmiu5 dx/dy/dz
    deriv[12] = miu_8_3_5 * const_2_C2_x;
    deriv[13] = temp * (temp_dmiu5_dy + const_2_C2_y * temp_miu5);
    deriv[14] = temp * (temp_dmiu5_dz + const_2_C2_z * temp_miu5);

    // dmiu6 dx/dy/dz
    deriv[15] = miu_8_3_6 * const_2_C2_x;
    deriv[16] = temp * (temp_dmiu6_dy + const_2_C2_y * temp_miu6);
    deriv[17] = temp * (temp_dmiu6_dz + const_2_C2_z * temp_miu6);

    value[0] = miu_8_3_1;
    value[1] = miu_8_3_2;
    value[2] = miu_8_3_3;
    value[3] = miu_8_3_4;
    value[4] = miu_8_3_5;
    value[5] = miu_8_3_6;
}


void calc_MCSH_8_4(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
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

    double lambda_x0_5 = lambda_x0_4 * lambda_x0;
    double lambda_y0_5 = lambda_y0_4 * lambda_y0;
    double lambda_z0_5 = lambda_z0_4 * lambda_z0;

    double lambda_x0_6 = lambda_x0_5 * lambda_x0;
    double lambda_y0_6 = lambda_y0_5 * lambda_y0;
    double lambda_z0_6 = lambda_z0_5 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = (2027025.0 * 15.0 / (2.0 * gamma)) - 2027025.0;
    double C4 = (2027025.0 * 45.0 / (4.0 * gamma * gamma)) - (2027025.0 * 3.0 / gamma) + 467775.0;
    double C5 = (2027025.0 * 15.0 / (8.0 * gamma * gamma * gamma)) - (2027025.0 * 3.0 / (4.0 * gamma * gamma)) + (467775.0 / (2.0 * gamma)) - 14175.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 2027025.0 * lambda_x0_6 + C3 * lambda_x0_4 + C4 * lambda_x0_sqr + C5;
    double temp_y = 2027025.0 * lambda_y0_6 + C3 * lambda_y0_4 + C4 * lambda_y0_sqr + C5;
    double temp_z = 2027025.0 * lambda_z0_6 + C3 * lambda_z0_4 + C4 * lambda_z0_sqr + C5;

    double temp_dx = lambda * (2027025.0 * 6.0 * lambda_x0_5 + 4.0 * C3 * lambda_x0_3 + 2.0 * C4 * lambda_x0);
    double temp_dy = lambda * (2027025.0 * 6.0 * lambda_y0_5 + 4.0 * C3 * lambda_y0_3 + 2.0 * C4 * lambda_y0);
    double temp_dz = lambda * (2027025.0 * 6.0 * lambda_z0_5 + 4.0 * C3 * lambda_z0_3 + 2.0 * C4 * lambda_z0);

    double miu_8_4_1 = temp * lambda_y0 * lambda_z0 * temp_x;
    double miu_8_4_2 = temp * lambda_x0 * lambda_z0 * temp_y;
    double miu_8_4_3 = temp * lambda_x0 * lambda_y0 * temp_z;

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

    value[0] = miu_8_4_1;
    value[1] = miu_8_4_2;
    value[2] = miu_8_4_3;
}

void calc_MCSH_8_5(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    // double lambda_sqr = lambda * lambda;

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

    double lambda_x0_5 = lambda_x0_4 * lambda_x0;
    double lambda_y0_5 = lambda_y0_4 * lambda_y0;
    double lambda_z0_5 = lambda_z0_4 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);

    double temp = C1 * exp( C2 * r0_sqr);

    double C3_1 = 3.0 / (2.0 * gamma);
    double temp_x_3 = lambda_x0_3 + C3_1 * lambda_x0;
    double temp_y_3 = lambda_y0_3 + C3_1 * lambda_y0;
    double temp_z_3 = lambda_z0_3 + C3_1 * lambda_z0;

    double temp_dx_3 = lambda * (3.0 * lambda_x0_sqr + C3_1);
    double temp_dy_3 = lambda * (3.0 * lambda_y0_sqr + C3_1);
    double temp_dz_3 = lambda * (3.0 * lambda_z0_sqr + C3_1);

    double C5_1 = 5.0 / gamma, C5_2 = 15.0 / (4.0 * gamma * gamma);
    double temp_x_5 = lambda_x0_5 + C5_1 * lambda_x0_3 + C5_2 * lambda_x0;
    double temp_y_5 = lambda_y0_5 + C5_1 * lambda_y0_3 + C5_2 * lambda_y0;
    double temp_z_5 = lambda_z0_5 + C5_1 * lambda_z0_3 + C5_2 * lambda_z0;

    double temp_dx_5 = lambda * (5.0 * lambda_x0_4 + 3.0 * C5_1 * lambda_x0_sqr + C5_2);
    double temp_dy_5 = lambda * (5.0 * lambda_y0_4 + 3.0 * C5_1 * lambda_y0_sqr + C5_2);
    double temp_dz_5 = lambda * (5.0 * lambda_z0_4 + 3.0 * C5_1 * lambda_z0_sqr + C5_2);



    double temp_term1_x = 2027025.0 * temp_x_5 - 1351350.0 * temp_x_3 + 155925.0 * lambda_x0;
    double temp_term1_y = 2027025.0 * temp_y_5 - 1351350.0 * temp_y_3 + 155925.0 * lambda_y0;
    double temp_term1_z = 2027025.0 * temp_z_5 - 1351350.0 * temp_z_3 + 155925.0 * lambda_z0;

    double temp_dterm1_dx = 2027025.0 * temp_dx_5 - 1351350.0 * temp_dx_3 + 155925.0 * lambda;
    double temp_dterm1_dy = 2027025.0 * temp_dy_5 - 1351350.0 * temp_dy_3 + 155925.0 * lambda;
    double temp_dterm1_dz = 2027025.0 * temp_dz_5 - 1351350.0 * temp_dz_3 + 155925.0 * lambda;

    double temp_term2_x = -405405.0 * temp_x_5 + 311850.0 * temp_x_3 - 42525.0 * lambda_x0;
    double temp_term2_y = -405405.0 * temp_y_5 + 311850.0 * temp_y_3 - 42525.0 * lambda_y0;
    double temp_term2_z = -405405.0 * temp_z_5 + 311850.0 * temp_z_3 - 42525.0 * lambda_z0;

    double temp_dterm2_dx = -405405.0 * temp_dx_5 + 311850.0 * temp_dx_3 - 42525.0 * lambda;
    double temp_dterm2_dy = -405405.0 * temp_dy_5 + 311850.0 * temp_dy_3 - 42525.0 * lambda;
    double temp_dterm2_dz = -405405.0 * temp_dz_5 + 311850.0 * temp_dz_3 - 42525.0 * lambda;

    double temp_miu1 = temp_y_3 * temp_term1_x + lambda_y0 * temp_term2_x;
    double temp_miu2 = temp_x_3 * temp_term1_y + lambda_x0 * temp_term2_y;
    double temp_miu3 = temp_z_3 * temp_term1_x + lambda_z0 * temp_term2_x;
    double temp_miu4 = temp_x_3 * temp_term1_z + lambda_x0 * temp_term2_z;
    double temp_miu5 = temp_z_3 * temp_term1_y + lambda_z0 * temp_term2_y;
    double temp_miu6 = temp_y_3 * temp_term1_z + lambda_y0 * temp_term2_z;

    double temp_dmiu1_dx = temp_y_3 * temp_dterm1_dx + lambda_y0 * temp_dterm2_dx;
    double temp_dmiu1_dy = temp_dy_3 * temp_term1_x + lambda * temp_term2_x;

    double temp_dmiu2_dx = temp_dx_3 * temp_term1_y + lambda * temp_term2_y;
    double temp_dmiu2_dy = temp_x_3 * temp_dterm1_dy + lambda_x0 * temp_dterm2_dy;

    double temp_dmiu3_dx = temp_z_3 * temp_dterm1_dx + lambda_z0 * temp_dterm2_dx;
    double temp_dmiu3_dz = temp_dz_3 * temp_term1_x + lambda * temp_term2_x;

    double temp_dmiu4_dx = temp_dx_3 * temp_term1_z + lambda * temp_term2_z;
    double temp_dmiu4_dz = temp_x_3 * temp_dterm1_dz + lambda_x0 * temp_dterm2_dz;

    double temp_dmiu5_dy = temp_z_3 * temp_dterm1_dy + lambda_z0 * temp_dterm2_dy;
    double temp_dmiu5_dz = temp_dz_3 * temp_term1_y + lambda * temp_term2_y;

    double temp_dmiu6_dy = temp_dy_3 * temp_term1_z + lambda * temp_term2_z;
    double temp_dmiu6_dz = temp_y_3 * temp_dterm1_dz + lambda_y0 * temp_dterm2_dz;

    double miu_8_5_1 = temp * temp_miu1;
    double miu_8_5_2 = temp * temp_miu2;
    double miu_8_5_3 = temp * temp_miu3;
    double miu_8_5_4 = temp * temp_miu4;
    double miu_8_5_5 = temp * temp_miu5;
    double miu_8_5_6 = temp * temp_miu6;

    //deriv consts
    double const_2_C2_x = 2.0 * C2 * x0; // 2 * C2 * x0
    double const_2_C2_y = 2.0 * C2 * y0;
    double const_2_C2_z = 2.0 * C2 * z0;

    // dmiu1 dx/dy/dz
    deriv[0] = temp * (temp_dmiu1_dx + const_2_C2_x * temp_miu1);
    deriv[1] = temp * (temp_dmiu1_dy + const_2_C2_y * temp_miu1);
    deriv[2] = miu_8_5_1 * const_2_C2_z;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * (temp_dmiu2_dx + const_2_C2_x * temp_miu2);
    deriv[4] = temp * (temp_dmiu2_dy + const_2_C2_y * temp_miu2);
    deriv[5] = miu_8_5_2 * const_2_C2_z;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * (temp_dmiu3_dx + const_2_C2_x * temp_miu3);
    deriv[7] = miu_8_5_3 * const_2_C2_y;
    deriv[8] = temp * (temp_dmiu3_dz + const_2_C2_z * temp_miu3);

    // dmiu4 dx/dy/dz
    deriv[9] = temp * (temp_dmiu4_dx + const_2_C2_x * temp_miu4);
    deriv[10] = miu_8_5_4 * const_2_C2_y;
    deriv[11] = temp * (temp_dmiu4_dz + const_2_C2_z * temp_miu4);

    // dmiu5 dx/dy/dz
    deriv[12] = miu_8_5_5 * const_2_C2_x;
    deriv[13] = temp * (temp_dmiu5_dy + const_2_C2_y * temp_miu5);
    deriv[14] = temp * (temp_dmiu5_dz + const_2_C2_z * temp_miu5);

    // dmiu6 dx/dy/dz
    deriv[15] = miu_8_5_6 * const_2_C2_x;
    deriv[16] = temp * (temp_dmiu6_dy + const_2_C2_y * temp_miu6);
    deriv[17] = temp * (temp_dmiu6_dz + const_2_C2_z * temp_miu6);

    value[0] = miu_8_5_1;
    value[1] = miu_8_5_2;
    value[2] = miu_8_5_3;
    value[3] = miu_8_5_4;
    value[4] = miu_8_5_5;
    value[5] = miu_8_5_6;
}

void calc_MCSH_8_6(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double lambda_sqr = lambda * lambda;
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

    double lambda_x0_5 = lambda_x0_4 * lambda_x0;
    double lambda_y0_5 = lambda_y0_4 * lambda_y0;
    double lambda_z0_5 = lambda_z0_4 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);

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

    double C5_1 = 5.0 / gamma, C5_2 = 15.0 / (4.0 * gamma * gamma);
    double temp_x_5 = lambda_x0_5 + C5_1 * lambda_x0_3 + C5_2 * lambda_x0;
    double temp_y_5 = lambda_y0_5 + C5_1 * lambda_y0_3 + C5_2 * lambda_y0;
    double temp_z_5 = lambda_z0_5 + C5_1 * lambda_z0_3 + C5_2 * lambda_z0;

    double temp_dx_5 = lambda * (5.0 * lambda_x0_4 + 3.0 * C5_1 * lambda_x0_sqr + C5_2);
    double temp_dy_5 = lambda * (5.0 * lambda_y0_4 + 3.0 * C5_1 * lambda_y0_sqr + C5_2);
    double temp_dz_5 = lambda * (5.0 * lambda_z0_4 + 3.0 * C5_1 * lambda_z0_sqr + C5_2);



    double temp_term1_x = 2027025.0 * temp_x_5 - 1351350.0 * temp_x_3 + 155925.0 * lambda_x0;
    double temp_term1_y = 2027025.0 * temp_y_5 - 1351350.0 * temp_y_3 + 155925.0 * lambda_y0;
    double temp_term1_z = 2027025.0 * temp_z_5 - 1351350.0 * temp_z_3 + 155925.0 * lambda_z0;

    double temp_dterm1_dx = 2027025.0 * temp_dx_5 - 1351350.0 * temp_dx_3 + 155925.0 * lambda;
    double temp_dterm1_dy = 2027025.0 * temp_dy_5 - 1351350.0 * temp_dy_3 + 155925.0 * lambda;
    double temp_dterm1_dz = 2027025.0 * temp_dz_5 - 1351350.0 * temp_dz_3 + 155925.0 * lambda;

    double temp_term2_x = -135135.0 * temp_x_5 + 103950.0 * temp_x_3 - 14175.0 * lambda_x0;
    double temp_term2_y = -135135.0 * temp_y_5 + 103950.0 * temp_y_3 - 14175.0 * lambda_y0;
    double temp_term2_z = -135135.0 * temp_z_5 + 103950.0 * temp_z_3 - 14175.0 * lambda_z0;

    double temp_dterm2_dx = -135135.0 * temp_dx_5 + 103950.0 * temp_dx_3 - 14175.0 * lambda;
    double temp_dterm2_dy = -135135.0 * temp_dy_5 + 103950.0 * temp_dy_3 - 14175.0 * lambda;
    double temp_dterm2_dz = -135135.0 * temp_dz_5 + 103950.0 * temp_dz_3 - 14175.0 * lambda;


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

    double miu_8_6_1 = temp * lambda_z0 * temp_miu1;
    double miu_8_6_2 = temp * lambda_z0 * temp_miu2;
    double miu_8_6_3 = temp * lambda_y0 * temp_miu3;
    double miu_8_6_4 = temp * lambda_y0 * temp_miu4;
    double miu_8_6_5 = temp * lambda_x0 * temp_miu5;
    double miu_8_6_6 = temp * lambda_x0 * temp_miu6;

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

    // dmiu2 dx/dy/dz
    deriv[3] = temp * lambda_z0 * (temp_dmiu2_dx + const_2_C2_x * temp_miu2);
    deriv[4] = temp * lambda_z0 * (temp_dmiu2_dy + const_2_C2_y * temp_miu2);
    deriv[5] = temp * temp_miu2 * lambda * const_1_p_2_C2_z2;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * lambda_y0 * (temp_dmiu3_dx + const_2_C2_x * temp_miu3);
    deriv[7] = temp * temp_miu3 * lambda * const_1_p_2_C2_y2;
    deriv[8] = temp * lambda_y0 * (temp_dmiu3_dz + const_2_C2_z * temp_miu3);

    // dmiu4 dx/dy/dz
    deriv[9] = temp * lambda_y0 * (temp_dmiu4_dx + const_2_C2_x * temp_miu4);
    deriv[10] = temp * temp_miu4 * lambda * const_1_p_2_C2_y2;
    deriv[11] = temp * lambda_y0 * (temp_dmiu4_dz + const_2_C2_z * temp_miu4);

    // dmiu5 dx/dy/dz
    deriv[12] = temp * temp_miu5 * lambda * const_1_p_2_C2_x2;
    deriv[13] = temp * lambda_x0 * (temp_dmiu5_dy + const_2_C2_y * temp_miu5);
    deriv[14] = temp * lambda_x0 * (temp_dmiu5_dz + const_2_C2_z * temp_miu5);

    // dmiu6 dx/dy/dz
    deriv[15] = temp * temp_miu6 * lambda * const_1_p_2_C2_x2;
    deriv[16] = temp * lambda_x0 * (temp_dmiu6_dy + const_2_C2_y * temp_miu6);
    deriv[17] = temp * lambda_x0 * (temp_dmiu6_dz + const_2_C2_z * temp_miu6);

    value[0] = miu_8_6_1;
    value[1] = miu_8_6_2;
    value[2] = miu_8_6_3;
    value[3] = miu_8_6_4;
    value[4] = miu_8_6_5;
    value[5] = miu_8_6_6;
}


void calc_MCSH_8_7(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double lambda_sqr = lambda * lambda;

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
    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x_2 = lambda_x0_sqr + (1.0/(2.0*gamma));
    double temp_y_2 = lambda_y0_sqr + (1.0/(2.0*gamma));
    double temp_z_2 = lambda_z0_sqr + (1.0/(2.0*gamma));

    double temp_dx_2 = 2.0 * lambda_sqr * x0;
    double temp_dy_2 = 2.0 * lambda_sqr * y0;
    double temp_dz_2 = 2.0 * lambda_sqr * z0;

    double C4_1 = 3.0 / gamma, C4_2 = 3.0 / (4.0 * gamma * gamma);
    double temp_x_4 = lambda_x0_4 + C4_1 * lambda_x0_sqr + C4_2;
    double temp_y_4 = lambda_y0_4 + C4_1 * lambda_y0_sqr + C4_2;
    double temp_z_4 = lambda_z0_4 + C4_1 * lambda_z0_sqr + C4_2;

    double temp_dx_4 = lambda * (4.0 * lambda_x0_3 + 2.0 * C4_1 * lambda_x0);
    double temp_dy_4 = lambda * (4.0 * lambda_y0_3 + 2.0 * C4_1 * lambda_y0);
    double temp_dz_4 = lambda * (4.0 * lambda_z0_3 + 2.0 * C4_1 * lambda_z0);


    double temp_term1_x = 2027025.0 * temp_x_4 - 810810.0 * temp_x_2 + 31185.0;
    double temp_term1_y = 2027025.0 * temp_y_4 - 810810.0 * temp_y_2 + 31185.0;
    // double temp_term1_z = 2027025.0 * temp_z_4 - 810810.0 * temp_z_2 + 31185.0;

    double temp_dterm1_dx = 2027025.0 * temp_dx_4 - 810810.0 * temp_dx_2;
    double temp_dterm1_dy = 2027025.0 * temp_dy_4 - 810810.0 * temp_dy_2;
    // double temp_dterm1_dz = 2027025.0 * temp_dz_4 - 810810.0 * temp_dz_2;

    double temp_term2_x = -810810.0 * temp_x_4 + 374220.0 * temp_x_2 - 17010.0;
    double temp_term2_y = -810810.0 * temp_y_4 + 374220.0 * temp_y_2 - 17010.0;
    // double temp_term2_z = -810810.0 * temp_z_4 + 374220.0 * temp_z_2 - 17010.0;

    double temp_dterm2_dx = -810810.0 * temp_dx_4 + 374220.0 * temp_dx_2;
    double temp_dterm2_dy = -810810.0 * temp_dy_4 + 374220.0 * temp_dy_2;
    // double temp_dterm2_dz = -810810.0 * temp_dz_4 + 374220.0 * temp_dz_2;

    double temp_term3_x = 31185.0 * temp_x_4 - 17010.0 * temp_x_2 + 945.0;
    double temp_term3_y = 31185.0 * temp_y_4 - 17010.0 * temp_y_2 + 945.0;
    // double temp_term3_z = 31185.0 * temp_z_4 - 17010.0 * temp_z_2 + 945.0;

    double temp_dterm3_dx = 31185.0 * temp_dx_4 - 17010.0 * temp_dx_2;
    double temp_dterm3_dy = 31185.0 * temp_dy_4 - 17010.0 * temp_dy_2;
    // double temp_dterm3_dz = 31185.0 * temp_dz_4 - 17010.0 * temp_dz_2;

    double temp_miu1 = temp_y_4 * temp_term1_x + temp_y_2 * temp_term2_x + temp_term3_x;
    double temp_miu2 = temp_z_4 * temp_term1_x + temp_z_2 * temp_term2_x + temp_term3_x;
    double temp_miu3 = temp_z_4 * temp_term1_y + temp_z_2 * temp_term2_y + temp_term3_y;

    double temp_dmiu1_dx = temp_y_4 * temp_dterm1_dx + temp_y_2 * temp_dterm2_dx + temp_dterm3_dx;
    double temp_dmiu1_dy = temp_dy_4 * temp_term1_x + temp_dy_2 * temp_term2_x;

    double temp_dmiu2_dx = temp_z_4 * temp_dterm1_dx + temp_z_2 * temp_dterm2_dx + temp_dterm3_dx;
    double temp_dmiu2_dz = temp_dz_4 * temp_term1_x + temp_dz_2 * temp_term2_x;

    double temp_dmiu3_dy = temp_z_4 * temp_dterm1_dy + temp_z_2 * temp_dterm2_dy + temp_dterm3_dy;
    double temp_dmiu3_dz = temp_dz_4 * temp_term1_y + temp_dz_2 * temp_term2_y;

    double miu_8_7_1 = temp * temp_miu1;
    double miu_8_7_2 = temp * temp_miu2;
    double miu_8_7_3 = temp * temp_miu3;

    //deriv consts
    double const_2_C2_x = 2.0 * C2 * x0; // 2 * C2 * x0
    double const_2_C2_y = 2.0 * C2 * y0;
    double const_2_C2_z = 2.0 * C2 * z0;

    // dmiu1 dx/dy/dz
    deriv[0] = temp * (temp_dmiu1_dx + const_2_C2_x * temp_miu1);
    deriv[1] = temp * (temp_dmiu1_dy + const_2_C2_y * temp_miu1);
    deriv[2] = miu_8_7_1 * const_2_C2_z;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * (temp_dmiu2_dx + const_2_C2_x * temp_miu2);
    deriv[4] = miu_8_7_2 * const_2_C2_y;
    deriv[5] = temp * (temp_dmiu2_dz + const_2_C2_z * temp_miu2);

    // dmiu3 dx/dy/dz
    deriv[6] = miu_8_7_3 * const_2_C2_x;
    deriv[7] = temp * (temp_dmiu3_dy + const_2_C2_y * temp_miu3);
    deriv[8] = temp * (temp_dmiu3_dz + const_2_C2_z * temp_miu3);

    value[0] = miu_8_7_1;
    value[1] = miu_8_7_2;
    value[2] = miu_8_7_3;
}


void calc_MCSH_8_8(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double lambda_sqr = lambda * lambda;
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

    double C4_1 = 3.0 / gamma, C4_2 = 3.0 / (4.0 * gamma * gamma);
    double temp_x_4 = lambda_x0_4 + C4_1 * lambda_x0_sqr + C4_2;
    double temp_y_4 = lambda_y0_4 + C4_1 * lambda_y0_sqr + C4_2;
    double temp_z_4 = lambda_z0_4 + C4_1 * lambda_z0_sqr + C4_2;

    double temp_dx_4 = lambda * (4.0 * lambda_x0_3 + 2.0 * C4_1 * lambda_x0);
    double temp_dy_4 = lambda * (4.0 * lambda_y0_3 + 2.0 * C4_1 * lambda_y0);
    double temp_dz_4 = lambda * (4.0 * lambda_z0_3 + 2.0 * C4_1 * lambda_z0);



    double temp_term1_x = 2027025.0 * temp_x_4 - 810810.0 * temp_x_2 + 31185.0;
    double temp_term1_y = 2027025.0 * temp_y_4 - 810810.0 * temp_y_2 + 31185.0;
    double temp_term1_z = 2027025.0 * temp_z_4 - 810810.0 * temp_z_2 + 31185.0;

    double temp_dterm1_dx = 2027025.0 * temp_dx_4 - 810810.0 * temp_dx_2;
    double temp_dterm1_dy = 2027025.0 * temp_dy_4 - 810810.0 * temp_dy_2;
    double temp_dterm1_dz = 2027025.0 * temp_dz_4 - 810810.0 * temp_dz_2;

    double temp_term2_x = -405405.0 * temp_x_4 + 187110.0 * temp_x_2 - 8505.0;
    double temp_term2_y = -405405.0 * temp_y_4 + 187110.0 * temp_y_2 - 8505.0;
    double temp_term2_z = -405405.0 * temp_z_4 + 187110.0 * temp_z_2 - 8505.0;

    double temp_dterm2_dx = -405405.0 * temp_dx_4 + 187110.0 * temp_dx_2;
    double temp_dterm2_dy = -405405.0 * temp_dy_4 + 187110.0 * temp_dy_2;
    double temp_dterm2_dz = -405405.0 * temp_dz_4 + 187110.0 * temp_dz_2;

    double temp_miu1 = temp_y_3 * temp_term1_x + lambda_y0 * temp_term2_x;
    double temp_miu2 = temp_x_3 * temp_term1_y + lambda_x0 * temp_term2_y;
    double temp_miu3 = temp_z_3 * temp_term1_x + lambda_z0 * temp_term2_x;
    double temp_miu4 = temp_x_3 * temp_term1_z + lambda_x0 * temp_term2_z;
    double temp_miu5 = temp_z_3 * temp_term1_y + lambda_z0 * temp_term2_y;
    double temp_miu6 = temp_y_3 * temp_term1_z + lambda_y0 * temp_term2_z;

    double temp_dmiu1_dx = temp_y_3 * temp_dterm1_dx + lambda_y0 * temp_dterm2_dx;
    double temp_dmiu1_dy = temp_dy_3 * temp_term1_x + lambda * temp_term2_x;

    double temp_dmiu2_dx = temp_dx_3 * temp_term1_y + lambda * temp_term2_y;
    double temp_dmiu2_dy = temp_x_3 * temp_dterm1_dy + lambda_x0 * temp_dterm2_dy;

    double temp_dmiu3_dx = temp_z_3 * temp_dterm1_dx + lambda_z0 * temp_dterm2_dx;
    double temp_dmiu3_dz = temp_dz_3 * temp_term1_x + lambda * temp_term2_x;

    double temp_dmiu4_dx = temp_dx_3 * temp_term1_z + lambda * temp_term2_z;
    double temp_dmiu4_dz = temp_x_3 * temp_dterm1_dz + lambda_x0 * temp_dterm2_dz;

    double temp_dmiu5_dy = temp_z_3 * temp_dterm1_dy + lambda_z0 * temp_dterm2_dy;
    double temp_dmiu5_dz = temp_dz_3 * temp_term1_y + lambda * temp_term2_y;

    double temp_dmiu6_dy = temp_dy_3 * temp_term1_z + lambda * temp_term2_z;
    double temp_dmiu6_dz = temp_y_3 * temp_dterm1_dz + lambda_y0 * temp_dterm2_dz;

    double miu_8_8_1 = temp * lambda_z0 * temp_miu1;
    double miu_8_8_2 = temp * lambda_z0 * temp_miu2;
    double miu_8_8_3 = temp * lambda_y0 * temp_miu3;
    double miu_8_8_4 = temp * lambda_y0 * temp_miu4;
    double miu_8_8_5 = temp * lambda_x0 * temp_miu5;
    double miu_8_8_6 = temp * lambda_x0 * temp_miu6;

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

    // dmiu2 dx/dy/dz
    deriv[3] = temp * lambda_z0 * (temp_dmiu2_dx + const_2_C2_x * temp_miu2);
    deriv[4] = temp * lambda_z0 * (temp_dmiu2_dy + const_2_C2_y * temp_miu2);
    deriv[5] = temp * temp_miu2 * lambda * const_1_p_2_C2_z2;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * lambda_y0 * (temp_dmiu3_dx + const_2_C2_x * temp_miu3);
    deriv[7] = temp * temp_miu3 * lambda * const_1_p_2_C2_y2;
    deriv[8] = temp * lambda_y0 * (temp_dmiu3_dz + const_2_C2_z * temp_miu3);

    // dmiu4 dx/dy/dz
    deriv[9] = temp * lambda_y0 * (temp_dmiu4_dx + const_2_C2_x * temp_miu4);
    deriv[10] = temp * temp_miu4 * lambda * const_1_p_2_C2_y2;
    deriv[11] = temp * lambda_y0 * (temp_dmiu4_dz + const_2_C2_z * temp_miu4);

    // dmiu5 dx/dy/dz
    deriv[12] = temp * temp_miu5 * lambda * const_1_p_2_C2_x2;
    deriv[13] = temp * lambda_x0 * (temp_dmiu5_dy + const_2_C2_y * temp_miu5);
    deriv[14] = temp * lambda_x0 * (temp_dmiu5_dz + const_2_C2_z * temp_miu5);

    // dmiu6 dx/dy/dz
    deriv[15] = temp * temp_miu6 * lambda * const_1_p_2_C2_x2;
    deriv[16] = temp * lambda_x0 * (temp_dmiu6_dy + const_2_C2_y * temp_miu6);
    deriv[17] = temp * lambda_x0 * (temp_dmiu6_dz + const_2_C2_z * temp_miu6);

    value[0] = miu_8_8_1;
    value[1] = miu_8_8_2;
    value[2] = miu_8_8_3;
    value[3] = miu_8_8_4;
    value[4] = miu_8_8_5;
    value[5] = miu_8_8_6;
}

void calc_MCSH_8_9(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double lambda_sqr = lambda * lambda;

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

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x_2 = lambda_x0_sqr + (1.0/(2.0*gamma));
    double temp_y_2 = lambda_y0_sqr + (1.0/(2.0*gamma));
    double temp_z_2 = lambda_z0_sqr + (1.0/(2.0*gamma));

    double temp_dx_2 = 2.0 * lambda_sqr * x0;
    double temp_dy_2 = 2.0 * lambda_sqr * y0;
    double temp_dz_2 = 2.0 * lambda_sqr * z0;

    double C4_1 = 3.0 / gamma, C4_2 = 3.0 / (4.0 * gamma * gamma);
    double temp_x_4 = lambda_x0_4 + C4_1 * lambda_x0_sqr + C4_2;
    double temp_y_4 = lambda_y0_4 + C4_1 * lambda_y0_sqr + C4_2;
    double temp_z_4 = lambda_z0_4 + C4_1 * lambda_z0_sqr + C4_2;

    double temp_dx_4 = lambda * (4.0 * lambda_x0_3 + 2.0 * C4_1 * lambda_x0);
    double temp_dy_4 = lambda * (4.0 * lambda_y0_3 + 2.0 * C4_1 * lambda_y0);
    double temp_dz_4 = lambda * (4.0 * lambda_z0_3 + 2.0 * C4_1 * lambda_z0);


    double temp_term1_x = 2027025.0 * temp_x_4 - 810810.0 * temp_x_2 + 31185.0;
    double temp_term1_y = 2027025.0 * temp_y_4 - 810810.0 * temp_y_2 + 31185.0;
    double temp_term1_z = 2027025.0 * temp_z_4 - 810810.0 * temp_z_2 + 31185.0;

    double temp_dterm1_dx = 2027025.0 * temp_dx_4 - 810810.0 * temp_dx_2;
    double temp_dterm1_dy = 2027025.0 * temp_dy_4 - 810810.0 * temp_dy_2;
    double temp_dterm1_dz = 2027025.0 * temp_dz_4 - 810810.0 * temp_dz_2;

    double temp_term2_x = -135135.0 * temp_x_4 + 62370.0 * temp_x_2 - 2835.0;
    double temp_term2_y = -135135.0 * temp_y_4 + 62370.0 * temp_y_2 - 2835.0;
    double temp_term2_z = -135135.0 * temp_z_4 + 62370.0 * temp_z_2 - 2835.0;

    double temp_dterm2_dx = -135135.0 * temp_dx_4 + 62370.0 * temp_dx_2;
    double temp_dterm2_dy = -135135.0 * temp_dy_4 + 62370.0 * temp_dy_2;
    double temp_dterm2_dz = -135135.0 * temp_dz_4 + 62370.0 * temp_dz_2;

    double temp_term3_x = 10395.0 * temp_x_4 - 5670.0 * temp_x_2 + 315.0;
    double temp_term3_y = 10395.0 * temp_y_4 - 5670.0 * temp_y_2 + 315.0;
    double temp_term3_z = 10395.0 * temp_z_4 - 5670.0 * temp_z_2 + 315.0;

    double temp_dterm3_dx = 10395.0 * temp_dx_4 - 5670.0 * temp_dx_2;
    double temp_dterm3_dy = 10395.0 * temp_dy_4 - 5670.0 * temp_dy_2;
    double temp_dterm3_dz = 10395.0 * temp_dz_4 - 5670.0 * temp_dz_2;

    double temp_miu1 = temp_y_2 * temp_z_2 * temp_term1_x + (temp_y_2 + temp_z_2) * temp_term2_x + temp_term3_x;
    double temp_miu2 = temp_x_2 * temp_z_2 * temp_term1_y + (temp_x_2 + temp_z_2) * temp_term2_y + temp_term3_y;
    double temp_miu3 = temp_x_2 * temp_y_2 * temp_term1_z + (temp_x_2 + temp_y_2) * temp_term2_z + temp_term3_z;

    double temp_dmiu1_dx = temp_y_2 * temp_z_2 * temp_dterm1_dx + (temp_y_2 + temp_z_2) * temp_dterm2_dx + temp_dterm3_dx;
    double temp_dmiu1_dy = temp_dy_2 * temp_z_2 * temp_term1_x + temp_dy_2 * temp_term2_x;
    double temp_dmiu1_dz = temp_y_2 * temp_dz_2 * temp_term1_x + temp_dz_2 * temp_term2_x;

    double temp_dmiu2_dx = temp_dx_2 * temp_z_2 * temp_term1_y + temp_dx_2 * temp_term2_y;
    double temp_dmiu2_dy = temp_x_2 * temp_z_2 * temp_dterm1_dy + (temp_x_2 + temp_z_2) * temp_dterm2_dy + temp_dterm3_dy;
    double temp_dmiu2_dz = temp_x_2 * temp_dz_2 * temp_term1_y + temp_dz_2 * temp_term2_y;

    double temp_dmiu3_dx = temp_dx_2 * temp_y_2 * temp_term1_z + temp_dx_2 * temp_term2_z;
    double temp_dmiu3_dy = temp_x_2 * temp_dy_2 * temp_term1_z + temp_dy_2 * temp_term2_z;
    double temp_dmiu3_dz = temp_x_2 * temp_y_2 * temp_dterm1_dz + (temp_x_2 + temp_y_2) * temp_dterm2_dz + temp_dterm3_dz;


    double miu_8_9_1 = temp * temp_miu1;
    double miu_8_9_2 = temp * temp_miu2;
    double miu_8_9_3 = temp * temp_miu3;

    //deriv consts
    double const_2_C2_x = 2.0 * C2 * x0; // 2 * C2 * x0
    double const_2_C2_y = 2.0 * C2 * y0;
    double const_2_C2_z = 2.0 * C2 * z0;

    // dmiu1 dx/dy/dz
    deriv[0] = temp * (temp_dmiu1_dx + const_2_C2_x * temp_miu1);
    deriv[1] = temp * (temp_dmiu1_dy + const_2_C2_y * temp_miu1);
    deriv[2] = temp * (temp_dmiu1_dz + const_2_C2_z * temp_miu1);

    // dmiu2 dx/dy/dz
    deriv[3] = temp * (temp_dmiu2_dx + const_2_C2_x * temp_miu2);
    deriv[4] = temp * (temp_dmiu2_dy + const_2_C2_y * temp_miu2);
    deriv[5] = temp * (temp_dmiu2_dz + const_2_C2_z * temp_miu2);

    // dmiu3 dx/dy/dz
    deriv[6] = temp * (temp_dmiu3_dx + const_2_C2_x * temp_miu3);
    deriv[7] = temp * (temp_dmiu3_dy + const_2_C2_y * temp_miu3);
    deriv[8] = temp * (temp_dmiu3_dz + const_2_C2_z * temp_miu3);

    value[0] = miu_8_9_1;
    value[1] = miu_8_9_2;
    value[2] = miu_8_9_3;
}

void calc_MCSH_8_10(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double lambda_sqr = lambda * lambda;

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


    double temp_term1_x = 2027025.0 * temp_x_3 - 405405.0 * lambda_x0;
    double temp_term1_y = 2027025.0 * temp_y_3 - 405405.0 * lambda_y0;
    // double temp_term1_z = 2027025.0 * temp_z_3 - 405405.0 * lambda_z0;

    double temp_dterm1_dx = 2027025.0 * temp_dx_3 - 405405.0 * lambda;
    double temp_dterm1_dy = 2027025.0 * temp_dy_3 - 405405.0 * lambda;
    // double temp_dterm1_dz = 2027025.0 * temp_dz_3 - 405405.0 * lambda;

    double temp_term2_x = -405405.0 * temp_x_3 + 93555.0 * lambda_x0;
    double temp_term2_y = -405405.0 * temp_y_3 + 93555.0 * lambda_y0;
    // double temp_term2_z = -405405.0 * temp_z_3 + 93555.0 * lambda_z0;

    double temp_dterm2_dx = -405405.0 * temp_dx_3 + 93555.0 * lambda;
    double temp_dterm2_dy = -405405.0 * temp_dy_3 + 93555.0 * lambda;
    // double temp_dterm2_dz = -405405.0 * temp_dz_3 + 93555.0 * lambda;

    double temp_term3_x = -135135.0 * temp_x_3 + 31185.0 * lambda_x0;
    double temp_term3_y = -135135.0 * temp_y_3 + 31185.0 * lambda_y0;
    // double temp_term3_z = -135135.0 * temp_z_3 + 31185.0 * lambda_z0;

    double temp_dterm3_dx = -135135.0 * temp_dx_3 + 31185.0 * lambda;
    double temp_dterm3_dy = -135135.0 * temp_dy_3 + 31185.0 * lambda;
    // double temp_dterm3_dz = -135135.0 * temp_dz_3 + 31185.0 * lambda;

    double temp_term4_x = 31185.0 * temp_x_3 - 8505.0 * lambda_x0;
    double temp_term4_y = 31185.0 * temp_y_3 - 8505.0 * lambda_y0;
    // double temp_term4_z = 31185.0 * temp_z_3 - 8505.0 * lambda_z0;

    double temp_dterm4_dx = 31185.0 * temp_dx_3 - 8505.0 * lambda;
    double temp_dterm4_dy = 31185.0 * temp_dy_3 - 8505.0 * lambda;
    // double temp_dterm4_dz = 31185.0 * temp_dz_3 - 8505.0 * lambda;



    double temp_miu1 = temp_y_3 * temp_z_2 * temp_term1_x + lambda_y0 * temp_z_2 * temp_term2_x + temp_y_3 * temp_term3_x + lambda_y0 * temp_term4_x;
    double temp_miu2 = temp_z_3 * temp_y_2 * temp_term1_x + lambda_z0 * temp_y_2 * temp_term2_x + temp_z_3 * temp_term3_x + lambda_z0 * temp_term4_x;
    double temp_miu3 = temp_z_3 * temp_x_2 * temp_term1_y + lambda_z0 * temp_x_2 * temp_term2_y + temp_z_3 * temp_term3_y + lambda_z0 * temp_term4_y;

    double temp_dmiu1_dx = temp_y_3 * temp_z_2 * temp_dterm1_dx + lambda_y0 * temp_z_2 * temp_dterm2_dx + temp_y_3 * temp_dterm3_dx + lambda_y0 * temp_dterm4_dx;
    double temp_dmiu1_dy = temp_dy_3 * temp_z_2 * temp_term1_x + lambda * temp_z_2 * temp_term2_x + temp_dy_3 * temp_term3_x + lambda * temp_term4_x;
    double temp_dmiu1_dz = temp_y_3 * temp_dz_2 * temp_term1_x + lambda_y0 * temp_dz_2 * temp_term2_x;

    double temp_dmiu2_dx = temp_z_3 * temp_y_2 * temp_dterm1_dx + lambda_z0 * temp_y_2 * temp_dterm2_dx + temp_z_3 * temp_dterm3_dx + lambda_z0 * temp_dterm4_dx;
    double temp_dmiu2_dy = temp_z_3 * temp_dy_2 * temp_term1_x + lambda_z0 * temp_dy_2 * temp_term2_x;
    double temp_dmiu2_dz = temp_dz_3 * temp_y_2 * temp_term1_x + lambda * temp_y_2 * temp_term2_x + temp_dz_3 * temp_term3_x + lambda * temp_term4_x;

    double temp_dmiu3_dx = temp_z_3 * temp_dx_2 * temp_term1_y + lambda_z0 * temp_dx_2 * temp_term2_y;
    double temp_dmiu3_dy = temp_z_3 * temp_x_2 * temp_dterm1_dy + lambda_z0 * temp_x_2 * temp_dterm2_dy + temp_z_3 * temp_dterm3_dy + lambda_z0 * temp_dterm4_dy;
    double temp_dmiu3_dz = temp_dz_3 * temp_x_2 * temp_term1_y + lambda * temp_x_2 * temp_term2_y + temp_dz_3 * temp_term3_y + lambda * temp_term4_y;

    double miu_8_10_1 = temp * temp_miu1;
    double miu_8_10_2 = temp * temp_miu2;
    double miu_8_10_3 = temp * temp_miu3;

    //deriv consts
    double const_2_C2_x = 2.0 * C2 * x0; // 2 * C2 * x0
    double const_2_C2_y = 2.0 * C2 * y0;
    double const_2_C2_z = 2.0 * C2 * z0;

    // dmiu1 dx/dy/dz
    deriv[0] = temp * (temp_dmiu1_dx + const_2_C2_x * temp_miu1);
    deriv[1] = temp * (temp_dmiu1_dy + const_2_C2_y * temp_miu1);
    deriv[2] = temp * (temp_dmiu1_dz + const_2_C2_z * temp_miu1);

    // dmiu2 dx/dy/dz
    deriv[3] = temp * (temp_dmiu2_dx + const_2_C2_x * temp_miu2);
    deriv[4] = temp * (temp_dmiu2_dy + const_2_C2_y * temp_miu2);
    deriv[5] = temp * (temp_dmiu2_dz + const_2_C2_z * temp_miu2);

    // dmiu3 dx/dy/dz
    deriv[6] = temp * (temp_dmiu3_dx + const_2_C2_x * temp_miu3);
    deriv[7] = temp * (temp_dmiu3_dy + const_2_C2_y * temp_miu3);
    deriv[8] = temp * (temp_dmiu3_dz + const_2_C2_z * temp_miu3);

    value[0] = miu_8_10_1;
    value[1] = miu_8_10_2;
    value[2] = miu_8_10_3;
}

void calc_MCSH_9_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    // double lambda_sqr = lambda * lambda;
    // double x0_sqr = x0*x0;
    // double y0_sqr = y0*y0;
    // double z0_sqr = z0*z0;

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

    double lambda_x0_6 = lambda_x0_5 * lambda_x0;
    double lambda_y0_6 = lambda_y0_5 * lambda_y0;
    double lambda_z0_6 = lambda_z0_5 * lambda_z0;

    double lambda_x0_7 = lambda_x0_6 * lambda_x0;
    double lambda_y0_7 = lambda_y0_6 * lambda_y0;
    double lambda_z0_7 = lambda_z0_6 * lambda_z0;

    double lambda_x0_8 = lambda_x0_7 * lambda_x0;
    double lambda_y0_8 = lambda_y0_7 * lambda_y0;
    double lambda_z0_8 = lambda_z0_7 * lambda_z0;

    double lambda_x0_9 = lambda_x0_8 * lambda_x0;
    double lambda_y0_9 = lambda_y0_8 * lambda_y0;
    double lambda_z0_9 = lambda_z0_8 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = (34459425.0 * 18.0 / gamma) - 72972900.0;
    double C4 = (34459425.0 * 189.0 / (2.0 * gamma * gamma)) - (72972900.0 * 21.0 / (2.0 * gamma)) + 51081030.0;
    double C5 = (34459425.0 * 315.0 / (2.0 * gamma * gamma * gamma)) - (72972900.0 * 105.0 / (4.0 * gamma * gamma)) + (51081030.0 * 5.0 / gamma) - 13097700.0;
    double C6 = (34459425.0 * 945.0 / (16.0 * gamma * gamma * gamma * gamma)) - (72972900.0 * 105.0 / (8.0 * gamma * gamma * gamma)) + (51081030.0 * 15.0 / (4.0 * gamma *gamma)) - (13097700.0 * 3.0 / (2.0 * gamma)) + 893025.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 34459425.0 * lambda_x0_9 + C3 * lambda_x0_7 + C4 * lambda_x0_5 + C5 * lambda_x0_3 + C6 * lambda_x0;
    double temp_y = 34459425.0 * lambda_y0_9 + C3 * lambda_y0_7 + C4 * lambda_y0_5 + C5 * lambda_y0_3 + C6 * lambda_y0;
    double temp_z = 34459425.0 * lambda_z0_9 + C3 * lambda_z0_7 + C4 * lambda_z0_5 + C5 * lambda_z0_3 + C6 * lambda_z0;

    double temp_dx = lambda * (34459425.0 * 9.0 * lambda_x0_8 + 7.0 * C3 * lambda_x0_6 + 5.0 * C4 * lambda_x0_4 + 3.0 * C5 * lambda_x0_sqr + C6);
    double temp_dy = lambda * (34459425.0 * 9.0 * lambda_y0_8 + 7.0 * C3 * lambda_y0_6 + 5.0 * C4 * lambda_y0_4 + 3.0 * C5 * lambda_y0_sqr + C6);
    double temp_dz = lambda * (34459425.0 * 9.0 * lambda_z0_8 + 7.0 * C3 * lambda_z0_6 + 5.0 * C4 * lambda_z0_4 + 3.0 * C5 * lambda_z0_sqr + C6);

    double miu_9_1_1 = temp * temp_x;
    double miu_9_1_2 = temp * temp_y;
    double miu_9_1_3 = temp * temp_z;

    //deriv consts
    double const_2_C2_x = 2.0 * C2 * x0;
    double const_2_C2_y = 2.0 * C2 * y0;
    double const_2_C2_z = 2.0 * C2 * z0;

    // dmiu1 dx/dy/dz
    deriv[0] = temp * (temp_dx + const_2_C2_x * temp_x);
    deriv[1] = miu_9_1_1 * const_2_C2_y;
    deriv[2] = miu_9_1_1 * const_2_C2_z;

    // dmiu2 dx/dy/dz
    deriv[3] = miu_9_1_2 * const_2_C2_x;
    deriv[4] = temp * (temp_dy + const_2_C2_y * temp_y);
    deriv[5] = miu_9_1_2 * const_2_C2_z;

    // dmiu3 dx/dy/dz
    deriv[6] = miu_9_1_3 * const_2_C2_x;
    deriv[7] = miu_9_1_3 * const_2_C2_y;
    deriv[8] = temp * (temp_dz + const_2_C2_z * temp_z);

    value[0] = miu_9_1_1;
    value[1] = miu_9_1_2;
    value[2] = miu_9_1_3;
}

void calc_MCSH_9_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
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

    double lambda_x0_5 = lambda_x0_4 * lambda_x0;
    double lambda_y0_5 = lambda_y0_4 * lambda_y0;
    double lambda_z0_5 = lambda_z0_4 * lambda_z0;

    double lambda_x0_6 = lambda_x0_5 * lambda_x0;
    double lambda_y0_6 = lambda_y0_5 * lambda_y0;
    double lambda_z0_6 = lambda_z0_5 * lambda_z0;

    double lambda_x0_7 = lambda_x0_6 * lambda_x0;
    double lambda_y0_7 = lambda_y0_6 * lambda_y0;
    double lambda_z0_7 = lambda_z0_6 * lambda_z0;

    double lambda_x0_8 = lambda_x0_7 * lambda_x0;
    double lambda_y0_8 = lambda_y0_7 * lambda_y0;
    double lambda_z0_8 = lambda_z0_7 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = (34459425.0 * 14.0 / gamma) - 56756700.0;
    double C4 = (34459425.0 * 105.0 / (2.0 * gamma * gamma)) - (56756700.0 * 15.0 / (2.0 * gamma)) + 28378350.0;
    double C5 = (34459425.0 * 105.0 / (2.0 * gamma * gamma * gamma)) - (56756700.0 * 45.20 / (4.0 * gamma * gamma)) + (28378350.0 * 3.0 / gamma) - 4365900.0;
    double C6 = (34459425.0 * 105.0 / (16.0 * gamma * gamma * gamma * gamma)) - (56756700.0 * 15.0 / (8.0 * gamma * gamma * gamma)) + (28378350.0 * 3.0 / (4.0 * gamma * gamma)) - (4365900.0 / (2.0 * gamma)) + 99225.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 34459425.0 * lambda_x0_8 + C3 * lambda_x0_6 + C4 * lambda_x0_4 + C5 * lambda_x0_sqr + C6;
    double temp_y = 34459425.0 * lambda_y0_8 + C3 * lambda_y0_6 + C4 * lambda_y0_4 + C5 * lambda_y0_sqr + C6;
    double temp_z = 34459425.0 * lambda_z0_8 + C3 * lambda_z0_6 + C4 * lambda_z0_4 + C5 * lambda_z0_sqr + C6;

    double temp_dx = lambda * (34459425.0 * 8.0 * lambda_x0_7 + 6.0 * C3 * lambda_x0_5 + 4.0 * C4 * lambda_x0_3 + 2.0 * C5 * lambda_x0);
    double temp_dy = lambda * (34459425.0 * 8.0 * lambda_y0_7 + 6.0 * C3 * lambda_y0_5 + 4.0 * C4 * lambda_y0_3 + 2.0 * C5 * lambda_y0);
    double temp_dz = lambda * (34459425.0 * 8.0 * lambda_z0_7 + 6.0 * C3 * lambda_z0_5 + 4.0 * C4 * lambda_z0_3 + 2.0 * C5 * lambda_z0);

    double miu_9_2_1 = temp * lambda_y0 * temp_x;
    double miu_9_2_2 = temp * lambda_x0 * temp_y;
    double miu_9_2_3 = temp * lambda_z0 * temp_x;
    double miu_9_2_4 = temp * lambda_x0 * temp_z;
    double miu_9_2_5 = temp * lambda_z0 * temp_y;
    double miu_9_2_6 = temp * lambda_y0 * temp_z;


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
    deriv[2] = miu_9_2_1 * const_2_C2_z;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * lambda * temp_y * const_1_p_2_C2_x2;
    deriv[4] = temp * lambda_x0 * (temp_dy + const_2_C2_y * temp_y);
    deriv[5] = miu_9_2_2 * const_2_C2_z;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * lambda_z0 * (temp_dx + const_2_C2_x * temp_x);
    deriv[7] = miu_9_2_3 * const_2_C2_y;
    deriv[8] = temp * lambda * temp_x * const_1_p_2_C2_z2;

    // dmiu4 dx/dy/dz
    deriv[9] = temp * lambda * temp_z * const_1_p_2_C2_x2;
    deriv[10] = miu_9_2_4 * const_2_C2_y;
    deriv[11] = temp * lambda_x0 * (temp_dz + const_2_C2_z * temp_z);

    // dmiu5 dx/dy/dz
    deriv[12] = miu_9_2_5 * const_2_C2_x;
    deriv[13] = temp * lambda_z0 * (temp_dy + const_2_C2_y * temp_y);
    deriv[14] = temp * lambda * temp_y * const_1_p_2_C2_z2;

    // dmiu6 dx/dy/dz
    deriv[15] = miu_9_2_6 * const_2_C2_x;
    deriv[16] = temp * lambda * temp_z * const_1_p_2_C2_y2;
    deriv[17] = temp * lambda_y0 * (temp_dz + const_2_C2_z * temp_z);

    value[0] = miu_9_2_1;
    value[1] = miu_9_2_2;
    value[2] = miu_9_2_3;
    value[3] = miu_9_2_4;
    value[4] = miu_9_2_5;
    value[5] = miu_9_2_6;
}

void calc_MCSH_9_3(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double lambda_sqr = lambda * lambda;

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

    double lambda_x0_5 = lambda_x0_4 * lambda_x0;
    double lambda_y0_5 = lambda_y0_4 * lambda_y0;
    double lambda_z0_5 = lambda_z0_4 * lambda_z0;

    double lambda_x0_6 = lambda_x0_5 * lambda_x0;
    double lambda_y0_6 = lambda_y0_5 * lambda_y0;
    double lambda_z0_6 = lambda_z0_5 * lambda_z0;

    double lambda_x0_7 = lambda_x0_6 * lambda_x0;
    double lambda_y0_7 = lambda_y0_6 * lambda_y0;
    double lambda_z0_7 = lambda_z0_6 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x_2 = lambda_x0_sqr + (1/(2*gamma));
    double temp_y_2 = lambda_y0_sqr + (1/(2*gamma));
    double temp_z_2 = lambda_z0_sqr + (1/(2*gamma));

    double temp_dx_2 = 2 * lambda_sqr * x0;
    double temp_dy_2 = 2 * lambda_sqr * y0;
    double temp_dz_2 = 2 * lambda_sqr * z0;

    double C3_1 = 3.0 / (2.0 * gamma);
    double temp_x_3 = lambda_x0_3 + C3_1 * lambda_x0;
    double temp_y_3 = lambda_y0_3 + C3_1 * lambda_y0;
    double temp_z_3 = lambda_z0_3 + C3_1 * lambda_z0;

    double temp_dx_3 = lambda * (3.0 * lambda_x0_sqr + C3_1);
    double temp_dy_3 = lambda * (3.0 * lambda_y0_sqr + C3_1);
    double temp_dz_3 = lambda * (3.0 * lambda_z0_sqr + C3_1);

    double C5_1 = 5.0 / gamma, C5_2 = 15.0 / (4.0 * gamma * gamma);
    double temp_x_5 = lambda_x0_5 + C5_1 * lambda_x0_3 + C5_2 * lambda_x0;
    double temp_y_5 = lambda_y0_5 + C5_1 * lambda_y0_3 + C5_2 * lambda_y0;
    double temp_z_5 = lambda_z0_5 + C5_1 * lambda_z0_3 + C5_2 * lambda_z0;

    double temp_dx_5 = lambda * (5.0 * lambda_x0_4 + 3.0 * C5_1 * lambda_x0_sqr + C5_2);
    double temp_dy_5 = lambda * (5.0 * lambda_y0_4 + 3.0 * C5_1 * lambda_y0_sqr + C5_2);
    double temp_dz_5 = lambda * (5.0 * lambda_z0_4 + 3.0 * C5_1 * lambda_z0_sqr + C5_2);

    double C7_1 = 21 / (2 * gamma), C7_2 = 105 / (4 * gamma * gamma), C7_3 = 105 / (8 * gamma * gamma * gamma);
    double temp_x_7 = lambda_x0_7 + C7_1 * lambda_x0_5 + C7_2 * lambda_x0_3 + C7_3 * lambda_x0;
    double temp_y_7 = lambda_y0_7 + C7_1 * lambda_y0_5 + C7_2 * lambda_y0_3 + C7_3 * lambda_y0;
    double temp_z_7 = lambda_z0_7 + C7_1 * lambda_z0_5 + C7_2 * lambda_z0_3 + C7_3 * lambda_z0;

    double temp_dx_7 = lambda * (7 * lambda_x0_6 + 5 * C7_1 * lambda_x0_4 + 3 * C7_2 * lambda_x0_sqr + C7_3);
    double temp_dy_7 = lambda * (7 * lambda_y0_6 + 5 * C7_1 * lambda_y0_4 + 3 * C7_2 * lambda_y0_sqr + C7_3);
    double temp_dz_7 = lambda * (7 * lambda_z0_6 + 5 * C7_1 * lambda_z0_4 + 3 * C7_2 * lambda_z0_sqr + C7_3);


    double temp_term1_x = 34459425.0 * temp_x_7 - 42567525.0 * temp_x_5 + 14189175.0 * temp_x_3 - 1091475.0 * lambda_x0;
    double temp_term1_y = 34459425.0 * temp_y_7 - 42567525.0 * temp_y_5 + 14189175.0 * temp_y_3 - 1091475.0 * lambda_y0;
    double temp_term1_z = 34459425.0 * temp_z_7 - 42567525.0 * temp_z_5 + 14189175.0 * temp_z_3 - 1091475.0 * lambda_z0;

    double temp_dterm1_dx = 34459425.0 * temp_dx_7 - 42567525.0 * temp_dx_5 + 14189175.0 * temp_dx_3 - 1091475.0 * lambda;
    double temp_dterm1_dy = 34459425.0 * temp_dy_7 - 42567525.0 * temp_dy_5 + 14189175.0 * temp_dy_3 - 1091475.0 * lambda;
    double temp_dterm1_dz = 34459425.0 * temp_dz_7 - 42567525.0 * temp_dz_5 + 14189175.0 * temp_dz_3 - 1091475.0 * lambda;

    double temp_term2_x = -2027025.0 * temp_x_7 + 2837835.0 * temp_x_5 - 1091475.0 * temp_x_3 + 99225.0 * lambda_x0;
    double temp_term2_y = -2027025.0 * temp_y_7 + 2837835.0 * temp_y_5 - 1091475.0 * temp_y_3 + 99225.0 * lambda_y0;
    double temp_term2_z = -2027025.0 * temp_z_7 + 2837835.0 * temp_z_5 - 1091475.0 * temp_z_3 + 99225.0 * lambda_z0;

    double temp_dterm2_dx = -2027025.0 * temp_dx_7 + 2837835.0 * temp_dx_5 - 1091475.0 * temp_dx_3 + 99225.0 * lambda;
    double temp_dterm2_dy = -2027025.0 * temp_dy_7 + 2837835.0 * temp_dy_5 - 1091475.0 * temp_dy_3 + 99225.0 * lambda;
    double temp_dterm2_dz = -2027025.0 * temp_dz_7 + 2837835.0 * temp_dz_5 - 1091475.0 * temp_dz_3 + 99225.0 * lambda;

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

    double miu_9_3_1 = temp * temp_miu1;
    double miu_9_3_2 = temp * temp_miu2;
    double miu_9_3_3 = temp * temp_miu3;
    double miu_9_3_4 = temp * temp_miu4;
    double miu_9_3_5 = temp * temp_miu5;
    double miu_9_3_6 = temp * temp_miu6;

    //deriv consts
    double const_2_C2_x = 2.0 * C2 * x0; // 2 * C2 * x0
    double const_2_C2_y = 2.0 * C2 * y0;
    double const_2_C2_z = 2.0 * C2 * z0;

    // dmiu1 dx/dy/dz
    deriv[0] = temp * (temp_dmiu1_dx + const_2_C2_x * temp_miu1);
    deriv[1] = temp * (temp_dmiu1_dy + const_2_C2_y * temp_miu1);
    deriv[2] = miu_9_3_1 * const_2_C2_z;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * (temp_dmiu2_dx + const_2_C2_x * temp_miu2);
    deriv[4] = temp * (temp_dmiu2_dy + const_2_C2_y * temp_miu2);
    deriv[5] = miu_9_3_2 * const_2_C2_z;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * (temp_dmiu3_dx + const_2_C2_x * temp_miu3);
    deriv[7] = miu_9_3_3 * const_2_C2_y;
    deriv[8] = temp * (temp_dmiu3_dz + const_2_C2_z * temp_miu3);

    // dmiu4 dx/dy/dz
    deriv[9] = temp * (temp_dmiu4_dx + const_2_C2_x * temp_miu4);
    deriv[10] = miu_9_3_4 * const_2_C2_y;
    deriv[11] = temp * (temp_dmiu4_dz + const_2_C2_z * temp_miu4);

    // dmiu5 dx/dy/dz
    deriv[12] = miu_9_3_5 * const_2_C2_x;
    deriv[13] = temp * (temp_dmiu5_dy + const_2_C2_y * temp_miu5);
    deriv[14] = temp * (temp_dmiu5_dz + const_2_C2_z * temp_miu5);

    // dmiu6 dx/dy/dz
    deriv[15] = miu_9_3_6 * const_2_C2_x;
    deriv[16] = temp * (temp_dmiu6_dy + const_2_C2_y * temp_miu6);
    deriv[17] = temp * (temp_dmiu6_dz + const_2_C2_z * temp_miu6);

    value[0] = miu_9_3_1;
    value[1] = miu_9_3_2;
    value[2] = miu_9_3_3;
    value[3] = miu_9_3_4;
    value[4] = miu_9_3_5;
    value[5] = miu_9_3_6;
}

void calc_MCSH_9_4(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
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

    double lambda_x0_5 = lambda_x0_4 * lambda_x0;
    double lambda_y0_5 = lambda_y0_4 * lambda_y0;
    double lambda_z0_5 = lambda_z0_4 * lambda_z0;

    double lambda_x0_6 = lambda_x0_5 * lambda_x0;
    double lambda_y0_6 = lambda_y0_5 * lambda_y0;
    double lambda_z0_6 = lambda_z0_5 * lambda_z0;

    double lambda_x0_7 = lambda_x0_6 * lambda_x0;
    double lambda_y0_7 = lambda_y0_6 * lambda_y0;
    double lambda_z0_7 = lambda_z0_6 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = (34459425.0 * 21.0 / (2.0 * gamma)) - 42567525.0;
    double C4 = (34459425.0 * 105.0 / (4.0 * gamma * gamma)) - (42567525.0 * 5.0 / gamma) + 14189175.0;
    double C5 = (34459425.0 * 105.0 / (8.0 * gamma * gamma * gamma)) - (42567525.0 * 15.0 / (4.0 * gamma * gamma)) + (14189175.0 * 3.0 / (2.0 * gamma)) - 1091475.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 34459425.0 * lambda_x0_7 + C3 * lambda_x0_5 + C4 * lambda_x0_3 + C5 * lambda_x0;
    double temp_y = 34459425.0 * lambda_y0_7 + C3 * lambda_y0_5 + C4 * lambda_y0_3 + C5 * lambda_y0;
    double temp_z = 34459425.0 * lambda_z0_7 + C3 * lambda_z0_5 + C4 * lambda_z0_3 + C5 * lambda_z0;

    double temp_dx = lambda * (34459425.0 * 7.0 * lambda_x0_6 + 5.0 * C3 * lambda_x0_4 + 3.0 * C4 * lambda_x0_sqr + C5);
    double temp_dy = lambda * (34459425.0 * 7.0 * lambda_y0_6 + 5.0 * C3 * lambda_y0_4 + 3.0 * C4 * lambda_y0_sqr + C5);
    double temp_dz = lambda * (34459425.0 * 7.0 * lambda_z0_6 + 5.0 * C3 * lambda_z0_4 + 3.0 * C4 * lambda_z0_sqr + C5);

    double miu_9_4_1 = temp * lambda_y0 * lambda_z0 * temp_x;
    double miu_9_4_2 = temp * lambda_x0 * lambda_z0 * temp_y;
    double miu_9_4_3 = temp * lambda_x0 * lambda_y0 * temp_z;

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

    value[0] = miu_9_4_1;
    value[1] = miu_9_4_2;
    value[2] = miu_9_4_3;
}


void calc_MCSH_9_5(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double lambda_sqr = lambda * lambda;

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

    double lambda_x0_5 = lambda_x0_4 * lambda_x0;
    double lambda_y0_5 = lambda_y0_4 * lambda_y0;
    double lambda_z0_5 = lambda_z0_4 * lambda_z0;

    double lambda_x0_6 = lambda_x0_5 * lambda_x0;
    double lambda_y0_6 = lambda_y0_5 * lambda_y0;
    double lambda_z0_6 = lambda_z0_5 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);

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

    double C4_1 = 3.0 / gamma, C4_2 = 3.0 / (4.0 * gamma * gamma);
    double temp_x_4 = lambda_x0_4 + C4_1 * lambda_x0_sqr + C4_2;
    double temp_y_4 = lambda_y0_4 + C4_1 * lambda_y0_sqr + C4_2;
    double temp_z_4 = lambda_z0_4 + C4_1 * lambda_z0_sqr + C4_2;

    double temp_dx_4 = lambda * (4.0 * lambda_x0_3 + 2.0 * C4_1 * lambda_x0);
    double temp_dy_4 = lambda * (4.0 * lambda_y0_3 + 2.0 * C4_1 * lambda_y0);
    double temp_dz_4 = lambda * (4.0 * lambda_z0_3 + 2.0 * C4_1 * lambda_z0);

    double C6_1 = 15.0 / (2.0 * gamma), C6_2 = 45.0 / (4.0 * gamma * gamma), C6_3 =  15.0 / (8.0 * gamma * gamma * gamma);
    double temp_x_6 = lambda_x0_6 + C6_1 * lambda_x0_4 + C6_2 * lambda_x0_sqr + C6_3;
    double temp_y_6 = lambda_y0_6 + C6_1 * lambda_y0_4 + C6_2 * lambda_y0_sqr + C6_3;
    double temp_z_6 = lambda_z0_6 + C6_1 * lambda_z0_4 + C6_2 * lambda_z0_sqr + C6_3;

    double temp_dx_6 = lambda * (6.0 * lambda_x0_5 + 4.0 * C6_1 * lambda_x0_3 + 2.0 * C6_2 * lambda_x0);
    double temp_dy_6 = lambda * (6.0 * lambda_y0_5 + 4.0 * C6_1 * lambda_y0_3 + 2.0 * C6_2 * lambda_y0);
    double temp_dz_6 = lambda * (6.0 * lambda_z0_5 + 4.0 * C6_1 * lambda_z0_3 + 2.0 * C6_2 * lambda_z0);


    double temp_term1_x = 34459425.0 * temp_x_6 - 30405375.0 * temp_x_4 + 6081075.0 * temp_x_2 - 155925.0;
    double temp_term1_y = 34459425.0 * temp_y_6 - 30405375.0 * temp_y_4 + 6081075.0 * temp_y_2 - 155925.0;
    double temp_term1_z = 34459425.0 * temp_z_6 - 30405375.0 * temp_z_4 + 6081075.0 * temp_z_2 - 155925.0;

    double temp_dterm1_dx = 34459425.0 * temp_dx_6 - 30405375.0 * temp_dx_4 + 6081075.0 * temp_dx_2;
    double temp_dterm1_dy = 34459425.0 * temp_dy_6 - 30405375.0 * temp_dy_4 + 6081075.0 * temp_dy_2;
    double temp_dterm1_dz = 34459425.0 * temp_dz_6 - 30405375.0 * temp_dz_4 + 6081075.0 * temp_dz_2;

    double temp_term2_x = -6081075.0 * temp_x_6 + 6081075.0 * temp_x_4 - 1403325.0 * temp_x_2 + 42525.0;
    double temp_term2_y = -6081075.0 * temp_y_6 + 6081075.0 * temp_y_4 - 1403325.0 * temp_x_2 + 42525.0;
    double temp_term2_z = -6081075.0 * temp_z_6 + 6081075.0 * temp_z_4 - 1403325.0 * temp_x_2 + 42525.0;

    double temp_dterm2_dx = -6081075.0 * temp_dx_6 + 6081075.0 * temp_dx_4 - 1403325.0 * temp_dx_2;
    double temp_dterm2_dy = -6081075.0 * temp_dy_6 + 6081075.0 * temp_dy_4 - 1403325.0 * temp_dy_2;
    double temp_dterm2_dz = -6081075.0 * temp_dz_6 + 6081075.0 * temp_dz_4 - 1403325.0 * temp_dz_2;

    double temp_miu1 = temp_y_3 * temp_term1_x + lambda_y0 * temp_term2_x;
    double temp_miu2 = temp_x_3 * temp_term1_y + lambda_x0 * temp_term2_y;
    double temp_miu3 = temp_z_3 * temp_term1_x + lambda_z0 * temp_term2_x;
    double temp_miu4 = temp_x_3 * temp_term1_z + lambda_x0 * temp_term2_z;
    double temp_miu5 = temp_z_3 * temp_term1_y + lambda_z0 * temp_term2_y;
    double temp_miu6 = temp_y_3 * temp_term1_z + lambda_y0 * temp_term2_z;

    double temp_dmiu1_dx = temp_y_3 * temp_dterm1_dx + lambda_y0 * temp_dterm2_dx;
    double temp_dmiu1_dy = temp_dy_3 * temp_term1_x + lambda * temp_term2_x;

    double temp_dmiu2_dx = temp_dx_3 * temp_term1_y + lambda * temp_term2_y;
    double temp_dmiu2_dy = temp_x_3 * temp_dterm1_dy + lambda_x0 * temp_dterm2_dy;

    double temp_dmiu3_dx = temp_z_3 * temp_dterm1_dx + lambda_z0 * temp_dterm2_dx;
    double temp_dmiu3_dz = temp_dz_3 * temp_term1_x + lambda * temp_term2_x;

    double temp_dmiu4_dx = temp_dx_3 * temp_term1_z + lambda * temp_term2_z;
    double temp_dmiu4_dz = temp_x_3 * temp_dterm1_dz + lambda_x0 * temp_dterm2_dz;

    double temp_dmiu5_dy = temp_z_3 * temp_dterm1_dy + lambda_z0 * temp_dterm2_dy;
    double temp_dmiu5_dz = temp_dz_3 * temp_term1_y + lambda * temp_term2_y;

    double temp_dmiu6_dy = temp_dy_3 * temp_term1_z + lambda * temp_term2_z;
    double temp_dmiu6_dz = temp_y_3 * temp_dterm1_dz + lambda_y0 * temp_dterm2_dz;

    double miu_9_5_1 = temp * temp_miu1;
    double miu_9_5_2 = temp * temp_miu2;
    double miu_9_5_3 = temp * temp_miu3;
    double miu_9_5_4 = temp * temp_miu4;
    double miu_9_5_5 = temp * temp_miu5;
    double miu_9_5_6 = temp * temp_miu6;

    //deriv consts
    double const_2_C2_x = 2.0 * C2 * x0; // 2 * C2 * x0
    double const_2_C2_y = 2.0 * C2 * y0;
    double const_2_C2_z = 2.0 * C2 * z0;

    // dmiu1 dx/dy/dz
    deriv[0] = temp * (temp_dmiu1_dx + const_2_C2_x * temp_miu1);
    deriv[1] = temp * (temp_dmiu1_dy + const_2_C2_y * temp_miu1);
    deriv[2] = miu_9_5_1 * const_2_C2_z;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * (temp_dmiu2_dx + const_2_C2_x * temp_miu2);
    deriv[4] = temp * (temp_dmiu2_dy + const_2_C2_y * temp_miu2);
    deriv[5] = miu_9_5_2 * const_2_C2_z;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * (temp_dmiu3_dx + const_2_C2_x * temp_miu3);
    deriv[7] = miu_9_5_3 * const_2_C2_y;
    deriv[8] = temp * (temp_dmiu3_dz + const_2_C2_z * temp_miu3);

    // dmiu4 dx/dy/dz
    deriv[9] = temp * (temp_dmiu4_dx + const_2_C2_x * temp_miu4);
    deriv[10] = miu_9_5_4 * const_2_C2_y;
    deriv[11] = temp * (temp_dmiu4_dz + const_2_C2_z * temp_miu4);

    // dmiu5 dx/dy/dz
    deriv[12] = miu_9_5_5 * const_2_C2_x;
    deriv[13] = temp * (temp_dmiu5_dy + const_2_C2_y * temp_miu5);
    deriv[14] = temp * (temp_dmiu5_dz + const_2_C2_z * temp_miu5);

    // dmiu6 dx/dy/dz
    deriv[15] = miu_9_5_6 * const_2_C2_x;
    deriv[16] = temp * (temp_dmiu6_dy + const_2_C2_y * temp_miu6);
    deriv[17] = temp * (temp_dmiu6_dz + const_2_C2_z * temp_miu6);

    value[0] = miu_9_5_1;
    value[1] = miu_9_5_2;
    value[2] = miu_9_5_3;
    value[3] = miu_9_5_4;
    value[4] = miu_9_5_5;
    value[5] = miu_9_5_6;
}


void calc_MCSH_9_6(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double lambda_sqr = lambda * lambda;
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

    double lambda_x0_5 = lambda_x0_4 * lambda_x0;
    double lambda_y0_5 = lambda_y0_4 * lambda_y0;
    double lambda_z0_5 = lambda_z0_4 * lambda_z0;

    double lambda_x0_6 = lambda_x0_5 * lambda_x0;
    double lambda_y0_6 = lambda_y0_5 * lambda_y0;
    double lambda_z0_6 = lambda_z0_5 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x_2 = lambda_x0_sqr + (1.0/(2.0*gamma));
    double temp_y_2 = lambda_y0_sqr + (1.0/(2.0*gamma));
    double temp_z_2 = lambda_z0_sqr + (1.0/(2.0*gamma));

    double temp_dx_2 = 2.0 * lambda_sqr * x0;
    double temp_dy_2 = 2.0 * lambda_sqr * y0;
    double temp_dz_2 = 2.0 * lambda_sqr * z0;

    double C4_1 = 3.0 / gamma, C4_2 = 3.0 / (4.0 * gamma * gamma);
    double temp_x_4 = lambda_x0_4 + C4_1 * lambda_x0_sqr + C4_2;
    double temp_y_4 = lambda_y0_4 + C4_1 * lambda_y0_sqr + C4_2;
    double temp_z_4 = lambda_z0_4 + C4_1 * lambda_z0_sqr + C4_2;

    double temp_dx_4 = lambda * (4.0 * lambda_x0_3 + 2.0 * C4_1 * lambda_x0);
    double temp_dy_4 = lambda * (4.0 * lambda_y0_3 + 2.0 * C4_1 * lambda_y0);
    double temp_dz_4 = lambda * (4.0 * lambda_z0_3 + 2.0 * C4_1 * lambda_z0);

    double C6_1 = 15.0 / (2.0 * gamma), C6_2 = 45.0 / (4.0 * gamma * gamma), C6_3 =  15.0 / (8.0 * gamma * gamma * gamma);
    double temp_x_6 = lambda_x0_6 + C6_1 * lambda_x0_4 + C6_2 * lambda_x0_sqr + C6_3;
    double temp_y_6 = lambda_y0_6 + C6_1 * lambda_y0_4 + C6_2 * lambda_y0_sqr + C6_3;
    double temp_z_6 = lambda_z0_6 + C6_1 * lambda_z0_4 + C6_2 * lambda_z0_sqr + C6_3;

    double temp_dx_6 = lambda * (6.0 * lambda_x0_5 + 4.0 * C6_1 * lambda_x0_3 + 2.0 * C6_2 * lambda_x0);
    double temp_dy_6 = lambda * (6.0 * lambda_y0_5 + 4.0 * C6_1 * lambda_y0_3 + 2.0 * C6_2 * lambda_y0);
    double temp_dz_6 = lambda * (6.0 * lambda_z0_5 + 4.0 * C6_1 * lambda_z0_3 + 2.0 * C6_2 * lambda_z0);


    double temp_term1_x = 34459425.0 * temp_x_6 - 30405375.0 * temp_x_4 + 6081075.0 * temp_x_2 - 155925.0;
    double temp_term1_y = 34459425.0 * temp_y_6 - 30405375.0 * temp_y_4 + 6081075.0 * temp_y_2 - 155925.0;
    double temp_term1_z = 34459425.0 * temp_z_6 - 30405375.0 * temp_z_4 + 6081075.0 * temp_z_2 - 155925.0;

    double temp_dterm1_dx = 34459425.0 * temp_dx_6 - 30405375.0 * temp_dx_4 + 6081075.0 * temp_dx_2;
    double temp_dterm1_dy = 34459425.0 * temp_dy_6 - 30405375.0 * temp_dy_4 + 6081075.0 * temp_dy_2;
    double temp_dterm1_dz = 34459425.0 * temp_dz_6 - 30405375.0 * temp_dz_4 + 6081075.0 * temp_dz_2;

    double temp_term2_x = -2027025.0 * temp_x_6 + 2027025.0 * temp_x_4 - 467775.0 * temp_x_2 + 14175.0;
    double temp_term2_y = -2027025.0 * temp_y_6 + 2027025.0 * temp_y_4 - 467775.0 * temp_x_2 + 14175.0;
    double temp_term2_z = -2027025.0 * temp_z_6 + 2027025.0 * temp_z_4 - 467775.0 * temp_x_2 + 14175.0;

    double temp_dterm2_dx = -2027025.0 * temp_dx_6 + 2027025.0 * temp_dx_4 - 467775.0 * temp_dx_2;
    double temp_dterm2_dy = -2027025.0 * temp_dy_6 + 2027025.0 * temp_dy_4 - 467775.0 * temp_dy_2;
    double temp_dterm2_dz = -2027025.0 * temp_dz_6 + 2027025.0 * temp_dz_4 - 467775.0 * temp_dz_2;

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

    double miu_9_6_1 = temp * lambda_z0 * temp_miu1;
    double miu_9_6_2 = temp * lambda_z0 * temp_miu2;
    double miu_9_6_3 = temp * lambda_y0 * temp_miu3;
    double miu_9_6_4 = temp * lambda_y0 * temp_miu4;
    double miu_9_6_5 = temp * lambda_x0 * temp_miu5;
    double miu_9_6_6 = temp * lambda_x0 * temp_miu6;

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

    // dmiu2 dx/dy/dz
    deriv[3] = temp * lambda_z0 * (temp_dmiu2_dx + const_2_C2_x * temp_miu2);
    deriv[4] = temp * lambda_z0 * (temp_dmiu2_dy + const_2_C2_y * temp_miu2);
    deriv[5] = temp * temp_miu2 * lambda * const_1_p_2_C2_z2;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * lambda_y0 * (temp_dmiu3_dx + const_2_C2_x * temp_miu3);
    deriv[7] = temp * temp_miu3 * lambda * const_1_p_2_C2_y2;
    deriv[8] = temp * lambda_y0 * (temp_dmiu3_dz + const_2_C2_z * temp_miu3);

    // dmiu4 dx/dy/dz
    deriv[9] = temp * lambda_y0 * (temp_dmiu4_dx + const_2_C2_x * temp_miu4);
    deriv[10] = temp * temp_miu4 * lambda * const_1_p_2_C2_y2;
    deriv[11] = temp * lambda_y0 * (temp_dmiu4_dz + const_2_C2_z * temp_miu4);

    // dmiu5 dx/dy/dz
    deriv[12] = temp * temp_miu5 * lambda * const_1_p_2_C2_x2;
    deriv[13] = temp * lambda_x0 * (temp_dmiu5_dy + const_2_C2_y * temp_miu5);
    deriv[14] = temp * lambda_x0 * (temp_dmiu5_dz + const_2_C2_z * temp_miu5);

    // dmiu6 dx/dy/dz
    deriv[15] = temp * temp_miu6 * lambda * const_1_p_2_C2_x2;
    deriv[16] = temp * lambda_x0 * (temp_dmiu6_dy + const_2_C2_y * temp_miu6);
    deriv[17] = temp * lambda_x0 * (temp_dmiu6_dz + const_2_C2_z * temp_miu6);

    value[0] = miu_9_6_1;
    value[1] = miu_9_6_2;
    value[2] = miu_9_6_3;
    value[3] = miu_9_6_4;
    value[4] = miu_9_6_5;
    value[5] = miu_9_6_6;
}

void calc_MCSH_9_7(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double lambda_sqr = lambda * lambda;
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

    double lambda_x0_4 = lambda_x0_3 * lambda_x0;
    double lambda_y0_4 = lambda_y0_3 * lambda_y0;
    double lambda_z0_4 = lambda_z0_3 * lambda_z0;

    double lambda_x0_5 = lambda_x0_4 * lambda_x0;
    double lambda_y0_5 = lambda_y0_4 * lambda_y0;
    double lambda_z0_5 = lambda_z0_4 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);

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

    double C4_1 = 3.0 / gamma, C4_2 = 3.0 / (4.0 * gamma * gamma);
    double temp_x_4 = lambda_x0_4 + C4_1 * lambda_x0_sqr + C4_2;
    double temp_y_4 = lambda_y0_4 + C4_1 * lambda_y0_sqr + C4_2;
    double temp_z_4 = lambda_z0_4 + C4_1 * lambda_z0_sqr + C4_2;

    double temp_dx_4 = lambda * (4.0 * lambda_x0_3 + 2.0 * C4_1 * lambda_x0);
    double temp_dy_4 = lambda * (4.0 * lambda_y0_3 + 2.0 * C4_1 * lambda_y0);
    double temp_dz_4 = lambda * (4.0 * lambda_z0_3 + 2.0 * C4_1 * lambda_z0);

    double C5_1 = 5.0 / gamma, C5_2 = 15.0 / (4.0 * gamma * gamma);
    double temp_x_5 = lambda_x0_5 + C5_1 * lambda_x0_3 + C5_2 * lambda_x0;
    double temp_y_5 = lambda_y0_5 + C5_1 * lambda_y0_3 + C5_2 * lambda_y0;
    double temp_z_5 = lambda_z0_5 + C5_1 * lambda_z0_3 + C5_2 * lambda_z0;

    double temp_dx_5 = lambda * (5.0 * lambda_x0_4 + 3.0 * C5_1 * lambda_x0_sqr + C5_2);
    double temp_dy_5 = lambda * (5.0 * lambda_y0_4 + 3.0 * C5_1 * lambda_y0_sqr + C5_2);
    double temp_dz_5 = lambda * (5.0 * lambda_z0_4 + 3.0 * C5_1 * lambda_z0_sqr + C5_2);


    double temp_term1_x = 34459425.0 * temp_x_5 - 20270250.0 * temp_x_3 + 2027025.0 * lambda_x0;
    double temp_term1_y = 34459425.0 * temp_y_5 - 20270250.0 * temp_y_3 + 2027025.0 * lambda_y0;
    double temp_term1_z = 34459425.0 * temp_z_5 - 20270250.0 * temp_z_3 + 2027025.0 * lambda_z0;

    double temp_dterm1_dx = 34459425.0 * temp_dx_5 - 20270250.0 * temp_dx_3 + 2027025.0 * lambda;
    double temp_dterm1_dy = 34459425.0 * temp_dy_5 - 20270250.0 * temp_dy_3 + 2027025.0 * lambda;
    double temp_dterm1_dz = 34459425.0 * temp_dz_5 - 20270250.0 * temp_dz_3 + 2027025.0 * lambda;

    double temp_term2_x = -12162150.0 * temp_x_5 + 8108100.0 * temp_x_3 - 935550.0 * lambda_x0;
    double temp_term2_y = -12162150.0 * temp_y_5 + 8108100.0 * temp_y_3 - 935550.0 * lambda_y0;
    double temp_term2_z = -12162150.0 * temp_z_5 + 8108100.0 * temp_z_3 - 935550.0 * lambda_z0;

    double temp_dterm2_dx = -12162150.0 * temp_dx_5 + 8108100.0 * temp_dx_3 - 935550.0 * lambda;
    double temp_dterm2_dy = -12162150.0 * temp_dy_5 + 8108100.0 * temp_dy_3 - 935550.0 * lambda;
    double temp_dterm2_dz = -12162150.0 * temp_dz_5 + 8108100.0 * temp_dz_3 - 935550.0 * lambda;

    double temp_term3_x = 405405.0 * temp_x_5 - 311850.0 * temp_x_3 + 425425.0 * lambda_x0;
    double temp_term3_y = 405405.0 * temp_y_5 - 311850.0 * temp_y_3 + 425425.0 * lambda_y0;
    double temp_term3_z = 405405.0 * temp_z_5 - 311850.0 * temp_z_3 + 425425.0 * lambda_z0;

    double temp_dterm3_dx = 405405.0 * temp_dx_5 - 311850.0 * temp_dx_3 + 425425.0 * lambda;
    double temp_dterm3_dy = 405405.0 * temp_dy_5 - 311850.0 * temp_dy_3 + 425425.0 * lambda;
    double temp_dterm3_dz = 405405.0 * temp_dz_5 - 311850.0 * temp_dz_3 + 425425.0 * lambda;

    double temp_miu1 = temp_y_4 * temp_term1_x + temp_y_2 * temp_term2_x + temp_term3_x;
    double temp_miu2 = temp_x_4 * temp_term1_y + temp_x_2 * temp_term2_y + temp_term3_y;
    double temp_miu3 = temp_z_4 * temp_term1_x + temp_z_2 * temp_term2_x + temp_term3_x;
    double temp_miu4 = temp_x_4 * temp_term1_z + temp_x_2 * temp_term2_z + temp_term3_z;
    double temp_miu5 = temp_z_4 * temp_term1_y + temp_z_2 * temp_term2_y + temp_term3_y;
    double temp_miu6 = temp_y_4 * temp_term1_z + temp_y_2 * temp_term2_z + temp_term3_z;

    double temp_dmiu1_dx = temp_y_4 * temp_dterm1_dx + temp_y_2 * temp_dterm2_dx + temp_dterm3_dx;
    double temp_dmiu1_dy = temp_dy_4 * temp_term1_x + temp_dy_2 * temp_term2_x;

    double temp_dmiu2_dx = temp_dx_4 * temp_term1_y + temp_dx_2 * temp_term2_y;
    double temp_dmiu2_dy = temp_x_4 * temp_dterm1_dy + temp_x_2 * temp_dterm2_dy + temp_dterm3_dy;

    double temp_dmiu3_dx = temp_z_4 * temp_dterm1_dx + temp_z_2 * temp_dterm2_dx + temp_dterm3_dx;
    double temp_dmiu3_dz = temp_dz_4 * temp_term1_x + temp_dz_2 * temp_term2_x;

    double temp_dmiu4_dx = temp_dx_4 * temp_term1_z + temp_dx_2 * temp_term2_z;
    double temp_dmiu4_dz = temp_x_4 * temp_dterm1_dz + temp_x_2 * temp_dterm2_dz + temp_dterm3_dz;

    double temp_dmiu5_dy = temp_z_4 * temp_dterm1_dy + temp_z_2 * temp_dterm2_dy + temp_dterm3_dy;
    double temp_dmiu5_dz = temp_dz_4 * temp_term1_y + temp_dz_2 * temp_term2_y;

    double temp_dmiu6_dy = temp_dy_4 * temp_term1_z + temp_dy_2 * temp_term2_z;
    double temp_dmiu6_dz = temp_y_4 * temp_dterm1_dz + temp_y_2 * temp_dterm2_dz + temp_dterm3_dz;

    double miu_9_7_1 = temp * temp_miu1;
    double miu_9_7_2 = temp * temp_miu2;
    double miu_9_7_3 = temp * temp_miu3;
    double miu_9_7_4 = temp * temp_miu4;
    double miu_9_7_5 = temp * temp_miu5;
    double miu_9_7_6 = temp * temp_miu6;

    //deriv consts
    double const_2_C2_x = 2.0 * C2 * x0; // 2 * C2 * x0
    double const_2_C2_y = 2.0 * C2 * y0;
    double const_2_C2_z = 2.0 * C2 * z0;

    // dmiu1 dx/dy/dz
    deriv[0] = temp * (temp_dmiu1_dx + const_2_C2_x * temp_miu1);
    deriv[1] = temp * (temp_dmiu1_dy + const_2_C2_y * temp_miu1);
    deriv[2] = miu_9_7_1 * const_2_C2_z;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * (temp_dmiu2_dx + const_2_C2_x * temp_miu2);
    deriv[4] = temp * (temp_dmiu2_dy + const_2_C2_y * temp_miu2);
    deriv[5] = miu_9_7_2 * const_2_C2_z;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * (temp_dmiu3_dx + const_2_C2_x * temp_miu3);
    deriv[7] = miu_9_7_3 * const_2_C2_y;
    deriv[8] = temp * (temp_dmiu3_dz + const_2_C2_z * temp_miu3);

    // dmiu4 dx/dy/dz
    deriv[9] = temp * (temp_dmiu4_dx + const_2_C2_x * temp_miu4);
    deriv[10] = miu_9_7_4 * const_2_C2_y;
    deriv[11] = temp * (temp_dmiu4_dz + const_2_C2_z * temp_miu4);

    // dmiu5 dx/dy/dz
    deriv[12] = miu_9_7_5 * const_2_C2_x;
    deriv[13] = temp * (temp_dmiu5_dy + const_2_C2_y * temp_miu5);
    deriv[14] = temp * (temp_dmiu5_dz + const_2_C2_z * temp_miu5);

    // dmiu6 dx/dy/dz
    deriv[15] = miu_9_7_6 * const_2_C2_x;
    deriv[16] = temp * (temp_dmiu6_dy + const_2_C2_y * temp_miu6);
    deriv[17] = temp * (temp_dmiu6_dz + const_2_C2_z * temp_miu6);

    value[0] = miu_9_7_1;
    value[1] = miu_9_7_2;
    value[2] = miu_9_7_3;
    value[3] = miu_9_7_4;
    value[4] = miu_9_7_5;
    value[5] = miu_9_7_6;
}

void calc_MCSH_9_8(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
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

    double lambda_x0_5 = lambda_x0_4 * lambda_x0;
    double lambda_y0_5 = lambda_y0_4 * lambda_y0;
    double lambda_z0_5 = lambda_z0_4 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);

    double temp = C1 * exp( C2 * r0_sqr);

    double C3_1 = 3.0 / (2.0 * gamma);
    double temp_x_3 = lambda_x0_3 + C3_1 * lambda_x0;
    double temp_y_3 = lambda_y0_3 + C3_1 * lambda_y0;
    double temp_z_3 = lambda_z0_3 + C3_1 * lambda_z0;

    double temp_dx_3 = lambda * (3.0 * lambda_x0_sqr + C3_1);
    double temp_dy_3 = lambda * (3.0 * lambda_y0_sqr + C3_1);
    double temp_dz_3 = lambda * (3.0 * lambda_z0_sqr + C3_1);

    double C5_1 = 5.0 / gamma, C5_2 = 15.0 / (4.0 * gamma * gamma);
    double temp_x_5 = lambda_x0_5 + C5_1 * lambda_x0_3 + C5_2 * lambda_x0;
    double temp_y_5 = lambda_y0_5 + C5_1 * lambda_y0_3 + C5_2 * lambda_y0;
    double temp_z_5 = lambda_z0_5 + C5_1 * lambda_z0_3 + C5_2 * lambda_z0;

    double temp_dx_5 = lambda * (5.0 * lambda_x0_4 + 3.0 * C5_1 * lambda_x0_sqr + C5_2);
    double temp_dy_5 = lambda * (5.0 * lambda_y0_4 + 3.0 * C5_1 * lambda_y0_sqr + C5_2);
    double temp_dz_5 = lambda * (5.0 * lambda_z0_4 + 3.0 * C5_1 * lambda_z0_sqr + C5_2);



    double temp_term1_x = 34459425.0 * temp_x_5 - 20270250.0 * temp_x_3 + 2027025.0 * lambda_x0;
    double temp_term1_y = 34459425.0 * temp_y_5 - 20270250.0 * temp_y_3 + 2027025.0 * lambda_y0;
    double temp_term1_z = 34459425.0 * temp_z_5 - 20270250.0 * temp_z_3 + 2027025.0 * lambda_z0;

    double temp_dterm1_dx = 34459425.0 * temp_dx_5 - 20270250.0 * temp_dx_3 + 2027025.0 * lambda;
    double temp_dterm1_dy = 34459425.0 * temp_dy_5 - 20270250.0 * temp_dy_3 + 2027025.0 * lambda;
    double temp_dterm1_dz = 34459425.0 * temp_dz_5 - 20270250.0 * temp_dz_3 + 2027025.0 * lambda;

    double temp_term2_x = -6081075.0 * temp_x_5 + 4054050.0 * temp_x_3 - 467775.0 * lambda_x0;
    double temp_term2_y = -6081075.0 * temp_y_5 + 4054050.0 * temp_y_3 - 467775.0 * lambda_x0;
    double temp_term2_z = -6081075.0 * temp_z_5 + 4054050.0 * temp_z_3 - 467775.0 * lambda_x0;

    double temp_dterm2_dx = -6081075.0 * temp_dx_5 + 4054050.0 * temp_dx_3 - 467775.0 * lambda;
    double temp_dterm2_dy = -6081075.0 * temp_dy_5 + 4054050.0 * temp_dy_3 - 467775.0 * lambda;
    double temp_dterm2_dz = -6081075.0 * temp_dz_5 + 4054050.0 * temp_dz_3 - 467775.0 * lambda;

    double temp_miu1 = temp_y_3 * temp_term1_x + lambda_y0 * temp_term2_x;
    double temp_miu2 = temp_x_3 * temp_term1_y + lambda_x0 * temp_term2_y;
    double temp_miu3 = temp_z_3 * temp_term1_x + lambda_z0 * temp_term2_x;
    double temp_miu4 = temp_x_3 * temp_term1_z + lambda_x0 * temp_term2_z;
    double temp_miu5 = temp_z_3 * temp_term1_y + lambda_z0 * temp_term2_y;
    double temp_miu6 = temp_y_3 * temp_term1_z + lambda_y0 * temp_term2_z;

    double temp_dmiu1_dx = temp_y_3 * temp_dterm1_dx + lambda_y0 * temp_dterm2_dx;
    double temp_dmiu1_dy = temp_dy_3 * temp_term1_x + lambda * temp_term2_x;

    double temp_dmiu2_dx = temp_dx_3 * temp_term1_y + lambda * temp_term2_y;
    double temp_dmiu2_dy = temp_x_3 * temp_dterm1_dy + lambda_x0 * temp_dterm2_dy;

    double temp_dmiu3_dx = temp_z_3 * temp_dterm1_dx + lambda_z0 * temp_dterm2_dx;
    double temp_dmiu3_dz = temp_dz_3 * temp_term1_x + lambda * temp_term2_x;

    double temp_dmiu4_dx = temp_dx_3 * temp_term1_z + lambda * temp_term2_z;
    double temp_dmiu4_dz = temp_x_3 * temp_dterm1_dz + lambda_x0 * temp_dterm2_dz;

    double temp_dmiu5_dy = temp_z_3 * temp_dterm1_dy + lambda_z0 * temp_dterm2_dy;
    double temp_dmiu5_dz = temp_dz_3 * temp_term1_y + lambda * temp_term2_y;

    double temp_dmiu6_dy = temp_dy_3 * temp_term1_z + lambda * temp_term2_z;
    double temp_dmiu6_dz = temp_y_3 * temp_dterm1_dz + lambda_y0 * temp_dterm2_dz;

    double miu_9_8_1 = temp * lambda_z0 * temp_miu1;
    double miu_9_8_2 = temp * lambda_z0 * temp_miu2;
    double miu_9_8_3 = temp * lambda_y0 * temp_miu3;
    double miu_9_8_4 = temp * lambda_y0 * temp_miu4;
    double miu_9_8_5 = temp * lambda_x0 * temp_miu5;
    double miu_9_8_6 = temp * lambda_x0 * temp_miu6;

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

    // dmiu2 dx/dy/dz
    deriv[3] = temp * lambda_z0 * (temp_dmiu2_dx + const_2_C2_x * temp_miu2);
    deriv[4] = temp * lambda_z0 * (temp_dmiu2_dy + const_2_C2_y * temp_miu2);
    deriv[5] = temp * temp_miu2 * lambda * const_1_p_2_C2_z2;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * lambda_y0 * (temp_dmiu3_dx + const_2_C2_x * temp_miu3);
    deriv[7] = temp * temp_miu3 * lambda * const_1_p_2_C2_y2;
    deriv[8] = temp * lambda_y0 * (temp_dmiu3_dz + const_2_C2_z * temp_miu3);

    // dmiu4 dx/dy/dz
    deriv[9] = temp * lambda_y0 * (temp_dmiu4_dx + const_2_C2_x * temp_miu4);
    deriv[10] = temp * temp_miu4 * lambda * const_1_p_2_C2_y2;
    deriv[11] = temp * lambda_y0 * (temp_dmiu4_dz + const_2_C2_z * temp_miu4);

    // dmiu5 dx/dy/dz
    deriv[12] = temp * temp_miu5 * lambda * const_1_p_2_C2_x2;
    deriv[13] = temp * lambda_x0 * (temp_dmiu5_dy + const_2_C2_y * temp_miu5);
    deriv[14] = temp * lambda_x0 * (temp_dmiu5_dz + const_2_C2_z * temp_miu5);

    // dmiu6 dx/dy/dz
    deriv[15] = temp * temp_miu6 * lambda * const_1_p_2_C2_x2;
    deriv[16] = temp * lambda_x0 * (temp_dmiu6_dy + const_2_C2_y * temp_miu6);
    deriv[17] = temp * lambda_x0 * (temp_dmiu6_dz + const_2_C2_z * temp_miu6);

    value[0] = miu_9_8_1;
    value[1] = miu_9_8_2;
    value[2] = miu_9_8_3;
    value[3] = miu_9_8_4;
    value[4] = miu_9_8_5;
    value[5] = miu_9_8_6;
}

void calc_MCSH_9_9(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double lambda_sqr = lambda * lambda;

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

    double lambda_x0_5 = lambda_x0_4 * lambda_x0;
    double lambda_y0_5 = lambda_y0_4 * lambda_y0;
    double lambda_z0_5 = lambda_z0_4 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);

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

    double C5_1 = 5.0 / gamma, C5_2 = 15.0 / (4.0 * gamma * gamma);
    double temp_x_5 = lambda_x0_5 + C5_1 * lambda_x0_3 + C5_2 * lambda_x0;
    double temp_y_5 = lambda_y0_5 + C5_1 * lambda_y0_3 + C5_2 * lambda_y0;
    double temp_z_5 = lambda_z0_5 + C5_1 * lambda_z0_3 + C5_2 * lambda_z0;

    double temp_dx_5 = lambda * (5.0 * lambda_x0_4 + 3.0 * C5_1 * lambda_x0_sqr + C5_2);
    double temp_dy_5 = lambda * (5.0 * lambda_y0_4 + 3.0 * C5_1 * lambda_y0_sqr + C5_2);
    double temp_dz_5 = lambda * (5.0 * lambda_z0_4 + 3.0 * C5_1 * lambda_z0_sqr + C5_2);


    double temp_term1_x = 34459425.0 * temp_x_5 - 20270250.0 * temp_x_3 + 2027025.0 * lambda_x0;
    double temp_term1_y = 34459425.0 * temp_y_5 - 20270250.0 * temp_y_3 + 2027025.0 * lambda_y0;
    double temp_term1_z = 34459425.0 * temp_z_5 - 20270250.0 * temp_z_3 + 2027025.0 * lambda_z0;

    double temp_dterm1_dx = 34459425.0 * temp_dx_5 - 20270250.0 * temp_dx_3 + 2027025.0 * lambda;
    double temp_dterm1_dy = 34459425.0 * temp_dy_5 - 20270250.0 * temp_dy_3 + 2027025.0 * lambda;
    double temp_dterm1_dz = 34459425.0 * temp_dz_5 - 20270250.0 * temp_dz_3 + 2027025.0 * lambda;

    double temp_term2_x = -2027025.0 * temp_x_5 + 1351350.0 * temp_x_3 - 155925.0 * lambda_x0;
    double temp_term2_y = -2027025.0 * temp_y_5 + 1351350.0 * temp_y_3 - 155925.0 * lambda_y0;
    double temp_term2_z = -2027025.0 * temp_z_5 + 1351350.0 * temp_z_3 - 155925.0 * lambda_z0;

    double temp_dterm2_dx = -2027025.0 * temp_dx_5 + 1351350.0 * temp_dx_3 - 155925.0 * lambda;
    double temp_dterm2_dy = -2027025.0 * temp_dy_5 + 1351350.0 * temp_dy_3 - 155925.0 * lambda;
    double temp_dterm2_dz = -2027025.0 * temp_dz_5 + 1351350.0 * temp_dz_3 - 155925.0 * lambda;

    double temp_term3_x = 135135.0 * temp_x_5 - 103950.0 * temp_x_3 + 14175.0 * lambda_x0;
    double temp_term3_y = 135135.0 * temp_y_5 - 103950.0 * temp_y_3 + 14175.0 * lambda_y0;
    double temp_term3_z = 135135.0 * temp_z_5 - 103950.0 * temp_z_3 + 14175.0 * lambda_z0;

    double temp_dterm3_dx = 135135.0 * temp_dx_5 - 103950.0 * temp_dx_3 + 14175.0 * lambda;
    double temp_dterm3_dy = 135135.0 * temp_dy_5 - 103950.0 * temp_dy_3 + 14175.0 * lambda;
    double temp_dterm3_dz = 135135.0 * temp_dz_5 - 103950.0 * temp_dz_3 + 14175.0 * lambda;

    double temp_miu1 = temp_y_2 * temp_z_2 * temp_term1_x + (temp_y_2 + temp_z_2) * temp_term2_x + temp_term3_x;
    double temp_miu2 = temp_x_2 * temp_z_2 * temp_term1_y + (temp_x_2 + temp_z_2) * temp_term2_y + temp_term3_y;
    double temp_miu3 = temp_x_2 * temp_y_2 * temp_term1_z + (temp_x_2 + temp_y_2) * temp_term2_z + temp_term3_z;

    double temp_dmiu1_dx = temp_y_2 * temp_z_2 * temp_dterm1_dx + (temp_y_2 + temp_z_2) * temp_dterm2_dx + temp_dterm3_dx;
    double temp_dmiu1_dy = temp_dy_2 * temp_z_2 * temp_term1_x + temp_dy_2 * temp_term2_x;
    double temp_dmiu1_dz = temp_y_2 * temp_dz_2 * temp_term1_x + temp_dz_2 * temp_term2_x;

    double temp_dmiu2_dx = temp_dx_2 * temp_z_2 * temp_term1_y + temp_dx_2 * temp_term2_y;
    double temp_dmiu2_dy = temp_x_2 * temp_z_2 * temp_dterm1_dy + (temp_x_2 + temp_z_2) * temp_dterm2_dy + temp_dterm3_dy;
    double temp_dmiu2_dz = temp_x_2 * temp_dz_2 * temp_term1_y + temp_dz_2 * temp_term2_y;

    double temp_dmiu3_dx = temp_dx_2 * temp_y_2 * temp_term1_z + temp_dx_2 * temp_term2_z;
    double temp_dmiu3_dy = temp_x_2 * temp_dy_2 * temp_term1_z + temp_dy_2 * temp_term2_z;
    double temp_dmiu3_dz = temp_x_2 * temp_y_2 * temp_dterm1_dz + (temp_x_2 + temp_y_2) * temp_dterm2_dz + temp_dterm3_dz;


    double miu_9_9_1 = temp * temp_miu1;
    double miu_9_9_2 = temp * temp_miu2;
    double miu_9_9_3 = temp * temp_miu3;

    //deriv consts
    double const_2_C2_x = 2.0 * C2 * x0; // 2 * C2 * x0
    double const_2_C2_y = 2.0 * C2 * y0;
    double const_2_C2_z = 2.0 * C2 * z0;

    // dmiu1 dx/dy/dz
    deriv[0] = temp * (temp_dmiu1_dx + const_2_C2_x * temp_miu1);
    deriv[1] = temp * (temp_dmiu1_dy + const_2_C2_y * temp_miu1);
    deriv[2] = temp * (temp_dmiu1_dz + const_2_C2_z * temp_miu1);

    // dmiu2 dx/dy/dz
    deriv[3] = temp * (temp_dmiu2_dx + const_2_C2_x * temp_miu2);
    deriv[4] = temp * (temp_dmiu2_dy + const_2_C2_y * temp_miu2);
    deriv[5] = temp * (temp_dmiu2_dz + const_2_C2_z * temp_miu2);

    // dmiu3 dx/dy/dz
    deriv[6] = temp * (temp_dmiu3_dx + const_2_C2_x * temp_miu3);
    deriv[7] = temp * (temp_dmiu3_dy + const_2_C2_y * temp_miu3);
    deriv[8] = temp * (temp_dmiu3_dz + const_2_C2_z * temp_miu3);

    value[0] = miu_9_9_1;
    value[1] = miu_9_9_2;
    value[2] = miu_9_9_3;
}


void calc_MCSH_9_10(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double lambda_sqr = lambda * lambda;
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

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x_2 = lambda_x0_sqr + (1.0/(2.0*gamma));
    double temp_y_2 = lambda_y0_sqr + (1.0/(2.0*gamma));
    double temp_z_2 = lambda_z0_sqr + (1.0/(2.0*gamma));

    double temp_dx_2 = 2.0 * lambda_sqr * x0;
    double temp_dy_2 = 2.0 * lambda_sqr * y0;
    double temp_dz_2 = 2.0 * lambda_sqr * z0;

    double C4_1 = 3.0 / gamma, C4_2 = 3.0 / (4.0 * gamma * gamma);
    double temp_x_4 = lambda_x0_4 + C4_1 * lambda_x0_sqr + C4_2;
    double temp_y_4 = lambda_y0_4 + C4_1 * lambda_y0_sqr + C4_2;
    double temp_z_4 = lambda_z0_4 + C4_1 * lambda_z0_sqr + C4_2;

    double temp_dx_4 = lambda * (4.0 * lambda_x0_3 + 2.0 * C4_1 * lambda_x0);
    double temp_dy_4 = lambda * (4.0 * lambda_y0_3 + 2.0 * C4_1 * lambda_y0);
    double temp_dz_4 = lambda * (4.0 * lambda_z0_3 + 2.0 * C4_1 * lambda_z0);



    double temp_term1_x = 34459425.0 * temp_x_4 - 12162150.0 * temp_x_2 + 405405.0;
    double temp_term1_y = 34459425.0 * temp_y_4 - 12162150.0 * temp_y_2 + 405405.0;

    double temp_dterm1_dx = 34459425.0 * temp_dx_4 - 12162150.0 * temp_dx_2;
    double temp_dterm1_dy = 34459425.0 * temp_dy_4 - 12162150.0 * temp_dy_2;

    double temp_term2_x = -12162150.0 * temp_x_4 + 4864860.0 * temp_x_2 - 187110.0;
    double temp_term2_y = -12162150.0 * temp_y_4 + 4864860.0 * temp_y_2 - 187110.0;

    double temp_dterm2_dx = -12162150.0 * temp_dx_4 + 4864860.0 * temp_dx_2;
    double temp_dterm2_dy = -12162150.0 * temp_dy_4 + 4864860.0 * temp_dy_2;

    double temp_term3_x = 405405.0 * temp_x_4 - 187110.0 * temp_x_2 + 8505.0;
    double temp_term3_y = 405405.0 * temp_y_4 - 187110.0 * temp_y_2 + 8505.0;

    double temp_dterm3_dx = 405405.0 * temp_dx_4 - 187110.0 * temp_dx_2;
    double temp_dterm3_dy = 405405.0 * temp_dy_4 - 187110.0 * temp_dy_2;


    double temp_miu1 = temp_y_4 * temp_term1_x + temp_y_2 * temp_term2_x + temp_term3_x;
    double temp_miu2 = temp_z_4 * temp_term1_x + temp_z_2 * temp_term2_x + temp_term3_x;
    double temp_miu3 = temp_z_4 * temp_term1_y + temp_z_2 * temp_term2_y + temp_term3_y;

    double temp_dmiu1_dx = temp_y_4 * temp_dterm1_dx + temp_y_2 * temp_dterm2_dx + temp_dterm3_dx;
    double temp_dmiu1_dy = temp_dy_4 * temp_term1_x + temp_dy_2 * temp_term2_x;

    double temp_dmiu2_dx = temp_z_4 * temp_dterm1_dx + temp_z_2 * temp_dterm2_dx + temp_dterm3_dx;
    double temp_dmiu2_dz = temp_dz_4 * temp_term1_x + temp_dz_2 * temp_term2_x;

    double temp_dmiu3_dy = temp_z_4 * temp_dterm1_dy + temp_z_2 * temp_dterm2_dy + temp_dterm3_dy;
    double temp_dmiu3_dz = temp_dz_4 * temp_term1_y + temp_dz_2 * temp_term2_y;


    double miu_9_10_1 = temp * lambda_z0 * temp_miu1;
    double miu_9_10_2 = temp * lambda_y0 * temp_miu2;
    double miu_9_10_3 = temp * lambda_x0 * temp_miu3;

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

    value[0] = miu_9_10_1;
    value[1] = miu_9_10_2;
    value[2] = miu_9_10_3;
}

void calc_MCSH_9_11(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double lambda_sqr = lambda * lambda;
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

    double lambda_x0_4 = lambda_x0_3 * lambda_x0;
    double lambda_y0_4 = lambda_y0_3 * lambda_y0;
    double lambda_z0_4 = lambda_z0_3 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);

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

    double C4_1 = 3.0 / gamma, C4_2 = 3.0 / (4.0 * gamma * gamma);
    double temp_x_4 = lambda_x0_4 + C4_1 * lambda_x0_sqr + C4_2;
    double temp_y_4 = lambda_y0_4 + C4_1 * lambda_y0_sqr + C4_2;
    double temp_z_4 = lambda_z0_4 + C4_1 * lambda_z0_sqr + C4_2;

    double temp_dx_4 = lambda * (4.0 * lambda_x0_3 + 2.0 * C4_1 * lambda_x0);
    double temp_dy_4 = lambda * (4.0 * lambda_y0_3 + 2.0 * C4_1 * lambda_y0);
    double temp_dz_4 = lambda * (4.0 * lambda_z0_3 + 2.0 * C4_1 * lambda_z0);


    double temp_term1_x = 34459425.0 * temp_x_4 - 12162150.0 * temp_x_2 + 405405.0;
    double temp_term1_y = 34459425.0 * temp_y_4 - 12162150.0 * temp_y_2 + 405405.0;
    double temp_term1_z = 34459425.0 * temp_z_4 - 12162150.0 * temp_z_2 + 405405.0;

    double temp_dterm1_dx = 34459425.0 * temp_dx_4 - 12162150.0 * temp_dx_2;
    double temp_dterm1_dy = 34459425.0 * temp_dy_4 - 12162150.0 * temp_dy_2;
    double temp_dterm1_dz = 34459425.0 * temp_dz_4 - 12162150.0 * temp_dz_2;

    double temp_term2_x = -2027025.0 * temp_x_4 + 810810.0 * temp_x_2 - 31185.0;
    double temp_term2_y = -2027025.0 * temp_y_4 + 810810.0 * temp_y_2 - 31185.0;
    double temp_term2_z = -2027025.0 * temp_z_4 + 810810.0 * temp_z_2 - 31185.0;

    double temp_dterm2_dx = -2027025.0 * temp_dx_4 + 810810.0 * temp_dx_2;
    double temp_dterm2_dy = -2027025.0 * temp_dy_4 + 810810.0 * temp_dy_2;
    double temp_dterm2_dz = -2027025.0 * temp_dz_4 + 810810.0 * temp_dz_2;

    double temp_term3_x = -6081075.0 * temp_x_4 + 2432430.0 * temp_x_2 - 93555.0;
    double temp_term3_y = -6081075.0 * temp_y_4 + 2432430.0 * temp_y_2 - 93555.0;
    double temp_term3_z = -6081075.0 * temp_z_4 + 2432430.0 * temp_z_2 - 93555.0;

    double temp_dterm3_dx = -6081075.0 * temp_dx_4 + 2432430.0 * temp_dx_2;
    double temp_dterm3_dy = -6081075.0 * temp_dy_4 + 2432430.0 * temp_dy_2;
    double temp_dterm3_dz = -6081075.0 * temp_dz_4 + 2432430.0 * temp_dz_2;

    double temp_term4_x = 405405.0 * temp_x_4 - 187110.0 * temp_x_2 + 8505.0;
    double temp_term4_y = 405405.0 * temp_y_4 - 187110.0 * temp_y_2 + 8505.0;
    double temp_term4_z = 405405.0 * temp_z_4 - 187110.0 * temp_z_2 + 8505.0;

    double temp_dterm4_dx = 405405.0 * temp_dx_4 - 187110.0 * temp_dx_2;
    double temp_dterm4_dy = 405405.0 * temp_dy_4 - 187110.0 * temp_dy_2;
    double temp_dterm4_dz = 405405.0 * temp_dz_4 - 187110.0 * temp_dz_2;


    double temp_miu1 = temp_y_3 * temp_z_2 * temp_term1_x + temp_y_3 * temp_term2_x + lambda_y0 * temp_z_2 * temp_term3_x + lambda_y0 * temp_term4_x;
    double temp_miu2 = temp_x_3 * temp_z_2 * temp_term1_y + temp_x_3 * temp_term2_y + lambda_x0 * temp_z_2 * temp_term3_y + lambda_x0 * temp_term4_y;
    double temp_miu3 = temp_z_3 * temp_y_2 * temp_term1_x + temp_z_3 * temp_term2_x + lambda_z0 * temp_y_2 * temp_term3_x + lambda_z0 * temp_term4_x;
    double temp_miu4 = temp_x_3 * temp_y_2 * temp_term1_z + temp_x_3 * temp_term2_z + lambda_x0 * temp_y_2 * temp_term3_z + lambda_x0 * temp_term4_z;
    double temp_miu5 = temp_z_3 * temp_x_2 * temp_term1_y + temp_z_3 * temp_term2_y + lambda_z0 * temp_x_2 * temp_term3_y + lambda_z0 * temp_term4_y;
    double temp_miu6 = temp_y_3 * temp_x_2 * temp_term1_z + temp_y_3 * temp_term2_z + lambda_y0 * temp_x_2 * temp_term3_z + lambda_y0 * temp_term4_z;

    double temp_dmiu1_dx = temp_y_3 * temp_z_2 * temp_dterm1_dx + temp_y_3 * temp_dterm2_dx + lambda_y0 * temp_z_2 * temp_dterm3_dx + lambda_y0 * temp_dterm4_dx;
    double temp_dmiu1_dy = temp_dy_3 * temp_z_2 * temp_term1_x + temp_dy_3 * temp_term2_x + lambda * temp_z_2 * temp_term3_x + lambda * temp_term4_x;
    double temp_dmiu1_dz = temp_y_3 * temp_dz_2 * temp_term1_x + lambda_y0 * temp_dz_2 * temp_term3_x;

    double temp_dmiu2_dx = temp_dx_3 * temp_z_2 * temp_term1_y + temp_dx_3 * temp_term2_y + lambda * temp_z_2 * temp_term3_y + lambda * temp_term4_y;
    double temp_dmiu2_dy = temp_x_3 * temp_z_2 * temp_dterm1_dy + temp_x_3 * temp_dterm2_dy + lambda_x0 * temp_z_2 * temp_dterm3_dy + lambda_x0 * temp_dterm4_dy;
    double temp_dmiu2_dz = temp_x_3 * temp_dz_2 * temp_term1_y + lambda_x0 * temp_dz_2 * temp_term3_y;

    double temp_dmiu3_dx = temp_z_3 * temp_y_2 * temp_dterm1_dx + temp_z_3 * temp_dterm2_dx + lambda_z0 * temp_y_2 * temp_dterm3_dx + lambda_z0 * temp_dterm4_dx;
    double temp_dmiu3_dy = temp_z_3 * temp_dy_2 * temp_term1_x + lambda_z0 * temp_dy_2 * temp_term3_x;
    double temp_dmiu3_dz = temp_dz_3 * temp_y_2 * temp_term1_x + temp_dz_3 * temp_term2_x + lambda * temp_y_2 * temp_term3_x + lambda * temp_term4_x;

    double temp_dmiu4_dx = temp_dx_3 * temp_y_2 * temp_term1_z + temp_dx_3 * temp_term2_z + lambda * temp_y_2 * temp_term3_z + lambda * temp_term4_z;
    double temp_dmiu4_dy = temp_x_3 * temp_dy_2 * temp_term1_z + lambda_x0 * temp_dy_2 * temp_term3_z;
    double temp_dmiu4_dz = temp_x_3 * temp_y_2 * temp_dterm1_dz + temp_x_3 * temp_dterm2_dz + lambda_x0 * temp_y_2 * temp_dterm3_dz + lambda_x0 * temp_dterm4_dz;

    double temp_dmiu5_dx = temp_z_3 * temp_dx_2 * temp_term1_y + lambda_z0 * temp_dx_2 * temp_term3_y;
    double temp_dmiu5_dy = temp_z_3 * temp_x_2 * temp_dterm1_dy + temp_z_3 * temp_dterm2_dy + lambda_z0 * temp_x_2 * temp_dterm3_dy + lambda_z0 * temp_dterm4_dy;
    double temp_dmiu5_dz = temp_dz_3 * temp_x_2 * temp_term1_y + temp_dz_3 * temp_term2_y + lambda * temp_x_2 * temp_term3_y + lambda * temp_term4_y;

    double temp_dmiu6_dx = temp_y_3 * temp_dx_2 * temp_term1_z + lambda_y0 * temp_dx_2 * temp_term3_z;
    double temp_dmiu6_dy = temp_dy_3 * temp_x_2 * temp_term1_z + temp_dy_3 * temp_term2_z + lambda * temp_x_2 * temp_term3_z + lambda * temp_term4_z;
    double temp_dmiu6_dz = temp_y_3 * temp_x_2 * temp_dterm1_dz + temp_y_3 * temp_dterm2_dz + lambda_y0 * temp_x_2 * temp_dterm3_dz + lambda_y0 * temp_dterm4_dz;

    double miu_9_11_1 = temp * temp_miu1;
    double miu_9_11_2 = temp * temp_miu2;
    double miu_9_11_3 = temp * temp_miu3;
    double miu_9_11_4 = temp * temp_miu4;
    double miu_9_11_5 = temp * temp_miu5;
    double miu_9_11_6 = temp * temp_miu6;

    //deriv consts
    double const_2_C2_x = 2.0 * C2 * x0; // 2 * C2 * x0
    double const_2_C2_y = 2.0 * C2 * y0;
    double const_2_C2_z = 2.0 * C2 * z0;

    // dmiu1 dx/dy/dz
    deriv[0] = temp * (temp_dmiu1_dx + const_2_C2_x * temp_miu1);
    deriv[1] = temp * (temp_dmiu1_dy + const_2_C2_y * temp_miu1);
    deriv[2] = temp * (temp_dmiu1_dz + const_2_C2_z * temp_miu1);

    // dmiu2 dx/dy/dz
    deriv[3] = temp * (temp_dmiu2_dx + const_2_C2_x * temp_miu2);
    deriv[4] = temp * (temp_dmiu2_dy + const_2_C2_y * temp_miu2);
    deriv[5] = temp * (temp_dmiu2_dz + const_2_C2_z * temp_miu2);

    // dmiu3 dx/dy/dz
    deriv[6] = temp * (temp_dmiu3_dx + const_2_C2_x * temp_miu3);
    deriv[7] = temp * (temp_dmiu3_dy + const_2_C2_y * temp_miu3);
    deriv[8] = temp * (temp_dmiu3_dz + const_2_C2_z * temp_miu3);

    // dmiu4 dx/dy/dz
    deriv[9] = temp * (temp_dmiu4_dx + const_2_C2_x * temp_miu4);
    deriv[10] = temp * (temp_dmiu4_dy + const_2_C2_y * temp_miu4);
    deriv[11] = temp * (temp_dmiu4_dz + const_2_C2_z * temp_miu4);

    // dmiu5 dx/dy/dz
    deriv[12] = temp * (temp_dmiu5_dx + const_2_C2_x * temp_miu5);
    deriv[13] = temp * (temp_dmiu5_dy + const_2_C2_y * temp_miu5);
    deriv[14] = temp * (temp_dmiu5_dz + const_2_C2_z * temp_miu5);

    // dmiu6 dx/dy/dz
    deriv[15] = temp * (temp_dmiu6_dx + const_2_C2_x * temp_miu6);
    deriv[16] = temp * (temp_dmiu6_dy + const_2_C2_y * temp_miu6);
    deriv[17] = temp * (temp_dmiu6_dz + const_2_C2_z * temp_miu6);

    value[0] = miu_9_11_1;
    value[1] = miu_9_11_2;
    value[2] = miu_9_11_3;
    value[3] = miu_9_11_4;
    value[4] = miu_9_11_5;
    value[5] = miu_9_11_6;
}


void calc_MCSH_9_12(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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

    double temp = C1 * exp( C2 * r0_sqr);

    double C3_1 = 3.0 / (2.0 * gamma);
    double temp_x_3 = lambda_x0_3 + C3_1 * lambda_x0;
    double temp_y_3 = lambda_y0_3 + C3_1 * lambda_y0;
    double temp_z_3 = lambda_z0_3 + C3_1 * lambda_z0;

    double temp_dx_3 = lambda * (3.0 * lambda_x0_sqr + C3_1);
    double temp_dy_3 = lambda * (3.0 * lambda_y0_sqr + C3_1);
    double temp_dz_3 = lambda * (3.0 * lambda_z0_sqr + C3_1);

    double t1 = 34459425.0 * temp_x_3 * temp_y_3 * temp_z_3;
    double t2 = -6081075.0 * temp_x_3 * temp_y_3 * lambda_z0;
    double t3 = -6081075.0 * temp_x_3 * lambda_y0 * temp_z_3;
    double t4 = 1216215.0 * temp_x_3 * lambda_y0 * lambda_z0;
    double t5 = -6081075.0 * lambda_x0 * temp_y_3 * temp_z_3;
    double t6 = 1216215.0 * lambda_x0 * temp_y_3 * lambda_z0;
    double t7 = 1216215.0 * lambda_x0 * lambda_y0 * temp_z_3;
    double t8 = -280665.0* lambda_x0 * lambda_y0 * lambda_z0;

    double sum_ts = t1 + t2 + t3 + t4 + t5 + t6 + t7 + t8;

    double sum_dtdx = temp_dx_3 * (34459425.0 * temp_y_3 * temp_z_3 - 6081075.0 * temp_y_3 * lambda_z0 - 6081075.0 * lambda_y0 * temp_z_3 + 1216215.0 * lambda_y0 * lambda_z0)
                       + lambda * (-6081075.0 * temp_y_3 * temp_z_3 + 1216215.0 * temp_y_3 * lambda_z0 + 1216215.0 * lambda_y0 * temp_z_3 - 280665.0 * lambda_y0 * lambda_z0);
    double sum_dtdy = temp_dy_3 * (34459425.0 * temp_x_3 * temp_z_3 - 6081075.0 * temp_x_3 * lambda_z0 - 6081075.0 * lambda_x0 * temp_z_3 + 1216215.0 * lambda_x0 * lambda_z0)
                       + lambda * (-6081075.0 * temp_x_3 * temp_z_3 + 1216215.0 * temp_x_3 * lambda_z0 + 1216215.0 * lambda_x0 * temp_z_3 - 280665.0 * lambda_x0 * lambda_z0);
    double sum_dtdz = temp_dz_3 * (34459425.0 * temp_x_3 * temp_y_3 - 6081075.0 * temp_x_3 * lambda_y0 - 6081075.0 * lambda_x0 * temp_y_3 + 1216215.0 * lambda_x0 * lambda_y0)
                       + lambda * (-6081075.0 * temp_x_3 * temp_y_3 + 1216215.0 * temp_x_3 * lambda_y0 + 1216215.0 * lambda_x0 * temp_y_3 - 280665.0 * lambda_x0 * lambda_y0);


    double m_9_12 = temp * sum_ts;

    deriv[0] = temp * (2.0 * C2 * x0 * sum_ts + sum_dtdx);
    deriv[1] = temp * (2.0 * C2 * y0 * sum_ts + sum_dtdy);
    deriv[2] = temp * (2.0 * C2 * z0 * sum_ts + sum_dtdz);

    value[0] = m_9_12;
}

















































void calc_MCSH_0_1_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double m_0_1 = C1 * exp( C2 * r0_sqr);

    value[0] = m_0_1;
}

void calc_MCSH_1_1_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double lambda = calc_lambda(alpha, beta);

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = lambda * x0;
    double temp_y = lambda * y0;
    double temp_z = lambda * z0;

    double miu_1_1_1 = temp * temp_x;
    double miu_1_1_2 = temp * temp_y;
    double miu_1_1_3 = temp * temp_z;

    value[0] = miu_1_1_1;
    value[1] = miu_1_1_2;
    value[2] = miu_1_1_3;

}

void calc_MCSH_2_1_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

    double lambda_x0 = x0 * lambda;
    double lambda_y0 = y0 * lambda;
    double lambda_z0 = z0 * lambda;

    double lambda_x0_sqr = lambda_x0 * lambda_x0;
    double lambda_y0_sqr = lambda_y0 * lambda_y0;
    double lambda_z0_sqr = lambda_z0 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = (3.0/(2.0*gamma)) - 1.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 3.0 * lambda_x0_sqr + C3;
    double temp_y = 3.0 * lambda_y0_sqr + C3;
    double temp_z = 3.0 * lambda_z0_sqr + C3;

    double miu_2_1_1 = temp * temp_x;
    double miu_2_1_2 = temp * temp_y;
    double miu_2_1_3 = temp * temp_z;

    value[0] = miu_2_1_1;
    value[1] = miu_2_1_2;
    value[2] = miu_2_1_3;
}

void calc_MCSH_2_2_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

    double temp = C1 * exp( C2 * r0_sqr) * lambda * lambda * 3.0;

    double miu_2_2_1 = temp * x0 * y0;
    double miu_2_2_2 = temp * x0 * z0;
    double miu_2_2_3 = temp * y0 * z0;

    value[0] = miu_2_2_1;
    value[1] = miu_2_2_2;
    value[2] = miu_2_2_3;
}

void calc_MCSH_3_1_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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
    double C3 = (45.0/(2.0*gamma)) - 9.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 15.0 * lambda_x0_3 + C3 * lambda_x0;
    double temp_y = 15.0 * lambda_y0_3 + C3 * lambda_y0;
    double temp_z = 15.0 * lambda_z0_3 + C3 * lambda_z0;

    double miu_3_1_1 = temp * temp_x;
    double miu_3_1_2 = temp * temp_y;
    double miu_3_1_3 = temp * temp_z;

    value[0] = miu_3_1_1;
    value[1] = miu_3_1_2;
    value[2] = miu_3_1_3;
}

void calc_MCSH_3_2_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

    double lambda_x0 = x0 * lambda;
    double lambda_y0 = y0 * lambda;
    double lambda_z0 = z0 * lambda;

    double lambda_x0_sqr = lambda_x0 * lambda_x0;
    double lambda_y0_sqr = lambda_y0 * lambda_y0;
    double lambda_z0_sqr = lambda_z0 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = (15.0/(2.0*gamma)) - 3.0;

    double temp = C1 * exp( C2 * r0_sqr);// * lambda;

    double temp_x = 15.0 * lambda_x0_sqr + C3;
    double temp_y = 15.0 * lambda_y0_sqr + C3;
    double temp_z = 15.0 * lambda_z0_sqr + C3;

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

void calc_MCSH_3_3_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double lambda = calc_lambda(alpha, beta);

    double temp =  C1 * exp( C2 * r0_sqr) * lambda * lambda * lambda * 15.0;
    double m_3_3 = temp * x0 * y0 * z0;

    value[0] = m_3_3;
}

void calc_MCSH_4_1_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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
    double C3 = (315.0/gamma) - 90.0;
    double C4 = (315.0/(4.0*gamma*gamma)) - (45.0/gamma) + 9.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 105.0 * lambda_x0_4 + C3 * lambda_x0_sqr + C4;
    double temp_y = 105.0 * lambda_y0_4 + C3 * lambda_y0_sqr + C4;
    double temp_z = 105.0 * lambda_z0_4 + C3 * lambda_z0_sqr + C4;

    double miu_4_1_1 = temp * temp_x;
    double miu_4_1_2 = temp * temp_y;
    double miu_4_1_3 = temp * temp_z;

    value[0] = miu_4_1_1;
    value[1] = miu_4_1_2;
    value[2] = miu_4_1_3;
}

void calc_MCSH_4_2_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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

    double C3 = (315.0/(2.0*gamma)) - 45.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 105.0 * lambda_x0_3 + C3 * lambda_x0;
    double temp_y = 105.0 * lambda_y0_3 + C3 * lambda_y0;
    double temp_z = 105.0 * lambda_z0_3 + C3 * lambda_z0;

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


void calc_MCSH_4_3_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    double lambda_sqr = lambda * lambda;
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

    double temp_term1_x = 105.0 * temp_x_2 - 15.0;
    double temp_term1_y = 105.0 * temp_y_2 - 15.0;
    // double temp_term1_z = 105.0 * temp_z_2 - 15.0;

    double temp_term2_x = -15.0 * temp_x_2 + 3.0;
    double temp_term2_y = -15.0 * temp_y_2 + 3.0;
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

void calc_MCSH_4_4_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

    double lambda_x0 = lambda * x0;
    double lambda_y0 = lambda * y0;
    double lambda_z0 = lambda * z0;

    double lambda_x0_sqr = lambda_x0 * lambda_x0;
    double lambda_y0_sqr = lambda_y0 * lambda_y0;
    double lambda_z0_sqr = lambda_z0 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = (105.0 / (2.0 * gamma)) - 15.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 105.0 * lambda_x0_sqr + C3;
    double temp_y = 105.0 * lambda_y0_sqr + C3;
    double temp_z = 105.0 * lambda_z0_sqr + C3;

    double miu_4_4_1 = temp * lambda_y0 * lambda_z0 * temp_x;
    double miu_4_4_2 = temp * lambda_x0 * lambda_z0 * temp_y;
    double miu_4_4_3 = temp * lambda_x0 * lambda_y0 * temp_z;

    value[0] = miu_4_4_1;
    value[1] = miu_4_4_2;
    value[2] = miu_4_4_3;

}


void calc_MCSH_5_1_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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
    double C3 = (4725.0 / gamma) - 1050.0;
    double C4 = (14175.0 / (4.0*gamma*gamma)) - (1575.0 / gamma) + 225.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 945.0 * lambda_x0_5 + C3 * lambda_x0_3 + C4 * lambda_x0;
    double temp_y = 945.0 * lambda_y0_5 + C3 * lambda_y0_3 + C4 * lambda_y0;
    double temp_z = 945.0 * lambda_z0_5 + C3 * lambda_z0_3 + C4 * lambda_z0;

    double miu_5_1_1 = temp * temp_x;
    double miu_5_1_2 = temp * temp_y;
    double miu_5_1_3 = temp * temp_z;

    value[0] = miu_5_1_1;
    value[1] = miu_5_1_2;
    value[2] = miu_5_1_3;
}

void calc_MCSH_5_2_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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
    double C3 = (945.0 * 3.0 / gamma) - 630.0;
    double C4 = (945.0 * 3.0 / (4.0 * gamma * gamma)) - (630.0 / (2.0 * gamma)) + 45.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 945.0 * lambda_x0_4 + C3 * lambda_x0_sqr + C4;
    double temp_y = 945.0 * lambda_y0_4 + C3 * lambda_y0_sqr + C4;
    double temp_z = 945.0 * lambda_z0_4 + C3 * lambda_z0_sqr + C4;

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

void calc_MCSH_5_3_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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

    double temp_term1_x = 945.0 * temp_x_3 - 315.0 * lambda_x0;
    double temp_term1_y = 945.0 * temp_y_3 - 315.0 * lambda_y0;
    double temp_term1_z = 945.0 * temp_z_3 - 315.0 * lambda_z0;

    double temp_term2_x = -105.0 * temp_x_3 + 45.0 * lambda_x0;
    double temp_term2_y = -105.0 * temp_y_3 + 45.0 * lambda_y0;
    double temp_term2_z = -105.0 * temp_z_3 + 45.0 * lambda_z0;

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

void calc_MCSH_5_4_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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
    double C3 = (2835.0 / (2.0 * gamma)) - 315.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 945.0 * lambda_x0_3 + C3 * lambda_x0;
    double temp_y = 945.0 * lambda_y0_3 + C3 * lambda_y0;
    double temp_z = 945.0 * lambda_z0_3 + C3 * lambda_z0;

    double miu_5_4_1 = temp * lambda_y0 * lambda_z0 * temp_x;
    double miu_5_4_2 = temp * lambda_x0 * lambda_z0 * temp_y;
    double miu_5_4_3 = temp * lambda_x0 * lambda_y0 * temp_z;

    value[0] = miu_5_4_1;
    value[1] = miu_5_4_2;
    value[2] = miu_5_4_3;

}

void calc_MCSH_5_5_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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

    double temp_miu1 = 945.0 * temp_x_2 * temp_y_2 - 105.0 * (temp_x_2 + temp_y_2) + 15.0;
    double temp_miu2 = 945.0 * temp_x_2 * temp_z_2 - 105.0 * (temp_x_2 + temp_z_2) + 15.0;
    double temp_miu3 = 945.0 * temp_y_2 * temp_z_2 - 105.0 * (temp_y_2 + temp_z_2) + 15.0;

    double miu_5_5_1 = temp * lambda_z0 * temp_miu1;
    double miu_5_5_2 = temp * lambda_y0 * temp_miu2;
    double miu_5_5_3 = temp * lambda_x0 * temp_miu3;

    value[0] = miu_5_5_1;
    value[1] = miu_5_5_2;
    value[2] = miu_5_5_3;
}

void calc_MCSH_6_1_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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

    double lambda_x0_6 = lambda_x0_5 * lambda_x0;
    double lambda_y0_6 = lambda_y0_5 * lambda_y0;
    double lambda_z0_6 = lambda_z0_5 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = (10395.0 * 15.0 / (2.0 * gamma)) - 14175.0;
    double C4 = (10395.0 * 45.0 / (4.0 * gamma * gamma)) - (14175.0 * 3.0 / gamma) + 4725.0;
    double C5 = (10395.0 * 15.0 / (8.0 * gamma * gamma * gamma)) - (14175.0 * 3.0 / (4.0 * gamma * gamma)) + (4725.0 / (2.0 * gamma)) - 225.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 10395.0 * lambda_x0_6 + C3 * lambda_x0_4 + C4 * lambda_x0_sqr + C5;
    double temp_y = 10395.0 * lambda_y0_6 + C3 * lambda_y0_4 + C4 * lambda_y0_sqr + C5;
    double temp_z = 10395.0 * lambda_z0_6 + C3 * lambda_z0_4 + C4 * lambda_z0_sqr + C5;

    double miu_6_1_1 = temp * temp_x;
    double miu_6_1_2 = temp * temp_y;
    double miu_6_1_3 = temp * temp_z;

    value[0] = miu_6_1_1;
    value[1] = miu_6_1_2;
    value[2] = miu_6_1_3;
}

void calc_MCSH_6_2_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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

    double lambda_x0_5 = lambda_x0_4 * lambda_x0;
    double lambda_y0_5 = lambda_y0_4 * lambda_y0;
    double lambda_z0_5 = lambda_z0_4 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = (10395.0 * 5.0 / gamma) - 9450.0;
    double C4 = (10395.0 * 15.0 / ( 4.0 * gamma * gamma)) - (9450.0 * 3.0 / (2.0 * gamma)) + 1575.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 10395.0 * lambda_x0_5 + C3 * lambda_x0_3 + C4 * lambda_x0;
    double temp_y = 10395.0 * lambda_y0_5 + C3 * lambda_y0_3 + C4 * lambda_y0;
    double temp_z = 10395.0 * lambda_z0_5 + C3 * lambda_z0_3 + C4 * lambda_z0;

    double miu_6_2_1 = temp * lambda_y0 * temp_x;
    double miu_6_2_2 = temp * lambda_x0 * temp_y;
    double miu_6_2_3 = temp * lambda_z0 * temp_x;
    double miu_6_2_4 = temp * lambda_x0 * temp_z;
    double miu_6_2_5 = temp * lambda_z0 * temp_y;
    double miu_6_2_6 = temp * lambda_y0 * temp_z;

    value[0] = miu_6_2_1;
    value[1] = miu_6_2_2;
    value[2] = miu_6_2_3;
    value[3] = miu_6_2_4;
    value[4] = miu_6_2_5;
    value[5] = miu_6_2_6;
}

void calc_MCSH_6_3_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x_2 = lambda_x0_sqr + (1.0/(2.0*gamma));
    double temp_y_2 = lambda_y0_sqr + (1.0/(2.0*gamma));
    double temp_z_2 = lambda_z0_sqr + (1.0/(2.0*gamma));

    double C4_1 = 3.0 / gamma, C4_2 = 3.0 / (4.0 * gamma * gamma);
    double temp_x_4 = lambda_x0_4 + C4_1 * lambda_x0_sqr + C4_2;
    double temp_y_4 = lambda_y0_4 + C4_1 * lambda_y0_sqr + C4_2;
    double temp_z_4 = lambda_z0_4 + C4_1 * lambda_z0_sqr + C4_2;

    double temp_term1_x = 10395.0 * temp_x_4 - 5670.0 * temp_x_2 + 315.0;
    double temp_term1_y = 10395.0 * temp_y_4 - 5670.0 * temp_y_2 + 315.0;
    double temp_term1_z = 10395.0 * temp_z_4 - 5670.0 * temp_z_2 + 315.0;

    double temp_term2_x = -945.0 * temp_x_4 + 630.0 * temp_x_2 - 45.0;
    double temp_term2_y = -945.0 * temp_y_4 + 630.0 * temp_y_2 - 45.0;
    double temp_term2_z = -945.0 * temp_z_4 + 630.0 * temp_z_2 - 45.0;

    double temp_miu1 = temp_y_2 * temp_term1_x + temp_term2_x;
    double temp_miu2 = temp_x_2 * temp_term1_y + temp_term2_y;
    double temp_miu3 = temp_z_2 * temp_term1_x + temp_term2_x;
    double temp_miu4 = temp_x_2 * temp_term1_z + temp_term2_z;
    double temp_miu5 = temp_z_2 * temp_term1_y + temp_term2_y;
    double temp_miu6 = temp_y_2 * temp_term1_z + temp_term2_z;

    double miu_6_3_1 = temp * temp_miu1;
    double miu_6_3_2 = temp * temp_miu2;
    double miu_6_3_3 = temp * temp_miu3;
    double miu_6_3_4 = temp * temp_miu4;
    double miu_6_3_5 = temp * temp_miu5;
    double miu_6_3_6 = temp * temp_miu6;

    value[0] = miu_6_3_1;
    value[1] = miu_6_3_2;
    value[2] = miu_6_3_3;
    value[3] = miu_6_3_4;
    value[4] = miu_6_3_5;
    value[5] = miu_6_3_6;
}


void calc_MCSH_6_4_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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
    double C3 = (10395.0 * 3.0 / gamma) - 5670.0;
    double C4 = (10395.0 * 3.0 / (4.0 * gamma * gamma)) - (5670.0 / (2.0 * gamma)) + 315.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 10395.0 * lambda_x0_4 + C3 * lambda_x0_sqr + C4;
    double temp_y = 10395.0 * lambda_y0_4 + C3 * lambda_y0_sqr + C4;
    double temp_z = 10395.0 * lambda_z0_4 + C3 * lambda_z0_sqr + C4;

    double miu_6_4_1 = temp * lambda_y0 * lambda_z0 * temp_x;
    double miu_6_4_2 = temp * lambda_x0 * lambda_z0 * temp_y;
    double miu_6_4_3 = temp * lambda_x0 * lambda_y0 * temp_z;

    value[0] = miu_6_4_1;
    value[1] = miu_6_4_2;
    value[2] = miu_6_4_3;
}

void calc_MCSH_6_5_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    // double lambda_sqr = lambda * lambda;

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

    double temp = C1 * exp( C2 * r0_sqr);

    double C3_1 = 3.0 / (2.0 * gamma);
    double temp_x_3 = lambda_x0_3 + C3_1 * lambda_x0;
    double temp_y_3 = lambda_y0_3 + C3_1 * lambda_y0;
    double temp_z_3 = lambda_z0_3 + C3_1 * lambda_z0;

    double temp_miu1 = 10395.0 * temp_x_3 * temp_y_3 - 2835.0 * (temp_x_3 * lambda_y0 + lambda_x0 * temp_y_3) + 945.0 * lambda_x0 * lambda_y0;
    double temp_miu2 = 10395.0 * temp_x_3 * temp_z_3 - 2835.0 * (temp_x_3 * lambda_z0 + lambda_x0 * temp_z_3) + 945.0 * lambda_x0 * lambda_z0;
    double temp_miu3 = 10395.0 * temp_y_3 * temp_z_3 - 2835.0 * (temp_y_3 * lambda_z0 + lambda_y0 * temp_z_3) + 945.0 * lambda_y0 * lambda_z0;

    double miu_6_5_1 = temp * temp_miu1;
    double miu_6_5_2 = temp * temp_miu2;
    double miu_6_5_3 = temp * temp_miu3;

    value[0] = miu_6_5_1;
    value[1] = miu_6_5_2;
    value[2] = miu_6_5_3;
}

void calc_MCSH_6_6_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x_2 = lambda_x0_sqr + (1.0/(2.0*gamma));
    double temp_y_2 = lambda_y0_sqr + (1.0/(2.0*gamma));
    double temp_z_2 = lambda_z0_sqr + (1.0/(2.0*gamma));

    double C3_1 = 3.0 / (2.0 * gamma);
    double temp_x_3 = lambda_x0_3 + C3_1 * lambda_x0;
    double temp_y_3 = lambda_y0_3 + C3_1 * lambda_y0;
    double temp_z_3 = lambda_z0_3 + C3_1 * lambda_z0;

    double temp_term1_x = 10395.0 * temp_x_3 - 2835.0 * lambda_x0;
    double temp_term1_y = 10395.0 * temp_y_3 - 2835.0 * lambda_y0;
    double temp_term1_z = 10395.0 * temp_z_3 - 2835.0 * lambda_z0;

    double temp_term2_x = -945.0 * temp_x_3 + 315.0 * lambda_x0;
    double temp_term2_y = -945.0 * temp_y_3 + 315.0 * lambda_y0;
    double temp_term2_z = -945.0 * temp_z_3 + 315.0 * lambda_z0;

    double temp_miu1 = temp_y_2 * temp_term1_x + temp_term2_x;
    double temp_miu2 = temp_x_2 * temp_term1_y + temp_term2_y;
    double temp_miu3 = temp_z_2 * temp_term1_x + temp_term2_x;
    double temp_miu4 = temp_x_2 * temp_term1_z + temp_term2_z;
    double temp_miu5 = temp_z_2 * temp_term1_y + temp_term2_y;
    double temp_miu6 = temp_y_2 * temp_term1_z + temp_term2_z;

    double miu_6_6_1 = temp * lambda_z0 * temp_miu1;
    double miu_6_6_2 = temp * lambda_z0 * temp_miu2;
    double miu_6_6_3 = temp * lambda_y0 * temp_miu3;
    double miu_6_6_4 = temp * lambda_y0 * temp_miu4;
    double miu_6_6_5 = temp * lambda_x0 * temp_miu5;
    double miu_6_6_6 = temp * lambda_x0 * temp_miu6;

    value[0] = miu_6_6_1;
    value[1] = miu_6_6_2;
    value[2] = miu_6_6_3;
    value[3] = miu_6_6_4;
    value[4] = miu_6_6_5;
    value[5] = miu_6_6_6;
}


void calc_MCSH_6_7_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

    double lambda_x0 = x0 * lambda;
    double lambda_y0 = y0 * lambda;
    double lambda_z0 = z0 * lambda;

    double lambda_x0_sqr = lambda_x0 * lambda_x0;
    double lambda_y0_sqr = lambda_y0 * lambda_y0;
    double lambda_z0_sqr = lambda_z0 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);

    double temp_x = lambda_x0_sqr + (1.0/(2.0*gamma));
    double temp_y = lambda_y0_sqr + (1.0/(2.0*gamma));
    double temp_z = lambda_z0_sqr + (1.0/(2.0*gamma));

    double t1 = 10395.0 * temp_x * temp_y * temp_z;
    double t2 = -945.0 * temp_x * temp_y;
    double t3 = -945.0 * temp_x * temp_z;
    double t4 = -945.0 * temp_y * temp_z;
    double t5 = 105.0 * temp_x;
    double t6 = 105.0 * temp_y;
    double t7 = 105.0 * temp_z;
    double t8 = -15.0;

    double sum_ts = t1 + t2 + t3 + t4 + t5 + t6 + t7 + t8;

    double temp =  C1 * exp( C2 * r0_sqr);
    double m_6_7 = temp * sum_ts;

    value[0] = m_6_7;
}

void calc_MCSH_7_1_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
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

    double lambda_x0_6 = lambda_x0_5 * lambda_x0;
    double lambda_y0_6 = lambda_y0_5 * lambda_y0;
    double lambda_z0_6 = lambda_z0_5 * lambda_z0;

    double lambda_x0_7 = lambda_x0_6 * lambda_x0;
    double lambda_y0_7 = lambda_y0_6 * lambda_y0;
    double lambda_z0_7 = lambda_z0_6 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = (135135.0 * 21.0 / (2.0 * gamma)) - 218295.0;
    double C4 = (135135.0 * 105.0 / (4.0 * gamma * gamma)) - (218295.0 * 5.0 / gamma) + 99225.0;
    double C5 = (135135.0 * 105.0 / (8.0 * gamma * gamma * gamma)) - (218295.0 * 15.0 / (4.0 * gamma * gamma)) + (99225.0 * 3.0 / (2.0 * gamma)) - 11025.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 135135.0 * lambda_x0_7 + C3 * lambda_x0_5 + C4 * lambda_x0_3 + C5 * lambda_x0;
    double temp_y = 135135.0 * lambda_y0_7 + C3 * lambda_y0_5 + C4 * lambda_y0_3 + C5 * lambda_y0;
    double temp_z = 135135.0 * lambda_z0_7 + C3 * lambda_z0_5 + C4 * lambda_z0_3 + C5 * lambda_z0;

    double miu_7_1_1 = temp * temp_x;
    double miu_7_1_2 = temp * temp_y;
    double miu_7_1_3 = temp * temp_z;

    value[0] = miu_7_1_1;
    value[1] = miu_7_1_2;
    value[2] = miu_7_1_3;
}

void calc_MCSH_7_2_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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

    double lambda_x0_5 = lambda_x0_4 * lambda_x0;
    double lambda_y0_5 = lambda_y0_4 * lambda_y0;
    double lambda_z0_5 = lambda_z0_4 * lambda_z0;

    double lambda_x0_6 = lambda_x0_5 * lambda_x0;
    double lambda_y0_6 = lambda_y0_5 * lambda_y0;
    double lambda_z0_6 = lambda_z0_5 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = (135135.0 * 15.0 / (2.0 * gamma)) - 155925.0;
    double C4 = (135135.0 * 45.0 / (4.0 * gamma * gamma)) - (155925.0 * 3.0 / gamma) + 42525.0;
    double C5 = (135135.0 * 15.0 / (8.0 * gamma * gamma * gamma)) - (155925.0 * 3.0 / (4.0 * gamma * gamma)) + (42525.0 / (2.0 * gamma)) - 1575.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 135135.0 * lambda_x0_6 + C3 * lambda_x0_4 + C4 * lambda_x0_sqr + C5;
    double temp_y = 135135.0 * lambda_y0_6 + C3 * lambda_y0_4 + C4 * lambda_y0_sqr + C5;
    double temp_z = 135135.0 * lambda_z0_6 + C3 * lambda_z0_4 + C4 * lambda_z0_sqr + C5;

    double miu_7_2_1 = temp * lambda_y0 * temp_x;
    double miu_7_2_2 = temp * lambda_x0 * temp_y;
    double miu_7_2_3 = temp * lambda_z0 * temp_x;
    double miu_7_2_4 = temp * lambda_x0 * temp_z;
    double miu_7_2_5 = temp * lambda_z0 * temp_y;
    double miu_7_2_6 = temp * lambda_y0 * temp_z;

    value[0] = miu_7_2_1;
    value[1] = miu_7_2_2;
    value[2] = miu_7_2_3;
    value[3] = miu_7_2_4;
    value[4] = miu_7_2_5;
    value[5] = miu_7_2_6;
}


void calc_MCSH_7_3_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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

    double lambda_x0_5 = lambda_x0_4 * lambda_x0;
    double lambda_y0_5 = lambda_y0_4 * lambda_y0;
    double lambda_z0_5 = lambda_z0_4 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x_2 = lambda_x0_sqr + (1.0/(2.0*gamma));
    double temp_y_2 = lambda_y0_sqr + (1.0/(2.0*gamma));
    double temp_z_2 = lambda_z0_sqr + (1.0/(2.0*gamma));

    double C3_1 = 3.0 / (2.0 * gamma);
    double temp_x_3 = lambda_x0_3 + C3_1 * lambda_x0;
    double temp_y_3 = lambda_y0_3 + C3_1 * lambda_y0;
    double temp_z_3 = lambda_z0_3 + C3_1 * lambda_z0;

    double C5_1 = 5.0 / gamma, C5_2 = 15.0 / (4.0 * gamma * gamma);
    double temp_x_5 = lambda_x0_5 + C5_1 * lambda_x0_3 + C5_2 * lambda_x0;
    double temp_y_5 = lambda_y0_5 + C5_1 * lambda_y0_3 + C5_2 * lambda_y0;
    double temp_z_5 = lambda_z0_5 + C5_1 * lambda_z0_3 + C5_2 * lambda_z0;

    double temp_term1_x = 135135.0 * temp_x_5 - 103950.0 * temp_x_3 + 14175.0 * lambda_x0;
    double temp_term1_y = 135135.0 * temp_y_5 - 103950.0 * temp_y_3 + 14175.0 * lambda_y0;
    double temp_term1_z = 135135.0 * temp_z_5 - 103950.0 * temp_z_3 + 14175.0 * lambda_z0;

    double temp_term2_x = -10395.0 * temp_x_5 + 9450.0 * temp_x_3 - 1575.0 * lambda_x0;
    double temp_term2_y = -10395.0 * temp_y_5 + 9450.0 * temp_y_3 - 1575.0 * lambda_y0;
    double temp_term2_z = -10395.0 * temp_z_5 + 9450.0 * temp_z_3 - 1575.0 * lambda_z0;

    double temp_miu1 = temp_y_2 * temp_term1_x + temp_term2_x;
    double temp_miu2 = temp_x_2 * temp_term1_y + temp_term2_y;
    double temp_miu3 = temp_z_2 * temp_term1_x + temp_term2_x;
    double temp_miu4 = temp_x_2 * temp_term1_z + temp_term2_z;
    double temp_miu5 = temp_z_2 * temp_term1_y + temp_term2_y;
    double temp_miu6 = temp_y_2 * temp_term1_z + temp_term2_z;

    double miu_7_3_1 = temp * temp_miu1;
    double miu_7_3_2 = temp * temp_miu2;
    double miu_7_3_3 = temp * temp_miu3;
    double miu_7_3_4 = temp * temp_miu4;
    double miu_7_3_5 = temp * temp_miu5;
    double miu_7_3_6 = temp * temp_miu6;

    value[0] = miu_7_3_1;
    value[1] = miu_7_3_2;
    value[2] = miu_7_3_3;
    value[3] = miu_7_3_4;
    value[4] = miu_7_3_5;
    value[5] = miu_7_3_6;
}


void calc_MCSH_7_4_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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

    double lambda_x0_5 = lambda_x0_4 * lambda_x0;
    double lambda_y0_5 = lambda_y0_4 * lambda_y0;
    double lambda_z0_5 = lambda_z0_4 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = (135135.0 * 5.0 / gamma) - 103950.0;
    double C4 = (135135.0 * 15.0 / (4.0 * gamma * gamma)) - (103950.0 * 3.0 / (2.0 * gamma)) + 14175.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 135135.0 * lambda_x0_5 + C3 * lambda_x0_3 + C4 * lambda_x0;
    double temp_y = 135135.0 * lambda_y0_5 + C3 * lambda_y0_3 + C4 * lambda_y0;
    double temp_z = 135135.0 * lambda_z0_5 + C3 * lambda_z0_3 + C4 * lambda_z0;

    double miu_7_4_1 = temp * lambda_y0 * lambda_z0 * temp_x;
    double miu_7_4_2 = temp * lambda_x0 * lambda_z0 * temp_y;
    double miu_7_4_3 = temp * lambda_x0 * lambda_y0 * temp_z;

    value[0] = miu_7_4_1;
    value[1] = miu_7_4_2;
    value[2] = miu_7_4_3;
}

void calc_MCSH_7_5_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x_2 = lambda_x0_sqr + (1.0/(2.0*gamma));
    double temp_y_2 = lambda_y0_sqr + (1.0/(2.0*gamma));
    double temp_z_2 = lambda_z0_sqr + (1.0/(2.0*gamma));

    double C3_1 = 3.0 / (2.0 * gamma);
    double temp_x_3 = lambda_x0_3 + C3_1 * lambda_x0;
    double temp_y_3 = lambda_y0_3 + C3_1 * lambda_y0;
    double temp_z_3 = lambda_z0_3 + C3_1 * lambda_z0;

    double C4_1 = 3.0 / gamma, C4_2 = 3.0 / (4.0 * gamma * gamma);
    double temp_x_4 = lambda_x0_4 + C4_1 * lambda_x0_sqr + C4_2;
    double temp_y_4 = lambda_y0_4 + C4_1 * lambda_y0_sqr + C4_2;
    double temp_z_4 = lambda_z0_4 + C4_1 * lambda_z0_sqr + C4_2;

    double temp_term1_x = 135135.0 * temp_x_4 - 62370.0 * temp_x_2 + 2835.0;
    double temp_term1_y = 135135.0 * temp_y_4 - 62370.0 * temp_y_2 + 2835.0;
    double temp_term1_z = 135135.0 * temp_z_4 - 62370.0 * temp_z_2 + 2835.0;

    double temp_term2_x = -31185.0 * temp_x_4 + 17010.0 * temp_x_2 - 945.0;
    double temp_term2_y = -31185.0 * temp_y_4 + 17010.0 * temp_y_2 - 945.0;
    double temp_term2_z = -31185.0 * temp_z_4 + 17010.0 * temp_z_2 - 945.0;

    double temp_miu1 = temp_y_3 * temp_term1_x + lambda_y0 * temp_term2_x;
    double temp_miu2 = temp_x_3 * temp_term1_y + lambda_x0 * temp_term2_y;
    double temp_miu3 = temp_z_3 * temp_term1_x + lambda_z0 * temp_term2_x;
    double temp_miu4 = temp_x_3 * temp_term1_z + lambda_x0 * temp_term2_z;
    double temp_miu5 = temp_z_3 * temp_term1_y + lambda_z0 * temp_term2_y;
    double temp_miu6 = temp_y_3 * temp_term1_z + lambda_y0 * temp_term2_z;

    double miu_7_5_1 = temp * temp_miu1;
    double miu_7_5_2 = temp * temp_miu2;
    double miu_7_5_3 = temp * temp_miu3;
    double miu_7_5_4 = temp * temp_miu4;
    double miu_7_5_5 = temp * temp_miu5;
    double miu_7_5_6 = temp * temp_miu6;

    value[0] = miu_7_5_1;
    value[1] = miu_7_5_2;
    value[2] = miu_7_5_3;
    value[3] = miu_7_5_4;
    value[4] = miu_7_5_5;
    value[5] = miu_7_5_6;
}

void calc_MCSH_7_6_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x_2 = lambda_x0_sqr + (1.0/(2.0*gamma));
    double temp_y_2 = lambda_y0_sqr + (1.0/(2.0*gamma));
    double temp_z_2 = lambda_z0_sqr + (1.0/(2.0*gamma));

    double C4_1 = 3.0 / gamma, C4_2 = 3.0 / (4.0 * gamma * gamma);
    double temp_x_4 = lambda_x0_4 + C4_1 * lambda_x0_sqr + C4_2;
    double temp_y_4 = lambda_y0_4 + C4_1 * lambda_y0_sqr + C4_2;
    double temp_z_4 = lambda_z0_4 + C4_1 * lambda_z0_sqr + C4_2;

    double temp_term1_x = 135135.0 * temp_x_4 - 62370.0 * temp_x_2 + 2835.0;
    double temp_term1_y = 135135.0 * temp_y_4 - 62370.0 * temp_y_2 + 2835.0;
    double temp_term1_z = 135135.0 * temp_z_4 - 62370.0 * temp_z_2 + 2835.0;

    double temp_term2_x = -10395.0 * temp_x_4 + 5670.0 * temp_x_2 - 315.0;
    double temp_term2_y = -10395.0 * temp_y_4 + 5670.0 * temp_y_2 - 315.0;
    double temp_term2_z = -10395.0 * temp_z_4 + 5670.0 * temp_z_2 - 315.0;

    double temp_miu1 = temp_y_2 * temp_term1_x + temp_term2_x;
    double temp_miu2 = temp_x_2 * temp_term1_y + temp_term2_y;
    double temp_miu3 = temp_z_2 * temp_term1_x + temp_term2_x;
    double temp_miu4 = temp_x_2 * temp_term1_z + temp_term2_z;
    double temp_miu5 = temp_z_2 * temp_term1_y + temp_term2_y;
    double temp_miu6 = temp_y_2 * temp_term1_z + temp_term2_z;

    double miu_7_6_1 = temp * lambda_z0 * temp_miu1;
    double miu_7_6_2 = temp * lambda_z0 * temp_miu2;
    double miu_7_6_3 = temp * lambda_y0 * temp_miu3;
    double miu_7_6_4 = temp * lambda_y0 * temp_miu4;
    double miu_7_6_5 = temp * lambda_x0 * temp_miu5;
    double miu_7_6_6 = temp * lambda_x0 * temp_miu6;

    value[0] = miu_7_6_1;
    value[1] = miu_7_6_2;
    value[2] = miu_7_6_3;
    value[3] = miu_7_6_4;
    value[4] = miu_7_6_5;
    value[5] = miu_7_6_6;
}


void calc_MCSH_7_7_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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

    double temp = C1 * exp( C2 * r0_sqr);

    double C3_1 = 3.0 / (2.0 * gamma);
    double temp_x_3 = lambda_x0_3 + C3_1 * lambda_x0;
    double temp_y_3 = lambda_y0_3 + C3_1 * lambda_y0;
    double temp_z_3 = lambda_z0_3 + C3_1 * lambda_z0;

    double temp_term1_x = 135135.0 * temp_x_3 - 31185.0 * lambda_x0;
    double temp_term1_y = 135135.0 * temp_y_3 - 31185.0 * lambda_y0;

    double temp_term2_x = -31185.0 * temp_x_3 + 8505.0 * lambda_x0;
    double temp_term2_y = -31185.0 * temp_y_3 + 8505.0 * lambda_y0;

    double temp_miu1 = temp_y_3 * temp_term1_x + lambda_y0 * temp_term2_x;
    double temp_miu2 = temp_z_3 * temp_term1_x + lambda_z0 * temp_term2_x;
    double temp_miu3 = temp_z_3 * temp_term1_y + lambda_z0 * temp_term2_y;

    double miu_7_7_1 = temp * lambda_z0 * temp_miu1;
    double miu_7_7_2 = temp * lambda_y0 * temp_miu2;
    double miu_7_7_3 = temp * lambda_x0 * temp_miu3;

    value[0] = miu_7_7_1;
    value[1] = miu_7_7_2;
    value[2] = miu_7_7_3;
}

void calc_MCSH_7_8_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x_2 = lambda_x0_sqr + (1.0/(2.0*gamma));
    double temp_y_2 = lambda_y0_sqr + (1.0/(2.0*gamma));
    double temp_z_2 = lambda_z0_sqr + (1.0/(2.0*gamma));

    double C3_1 = 3.0 / (2.0 * gamma);
    double temp_x_3 = lambda_x0_3 + C3_1 * lambda_x0;
    double temp_y_3 = lambda_y0_3 + C3_1 * lambda_y0;
    double temp_z_3 = lambda_z0_3 + C3_1 * lambda_z0;

    double temp_term1_x = 135135.0 * temp_x_3 - 31185.0 * lambda_x0;
    double temp_term1_y = 135135.0 * temp_y_3 - 31185.0 * lambda_y0;
    double temp_term1_z = 135135.0 * temp_z_3 - 31185.0 * lambda_z0;

    double temp_term2_x = -10395.0 * temp_x_3 + 2835.0 * lambda_x0;
    double temp_term2_y = -10395.0 * temp_y_3 + 2835.0 * lambda_y0;
    double temp_term2_z = -10395.0 * temp_z_3 + 2835.0 * lambda_z0;

    double temp_term3_x = 945.0 * temp_x_3 - 315.0 * lambda_x0;
    double temp_term3_y = 945.0 * temp_y_3 - 315.0 * lambda_y0;
    double temp_term3_z = 945.0 * temp_z_3 - 315.0 * lambda_z0;

    double temp_miu1 = temp_y_2 * temp_z_2 * temp_term1_x + (temp_y_2 + temp_z_2) * temp_term2_x + temp_term3_x;
    double temp_miu2 = temp_x_2 * temp_z_2 * temp_term1_y + (temp_x_2 + temp_z_2) * temp_term2_y + temp_term3_y;
    double temp_miu3 = temp_x_2 * temp_y_2 * temp_term1_z + (temp_x_2 + temp_y_2) * temp_term2_z + temp_term3_z;

    double miu_7_8_1 = temp * temp_miu1;
    double miu_7_8_2 = temp * temp_miu2;
    double miu_7_8_3 = temp * temp_miu3;

    value[0] = miu_7_8_1;
    value[1] = miu_7_8_2;
    value[2] = miu_7_8_3;
}

void calc_MCSH_8_1_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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

    double lambda_x0_6 = lambda_x0_5 * lambda_x0;
    double lambda_y0_6 = lambda_y0_5 * lambda_y0;
    double lambda_z0_6 = lambda_z0_5 * lambda_z0;

    double lambda_x0_7 = lambda_x0_6 * lambda_x0;
    double lambda_y0_7 = lambda_y0_6 * lambda_y0;
    double lambda_z0_7 = lambda_z0_6 * lambda_z0;

    double lambda_x0_8 = lambda_x0_7 * lambda_x0;
    double lambda_y0_8 = lambda_y0_7 * lambda_y0;
    double lambda_z0_8 = lambda_z0_7 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = (2027025.0 * 14.0 / gamma) - 3783780.0;
    double C4 = (2027025.0 * 105.0 / (2.0 * gamma * gamma)) - (3783780.0 * 15.0 / (2.0 * gamma)) + 2182950.0;
    double C5 = (2027025.0 * 105.0 / (2.0 * gamma * gamma * gamma)) - (3783780.0 * 45.0 / (4.0 * gamma * gamma)) + (2182950.0 * 3.0 / gamma) - 396900.0;
    double C6 = (2027025.0 * 105.0 / (16.0 * gamma * gamma * gamma * gamma)) - (3783780.0 * 15.0 / (8.0 * gamma * gamma * gamma)) + (2182950.0 * 3.0 / (4.0 * gamma *gamma)) - (396900.0 / (2.0 * gamma)) + 11025.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 2027025.0 * lambda_x0_8 + C3 * lambda_x0_6 + C4 * lambda_x0_4 + C5 * lambda_x0_sqr + C6;
    double temp_y = 2027025.0 * lambda_y0_8 + C3 * lambda_y0_6 + C4 * lambda_y0_4 + C5 * lambda_y0_sqr + C6;
    double temp_z = 2027025.0 * lambda_z0_8 + C3 * lambda_z0_6 + C4 * lambda_z0_4 + C5 * lambda_z0_sqr + C6;

    double miu_8_1_1 = temp * temp_x;
    double miu_8_1_2 = temp * temp_y;
    double miu_8_1_3 = temp * temp_z;

    value[0] = miu_8_1_1;
    value[1] = miu_8_1_2;
    value[2] = miu_8_1_3;
}

void calc_MCSH_8_2_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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

    double lambda_x0_5 = lambda_x0_4 * lambda_x0;
    double lambda_y0_5 = lambda_y0_4 * lambda_y0;
    double lambda_z0_5 = lambda_z0_4 * lambda_z0;

    double lambda_x0_6 = lambda_x0_5 * lambda_x0;
    double lambda_y0_6 = lambda_y0_5 * lambda_y0;
    double lambda_z0_6 = lambda_z0_5 * lambda_z0;

    double lambda_x0_7 = lambda_x0_6 * lambda_x0;
    double lambda_y0_7 = lambda_y0_6 * lambda_y0;
    double lambda_z0_7 = lambda_z0_6 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = (2027025.0 * 21.0 / (2.0 * gamma)) - 2837835.0;
    double C4 = (2027025.0 * 105.0 / (4.0 * gamma * gamma)) - (2837835.0 * 5.0 / gamma) + 1091475.0;
    double C5 = (2027025.0 * 105.0 / (8.0 * gamma * gamma * gamma)) - (2837835.0 * 15.0 / (4.0 * gamma * gamma)) + (1091475.0 * 3.0 / (2.0 * gamma)) - 99225.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 2027025.0 * lambda_x0_7 + C3 * lambda_x0_5 + C4 * lambda_x0_3 + C5 * lambda_x0;
    double temp_y = 2027025.0 * lambda_y0_7 + C3 * lambda_y0_5 + C4 * lambda_y0_3 + C5 * lambda_y0;
    double temp_z = 2027025.0 * lambda_z0_7 + C3 * lambda_z0_5 + C4 * lambda_z0_3 + C5 * lambda_z0;

    double miu_8_2_1 = temp * lambda_y0 * temp_x;
    double miu_8_2_2 = temp * lambda_x0 * temp_y;
    double miu_8_2_3 = temp * lambda_z0 * temp_x;
    double miu_8_2_4 = temp * lambda_x0 * temp_z;
    double miu_8_2_5 = temp * lambda_z0 * temp_y;
    double miu_8_2_6 = temp * lambda_y0 * temp_z;

    value[0] = miu_8_2_1;
    value[1] = miu_8_2_2;
    value[2] = miu_8_2_3;
    value[3] = miu_8_2_4;
    value[4] = miu_8_2_5;
    value[5] = miu_8_2_6;
}


void calc_MCSH_8_3_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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

    double lambda_x0_5 = lambda_x0_4 * lambda_x0;
    double lambda_y0_5 = lambda_y0_4 * lambda_y0;
    double lambda_z0_5 = lambda_z0_4 * lambda_z0;

    double lambda_x0_6 = lambda_x0_5 * lambda_x0;
    double lambda_y0_6 = lambda_y0_5 * lambda_y0;
    double lambda_z0_6 = lambda_z0_5 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x_2 = lambda_x0_sqr + (1.0/(2.0*gamma));
    double temp_y_2 = lambda_y0_sqr + (1.0/(2.0*gamma));
    double temp_z_2 = lambda_z0_sqr + (1.0/(2.0*gamma));

    double C4_1 = 3.0 / gamma, C4_2 = 3.0 / (4.0 * gamma * gamma);
    double temp_x_4 = lambda_x0_4 + C4_1 * lambda_x0_sqr + C4_2;
    double temp_y_4 = lambda_y0_4 + C4_1 * lambda_y0_sqr + C4_2;
    double temp_z_4 = lambda_z0_4 + C4_1 * lambda_z0_sqr + C4_2;

    double C6_1 = 15.0 / (2.0 * gamma), C6_2 = 45.0 / (4.0 * gamma * gamma), C6_3 =  15.0 / (8.0 * gamma * gamma * gamma);
    double temp_x_6 = lambda_x0_6 + C6_1 * lambda_x0_4 + C6_2 * lambda_x0_sqr + C6_3;
    double temp_y_6 = lambda_y0_6 + C6_1 * lambda_y0_4 + C6_2 * lambda_y0_sqr + C6_3;
    double temp_z_6 = lambda_z0_6 + C6_1 * lambda_z0_4 + C6_2 * lambda_z0_sqr + C6_3;

    double temp_term1_x = 2027025.0 * temp_x_6 - 2027025.0 * temp_x_4 + 467775.0 * temp_x_2 - 14175.0;
    double temp_term1_y = 2027025.0 * temp_y_6 - 2027025.0 * temp_y_4 + 467775.0 * temp_y_2 - 14175.0;
    double temp_term1_z = 2027025.0 * temp_z_6 - 2027025.0 * temp_z_4 + 467775.0 * temp_z_2 - 14175.0;

    double temp_term2_x = -135135.0 * temp_x_6 + 155925.0 * temp_x_4 - 42525.0 * temp_x_2 + 1575.0;
    double temp_term2_y = -135135.0 * temp_y_6 + 155925.0 * temp_y_4 - 42525.0 * temp_y_2 + 1575.0;
    double temp_term2_z = -135135.0 * temp_z_6 + 155925.0 * temp_z_4 - 42525.0 * temp_z_2 + 1575.0;

    double temp_miu1 = temp_y_2 * temp_term1_x + temp_term2_x;
    double temp_miu2 = temp_x_2 * temp_term1_y + temp_term2_y;
    double temp_miu3 = temp_z_2 * temp_term1_x + temp_term2_x;
    double temp_miu4 = temp_x_2 * temp_term1_z + temp_term2_z;
    double temp_miu5 = temp_z_2 * temp_term1_y + temp_term2_y;
    double temp_miu6 = temp_y_2 * temp_term1_z + temp_term2_z;

    double miu_8_3_1 = temp * temp_miu1;
    double miu_8_3_2 = temp * temp_miu2;
    double miu_8_3_3 = temp * temp_miu3;
    double miu_8_3_4 = temp * temp_miu4;
    double miu_8_3_5 = temp * temp_miu5;
    double miu_8_3_6 = temp * temp_miu6;

    value[0] = miu_8_3_1;
    value[1] = miu_8_3_2;
    value[2] = miu_8_3_3;
    value[3] = miu_8_3_4;
    value[4] = miu_8_3_5;
    value[5] = miu_8_3_6;
}


void calc_MCSH_8_4_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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

    double lambda_x0_5 = lambda_x0_4 * lambda_x0;
    double lambda_y0_5 = lambda_y0_4 * lambda_y0;
    double lambda_z0_5 = lambda_z0_4 * lambda_z0;

    double lambda_x0_6 = lambda_x0_5 * lambda_x0;
    double lambda_y0_6 = lambda_y0_5 * lambda_y0;
    double lambda_z0_6 = lambda_z0_5 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = (2027025.0 * 15.0 / (2.0 * gamma)) - 2027025.0;
    double C4 = (2027025.0 * 45.0 / (4.0 * gamma * gamma)) - (2027025.0 * 3.0 / gamma) + 467775.0;
    double C5 = (2027025.0 * 15.0 / (8.0 * gamma * gamma * gamma)) - (2027025.0 * 3.0 / (4.0 * gamma * gamma)) + (467775.0 / (2.0 * gamma)) - 14175.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 2027025.0 * lambda_x0_6 + C3 * lambda_x0_4 + C4 * lambda_x0_sqr + C5;
    double temp_y = 2027025.0 * lambda_y0_6 + C3 * lambda_y0_4 + C4 * lambda_y0_sqr + C5;
    double temp_z = 2027025.0 * lambda_z0_6 + C3 * lambda_z0_4 + C4 * lambda_z0_sqr + C5;

    double miu_8_4_1 = temp * lambda_y0 * lambda_z0 * temp_x;
    double miu_8_4_2 = temp * lambda_x0 * lambda_z0 * temp_y;
    double miu_8_4_3 = temp * lambda_x0 * lambda_y0 * temp_z;

    value[0] = miu_8_4_1;
    value[1] = miu_8_4_2;
    value[2] = miu_8_4_3;
}

void calc_MCSH_8_5_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);
    // double lambda_sqr = lambda * lambda;

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

    double lambda_x0_5 = lambda_x0_4 * lambda_x0;
    double lambda_y0_5 = lambda_y0_4 * lambda_y0;
    double lambda_z0_5 = lambda_z0_4 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);

    double temp = C1 * exp( C2 * r0_sqr);

    double C3_1 = 3.0 / (2.0 * gamma);
    double temp_x_3 = lambda_x0_3 + C3_1 * lambda_x0;
    double temp_y_3 = lambda_y0_3 + C3_1 * lambda_y0;
    double temp_z_3 = lambda_z0_3 + C3_1 * lambda_z0;

    double C5_1 = 5.0 / gamma, C5_2 = 15.0 / (4.0 * gamma * gamma);
    double temp_x_5 = lambda_x0_5 + C5_1 * lambda_x0_3 + C5_2 * lambda_x0;
    double temp_y_5 = lambda_y0_5 + C5_1 * lambda_y0_3 + C5_2 * lambda_y0;
    double temp_z_5 = lambda_z0_5 + C5_1 * lambda_z0_3 + C5_2 * lambda_z0;

    double temp_term1_x = 2027025.0 * temp_x_5 - 1351350.0 * temp_x_3 + 155925.0 * lambda_x0;
    double temp_term1_y = 2027025.0 * temp_y_5 - 1351350.0 * temp_y_3 + 155925.0 * lambda_y0;
    double temp_term1_z = 2027025.0 * temp_z_5 - 1351350.0 * temp_z_3 + 155925.0 * lambda_z0;

    double temp_term2_x = -405405.0 * temp_x_5 + 311850.0 * temp_x_3 - 42525.0 * lambda_x0;
    double temp_term2_y = -405405.0 * temp_y_5 + 311850.0 * temp_y_3 - 42525.0 * lambda_y0;
    double temp_term2_z = -405405.0 * temp_z_5 + 311850.0 * temp_z_3 - 42525.0 * lambda_z0;

    double temp_miu1 = temp_y_3 * temp_term1_x + lambda_y0 * temp_term2_x;
    double temp_miu2 = temp_x_3 * temp_term1_y + lambda_x0 * temp_term2_y;
    double temp_miu3 = temp_z_3 * temp_term1_x + lambda_z0 * temp_term2_x;
    double temp_miu4 = temp_x_3 * temp_term1_z + lambda_x0 * temp_term2_z;
    double temp_miu5 = temp_z_3 * temp_term1_y + lambda_z0 * temp_term2_y;
    double temp_miu6 = temp_y_3 * temp_term1_z + lambda_y0 * temp_term2_z;

    double miu_8_5_1 = temp * temp_miu1;
    double miu_8_5_2 = temp * temp_miu2;
    double miu_8_5_3 = temp * temp_miu3;
    double miu_8_5_4 = temp * temp_miu4;
    double miu_8_5_5 = temp * temp_miu5;
    double miu_8_5_6 = temp * temp_miu6;

    value[0] = miu_8_5_1;
    value[1] = miu_8_5_2;
    value[2] = miu_8_5_3;
    value[3] = miu_8_5_4;
    value[4] = miu_8_5_5;
    value[5] = miu_8_5_6;
}

void calc_MCSH_8_6_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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

    double lambda_x0_5 = lambda_x0_4 * lambda_x0;
    double lambda_y0_5 = lambda_y0_4 * lambda_y0;
    double lambda_z0_5 = lambda_z0_4 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x_2 = lambda_x0_sqr + (1.0/(2.0*gamma));
    double temp_y_2 = lambda_y0_sqr + (1.0/(2.0*gamma));
    double temp_z_2 = lambda_z0_sqr + (1.0/(2.0*gamma));

    double C3_1 = 3.0 / (2.0 * gamma);
    double temp_x_3 = lambda_x0_3 + C3_1 * lambda_x0;
    double temp_y_3 = lambda_y0_3 + C3_1 * lambda_y0;
    double temp_z_3 = lambda_z0_3 + C3_1 * lambda_z0;

    double C5_1 = 5.0 / gamma, C5_2 = 15.0 / (4.0 * gamma * gamma);
    double temp_x_5 = lambda_x0_5 + C5_1 * lambda_x0_3 + C5_2 * lambda_x0;
    double temp_y_5 = lambda_y0_5 + C5_1 * lambda_y0_3 + C5_2 * lambda_y0;
    double temp_z_5 = lambda_z0_5 + C5_1 * lambda_z0_3 + C5_2 * lambda_z0;

    double temp_term1_x = 2027025.0 * temp_x_5 - 1351350.0 * temp_x_3 + 155925.0 * lambda_x0;
    double temp_term1_y = 2027025.0 * temp_y_5 - 1351350.0 * temp_y_3 + 155925.0 * lambda_y0;
    double temp_term1_z = 2027025.0 * temp_z_5 - 1351350.0 * temp_z_3 + 155925.0 * lambda_z0;

    double temp_term2_x = -135135.0 * temp_x_5 + 103950.0 * temp_x_3 - 14175.0 * lambda_x0;
    double temp_term2_y = -135135.0 * temp_y_5 + 103950.0 * temp_y_3 - 14175.0 * lambda_y0;
    double temp_term2_z = -135135.0 * temp_z_5 + 103950.0 * temp_z_3 - 14175.0 * lambda_z0;

    double temp_miu1 = temp_y_2 * temp_term1_x + temp_term2_x;
    double temp_miu2 = temp_x_2 * temp_term1_y + temp_term2_y;
    double temp_miu3 = temp_z_2 * temp_term1_x + temp_term2_x;
    double temp_miu4 = temp_x_2 * temp_term1_z + temp_term2_z;
    double temp_miu5 = temp_z_2 * temp_term1_y + temp_term2_y;
    double temp_miu6 = temp_y_2 * temp_term1_z + temp_term2_z;

    double miu_8_6_1 = temp * lambda_z0 * temp_miu1;
    double miu_8_6_2 = temp * lambda_z0 * temp_miu2;
    double miu_8_6_3 = temp * lambda_y0 * temp_miu3;
    double miu_8_6_4 = temp * lambda_y0 * temp_miu4;
    double miu_8_6_5 = temp * lambda_x0 * temp_miu5;
    double miu_8_6_6 = temp * lambda_x0 * temp_miu6;

    value[0] = miu_8_6_1;
    value[1] = miu_8_6_2;
    value[2] = miu_8_6_3;
    value[3] = miu_8_6_4;
    value[4] = miu_8_6_5;
    value[5] = miu_8_6_6;
}


void calc_MCSH_8_7_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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
    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x_2 = lambda_x0_sqr + (1.0/(2.0*gamma));
    double temp_y_2 = lambda_y0_sqr + (1.0/(2.0*gamma));
    double temp_z_2 = lambda_z0_sqr + (1.0/(2.0*gamma));

    double C4_1 = 3.0 / gamma, C4_2 = 3.0 / (4.0 * gamma * gamma);
    double temp_x_4 = lambda_x0_4 + C4_1 * lambda_x0_sqr + C4_2;
    double temp_y_4 = lambda_y0_4 + C4_1 * lambda_y0_sqr + C4_2;
    double temp_z_4 = lambda_z0_4 + C4_1 * lambda_z0_sqr + C4_2;

    double temp_term1_x = 2027025.0 * temp_x_4 - 810810.0 * temp_x_2 + 31185.0;
    double temp_term1_y = 2027025.0 * temp_y_4 - 810810.0 * temp_y_2 + 31185.0;
    // double temp_term1_z = 2027025.0 * temp_z_4 - 810810.0 * temp_z_2 + 31185.0;

    double temp_term2_x = -810810.0 * temp_x_4 + 374220.0 * temp_x_2 - 17010.0;
    double temp_term2_y = -810810.0 * temp_y_4 + 374220.0 * temp_y_2 - 17010.0;
    // double temp_term2_z = -810810.0 * temp_z_4 + 374220.0 * temp_z_2 - 17010.0;

    double temp_term3_x = 31185.0 * temp_x_4 - 17010.0 * temp_x_2 + 945.0;
    double temp_term3_y = 31185.0 * temp_y_4 - 17010.0 * temp_y_2 + 945.0;
    // double temp_term3_z = 31185.0 * temp_z_4 - 17010.0 * temp_z_2 + 945.0;

    double temp_miu1 = temp_y_4 * temp_term1_x + temp_y_2 * temp_term2_x + temp_term3_x;
    double temp_miu2 = temp_z_4 * temp_term1_x + temp_z_2 * temp_term2_x + temp_term3_x;
    double temp_miu3 = temp_z_4 * temp_term1_y + temp_z_2 * temp_term2_y + temp_term3_y;

    double miu_8_7_1 = temp * temp_miu1;
    double miu_8_7_2 = temp * temp_miu2;
    double miu_8_7_3 = temp * temp_miu3;

    value[0] = miu_8_7_1;
    value[1] = miu_8_7_2;
    value[2] = miu_8_7_3;
}


void calc_MCSH_8_8_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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
    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x_2 = lambda_x0_sqr + (1.0/(2.0*gamma));
    double temp_y_2 = lambda_y0_sqr + (1.0/(2.0*gamma));
    double temp_z_2 = lambda_z0_sqr + (1.0/(2.0*gamma));

    double C3_1 = 3.0 / (2.0 * gamma);
    double temp_x_3 = lambda_x0_3 + C3_1 * lambda_x0;
    double temp_y_3 = lambda_y0_3 + C3_1 * lambda_y0;
    double temp_z_3 = lambda_z0_3 + C3_1 * lambda_z0;

    double C4_1 = 3.0 / gamma, C4_2 = 3.0 / (4.0 * gamma * gamma);
    double temp_x_4 = lambda_x0_4 + C4_1 * lambda_x0_sqr + C4_2;
    double temp_y_4 = lambda_y0_4 + C4_1 * lambda_y0_sqr + C4_2;
    double temp_z_4 = lambda_z0_4 + C4_1 * lambda_z0_sqr + C4_2;

    double temp_term1_x = 2027025.0 * temp_x_4 - 810810.0 * temp_x_2 + 31185.0;
    double temp_term1_y = 2027025.0 * temp_y_4 - 810810.0 * temp_y_2 + 31185.0;
    double temp_term1_z = 2027025.0 * temp_z_4 - 810810.0 * temp_z_2 + 31185.0;

    double temp_term2_x = -405405.0 * temp_x_4 + 187110.0 * temp_x_2 - 8505.0;
    double temp_term2_y = -405405.0 * temp_y_4 + 187110.0 * temp_y_2 - 8505.0;
    double temp_term2_z = -405405.0 * temp_z_4 + 187110.0 * temp_z_2 - 8505.0;

    double temp_miu1 = temp_y_3 * temp_term1_x + lambda_y0 * temp_term2_x;
    double temp_miu2 = temp_x_3 * temp_term1_y + lambda_x0 * temp_term2_y;
    double temp_miu3 = temp_z_3 * temp_term1_x + lambda_z0 * temp_term2_x;
    double temp_miu4 = temp_x_3 * temp_term1_z + lambda_x0 * temp_term2_z;
    double temp_miu5 = temp_z_3 * temp_term1_y + lambda_z0 * temp_term2_y;
    double temp_miu6 = temp_y_3 * temp_term1_z + lambda_y0 * temp_term2_z;

    double miu_8_8_1 = temp * lambda_z0 * temp_miu1;
    double miu_8_8_2 = temp * lambda_z0 * temp_miu2;
    double miu_8_8_3 = temp * lambda_y0 * temp_miu3;
    double miu_8_8_4 = temp * lambda_y0 * temp_miu4;
    double miu_8_8_5 = temp * lambda_x0 * temp_miu5;
    double miu_8_8_6 = temp * lambda_x0 * temp_miu6;

    value[0] = miu_8_8_1;
    value[1] = miu_8_8_2;
    value[2] = miu_8_8_3;
    value[3] = miu_8_8_4;
    value[4] = miu_8_8_5;
    value[5] = miu_8_8_6;
}

void calc_MCSH_8_9_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x_2 = lambda_x0_sqr + (1.0/(2.0*gamma));
    double temp_y_2 = lambda_y0_sqr + (1.0/(2.0*gamma));
    double temp_z_2 = lambda_z0_sqr + (1.0/(2.0*gamma));

    double C4_1 = 3.0 / gamma, C4_2 = 3.0 / (4.0 * gamma * gamma);
    double temp_x_4 = lambda_x0_4 + C4_1 * lambda_x0_sqr + C4_2;
    double temp_y_4 = lambda_y0_4 + C4_1 * lambda_y0_sqr + C4_2;
    double temp_z_4 = lambda_z0_4 + C4_1 * lambda_z0_sqr + C4_2;

    double temp_term1_x = 2027025.0 * temp_x_4 - 810810.0 * temp_x_2 + 31185.0;
    double temp_term1_y = 2027025.0 * temp_y_4 - 810810.0 * temp_y_2 + 31185.0;
    double temp_term1_z = 2027025.0 * temp_z_4 - 810810.0 * temp_z_2 + 31185.0;

    double temp_term2_x = -135135.0 * temp_x_4 + 62370.0 * temp_x_2 - 2835.0;
    double temp_term2_y = -135135.0 * temp_y_4 + 62370.0 * temp_y_2 - 2835.0;
    double temp_term2_z = -135135.0 * temp_z_4 + 62370.0 * temp_z_2 - 2835.0;

    double temp_term3_x = 10395.0 * temp_x_4 - 5670.0 * temp_x_2 + 315.0;
    double temp_term3_y = 10395.0 * temp_y_4 - 5670.0 * temp_y_2 + 315.0;
    double temp_term3_z = 10395.0 * temp_z_4 - 5670.0 * temp_z_2 + 315.0;

    double temp_miu1 = temp_y_2 * temp_z_2 * temp_term1_x + (temp_y_2 + temp_z_2) * temp_term2_x + temp_term3_x;
    double temp_miu2 = temp_x_2 * temp_z_2 * temp_term1_y + (temp_x_2 + temp_z_2) * temp_term2_y + temp_term3_y;
    double temp_miu3 = temp_x_2 * temp_y_2 * temp_term1_z + (temp_x_2 + temp_y_2) * temp_term2_z + temp_term3_z;

    double miu_8_9_1 = temp * temp_miu1;
    double miu_8_9_2 = temp * temp_miu2;
    double miu_8_9_3 = temp * temp_miu3;

    value[0] = miu_8_9_1;
    value[1] = miu_8_9_2;
    value[2] = miu_8_9_3;
}

void calc_MCSH_8_10_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x_2 = lambda_x0_sqr + (1.0/(2.0*gamma));
    double temp_y_2 = lambda_y0_sqr + (1.0/(2.0*gamma));
    double temp_z_2 = lambda_z0_sqr + (1.0/(2.0*gamma));

    double C3_1 = 3.0 / (2.0 * gamma);
    double temp_x_3 = lambda_x0_3 + C3_1 * lambda_x0;
    double temp_y_3 = lambda_y0_3 + C3_1 * lambda_y0;
    double temp_z_3 = lambda_z0_3 + C3_1 * lambda_z0;

    double temp_term1_x = 2027025.0 * temp_x_3 - 405405.0 * lambda_x0;
    double temp_term1_y = 2027025.0 * temp_y_3 - 405405.0 * lambda_y0;
    // double temp_term1_z = 2027025.0 * temp_z_3 - 405405.0 * lambda_z0;

    double temp_term2_x = -405405.0 * temp_x_3 + 93555.0 * lambda_x0;
    double temp_term2_y = -405405.0 * temp_y_3 + 93555.0 * lambda_y0;
    // double temp_term2_z = -405405.0 * temp_z_3 + 93555.0 * lambda_z0;

    double temp_term3_x = -135135.0 * temp_x_3 + 31185.0 * lambda_x0;
    double temp_term3_y = -135135.0 * temp_y_3 + 31185.0 * lambda_y0;
    // double temp_term3_z = -135135.0 * temp_z_3 + 31185.0 * lambda_z0;

    double temp_term4_x = 31185.0 * temp_x_3 - 8505.0 * lambda_x0;
    double temp_term4_y = 31185.0 * temp_y_3 - 8505.0 * lambda_y0;
    // double temp_term4_z = 31185.0 * temp_z_3 - 8505.0 * lambda_z0;

    double temp_miu1 = temp_y_3 * temp_z_2 * temp_term1_x + lambda_y0 * temp_z_2 * temp_term2_x + temp_y_3 * temp_term3_x + lambda_y0 * temp_term4_x;
    double temp_miu2 = temp_z_3 * temp_y_2 * temp_term1_x + lambda_z0 * temp_y_2 * temp_term2_x + temp_z_3 * temp_term3_x + lambda_z0 * temp_term4_x;
    double temp_miu3 = temp_z_3 * temp_x_2 * temp_term1_y + lambda_z0 * temp_x_2 * temp_term2_y + temp_z_3 * temp_term3_y + lambda_z0 * temp_term4_y;

    double miu_8_10_1 = temp * temp_miu1;
    double miu_8_10_2 = temp * temp_miu2;
    double miu_8_10_3 = temp * temp_miu3;

    value[0] = miu_8_10_1;
    value[1] = miu_8_10_2;
    value[2] = miu_8_10_3;
}

void calc_MCSH_9_1_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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

    double lambda_x0_6 = lambda_x0_5 * lambda_x0;
    double lambda_y0_6 = lambda_y0_5 * lambda_y0;
    double lambda_z0_6 = lambda_z0_5 * lambda_z0;

    double lambda_x0_7 = lambda_x0_6 * lambda_x0;
    double lambda_y0_7 = lambda_y0_6 * lambda_y0;
    double lambda_z0_7 = lambda_z0_6 * lambda_z0;

    double lambda_x0_8 = lambda_x0_7 * lambda_x0;
    double lambda_y0_8 = lambda_y0_7 * lambda_y0;
    double lambda_z0_8 = lambda_z0_7 * lambda_z0;

    double lambda_x0_9 = lambda_x0_8 * lambda_x0;
    double lambda_y0_9 = lambda_y0_8 * lambda_y0;
    double lambda_z0_9 = lambda_z0_8 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = (34459425.0 * 18.0 / gamma) - 72972900.0;
    double C4 = (34459425.0 * 189.0 / (2.0 * gamma * gamma)) - (72972900.0 * 21.0 / (2.0 * gamma)) + 51081030.0;
    double C5 = (34459425.0 * 315.0 / (2.0 * gamma * gamma * gamma)) - (72972900.0 * 105.0 / (4.0 * gamma * gamma)) + (51081030.0 * 5.0 / gamma) - 13097700.0;
    double C6 = (34459425.0 * 945.0 / (16.0 * gamma * gamma * gamma * gamma)) - (72972900.0 * 105.0 / (8.0 * gamma * gamma * gamma)) + (51081030.0 * 15.0 / (4.0 * gamma *gamma)) - (13097700.0 * 3.0 / (2.0 * gamma)) + 893025.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 34459425.0 * lambda_x0_9 + C3 * lambda_x0_7 + C4 * lambda_x0_5 + C5 * lambda_x0_3 + C6 * lambda_x0;
    double temp_y = 34459425.0 * lambda_y0_9 + C3 * lambda_y0_7 + C4 * lambda_y0_5 + C5 * lambda_y0_3 + C6 * lambda_y0;
    double temp_z = 34459425.0 * lambda_z0_9 + C3 * lambda_z0_7 + C4 * lambda_z0_5 + C5 * lambda_z0_3 + C6 * lambda_z0;

    double miu_9_1_1 = temp * temp_x;
    double miu_9_1_2 = temp * temp_y;
    double miu_9_1_3 = temp * temp_z;

    value[0] = miu_9_1_1;
    value[1] = miu_9_1_2;
    value[2] = miu_9_1_3;
}

void calc_MCSH_9_2_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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

    double lambda_x0_5 = lambda_x0_4 * lambda_x0;
    double lambda_y0_5 = lambda_y0_4 * lambda_y0;
    double lambda_z0_5 = lambda_z0_4 * lambda_z0;

    double lambda_x0_6 = lambda_x0_5 * lambda_x0;
    double lambda_y0_6 = lambda_y0_5 * lambda_y0;
    double lambda_z0_6 = lambda_z0_5 * lambda_z0;

    double lambda_x0_7 = lambda_x0_6 * lambda_x0;
    double lambda_y0_7 = lambda_y0_6 * lambda_y0;
    double lambda_z0_7 = lambda_z0_6 * lambda_z0;

    double lambda_x0_8 = lambda_x0_7 * lambda_x0;
    double lambda_y0_8 = lambda_y0_7 * lambda_y0;
    double lambda_z0_8 = lambda_z0_7 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = (34459425.0 * 14.0 / gamma) - 56756700.0;
    double C4 = (34459425.0 * 105.0 / (2.0 * gamma * gamma)) - (56756700.0 * 15.0 / (2.0 * gamma)) + 28378350.0;
    double C5 = (34459425.0 * 105.0 / (2.0 * gamma * gamma * gamma)) - (56756700.0 * 45.20 / (4.0 * gamma * gamma)) + (28378350.0 * 3.0 / gamma) - 4365900.0;
    double C6 = (34459425.0 * 105.0 / (16.0 * gamma * gamma * gamma * gamma)) - (56756700.0 * 15.0 / (8.0 * gamma * gamma * gamma)) + (28378350.0 * 3.0 / (4.0 * gamma * gamma)) - (4365900.0 / (2.0 * gamma)) + 99225.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 34459425.0 * lambda_x0_8 + C3 * lambda_x0_6 + C4 * lambda_x0_4 + C5 * lambda_x0_sqr + C6;
    double temp_y = 34459425.0 * lambda_y0_8 + C3 * lambda_y0_6 + C4 * lambda_y0_4 + C5 * lambda_y0_sqr + C6;
    double temp_z = 34459425.0 * lambda_z0_8 + C3 * lambda_z0_6 + C4 * lambda_z0_4 + C5 * lambda_z0_sqr + C6;

    double miu_9_2_1 = temp * lambda_y0 * temp_x;
    double miu_9_2_2 = temp * lambda_x0 * temp_y;
    double miu_9_2_3 = temp * lambda_z0 * temp_x;
    double miu_9_2_4 = temp * lambda_x0 * temp_z;
    double miu_9_2_5 = temp * lambda_z0 * temp_y;
    double miu_9_2_6 = temp * lambda_y0 * temp_z;

    value[0] = miu_9_2_1;
    value[1] = miu_9_2_2;
    value[2] = miu_9_2_3;
    value[3] = miu_9_2_4;
    value[4] = miu_9_2_5;
    value[5] = miu_9_2_6;
}

void calc_MCSH_9_3_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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

    double lambda_x0_5 = lambda_x0_4 * lambda_x0;
    double lambda_y0_5 = lambda_y0_4 * lambda_y0;
    double lambda_z0_5 = lambda_z0_4 * lambda_z0;

    double lambda_x0_6 = lambda_x0_5 * lambda_x0;
    double lambda_y0_6 = lambda_y0_5 * lambda_y0;
    double lambda_z0_6 = lambda_z0_5 * lambda_z0;

    double lambda_x0_7 = lambda_x0_6 * lambda_x0;
    double lambda_y0_7 = lambda_y0_6 * lambda_y0;
    double lambda_z0_7 = lambda_z0_6 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x_2 = lambda_x0_sqr + (1/(2*gamma));
    double temp_y_2 = lambda_y0_sqr + (1/(2*gamma));
    double temp_z_2 = lambda_z0_sqr + (1/(2*gamma));

    double C3_1 = 3.0 / (2.0 * gamma);
    double temp_x_3 = lambda_x0_3 + C3_1 * lambda_x0;
    double temp_y_3 = lambda_y0_3 + C3_1 * lambda_y0;
    double temp_z_3 = lambda_z0_3 + C3_1 * lambda_z0;

    double C5_1 = 5.0 / gamma, C5_2 = 15.0 / (4.0 * gamma * gamma);
    double temp_x_5 = lambda_x0_5 + C5_1 * lambda_x0_3 + C5_2 * lambda_x0;
    double temp_y_5 = lambda_y0_5 + C5_1 * lambda_y0_3 + C5_2 * lambda_y0;
    double temp_z_5 = lambda_z0_5 + C5_1 * lambda_z0_3 + C5_2 * lambda_z0;

    double C7_1 = 21 / (2 * gamma), C7_2 = 105 / (4 * gamma * gamma), C7_3 = 105 / (8 * gamma * gamma * gamma);
    double temp_x_7 = lambda_x0_7 + C7_1 * lambda_x0_5 + C7_2 * lambda_x0_3 + C7_3 * lambda_x0;
    double temp_y_7 = lambda_y0_7 + C7_1 * lambda_y0_5 + C7_2 * lambda_y0_3 + C7_3 * lambda_y0;
    double temp_z_7 = lambda_z0_7 + C7_1 * lambda_z0_5 + C7_2 * lambda_z0_3 + C7_3 * lambda_z0;

    double temp_term1_x = 34459425.0 * temp_x_7 - 42567525.0 * temp_x_5 + 14189175.0 * temp_x_3 - 1091475.0 * lambda_x0;
    double temp_term1_y = 34459425.0 * temp_y_7 - 42567525.0 * temp_y_5 + 14189175.0 * temp_y_3 - 1091475.0 * lambda_y0;
    double temp_term1_z = 34459425.0 * temp_z_7 - 42567525.0 * temp_z_5 + 14189175.0 * temp_z_3 - 1091475.0 * lambda_z0;

    double temp_term2_x = -2027025.0 * temp_x_7 + 2837835.0 * temp_x_5 - 1091475.0 * temp_x_3 + 99225.0 * lambda_x0;
    double temp_term2_y = -2027025.0 * temp_y_7 + 2837835.0 * temp_y_5 - 1091475.0 * temp_y_3 + 99225.0 * lambda_y0;
    double temp_term2_z = -2027025.0 * temp_z_7 + 2837835.0 * temp_z_5 - 1091475.0 * temp_z_3 + 99225.0 * lambda_z0;

    double temp_miu1 = temp_y_2 * temp_term1_x + temp_term2_x;
    double temp_miu2 = temp_x_2 * temp_term1_y + temp_term2_y;
    double temp_miu3 = temp_z_2 * temp_term1_x + temp_term2_x;
    double temp_miu4 = temp_x_2 * temp_term1_z + temp_term2_z;
    double temp_miu5 = temp_z_2 * temp_term1_y + temp_term2_y;
    double temp_miu6 = temp_y_2 * temp_term1_z + temp_term2_z;

    double miu_9_3_1 = temp * temp_miu1;
    double miu_9_3_2 = temp * temp_miu2;
    double miu_9_3_3 = temp * temp_miu3;
    double miu_9_3_4 = temp * temp_miu4;
    double miu_9_3_5 = temp * temp_miu5;
    double miu_9_3_6 = temp * temp_miu6;

    value[0] = miu_9_3_1;
    value[1] = miu_9_3_2;
    value[2] = miu_9_3_3;
    value[3] = miu_9_3_4;
    value[4] = miu_9_3_5;
    value[5] = miu_9_3_6;
}

void calc_MCSH_9_4_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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

    double lambda_x0_5 = lambda_x0_4 * lambda_x0;
    double lambda_y0_5 = lambda_y0_4 * lambda_y0;
    double lambda_z0_5 = lambda_z0_4 * lambda_z0;

    double lambda_x0_6 = lambda_x0_5 * lambda_x0;
    double lambda_y0_6 = lambda_y0_5 * lambda_y0;
    double lambda_z0_6 = lambda_z0_5 * lambda_z0;

    double lambda_x0_7 = lambda_x0_6 * lambda_x0;
    double lambda_y0_7 = lambda_y0_6 * lambda_y0;
    double lambda_z0_7 = lambda_z0_6 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);
    double C3 = (34459425.0 * 21.0 / (2.0 * gamma)) - 42567525.0;
    double C4 = (34459425.0 * 105.0 / (4.0 * gamma * gamma)) - (42567525.0 * 5.0 / gamma) + 14189175.0;
    double C5 = (34459425.0 * 105.0 / (8.0 * gamma * gamma * gamma)) - (42567525.0 * 15.0 / (4.0 * gamma * gamma)) + (14189175.0 * 3.0 / (2.0 * gamma)) - 1091475.0;

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x = 34459425.0 * lambda_x0_7 + C3 * lambda_x0_5 + C4 * lambda_x0_3 + C5 * lambda_x0;
    double temp_y = 34459425.0 * lambda_y0_7 + C3 * lambda_y0_5 + C4 * lambda_y0_3 + C5 * lambda_y0;
    double temp_z = 34459425.0 * lambda_z0_7 + C3 * lambda_z0_5 + C4 * lambda_z0_3 + C5 * lambda_z0;

    double miu_9_4_1 = temp * lambda_y0 * lambda_z0 * temp_x;
    double miu_9_4_2 = temp * lambda_x0 * lambda_z0 * temp_y;
    double miu_9_4_3 = temp * lambda_x0 * lambda_y0 * temp_z;

    value[0] = miu_9_4_1;
    value[1] = miu_9_4_2;
    value[2] = miu_9_4_3;
}


void calc_MCSH_9_5_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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

    double lambda_x0_5 = lambda_x0_4 * lambda_x0;
    double lambda_y0_5 = lambda_y0_4 * lambda_y0;
    double lambda_z0_5 = lambda_z0_4 * lambda_z0;

    double lambda_x0_6 = lambda_x0_5 * lambda_x0;
    double lambda_y0_6 = lambda_y0_5 * lambda_y0;
    double lambda_z0_6 = lambda_z0_5 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x_2 = lambda_x0_sqr + (1.0/(2.0*gamma));
    double temp_y_2 = lambda_y0_sqr + (1.0/(2.0*gamma));
    double temp_z_2 = lambda_z0_sqr + (1.0/(2.0*gamma));

    double C3_1 = 3.0 / (2.0 * gamma);
    double temp_x_3 = lambda_x0_3 + C3_1 * lambda_x0;
    double temp_y_3 = lambda_y0_3 + C3_1 * lambda_y0;
    double temp_z_3 = lambda_z0_3 + C3_1 * lambda_z0;

    double C4_1 = 3.0 / gamma, C4_2 = 3.0 / (4.0 * gamma * gamma);
    double temp_x_4 = lambda_x0_4 + C4_1 * lambda_x0_sqr + C4_2;
    double temp_y_4 = lambda_y0_4 + C4_1 * lambda_y0_sqr + C4_2;
    double temp_z_4 = lambda_z0_4 + C4_1 * lambda_z0_sqr + C4_2;

    double C6_1 = 15.0 / (2.0 * gamma), C6_2 = 45.0 / (4.0 * gamma * gamma), C6_3 =  15.0 / (8.0 * gamma * gamma * gamma);
    double temp_x_6 = lambda_x0_6 + C6_1 * lambda_x0_4 + C6_2 * lambda_x0_sqr + C6_3;
    double temp_y_6 = lambda_y0_6 + C6_1 * lambda_y0_4 + C6_2 * lambda_y0_sqr + C6_3;
    double temp_z_6 = lambda_z0_6 + C6_1 * lambda_z0_4 + C6_2 * lambda_z0_sqr + C6_3;

    double temp_term1_x = 34459425.0 * temp_x_6 - 30405375.0 * temp_x_4 + 6081075.0 * temp_x_2 - 155925.0;
    double temp_term1_y = 34459425.0 * temp_y_6 - 30405375.0 * temp_y_4 + 6081075.0 * temp_y_2 - 155925.0;
    double temp_term1_z = 34459425.0 * temp_z_6 - 30405375.0 * temp_z_4 + 6081075.0 * temp_z_2 - 155925.0;

    double temp_term2_x = -6081075.0 * temp_x_6 + 6081075.0 * temp_x_4 - 1403325.0 * temp_x_2 + 42525.0;
    double temp_term2_y = -6081075.0 * temp_y_6 + 6081075.0 * temp_y_4 - 1403325.0 * temp_x_2 + 42525.0;
    double temp_term2_z = -6081075.0 * temp_z_6 + 6081075.0 * temp_z_4 - 1403325.0 * temp_x_2 + 42525.0;

    double temp_miu1 = temp_y_3 * temp_term1_x + lambda_y0 * temp_term2_x;
    double temp_miu2 = temp_x_3 * temp_term1_y + lambda_x0 * temp_term2_y;
    double temp_miu3 = temp_z_3 * temp_term1_x + lambda_z0 * temp_term2_x;
    double temp_miu4 = temp_x_3 * temp_term1_z + lambda_x0 * temp_term2_z;
    double temp_miu5 = temp_z_3 * temp_term1_y + lambda_z0 * temp_term2_y;
    double temp_miu6 = temp_y_3 * temp_term1_z + lambda_y0 * temp_term2_z;

    double miu_9_5_1 = temp * temp_miu1;
    double miu_9_5_2 = temp * temp_miu2;
    double miu_9_5_3 = temp * temp_miu3;
    double miu_9_5_4 = temp * temp_miu4;
    double miu_9_5_5 = temp * temp_miu5;
    double miu_9_5_6 = temp * temp_miu6;

    value[0] = miu_9_5_1;
    value[1] = miu_9_5_2;
    value[2] = miu_9_5_3;
    value[3] = miu_9_5_4;
    value[4] = miu_9_5_5;
    value[5] = miu_9_5_6;
}


void calc_MCSH_9_6_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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

    double lambda_x0_5 = lambda_x0_4 * lambda_x0;
    double lambda_y0_5 = lambda_y0_4 * lambda_y0;
    double lambda_z0_5 = lambda_z0_4 * lambda_z0;

    double lambda_x0_6 = lambda_x0_5 * lambda_x0;
    double lambda_y0_6 = lambda_y0_5 * lambda_y0;
    double lambda_z0_6 = lambda_z0_5 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x_2 = lambda_x0_sqr + (1.0/(2.0*gamma));
    double temp_y_2 = lambda_y0_sqr + (1.0/(2.0*gamma));
    double temp_z_2 = lambda_z0_sqr + (1.0/(2.0*gamma));

    double C4_1 = 3.0 / gamma, C4_2 = 3.0 / (4.0 * gamma * gamma);
    double temp_x_4 = lambda_x0_4 + C4_1 * lambda_x0_sqr + C4_2;
    double temp_y_4 = lambda_y0_4 + C4_1 * lambda_y0_sqr + C4_2;
    double temp_z_4 = lambda_z0_4 + C4_1 * lambda_z0_sqr + C4_2;

    double C6_1 = 15.0 / (2.0 * gamma), C6_2 = 45.0 / (4.0 * gamma * gamma), C6_3 =  15.0 / (8.0 * gamma * gamma * gamma);
    double temp_x_6 = lambda_x0_6 + C6_1 * lambda_x0_4 + C6_2 * lambda_x0_sqr + C6_3;
    double temp_y_6 = lambda_y0_6 + C6_1 * lambda_y0_4 + C6_2 * lambda_y0_sqr + C6_3;
    double temp_z_6 = lambda_z0_6 + C6_1 * lambda_z0_4 + C6_2 * lambda_z0_sqr + C6_3;

    double temp_term1_x = 34459425.0 * temp_x_6 - 30405375.0 * temp_x_4 + 6081075.0 * temp_x_2 - 155925.0;
    double temp_term1_y = 34459425.0 * temp_y_6 - 30405375.0 * temp_y_4 + 6081075.0 * temp_y_2 - 155925.0;
    double temp_term1_z = 34459425.0 * temp_z_6 - 30405375.0 * temp_z_4 + 6081075.0 * temp_z_2 - 155925.0;

    double temp_term2_x = -2027025.0 * temp_x_6 + 2027025.0 * temp_x_4 - 467775.0 * temp_x_2 + 14175.0;
    double temp_term2_y = -2027025.0 * temp_y_6 + 2027025.0 * temp_y_4 - 467775.0 * temp_x_2 + 14175.0;
    double temp_term2_z = -2027025.0 * temp_z_6 + 2027025.0 * temp_z_4 - 467775.0 * temp_x_2 + 14175.0;

    double temp_miu1 = temp_y_2 * temp_term1_x + temp_term2_x;
    double temp_miu2 = temp_x_2 * temp_term1_y + temp_term2_y;
    double temp_miu3 = temp_z_2 * temp_term1_x + temp_term2_x;
    double temp_miu4 = temp_x_2 * temp_term1_z + temp_term2_z;
    double temp_miu5 = temp_z_2 * temp_term1_y + temp_term2_y;
    double temp_miu6 = temp_y_2 * temp_term1_z + temp_term2_z;

    double miu_9_6_1 = temp * lambda_z0 * temp_miu1;
    double miu_9_6_2 = temp * lambda_z0 * temp_miu2;
    double miu_9_6_3 = temp * lambda_y0 * temp_miu3;
    double miu_9_6_4 = temp * lambda_y0 * temp_miu4;
    double miu_9_6_5 = temp * lambda_x0 * temp_miu5;
    double miu_9_6_6 = temp * lambda_x0 * temp_miu6;

    value[0] = miu_9_6_1;
    value[1] = miu_9_6_2;
    value[2] = miu_9_6_3;
    value[3] = miu_9_6_4;
    value[4] = miu_9_6_5;
    value[5] = miu_9_6_6;
}

void calc_MCSH_9_7_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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

    double lambda_x0_5 = lambda_x0_4 * lambda_x0;
    double lambda_y0_5 = lambda_y0_4 * lambda_y0;
    double lambda_z0_5 = lambda_z0_4 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x_2 = lambda_x0_sqr + (1.0/(2.0*gamma));
    double temp_y_2 = lambda_y0_sqr + (1.0/(2.0*gamma));
    double temp_z_2 = lambda_z0_sqr + (1.0/(2.0*gamma));

    double C3_1 = 3.0 / (2.0 * gamma);
    double temp_x_3 = lambda_x0_3 + C3_1 * lambda_x0;
    double temp_y_3 = lambda_y0_3 + C3_1 * lambda_y0;
    double temp_z_3 = lambda_z0_3 + C3_1 * lambda_z0;

    double C4_1 = 3.0 / gamma, C4_2 = 3.0 / (4.0 * gamma * gamma);
    double temp_x_4 = lambda_x0_4 + C4_1 * lambda_x0_sqr + C4_2;
    double temp_y_4 = lambda_y0_4 + C4_1 * lambda_y0_sqr + C4_2;
    double temp_z_4 = lambda_z0_4 + C4_1 * lambda_z0_sqr + C4_2;

    double C5_1 = 5.0 / gamma, C5_2 = 15.0 / (4.0 * gamma * gamma);
    double temp_x_5 = lambda_x0_5 + C5_1 * lambda_x0_3 + C5_2 * lambda_x0;
    double temp_y_5 = lambda_y0_5 + C5_1 * lambda_y0_3 + C5_2 * lambda_y0;
    double temp_z_5 = lambda_z0_5 + C5_1 * lambda_z0_3 + C5_2 * lambda_z0;

    double temp_term1_x = 34459425.0 * temp_x_5 - 20270250.0 * temp_x_3 + 2027025.0 * lambda_x0;
    double temp_term1_y = 34459425.0 * temp_y_5 - 20270250.0 * temp_y_3 + 2027025.0 * lambda_y0;
    double temp_term1_z = 34459425.0 * temp_z_5 - 20270250.0 * temp_z_3 + 2027025.0 * lambda_z0;

    double temp_term2_x = -12162150.0 * temp_x_5 + 8108100.0 * temp_x_3 - 935550.0 * lambda_x0;
    double temp_term2_y = -12162150.0 * temp_y_5 + 8108100.0 * temp_y_3 - 935550.0 * lambda_y0;
    double temp_term2_z = -12162150.0 * temp_z_5 + 8108100.0 * temp_z_3 - 935550.0 * lambda_z0;

    double temp_term3_x = 405405.0 * temp_x_5 - 311850.0 * temp_x_3 + 425425.0 * lambda_x0;
    double temp_term3_y = 405405.0 * temp_y_5 - 311850.0 * temp_y_3 + 425425.0 * lambda_y0;
    double temp_term3_z = 405405.0 * temp_z_5 - 311850.0 * temp_z_3 + 425425.0 * lambda_z0;

    double temp_miu1 = temp_y_4 * temp_term1_x + temp_y_2 * temp_term2_x + temp_term3_x;
    double temp_miu2 = temp_x_4 * temp_term1_y + temp_x_2 * temp_term2_y + temp_term3_y;
    double temp_miu3 = temp_z_4 * temp_term1_x + temp_z_2 * temp_term2_x + temp_term3_x;
    double temp_miu4 = temp_x_4 * temp_term1_z + temp_x_2 * temp_term2_z + temp_term3_z;
    double temp_miu5 = temp_z_4 * temp_term1_y + temp_z_2 * temp_term2_y + temp_term3_y;
    double temp_miu6 = temp_y_4 * temp_term1_z + temp_y_2 * temp_term2_z + temp_term3_z;

    double miu_9_7_1 = temp * temp_miu1;
    double miu_9_7_2 = temp * temp_miu2;
    double miu_9_7_3 = temp * temp_miu3;
    double miu_9_7_4 = temp * temp_miu4;
    double miu_9_7_5 = temp * temp_miu5;
    double miu_9_7_6 = temp * temp_miu6;

    value[0] = miu_9_7_1;
    value[1] = miu_9_7_2;
    value[2] = miu_9_7_3;
    value[3] = miu_9_7_4;
    value[4] = miu_9_7_5;
    value[5] = miu_9_7_6;
}

void calc_MCSH_9_8_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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

    double lambda_x0_5 = lambda_x0_4 * lambda_x0;
    double lambda_y0_5 = lambda_y0_4 * lambda_y0;
    double lambda_z0_5 = lambda_z0_4 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);

    double temp = C1 * exp( C2 * r0_sqr);

    double C3_1 = 3.0 / (2.0 * gamma);
    double temp_x_3 = lambda_x0_3 + C3_1 * lambda_x0;
    double temp_y_3 = lambda_y0_3 + C3_1 * lambda_y0;
    double temp_z_3 = lambda_z0_3 + C3_1 * lambda_z0;

    double C5_1 = 5.0 / gamma, C5_2 = 15.0 / (4.0 * gamma * gamma);
    double temp_x_5 = lambda_x0_5 + C5_1 * lambda_x0_3 + C5_2 * lambda_x0;
    double temp_y_5 = lambda_y0_5 + C5_1 * lambda_y0_3 + C5_2 * lambda_y0;
    double temp_z_5 = lambda_z0_5 + C5_1 * lambda_z0_3 + C5_2 * lambda_z0;

    double temp_term1_x = 34459425.0 * temp_x_5 - 20270250.0 * temp_x_3 + 2027025.0 * lambda_x0;
    double temp_term1_y = 34459425.0 * temp_y_5 - 20270250.0 * temp_y_3 + 2027025.0 * lambda_y0;
    double temp_term1_z = 34459425.0 * temp_z_5 - 20270250.0 * temp_z_3 + 2027025.0 * lambda_z0;

    double temp_term2_x = -6081075.0 * temp_x_5 + 4054050.0 * temp_x_3 - 467775.0 * lambda_x0;
    double temp_term2_y = -6081075.0 * temp_y_5 + 4054050.0 * temp_y_3 - 467775.0 * lambda_x0;
    double temp_term2_z = -6081075.0 * temp_z_5 + 4054050.0 * temp_z_3 - 467775.0 * lambda_x0;

    double temp_miu1 = temp_y_3 * temp_term1_x + lambda_y0 * temp_term2_x;
    double temp_miu2 = temp_x_3 * temp_term1_y + lambda_x0 * temp_term2_y;
    double temp_miu3 = temp_z_3 * temp_term1_x + lambda_z0 * temp_term2_x;
    double temp_miu4 = temp_x_3 * temp_term1_z + lambda_x0 * temp_term2_z;
    double temp_miu5 = temp_z_3 * temp_term1_y + lambda_z0 * temp_term2_y;
    double temp_miu6 = temp_y_3 * temp_term1_z + lambda_y0 * temp_term2_z;

    double miu_9_8_1 = temp * lambda_z0 * temp_miu1;
    double miu_9_8_2 = temp * lambda_z0 * temp_miu2;
    double miu_9_8_3 = temp * lambda_y0 * temp_miu3;
    double miu_9_8_4 = temp * lambda_y0 * temp_miu4;
    double miu_9_8_5 = temp * lambda_x0 * temp_miu5;
    double miu_9_8_6 = temp * lambda_x0 * temp_miu6;

    value[0] = miu_9_8_1;
    value[1] = miu_9_8_2;
    value[2] = miu_9_8_3;
    value[3] = miu_9_8_4;
    value[4] = miu_9_8_5;
    value[5] = miu_9_8_6;
}

void calc_MCSH_9_9_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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

    double lambda_x0_5 = lambda_x0_4 * lambda_x0;
    double lambda_y0_5 = lambda_y0_4 * lambda_y0;
    double lambda_z0_5 = lambda_z0_4 * lambda_z0;

    double gamma = calc_gamma(alpha, beta);

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x_2 = lambda_x0_sqr + (1.0/(2.0*gamma));
    double temp_y_2 = lambda_y0_sqr + (1.0/(2.0*gamma));
    double temp_z_2 = lambda_z0_sqr + (1.0/(2.0*gamma));

    double C3_1 = 3.0 / (2.0 * gamma);
    double temp_x_3 = lambda_x0_3 + C3_1 * lambda_x0;
    double temp_y_3 = lambda_y0_3 + C3_1 * lambda_y0;
    double temp_z_3 = lambda_z0_3 + C3_1 * lambda_z0;

    double C5_1 = 5.0 / gamma, C5_2 = 15.0 / (4.0 * gamma * gamma);
    double temp_x_5 = lambda_x0_5 + C5_1 * lambda_x0_3 + C5_2 * lambda_x0;
    double temp_y_5 = lambda_y0_5 + C5_1 * lambda_y0_3 + C5_2 * lambda_y0;
    double temp_z_5 = lambda_z0_5 + C5_1 * lambda_z0_3 + C5_2 * lambda_z0;


    double temp_term1_x = 34459425.0 * temp_x_5 - 20270250.0 * temp_x_3 + 2027025.0 * lambda_x0;
    double temp_term1_y = 34459425.0 * temp_y_5 - 20270250.0 * temp_y_3 + 2027025.0 * lambda_y0;
    double temp_term1_z = 34459425.0 * temp_z_5 - 20270250.0 * temp_z_3 + 2027025.0 * lambda_z0;

    double temp_term2_x = -2027025.0 * temp_x_5 + 1351350.0 * temp_x_3 - 155925.0 * lambda_x0;
    double temp_term2_y = -2027025.0 * temp_y_5 + 1351350.0 * temp_y_3 - 155925.0 * lambda_y0;
    double temp_term2_z = -2027025.0 * temp_z_5 + 1351350.0 * temp_z_3 - 155925.0 * lambda_z0;

    double temp_term3_x = 135135.0 * temp_x_5 - 103950.0 * temp_x_3 + 14175.0 * lambda_x0;
    double temp_term3_y = 135135.0 * temp_y_5 - 103950.0 * temp_y_3 + 14175.0 * lambda_y0;
    double temp_term3_z = 135135.0 * temp_z_5 - 103950.0 * temp_z_3 + 14175.0 * lambda_z0;

    double temp_miu1 = temp_y_2 * temp_z_2 * temp_term1_x + (temp_y_2 + temp_z_2) * temp_term2_x + temp_term3_x;
    double temp_miu2 = temp_x_2 * temp_z_2 * temp_term1_y + (temp_x_2 + temp_z_2) * temp_term2_y + temp_term3_y;
    double temp_miu3 = temp_x_2 * temp_y_2 * temp_term1_z + (temp_x_2 + temp_y_2) * temp_term2_z + temp_term3_z;

    double miu_9_9_1 = temp * temp_miu1;
    double miu_9_9_2 = temp * temp_miu2;
    double miu_9_9_3 = temp * temp_miu3;

    value[0] = miu_9_9_1;
    value[1] = miu_9_9_2;
    value[2] = miu_9_9_3;
}


void calc_MCSH_9_10_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x_2 = lambda_x0_sqr + (1.0/(2.0*gamma));
    double temp_y_2 = lambda_y0_sqr + (1.0/(2.0*gamma));
    double temp_z_2 = lambda_z0_sqr + (1.0/(2.0*gamma));

    double C4_1 = 3.0 / gamma, C4_2 = 3.0 / (4.0 * gamma * gamma);
    double temp_x_4 = lambda_x0_4 + C4_1 * lambda_x0_sqr + C4_2;
    double temp_y_4 = lambda_y0_4 + C4_1 * lambda_y0_sqr + C4_2;
    double temp_z_4 = lambda_z0_4 + C4_1 * lambda_z0_sqr + C4_2;

    double temp_term1_x = 34459425.0 * temp_x_4 - 12162150.0 * temp_x_2 + 405405.0;
    double temp_term1_y = 34459425.0 * temp_y_4 - 12162150.0 * temp_y_2 + 405405.0;

    double temp_term2_x = -12162150.0 * temp_x_4 + 4864860.0 * temp_x_2 - 187110.0;
    double temp_term2_y = -12162150.0 * temp_y_4 + 4864860.0 * temp_y_2 - 187110.0;

    double temp_term3_x = 405405.0 * temp_x_4 - 187110.0 * temp_x_2 + 8505.0;
    double temp_term3_y = 405405.0 * temp_y_4 - 187110.0 * temp_y_2 + 8505.0;

    double temp_miu1 = temp_y_4 * temp_term1_x + temp_y_2 * temp_term2_x + temp_term3_x;
    double temp_miu2 = temp_z_4 * temp_term1_x + temp_z_2 * temp_term2_x + temp_term3_x;
    double temp_miu3 = temp_z_4 * temp_term1_y + temp_z_2 * temp_term2_y + temp_term3_y;

    double miu_9_10_1 = temp * lambda_z0 * temp_miu1;
    double miu_9_10_2 = temp * lambda_y0 * temp_miu2;
    double miu_9_10_3 = temp * lambda_x0 * temp_miu3;

    value[0] = miu_9_10_1;
    value[1] = miu_9_10_2;
    value[2] = miu_9_10_3;
}

void calc_MCSH_9_11_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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

    double temp = C1 * exp( C2 * r0_sqr);

    double temp_x_2 = lambda_x0_sqr + (1.0/(2.0*gamma));
    double temp_y_2 = lambda_y0_sqr + (1.0/(2.0*gamma));
    double temp_z_2 = lambda_z0_sqr + (1.0/(2.0*gamma));

    double C3_1 = 3.0 / (2.0 * gamma);
    double temp_x_3 = lambda_x0_3 + C3_1 * lambda_x0;
    double temp_y_3 = lambda_y0_3 + C3_1 * lambda_y0;
    double temp_z_3 = lambda_z0_3 + C3_1 * lambda_z0;

    double C4_1 = 3.0 / gamma, C4_2 = 3.0 / (4.0 * gamma * gamma);
    double temp_x_4 = lambda_x0_4 + C4_1 * lambda_x0_sqr + C4_2;
    double temp_y_4 = lambda_y0_4 + C4_1 * lambda_y0_sqr + C4_2;
    double temp_z_4 = lambda_z0_4 + C4_1 * lambda_z0_sqr + C4_2;

    double temp_term1_x = 34459425.0 * temp_x_4 - 12162150.0 * temp_x_2 + 405405.0;
    double temp_term1_y = 34459425.0 * temp_y_4 - 12162150.0 * temp_y_2 + 405405.0;
    double temp_term1_z = 34459425.0 * temp_z_4 - 12162150.0 * temp_z_2 + 405405.0;

    double temp_term2_x = -2027025.0 * temp_x_4 + 810810.0 * temp_x_2 - 31185.0;
    double temp_term2_y = -2027025.0 * temp_y_4 + 810810.0 * temp_y_2 - 31185.0;
    double temp_term2_z = -2027025.0 * temp_z_4 + 810810.0 * temp_z_2 - 31185.0;

    double temp_term3_x = -6081075.0 * temp_x_4 + 2432430.0 * temp_x_2 - 93555.0;
    double temp_term3_y = -6081075.0 * temp_y_4 + 2432430.0 * temp_y_2 - 93555.0;
    double temp_term3_z = -6081075.0 * temp_z_4 + 2432430.0 * temp_z_2 - 93555.0;

    double temp_term4_x = 405405.0 * temp_x_4 - 187110.0 * temp_x_2 + 8505.0;
    double temp_term4_y = 405405.0 * temp_y_4 - 187110.0 * temp_y_2 + 8505.0;
    double temp_term4_z = 405405.0 * temp_z_4 - 187110.0 * temp_z_2 + 8505.0;


    double temp_miu1 = temp_y_3 * temp_z_2 * temp_term1_x + temp_y_3 * temp_term2_x + lambda_y0 * temp_z_2 * temp_term3_x + lambda_y0 * temp_term4_x;
    double temp_miu2 = temp_x_3 * temp_z_2 * temp_term1_y + temp_x_3 * temp_term2_y + lambda_x0 * temp_z_2 * temp_term3_y + lambda_x0 * temp_term4_y;
    double temp_miu3 = temp_z_3 * temp_y_2 * temp_term1_x + temp_z_3 * temp_term2_x + lambda_z0 * temp_y_2 * temp_term3_x + lambda_z0 * temp_term4_x;
    double temp_miu4 = temp_x_3 * temp_y_2 * temp_term1_z + temp_x_3 * temp_term2_z + lambda_x0 * temp_y_2 * temp_term3_z + lambda_x0 * temp_term4_z;
    double temp_miu5 = temp_z_3 * temp_x_2 * temp_term1_y + temp_z_3 * temp_term2_y + lambda_z0 * temp_x_2 * temp_term3_y + lambda_z0 * temp_term4_y;
    double temp_miu6 = temp_y_3 * temp_x_2 * temp_term1_z + temp_y_3 * temp_term2_z + lambda_y0 * temp_x_2 * temp_term3_z + lambda_y0 * temp_term4_z;

    double miu_9_11_1 = temp * temp_miu1;
    double miu_9_11_2 = temp * temp_miu2;
    double miu_9_11_3 = temp * temp_miu3;
    double miu_9_11_4 = temp * temp_miu4;
    double miu_9_11_5 = temp * temp_miu5;
    double miu_9_11_6 = temp * temp_miu6;

    value[0] = miu_9_11_1;
    value[1] = miu_9_11_2;
    value[2] = miu_9_11_3;
    value[3] = miu_9_11_4;
    value[4] = miu_9_11_5;
    value[5] = miu_9_11_6;
}


void calc_MCSH_9_12_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);

    double lambda = calc_lambda(alpha, beta);

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

    double temp = C1 * exp( C2 * r0_sqr);

    double C3_1 = 3.0 / (2.0 * gamma);
    double temp_x_3 = lambda_x0_3 + C3_1 * lambda_x0;
    double temp_y_3 = lambda_y0_3 + C3_1 * lambda_y0;
    double temp_z_3 = lambda_z0_3 + C3_1 * lambda_z0;

    double t1 = 34459425.0 * temp_x_3 * temp_y_3 * temp_z_3;
    double t2 = -6081075.0 * temp_x_3 * temp_y_3 * lambda_z0;
    double t3 = -6081075.0 * temp_x_3 * lambda_y0 * temp_z_3;
    double t4 = 1216215.0 * temp_x_3 * lambda_y0 * lambda_z0;
    double t5 = -6081075.0 * lambda_x0 * temp_y_3 * temp_z_3;
    double t6 = 1216215.0 * lambda_x0 * temp_y_3 * lambda_z0;
    double t7 = 1216215.0 * lambda_x0 * lambda_y0 * temp_z_3;
    double t8 = -280665.0* lambda_x0 * lambda_y0 * lambda_z0;

    double sum_ts = t1 + t2 + t3 + t4 + t5 + t6 + t7 + t8;

    double m_9_12 = temp * sum_ts;

    value[0] = m_9_12;
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
    } else if (mcsh_order == 6) {
        if (group_num == 1) {
            result = calc_MCSH_6_1;
        } else if (group_num == 2){
            result = calc_MCSH_6_2;
        } else if (group_num == 3){
            result = calc_MCSH_6_3;
        } else if (group_num == 4){
            result = calc_MCSH_6_4;
        } else if (group_num == 5){
            result = calc_MCSH_6_5;
        } else if (group_num == 6){
            result = calc_MCSH_6_6;
        } else if (group_num == 7){
            result = calc_MCSH_6_7;
        }
    } else if (mcsh_order == 7) {
        if (group_num == 1) {
            result = calc_MCSH_7_1;
        } else if (group_num == 2){
            result = calc_MCSH_7_2;
        } else if (group_num == 3){
            result = calc_MCSH_7_3;
        } else if (group_num == 4){
            result = calc_MCSH_7_4;
        } else if (group_num == 5){
            result = calc_MCSH_7_5;
        } else if (group_num == 6){
            result = calc_MCSH_7_6;
        } else if (group_num == 7){
            result = calc_MCSH_7_7;
        } else if (group_num == 8){
            result = calc_MCSH_7_8;
        }
    } else if (mcsh_order == 8) {
        if (group_num == 1) {
            result = calc_MCSH_8_1;
        } else if (group_num == 2){
            result = calc_MCSH_8_2;
        } else if (group_num == 3){
            result = calc_MCSH_8_3;
        } else if (group_num == 4){
            result = calc_MCSH_8_4;
        } else if (group_num == 5){
            result = calc_MCSH_8_5;
        } else if (group_num == 6){
            result = calc_MCSH_8_6;
        } else if (group_num == 7){
            result = calc_MCSH_8_7;
        } else if (group_num == 8){
            result = calc_MCSH_8_8;
        } else if (group_num == 9){
            result = calc_MCSH_8_9;
        } else if (group_num == 10){
            result = calc_MCSH_8_10;
        }
    } else if (mcsh_order == 9) {
        if (group_num == 1) {
            result = calc_MCSH_9_1;
        } else if (group_num == 2){
            result = calc_MCSH_9_2;
        } else if (group_num == 3){
            result = calc_MCSH_9_3;
        } else if (group_num == 4){
            result = calc_MCSH_9_4;
        } else if (group_num == 5){
            result = calc_MCSH_9_5;
        } else if (group_num == 6){
            result = calc_MCSH_9_6;
        } else if (group_num == 7){
            result = calc_MCSH_9_7;
        } else if (group_num == 8){
            result = calc_MCSH_9_8;
        } else if (group_num == 9){
            result = calc_MCSH_9_9;
        } else if (group_num == 10){
            result = calc_MCSH_9_10;
        } else if (group_num == 11){
            result = calc_MCSH_9_11;
        } else if (group_num == 12){
            result = calc_MCSH_9_12;
        }
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
    } else if (mcsh_order == 9) {
        if (group_num == 1) {
            result = calc_MCSH_9_1_noderiv;
        } else if (group_num == 2){
            result = calc_MCSH_9_2_noderiv;
        } else if (group_num == 3){
            result = calc_MCSH_9_3_noderiv;
        } else if (group_num == 4){
            result = calc_MCSH_9_4_noderiv;
        } else if (group_num == 5){
            result = calc_MCSH_9_5_noderiv;
        } else if (group_num == 6){
            result = calc_MCSH_9_6_noderiv;
        } else if (group_num == 7){
            result = calc_MCSH_9_7_noderiv;
        } else if (group_num == 8){
            result = calc_MCSH_9_8_noderiv;
        } else if (group_num == 9){
            result = calc_MCSH_9_9_noderiv;
        } else if (group_num == 10){
            result = calc_MCSH_9_10_noderiv;
        } else if (group_num == 11){
            result = calc_MCSH_9_11_noderiv;
        } else if (group_num == 12){
            result = calc_MCSH_9_12_noderiv;
        }
    }

    return result;
}
