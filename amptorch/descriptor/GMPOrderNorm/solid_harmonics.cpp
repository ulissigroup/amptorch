#include <math.h>
#include "solid_harmonics.h"

SolidGMPFunctionNoderiv get_solid_mcsh_function_noderiv(int mcsh_order, int group_num)
{
    SolidGMPFunctionNoderiv result;

    if (mcsh_order == 0) {
        if (group_num == 1) {
            result = calc_solid_MCSH_0_1_noderiv;
        }
    } else if (mcsh_order == 1) {
        if (group_num == 1) {
            result = calc_solid_MCSH_1_1_noderiv;
        }
    } else if (mcsh_order == 2) {
        if (group_num == 1) {
            result = calc_solid_MCSH_2_1_noderiv;
        } else if (group_num == 2){
            result = calc_solid_MCSH_2_2_noderiv;
        }
    } else if (mcsh_order == 3) {
        if (group_num == 1) {
            result = calc_solid_MCSH_3_1_noderiv;
        } else if (group_num == 2){
            result = calc_solid_MCSH_3_2_noderiv;
        } else if (group_num == 3){
            result = calc_solid_MCSH_3_3_noderiv;
        }
    } else if (mcsh_order == 4) {
        if (group_num == 1) {
            result = calc_solid_MCSH_4_1_noderiv;
        } else if (group_num == 2){
            result = calc_solid_MCSH_4_2_noderiv;
        } else if (group_num == 3){
            result = calc_solid_MCSH_4_3_noderiv;
        } else if (group_num == 4){
            result = calc_solid_MCSH_4_4_noderiv;
        }
    } else if (mcsh_order == 5) {
        if (group_num == 1) {
            result = calc_solid_MCSH_5_1_noderiv;
        } else if (group_num == 2){
            result = calc_solid_MCSH_5_2_noderiv;
        } else if (group_num == 3){
            result = calc_solid_MCSH_5_3_noderiv;
        } else if (group_num == 4){
            result = calc_solid_MCSH_5_4_noderiv;
        } else if (group_num == 5){
            result = calc_solid_MCSH_5_5_noderiv;
        }
    } else if (mcsh_order == 6) {
        if (group_num == 1) {
            result = calc_solid_MCSH_6_1_noderiv;
        } else if (group_num == 2){
            result = calc_solid_MCSH_6_2_noderiv;
        } else if (group_num == 3){
            result = calc_solid_MCSH_6_3_noderiv;
        } else if (group_num == 4){
            result = calc_solid_MCSH_6_4_noderiv;
        } else if (group_num == 5){
            result = calc_solid_MCSH_6_5_noderiv;
        } else if (group_num == 6){
            result = calc_solid_MCSH_6_6_noderiv;
        } else if (group_num == 7){
            result = calc_solid_MCSH_6_7_noderiv;
        }
    } else if (mcsh_order == 7) {
        if (group_num == 1) {
            result = calc_solid_MCSH_7_1_noderiv;
        } else if (group_num == 2){
            result = calc_solid_MCSH_7_2_noderiv;
        } else if (group_num == 3){
            result = calc_solid_MCSH_7_3_noderiv;
        } else if (group_num == 4){
            result = calc_solid_MCSH_7_4_noderiv;
        } else if (group_num == 5){
            result = calc_solid_MCSH_7_5_noderiv;
        } else if (group_num == 6){
            result = calc_solid_MCSH_7_6_noderiv;
        } else if (group_num == 7){
            result = calc_solid_MCSH_7_7_noderiv;
        } else if (group_num == 8){
            result = calc_solid_MCSH_7_8_noderiv;
        }
    } else if (mcsh_order == 8) {
        if (group_num == 1) {
            result = calc_solid_MCSH_8_1_noderiv;
        } else if (group_num == 2){
            result = calc_solid_MCSH_8_2_noderiv;
        } else if (group_num == 3){
            result = calc_solid_MCSH_8_3_noderiv;
        } else if (group_num == 4){
            result = calc_solid_MCSH_8_4_noderiv;
        } else if (group_num == 5){
            result = calc_solid_MCSH_8_5_noderiv;
        } else if (group_num == 6){
            result = calc_solid_MCSH_8_6_noderiv;
        } else if (group_num == 7){
            result = calc_solid_MCSH_8_7_noderiv;
        } else if (group_num == 8){
            result = calc_solid_MCSH_8_8_noderiv;
        } else if (group_num == 9){
            result = calc_solid_MCSH_8_9_noderiv;
        } else if (group_num == 10){
            result = calc_solid_MCSH_8_10_noderiv;
        }
    } else if (mcsh_order == 9) {
        if (group_num == 1) {
            result = calc_solid_MCSH_9_1_noderiv;
        } else if (group_num == 2){
            result = calc_solid_MCSH_9_2_noderiv;
        } else if (group_num == 3){
            result = calc_solid_MCSH_9_3_noderiv;
        } else if (group_num == 4){
            result = calc_solid_MCSH_9_4_noderiv;
        } else if (group_num == 5){
            result = calc_solid_MCSH_9_5_noderiv;
        } else if (group_num == 6){
            result = calc_solid_MCSH_9_6_noderiv;
        } else if (group_num == 7){
            result = calc_solid_MCSH_9_7_noderiv;
        } else if (group_num == 8){
            result = calc_solid_MCSH_9_8_noderiv;
        } else if (group_num == 9){
            result = calc_solid_MCSH_9_9_noderiv;
        } else if (group_num == 10){
            result = calc_solid_MCSH_9_10_noderiv;
        } else if (group_num == 11){
            result = calc_solid_MCSH_9_11_noderiv;
        } else if (group_num == 12){
            result = calc_solid_MCSH_9_12_noderiv;
        }
    }

    return result;
}


void calc_solid_MCSH_0_1_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double m_0_1 = C1 * exp( C2 * r0_sqr);

    value[0] = m_0_1;
}

void calc_solid_MCSH_1_1_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
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

void calc_solid_MCSH_2_1_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double term_x = (2.0 * P2x) - (P2y + P2z);
    double term_y = (2.0 * P2y) - (P2x + P2z);
    double term_z = (2.0 * P2z) - (P2x + P2y);

    double miu_1 = temp * term_x;
    double miu_2 = temp * term_y;
    double miu_3 = temp * term_z;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
}

void calc_solid_MCSH_2_2_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
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

void calc_solid_MCSH_3_1_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

    double P1x = P1(lambda, x0, gamma);
    double P1y = P1(lambda, y0, gamma);
    double P1z = P1(lambda, z0, gamma);

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double P3x = P3(lambda, x0, gamma);
    double P3y = P3(lambda, y0, gamma);
    double P3z = P3(lambda, z0, gamma);

    double term_x = (6.0 * P3x) - (9.0 * P1x * (P2y + P2z));
    double term_y = (6.0 * P3y) - (9.0 * P1y * (P2x + P2z));
    double term_z = (6.0 * P3z) - (9.0 * P1z * (P2x + P2y));

    double miu_1 = temp * term_x;
    double miu_2 = temp * term_y;
    double miu_3 = temp * term_z;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
}

void calc_solid_MCSH_3_2_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

    double P1x = P1(lambda, x0, gamma);
    double P1y = P1(lambda, y0, gamma);
    double P1z = P1(lambda, z0, gamma);

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double P3x = P3(lambda, x0, gamma);
    double P3y = P3(lambda, y0, gamma);
    double P3z = P3(lambda, z0, gamma);

    double term_1 = (12.0 * P2x * P1y) - (3.0 * P3y) - (3.0 * P1y * P2z);
    double term_2 = (12.0 * P2y * P1x) - (3.0 * P3x) - (3.0 * P1x * P2z);
    double term_3 = (12.0 * P2x * P1z) - (3.0 * P3z) - (3.0 * P1z * P2y);
    double term_4 = (12.0 * P2z * P1x) - (3.0 * P3x) - (3.0 * P1x * P2y);
    double term_5 = (12.0 * P2y * P1z) - (3.0 * P3z) - (3.0 * P1z * P2x);
    double term_6 = (12.0 * P2z * P1y) - (3.0 * P3y) - (3.0 * P1y * P2x);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;
}

void calc_solid_MCSH_3_3_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double lambda = calc_lambda(alpha, beta);

    double temp =  C1 * exp( C2 * r0_sqr) * lambda * lambda * lambda * 15.0 ;
    double m_3_3 = temp * x0 * y0 * z0;

    value[0] = m_3_3;
}

void calc_solid_MCSH_4_1_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double P4x = P4(lambda, x0, gamma);
    double P4y = P4(lambda, y0, gamma);
    double P4z = P4(lambda, z0, gamma);

    double term_x = (24.0 * P4x) + (9.0 * (P4y + P4z)) - (72.0 * P2x * (P2y + P2z)) + (18.0 * P2y * P2z);
    double term_y = (24.0 * P4y) + (9.0 * (P4x + P4z)) - (72.0 * P2y * (P2x + P2z)) + (18.0 * P2x * P2z);
    double term_z = (24.0 * P4z) + (9.0 * (P4x + P4y)) - (72.0 * P2z * (P2x + P2y)) + (18.0 * P2x * P2y);

    double miu_1 = temp * term_x;
    double miu_2 = temp * term_y;
    double miu_3 = temp * term_z;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
}

void calc_solid_MCSH_4_2_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

    double P1x = P1(lambda, x0, gamma);
    double P1y = P1(lambda, y0, gamma);
    double P1z = P1(lambda, z0, gamma);

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double P3x = P3(lambda, x0, gamma);
    double P3y = P3(lambda, y0, gamma);
    double P3z = P3(lambda, z0, gamma);

    double term_1 = (60.0 * P3x * P1y) - (45.0 * (P1x * P3y + P1x * P1y * P2z));
    double term_2 = (60.0 * P3y * P1x) - (45.0 * (P1y * P3x + P1y * P1x * P2z));
    double term_3 = (60.0 * P3x * P1z) - (45.0 * (P1x * P3z + P1x * P1z * P2y));
    double term_4 = (60.0 * P3z * P1x) - (45.0 * (P1z * P3x + P1z * P1x * P2y));
    double term_5 = (60.0 * P3y * P1z) - (45.0 * (P1y * P3z + P1y * P1z * P2x));
    double term_6 = (60.0 * P3z * P1y) - (45.0 * (P1z * P3y + P1z * P1y * P2x));

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;
}



void calc_solid_MCSH_4_3_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double P4x = P4(lambda, x0, gamma);
    double P4y = P4(lambda, y0, gamma);
    double P4z = P4(lambda, z0, gamma);

    double term_1 = (81.0 * P2x * P2y) - (12.0 * (P4x + P4y)) + (3.0 * P4z) - (9.0 * P2z * (P2x + P2y));
    double term_2 = (81.0 * P2x * P2z) - (12.0 * (P4x + P4z)) + (3.0 * P4y) - (9.0 * P2y * (P2x + P2z));
    double term_3 = (81.0 * P2y * P2z) - (12.0 * (P4y + P4z)) + (3.0 * P4x) - (9.0 * P2x * (P2y + P2z));

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
}

void calc_solid_MCSH_4_4_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

    double P1x = P1(lambda, x0, gamma);
    double P1y = P1(lambda, y0, gamma);
    double P1z = P1(lambda, z0, gamma);

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double P3x = P3(lambda, x0, gamma);
    double P3y = P3(lambda, y0, gamma);
    double P3z = P3(lambda, z0, gamma);

    double term_1 = 90.0 * P2x * P1y * P1z - 15.0 * (P3y * P1z + P1y * P3z);
    double term_2 = 90.0 * P2y * P1x * P1z - 15.0 * (P3x * P1z + P1x * P3z);
    double term_3 = 90.0 * P2z * P1x * P1y - 15.0 * (P3x * P1y + P1x * P3y);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
}

void calc_solid_MCSH_5_1_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double term_x = (120.0 * P5x) - (600.0 * P3x * (P2y + P2z)) + (225.0 * P1x * (P4y + P4z)) + (450.0 * P1x * P2y * P2z);
    double term_y = (120.0 * P5y) - (600.0 * P3y * (P2x + P2z)) + (225.0 * P1y * (P4x + P4z)) + (450.0 * P1y * P2x * P2z);
    double term_z = (120.0 * P5z) - (600.0 * P3z * (P2x + P2y)) + (225.0 * P1z * (P4x + P4y)) + (450.0 * P1z * P2x * P2y);

    double miu_1 = temp * term_x;
    double miu_2 = temp * term_y;
    double miu_3 = temp * term_z;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
}

void calc_solid_MCSH_5_2_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double term_1 = (360.0 * P4x * P1y) - (540.0 * P2x * (P3y + P1y * P2z)) + (45.0 * (P5y + P1y * P4z)) + (90.0 * P3y * P2z);
    double term_2 = (360.0 * P4y * P1x) - (540.0 * P2y * (P3x + P1x * P2z)) + (45.0 * (P5x + P1x * P4z)) + (90.0 * P3x * P2z);
    double term_3 = (360.0 * P4x * P1z) - (540.0 * P2x * (P3z + P1z * P2y)) + (45.0 * (P5z + P1z * P4y)) + (90.0 * P3z * P2y);
    double term_4 = (360.0 * P4z * P1x) - (540.0 * P2z * (P3x + P1x * P2y)) + (45.0 * (P5x + P1x * P4y)) + (90.0 * P3x * P2y);
    double term_5 = (360.0 * P4y * P1z) - (540.0 * P2y * (P3z + P1z * P2x)) + (45.0 * (P5z + P1z * P4x)) + (90.0 * P3z * P2x);
    double term_6 = (360.0 * P4z * P1y) - (540.0 * P2z * (P3y + P1y * P2x)) + (45.0 * (P5y + P1y * P4x)) + (90.0 * P3y * P2x);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;
}



void calc_solid_MCSH_5_3_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    // double r0_sqr = x0*x0 + y0*y0 + z0*z0;
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double term_1 = (615.0 * P3x * P2y) - (60.0 * P5x) - (15.0 * P3x * P2z) - (270.0 * P1x * P4y) + (45.0 * P1x * P4z) - (225.0 * P1x * P2y * P2z);
    double term_2 = (615.0 * P3y * P2x) - (60.0 * P5y) - (15.0 * P3y * P2z) - (270.0 * P1y * P4x) + (45.0 * P1y * P4z) - (225.0 * P1y * P2x * P2z);
    double term_3 = (615.0 * P3x * P2z) - (60.0 * P5x) - (15.0 * P3x * P2y) - (270.0 * P1x * P4z) + (45.0 * P1x * P4y) - (225.0 * P1x * P2z * P2y);
    double term_4 = (615.0 * P3z * P2x) - (60.0 * P5z) - (15.0 * P3z * P2y) - (270.0 * P1z * P4x) + (45.0 * P1z * P4y) - (225.0 * P1z * P2x * P2y);
    double term_5 = (615.0 * P3y * P2z) - (60.0 * P5y) - (15.0 * P3y * P2x) - (270.0 * P1y * P4z) + (45.0 * P1y * P4x) - (225.0 * P1y * P2z * P2x);
    double term_6 = (615.0 * P3z * P2y) - (60.0 * P5z) - (15.0 * P3z * P2x) - (270.0 * P1z * P4y) + (45.0 * P1z * P4x) - (225.0 * P1z * P2y * P2x);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;
}

void calc_solid_MCSH_5_4_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

    double P1x = P1(lambda, x0, gamma);
    double P1y = P1(lambda, y0, gamma);
    double P1z = P1(lambda, z0, gamma);

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double P3x = P3(lambda, x0, gamma);
    double P3y = P3(lambda, y0, gamma);
    double P3z = P3(lambda, z0, gamma);

    double term_1 = (630.0 * P3x * P1y * P1z) - (315 * P1x * (P3y * P1z + P1y * P3z));
    double term_2 = (630.0 * P3y * P1x * P1z) - (315 * P1y * (P3x * P1z + P1x * P3z));
    double term_3 = (630.0 * P3z * P1x * P1y) - (315 * P1z * (P3x * P1y + P1x * P3y));

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
}

void calc_solid_MCSH_5_5_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double term_1 = (765.0 * P2x * P2y * P1z) - (90.0 * P1z * (P4x + P4y)) - (75.0 * P3z * (P2x + P2y)) + (15.0 * P5z);
    double term_2 = (765.0 * P2x * P2z * P1y) - (90.0 * P1y * (P4x + P4z)) - (75.0 * P3y * (P2x + P2z)) + (15.0 * P5y);
    double term_3 = (765.0 * P2y * P2z * P1x) - (90.0 * P1x * (P4y + P4z)) - (75.0 * P3x * (P2y + P2z)) + (15.0 * P5x);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
}

void calc_solid_MCSH_6_1_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double P4x = P4(lambda, x0, gamma);
    double P4y = P4(lambda, y0, gamma);
    double P4z = P4(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double term_1 = (720.0 * P6x) - (5400.0 * P4x * P2y) - (5400.0 * P4x * P2z) + (4050.0 * P2x * P4y) + (8100.0 * P2x * P2y * P2z) + (4050.0 * P2x * P4z) - (225.0 * P6y) - (675.0 * P4y * P2z) - (675.0 * P2y * P4z) - (225.0 * P6z);
    double term_2 = (720.0 * P6y) - (5400.0 * P4y * P2x) - (5400.0 * P4y * P2z) + (4050.0 * P2y * P4x) + (8100.0 * P2y * P2x * P2z) + (4050.0 * P2y * P4z) - (225.0 * P6x) - (675.0 * P4x * P2z) - (675.0 * P2x * P4z) - (225.0 * P6z);
    double term_3 = (720.0 * P6z) - (5400.0 * P4z * P2x) - (5400.0 * P4z * P2y) + (4050.0 * P2z * P4x) + (8100.0 * P2z * P2x * P2y) + (4050.0 * P2z * P4y) - (225.0 * P6x) - (675.0 * P4x * P2y) - (675.0 * P2x * P4y) - (225.0 * P6y);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;

}

void calc_solid_MCSH_6_2_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double term_1 = (2520.0 * P5x * P1y) - (6300.0 * P3x * P3y) - (6300.0 * P3x * P1y * P2z) + (1575.0 * P1x * P5y) + (3150.0 * P1x * P3y * P2z) + (1575.0 * P1x * P1y * P4z);
    double term_2 = (2520.0 * P5y * P1x) - (6300.0 * P3y * P3x) - (6300.0 * P3y * P1x * P2z) + (1575.0 * P1y * P5x) + (3150.0 * P1y * P3x * P2z) + (1575.0 * P1y * P1x * P4z);
    double term_3 = (2520.0 * P5x * P1z) - (6300.0 * P3x * P3z) - (6300.0 * P3x * P1z * P2y) + (1575.0 * P1x * P5z) + (3150.0 * P1x * P3z * P2y) + (1575.0 * P1x * P1z * P4y);
    double term_4 = (2520.0 * P5z * P1x) - (6300.0 * P3z * P3x) - (6300.0 * P3z * P1x * P2y) + (1575.0 * P1z * P5x) + (3150.0 * P1z * P3x * P2y) + (1575.0 * P1z * P1x * P4y);
    double term_5 = (2520.0 * P5y * P1z) - (6300.0 * P3y * P3z) - (6300.0 * P3y * P1z * P2x) + (1575.0 * P1y * P5z) + (3150.0 * P1y * P3z * P2x) + (1575.0 * P1y * P1z * P4x);
    double term_6 = (2520.0 * P5z * P1y) - (6300.0 * P3z * P3y) - (6300.0 * P3z * P1y * P2x) + (1575.0 * P1z * P5y) + (3150.0 * P1z * P3y * P2x) + (1575.0 * P1z * P1y * P4x);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;

}

void calc_solid_MCSH_6_3_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double P4x = P4(lambda, x0, gamma);
    double P4y = P4(lambda, y0, gamma);
    double P4z = P4(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double term_1 = (5220.0 * P4x * P2y) - (360.0 * P6x) + (180.0 * P4x * P2z) - (4545.0 * P2x * P4y) - (4050.0 * P2x * P2y * P2z) + (495.0 * P2x * P4z) + (270.0 * P6y) + (495.0 * P4y * P2z) + (180.0 * P2y * P4z) - (45.0 * P6z);
    double term_2 = (5220.0 * P4y * P2x) - (360.0 * P6y) + (180.0 * P4y * P2z) - (4545.0 * P2y * P4x) - (4050.0 * P2y * P2x * P2z) + (495.0 * P2y * P4z) + (270.0 * P6x) + (495.0 * P4x * P2z) + (180.0 * P2x * P4z) - (45.0 * P6z);
    double term_3 = (5220.0 * P4x * P2z) - (360.0 * P6x) + (180.0 * P4x * P2y) - (4545.0 * P2x * P4z) - (4050.0 * P2x * P2z * P2y) + (495.0 * P2x * P4y) + (270.0 * P6z) + (495.0 * P4z * P2y) + (180.0 * P2z * P4y) - (45.0 * P6y);
    double term_4 = (5220.0 * P4z * P2x) - (360.0 * P6z) + (180.0 * P4z * P2y) - (4545.0 * P2z * P4x) - (4050.0 * P2z * P2x * P2y) + (495.0 * P2z * P4y) + (270.0 * P6x) + (495.0 * P4x * P2y) + (180.0 * P2x * P4y) - (45.0 * P6y);
    double term_5 = (5220.0 * P4y * P2z) - (360.0 * P6y) + (180.0 * P4y * P2x) - (4545.0 * P2y * P4z) - (4050.0 * P2y * P2z * P2x) + (495.0 * P2y * P4x) + (270.0 * P6z) + (495.0 * P4z * P2x) + (180.0 * P2z * P4x) - (45.0 * P6x);
    double term_6 = (5220.0 * P4z * P2y) - (360.0 * P6z) + (180.0 * P4z * P2x) - (4545.0 * P2z * P4y) - (4050.0 * P2z * P2y * P2x) + (495.0 * P2z * P4x) + (270.0 * P6y) + (495.0 * P4y * P2x) + (180.0 * P2y * P4x) - (45.0 * P6x);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;

}

void calc_solid_MCSH_6_4_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double term_1 = (5040.0 * P4x * P1y * P1z) - (5040.0 * P2x * P3y * P1z) - (5040.0 * P2x * P1y * P3z) + (315.0 * P5y * P1z) + (630.0 * P3y * P3z) + (315.0 * P1y * P5z);
    double term_2 = (5040.0 * P4y * P1x * P1z) - (5040.0 * P2y * P3x * P1z) - (5040.0 * P2y * P1x * P3z) + (315.0 * P5x * P1z) + (630.0 * P3x * P3z) + (315.0 * P1x * P5z);
    double term_3 = (5040.0 * P4z * P1x * P1y) - (5040.0 * P2z * P3x * P1y) - (5040.0 * P2z * P1x * P3y) + (315.0 * P5x * P1y) + (630.0 * P3x * P3y) + (315.0 * P1x * P5y);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;

}

void calc_solid_MCSH_6_5_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double term_1 = (6615.0 * P3x * P3y) - (1890.0 * P5x * P1y) - (945.0 * P3x * P1y * P2z) - (1890.0 * P1x * P5y) - (945.0 * P1x * P3y * P2z) + (945.0 * P1x * P1y * P4z);
    double term_2 = (6615.0 * P3x * P3z) - (1890.0 * P5x * P1z) - (945.0 * P3x * P1z * P2y) - (1890.0 * P1x * P5z) - (945.0 * P1x * P3z * P2y) + (945.0 * P1x * P1z * P4y);
    double term_3 = (6615.0 * P3y * P3z) - (1890.0 * P5y * P1z) - (945.0 * P3y * P1z * P2x) - (1890.0 * P1y * P5z) - (945.0 * P1y * P3z * P2x) + (945.0 * P1y * P1z * P4x);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;

}

void calc_solid_MCSH_6_6_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double term_1 = (7245.0 * P3x * P2y * P1z) - (630.0 * P5x * P1z) - (315.0 * P3x * P3z) - (2520.0 * P1x * P4y * P1z) - (2205.0 * P1x * P2y * P3z) + (315.0 * P1x * P5z);
    double term_2 = (7245.0 * P3y * P2x * P1z) - (630.0 * P5y * P1z) - (315.0 * P3y * P3z) - (2520.0 * P1y * P4x * P1z) - (2205.0 * P1y * P2x * P3z) + (315.0 * P1y * P5z);
    double term_3 = (7245.0 * P3x * P2z * P1y) - (630.0 * P5x * P1y) - (315.0 * P3x * P3y) - (2520.0 * P1x * P4z * P1y) - (2205.0 * P1x * P2z * P3y) + (315.0 * P1x * P5y);
    double term_4 = (7245.0 * P3z * P2x * P1y) - (630.0 * P5z * P1y) - (315.0 * P3z * P3y) - (2520.0 * P1z * P4x * P1y) - (2205.0 * P1z * P2x * P3y) + (315.0 * P1z * P5y);
    double term_5 = (7245.0 * P3y * P2z * P1x) - (630.0 * P5y * P1x) - (315.0 * P3y * P3x) - (2520.0 * P1y * P4z * P1x) - (2205.0 * P1y * P2z * P3x) + (315.0 * P1y * P5x);
    double term_6 = (7245.0 * P3z * P2y * P1x) - (630.0 * P5z * P1x) - (315.0 * P3z * P3x) - (2520.0 * P1z * P4y * P1x) - (2205.0 * P1z * P2y * P3x) + (315.0 * P1z * P5x);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;

}

void calc_solid_MCSH_6_7_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double P4x = P4(lambda, x0, gamma);
    double P4y = P4(lambda, y0, gamma);
    double P4z = P4(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double term = (8100.0 * P2x * P2y * P2z) - (675.0 * P4x * P2y) - (675.0 * P2x * P4y) - (675.0 * P4x * P2z) - (675.0 * P2x * P4z) - (675.0 * P4y * P2z) - (675.0 * P2y * P4z) + (90.0 * P6x) + (90.0 * P6y) + (90.0 * P6z);

    double m = temp * term;

    value[0] = m;

}

void calc_solid_MCSH_7_1_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double term_1 = (5040.0 * P7x) - (52920.0 * P5x * P2y) - (52920.0 * P5x * P2z) + (66150.0 * P3x * P4y) + (132300.0 * P3x * P2y * P2z) + (66150.0 * P3x * P4z) - (11025.0 * P1x * P6y) - (33075.0 * P1x * P4y * P2z) - (33075.0 * P1x * P2y * P4z) - (11025.0 * P1x * P6z);
    double term_2 = (5040.0 * P7y) - (52920.0 * P5y * P2x) - (52920.0 * P5y * P2z) + (66150.0 * P3y * P4x) + (132300.0 * P3y * P2x * P2z) + (66150.0 * P3y * P4z) - (11025.0 * P1y * P6x) - (33075.0 * P1y * P4x * P2z) - (33075.0 * P1y * P2x * P4z) - (11025.0 * P1y * P6z);
    double term_3 = (5040.0 * P7z) - (52920.0 * P5z * P2x) - (52920.0 * P5z * P2y) + (66150.0 * P3z * P4x) + (132300.0 * P3z * P2x * P2y) + (66150.0 * P3z * P4y) - (11025.0 * P1z * P6x) - (33075.0 * P1z * P4x * P2y) - (33075.0 * P1z * P2x * P4y) - (11025.0 * P1z * P6y);
    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
}

void calc_solid_MCSH_7_2_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double term_1 = (20160.0 * P6x * P1y) - (75600.0 * P4x * P3y) - (75600.0 * P4x * P1y * P2z) + (37800.0 * P2x * P5y) + (75600.0 * P2x * P3y * P2z) + (37800.0 * P2x * P1y * P4z) - (1575.0 * P7y) - (4725.0 * P5y * P2z) - (4725.0 * P3y * P4z) - (1575.0 * P1y * P6z);
    double term_2 = (20160.0 * P6y * P1x) - (75600.0 * P4y * P3x) - (75600.0 * P4y * P1x * P2z) + (37800.0 * P2y * P5x) + (75600.0 * P2y * P3x * P2z) + (37800.0 * P2y * P1x * P4z) - (1575.0 * P7x) - (4725.0 * P5x * P2z) - (4725.0 * P3x * P4z) - (1575.0 * P1x * P6z);
    double term_3 = (20160.0 * P6x * P1z) - (75600.0 * P4x * P3z) - (75600.0 * P4x * P1z * P2y) + (37800.0 * P2x * P5z) + (75600.0 * P2x * P3z * P2y) + (37800.0 * P2x * P1z * P4y) - (1575.0 * P7z) - (4725.0 * P5z * P2y) - (4725.0 * P3z * P4y) - (1575.0 * P1z * P6y);
    double term_4 = (20160.0 * P6z * P1x) - (75600.0 * P4z * P3x) - (75600.0 * P4z * P1x * P2y) + (37800.0 * P2z * P5x) + (75600.0 * P2z * P3x * P2y) + (37800.0 * P2z * P1x * P4y) - (1575.0 * P7x) - (4725.0 * P5x * P2y) - (4725.0 * P3x * P4y) - (1575.0 * P1x * P6y);
    double term_5 = (20160.0 * P6y * P1z) - (75600.0 * P4y * P3z) - (75600.0 * P4y * P1z * P2x) + (37800.0 * P2y * P5z) + (75600.0 * P2y * P3z * P2x) + (37800.0 * P2y * P1z * P4x) - (1575.0 * P7z) - (4725.0 * P5z * P2x) - (4725.0 * P3z * P4x) - (1575.0 * P1z * P6x);
    double term_6 = (20160.0 * P6z * P1y) - (75600.0 * P4z * P3y) - (75600.0 * P4z * P1y * P2x) + (37800.0 * P2z * P5y) + (75600.0 * P2z * P3y * P2x) + (37800.0 * P2z * P1y * P4x) - (1575.0 * P7y) - (4725.0 * P5y * P2x) - (4725.0 * P3y * P4x) - (1575.0 * P1y * P6x);
    
    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;
}

void calc_solid_MCSH_7_3_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double term_1 = (49140.0 * P5x * P2y) - (2520.0 * P7x) + (3780.0 * P5x * P2z) - (70875.0 * P3x * P4y) - (66150.0 * P3x * P2y * P2z) + (4725.0 * P3x * P4z) + (12600.0 * P1x * P6y) + (23625.0 * P1x * P4y * P2z) + (9450.0 * P1x * P2y * P4z) - (1575.0 * P1x * P6z);
    double term_2 = (49140.0 * P5y * P2x) - (2520.0 * P7y) + (3780.0 * P5y * P2z) - (70875.0 * P3y * P4x) - (66150.0 * P3y * P2x * P2z) + (4725.0 * P3y * P4z) + (12600.0 * P1y * P6x) + (23625.0 * P1y * P4x * P2z) + (9450.0 * P1y * P2x * P4z) - (1575.0 * P1y * P6z);
    double term_3 = (49140.0 * P5x * P2z) - (2520.0 * P7x) + (3780.0 * P5x * P2y) - (70875.0 * P3x * P4z) - (66150.0 * P3x * P2z * P2y) + (4725.0 * P3x * P4y) + (12600.0 * P1x * P6z) + (23625.0 * P1x * P4z * P2y) + (9450.0 * P1x * P2z * P4y) - (1575.0 * P1x * P6y);
    double term_4 = (49140.0 * P5z * P2x) - (2520.0 * P7z) + (3780.0 * P5z * P2y) - (70875.0 * P3z * P4x) - (66150.0 * P3z * P2x * P2y) + (4725.0 * P3z * P4y) + (12600.0 * P1z * P6x) + (23625.0 * P1z * P4x * P2y) + (9450.0 * P1z * P2x * P4y) - (1575.0 * P1z * P6y);
    double term_5 = (49140.0 * P5y * P2z) - (2520.0 * P7y) + (3780.0 * P5y * P2x) - (70875.0 * P3y * P4z) - (66150.0 * P3y * P2z * P2x) + (4725.0 * P3y * P4x) + (12600.0 * P1y * P6z) + (23625.0 * P1y * P4z * P2x) + (9450.0 * P1y * P2z * P4x) - (1575.0 * P1y * P6x);
    double term_6 = (49140.0 * P5z * P2y) - (2520.0 * P7z) + (3780.0 * P5z * P2x) - (70875.0 * P3z * P4y) - (66150.0 * P3z * P2y * P2x) + (4725.0 * P3z * P4x) + (12600.0 * P1z * P6y) + (23625.0 * P1z * P4y * P2x) + (9450.0 * P1z * P2y * P4x) - (1575.0 * P1z * P6x);
    
    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;
}

void calc_solid_MCSH_7_4_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    // double term_1 = (45360.0 * P5x * P1y * P1z) - (75600.0 * P3x * (P3y * P1z + P1y * P3z)) + (14175.0 * P1x * (P5y * P1z + P1y * P5z)) + (28350.0 * P1x * P3y * P3z);
    // double term_2 = (45360.0 * P5y * P1x * P1z) - (75600.0 * P3y * (P3x * P1z + P1x * P3z)) + (14175.0 * P1y * (P5x * P1z + P1x * P5z)) + (28350.0 * P1y * P3x * P3z);
    // double term_3 = (45360.0 * P5z * P1x * P1y) - (75600.0 * P3z * (P3x * P1y + P1x * P3y)) + (14175.0 * P1z * (P5x * P1y + P1x * P5y)) + (28350.0 * P1z * P3x * P3y);
    double term_1 = (45360.0 * P5x * P1y * P1z) - (75600.0 * P3x * P3y * P1z) - (75600.0 * P3x * P1y * P3z) + (14175.0 * P1x * P5y * P1z) + (28350.0 * P1x * P3y * P3z) + (14175.0 * P1x * P1y * P5z);
    double term_2 = (45360.0 * P5y * P1x * P1z) - (75600.0 * P3y * P3x * P1z) - (75600.0 * P3y * P1x * P3z) + (14175.0 * P1y * P5x * P1z) + (28350.0 * P1y * P3x * P3z) + (14175.0 * P1y * P1x * P5z);
    double term_3 = (45360.0 * P5z * P1x * P1y) - (75600.0 * P3z * P3x * P1y) - (75600.0 * P3z * P1x * P3y) + (14175.0 * P1z * P5x * P1y) + (28350.0 * P1z * P3x * P3y) + (14175.0 * P1z * P1x * P5y);
    
    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
}

void calc_solid_MCSH_7_5_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double term_1 = (75600.0 * P4x * P3y) - (15120.0 * P6x * P1y) - (42525.0 * P2x * P5y) - (28350.0 * P2x * P3y * P2z) + (14175.0 * P2x * P1y * P4z) + (1890.0 * P7y) + (2835.0 * P5y * P2z) - (945.0 * P1y * P6z);
    double term_2 = (75600.0 * P4y * P3x) - (15120.0 * P6y * P1x) - (42525.0 * P2y * P5x) - (28350.0 * P2y * P3x * P2z) + (14175.0 * P2y * P1x * P4z) + (1890.0 * P7x) + (2835.0 * P5x * P2z) - (945.0 * P1x * P6z);
    double term_3 = (75600.0 * P4x * P3z) - (15120.0 * P6x * P1z) - (42525.0 * P2x * P5z) - (28350.0 * P2x * P3z * P2y) + (14175.0 * P2x * P1z * P4y) + (1890.0 * P7z) + (2835.0 * P5z * P2y) - (945.0 * P1z * P6y);
    double term_4 = (75600.0 * P4z * P3x) - (15120.0 * P6z * P1x) - (42525.0 * P2z * P5x) - (28350.0 * P2z * P3x * P2y) + (14175.0 * P2z * P1x * P4y) + (1890.0 * P7x) + (2835.0 * P5x * P2y) - (945.0 * P1x * P6y);
    double term_5 = (75600.0 * P4y * P3z) - (15120.0 * P6y * P1z) - (42525.0 * P2y * P5z) - (28350.0 * P2y * P3z * P2x) + (14175.0 * P2y * P1z * P4x) + (1890.0 * P7z) + (2835.0 * P5z * P2x) - (945.0 * P1z * P6x);
    double term_6 = (75600.0 * P4z * P3y) - (15120.0 * P6z * P1y) - (42525.0 * P2z * P5y) - (28350.0 * P2z * P3y * P2x) + (14175.0 * P2z * P1y * P4x) + (1890.0 * P7y) + (2835.0 * P5y * P2x) - (945.0 * P1y * P6x);
    
    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;
}

void calc_solid_MCSH_7_6_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double term_1 = (75600.0 * P4x * P2y * P1z) - (5040.0 * P6x * P1z) - (51975.0 * P2x * P4y * P1z) - (47250.0 * P2x * P2y * P3z) + (4725.0 * P2x * P5z) + (2520.0 * P6y * P1z) + (4725.0 * P4y * P3z) + (1890.0 * P2y * P5z) - (315.0 * P7z);
    double term_2 = (75600.0 * P4y * P2x * P1z) - (5040.0 * P6y * P1z) - (51975.0 * P2y * P4x * P1z) - (47250.0 * P2y * P2x * P3z) + (4725.0 * P2y * P5z) + (2520.0 * P6x * P1z) + (4725.0 * P4x * P3z) + (1890.0 * P2x * P5z) - (315.0 * P7z);
    double term_3 = (75600.0 * P4x * P2z * P1y) - (5040.0 * P6x * P1y) - (51975.0 * P2x * P4z * P1y) - (47250.0 * P2x * P2z * P3y) + (4725.0 * P2x * P5y) + (2520.0 * P6z * P1y) + (4725.0 * P4z * P3y) + (1890.0 * P2z * P5y) - (315.0 * P7y);
    double term_4 = (75600.0 * P4z * P2x * P1y) - (5040.0 * P6z * P1y) - (51975.0 * P2z * P4x * P1y) - (47250.0 * P2z * P2x * P3y) + (4725.0 * P2z * P5y) + (2520.0 * P6x * P1y) + (4725.0 * P4x * P3y) + (1890.0 * P2x * P5y) - (315.0 * P7y);
    double term_5 = (75600.0 * P4y * P2z * P1x) - (5040.0 * P6y * P1x) - (51975.0 * P2y * P4z * P1x) - (47250.0 * P2y * P2z * P3x) + (4725.0 * P2y * P5x) + (2520.0 * P6z * P1x) + (4725.0 * P4z * P3x) + (1890.0 * P2z * P5x) - (315.0 * P7x);
    double term_6 = (75600.0 * P4z * P2y * P1x) - (5040.0 * P6z * P1x) - (51975.0 * P2z * P4y * P1x) - (47250.0 * P2z * P2y * P3x) + (4725.0 * P2z * P5x) + (2520.0 * P6y * P1x) + (4725.0 * P4y * P3x) + (1890.0 * P2y * P5x) - (315.0 * P7x);  
    
    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;
}


void calc_solid_MCSH_7_7_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    // double term_1 = (89775.0 * P3x * P3y * P1z) - (22680.0 * P1z * (P5x * P1y + P1x * P5y)) - (14175.0 * P3z * (P3x * P1y + P1x * P3y)) + (8505.0 * P1x * P1y * P5z);
    // double term_2 = (89775.0 * P3y * P3x * P1z) - (22680.0 * P1z * (P5y * P1x + P1y * P5x)) - (14175.0 * P3z * (P3y * P1x + P1y * P3x)) + (8505.0 * P1y * P1x * P5z);
    // double term_3 = (89775.0 * P3z * P3x * P1y) - (22680.0 * P1y * (P5z * P1x + P1z * P5x)) - (14175.0 * P3y * (P3z * P1x + P1z * P3x)) + (8505.0 * P1z * P1x * P5y);
    double term_1 = (89775.0 * P3x * P3y * P1z) - (22680.0 * P5x * P1y * P1z) - (14175.0 * P3x * P1y * P3z) - (22680.0 * P1x * P5y * P1z) - (14175.0 * P1x * P3y * P3z) + (8505.0 * P1x * P1y * P5z);
    double term_2 = (89775.0 * P3x * P3z * P1y) - (22680.0 * P5x * P1z * P1y) - (14175.0 * P3x * P1z * P3y) - (22680.0 * P1x * P5z * P1y) - (14175.0 * P1x * P3z * P3y) + (8505.0 * P1x * P1z * P5y);
    double term_3 = (89775.0 * P3y * P3z * P1x) - (22680.0 * P5y * P1z * P1x) - (14175.0 * P3y * P1z * P3x) - (22680.0 * P1y * P5z * P1x) - (14175.0 * P1y * P3z * P3x) + (8505.0 * P1y * P1z * P5x);
    
    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
}

void calc_solid_MCSH_7_8_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    // double term_1 = (94500.0 * P3x * P2y * P2z) - (6615.0 * P5x * (P2y + P2z)) - (4725.0 * P3x * (P4y + P4z)) + (630.0 * P7x) - (23625.0 * P1x * (P4y * P2z + P2y * P4z)) + (2520.0 * P1x * (P6y + P6z));
    // double term_2 = (94500.0 * P3y * P2x * P2z) - (6615.0 * P5y * (P2x + P2z)) - (4725.0 * P3y * (P4x + P4z)) + (630.0 * P7y) - (23625.0 * P1y * (P4x * P2z + P2x * P4z)) + (2520.0 * P1y * (P6x + P6z));
    // double term_3 = (94500.0 * P3z * P2x * P2y) - (6615.0 * P5z * (P2x + P2y)) - (4725.0 * P3z * (P4x + P4y)) + (630.0 * P7z) - (23625.0 * P1z * (P4x * P2y + P2x * P4y)) + (2520.0 * P1z * (P6x + P6y));   
    double term_1 = (94500.0 * P3x * P2y * P2z) - (6615.0 * P5x * P2z) - (4725.0 * P3x * P4z) - (6615.0 * P5x * P2y) - (4725.0 * P3x * P4y) + (630.0 * P7x) - (23625.0 * P1x * P4y * P2z) - (23625.0 * P1x * P2y * P4z) + (2520.0 * P1x * P6z) + (2520.0 * P1x * P6y);
    double term_2 = (94500.0 * P3y * P2x * P2z) - (6615.0 * P5y * P2z) - (4725.0 * P3y * P4z) - (6615.0 * P5y * P2x) - (4725.0 * P3y * P4x) + (630.0 * P7y) - (23625.0 * P1y * P4x * P2z) - (23625.0 * P1y * P2x * P4z) + (2520.0 * P1y * P6z) + (2520.0 * P1y * P6x);
    double term_3 = (94500.0 * P3z * P2x * P2y) - (6615.0 * P5z * P2y) - (4725.0 * P3z * P4y) - (6615.0 * P5z * P2x) - (4725.0 * P3z * P4x) + (630.0 * P7z) - (23625.0 * P1z * P4x * P2y) - (23625.0 * P1z * P2x * P4y) + (2520.0 * P1z * P6y) + (2520.0 * P1z * P6x);
    
    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
}


void calc_solid_MCSH_8_1_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double term_1 = (40320.0 * P8x) - (564480.0 * P6x * P2y) - (564480.0 * P6x * P2z) + (1058400.0 * P4x * P4y) + (2116800.0 * P4x * P2y * P2z) + (1058400.0 * P4x * P4z) - (352800.0 * P2x * P6y) - (1058400.0 * P2x * P4y * P2z) - (1058400.0 * P2x * P2y * P4z) - (352800.0 * P2x * P6z) + (11025.0 * P8y) + (44100.0 * P6y * P2z) + (66150.0 * P4y * P4z) + (44100.0 * P2y * P6z) + (11025.0 * P8z);
    double term_2 = (40320.0 * P8y) - (564480.0 * P6y * P2x) - (564480.0 * P6y * P2z) + (1058400.0 * P4y * P4x) + (2116800.0 * P4y * P2x * P2z) + (1058400.0 * P4y * P4z) - (352800.0 * P2y * P6x) - (1058400.0 * P2y * P4x * P2z) - (1058400.0 * P2y * P2x * P4z) - (352800.0 * P2y * P6z) + (11025.0 * P8x) + (44100.0 * P6x * P2z) + (66150.0 * P4x * P4z) + (44100.0 * P2x * P6z) + (11025.0 * P8z);
    double term_3 = (40320.0 * P8z) - (564480.0 * P6z * P2x) - (564480.0 * P6z * P2y) + (1058400.0 * P4z * P4x) + (2116800.0 * P4z * P2x * P2y) + (1058400.0 * P4z * P4y) - (352800.0 * P2z * P6x) - (1058400.0 * P2z * P4x * P2y) - (1058400.0 * P2z * P2x * P4y) - (352800.0 * P2z * P6y) + (11025.0 * P8x) + (44100.0 * P6x * P2y) + (66150.0 * P4x * P4y) + (44100.0 * P2x * P6y) + (11025.0 * P8y);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;

}

void calc_solid_MCSH_8_2_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double term_1 = (181440.0 * P7x * P1y) - (952560.0 * P5x * P3y) - (952560.0 * P5x * P1y * P2z) + (793800.0 * P3x * P5y) + (1587600.0 * P3x * P3y * P2z) + (793800.0 * P3x * P1y * P4z) - (99225.0 * P1x * P7y) - (297675.0 * P1x * P5y * P2z) - (297675.0 * P1x * P3y * P4z) - (99225.0 * P1x * P1y * P6z);
    double term_2 = (181440.0 * P7y * P1x) - (952560.0 * P5y * P3x) - (952560.0 * P5y * P1x * P2z) + (793800.0 * P3y * P5x) + (1587600.0 * P3y * P3x * P2z) + (793800.0 * P3y * P1x * P4z) - (99225.0 * P1y * P7x) - (297675.0 * P1y * P5x * P2z) - (297675.0 * P1y * P3x * P4z) - (99225.0 * P1y * P1x * P6z);
    double term_3 = (181440.0 * P7x * P1z) - (952560.0 * P5x * P3z) - (952560.0 * P5x * P1z * P2y) + (793800.0 * P3x * P5z) + (1587600.0 * P3x * P3z * P2y) + (793800.0 * P3x * P1z * P4y) - (99225.0 * P1x * P7z) - (297675.0 * P1x * P5z * P2y) - (297675.0 * P1x * P3z * P4y) - (99225.0 * P1x * P1z * P6y);
    double term_4 = (181440.0 * P7z * P1x) - (952560.0 * P5z * P3x) - (952560.0 * P5z * P1x * P2y) + (793800.0 * P3z * P5x) + (1587600.0 * P3z * P3x * P2y) + (793800.0 * P3z * P1x * P4y) - (99225.0 * P1z * P7x) - (297675.0 * P1z * P5x * P2y) - (297675.0 * P1z * P3x * P4y) - (99225.0 * P1z * P1x * P6y);
    double term_5 = (181440.0 * P7y * P1z) - (952560.0 * P5y * P3z) - (952560.0 * P5y * P1z * P2x) + (793800.0 * P3y * P5z) + (1587600.0 * P3y * P3z * P2x) + (793800.0 * P3y * P1z * P4x) - (99225.0 * P1y * P7z) - (297675.0 * P1y * P5z * P2x) - (297675.0 * P1y * P3z * P4x) - (99225.0 * P1y * P1z * P6x);
    double term_6 = (181440.0 * P7z * P1y) - (952560.0 * P5z * P3y) - (952560.0 * P5z * P1y * P2x) + (793800.0 * P3z * P5y) + (1587600.0 * P3z * P3y * P2x) + (793800.0 * P3z * P1y * P4x) - (99225.0 * P1z * P7y) - (297675.0 * P1z * P5y * P2x) - (297675.0 * P1z * P3y * P4x) - (99225.0 * P1z * P1y * P6x);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;

}

void calc_solid_MCSH_8_3_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double term_1 = (509040.0 * P6x * P2y) - (20160.0 * P8x) + (55440.0 * P6x * P2z) - (1096200.0 * P4x * P4y) - (1058400.0 * P4x * P2y * P2z) + (37800.0 * P4x * P4z) + (389025.0 * P2x * P6y) + (741825.0 * P2x * P4y * P2z) + (316575.0 * P2x * P2y * P4z) - (36225.0 * P2x * P6z) - (12600.0 * P8y) - (36225.0 * P6y * P2z) - (33075.0 * P4y * P4z) - (7875.0 * P2y * P6z) + (1575.0 * P8z);
    double term_2 = (509040.0 * P6y * P2x) - (20160.0 * P8y) + (55440.0 * P6y * P2z) - (1096200.0 * P4y * P4x) - (1058400.0 * P4y * P2x * P2z) + (37800.0 * P4y * P4z) + (389025.0 * P2y * P6x) + (741825.0 * P2y * P4x * P2z) + (316575.0 * P2y * P2x * P4z) - (36225.0 * P2y * P6z) - (12600.0 * P8x) - (36225.0 * P6x * P2z) - (33075.0 * P4x * P4z) - (7875.0 * P2x * P6z) + (1575.0 * P8z);
    double term_3 = (509040.0 * P6x * P2z) - (20160.0 * P8x) + (55440.0 * P6x * P2y) - (1096200.0 * P4x * P4z) - (1058400.0 * P4x * P2z * P2y) + (37800.0 * P4x * P4y) + (389025.0 * P2x * P6z) + (741825.0 * P2x * P4z * P2y) + (316575.0 * P2x * P2z * P4y) - (36225.0 * P2x * P6y) - (12600.0 * P8z) - (36225.0 * P6z * P2y) - (33075.0 * P4z * P4y) - (7875.0 * P2z * P6y) + (1575.0 * P8y);
    double term_4 = (509040.0 * P6z * P2x) - (20160.0 * P8z) + (55440.0 * P6z * P2y) - (1096200.0 * P4z * P4x) - (1058400.0 * P4z * P2x * P2y) + (37800.0 * P4z * P4y) + (389025.0 * P2z * P6x) + (741825.0 * P2z * P4x * P2y) + (316575.0 * P2z * P2x * P4y) - (36225.0 * P2z * P6y) - (12600.0 * P8x) - (36225.0 * P6x * P2y) - (33075.0 * P4x * P4y) - (7875.0 * P2x * P6y) + (1575.0 * P8y);
    double term_5 = (509040.0 * P6y * P2z) - (20160.0 * P8y) + (55440.0 * P6y * P2x) - (1096200.0 * P4y * P4z) - (1058400.0 * P4y * P2z * P2x) + (37800.0 * P4y * P4x) + (389025.0 * P2y * P6z) + (741825.0 * P2y * P4z * P2x) + (316575.0 * P2y * P2z * P4x) - (36225.0 * P2y * P6x) - (12600.0 * P8z) - (36225.0 * P6z * P2x) - (33075.0 * P4z * P4x) - (7875.0 * P2z * P6x) + (1575.0 * P8x);
    double term_6 = (509040.0 * P6z * P2y) - (20160.0 * P8z) + (55440.0 * P6z * P2x) - (1096200.0 * P4z * P4y) - (1058400.0 * P4z * P2y * P2x) + (37800.0 * P4z * P4x) + (389025.0 * P2z * P6y) + (741825.0 * P2z * P4y * P2x) + (316575.0 * P2z * P2y * P4x) - (36225.0 * P2z * P6x) - (12600.0 * P8y) - (36225.0 * P6y * P2x) - (33075.0 * P4y * P4x) - (7875.0 * P2y * P6x) + (1575.0 * P8x);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;

}

void calc_solid_MCSH_8_4_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double term_1 = (453600.0 * P6x * P1y * P1z) - (1134000.0 * P4x * P3y * P1z) - (1134000.0 * P4x * P1y * P3z) + (425250.0 * P2x * P5y * P1z) + (850500.0 * P2x * P3y * P3z) + (425250.0 * P2x * P1y * P5z) - (14175.0 * P7y * P1z) - (42525.0 * P5y * P3z) - (42525.0 * P3y * P5z) - (14175.0 * P1y * P7z);
    double term_2 = (453600.0 * P6y * P1x * P1z) - (1134000.0 * P4y * P3x * P1z) - (1134000.0 * P4y * P1x * P3z) + (425250.0 * P2y * P5x * P1z) + (850500.0 * P2y * P3x * P3z) + (425250.0 * P2y * P1x * P5z) - (14175.0 * P7x * P1z) - (42525.0 * P5x * P3z) - (42525.0 * P3x * P5z) - (14175.0 * P1x * P7z);
    double term_3 = (453600.0 * P6z * P1x * P1y) - (1134000.0 * P4z * P3x * P1y) - (1134000.0 * P4z * P1x * P3y) + (425250.0 * P2z * P5x * P1y) + (850500.0 * P2z * P3x * P3y) + (425250.0 * P2z * P1x * P5y) - (14175.0 * P7x * P1y) - (42525.0 * P5x * P3y) - (42525.0 * P3x * P5y) - (14175.0 * P1x * P7y);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;

}

void calc_solid_MCSH_8_5_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double term_1 = (922320.0 * P5x * P3y) - (136080.0 * P7x * P1y) + (90720.0 * P5x * P1y * P2z) - (855225.0 * P3x * P5y) - (670950.0 * P3x * P3y * P2z) + (184275.0 * P3x * P1y * P4z) + (113400.0 * P1x * P7y) + (184275.0 * P1x * P5y * P2z) + (28350.0 * P1x * P3y * P4z) - (42525.0 * P1x * P1y * P6z);
    double term_2 = (922320.0 * P5y * P3x) - (136080.0 * P7y * P1x) + (90720.0 * P5y * P1x * P2z) - (855225.0 * P3y * P5x) - (670950.0 * P3y * P3x * P2z) + (184275.0 * P3y * P1x * P4z) + (113400.0 * P1y * P7x) + (184275.0 * P1y * P5x * P2z) + (28350.0 * P1y * P3x * P4z) - (42525.0 * P1y * P1x * P6z);
    double term_3 = (922320.0 * P5x * P3z) - (136080.0 * P7x * P1z) + (90720.0 * P5x * P1z * P2y) - (855225.0 * P3x * P5z) - (670950.0 * P3x * P3z * P2y) + (184275.0 * P3x * P1z * P4y) + (113400.0 * P1x * P7z) + (184275.0 * P1x * P5z * P2y) + (28350.0 * P1x * P3z * P4y) - (42525.0 * P1x * P1z * P6y);
    double term_4 = (922320.0 * P5z * P3x) - (136080.0 * P7z * P1x) + (90720.0 * P5z * P1x * P2y) - (855225.0 * P3z * P5x) - (670950.0 * P3z * P3x * P2y) + (184275.0 * P3z * P1x * P4y) + (113400.0 * P1z * P7x) + (184275.0 * P1z * P5x * P2y) + (28350.0 * P1z * P3x * P4y) - (42525.0 * P1z * P1x * P6y);
    double term_5 = (922320.0 * P5y * P3z) - (136080.0 * P7y * P1z) + (90720.0 * P5y * P1z * P2x) - (855225.0 * P3y * P5z) - (670950.0 * P3y * P3z * P2x) + (184275.0 * P3y * P1z * P4x) + (113400.0 * P1y * P7z) + (184275.0 * P1y * P5z * P2x) + (28350.0 * P1y * P3z * P4x) - (42525.0 * P1y * P1z * P6x);
    double term_6 = (922320.0 * P5z * P3y) - (136080.0 * P7z * P1y) + (90720.0 * P5z * P1y * P2x) - (855225.0 * P3z * P5y) - (670950.0 * P3z * P3y * P2x) + (184275.0 * P3z * P1y * P4x) + (113400.0 * P1z * P7y) + (184275.0 * P1z * P5y * P2x) + (28350.0 * P1z * P3y * P4x) - (42525.0 * P1z * P1y * P6x);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;

}

void calc_solid_MCSH_8_6_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double term_1 = (861840.0 * P5x * P2y * P1z) - (45360.0 * P7x * P1z) + (30240.0 * P5x * P3z) - (978075.0 * P3x * P4y * P1z) - (916650.0 * P3x * P2y * P3z) + (61425.0 * P3x * P5z) + (141750.0 * P1x * P6y * P1z) + (269325.0 * P1x * P4y * P3z) + (113400.0 * P1x * P2y * P5z) - (14175.0 * P1x * P7z);
    double term_2 = (861840.0 * P5y * P2x * P1z) - (45360.0 * P7y * P1z) + (30240.0 * P5y * P3z) - (978075.0 * P3y * P4x * P1z) - (916650.0 * P3y * P2x * P3z) + (61425.0 * P3y * P5z) + (141750.0 * P1y * P6x * P1z) + (269325.0 * P1y * P4x * P3z) + (113400.0 * P1y * P2x * P5z) - (14175.0 * P1y * P7z);
    double term_3 = (861840.0 * P5x * P2z * P1y) - (45360.0 * P7x * P1y) + (30240.0 * P5x * P3y) - (978075.0 * P3x * P4z * P1y) - (916650.0 * P3x * P2z * P3y) + (61425.0 * P3x * P5y) + (141750.0 * P1x * P6z * P1y) + (269325.0 * P1x * P4z * P3y) + (113400.0 * P1x * P2z * P5y) - (14175.0 * P1x * P7y);
    double term_4 = (861840.0 * P5z * P2x * P1y) - (45360.0 * P7z * P1y) + (30240.0 * P5z * P3y) - (978075.0 * P3z * P4x * P1y) - (916650.0 * P3z * P2x * P3y) + (61425.0 * P3z * P5y) + (141750.0 * P1z * P6x * P1y) + (269325.0 * P1z * P4x * P3y) + (113400.0 * P1z * P2x * P5y) - (14175.0 * P1z * P7y);
    double term_5 = (861840.0 * P5y * P2z * P1x) - (45360.0 * P7y * P1x) + (30240.0 * P5y * P3x) - (978075.0 * P3y * P4z * P1x) - (916650.0 * P3y * P2z * P3x) + (61425.0 * P3y * P5x) + (141750.0 * P1y * P6z * P1x) + (269325.0 * P1y * P4z * P3x) + (113400.0 * P1y * P2z * P5x) - (14175.0 * P1y * P7x);
    double term_6 = (861840.0 * P5z * P2y * P1x) - (45360.0 * P7z * P1x) + (30240.0 * P5z * P3x) - (978075.0 * P3z * P4y * P1x) - (916650.0 * P3z * P2y * P3x) + (61425.0 * P3z * P5x) + (141750.0 * P1z * P6y * P1x) + (269325.0 * P1z * P4y * P3x) + (113400.0 * P1z * P2y * P5x) - (14175.0 * P1z * P7x);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;

}

void calc_solid_MCSH_8_7_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double term_1 = (1119825.0 * P4x * P4y) - (438480.0 * P6x * P2y) - (141750.0 * P4x * P2y * P2z) + (15120.0 * P8x) + (15120.0 * P6x * P2z) - (14175.0 * P4x * P4z) - (438480.0 * P2x * P6y) - (141750.0 * P2x * P4y * P2z) + (283500.0 * P2x * P2y * P4z) - (13230.0 * P2x * P6z) + (15120.0 * P8y) + (15120.0 * P6y * P2z) - (14175.0 * P4y * P4z) - (13230.0 * P2y * P6z) + (945.0 * P8z);
    double term_2 = (1119825.0 * P4x * P4z) - (438480.0 * P6x * P2z) - (141750.0 * P4x * P2z * P2y) + (15120.0 * P8x) + (15120.0 * P6x * P2y) - (14175.0 * P4x * P4y) - (438480.0 * P2x * P6z) - (141750.0 * P2x * P4z * P2y) + (283500.0 * P2x * P2z * P4y) - (13230.0 * P2x * P6y) + (15120.0 * P8z) + (15120.0 * P6z * P2y) - (14175.0 * P4z * P4y) - (13230.0 * P2z * P6y) + (945.0 * P8y);
    double term_3 = (1119825.0 * P4y * P4z) - (438480.0 * P6y * P2z) - (141750.0 * P4y * P2z * P2x) + (15120.0 * P8y) + (15120.0 * P6y * P2x) - (14175.0 * P4y * P4x) - (438480.0 * P2y * P6z) - (141750.0 * P2y * P4z * P2x) + (283500.0 * P2y * P2z * P4x) - (13230.0 * P2y * P6x) + (15120.0 * P8z) + (15120.0 * P6z * P2x) - (14175.0 * P4z * P4x) - (13230.0 * P2z * P6x) + (945.0 * P8x);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;

}

void calc_solid_MCSH_8_8_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double term_1 = (1190700.0 * P4x * P3y * P1z) - (226800.0 * P6x * P1y * P1z) - (56700.0 * P4x * P1y * P3z) - (586845.0 * P2x * P5y * P1z) - (425250.0 * P2x * P3y * P3z) + (161595.0 * P2x * P1y * P5z) + (22680.0 * P7y * P1z) + (36855.0 * P5y * P3z) + (5670.0 * P3y * P5z) - (8505.0 * P1y * P7z);
    double term_2 = (1190700.0 * P4y * P3x * P1z) - (226800.0 * P6y * P1x * P1z) - (56700.0 * P4y * P1x * P3z) - (586845.0 * P2y * P5x * P1z) - (425250.0 * P2y * P3x * P3z) + (161595.0 * P2y * P1x * P5z) + (22680.0 * P7x * P1z) + (36855.0 * P5x * P3z) + (5670.0 * P3x * P5z) - (8505.0 * P1x * P7z);
    double term_3 = (1190700.0 * P4x * P3z * P1y) - (226800.0 * P6x * P1z * P1y) - (56700.0 * P4x * P1z * P3y) - (586845.0 * P2x * P5z * P1y) - (425250.0 * P2x * P3z * P3y) + (161595.0 * P2x * P1z * P5y) + (22680.0 * P7z * P1y) + (36855.0 * P5z * P3y) + (5670.0 * P3z * P5y) - (8505.0 * P1z * P7y);
    double term_4 = (1190700.0 * P4z * P3x * P1y) - (226800.0 * P6z * P1x * P1y) - (56700.0 * P4z * P1x * P3y) - (586845.0 * P2z * P5x * P1y) - (425250.0 * P2z * P3x * P3y) + (161595.0 * P2z * P1x * P5y) + (22680.0 * P7x * P1y) + (36855.0 * P5x * P3y) + (5670.0 * P3x * P5y) - (8505.0 * P1x * P7y);
    double term_5 = (1190700.0 * P4y * P3z * P1x) - (226800.0 * P6y * P1z * P1x) - (56700.0 * P4y * P1z * P3x) - (586845.0 * P2y * P5z * P1x) - (425250.0 * P2y * P3z * P3x) + (161595.0 * P2y * P1z * P5x) + (22680.0 * P7z * P1x) + (36855.0 * P5z * P3x) + (5670.0 * P3z * P5x) - (8505.0 * P1z * P7x);
    double term_6 = (1190700.0 * P4z * P3y * P1x) - (226800.0 * P6z * P1y * P1x) - (56700.0 * P4z * P1y * P3x) - (586845.0 * P2z * P5y * P1x) - (425250.0 * P2z * P3y * P3x) + (161595.0 * P2z * P1y * P5x) + (22680.0 * P7y * P1x) + (36855.0 * P5y * P3x) + (5670.0 * P3y * P5x) - (8505.0 * P1y * P7x);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;

}

void calc_solid_MCSH_8_9_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double term_1 = (1200150.0 * P4x * P2y * P2z) - (70560.0 * P6x * P2y) - (23625.0 * P4x * P4y) - (70560.0 * P6x * P2z) - (23625.0 * P4x * P4z) + (5040.0 * P8x) - (600075.0 * P2x * P4y * P2z) - (600075.0 * P2x * P2y * P4z) + (49455.0 * P2x * P6y) + (49455.0 * P2x * P6z) + (21105.0 * P6y * P2z) + (47250.0 * P4y * P4z) + (21105.0 * P2y * P6z) - (2520.0 * P8y) - (2520.0 * P8z);
    double term_2 = (1200150.0 * P4y * P2x * P2z) - (70560.0 * P6y * P2x) - (23625.0 * P4y * P4x) - (70560.0 * P6y * P2z) - (23625.0 * P4y * P4z) + (5040.0 * P8y) - (600075.0 * P2y * P4x * P2z) - (600075.0 * P2y * P2x * P4z) + (49455.0 * P2y * P6x) + (49455.0 * P2y * P6z) + (21105.0 * P6x * P2z) + (47250.0 * P4x * P4z) + (21105.0 * P2x * P6z) - (2520.0 * P8x) - (2520.0 * P8z);
    double term_3 = (1200150.0 * P4z * P2x * P2y) - (70560.0 * P6z * P2x) - (23625.0 * P4z * P4x) - (70560.0 * P6z * P2y) - (23625.0 * P4z * P4y) + (5040.0 * P8z) - (600075.0 * P2z * P4x * P2y) - (600075.0 * P2z * P2x * P4y) + (49455.0 * P2z * P6x) + (49455.0 * P2z * P6y) + (21105.0 * P6x * P2y) + (47250.0 * P4x * P4y) + (21105.0 * P2x * P6y) - (2520.0 * P8x) - (2520.0 * P8y);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;

}

void calc_solid_MCSH_8_10_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double term_1 = (1341900.0 * P3x * P3y * P2z) - (67095.0 * P5x * P3y) - (67095.0 * P3x * P5y) - (274995.0 * P5x * P1y * P2z) - (212625.0 * P3x * P1y * P4z) + (22680.0 * P7x * P1y) - (274995.0 * P1x * P5y * P2z) - (212625.0 * P1x * P3y * P4z) + (22680.0 * P1x * P7y) + (85050.0 * P1x * P1y * P6z);
    double term_2 = (1341900.0 * P3x * P3z * P2y) - (67095.0 * P5x * P3z) - (67095.0 * P3x * P5z) - (274995.0 * P5x * P1z * P2y) - (212625.0 * P3x * P1z * P4y) + (22680.0 * P7x * P1z) - (274995.0 * P1x * P5z * P2y) - (212625.0 * P1x * P3z * P4y) + (22680.0 * P1x * P7z) + (85050.0 * P1x * P1z * P6y);
    double term_3 = (1341900.0 * P3y * P3z * P2x) - (67095.0 * P5y * P3z) - (67095.0 * P3y * P5z) - (274995.0 * P5y * P1z * P2x) - (212625.0 * P3y * P1z * P4x) + (22680.0 * P7y * P1z) - (274995.0 * P1y * P5z * P2x) - (212625.0 * P1y * P3z * P4x) + (22680.0 * P1y * P7z) + (85050.0 * P1y * P1z * P6x);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;

}


void calc_solid_MCSH_9_1_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double P8x = P8(lambda, x0, gamma);
    double P8y = P8(lambda, y0, gamma);
    double P8z = P8(lambda, z0, gamma);

    double P9x = P9(lambda, x0, gamma);
    double P9y = P9(lambda, y0, gamma);
    double P9z = P9(lambda, z0, gamma);

    double term_1 = (362880.0 * P9x) - (6531840.0 * P7x * P2y) - (6531840.0 * P7x * P2z) + (17146080.0 * P5x * P4y) + (34292160.0 * P5x * P2y * P2z) + (17146080.0 * P5x * P4z) - (9525600.0 * P3x * P6y) - (28576800.0 * P3x * P4y * P2z) - (28576800.0 * P3x * P2y * P4z) - (9525600.0 * P3x * P6z) + (893025.0 * P1x * P8y) + (3572100.0 * P1x * P6y * P2z) + (5358150.0 * P1x * P4y * P4z) + (3572100.0 * P1x * P2y * P6z) + (893025.0 * P1x * P8z);
    double term_2 = (362880.0 * P9y) - (6531840.0 * P7y * P2x) - (6531840.0 * P7y * P2z) + (17146080.0 * P5y * P4x) + (34292160.0 * P5y * P2x * P2z) + (17146080.0 * P5y * P4z) - (9525600.0 * P3y * P6x) - (28576800.0 * P3y * P4x * P2z) - (28576800.0 * P3y * P2x * P4z) - (9525600.0 * P3y * P6z) + (893025.0 * P1y * P8x) + (3572100.0 * P1y * P6x * P2z) + (5358150.0 * P1y * P4x * P4z) + (3572100.0 * P1y * P2x * P6z) + (893025.0 * P1y * P8z);
    double term_3 = (362880.0 * P9z) - (6531840.0 * P7z * P2x) - (6531840.0 * P7z * P2y) + (17146080.0 * P5z * P4x) + (34292160.0 * P5z * P2x * P2y) + (17146080.0 * P5z * P4y) - (9525600.0 * P3z * P6x) - (28576800.0 * P3z * P4x * P2y) - (28576800.0 * P3z * P2x * P4y) - (9525600.0 * P3z * P6y) + (893025.0 * P1z * P8x) + (3572100.0 * P1z * P6x * P2y) + (5358150.0 * P1z * P4x * P4y) + (3572100.0 * P1z * P2x * P6y) + (893025.0 * P1z * P8y);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;

}

void calc_solid_MCSH_9_2_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double P8x = P8(lambda, x0, gamma);
    double P8y = P8(lambda, y0, gamma);
    double P8z = P8(lambda, z0, gamma);

    double P9x = P9(lambda, x0, gamma);
    double P9y = P9(lambda, y0, gamma);
    double P9z = P9(lambda, z0, gamma);

    double term_1 = (1814400.0 * P8x * P1y) - (12700800.0 * P6x * P3y) - (12700800.0 * P6x * P1y * P2z) + (15876000.0 * P4x * P5y) + (31752000.0 * P4x * P3y * P2z) + (15876000.0 * P4x * P1y * P4z) - (3969000.0 * P2x * P7y) - (11907000.0 * P2x * P5y * P2z) - (11907000.0 * P2x * P3y * P4z) - (3969000.0 * P2x * P1y * P6z) + (99225.0 * P9y) + (396900.0 * P7y * P2z) + (595350.0 * P5y * P4z) + (396900.0 * P3y * P6z) + (99225.0 * P1y * P8z);
    double term_2 = (1814400.0 * P8y * P1x) - (12700800.0 * P6y * P3x) - (12700800.0 * P6y * P1x * P2z) + (15876000.0 * P4y * P5x) + (31752000.0 * P4y * P3x * P2z) + (15876000.0 * P4y * P1x * P4z) - (3969000.0 * P2y * P7x) - (11907000.0 * P2y * P5x * P2z) - (11907000.0 * P2y * P3x * P4z) - (3969000.0 * P2y * P1x * P6z) + (99225.0 * P9x) + (396900.0 * P7x * P2z) + (595350.0 * P5x * P4z) + (396900.0 * P3x * P6z) + (99225.0 * P1x * P8z);
    double term_3 = (1814400.0 * P8x * P1z) - (12700800.0 * P6x * P3z) - (12700800.0 * P6x * P1z * P2y) + (15876000.0 * P4x * P5z) + (31752000.0 * P4x * P3z * P2y) + (15876000.0 * P4x * P1z * P4y) - (3969000.0 * P2x * P7z) - (11907000.0 * P2x * P5z * P2y) - (11907000.0 * P2x * P3z * P4y) - (3969000.0 * P2x * P1z * P6y) + (99225.0 * P9z) + (396900.0 * P7z * P2y) + (595350.0 * P5z * P4y) + (396900.0 * P3z * P6y) + (99225.0 * P1z * P8y);
    double term_4 = (1814400.0 * P8z * P1x) - (12700800.0 * P6z * P3x) - (12700800.0 * P6z * P1x * P2y) + (15876000.0 * P4z * P5x) + (31752000.0 * P4z * P3x * P2y) + (15876000.0 * P4z * P1x * P4y) - (3969000.0 * P2z * P7x) - (11907000.0 * P2z * P5x * P2y) - (11907000.0 * P2z * P3x * P4y) - (3969000.0 * P2z * P1x * P6y) + (99225.0 * P9x) + (396900.0 * P7x * P2y) + (595350.0 * P5x * P4y) + (396900.0 * P3x * P6y) + (99225.0 * P1x * P8y);
    double term_5 = (1814400.0 * P8y * P1z) - (12700800.0 * P6y * P3z) - (12700800.0 * P6y * P1z * P2x) + (15876000.0 * P4y * P5z) + (31752000.0 * P4y * P3z * P2x) + (15876000.0 * P4y * P1z * P4x) - (3969000.0 * P2y * P7z) - (11907000.0 * P2y * P5z * P2x) - (11907000.0 * P2y * P3z * P4x) - (3969000.0 * P2y * P1z * P6x) + (99225.0 * P9z) + (396900.0 * P7z * P2x) + (595350.0 * P5z * P4x) + (396900.0 * P3z * P6x) + (99225.0 * P1z * P8x);
    double term_6 = (1814400.0 * P8z * P1y) - (12700800.0 * P6z * P3y) - (12700800.0 * P6z * P1y * P2x) + (15876000.0 * P4z * P5y) + (31752000.0 * P4z * P3y * P2x) + (15876000.0 * P4z * P1y * P4x) - (3969000.0 * P2z * P7y) - (11907000.0 * P2z * P5y * P2x) - (11907000.0 * P2z * P3y * P4x) - (3969000.0 * P2z * P1y * P6x) + (99225.0 * P9y) + (396900.0 * P7y * P2x) + (595350.0 * P5y * P4x) + (396900.0 * P3y * P6x) + (99225.0 * P1y * P8x);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;

}

void calc_solid_MCSH_9_3_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double P8x = P8(lambda, x0, gamma);
    double P8y = P8(lambda, y0, gamma);
    double P8z = P8(lambda, z0, gamma);

    double P9x = P9(lambda, x0, gamma);
    double P9y = P9(lambda, y0, gamma);
    double P9z = P9(lambda, z0, gamma);

    double term_1 = (5760720.0 * P7x * P2y) - (181440.0 * P9x) + (771120.0 * P7x * P2z) - (17304840.0 * P5x * P4y) - (17146080.0 * P5x * P2y * P2z) + (158760.0 * P5x * P4z) + (10220175.0 * P3x * P6y) + (19745775.0 * P3x * P4y * P2z) + (8831025.0 * P3x * P2y * P4z) - (694575.0 * P3x * P6z) - (992250.0 * P1x * P8y) - (2877525.0 * P1x * P6y * P2z) - (2679075.0 * P1x * P4y * P4z) - (694575.0 * P1x * P2y * P6z) + (99225.0 * P1x * P8z);
    double term_2 = (5760720.0 * P7y * P2x) - (181440.0 * P9y) + (771120.0 * P7y * P2z) - (17304840.0 * P5y * P4x) - (17146080.0 * P5y * P2x * P2z) + (158760.0 * P5y * P4z) + (10220175.0 * P3y * P6x) + (19745775.0 * P3y * P4x * P2z) + (8831025.0 * P3y * P2x * P4z) - (694575.0 * P3y * P6z) - (992250.0 * P1y * P8x) - (2877525.0 * P1y * P6x * P2z) - (2679075.0 * P1y * P4x * P4z) - (694575.0 * P1y * P2x * P6z) + (99225.0 * P1y * P8z);
    double term_3 = (5760720.0 * P7x * P2z) - (181440.0 * P9x) + (771120.0 * P7x * P2y) - (17304840.0 * P5x * P4z) - (17146080.0 * P5x * P2z * P2y) + (158760.0 * P5x * P4y) + (10220175.0 * P3x * P6z) + (19745775.0 * P3x * P4z * P2y) + (8831025.0 * P3x * P2z * P4y) - (694575.0 * P3x * P6y) - (992250.0 * P1x * P8z) - (2877525.0 * P1x * P6z * P2y) - (2679075.0 * P1x * P4z * P4y) - (694575.0 * P1x * P2z * P6y) + (99225.0 * P1x * P8y);
    double term_4 = (5760720.0 * P7z * P2x) - (181440.0 * P9z) + (771120.0 * P7z * P2y) - (17304840.0 * P5z * P4x) - (17146080.0 * P5z * P2x * P2y) + (158760.0 * P5z * P4y) + (10220175.0 * P3z * P6x) + (19745775.0 * P3z * P4x * P2y) + (8831025.0 * P3z * P2x * P4y) - (694575.0 * P3z * P6y) - (992250.0 * P1z * P8x) - (2877525.0 * P1z * P6x * P2y) - (2679075.0 * P1z * P4x * P4y) - (694575.0 * P1z * P2x * P6y) + (99225.0 * P1z * P8y);
    double term_5 = (5760720.0 * P7y * P2z) - (181440.0 * P9y) + (771120.0 * P7y * P2x) - (17304840.0 * P5y * P4z) - (17146080.0 * P5y * P2z * P2x) + (158760.0 * P5y * P4x) + (10220175.0 * P3y * P6z) + (19745775.0 * P3y * P4z * P2x) + (8831025.0 * P3y * P2z * P4x) - (694575.0 * P3y * P6x) - (992250.0 * P1y * P8z) - (2877525.0 * P1y * P6z * P2x) - (2679075.0 * P1y * P4z * P4x) - (694575.0 * P1y * P2z * P6x) + (99225.0 * P1y * P8x);
    double term_6 = (5760720.0 * P7z * P2y) - (181440.0 * P9z) + (771120.0 * P7z * P2x) - (17304840.0 * P5z * P4y) - (17146080.0 * P5z * P2y * P2x) + (158760.0 * P5z * P4x) + (10220175.0 * P3z * P6y) + (19745775.0 * P3z * P4y * P2x) + (8831025.0 * P3z * P2y * P4x) - (694575.0 * P3z * P6x) - (992250.0 * P1z * P8y) - (2877525.0 * P1z * P6y * P2x) - (2679075.0 * P1z * P4y * P4x) - (694575.0 * P1z * P2y * P6x) + (99225.0 * P1z * P8x);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;

}

void calc_solid_MCSH_9_4_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double term_1 = (4989600.0 * P7x * P1y * P1z) - (17463600.0 * P5x * P3y * P1z) - (17463600.0 * P5x * P1y * P3z) + (10914750.0 * P3x * P5y * P1z) + (21829500.0 * P3x * P3y * P3z) + (10914750.0 * P3x * P1y * P5z) - (1091475.0 * P1x * P7y * P1z) - (3274425.0 * P1x * P5y * P3z) - (3274425.0 * P1x * P3y * P5z) - (1091475.0 * P1x * P1y * P7z);
    double term_2 = (4989600.0 * P7y * P1x * P1z) - (17463600.0 * P5y * P3x * P1z) - (17463600.0 * P5y * P1x * P3z) + (10914750.0 * P3y * P5x * P1z) + (21829500.0 * P3y * P3x * P3z) + (10914750.0 * P3y * P1x * P5z) - (1091475.0 * P1y * P7x * P1z) - (3274425.0 * P1y * P5x * P3z) - (3274425.0 * P1y * P3x * P5z) - (1091475.0 * P1y * P1x * P7z);
    double term_3 = (4989600.0 * P7z * P1x * P1y) - (17463600.0 * P5z * P3x * P1y) - (17463600.0 * P5z * P1x * P3y) + (10914750.0 * P3z * P5x * P1y) + (21829500.0 * P3z * P3x * P3y) + (10914750.0 * P3z * P1x * P5y) - (1091475.0 * P1z * P7x * P1y) - (3274425.0 * P1z * P5x * P3y) - (3274425.0 * P1z * P3x * P5y) - (1091475.0 * P1z * P1x * P7y);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;

}

void calc_solid_MCSH_9_5_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double P8x = P8(lambda, x0, gamma);
    double P8y = P8(lambda, y0, gamma);
    double P8z = P8(lambda, z0, gamma);

    double P9x = P9(lambda, x0, gamma);
    double P9y = P9(lambda, y0, gamma);
    double P9z = P9(lambda, z0, gamma);

    double term_1 = (12020400.0 * P6x * P3y) - (1360800.0 * P8x * P1y) + (2041200.0 * P6x * P1y * P2z) - (16584750.0 * P4x * P5y) - (14458500.0 * P4x * P3y * P2z) + (2126250.0 * P4x * P1y * P4z) + (4380075.0 * P2x * P7y) + (7526925.0 * P2x * P5y * P2z) + (1913625.0 * P2x * P3y * P4z) - (1233225.0 * P2x * P1y * P6z) - (113400.0 * P9y) - (297675.0 * P7y * P2z) - (212625.0 * P5y * P4z) + (14175.0 * P3y * P6z) + (42525.0 * P1y * P8z);
    double term_2 = (12020400.0 * P6y * P3x) - (1360800.0 * P8y * P1x) + (2041200.0 * P6y * P1x * P2z) - (16584750.0 * P4y * P5x) - (14458500.0 * P4y * P3x * P2z) + (2126250.0 * P4y * P1x * P4z) + (4380075.0 * P2y * P7x) + (7526925.0 * P2y * P5x * P2z) + (1913625.0 * P2y * P3x * P4z) - (1233225.0 * P2y * P1x * P6z) - (113400.0 * P9x) - (297675.0 * P7x * P2z) - (212625.0 * P5x * P4z) + (14175.0 * P3x * P6z) + (42525.0 * P1x * P8z);
    double term_3 = (12020400.0 * P6x * P3z) - (1360800.0 * P8x * P1z) + (2041200.0 * P6x * P1z * P2y) - (16584750.0 * P4x * P5z) - (14458500.0 * P4x * P3z * P2y) + (2126250.0 * P4x * P1z * P4y) + (4380075.0 * P2x * P7z) + (7526925.0 * P2x * P5z * P2y) + (1913625.0 * P2x * P3z * P4y) - (1233225.0 * P2x * P1z * P6y) - (113400.0 * P9z) - (297675.0 * P7z * P2y) - (212625.0 * P5z * P4y) + (14175.0 * P3z * P6y) + (42525.0 * P1z * P8y);
    double term_4 = (12020400.0 * P6z * P3x) - (1360800.0 * P8z * P1x) + (2041200.0 * P6z * P1x * P2y) - (16584750.0 * P4z * P5x) - (14458500.0 * P4z * P3x * P2y) + (2126250.0 * P4z * P1x * P4y) + (4380075.0 * P2z * P7x) + (7526925.0 * P2z * P5x * P2y) + (1913625.0 * P2z * P3x * P4y) - (1233225.0 * P2z * P1x * P6y) - (113400.0 * P9x) - (297675.0 * P7x * P2y) - (212625.0 * P5x * P4y) + (14175.0 * P3x * P6y) + (42525.0 * P1x * P8y);
    double term_5 = (12020400.0 * P6y * P3z) - (1360800.0 * P8y * P1z) + (2041200.0 * P6y * P1z * P2x) - (16584750.0 * P4y * P5z) - (14458500.0 * P4y * P3z * P2x) + (2126250.0 * P4y * P1z * P4x) + (4380075.0 * P2y * P7z) + (7526925.0 * P2y * P5z * P2x) + (1913625.0 * P2y * P3z * P4x) - (1233225.0 * P2y * P1z * P6x) - (113400.0 * P9z) - (297675.0 * P7z * P2x) - (212625.0 * P5z * P4x) + (14175.0 * P3z * P6x) + (42525.0 * P1z * P8x);
    double term_6 = (12020400.0 * P6z * P3y) - (1360800.0 * P8z * P1y) + (2041200.0 * P6z * P1y * P2x) - (16584750.0 * P4z * P5y) - (14458500.0 * P4z * P3y * P2x) + (2126250.0 * P4z * P1y * P4x) + (4380075.0 * P2z * P7y) + (7526925.0 * P2z * P5y * P2x) + (1913625.0 * P2z * P3y * P4x) - (1233225.0 * P2z * P1y * P6x) - (113400.0 * P9y) - (297675.0 * P7y * P2x) - (212625.0 * P5y * P4x) + (14175.0 * P3y * P6x) + (42525.0 * P1y * P8x);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;

}

void calc_solid_MCSH_9_6_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double P8x = P8(lambda, x0, gamma);
    double P8y = P8(lambda, y0, gamma);
    double P8z = P8(lambda, z0, gamma);

    double P9x = P9(lambda, x0, gamma);
    double P9y = P9(lambda, y0, gamma);
    double P9z = P9(lambda, z0, gamma);

    double term_1 = (10659600.0 * P6x * P2y * P1z) - (453600.0 * P8x * P1z) + (680400.0 * P6x * P3z) - (18002250.0 * P4x * P4y * P1z) - (17293500.0 * P4x * P2y * P3z) + (708750.0 * P4x * P5z) + (5202225.0 * P2x * P6y * P1z) + (9993375.0 * P2x * P4y * P3z) + (4380075.0 * P2x * P2y * P5z) - (411075.0 * P2x * P7z) - (141750.0 * P8y * P1z) - (411075.0 * P6y * P3z) - (382725.0 * P4y * P5z) - (99225.0 * P2y * P7z) + (14175.0 * P9z);
    double term_2 = (10659600.0 * P6y * P2x * P1z) - (453600.0 * P8y * P1z) + (680400.0 * P6y * P3z) - (18002250.0 * P4y * P4x * P1z) - (17293500.0 * P4y * P2x * P3z) + (708750.0 * P4y * P5z) + (5202225.0 * P2y * P6x * P1z) + (9993375.0 * P2y * P4x * P3z) + (4380075.0 * P2y * P2x * P5z) - (411075.0 * P2y * P7z) - (141750.0 * P8x * P1z) - (411075.0 * P6x * P3z) - (382725.0 * P4x * P5z) - (99225.0 * P2x * P7z) + (14175.0 * P9z);
    double term_3 = (10659600.0 * P6x * P2z * P1y) - (453600.0 * P8x * P1y) + (680400.0 * P6x * P3y) - (18002250.0 * P4x * P4z * P1y) - (17293500.0 * P4x * P2z * P3y) + (708750.0 * P4x * P5y) + (5202225.0 * P2x * P6z * P1y) + (9993375.0 * P2x * P4z * P3y) + (4380075.0 * P2x * P2z * P5y) - (411075.0 * P2x * P7y) - (141750.0 * P8z * P1y) - (411075.0 * P6z * P3y) - (382725.0 * P4z * P5y) - (99225.0 * P2z * P7y) + (14175.0 * P9y);
    double term_4 = (10659600.0 * P6z * P2x * P1y) - (453600.0 * P8z * P1y) + (680400.0 * P6z * P3y) - (18002250.0 * P4z * P4x * P1y) - (17293500.0 * P4z * P2x * P3y) + (708750.0 * P4z * P5y) + (5202225.0 * P2z * P6x * P1y) + (9993375.0 * P2z * P4x * P3y) + (4380075.0 * P2z * P2x * P5y) - (411075.0 * P2z * P7y) - (141750.0 * P8x * P1y) - (411075.0 * P6x * P3y) - (382725.0 * P4x * P5y) - (99225.0 * P2x * P7y) + (14175.0 * P9y);
    double term_5 = (10659600.0 * P6y * P2z * P1x) - (453600.0 * P8y * P1x) + (680400.0 * P6y * P3x) - (18002250.0 * P4y * P4z * P1x) - (17293500.0 * P4y * P2z * P3x) + (708750.0 * P4y * P5x) + (5202225.0 * P2y * P6z * P1x) + (9993375.0 * P2y * P4z * P3x) + (4380075.0 * P2y * P2z * P5x) - (411075.0 * P2y * P7x) - (141750.0 * P8z * P1x) - (411075.0 * P6z * P3x) - (382725.0 * P4z * P5x) - (99225.0 * P2z * P7x) + (14175.0 * P9x);
    double term_6 = (10659600.0 * P6z * P2y * P1x) - (453600.0 * P8z * P1x) + (680400.0 * P6z * P3x) - (18002250.0 * P4z * P4y * P1x) - (17293500.0 * P4z * P2y * P3x) + (708750.0 * P4z * P5x) + (5202225.0 * P2z * P6y * P1x) + (9993375.0 * P2z * P4y * P3x) + (4380075.0 * P2z * P2y * P5x) - (411075.0 * P2z * P7x) - (141750.0 * P8y * P1x) - (411075.0 * P6y * P3x) - (382725.0 * P4y * P5x) - (99225.0 * P2y * P7x) + (14175.0 * P9x);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;

}

void calc_solid_MCSH_9_7_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double P8x = P8(lambda, x0, gamma);
    double P8y = P8(lambda, y0, gamma);
    double P8z = P8(lambda, z0, gamma);

    double P9x = P9(lambda, x0, gamma);
    double P9y = P9(lambda, y0, gamma);
    double P9z = P9(lambda, z0, gamma);

    double term_1 = (19486005.0 * P5x * P4y) - (3412640.0 * P7x * P2y) + (5292210.0 * P5x * P2y * P2z) + (518980.0 * P9x) + (1576960.0 * P7x * P2z) + (2022405.0 * P5x * P4z) - (9524900.0 * P3x * P6y) - (1443750.0 * P3x * P4y * P2z) + (9471000.0 * P3x * P2y * P4z) + (1389850.0 * P3x * P6z) + (1516900.0 * P1x * P8y) + (2949100.0 * P1x * P6y * P2z) + (1772925.0 * P1x * P4y * P4z) + (766150.0 * P1x * P2y * P6z) + (425425.0 * P1x * P8z);
    double term_2 = (19486005.0 * P5y * P4x) - (3412640.0 * P7y * P2x) + (5292210.0 * P5y * P2x * P2z) + (518980.0 * P9y) + (1576960.0 * P7y * P2z) + (2022405.0 * P5y * P4z) - (9524900.0 * P3y * P6x) - (1443750.0 * P3y * P4x * P2z) + (9471000.0 * P3y * P2x * P4z) + (1389850.0 * P3y * P6z) + (1516900.0 * P1y * P8x) + (2949100.0 * P1y * P6x * P2z) + (1772925.0 * P1y * P4x * P4z) + (766150.0 * P1y * P2x * P6z) + (425425.0 * P1y * P8z);
    double term_3 = (19486005.0 * P5x * P4z) - (3412640.0 * P7x * P2z) + (5292210.0 * P5x * P2z * P2y) + (518980.0 * P9x) + (1576960.0 * P7x * P2y) + (2022405.0 * P5x * P4y) - (9524900.0 * P3x * P6z) - (1443750.0 * P3x * P4z * P2y) + (9471000.0 * P3x * P2z * P4y) + (1389850.0 * P3x * P6y) + (1516900.0 * P1x * P8z) + (2949100.0 * P1x * P6z * P2y) + (1772925.0 * P1x * P4z * P4y) + (766150.0 * P1x * P2z * P6y) + (425425.0 * P1x * P8y);
    double term_4 = (19486005.0 * P5z * P4x) - (3412640.0 * P7z * P2x) + (5292210.0 * P5z * P2x * P2y) + (518980.0 * P9z) + (1576960.0 * P7z * P2y) + (2022405.0 * P5z * P4y) - (9524900.0 * P3z * P6x) - (1443750.0 * P3z * P4x * P2y) + (9471000.0 * P3z * P2x * P4y) + (1389850.0 * P3z * P6y) + (1516900.0 * P1z * P8x) + (2949100.0 * P1z * P6x * P2y) + (1772925.0 * P1z * P4x * P4y) + (766150.0 * P1z * P2x * P6y) + (425425.0 * P1z * P8y);
    double term_5 = (19486005.0 * P5y * P4z) - (3412640.0 * P7y * P2z) + (5292210.0 * P5y * P2z * P2x) + (518980.0 * P9y) + (1576960.0 * P7y * P2x) + (2022405.0 * P5y * P4x) - (9524900.0 * P3y * P6z) - (1443750.0 * P3y * P4z * P2x) + (9471000.0 * P3y * P2z * P4x) + (1389850.0 * P3y * P6x) + (1516900.0 * P1y * P8z) + (2949100.0 * P1y * P6z * P2x) + (1772925.0 * P1y * P4z * P4x) + (766150.0 * P1y * P2z * P6x) + (425425.0 * P1y * P8x);
    double term_6 = (19486005.0 * P5z * P4y) - (3412640.0 * P7z * P2y) + (5292210.0 * P5z * P2y * P2x) + (518980.0 * P9z) + (1576960.0 * P7z * P2x) + (2022405.0 * P5z * P4x) - (9524900.0 * P3z * P6y) - (1443750.0 * P3z * P4y * P2x) + (9471000.0 * P3z * P2y * P4x) + (1389850.0 * P3z * P6x) + (1516900.0 * P1z * P8y) + (2949100.0 * P1z * P6y * P2x) + (1772925.0 * P1z * P4y * P4x) + (766150.0 * P1z * P2y * P6x) + (425425.0 * P1z * P8x);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;

}

void calc_solid_MCSH_9_8_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double term_1 = (16839900.0 * P5x * P3y * P1z) - (2494800.0 * P7x * P1y * P1z) + (623700.0 * P5x * P1y * P3z) - (13565475.0 * P3x * P5y * P1z) - (10914750.0 * P3x * P3y * P3z) + (2650725.0 * P3x * P1y * P5z) + (1559250.0 * P1x * P7y * P1z) + (2650725.0 * P1x * P5y * P3z) + (623700.0 * P1x * P3y * P5z) - (467775.0 * P1x * P1y * P7z);
    double term_2 = (16839900.0 * P5y * P3x * P1z) - (2494800.0 * P7y * P1x * P1z) + (623700.0 * P5y * P1x * P3z) - (13565475.0 * P3y * P5x * P1z) - (10914750.0 * P3y * P3x * P3z) + (2650725.0 * P3y * P1x * P5z) + (1559250.0 * P1y * P7x * P1z) + (2650725.0 * P1y * P5x * P3z) + (623700.0 * P1y * P3x * P5z) - (467775.0 * P1y * P1x * P7z);
    double term_3 = (16839900.0 * P5x * P3z * P1y) - (2494800.0 * P7x * P1z * P1y) + (623700.0 * P5x * P1z * P3y) - (13565475.0 * P3x * P5z * P1y) - (10914750.0 * P3x * P3z * P3y) + (2650725.0 * P3x * P1z * P5y) + (1559250.0 * P1x * P7z * P1y) + (2650725.0 * P1x * P5z * P3y) + (623700.0 * P1x * P3z * P5y) - (467775.0 * P1x * P1z * P7y);
    double term_4 = (16839900.0 * P5z * P3x * P1y) - (2494800.0 * P7z * P1x * P1y) + (623700.0 * P5z * P1x * P3y) - (13565475.0 * P3z * P5x * P1y) - (10914750.0 * P3z * P3x * P3y) + (2650725.0 * P3z * P1x * P5y) + (1559250.0 * P1z * P7x * P1y) + (2650725.0 * P1z * P5x * P3y) + (623700.0 * P1z * P3x * P5y) - (467775.0 * P1z * P1x * P7y);
    double term_5 = (16839900.0 * P5y * P3z * P1x) - (2494800.0 * P7y * P1z * P1x) + (623700.0 * P5y * P1z * P3x) - (13565475.0 * P3y * P5z * P1x) - (10914750.0 * P3y * P3z * P3x) + (2650725.0 * P3y * P1z * P5x) + (1559250.0 * P1y * P7z * P1x) + (2650725.0 * P1y * P5z * P3x) + (623700.0 * P1y * P3z * P5x) - (467775.0 * P1y * P1z * P7x);
    double term_6 = (16839900.0 * P5z * P3y * P1x) - (2494800.0 * P7z * P1y * P1x) + (623700.0 * P5z * P1y * P3x) - (13565475.0 * P3z * P5y * P1x) - (10914750.0 * P3z * P3y * P3x) + (2650725.0 * P3z * P1y * P5x) + (1559250.0 * P1z * P7y * P1x) + (2650725.0 * P1z * P5y * P3x) + (623700.0 * P1z * P3y * P5x) - (467775.0 * P1z * P1y * P7x);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;

}

void calc_solid_MCSH_9_9_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double P8x = P8(lambda, x0, gamma);
    double P8y = P8(lambda, y0, gamma);
    double P8z = P8(lambda, z0, gamma);

    double P9x = P9(lambda, x0, gamma);
    double P9y = P9(lambda, y0, gamma);
    double P9z = P9(lambda, z0, gamma);

    double term_1 = (16448670.0 * P5x * P2y * P2z) - (816480.0 * P7x * P2y) + (116235.0 * P5x * P4y) - (816480.0 * P7x * P2z) + (116235.0 * P5x * P4z) + (45360.0 * P9x) - (13707225.0 * P3x * P4y * P2z) - (13707225.0 * P3x * P2y * P4z) + (836325.0 * P3x * P6y) + (836325.0 * P3x * P6z) + (1460025.0 * P1x * P6y * P2z) + (3203550.0 * P1x * P4y * P4z) + (1460025.0 * P1x * P2y * P6z) - (141750.0 * P1x * P8y) - (141750.0 * P1x * P8z);
    double term_2 = (16448670.0 * P5y * P2x * P2z) - (816480.0 * P7y * P2x) + (116235.0 * P5y * P4x) - (816480.0 * P7y * P2z) + (116235.0 * P5y * P4z) + (45360.0 * P9y) - (13707225.0 * P3y * P4x * P2z) - (13707225.0 * P3y * P2x * P4z) + (836325.0 * P3y * P6x) + (836325.0 * P3y * P6z) + (1460025.0 * P1y * P6x * P2z) + (3203550.0 * P1y * P4x * P4z) + (1460025.0 * P1y * P2x * P6z) - (141750.0 * P1y * P8x) - (141750.0 * P1y * P8z);
    double term_3 = (16448670.0 * P5z * P2x * P2y) - (816480.0 * P7z * P2x) + (116235.0 * P5z * P4x) - (816480.0 * P7z * P2y) + (116235.0 * P5z * P4y) + (45360.0 * P9z) - (13707225.0 * P3z * P4x * P2y) - (13707225.0 * P3z * P2x * P4y) + (836325.0 * P3z * P6x) + (836325.0 * P3z * P6y) + (1460025.0 * P1z * P6x * P2y) + (3203550.0 * P1z * P4x * P4y) + (1460025.0 * P1z * P2x * P6y) - (141750.0 * P1z * P8x) - (141750.0 * P1z * P8y);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;

}

void calc_solid_MCSH_9_10_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double P8x = P8(lambda, x0, gamma);
    double P8y = P8(lambda, y0, gamma);
    double P8z = P8(lambda, z0, gamma);

    double P9x = P9(lambda, x0, gamma);
    double P9y = P9(lambda, y0, gamma);
    double P9z = P9(lambda, z0, gamma);

    double term_1 = (19604025.0 * P4x * P4y * P1z) - (7200900.0 * P6x * P2y * P1z) - (3203550.0 * P4x * P2y * P3z) + (226800.0 * P8x * P1z) + (283500.0 * P6x * P3z) - (104895.0 * P4x * P5z) - (7200900.0 * P2x * P6y * P1z) - (3203550.0 * P2x * P4y * P3z) + (3844260.0 * P2x * P2y * P5z) - (153090.0 * P2x * P7z) + (226800.0 * P8y * P1z) + (283500.0 * P6y * P3z) - (104895.0 * P4y * P5z) - (153090.0 * P2y * P7z) + (8505.0 * P9z);
    double term_2 = (19604025.0 * P4x * P4z * P1y) - (7200900.0 * P6x * P2z * P1y) - (3203550.0 * P4x * P2z * P3y) + (226800.0 * P8x * P1y) + (283500.0 * P6x * P3y) - (104895.0 * P4x * P5y) - (7200900.0 * P2x * P6z * P1y) - (3203550.0 * P2x * P4z * P3y) + (3844260.0 * P2x * P2z * P5y) - (153090.0 * P2x * P7y) + (226800.0 * P8z * P1y) + (283500.0 * P6z * P3y) - (104895.0 * P4z * P5y) - (153090.0 * P2z * P7y) + (8505.0 * P9y);
    double term_3 = (19604025.0 * P4y * P4z * P1x) - (7200900.0 * P6y * P2z * P1x) - (3203550.0 * P4y * P2z * P3x) + (226800.0 * P8y * P1x) + (283500.0 * P6y * P3x) - (104895.0 * P4y * P5x) - (7200900.0 * P2y * P6z * P1x) - (3203550.0 * P2y * P4z * P3x) + (3844260.0 * P2y * P2z * P5x) - (153090.0 * P2y * P7x) + (226800.0 * P8z * P1x) + (283500.0 * P6z * P3x) - (104895.0 * P4z * P5x) - (153090.0 * P2z * P7x) + (8505.0 * P9x);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;

}

void calc_solid_MCSH_9_11_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double P8x = P8(lambda, x0, gamma);
    double P8y = P8(lambda, y0, gamma);
    double P8z = P8(lambda, z0, gamma);

    double P9x = P9(lambda, x0, gamma);
    double P9y = P9(lambda, y0, gamma);
    double P9z = P9(lambda, z0, gamma);

    double term_1 = (20497050.0 * P4x * P3y * P2z) - (963900.0 * P6x * P3y) - (603855.0 * P4x * P5y) - (3458700.0 * P6x * P1y * P2z) - (1601775.0 * P4x * P1y * P4z) + (226800.0 * P8x * P1y) - (8224335.0 * P2x * P5y * P2z) - (6789825.0 * P2x * P3y * P4z) + (564165.0 * P2x * P7y) + (1998675.0 * P2x * P1y * P6z) + (252315.0 * P7y * P2z) + (487620.0 * P5y * P4z) + (127575.0 * P3y * P6z) - (22680.0 * P9y) - (85050.0 * P1y * P8z);
    double term_2 = (20497050.0 * P4y * P3x * P2z) - (963900.0 * P6y * P3x) - (603855.0 * P4y * P5x) - (3458700.0 * P6y * P1x * P2z) - (1601775.0 * P4y * P1x * P4z) + (226800.0 * P8y * P1x) - (8224335.0 * P2y * P5x * P2z) - (6789825.0 * P2y * P3x * P4z) + (564165.0 * P2y * P7x) + (1998675.0 * P2y * P1x * P6z) + (252315.0 * P7x * P2z) + (487620.0 * P5x * P4z) + (127575.0 * P3x * P6z) - (22680.0 * P9x) - (85050.0 * P1x * P8z);
    double term_3 = (20497050.0 * P4x * P3z * P2y) - (963900.0 * P6x * P3z) - (603855.0 * P4x * P5z) - (3458700.0 * P6x * P1z * P2y) - (1601775.0 * P4x * P1z * P4y) + (226800.0 * P8x * P1z) - (8224335.0 * P2x * P5z * P2y) - (6789825.0 * P2x * P3z * P4y) + (564165.0 * P2x * P7z) + (1998675.0 * P2x * P1z * P6y) + (252315.0 * P7z * P2y) + (487620.0 * P5z * P4y) + (127575.0 * P3z * P6y) - (22680.0 * P9z) - (85050.0 * P1z * P8y);
    double term_4 = (20497050.0 * P4z * P3x * P2y) - (963900.0 * P6z * P3x) - (603855.0 * P4z * P5x) - (3458700.0 * P6z * P1x * P2y) - (1601775.0 * P4z * P1x * P4y) + (226800.0 * P8z * P1x) - (8224335.0 * P2z * P5x * P2y) - (6789825.0 * P2z * P3x * P4y) + (564165.0 * P2z * P7x) + (1998675.0 * P2z * P1x * P6y) + (252315.0 * P7x * P2y) + (487620.0 * P5x * P4y) + (127575.0 * P3x * P6y) - (22680.0 * P9x) - (85050.0 * P1x * P8y);
    double term_5 = (20497050.0 * P4y * P3z * P2x) - (963900.0 * P6y * P3z) - (603855.0 * P4y * P5z) - (3458700.0 * P6y * P1z * P2x) - (1601775.0 * P4y * P1z * P4x) + (226800.0 * P8y * P1z) - (8224335.0 * P2y * P5z * P2x) - (6789825.0 * P2y * P3z * P4x) + (564165.0 * P2y * P7z) + (1998675.0 * P2y * P1z * P6x) + (252315.0 * P7z * P2x) + (487620.0 * P5z * P4x) + (127575.0 * P3z * P6x) - (22680.0 * P9z) - (85050.0 * P1z * P8x);
    double term_6 = (20497050.0 * P4z * P3y * P2x) - (963900.0 * P6z * P3y) - (603855.0 * P4z * P5y) - (3458700.0 * P6z * P1y * P2x) - (1601775.0 * P4z * P1y * P4x) + (226800.0 * P8z * P1y) - (8224335.0 * P2z * P5y * P2x) - (6789825.0 * P2z * P3y * P4x) + (564165.0 * P2z * P7y) + (1998675.0 * P2z * P1y * P6x) + (252315.0 * P7y * P2x) + (487620.0 * P5y * P4x) + (127575.0 * P3y * P6x) - (22680.0 * P9y) - (85050.0 * P1y * P8x);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;

}

void calc_solid_MCSH_9_12_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double term = (21829320.0 * P3x * P3y * P3z) - (3274515.0 * P5x * P3y * P1z) - (3274605.0 * P3x * P5y * P1z) - (3274425.0 * P5x * P1y * P3z) - (3274425.0 * P3x * P1y * P5z) + (935550.0 * P7x * P1y * P1z) - (3274605.0 * P1x * P5y * P3z) - (3274515.0 * P1x * P3y * P5z) + (935460.0 * P1x * P7y * P1z) + (935550.0 * P1x * P1y * P7z);

    double m = temp * term;

    value[0] = m;

}



SolidGMPFunction get_solid_mcsh_function(int mcsh_order, int group_num)
{
    SolidGMPFunction result;

    if (mcsh_order == 0) {
        if (group_num == 1) {
            result = calc_solid_MCSH_0_1;
        }
    } else if (mcsh_order == 1) {
        if (group_num == 1) {
            result = calc_solid_MCSH_1_1;
        }
    } else if (mcsh_order == 2) {
        if (group_num == 1) {
            result = calc_solid_MCSH_2_1;
        } else if (group_num == 2){
            result = calc_solid_MCSH_2_2;
        }
    } else if (mcsh_order == 3) {
        if (group_num == 1) {
            result = calc_solid_MCSH_3_1;
        } else if (group_num == 2){
            result = calc_solid_MCSH_3_2;
        } else if (group_num == 3){
            result = calc_solid_MCSH_3_3;
        }
    } else if (mcsh_order == 4) {
        if (group_num == 1) {
            result = calc_solid_MCSH_4_1;
        } else if (group_num == 2){
            result = calc_solid_MCSH_4_2;
        } else if (group_num == 3){
            result = calc_solid_MCSH_4_3;
        } else if (group_num == 4){
            result = calc_solid_MCSH_4_4;
        }
    } else if (mcsh_order == 5) {
        if (group_num == 1) {
            result = calc_solid_MCSH_5_1;
        } else if (group_num == 2){
            result = calc_solid_MCSH_5_2;
        } else if (group_num == 3){
            result = calc_solid_MCSH_5_3;
        } else if (group_num == 4){
            result = calc_solid_MCSH_5_4;
        } else if (group_num == 5){
            result = calc_solid_MCSH_5_5;
        }
    } else if (mcsh_order == 6) {
        if (group_num == 1) {
            result = calc_solid_MCSH_6_1;
        } else if (group_num == 2){
            result = calc_solid_MCSH_6_2;
        } else if (group_num == 3){
            result = calc_solid_MCSH_6_3;
        } else if (group_num == 4){
            result = calc_solid_MCSH_6_4;
        } else if (group_num == 5){
            result = calc_solid_MCSH_6_5;
        } else if (group_num == 6){
            result = calc_solid_MCSH_6_6;
        } else if (group_num == 7){
            result = calc_solid_MCSH_6_7;
        }
    } else if (mcsh_order == 7) {
        if (group_num == 1) {
            result = calc_solid_MCSH_7_1;
        } else if (group_num == 2){
            result = calc_solid_MCSH_7_2;
        } else if (group_num == 3){
            result = calc_solid_MCSH_7_3;
        } else if (group_num == 4){
            result = calc_solid_MCSH_7_4;
        } else if (group_num == 5){
            result = calc_solid_MCSH_7_5;
        } else if (group_num == 6){
            result = calc_solid_MCSH_7_6;
        } else if (group_num == 7){
            result = calc_solid_MCSH_7_7;
        } else if (group_num == 8){
            result = calc_solid_MCSH_7_8;
        }
    } else if (mcsh_order == 8) {
        if (group_num == 1) {
            result = calc_solid_MCSH_8_1;
        } else if (group_num == 2){
            result = calc_solid_MCSH_8_2;
        } else if (group_num == 3){
            result = calc_solid_MCSH_8_3;
        } else if (group_num == 4){
            result = calc_solid_MCSH_8_4;
        } else if (group_num == 5){
            result = calc_solid_MCSH_8_5;
        } else if (group_num == 6){
            result = calc_solid_MCSH_8_6;
        } else if (group_num == 7){
            result = calc_solid_MCSH_8_7;
        } else if (group_num == 8){
            result = calc_solid_MCSH_8_8;
        } else if (group_num == 9){
            result = calc_solid_MCSH_8_9;
        } else if (group_num == 10){
            result = calc_solid_MCSH_8_10;
        }
    } else if (mcsh_order == 9) {
        if (group_num == 1) {
            result = calc_solid_MCSH_9_1;
        } else if (group_num == 2){
            result = calc_solid_MCSH_9_2;
        } else if (group_num == 3){
            result = calc_solid_MCSH_9_3;
        } else if (group_num == 4){
            result = calc_solid_MCSH_9_4;
        } else if (group_num == 5){
            result = calc_solid_MCSH_9_5;
        } else if (group_num == 6){
            result = calc_solid_MCSH_9_6;
        } else if (group_num == 7){
            result = calc_solid_MCSH_9_7;
        } else if (group_num == 8){
            result = calc_solid_MCSH_9_8;
        } else if (group_num == 9){
            result = calc_solid_MCSH_9_9;
        } else if (group_num == 10){
            result = calc_solid_MCSH_9_10;
        } else if (group_num == 11){
            result = calc_solid_MCSH_9_11;
        } else if (group_num == 12){
            result = calc_solid_MCSH_9_12;
        }
    }

    return result;
}




void calc_solid_MCSH_0_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
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

void calc_solid_MCSH_1_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

    double P1x = P1(lambda, x0, gamma);
    double P1y = P1(lambda, y0, gamma);
    double P1z = P1(lambda, z0, gamma);

    double term_1 = (1.0 * P1x);
    double term_2 = (1.0 * P1y);
    double term_3 = (1.0 * P1z);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dterm_1_dx = (1.0 * dP1x_exp);
    double dterm_1_dy = (1.0 * P1x * dP0y_exp);
    double dterm_1_dz = (1.0 * P1x * dP0z_exp);

    double dterm_2_dx = (1.0 * dP0x_exp * P1y);
    double dterm_2_dy = (1.0 * dP1y_exp);
    double dterm_2_dz = (1.0 * P1y * dP0z_exp);

    double dterm_3_dx = (1.0 * dP0x_exp * P1z);
    double dterm_3_dy = (1.0 * dP0y_exp * P1z);
    double dterm_3_dz = (1.0 * dP1z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

}


void calc_solid_MCSH_2_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double term_1 = (2.0 * P2x) - (1.0 * P2y) - (1.0 * P2z);
    double term_2 = (2.0 * P2y) - (1.0 * P2x) - (1.0 * P2z);
    double term_3 = (2.0 * P2z) - (1.0 * P2x) - (1.0 * P2y);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dterm_1_dx = (2.0 * dP2x_exp) - (1.0 * dP0x_exp * P2y) - (1.0 * dP0x_exp * P2z);
    double dterm_1_dy = (2.0 * P2x * dP0y_exp) - (1.0 * dP2y_exp) - (1.0 * dP0y_exp * P2z);
    double dterm_1_dz = (2.0 * P2x * dP0z_exp) - (1.0 * P2y * dP0z_exp) - (1.0 * dP2z_exp);

    double dterm_2_dx = (2.0 * dP0x_exp * P2y) - (1.0 * dP2x_exp) - (1.0 * dP0x_exp * P2z);
    double dterm_2_dy = (2.0 * dP2y_exp) - (1.0 * P2x * dP0y_exp) - (1.0 * dP0y_exp * P2z);
    double dterm_2_dz = (2.0 * P2y * dP0z_exp) - (1.0 * P2x * dP0z_exp) - (1.0 * dP2z_exp);

    double dterm_3_dx = (2.0 * dP0x_exp * P2z) - (1.0 * dP2x_exp) - (1.0 * dP0x_exp * P2y);
    double dterm_3_dy = (2.0 * dP0y_exp * P2z) - (1.0 * P2x * dP0y_exp) - (1.0 * dP2y_exp);
    double dterm_3_dz = (2.0 * dP2z_exp) - (1.0 * P2x * dP0z_exp) - (1.0 * P2y * dP0z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

}

void calc_solid_MCSH_2_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

    double P1x = P1(lambda, x0, gamma);
    double P1y = P1(lambda, y0, gamma);
    double P1z = P1(lambda, z0, gamma);

    double term_1 = (3.0 * P1x * P1y);
    double term_2 = (3.0 * P1x * P1z);
    double term_3 = (3.0 * P1y * P1z);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dterm_1_dx = (3.0 * dP1x_exp * P1y);
    double dterm_1_dy = (3.0 * P1x * dP1y_exp);
    double dterm_1_dz = (3.0 * P1x * P1y * dP0z_exp);

    double dterm_2_dx = (3.0 * dP1x_exp * P1z);
    double dterm_2_dy = (3.0 * P1x * dP0y_exp * P1z);
    double dterm_2_dz = (3.0 * P1x * dP1z_exp);

    double dterm_3_dx = (3.0 * dP0x_exp * P1y * P1z);
    double dterm_3_dy = (3.0 * dP1y_exp * P1z);
    double dterm_3_dz = (3.0 * P1y * dP1z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

}


void calc_solid_MCSH_3_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

    double P1x = P1(lambda, x0, gamma);
    double P1y = P1(lambda, y0, gamma);
    double P1z = P1(lambda, z0, gamma);

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double P3x = P3(lambda, x0, gamma);
    double P3y = P3(lambda, y0, gamma);
    double P3z = P3(lambda, z0, gamma);

    double term_1 = (6.0 * P3x) - (9.0 * P1x * P2y) - (9.0 * P1x * P2z);
    double term_2 = (6.0 * P3y) - (9.0 * P2x * P1y) - (9.0 * P1y * P2z);
    double term_3 = (6.0 * P3z) - (9.0 * P2x * P1z) - (9.0 * P2y * P1z);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP3x_exp = dP3_exp(P3x, C2, lambda, x0, gamma);
    double dP3y_exp = dP3_exp(P3y, C2, lambda, y0, gamma);
    double dP3z_exp = dP3_exp(P3z, C2, lambda, z0, gamma);

    double dterm_1_dx = (6.0 * dP3x_exp) - (9.0 * dP1x_exp * P2y) - (9.0 * dP1x_exp * P2z);
    double dterm_1_dy = (6.0 * P3x * dP0y_exp) - (9.0 * P1x * dP2y_exp) - (9.0 * P1x * dP0y_exp * P2z);
    double dterm_1_dz = (6.0 * P3x * dP0z_exp) - (9.0 * P1x * P2y * dP0z_exp) - (9.0 * P1x * dP2z_exp);

    double dterm_2_dx = (6.0 * dP0x_exp * P3y) - (9.0 * dP2x_exp * P1y) - (9.0 * dP0x_exp * P1y * P2z);
    double dterm_2_dy = (6.0 * dP3y_exp) - (9.0 * P2x * dP1y_exp) - (9.0 * dP1y_exp * P2z);
    double dterm_2_dz = (6.0 * P3y * dP0z_exp) - (9.0 * P2x * P1y * dP0z_exp) - (9.0 * P1y * dP2z_exp);

    double dterm_3_dx = (6.0 * dP0x_exp * P3z) - (9.0 * dP2x_exp * P1z) - (9.0 * dP0x_exp * P2y * P1z);
    double dterm_3_dy = (6.0 * dP0y_exp * P3z) - (9.0 * P2x * dP0y_exp * P1z) - (9.0 * dP2y_exp * P1z);
    double dterm_3_dz = (6.0 * dP3z_exp) - (9.0 * P2x * dP1z_exp) - (9.0 * P2y * dP1z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

}

void calc_solid_MCSH_3_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

    double P1x = P1(lambda, x0, gamma);
    double P1y = P1(lambda, y0, gamma);
    double P1z = P1(lambda, z0, gamma);

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double P3x = P3(lambda, x0, gamma);
    double P3y = P3(lambda, y0, gamma);
    double P3z = P3(lambda, z0, gamma);

    double term_1 = (12.0 * P2x * P1y) - (3.0 * P3y) - (3.0 * P1y * P2z);
    double term_2 = (12.0 * P1x * P2y) - (3.0 * P3x) - (3.0 * P1x * P2z);
    double term_3 = (12.0 * P2x * P1z) - (3.0 * P3z) - (3.0 * P2y * P1z);
    double term_4 = (12.0 * P1x * P2z) - (3.0 * P3x) - (3.0 * P1x * P2y);
    double term_5 = (12.0 * P2y * P1z) - (3.0 * P3z) - (3.0 * P2x * P1z);
    double term_6 = (12.0 * P1y * P2z) - (3.0 * P3y) - (3.0 * P2x * P1y);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP3x_exp = dP3_exp(P3x, C2, lambda, x0, gamma);
    double dP3y_exp = dP3_exp(P3y, C2, lambda, y0, gamma);
    double dP3z_exp = dP3_exp(P3z, C2, lambda, z0, gamma);

    double dterm_1_dx = (12.0 * dP2x_exp * P1y) - (3.0 * dP0x_exp * P3y) - (3.0 * dP0x_exp * P1y * P2z);
    double dterm_1_dy = (12.0 * P2x * dP1y_exp) - (3.0 * dP3y_exp) - (3.0 * dP1y_exp * P2z);
    double dterm_1_dz = (12.0 * P2x * P1y * dP0z_exp) - (3.0 * P3y * dP0z_exp) - (3.0 * P1y * dP2z_exp);

    double dterm_2_dx = (12.0 * dP1x_exp * P2y) - (3.0 * dP3x_exp) - (3.0 * dP1x_exp * P2z);
    double dterm_2_dy = (12.0 * P1x * dP2y_exp) - (3.0 * P3x * dP0y_exp) - (3.0 * P1x * dP0y_exp * P2z);
    double dterm_2_dz = (12.0 * P1x * P2y * dP0z_exp) - (3.0 * P3x * dP0z_exp) - (3.0 * P1x * dP2z_exp);

    double dterm_3_dx = (12.0 * dP2x_exp * P1z) - (3.0 * dP0x_exp * P3z) - (3.0 * dP0x_exp * P2y * P1z);
    double dterm_3_dy = (12.0 * P2x * dP0y_exp * P1z) - (3.0 * dP0y_exp * P3z) - (3.0 * dP2y_exp * P1z);
    double dterm_3_dz = (12.0 * P2x * dP1z_exp) - (3.0 * dP3z_exp) - (3.0 * P2y * dP1z_exp);

    double dterm_4_dx = (12.0 * dP1x_exp * P2z) - (3.0 * dP3x_exp) - (3.0 * dP1x_exp * P2y);
    double dterm_4_dy = (12.0 * P1x * dP0y_exp * P2z) - (3.0 * P3x * dP0y_exp) - (3.0 * P1x * dP2y_exp);
    double dterm_4_dz = (12.0 * P1x * dP2z_exp) - (3.0 * P3x * dP0z_exp) - (3.0 * P1x * P2y * dP0z_exp);

    double dterm_5_dx = (12.0 * dP0x_exp * P2y * P1z) - (3.0 * dP0x_exp * P3z) - (3.0 * dP2x_exp * P1z);
    double dterm_5_dy = (12.0 * dP2y_exp * P1z) - (3.0 * dP0y_exp * P3z) - (3.0 * P2x * dP0y_exp * P1z);
    double dterm_5_dz = (12.0 * P2y * dP1z_exp) - (3.0 * dP3z_exp) - (3.0 * P2x * dP1z_exp);

    double dterm_6_dx = (12.0 * dP0x_exp * P1y * P2z) - (3.0 * dP0x_exp * P3y) - (3.0 * dP2x_exp * P1y);
    double dterm_6_dy = (12.0 * dP1y_exp * P2z) - (3.0 * dP3y_exp) - (3.0 * P2x * dP1y_exp);
    double dterm_6_dz = (12.0 * P1y * dP2z_exp) - (3.0 * P3y * dP0z_exp) - (3.0 * P2x * P1y * dP0z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

    // dmiu4 dx/dy/dz
    deriv[9] = temp * dterm_4_dx;
    deriv[10] = temp * dterm_4_dy;
    deriv[11] = temp * dterm_4_dz;

    // dmiu5 dx/dy/dz
    deriv[12] = temp * dterm_5_dx;
    deriv[13] = temp * dterm_5_dy;
    deriv[14] = temp * dterm_5_dz;

    // dmiu6 dx/dy/dz
    deriv[15] = temp * dterm_6_dx;
    deriv[16] = temp * dterm_6_dy;
    deriv[17] = temp * dterm_6_dz;

}

void calc_solid_MCSH_3_3(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

    double P1x = P1(lambda, x0, gamma);
    double P1y = P1(lambda, y0, gamma);
    double P1z = P1(lambda, z0, gamma);

    double term_1 = (15.0 * P1x * P1y * P1z);

    double m = temp * term_1;

    value[0] = m;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dterm_1_dx = (15.0 * dP1x_exp * P1y * P1z);
    double dterm_1_dy = (15.0 * P1x * dP1y_exp * P1z);
    double dterm_1_dz = (15.0 * P1x * P1y * dP1z_exp);



    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

}


void calc_solid_MCSH_4_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double P4x = P4(lambda, x0, gamma);
    double P4y = P4(lambda, y0, gamma);
    double P4z = P4(lambda, z0, gamma);

    double term_1 = (24.0 * P4x) - (72.0 * P2x * P2y) - (72.0 * P2x * P2z) + (9.0 * P4y) + (18.0 * P2y * P2z) + (9.0 * P4z);
    double term_2 = (24.0 * P4y) - (72.0 * P2x * P2y) - (72.0 * P2y * P2z) + (9.0 * P4x) + (18.0 * P2x * P2z) + (9.0 * P4z);
    double term_3 = (24.0 * P4z) - (72.0 * P2x * P2z) - (72.0 * P2y * P2z) + (9.0 * P4x) + (18.0 * P2x * P2y) + (9.0 * P4y);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP4x_exp = dP4_exp(P4x, C2, lambda, x0, gamma);
    double dP4y_exp = dP4_exp(P4y, C2, lambda, y0, gamma);
    double dP4z_exp = dP4_exp(P4z, C2, lambda, z0, gamma);

    double dterm_1_dx = (24.0 * dP4x_exp) - (72.0 * dP2x_exp * P2y) - (72.0 * dP2x_exp * P2z) + (9.0 * dP0x_exp * P4y) + (18.0 * dP0x_exp * P2y * P2z) + (9.0 * dP0x_exp * P4z);
    double dterm_1_dy = (24.0 * P4x * dP0y_exp) - (72.0 * P2x * dP2y_exp) - (72.0 * P2x * dP0y_exp * P2z) + (9.0 * dP4y_exp) + (18.0 * dP2y_exp * P2z) + (9.0 * dP0y_exp * P4z);
    double dterm_1_dz = (24.0 * P4x * dP0z_exp) - (72.0 * P2x * P2y * dP0z_exp) - (72.0 * P2x * dP2z_exp) + (9.0 * P4y * dP0z_exp) + (18.0 * P2y * dP2z_exp) + (9.0 * dP4z_exp);

    double dterm_2_dx = (24.0 * dP0x_exp * P4y) - (72.0 * dP2x_exp * P2y) - (72.0 * dP0x_exp * P2y * P2z) + (9.0 * dP4x_exp) + (18.0 * dP2x_exp * P2z) + (9.0 * dP0x_exp * P4z);
    double dterm_2_dy = (24.0 * dP4y_exp) - (72.0 * P2x * dP2y_exp) - (72.0 * dP2y_exp * P2z) + (9.0 * P4x * dP0y_exp) + (18.0 * P2x * dP0y_exp * P2z) + (9.0 * dP0y_exp * P4z);
    double dterm_2_dz = (24.0 * P4y * dP0z_exp) - (72.0 * P2x * P2y * dP0z_exp) - (72.0 * P2y * dP2z_exp) + (9.0 * P4x * dP0z_exp) + (18.0 * P2x * dP2z_exp) + (9.0 * dP4z_exp);

    double dterm_3_dx = (24.0 * dP0x_exp * P4z) - (72.0 * dP2x_exp * P2z) - (72.0 * dP0x_exp * P2y * P2z) + (9.0 * dP4x_exp) + (18.0 * dP2x_exp * P2y) + (9.0 * dP0x_exp * P4y);
    double dterm_3_dy = (24.0 * dP0y_exp * P4z) - (72.0 * P2x * dP0y_exp * P2z) - (72.0 * dP2y_exp * P2z) + (9.0 * P4x * dP0y_exp) + (18.0 * P2x * dP2y_exp) + (9.0 * dP4y_exp);
    double dterm_3_dz = (24.0 * dP4z_exp) - (72.0 * P2x * dP2z_exp) - (72.0 * P2y * dP2z_exp) + (9.0 * P4x * dP0z_exp) + (18.0 * P2x * P2y * dP0z_exp) + (9.0 * P4y * dP0z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

}

void calc_solid_MCSH_4_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

    double P1x = P1(lambda, x0, gamma);
    double P1y = P1(lambda, y0, gamma);
    double P1z = P1(lambda, z0, gamma);

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double P3x = P3(lambda, x0, gamma);
    double P3y = P3(lambda, y0, gamma);
    double P3z = P3(lambda, z0, gamma);

    double term_1 = (60.0 * P3x * P1y) - (45.0 * P1x * P3y) - (45.0 * P1x * P1y * P2z);
    double term_2 = (60.0 * P1x * P3y) - (45.0 * P3x * P1y) - (45.0 * P1x * P1y * P2z);
    double term_3 = (60.0 * P3x * P1z) - (45.0 * P1x * P3z) - (45.0 * P1x * P2y * P1z);
    double term_4 = (60.0 * P1x * P3z) - (45.0 * P3x * P1z) - (45.0 * P1x * P2y * P1z);
    double term_5 = (60.0 * P3y * P1z) - (45.0 * P1y * P3z) - (45.0 * P2x * P1y * P1z);
    double term_6 = (60.0 * P1y * P3z) - (45.0 * P3y * P1z) - (45.0 * P2x * P1y * P1z);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP3x_exp = dP3_exp(P3x, C2, lambda, x0, gamma);
    double dP3y_exp = dP3_exp(P3y, C2, lambda, y0, gamma);
    double dP3z_exp = dP3_exp(P3z, C2, lambda, z0, gamma);

    double dterm_1_dx = (60.0 * dP3x_exp * P1y) - (45.0 * dP1x_exp * P3y) - (45.0 * dP1x_exp * P1y * P2z);
    double dterm_1_dy = (60.0 * P3x * dP1y_exp) - (45.0 * P1x * dP3y_exp) - (45.0 * P1x * dP1y_exp * P2z);
    double dterm_1_dz = (60.0 * P3x * P1y * dP0z_exp) - (45.0 * P1x * P3y * dP0z_exp) - (45.0 * P1x * P1y * dP2z_exp);

    double dterm_2_dx = (60.0 * dP1x_exp * P3y) - (45.0 * dP3x_exp * P1y) - (45.0 * dP1x_exp * P1y * P2z);
    double dterm_2_dy = (60.0 * P1x * dP3y_exp) - (45.0 * P3x * dP1y_exp) - (45.0 * P1x * dP1y_exp * P2z);
    double dterm_2_dz = (60.0 * P1x * P3y * dP0z_exp) - (45.0 * P3x * P1y * dP0z_exp) - (45.0 * P1x * P1y * dP2z_exp);

    double dterm_3_dx = (60.0 * dP3x_exp * P1z) - (45.0 * dP1x_exp * P3z) - (45.0 * dP1x_exp * P2y * P1z);
    double dterm_3_dy = (60.0 * P3x * dP0y_exp * P1z) - (45.0 * P1x * dP0y_exp * P3z) - (45.0 * P1x * dP2y_exp * P1z);
    double dterm_3_dz = (60.0 * P3x * dP1z_exp) - (45.0 * P1x * dP3z_exp) - (45.0 * P1x * P2y * dP1z_exp);

    double dterm_4_dx = (60.0 * dP1x_exp * P3z) - (45.0 * dP3x_exp * P1z) - (45.0 * dP1x_exp * P2y * P1z);
    double dterm_4_dy = (60.0 * P1x * dP0y_exp * P3z) - (45.0 * P3x * dP0y_exp * P1z) - (45.0 * P1x * dP2y_exp * P1z);
    double dterm_4_dz = (60.0 * P1x * dP3z_exp) - (45.0 * P3x * dP1z_exp) - (45.0 * P1x * P2y * dP1z_exp);

    double dterm_5_dx = (60.0 * dP0x_exp * P3y * P1z) - (45.0 * dP0x_exp * P1y * P3z) - (45.0 * dP2x_exp * P1y * P1z);
    double dterm_5_dy = (60.0 * dP3y_exp * P1z) - (45.0 * dP1y_exp * P3z) - (45.0 * P2x * dP1y_exp * P1z);
    double dterm_5_dz = (60.0 * P3y * dP1z_exp) - (45.0 * P1y * dP3z_exp) - (45.0 * P2x * P1y * dP1z_exp);

    double dterm_6_dx = (60.0 * dP0x_exp * P1y * P3z) - (45.0 * dP0x_exp * P3y * P1z) - (45.0 * dP2x_exp * P1y * P1z);
    double dterm_6_dy = (60.0 * dP1y_exp * P3z) - (45.0 * dP3y_exp * P1z) - (45.0 * P2x * dP1y_exp * P1z);
    double dterm_6_dz = (60.0 * P1y * dP3z_exp) - (45.0 * P3y * dP1z_exp) - (45.0 * P2x * P1y * dP1z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

    // dmiu4 dx/dy/dz
    deriv[9] = temp * dterm_4_dx;
    deriv[10] = temp * dterm_4_dy;
    deriv[11] = temp * dterm_4_dz;

    // dmiu5 dx/dy/dz
    deriv[12] = temp * dterm_5_dx;
    deriv[13] = temp * dterm_5_dy;
    deriv[14] = temp * dterm_5_dz;

    // dmiu6 dx/dy/dz
    deriv[15] = temp * dterm_6_dx;
    deriv[16] = temp * dterm_6_dy;
    deriv[17] = temp * dterm_6_dz;

}

void calc_solid_MCSH_4_3(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double P4x = P4(lambda, x0, gamma);
    double P4y = P4(lambda, y0, gamma);
    double P4z = P4(lambda, z0, gamma);

    double term_1 = (81.0 * P2x * P2y) - (12.0 * P4x) - (9.0 * P2x * P2z) - (12.0 * P4y) - (9.0 * P2y * P2z) + (3.0 * P4z);
    double term_2 = (81.0 * P2x * P2z) - (12.0 * P4x) - (9.0 * P2x * P2y) - (12.0 * P4z) - (9.0 * P2y * P2z) + (3.0 * P4y);
    double term_3 = (81.0 * P2y * P2z) - (12.0 * P4y) - (9.0 * P2x * P2y) - (12.0 * P4z) - (9.0 * P2x * P2z) + (3.0 * P4x);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP4x_exp = dP4_exp(P4x, C2, lambda, x0, gamma);
    double dP4y_exp = dP4_exp(P4y, C2, lambda, y0, gamma);
    double dP4z_exp = dP4_exp(P4z, C2, lambda, z0, gamma);

    double dterm_1_dx = (81.0 * dP2x_exp * P2y) - (12.0 * dP4x_exp) - (9.0 * dP2x_exp * P2z) - (12.0 * dP0x_exp * P4y) - (9.0 * dP0x_exp * P2y * P2z) + (3.0 * dP0x_exp * P4z);
    double dterm_1_dy = (81.0 * P2x * dP2y_exp) - (12.0 * P4x * dP0y_exp) - (9.0 * P2x * dP0y_exp * P2z) - (12.0 * dP4y_exp) - (9.0 * dP2y_exp * P2z) + (3.0 * dP0y_exp * P4z);
    double dterm_1_dz = (81.0 * P2x * P2y * dP0z_exp) - (12.0 * P4x * dP0z_exp) - (9.0 * P2x * dP2z_exp) - (12.0 * P4y * dP0z_exp) - (9.0 * P2y * dP2z_exp) + (3.0 * dP4z_exp);

    double dterm_2_dx = (81.0 * dP2x_exp * P2z) - (12.0 * dP4x_exp) - (9.0 * dP2x_exp * P2y) - (12.0 * dP0x_exp * P4z) - (9.0 * dP0x_exp * P2y * P2z) + (3.0 * dP0x_exp * P4y);
    double dterm_2_dy = (81.0 * P2x * dP0y_exp * P2z) - (12.0 * P4x * dP0y_exp) - (9.0 * P2x * dP2y_exp) - (12.0 * dP0y_exp * P4z) - (9.0 * dP2y_exp * P2z) + (3.0 * dP4y_exp);
    double dterm_2_dz = (81.0 * P2x * dP2z_exp) - (12.0 * P4x * dP0z_exp) - (9.0 * P2x * P2y * dP0z_exp) - (12.0 * dP4z_exp) - (9.0 * P2y * dP2z_exp) + (3.0 * P4y * dP0z_exp);

    double dterm_3_dx = (81.0 * dP0x_exp * P2y * P2z) - (12.0 * dP0x_exp * P4y) - (9.0 * dP2x_exp * P2y) - (12.0 * dP0x_exp * P4z) - (9.0 * dP2x_exp * P2z) + (3.0 * dP4x_exp);
    double dterm_3_dy = (81.0 * dP2y_exp * P2z) - (12.0 * dP4y_exp) - (9.0 * P2x * dP2y_exp) - (12.0 * dP0y_exp * P4z) - (9.0 * P2x * dP0y_exp * P2z) + (3.0 * P4x * dP0y_exp);
    double dterm_3_dz = (81.0 * P2y * dP2z_exp) - (12.0 * P4y * dP0z_exp) - (9.0 * P2x * P2y * dP0z_exp) - (12.0 * dP4z_exp) - (9.0 * P2x * dP2z_exp) + (3.0 * P4x * dP0z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

}

void calc_solid_MCSH_4_4(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

    double P1x = P1(lambda, x0, gamma);
    double P1y = P1(lambda, y0, gamma);
    double P1z = P1(lambda, z0, gamma);

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double P3x = P3(lambda, x0, gamma);
    double P3y = P3(lambda, y0, gamma);
    double P3z = P3(lambda, z0, gamma);

    double term_1 = (90.0 * P2x * P1y * P1z) - (15.0 * P3y * P1z) - (15.0 * P1y * P3z);
    double term_2 = (90.0 * P1x * P2y * P1z) - (15.0 * P3x * P1z) - (15.0 * P1x * P3z);
    double term_3 = (90.0 * P1x * P1y * P2z) - (15.0 * P3x * P1y) - (15.0 * P1x * P3y);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP3x_exp = dP3_exp(P3x, C2, lambda, x0, gamma);
    double dP3y_exp = dP3_exp(P3y, C2, lambda, y0, gamma);
    double dP3z_exp = dP3_exp(P3z, C2, lambda, z0, gamma);

    double dterm_1_dx = (90.0 * dP2x_exp * P1y * P1z) - (15.0 * dP0x_exp * P3y * P1z) - (15.0 * dP0x_exp * P1y * P3z);
    double dterm_1_dy = (90.0 * P2x * dP1y_exp * P1z) - (15.0 * dP3y_exp * P1z) - (15.0 * dP1y_exp * P3z);
    double dterm_1_dz = (90.0 * P2x * P1y * dP1z_exp) - (15.0 * P3y * dP1z_exp) - (15.0 * P1y * dP3z_exp);

    double dterm_2_dx = (90.0 * dP1x_exp * P2y * P1z) - (15.0 * dP3x_exp * P1z) - (15.0 * dP1x_exp * P3z);
    double dterm_2_dy = (90.0 * P1x * dP2y_exp * P1z) - (15.0 * P3x * dP0y_exp * P1z) - (15.0 * P1x * dP0y_exp * P3z);
    double dterm_2_dz = (90.0 * P1x * P2y * dP1z_exp) - (15.0 * P3x * dP1z_exp) - (15.0 * P1x * dP3z_exp);

    double dterm_3_dx = (90.0 * dP1x_exp * P1y * P2z) - (15.0 * dP3x_exp * P1y) - (15.0 * dP1x_exp * P3y);
    double dterm_3_dy = (90.0 * P1x * dP1y_exp * P2z) - (15.0 * P3x * dP1y_exp) - (15.0 * P1x * dP3y_exp);
    double dterm_3_dz = (90.0 * P1x * P1y * dP2z_exp) - (15.0 * P3x * P1y * dP0z_exp) - (15.0 * P1x * P3y * dP0z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

}


void calc_solid_MCSH_5_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double term_1 = (120.0 * P5x) - (600.0 * P3x * P2y) - (600.0 * P3x * P2z) + (225.0 * P1x * P4y) + (450.0 * P1x * P2y * P2z) + (225.0 * P1x * P4z);
    double term_2 = (120.0 * P5y) - (600.0 * P2x * P3y) - (600.0 * P3y * P2z) + (225.0 * P4x * P1y) + (450.0 * P2x * P1y * P2z) + (225.0 * P1y * P4z);
    double term_3 = (120.0 * P5z) - (600.0 * P2x * P3z) - (600.0 * P2y * P3z) + (225.0 * P4x * P1z) + (450.0 * P2x * P2y * P1z) + (225.0 * P4y * P1z);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP3x_exp = dP3_exp(P3x, C2, lambda, x0, gamma);
    double dP3y_exp = dP3_exp(P3y, C2, lambda, y0, gamma);
    double dP3z_exp = dP3_exp(P3z, C2, lambda, z0, gamma);

    double dP4x_exp = dP4_exp(P4x, C2, lambda, x0, gamma);
    double dP4y_exp = dP4_exp(P4y, C2, lambda, y0, gamma);
    double dP4z_exp = dP4_exp(P4z, C2, lambda, z0, gamma);

    double dP5x_exp = dP5_exp(P5x, C2, lambda, x0, gamma);
    double dP5y_exp = dP5_exp(P5y, C2, lambda, y0, gamma);
    double dP5z_exp = dP5_exp(P5z, C2, lambda, z0, gamma);

    double dterm_1_dx = (120.0 * dP5x_exp) - (600.0 * dP3x_exp * P2y) - (600.0 * dP3x_exp * P2z) + (225.0 * dP1x_exp * P4y) + (450.0 * dP1x_exp * P2y * P2z) + (225.0 * dP1x_exp * P4z);
    double dterm_1_dy = (120.0 * P5x * dP0y_exp) - (600.0 * P3x * dP2y_exp) - (600.0 * P3x * dP0y_exp * P2z) + (225.0 * P1x * dP4y_exp) + (450.0 * P1x * dP2y_exp * P2z) + (225.0 * P1x * dP0y_exp * P4z);
    double dterm_1_dz = (120.0 * P5x * dP0z_exp) - (600.0 * P3x * P2y * dP0z_exp) - (600.0 * P3x * dP2z_exp) + (225.0 * P1x * P4y * dP0z_exp) + (450.0 * P1x * P2y * dP2z_exp) + (225.0 * P1x * dP4z_exp);

    double dterm_2_dx = (120.0 * dP0x_exp * P5y) - (600.0 * dP2x_exp * P3y) - (600.0 * dP0x_exp * P3y * P2z) + (225.0 * dP4x_exp * P1y) + (450.0 * dP2x_exp * P1y * P2z) + (225.0 * dP0x_exp * P1y * P4z);
    double dterm_2_dy = (120.0 * dP5y_exp) - (600.0 * P2x * dP3y_exp) - (600.0 * dP3y_exp * P2z) + (225.0 * P4x * dP1y_exp) + (450.0 * P2x * dP1y_exp * P2z) + (225.0 * dP1y_exp * P4z);
    double dterm_2_dz = (120.0 * P5y * dP0z_exp) - (600.0 * P2x * P3y * dP0z_exp) - (600.0 * P3y * dP2z_exp) + (225.0 * P4x * P1y * dP0z_exp) + (450.0 * P2x * P1y * dP2z_exp) + (225.0 * P1y * dP4z_exp);

    double dterm_3_dx = (120.0 * dP0x_exp * P5z) - (600.0 * dP2x_exp * P3z) - (600.0 * dP0x_exp * P2y * P3z) + (225.0 * dP4x_exp * P1z) + (450.0 * dP2x_exp * P2y * P1z) + (225.0 * dP0x_exp * P4y * P1z);
    double dterm_3_dy = (120.0 * dP0y_exp * P5z) - (600.0 * P2x * dP0y_exp * P3z) - (600.0 * dP2y_exp * P3z) + (225.0 * P4x * dP0y_exp * P1z) + (450.0 * P2x * dP2y_exp * P1z) + (225.0 * dP4y_exp * P1z);
    double dterm_3_dz = (120.0 * dP5z_exp) - (600.0 * P2x * dP3z_exp) - (600.0 * P2y * dP3z_exp) + (225.0 * P4x * dP1z_exp) + (450.0 * P2x * P2y * dP1z_exp) + (225.0 * P4y * dP1z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

}

void calc_solid_MCSH_5_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double term_1 = (360.0 * P4x * P1y) - (540.0 * P2x * P3y) - (540.0 * P2x * P1y * P2z) + (45.0 * P5y) + (90.0 * P3y * P2z) + (45.0 * P1y * P4z);
    double term_2 = (360.0 * P1x * P4y) - (540.0 * P3x * P2y) - (540.0 * P1x * P2y * P2z) + (45.0 * P5x) + (90.0 * P3x * P2z) + (45.0 * P1x * P4z);
    double term_3 = (360.0 * P4x * P1z) - (540.0 * P2x * P3z) - (540.0 * P2x * P2y * P1z) + (45.0 * P5z) + (90.0 * P2y * P3z) + (45.0 * P4y * P1z);
    double term_4 = (360.0 * P1x * P4z) - (540.0 * P3x * P2z) - (540.0 * P1x * P2y * P2z) + (45.0 * P5x) + (90.0 * P3x * P2y) + (45.0 * P1x * P4y);
    double term_5 = (360.0 * P4y * P1z) - (540.0 * P2y * P3z) - (540.0 * P2x * P2y * P1z) + (45.0 * P5z) + (90.0 * P2x * P3z) + (45.0 * P4x * P1z);
    double term_6 = (360.0 * P1y * P4z) - (540.0 * P3y * P2z) - (540.0 * P2x * P1y * P2z) + (45.0 * P5y) + (90.0 * P2x * P3y) + (45.0 * P4x * P1y);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP3x_exp = dP3_exp(P3x, C2, lambda, x0, gamma);
    double dP3y_exp = dP3_exp(P3y, C2, lambda, y0, gamma);
    double dP3z_exp = dP3_exp(P3z, C2, lambda, z0, gamma);

    double dP4x_exp = dP4_exp(P4x, C2, lambda, x0, gamma);
    double dP4y_exp = dP4_exp(P4y, C2, lambda, y0, gamma);
    double dP4z_exp = dP4_exp(P4z, C2, lambda, z0, gamma);

    double dP5x_exp = dP5_exp(P5x, C2, lambda, x0, gamma);
    double dP5y_exp = dP5_exp(P5y, C2, lambda, y0, gamma);
    double dP5z_exp = dP5_exp(P5z, C2, lambda, z0, gamma);

    double dterm_1_dx = (360.0 * dP4x_exp * P1y) - (540.0 * dP2x_exp * P3y) - (540.0 * dP2x_exp * P1y * P2z) + (45.0 * dP0x_exp * P5y) + (90.0 * dP0x_exp * P3y * P2z) + (45.0 * dP0x_exp * P1y * P4z);
    double dterm_1_dy = (360.0 * P4x * dP1y_exp) - (540.0 * P2x * dP3y_exp) - (540.0 * P2x * dP1y_exp * P2z) + (45.0 * dP5y_exp) + (90.0 * dP3y_exp * P2z) + (45.0 * dP1y_exp * P4z);
    double dterm_1_dz = (360.0 * P4x * P1y * dP0z_exp) - (540.0 * P2x * P3y * dP0z_exp) - (540.0 * P2x * P1y * dP2z_exp) + (45.0 * P5y * dP0z_exp) + (90.0 * P3y * dP2z_exp) + (45.0 * P1y * dP4z_exp);

    double dterm_2_dx = (360.0 * dP1x_exp * P4y) - (540.0 * dP3x_exp * P2y) - (540.0 * dP1x_exp * P2y * P2z) + (45.0 * dP5x_exp) + (90.0 * dP3x_exp * P2z) + (45.0 * dP1x_exp * P4z);
    double dterm_2_dy = (360.0 * P1x * dP4y_exp) - (540.0 * P3x * dP2y_exp) - (540.0 * P1x * dP2y_exp * P2z) + (45.0 * P5x * dP0y_exp) + (90.0 * P3x * dP0y_exp * P2z) + (45.0 * P1x * dP0y_exp * P4z);
    double dterm_2_dz = (360.0 * P1x * P4y * dP0z_exp) - (540.0 * P3x * P2y * dP0z_exp) - (540.0 * P1x * P2y * dP2z_exp) + (45.0 * P5x * dP0z_exp) + (90.0 * P3x * dP2z_exp) + (45.0 * P1x * dP4z_exp);

    double dterm_3_dx = (360.0 * dP4x_exp * P1z) - (540.0 * dP2x_exp * P3z) - (540.0 * dP2x_exp * P2y * P1z) + (45.0 * dP0x_exp * P5z) + (90.0 * dP0x_exp * P2y * P3z) + (45.0 * dP0x_exp * P4y * P1z);
    double dterm_3_dy = (360.0 * P4x * dP0y_exp * P1z) - (540.0 * P2x * dP0y_exp * P3z) - (540.0 * P2x * dP2y_exp * P1z) + (45.0 * dP0y_exp * P5z) + (90.0 * dP2y_exp * P3z) + (45.0 * dP4y_exp * P1z);
    double dterm_3_dz = (360.0 * P4x * dP1z_exp) - (540.0 * P2x * dP3z_exp) - (540.0 * P2x * P2y * dP1z_exp) + (45.0 * dP5z_exp) + (90.0 * P2y * dP3z_exp) + (45.0 * P4y * dP1z_exp);

    double dterm_4_dx = (360.0 * dP1x_exp * P4z) - (540.0 * dP3x_exp * P2z) - (540.0 * dP1x_exp * P2y * P2z) + (45.0 * dP5x_exp) + (90.0 * dP3x_exp * P2y) + (45.0 * dP1x_exp * P4y);
    double dterm_4_dy = (360.0 * P1x * dP0y_exp * P4z) - (540.0 * P3x * dP0y_exp * P2z) - (540.0 * P1x * dP2y_exp * P2z) + (45.0 * P5x * dP0y_exp) + (90.0 * P3x * dP2y_exp) + (45.0 * P1x * dP4y_exp);
    double dterm_4_dz = (360.0 * P1x * dP4z_exp) - (540.0 * P3x * dP2z_exp) - (540.0 * P1x * P2y * dP2z_exp) + (45.0 * P5x * dP0z_exp) + (90.0 * P3x * P2y * dP0z_exp) + (45.0 * P1x * P4y * dP0z_exp);

    double dterm_5_dx = (360.0 * dP0x_exp * P4y * P1z) - (540.0 * dP0x_exp * P2y * P3z) - (540.0 * dP2x_exp * P2y * P1z) + (45.0 * dP0x_exp * P5z) + (90.0 * dP2x_exp * P3z) + (45.0 * dP4x_exp * P1z);
    double dterm_5_dy = (360.0 * dP4y_exp * P1z) - (540.0 * dP2y_exp * P3z) - (540.0 * P2x * dP2y_exp * P1z) + (45.0 * dP0y_exp * P5z) + (90.0 * P2x * dP0y_exp * P3z) + (45.0 * P4x * dP0y_exp * P1z);
    double dterm_5_dz = (360.0 * P4y * dP1z_exp) - (540.0 * P2y * dP3z_exp) - (540.0 * P2x * P2y * dP1z_exp) + (45.0 * dP5z_exp) + (90.0 * P2x * dP3z_exp) + (45.0 * P4x * dP1z_exp);

    double dterm_6_dx = (360.0 * dP0x_exp * P1y * P4z) - (540.0 * dP0x_exp * P3y * P2z) - (540.0 * dP2x_exp * P1y * P2z) + (45.0 * dP0x_exp * P5y) + (90.0 * dP2x_exp * P3y) + (45.0 * dP4x_exp * P1y);
    double dterm_6_dy = (360.0 * dP1y_exp * P4z) - (540.0 * dP3y_exp * P2z) - (540.0 * P2x * dP1y_exp * P2z) + (45.0 * dP5y_exp) + (90.0 * P2x * dP3y_exp) + (45.0 * P4x * dP1y_exp);
    double dterm_6_dz = (360.0 * P1y * dP4z_exp) - (540.0 * P3y * dP2z_exp) - (540.0 * P2x * P1y * dP2z_exp) + (45.0 * P5y * dP0z_exp) + (90.0 * P2x * P3y * dP0z_exp) + (45.0 * P4x * P1y * dP0z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

    // dmiu4 dx/dy/dz
    deriv[9] = temp * dterm_4_dx;
    deriv[10] = temp * dterm_4_dy;
    deriv[11] = temp * dterm_4_dz;

    // dmiu5 dx/dy/dz
    deriv[12] = temp * dterm_5_dx;
    deriv[13] = temp * dterm_5_dy;
    deriv[14] = temp * dterm_5_dz;

    // dmiu6 dx/dy/dz
    deriv[15] = temp * dterm_6_dx;
    deriv[16] = temp * dterm_6_dy;
    deriv[17] = temp * dterm_6_dz;

}

void calc_solid_MCSH_5_3(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double term_1 = (615.0 * P3x * P2y) - (60.0 * P5x) - (15.0 * P3x * P2z) - (270.0 * P1x * P4y) - (225.0 * P1x * P2y * P2z) + (45.0 * P1x * P4z);
    double term_2 = (615.0 * P2x * P3y) - (60.0 * P5y) - (15.0 * P3y * P2z) - (270.0 * P4x * P1y) - (225.0 * P2x * P1y * P2z) + (45.0 * P1y * P4z);
    double term_3 = (615.0 * P3x * P2z) - (60.0 * P5x) - (15.0 * P3x * P2y) - (270.0 * P1x * P4z) - (225.0 * P1x * P2y * P2z) + (45.0 * P1x * P4y);
    double term_4 = (615.0 * P2x * P3z) - (60.0 * P5z) - (15.0 * P2y * P3z) - (270.0 * P4x * P1z) - (225.0 * P2x * P2y * P1z) + (45.0 * P4y * P1z);
    double term_5 = (615.0 * P3y * P2z) - (60.0 * P5y) - (15.0 * P2x * P3y) - (270.0 * P1y * P4z) - (225.0 * P2x * P1y * P2z) + (45.0 * P4x * P1y);
    double term_6 = (615.0 * P2y * P3z) - (60.0 * P5z) - (15.0 * P2x * P3z) - (270.0 * P4y * P1z) - (225.0 * P2x * P2y * P1z) + (45.0 * P4x * P1z);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP3x_exp = dP3_exp(P3x, C2, lambda, x0, gamma);
    double dP3y_exp = dP3_exp(P3y, C2, lambda, y0, gamma);
    double dP3z_exp = dP3_exp(P3z, C2, lambda, z0, gamma);

    double dP4x_exp = dP4_exp(P4x, C2, lambda, x0, gamma);
    double dP4y_exp = dP4_exp(P4y, C2, lambda, y0, gamma);
    double dP4z_exp = dP4_exp(P4z, C2, lambda, z0, gamma);

    double dP5x_exp = dP5_exp(P5x, C2, lambda, x0, gamma);
    double dP5y_exp = dP5_exp(P5y, C2, lambda, y0, gamma);
    double dP5z_exp = dP5_exp(P5z, C2, lambda, z0, gamma);

    double dterm_1_dx = (615.0 * dP3x_exp * P2y) - (60.0 * dP5x_exp) - (15.0 * dP3x_exp * P2z) - (270.0 * dP1x_exp * P4y) - (225.0 * dP1x_exp * P2y * P2z) + (45.0 * dP1x_exp * P4z);
    double dterm_1_dy = (615.0 * P3x * dP2y_exp) - (60.0 * P5x * dP0y_exp) - (15.0 * P3x * dP0y_exp * P2z) - (270.0 * P1x * dP4y_exp) - (225.0 * P1x * dP2y_exp * P2z) + (45.0 * P1x * dP0y_exp * P4z);
    double dterm_1_dz = (615.0 * P3x * P2y * dP0z_exp) - (60.0 * P5x * dP0z_exp) - (15.0 * P3x * dP2z_exp) - (270.0 * P1x * P4y * dP0z_exp) - (225.0 * P1x * P2y * dP2z_exp) + (45.0 * P1x * dP4z_exp);

    double dterm_2_dx = (615.0 * dP2x_exp * P3y) - (60.0 * dP0x_exp * P5y) - (15.0 * dP0x_exp * P3y * P2z) - (270.0 * dP4x_exp * P1y) - (225.0 * dP2x_exp * P1y * P2z) + (45.0 * dP0x_exp * P1y * P4z);
    double dterm_2_dy = (615.0 * P2x * dP3y_exp) - (60.0 * dP5y_exp) - (15.0 * dP3y_exp * P2z) - (270.0 * P4x * dP1y_exp) - (225.0 * P2x * dP1y_exp * P2z) + (45.0 * dP1y_exp * P4z);
    double dterm_2_dz = (615.0 * P2x * P3y * dP0z_exp) - (60.0 * P5y * dP0z_exp) - (15.0 * P3y * dP2z_exp) - (270.0 * P4x * P1y * dP0z_exp) - (225.0 * P2x * P1y * dP2z_exp) + (45.0 * P1y * dP4z_exp);

    double dterm_3_dx = (615.0 * dP3x_exp * P2z) - (60.0 * dP5x_exp) - (15.0 * dP3x_exp * P2y) - (270.0 * dP1x_exp * P4z) - (225.0 * dP1x_exp * P2y * P2z) + (45.0 * dP1x_exp * P4y);
    double dterm_3_dy = (615.0 * P3x * dP0y_exp * P2z) - (60.0 * P5x * dP0y_exp) - (15.0 * P3x * dP2y_exp) - (270.0 * P1x * dP0y_exp * P4z) - (225.0 * P1x * dP2y_exp * P2z) + (45.0 * P1x * dP4y_exp);
    double dterm_3_dz = (615.0 * P3x * dP2z_exp) - (60.0 * P5x * dP0z_exp) - (15.0 * P3x * P2y * dP0z_exp) - (270.0 * P1x * dP4z_exp) - (225.0 * P1x * P2y * dP2z_exp) + (45.0 * P1x * P4y * dP0z_exp);

    double dterm_4_dx = (615.0 * dP2x_exp * P3z) - (60.0 * dP0x_exp * P5z) - (15.0 * dP0x_exp * P2y * P3z) - (270.0 * dP4x_exp * P1z) - (225.0 * dP2x_exp * P2y * P1z) + (45.0 * dP0x_exp * P4y * P1z);
    double dterm_4_dy = (615.0 * P2x * dP0y_exp * P3z) - (60.0 * dP0y_exp * P5z) - (15.0 * dP2y_exp * P3z) - (270.0 * P4x * dP0y_exp * P1z) - (225.0 * P2x * dP2y_exp * P1z) + (45.0 * dP4y_exp * P1z);
    double dterm_4_dz = (615.0 * P2x * dP3z_exp) - (60.0 * dP5z_exp) - (15.0 * P2y * dP3z_exp) - (270.0 * P4x * dP1z_exp) - (225.0 * P2x * P2y * dP1z_exp) + (45.0 * P4y * dP1z_exp);

    double dterm_5_dx = (615.0 * dP0x_exp * P3y * P2z) - (60.0 * dP0x_exp * P5y) - (15.0 * dP2x_exp * P3y) - (270.0 * dP0x_exp * P1y * P4z) - (225.0 * dP2x_exp * P1y * P2z) + (45.0 * dP4x_exp * P1y);
    double dterm_5_dy = (615.0 * dP3y_exp * P2z) - (60.0 * dP5y_exp) - (15.0 * P2x * dP3y_exp) - (270.0 * dP1y_exp * P4z) - (225.0 * P2x * dP1y_exp * P2z) + (45.0 * P4x * dP1y_exp);
    double dterm_5_dz = (615.0 * P3y * dP2z_exp) - (60.0 * P5y * dP0z_exp) - (15.0 * P2x * P3y * dP0z_exp) - (270.0 * P1y * dP4z_exp) - (225.0 * P2x * P1y * dP2z_exp) + (45.0 * P4x * P1y * dP0z_exp);

    double dterm_6_dx = (615.0 * dP0x_exp * P2y * P3z) - (60.0 * dP0x_exp * P5z) - (15.0 * dP2x_exp * P3z) - (270.0 * dP0x_exp * P4y * P1z) - (225.0 * dP2x_exp * P2y * P1z) + (45.0 * dP4x_exp * P1z);
    double dterm_6_dy = (615.0 * dP2y_exp * P3z) - (60.0 * dP0y_exp * P5z) - (15.0 * P2x * dP0y_exp * P3z) - (270.0 * dP4y_exp * P1z) - (225.0 * P2x * dP2y_exp * P1z) + (45.0 * P4x * dP0y_exp * P1z);
    double dterm_6_dz = (615.0 * P2y * dP3z_exp) - (60.0 * dP5z_exp) - (15.0 * P2x * dP3z_exp) - (270.0 * P4y * dP1z_exp) - (225.0 * P2x * P2y * dP1z_exp) + (45.0 * P4x * dP1z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

    // dmiu4 dx/dy/dz
    deriv[9] = temp * dterm_4_dx;
    deriv[10] = temp * dterm_4_dy;
    deriv[11] = temp * dterm_4_dz;

    // dmiu5 dx/dy/dz
    deriv[12] = temp * dterm_5_dx;
    deriv[13] = temp * dterm_5_dy;
    deriv[14] = temp * dterm_5_dz;

    // dmiu6 dx/dy/dz
    deriv[15] = temp * dterm_6_dx;
    deriv[16] = temp * dterm_6_dy;
    deriv[17] = temp * dterm_6_dz;

}

void calc_solid_MCSH_5_4(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

    double P1x = P1(lambda, x0, gamma);
    double P1y = P1(lambda, y0, gamma);
    double P1z = P1(lambda, z0, gamma);

    double P3x = P3(lambda, x0, gamma);
    double P3y = P3(lambda, y0, gamma);
    double P3z = P3(lambda, z0, gamma);

    double term_1 = (630.0 * P3x * P1y * P1z) - (315.0 * P1x * P3y * P1z) - (315.0 * P1x * P1y * P3z);
    double term_2 = (630.0 * P1x * P3y * P1z) - (315.0 * P3x * P1y * P1z) - (315.0 * P1x * P1y * P3z);
    double term_3 = (630.0 * P1x * P1y * P3z) - (315.0 * P3x * P1y * P1z) - (315.0 * P1x * P3y * P1z);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dP3x_exp = dP3_exp(P3x, C2, lambda, x0, gamma);
    double dP3y_exp = dP3_exp(P3y, C2, lambda, y0, gamma);
    double dP3z_exp = dP3_exp(P3z, C2, lambda, z0, gamma);

    double dterm_1_dx = (630.0 * dP3x_exp * P1y * P1z) - (315.0 * dP1x_exp * P3y * P1z) - (315.0 * dP1x_exp * P1y * P3z);
    double dterm_1_dy = (630.0 * P3x * dP1y_exp * P1z) - (315.0 * P1x * dP3y_exp * P1z) - (315.0 * P1x * dP1y_exp * P3z);
    double dterm_1_dz = (630.0 * P3x * P1y * dP1z_exp) - (315.0 * P1x * P3y * dP1z_exp) - (315.0 * P1x * P1y * dP3z_exp);

    double dterm_2_dx = (630.0 * dP1x_exp * P3y * P1z) - (315.0 * dP3x_exp * P1y * P1z) - (315.0 * dP1x_exp * P1y * P3z);
    double dterm_2_dy = (630.0 * P1x * dP3y_exp * P1z) - (315.0 * P3x * dP1y_exp * P1z) - (315.0 * P1x * dP1y_exp * P3z);
    double dterm_2_dz = (630.0 * P1x * P3y * dP1z_exp) - (315.0 * P3x * P1y * dP1z_exp) - (315.0 * P1x * P1y * dP3z_exp);

    double dterm_3_dx = (630.0 * dP1x_exp * P1y * P3z) - (315.0 * dP3x_exp * P1y * P1z) - (315.0 * dP1x_exp * P3y * P1z);
    double dterm_3_dy = (630.0 * P1x * dP1y_exp * P3z) - (315.0 * P3x * dP1y_exp * P1z) - (315.0 * P1x * dP3y_exp * P1z);
    double dterm_3_dz = (630.0 * P1x * P1y * dP3z_exp) - (315.0 * P3x * P1y * dP1z_exp) - (315.0 * P1x * P3y * dP1z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

}

void calc_solid_MCSH_5_5(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double term_1 = (765.0 * P2x * P2y * P1z) - (90.0 * P4x * P1z) - (75.0 * P2x * P3z) - (90.0 * P4y * P1z) - (75.0 * P2y * P3z) + (15.0 * P5z);
    double term_2 = (765.0 * P2x * P1y * P2z) - (90.0 * P4x * P1y) - (75.0 * P2x * P3y) - (90.0 * P1y * P4z) - (75.0 * P3y * P2z) + (15.0 * P5y);
    double term_3 = (765.0 * P1x * P2y * P2z) - (90.0 * P1x * P4y) - (75.0 * P3x * P2y) - (90.0 * P1x * P4z) - (75.0 * P3x * P2z) + (15.0 * P5x);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP3x_exp = dP3_exp(P3x, C2, lambda, x0, gamma);
    double dP3y_exp = dP3_exp(P3y, C2, lambda, y0, gamma);
    double dP3z_exp = dP3_exp(P3z, C2, lambda, z0, gamma);

    double dP4x_exp = dP4_exp(P4x, C2, lambda, x0, gamma);
    double dP4y_exp = dP4_exp(P4y, C2, lambda, y0, gamma);
    double dP4z_exp = dP4_exp(P4z, C2, lambda, z0, gamma);

    double dP5x_exp = dP5_exp(P5x, C2, lambda, x0, gamma);
    double dP5y_exp = dP5_exp(P5y, C2, lambda, y0, gamma);
    double dP5z_exp = dP5_exp(P5z, C2, lambda, z0, gamma);

    double dterm_1_dx = (765.0 * dP2x_exp * P2y * P1z) - (90.0 * dP4x_exp * P1z) - (75.0 * dP2x_exp * P3z) - (90.0 * dP0x_exp * P4y * P1z) - (75.0 * dP0x_exp * P2y * P3z) + (15.0 * dP0x_exp * P5z);
    double dterm_1_dy = (765.0 * P2x * dP2y_exp * P1z) - (90.0 * P4x * dP0y_exp * P1z) - (75.0 * P2x * dP0y_exp * P3z) - (90.0 * dP4y_exp * P1z) - (75.0 * dP2y_exp * P3z) + (15.0 * dP0y_exp * P5z);
    double dterm_1_dz = (765.0 * P2x * P2y * dP1z_exp) - (90.0 * P4x * dP1z_exp) - (75.0 * P2x * dP3z_exp) - (90.0 * P4y * dP1z_exp) - (75.0 * P2y * dP3z_exp) + (15.0 * dP5z_exp);

    double dterm_2_dx = (765.0 * dP2x_exp * P1y * P2z) - (90.0 * dP4x_exp * P1y) - (75.0 * dP2x_exp * P3y) - (90.0 * dP0x_exp * P1y * P4z) - (75.0 * dP0x_exp * P3y * P2z) + (15.0 * dP0x_exp * P5y);
    double dterm_2_dy = (765.0 * P2x * dP1y_exp * P2z) - (90.0 * P4x * dP1y_exp) - (75.0 * P2x * dP3y_exp) - (90.0 * dP1y_exp * P4z) - (75.0 * dP3y_exp * P2z) + (15.0 * dP5y_exp);
    double dterm_2_dz = (765.0 * P2x * P1y * dP2z_exp) - (90.0 * P4x * P1y * dP0z_exp) - (75.0 * P2x * P3y * dP0z_exp) - (90.0 * P1y * dP4z_exp) - (75.0 * P3y * dP2z_exp) + (15.0 * P5y * dP0z_exp);

    double dterm_3_dx = (765.0 * dP1x_exp * P2y * P2z) - (90.0 * dP1x_exp * P4y) - (75.0 * dP3x_exp * P2y) - (90.0 * dP1x_exp * P4z) - (75.0 * dP3x_exp * P2z) + (15.0 * dP5x_exp);
    double dterm_3_dy = (765.0 * P1x * dP2y_exp * P2z) - (90.0 * P1x * dP4y_exp) - (75.0 * P3x * dP2y_exp) - (90.0 * P1x * dP0y_exp * P4z) - (75.0 * P3x * dP0y_exp * P2z) + (15.0 * P5x * dP0y_exp);
    double dterm_3_dz = (765.0 * P1x * P2y * dP2z_exp) - (90.0 * P1x * P4y * dP0z_exp) - (75.0 * P3x * P2y * dP0z_exp) - (90.0 * P1x * dP4z_exp) - (75.0 * P3x * dP2z_exp) + (15.0 * P5x * dP0z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

}


void calc_solid_MCSH_6_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double P4x = P4(lambda, x0, gamma);
    double P4y = P4(lambda, y0, gamma);
    double P4z = P4(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double term_1 = (720.0 * P6x) - (5400.0 * P4x * P2y) - (5400.0 * P4x * P2z) + (4050.0 * P2x * P4y) + (8100.0 * P2x * P2y * P2z) + (4050.0 * P2x * P4z) - (225.0 * P6y) - (675.0 * P4y * P2z) - (675.0 * P2y * P4z) - (225.0 * P6z);
    double term_2 = (720.0 * P6y) - (5400.0 * P2x * P4y) - (5400.0 * P4y * P2z) + (4050.0 * P4x * P2y) + (8100.0 * P2x * P2y * P2z) + (4050.0 * P2y * P4z) - (225.0 * P6x) - (675.0 * P4x * P2z) - (675.0 * P2x * P4z) - (225.0 * P6z);
    double term_3 = (720.0 * P6z) - (5400.0 * P2x * P4z) - (5400.0 * P2y * P4z) + (4050.0 * P4x * P2z) + (8100.0 * P2x * P2y * P2z) + (4050.0 * P4y * P2z) - (225.0 * P6x) - (675.0 * P4x * P2y) - (675.0 * P2x * P4y) - (225.0 * P6y);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP4x_exp = dP4_exp(P4x, C2, lambda, x0, gamma);
    double dP4y_exp = dP4_exp(P4y, C2, lambda, y0, gamma);
    double dP4z_exp = dP4_exp(P4z, C2, lambda, z0, gamma);

    double dP6x_exp = dP6_exp(P6x, C2, lambda, x0, gamma);
    double dP6y_exp = dP6_exp(P6y, C2, lambda, y0, gamma);
    double dP6z_exp = dP6_exp(P6z, C2, lambda, z0, gamma);

    double dterm_1_dx = (720.0 * dP6x_exp) - (5400.0 * dP4x_exp * P2y) - (5400.0 * dP4x_exp * P2z) + (4050.0 * dP2x_exp * P4y) + (8100.0 * dP2x_exp * P2y * P2z) + (4050.0 * dP2x_exp * P4z) - (225.0 * dP0x_exp * P6y) - (675.0 * dP0x_exp * P4y * P2z) - (675.0 * dP0x_exp * P2y * P4z) - (225.0 * dP0x_exp * P6z);
    double dterm_1_dy = (720.0 * P6x * dP0y_exp) - (5400.0 * P4x * dP2y_exp) - (5400.0 * P4x * dP0y_exp * P2z) + (4050.0 * P2x * dP4y_exp) + (8100.0 * P2x * dP2y_exp * P2z) + (4050.0 * P2x * dP0y_exp * P4z) - (225.0 * dP6y_exp) - (675.0 * dP4y_exp * P2z) - (675.0 * dP2y_exp * P4z) - (225.0 * dP0y_exp * P6z);
    double dterm_1_dz = (720.0 * P6x * dP0z_exp) - (5400.0 * P4x * P2y * dP0z_exp) - (5400.0 * P4x * dP2z_exp) + (4050.0 * P2x * P4y * dP0z_exp) + (8100.0 * P2x * P2y * dP2z_exp) + (4050.0 * P2x * dP4z_exp) - (225.0 * P6y * dP0z_exp) - (675.0 * P4y * dP2z_exp) - (675.0 * P2y * dP4z_exp) - (225.0 * dP6z_exp);

    double dterm_2_dx = (720.0 * dP0x_exp * P6y) - (5400.0 * dP2x_exp * P4y) - (5400.0 * dP0x_exp * P4y * P2z) + (4050.0 * dP4x_exp * P2y) + (8100.0 * dP2x_exp * P2y * P2z) + (4050.0 * dP0x_exp * P2y * P4z) - (225.0 * dP6x_exp) - (675.0 * dP4x_exp * P2z) - (675.0 * dP2x_exp * P4z) - (225.0 * dP0x_exp * P6z);
    double dterm_2_dy = (720.0 * dP6y_exp) - (5400.0 * P2x * dP4y_exp) - (5400.0 * dP4y_exp * P2z) + (4050.0 * P4x * dP2y_exp) + (8100.0 * P2x * dP2y_exp * P2z) + (4050.0 * dP2y_exp * P4z) - (225.0 * P6x * dP0y_exp) - (675.0 * P4x * dP0y_exp * P2z) - (675.0 * P2x * dP0y_exp * P4z) - (225.0 * dP0y_exp * P6z);
    double dterm_2_dz = (720.0 * P6y * dP0z_exp) - (5400.0 * P2x * P4y * dP0z_exp) - (5400.0 * P4y * dP2z_exp) + (4050.0 * P4x * P2y * dP0z_exp) + (8100.0 * P2x * P2y * dP2z_exp) + (4050.0 * P2y * dP4z_exp) - (225.0 * P6x * dP0z_exp) - (675.0 * P4x * dP2z_exp) - (675.0 * P2x * dP4z_exp) - (225.0 * dP6z_exp);

    double dterm_3_dx = (720.0 * dP0x_exp * P6z) - (5400.0 * dP2x_exp * P4z) - (5400.0 * dP0x_exp * P2y * P4z) + (4050.0 * dP4x_exp * P2z) + (8100.0 * dP2x_exp * P2y * P2z) + (4050.0 * dP0x_exp * P4y * P2z) - (225.0 * dP6x_exp) - (675.0 * dP4x_exp * P2y) - (675.0 * dP2x_exp * P4y) - (225.0 * dP0x_exp * P6y);
    double dterm_3_dy = (720.0 * dP0y_exp * P6z) - (5400.0 * P2x * dP0y_exp * P4z) - (5400.0 * dP2y_exp * P4z) + (4050.0 * P4x * dP0y_exp * P2z) + (8100.0 * P2x * dP2y_exp * P2z) + (4050.0 * dP4y_exp * P2z) - (225.0 * P6x * dP0y_exp) - (675.0 * P4x * dP2y_exp) - (675.0 * P2x * dP4y_exp) - (225.0 * dP6y_exp);
    double dterm_3_dz = (720.0 * dP6z_exp) - (5400.0 * P2x * dP4z_exp) - (5400.0 * P2y * dP4z_exp) + (4050.0 * P4x * dP2z_exp) + (8100.0 * P2x * P2y * dP2z_exp) + (4050.0 * P4y * dP2z_exp) - (225.0 * P6x * dP0z_exp) - (675.0 * P4x * P2y * dP0z_exp) - (675.0 * P2x * P4y * dP0z_exp) - (225.0 * P6y * dP0z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

}

void calc_solid_MCSH_6_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double term_1 = (2520.0 * P5x * P1y) - (6300.0 * P3x * P3y) - (6300.0 * P3x * P1y * P2z) + (1575.0 * P1x * P5y) + (3150.0 * P1x * P3y * P2z) + (1575.0 * P1x * P1y * P4z);
    double term_2 = (2520.0 * P1x * P5y) - (6300.0 * P3x * P3y) - (6300.0 * P1x * P3y * P2z) + (1575.0 * P5x * P1y) + (3150.0 * P3x * P1y * P2z) + (1575.0 * P1x * P1y * P4z);
    double term_3 = (2520.0 * P5x * P1z) - (6300.0 * P3x * P3z) - (6300.0 * P3x * P2y * P1z) + (1575.0 * P1x * P5z) + (3150.0 * P1x * P2y * P3z) + (1575.0 * P1x * P4y * P1z);
    double term_4 = (2520.0 * P1x * P5z) - (6300.0 * P3x * P3z) - (6300.0 * P1x * P2y * P3z) + (1575.0 * P5x * P1z) + (3150.0 * P3x * P2y * P1z) + (1575.0 * P1x * P4y * P1z);
    double term_5 = (2520.0 * P5y * P1z) - (6300.0 * P3y * P3z) - (6300.0 * P2x * P3y * P1z) + (1575.0 * P1y * P5z) + (3150.0 * P2x * P1y * P3z) + (1575.0 * P4x * P1y * P1z);
    double term_6 = (2520.0 * P1y * P5z) - (6300.0 * P3y * P3z) - (6300.0 * P2x * P1y * P3z) + (1575.0 * P5y * P1z) + (3150.0 * P2x * P3y * P1z) + (1575.0 * P4x * P1y * P1z);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP3x_exp = dP3_exp(P3x, C2, lambda, x0, gamma);
    double dP3y_exp = dP3_exp(P3y, C2, lambda, y0, gamma);
    double dP3z_exp = dP3_exp(P3z, C2, lambda, z0, gamma);

    double dP4x_exp = dP4_exp(P4x, C2, lambda, x0, gamma);
    double dP4y_exp = dP4_exp(P4y, C2, lambda, y0, gamma);
    double dP4z_exp = dP4_exp(P4z, C2, lambda, z0, gamma);

    double dP5x_exp = dP5_exp(P5x, C2, lambda, x0, gamma);
    double dP5y_exp = dP5_exp(P5y, C2, lambda, y0, gamma);
    double dP5z_exp = dP5_exp(P5z, C2, lambda, z0, gamma);

    double dterm_1_dx = (2520.0 * dP5x_exp * P1y) - (6300.0 * dP3x_exp * P3y) - (6300.0 * dP3x_exp * P1y * P2z) + (1575.0 * dP1x_exp * P5y) + (3150.0 * dP1x_exp * P3y * P2z) + (1575.0 * dP1x_exp * P1y * P4z);
    double dterm_1_dy = (2520.0 * P5x * dP1y_exp) - (6300.0 * P3x * dP3y_exp) - (6300.0 * P3x * dP1y_exp * P2z) + (1575.0 * P1x * dP5y_exp) + (3150.0 * P1x * dP3y_exp * P2z) + (1575.0 * P1x * dP1y_exp * P4z);
    double dterm_1_dz = (2520.0 * P5x * P1y * dP0z_exp) - (6300.0 * P3x * P3y * dP0z_exp) - (6300.0 * P3x * P1y * dP2z_exp) + (1575.0 * P1x * P5y * dP0z_exp) + (3150.0 * P1x * P3y * dP2z_exp) + (1575.0 * P1x * P1y * dP4z_exp);

    double dterm_2_dx = (2520.0 * dP1x_exp * P5y) - (6300.0 * dP3x_exp * P3y) - (6300.0 * dP1x_exp * P3y * P2z) + (1575.0 * dP5x_exp * P1y) + (3150.0 * dP3x_exp * P1y * P2z) + (1575.0 * dP1x_exp * P1y * P4z);
    double dterm_2_dy = (2520.0 * P1x * dP5y_exp) - (6300.0 * P3x * dP3y_exp) - (6300.0 * P1x * dP3y_exp * P2z) + (1575.0 * P5x * dP1y_exp) + (3150.0 * P3x * dP1y_exp * P2z) + (1575.0 * P1x * dP1y_exp * P4z);
    double dterm_2_dz = (2520.0 * P1x * P5y * dP0z_exp) - (6300.0 * P3x * P3y * dP0z_exp) - (6300.0 * P1x * P3y * dP2z_exp) + (1575.0 * P5x * P1y * dP0z_exp) + (3150.0 * P3x * P1y * dP2z_exp) + (1575.0 * P1x * P1y * dP4z_exp);

    double dterm_3_dx = (2520.0 * dP5x_exp * P1z) - (6300.0 * dP3x_exp * P3z) - (6300.0 * dP3x_exp * P2y * P1z) + (1575.0 * dP1x_exp * P5z) + (3150.0 * dP1x_exp * P2y * P3z) + (1575.0 * dP1x_exp * P4y * P1z);
    double dterm_3_dy = (2520.0 * P5x * dP0y_exp * P1z) - (6300.0 * P3x * dP0y_exp * P3z) - (6300.0 * P3x * dP2y_exp * P1z) + (1575.0 * P1x * dP0y_exp * P5z) + (3150.0 * P1x * dP2y_exp * P3z) + (1575.0 * P1x * dP4y_exp * P1z);
    double dterm_3_dz = (2520.0 * P5x * dP1z_exp) - (6300.0 * P3x * dP3z_exp) - (6300.0 * P3x * P2y * dP1z_exp) + (1575.0 * P1x * dP5z_exp) + (3150.0 * P1x * P2y * dP3z_exp) + (1575.0 * P1x * P4y * dP1z_exp);

    double dterm_4_dx = (2520.0 * dP1x_exp * P5z) - (6300.0 * dP3x_exp * P3z) - (6300.0 * dP1x_exp * P2y * P3z) + (1575.0 * dP5x_exp * P1z) + (3150.0 * dP3x_exp * P2y * P1z) + (1575.0 * dP1x_exp * P4y * P1z);
    double dterm_4_dy = (2520.0 * P1x * dP0y_exp * P5z) - (6300.0 * P3x * dP0y_exp * P3z) - (6300.0 * P1x * dP2y_exp * P3z) + (1575.0 * P5x * dP0y_exp * P1z) + (3150.0 * P3x * dP2y_exp * P1z) + (1575.0 * P1x * dP4y_exp * P1z);
    double dterm_4_dz = (2520.0 * P1x * dP5z_exp) - (6300.0 * P3x * dP3z_exp) - (6300.0 * P1x * P2y * dP3z_exp) + (1575.0 * P5x * dP1z_exp) + (3150.0 * P3x * P2y * dP1z_exp) + (1575.0 * P1x * P4y * dP1z_exp);

    double dterm_5_dx = (2520.0 * dP0x_exp * P5y * P1z) - (6300.0 * dP0x_exp * P3y * P3z) - (6300.0 * dP2x_exp * P3y * P1z) + (1575.0 * dP0x_exp * P1y * P5z) + (3150.0 * dP2x_exp * P1y * P3z) + (1575.0 * dP4x_exp * P1y * P1z);
    double dterm_5_dy = (2520.0 * dP5y_exp * P1z) - (6300.0 * dP3y_exp * P3z) - (6300.0 * P2x * dP3y_exp * P1z) + (1575.0 * dP1y_exp * P5z) + (3150.0 * P2x * dP1y_exp * P3z) + (1575.0 * P4x * dP1y_exp * P1z);
    double dterm_5_dz = (2520.0 * P5y * dP1z_exp) - (6300.0 * P3y * dP3z_exp) - (6300.0 * P2x * P3y * dP1z_exp) + (1575.0 * P1y * dP5z_exp) + (3150.0 * P2x * P1y * dP3z_exp) + (1575.0 * P4x * P1y * dP1z_exp);

    double dterm_6_dx = (2520.0 * dP0x_exp * P1y * P5z) - (6300.0 * dP0x_exp * P3y * P3z) - (6300.0 * dP2x_exp * P1y * P3z) + (1575.0 * dP0x_exp * P5y * P1z) + (3150.0 * dP2x_exp * P3y * P1z) + (1575.0 * dP4x_exp * P1y * P1z);
    double dterm_6_dy = (2520.0 * dP1y_exp * P5z) - (6300.0 * dP3y_exp * P3z) - (6300.0 * P2x * dP1y_exp * P3z) + (1575.0 * dP5y_exp * P1z) + (3150.0 * P2x * dP3y_exp * P1z) + (1575.0 * P4x * dP1y_exp * P1z);
    double dterm_6_dz = (2520.0 * P1y * dP5z_exp) - (6300.0 * P3y * dP3z_exp) - (6300.0 * P2x * P1y * dP3z_exp) + (1575.0 * P5y * dP1z_exp) + (3150.0 * P2x * P3y * dP1z_exp) + (1575.0 * P4x * P1y * dP1z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

    // dmiu4 dx/dy/dz
    deriv[9] = temp * dterm_4_dx;
    deriv[10] = temp * dterm_4_dy;
    deriv[11] = temp * dterm_4_dz;

    // dmiu5 dx/dy/dz
    deriv[12] = temp * dterm_5_dx;
    deriv[13] = temp * dterm_5_dy;
    deriv[14] = temp * dterm_5_dz;

    // dmiu6 dx/dy/dz
    deriv[15] = temp * dterm_6_dx;
    deriv[16] = temp * dterm_6_dy;
    deriv[17] = temp * dterm_6_dz;

}

void calc_solid_MCSH_6_3(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double P4x = P4(lambda, x0, gamma);
    double P4y = P4(lambda, y0, gamma);
    double P4z = P4(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double term_1 = (5220.0 * P4x * P2y) - (360.0 * P6x) + (180.0 * P4x * P2z) - (4545.0 * P2x * P4y) - (4050.0 * P2x * P2y * P2z) + (495.0 * P2x * P4z) + (270.0 * P6y) + (495.0 * P4y * P2z) + (180.0 * P2y * P4z) - (45.0 * P6z);
    double term_2 = (5220.0 * P2x * P4y) - (360.0 * P6y) + (180.0 * P4y * P2z) - (4545.0 * P4x * P2y) - (4050.0 * P2x * P2y * P2z) + (495.0 * P2y * P4z) + (270.0 * P6x) + (495.0 * P4x * P2z) + (180.0 * P2x * P4z) - (45.0 * P6z);
    double term_3 = (5220.0 * P4x * P2z) - (360.0 * P6x) + (180.0 * P4x * P2y) - (4545.0 * P2x * P4z) - (4050.0 * P2x * P2y * P2z) + (495.0 * P2x * P4y) + (270.0 * P6z) + (495.0 * P2y * P4z) + (180.0 * P4y * P2z) - (45.0 * P6y);
    double term_4 = (5220.0 * P2x * P4z) - (360.0 * P6z) + (180.0 * P2y * P4z) - (4545.0 * P4x * P2z) - (4050.0 * P2x * P2y * P2z) + (495.0 * P4y * P2z) + (270.0 * P6x) + (495.0 * P4x * P2y) + (180.0 * P2x * P4y) - (45.0 * P6y);
    double term_5 = (5220.0 * P4y * P2z) - (360.0 * P6y) + (180.0 * P2x * P4y) - (4545.0 * P2y * P4z) - (4050.0 * P2x * P2y * P2z) + (495.0 * P4x * P2y) + (270.0 * P6z) + (495.0 * P2x * P4z) + (180.0 * P4x * P2z) - (45.0 * P6x);
    double term_6 = (5220.0 * P2y * P4z) - (360.0 * P6z) + (180.0 * P2x * P4z) - (4545.0 * P4y * P2z) - (4050.0 * P2x * P2y * P2z) + (495.0 * P4x * P2z) + (270.0 * P6y) + (495.0 * P2x * P4y) + (180.0 * P4x * P2y) - (45.0 * P6x);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP4x_exp = dP4_exp(P4x, C2, lambda, x0, gamma);
    double dP4y_exp = dP4_exp(P4y, C2, lambda, y0, gamma);
    double dP4z_exp = dP4_exp(P4z, C2, lambda, z0, gamma);

    double dP6x_exp = dP6_exp(P6x, C2, lambda, x0, gamma);
    double dP6y_exp = dP6_exp(P6y, C2, lambda, y0, gamma);
    double dP6z_exp = dP6_exp(P6z, C2, lambda, z0, gamma);

    double dterm_1_dx = (5220.0 * dP4x_exp * P2y) - (360.0 * dP6x_exp) + (180.0 * dP4x_exp * P2z) - (4545.0 * dP2x_exp * P4y) - (4050.0 * dP2x_exp * P2y * P2z) + (495.0 * dP2x_exp * P4z) + (270.0 * dP0x_exp * P6y) + (495.0 * dP0x_exp * P4y * P2z) + (180.0 * dP0x_exp * P2y * P4z) - (45.0 * dP0x_exp * P6z);
    double dterm_1_dy = (5220.0 * P4x * dP2y_exp) - (360.0 * P6x * dP0y_exp) + (180.0 * P4x * dP0y_exp * P2z) - (4545.0 * P2x * dP4y_exp) - (4050.0 * P2x * dP2y_exp * P2z) + (495.0 * P2x * dP0y_exp * P4z) + (270.0 * dP6y_exp) + (495.0 * dP4y_exp * P2z) + (180.0 * dP2y_exp * P4z) - (45.0 * dP0y_exp * P6z);
    double dterm_1_dz = (5220.0 * P4x * P2y * dP0z_exp) - (360.0 * P6x * dP0z_exp) + (180.0 * P4x * dP2z_exp) - (4545.0 * P2x * P4y * dP0z_exp) - (4050.0 * P2x * P2y * dP2z_exp) + (495.0 * P2x * dP4z_exp) + (270.0 * P6y * dP0z_exp) + (495.0 * P4y * dP2z_exp) + (180.0 * P2y * dP4z_exp) - (45.0 * dP6z_exp);

    double dterm_2_dx = (5220.0 * dP2x_exp * P4y) - (360.0 * dP0x_exp * P6y) + (180.0 * dP0x_exp * P4y * P2z) - (4545.0 * dP4x_exp * P2y) - (4050.0 * dP2x_exp * P2y * P2z) + (495.0 * dP0x_exp * P2y * P4z) + (270.0 * dP6x_exp) + (495.0 * dP4x_exp * P2z) + (180.0 * dP2x_exp * P4z) - (45.0 * dP0x_exp * P6z);
    double dterm_2_dy = (5220.0 * P2x * dP4y_exp) - (360.0 * dP6y_exp) + (180.0 * dP4y_exp * P2z) - (4545.0 * P4x * dP2y_exp) - (4050.0 * P2x * dP2y_exp * P2z) + (495.0 * dP2y_exp * P4z) + (270.0 * P6x * dP0y_exp) + (495.0 * P4x * dP0y_exp * P2z) + (180.0 * P2x * dP0y_exp * P4z) - (45.0 * dP0y_exp * P6z);
    double dterm_2_dz = (5220.0 * P2x * P4y * dP0z_exp) - (360.0 * P6y * dP0z_exp) + (180.0 * P4y * dP2z_exp) - (4545.0 * P4x * P2y * dP0z_exp) - (4050.0 * P2x * P2y * dP2z_exp) + (495.0 * P2y * dP4z_exp) + (270.0 * P6x * dP0z_exp) + (495.0 * P4x * dP2z_exp) + (180.0 * P2x * dP4z_exp) - (45.0 * dP6z_exp);

    double dterm_3_dx = (5220.0 * dP4x_exp * P2z) - (360.0 * dP6x_exp) + (180.0 * dP4x_exp * P2y) - (4545.0 * dP2x_exp * P4z) - (4050.0 * dP2x_exp * P2y * P2z) + (495.0 * dP2x_exp * P4y) + (270.0 * dP0x_exp * P6z) + (495.0 * dP0x_exp * P2y * P4z) + (180.0 * dP0x_exp * P4y * P2z) - (45.0 * dP0x_exp * P6y);
    double dterm_3_dy = (5220.0 * P4x * dP0y_exp * P2z) - (360.0 * P6x * dP0y_exp) + (180.0 * P4x * dP2y_exp) - (4545.0 * P2x * dP0y_exp * P4z) - (4050.0 * P2x * dP2y_exp * P2z) + (495.0 * P2x * dP4y_exp) + (270.0 * dP0y_exp * P6z) + (495.0 * dP2y_exp * P4z) + (180.0 * dP4y_exp * P2z) - (45.0 * dP6y_exp);
    double dterm_3_dz = (5220.0 * P4x * dP2z_exp) - (360.0 * P6x * dP0z_exp) + (180.0 * P4x * P2y * dP0z_exp) - (4545.0 * P2x * dP4z_exp) - (4050.0 * P2x * P2y * dP2z_exp) + (495.0 * P2x * P4y * dP0z_exp) + (270.0 * dP6z_exp) + (495.0 * P2y * dP4z_exp) + (180.0 * P4y * dP2z_exp) - (45.0 * P6y * dP0z_exp);

    double dterm_4_dx = (5220.0 * dP2x_exp * P4z) - (360.0 * dP0x_exp * P6z) + (180.0 * dP0x_exp * P2y * P4z) - (4545.0 * dP4x_exp * P2z) - (4050.0 * dP2x_exp * P2y * P2z) + (495.0 * dP0x_exp * P4y * P2z) + (270.0 * dP6x_exp) + (495.0 * dP4x_exp * P2y) + (180.0 * dP2x_exp * P4y) - (45.0 * dP0x_exp * P6y);
    double dterm_4_dy = (5220.0 * P2x * dP0y_exp * P4z) - (360.0 * dP0y_exp * P6z) + (180.0 * dP2y_exp * P4z) - (4545.0 * P4x * dP0y_exp * P2z) - (4050.0 * P2x * dP2y_exp * P2z) + (495.0 * dP4y_exp * P2z) + (270.0 * P6x * dP0y_exp) + (495.0 * P4x * dP2y_exp) + (180.0 * P2x * dP4y_exp) - (45.0 * dP6y_exp);
    double dterm_4_dz = (5220.0 * P2x * dP4z_exp) - (360.0 * dP6z_exp) + (180.0 * P2y * dP4z_exp) - (4545.0 * P4x * dP2z_exp) - (4050.0 * P2x * P2y * dP2z_exp) + (495.0 * P4y * dP2z_exp) + (270.0 * P6x * dP0z_exp) + (495.0 * P4x * P2y * dP0z_exp) + (180.0 * P2x * P4y * dP0z_exp) - (45.0 * P6y * dP0z_exp);

    double dterm_5_dx = (5220.0 * dP0x_exp * P4y * P2z) - (360.0 * dP0x_exp * P6y) + (180.0 * dP2x_exp * P4y) - (4545.0 * dP0x_exp * P2y * P4z) - (4050.0 * dP2x_exp * P2y * P2z) + (495.0 * dP4x_exp * P2y) + (270.0 * dP0x_exp * P6z) + (495.0 * dP2x_exp * P4z) + (180.0 * dP4x_exp * P2z) - (45.0 * dP6x_exp);
    double dterm_5_dy = (5220.0 * dP4y_exp * P2z) - (360.0 * dP6y_exp) + (180.0 * P2x * dP4y_exp) - (4545.0 * dP2y_exp * P4z) - (4050.0 * P2x * dP2y_exp * P2z) + (495.0 * P4x * dP2y_exp) + (270.0 * dP0y_exp * P6z) + (495.0 * P2x * dP0y_exp * P4z) + (180.0 * P4x * dP0y_exp * P2z) - (45.0 * P6x * dP0y_exp);
    double dterm_5_dz = (5220.0 * P4y * dP2z_exp) - (360.0 * P6y * dP0z_exp) + (180.0 * P2x * P4y * dP0z_exp) - (4545.0 * P2y * dP4z_exp) - (4050.0 * P2x * P2y * dP2z_exp) + (495.0 * P4x * P2y * dP0z_exp) + (270.0 * dP6z_exp) + (495.0 * P2x * dP4z_exp) + (180.0 * P4x * dP2z_exp) - (45.0 * P6x * dP0z_exp);

    double dterm_6_dx = (5220.0 * dP0x_exp * P2y * P4z) - (360.0 * dP0x_exp * P6z) + (180.0 * dP2x_exp * P4z) - (4545.0 * dP0x_exp * P4y * P2z) - (4050.0 * dP2x_exp * P2y * P2z) + (495.0 * dP4x_exp * P2z) + (270.0 * dP0x_exp * P6y) + (495.0 * dP2x_exp * P4y) + (180.0 * dP4x_exp * P2y) - (45.0 * dP6x_exp);
    double dterm_6_dy = (5220.0 * dP2y_exp * P4z) - (360.0 * dP0y_exp * P6z) + (180.0 * P2x * dP0y_exp * P4z) - (4545.0 * dP4y_exp * P2z) - (4050.0 * P2x * dP2y_exp * P2z) + (495.0 * P4x * dP0y_exp * P2z) + (270.0 * dP6y_exp) + (495.0 * P2x * dP4y_exp) + (180.0 * P4x * dP2y_exp) - (45.0 * P6x * dP0y_exp);
    double dterm_6_dz = (5220.0 * P2y * dP4z_exp) - (360.0 * dP6z_exp) + (180.0 * P2x * dP4z_exp) - (4545.0 * P4y * dP2z_exp) - (4050.0 * P2x * P2y * dP2z_exp) + (495.0 * P4x * dP2z_exp) + (270.0 * P6y * dP0z_exp) + (495.0 * P2x * P4y * dP0z_exp) + (180.0 * P4x * P2y * dP0z_exp) - (45.0 * P6x * dP0z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

    // dmiu4 dx/dy/dz
    deriv[9] = temp * dterm_4_dx;
    deriv[10] = temp * dterm_4_dy;
    deriv[11] = temp * dterm_4_dz;

    // dmiu5 dx/dy/dz
    deriv[12] = temp * dterm_5_dx;
    deriv[13] = temp * dterm_5_dy;
    deriv[14] = temp * dterm_5_dz;

    // dmiu6 dx/dy/dz
    deriv[15] = temp * dterm_6_dx;
    deriv[16] = temp * dterm_6_dy;
    deriv[17] = temp * dterm_6_dz;

}

void calc_solid_MCSH_6_4(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double term_1 = (5040.0 * P4x * P1y * P1z) - (5040.0 * P2x * P3y * P1z) - (5040.0 * P2x * P1y * P3z) + (315.0 * P5y * P1z) + (630.0 * P3y * P3z) + (315.0 * P1y * P5z);
    double term_2 = (5040.0 * P1x * P4y * P1z) - (5040.0 * P3x * P2y * P1z) - (5040.0 * P1x * P2y * P3z) + (315.0 * P5x * P1z) + (630.0 * P3x * P3z) + (315.0 * P1x * P5z);
    double term_3 = (5040.0 * P1x * P1y * P4z) - (5040.0 * P3x * P1y * P2z) - (5040.0 * P1x * P3y * P2z) + (315.0 * P5x * P1y) + (630.0 * P3x * P3y) + (315.0 * P1x * P5y);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP3x_exp = dP3_exp(P3x, C2, lambda, x0, gamma);
    double dP3y_exp = dP3_exp(P3y, C2, lambda, y0, gamma);
    double dP3z_exp = dP3_exp(P3z, C2, lambda, z0, gamma);

    double dP4x_exp = dP4_exp(P4x, C2, lambda, x0, gamma);
    double dP4y_exp = dP4_exp(P4y, C2, lambda, y0, gamma);
    double dP4z_exp = dP4_exp(P4z, C2, lambda, z0, gamma);

    double dP5x_exp = dP5_exp(P5x, C2, lambda, x0, gamma);
    double dP5y_exp = dP5_exp(P5y, C2, lambda, y0, gamma);
    double dP5z_exp = dP5_exp(P5z, C2, lambda, z0, gamma);

    double dterm_1_dx = (5040.0 * dP4x_exp * P1y * P1z) - (5040.0 * dP2x_exp * P3y * P1z) - (5040.0 * dP2x_exp * P1y * P3z) + (315.0 * dP0x_exp * P5y * P1z) + (630.0 * dP0x_exp * P3y * P3z) + (315.0 * dP0x_exp * P1y * P5z);
    double dterm_1_dy = (5040.0 * P4x * dP1y_exp * P1z) - (5040.0 * P2x * dP3y_exp * P1z) - (5040.0 * P2x * dP1y_exp * P3z) + (315.0 * dP5y_exp * P1z) + (630.0 * dP3y_exp * P3z) + (315.0 * dP1y_exp * P5z);
    double dterm_1_dz = (5040.0 * P4x * P1y * dP1z_exp) - (5040.0 * P2x * P3y * dP1z_exp) - (5040.0 * P2x * P1y * dP3z_exp) + (315.0 * P5y * dP1z_exp) + (630.0 * P3y * dP3z_exp) + (315.0 * P1y * dP5z_exp);

    double dterm_2_dx = (5040.0 * dP1x_exp * P4y * P1z) - (5040.0 * dP3x_exp * P2y * P1z) - (5040.0 * dP1x_exp * P2y * P3z) + (315.0 * dP5x_exp * P1z) + (630.0 * dP3x_exp * P3z) + (315.0 * dP1x_exp * P5z);
    double dterm_2_dy = (5040.0 * P1x * dP4y_exp * P1z) - (5040.0 * P3x * dP2y_exp * P1z) - (5040.0 * P1x * dP2y_exp * P3z) + (315.0 * P5x * dP0y_exp * P1z) + (630.0 * P3x * dP0y_exp * P3z) + (315.0 * P1x * dP0y_exp * P5z);
    double dterm_2_dz = (5040.0 * P1x * P4y * dP1z_exp) - (5040.0 * P3x * P2y * dP1z_exp) - (5040.0 * P1x * P2y * dP3z_exp) + (315.0 * P5x * dP1z_exp) + (630.0 * P3x * dP3z_exp) + (315.0 * P1x * dP5z_exp);

    double dterm_3_dx = (5040.0 * dP1x_exp * P1y * P4z) - (5040.0 * dP3x_exp * P1y * P2z) - (5040.0 * dP1x_exp * P3y * P2z) + (315.0 * dP5x_exp * P1y) + (630.0 * dP3x_exp * P3y) + (315.0 * dP1x_exp * P5y);
    double dterm_3_dy = (5040.0 * P1x * dP1y_exp * P4z) - (5040.0 * P3x * dP1y_exp * P2z) - (5040.0 * P1x * dP3y_exp * P2z) + (315.0 * P5x * dP1y_exp) + (630.0 * P3x * dP3y_exp) + (315.0 * P1x * dP5y_exp);
    double dterm_3_dz = (5040.0 * P1x * P1y * dP4z_exp) - (5040.0 * P3x * P1y * dP2z_exp) - (5040.0 * P1x * P3y * dP2z_exp) + (315.0 * P5x * P1y * dP0z_exp) + (630.0 * P3x * P3y * dP0z_exp) + (315.0 * P1x * P5y * dP0z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

}

void calc_solid_MCSH_6_5(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double term_1 = (6615.0 * P3x * P3y) - (1890.0 * P5x * P1y) - (945.0 * P3x * P1y * P2z) - (1890.0 * P1x * P5y) - (945.0 * P1x * P3y * P2z) + (945.0 * P1x * P1y * P4z);
    double term_2 = (6615.0 * P3x * P3z) - (1890.0 * P5x * P1z) - (945.0 * P3x * P2y * P1z) - (1890.0 * P1x * P5z) - (945.0 * P1x * P2y * P3z) + (945.0 * P1x * P4y * P1z);
    double term_3 = (6615.0 * P3y * P3z) - (1890.0 * P5y * P1z) - (945.0 * P2x * P3y * P1z) - (1890.0 * P1y * P5z) - (945.0 * P2x * P1y * P3z) + (945.0 * P4x * P1y * P1z);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP3x_exp = dP3_exp(P3x, C2, lambda, x0, gamma);
    double dP3y_exp = dP3_exp(P3y, C2, lambda, y0, gamma);
    double dP3z_exp = dP3_exp(P3z, C2, lambda, z0, gamma);

    double dP4x_exp = dP4_exp(P4x, C2, lambda, x0, gamma);
    double dP4y_exp = dP4_exp(P4y, C2, lambda, y0, gamma);
    double dP4z_exp = dP4_exp(P4z, C2, lambda, z0, gamma);

    double dP5x_exp = dP5_exp(P5x, C2, lambda, x0, gamma);
    double dP5y_exp = dP5_exp(P5y, C2, lambda, y0, gamma);
    double dP5z_exp = dP5_exp(P5z, C2, lambda, z0, gamma);

    double dterm_1_dx = (6615.0 * dP3x_exp * P3y) - (1890.0 * dP5x_exp * P1y) - (945.0 * dP3x_exp * P1y * P2z) - (1890.0 * dP1x_exp * P5y) - (945.0 * dP1x_exp * P3y * P2z) + (945.0 * dP1x_exp * P1y * P4z);
    double dterm_1_dy = (6615.0 * P3x * dP3y_exp) - (1890.0 * P5x * dP1y_exp) - (945.0 * P3x * dP1y_exp * P2z) - (1890.0 * P1x * dP5y_exp) - (945.0 * P1x * dP3y_exp * P2z) + (945.0 * P1x * dP1y_exp * P4z);
    double dterm_1_dz = (6615.0 * P3x * P3y * dP0z_exp) - (1890.0 * P5x * P1y * dP0z_exp) - (945.0 * P3x * P1y * dP2z_exp) - (1890.0 * P1x * P5y * dP0z_exp) - (945.0 * P1x * P3y * dP2z_exp) + (945.0 * P1x * P1y * dP4z_exp);

    double dterm_2_dx = (6615.0 * dP3x_exp * P3z) - (1890.0 * dP5x_exp * P1z) - (945.0 * dP3x_exp * P2y * P1z) - (1890.0 * dP1x_exp * P5z) - (945.0 * dP1x_exp * P2y * P3z) + (945.0 * dP1x_exp * P4y * P1z);
    double dterm_2_dy = (6615.0 * P3x * dP0y_exp * P3z) - (1890.0 * P5x * dP0y_exp * P1z) - (945.0 * P3x * dP2y_exp * P1z) - (1890.0 * P1x * dP0y_exp * P5z) - (945.0 * P1x * dP2y_exp * P3z) + (945.0 * P1x * dP4y_exp * P1z);
    double dterm_2_dz = (6615.0 * P3x * dP3z_exp) - (1890.0 * P5x * dP1z_exp) - (945.0 * P3x * P2y * dP1z_exp) - (1890.0 * P1x * dP5z_exp) - (945.0 * P1x * P2y * dP3z_exp) + (945.0 * P1x * P4y * dP1z_exp);

    double dterm_3_dx = (6615.0 * dP0x_exp * P3y * P3z) - (1890.0 * dP0x_exp * P5y * P1z) - (945.0 * dP2x_exp * P3y * P1z) - (1890.0 * dP0x_exp * P1y * P5z) - (945.0 * dP2x_exp * P1y * P3z) + (945.0 * dP4x_exp * P1y * P1z);
    double dterm_3_dy = (6615.0 * dP3y_exp * P3z) - (1890.0 * dP5y_exp * P1z) - (945.0 * P2x * dP3y_exp * P1z) - (1890.0 * dP1y_exp * P5z) - (945.0 * P2x * dP1y_exp * P3z) + (945.0 * P4x * dP1y_exp * P1z);
    double dterm_3_dz = (6615.0 * P3y * dP3z_exp) - (1890.0 * P5y * dP1z_exp) - (945.0 * P2x * P3y * dP1z_exp) - (1890.0 * P1y * dP5z_exp) - (945.0 * P2x * P1y * dP3z_exp) + (945.0 * P4x * P1y * dP1z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

}

void calc_solid_MCSH_6_6(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double term_1 = (7245.0 * P3x * P2y * P1z) - (630.0 * P5x * P1z) - (315.0 * P3x * P3z) - (2520.0 * P1x * P4y * P1z) - (2205.0 * P1x * P2y * P3z) + (315.0 * P1x * P5z);
    double term_2 = (7245.0 * P2x * P3y * P1z) - (630.0 * P5y * P1z) - (315.0 * P3y * P3z) - (2520.0 * P4x * P1y * P1z) - (2205.0 * P2x * P1y * P3z) + (315.0 * P1y * P5z);
    double term_3 = (7245.0 * P3x * P1y * P2z) - (630.0 * P5x * P1y) - (315.0 * P3x * P3y) - (2520.0 * P1x * P1y * P4z) - (2205.0 * P1x * P3y * P2z) + (315.0 * P1x * P5y);
    double term_4 = (7245.0 * P2x * P1y * P3z) - (630.0 * P1y * P5z) - (315.0 * P3y * P3z) - (2520.0 * P4x * P1y * P1z) - (2205.0 * P2x * P3y * P1z) + (315.0 * P5y * P1z);
    double term_5 = (7245.0 * P1x * P3y * P2z) - (630.0 * P1x * P5y) - (315.0 * P3x * P3y) - (2520.0 * P1x * P1y * P4z) - (2205.0 * P3x * P1y * P2z) + (315.0 * P5x * P1y);
    double term_6 = (7245.0 * P1x * P2y * P3z) - (630.0 * P1x * P5z) - (315.0 * P3x * P3z) - (2520.0 * P1x * P4y * P1z) - (2205.0 * P3x * P2y * P1z) + (315.0 * P5x * P1z);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP3x_exp = dP3_exp(P3x, C2, lambda, x0, gamma);
    double dP3y_exp = dP3_exp(P3y, C2, lambda, y0, gamma);
    double dP3z_exp = dP3_exp(P3z, C2, lambda, z0, gamma);

    double dP4x_exp = dP4_exp(P4x, C2, lambda, x0, gamma);
    double dP4y_exp = dP4_exp(P4y, C2, lambda, y0, gamma);
    double dP4z_exp = dP4_exp(P4z, C2, lambda, z0, gamma);

    double dP5x_exp = dP5_exp(P5x, C2, lambda, x0, gamma);
    double dP5y_exp = dP5_exp(P5y, C2, lambda, y0, gamma);
    double dP5z_exp = dP5_exp(P5z, C2, lambda, z0, gamma);

    double dterm_1_dx = (7245.0 * dP3x_exp * P2y * P1z) - (630.0 * dP5x_exp * P1z) - (315.0 * dP3x_exp * P3z) - (2520.0 * dP1x_exp * P4y * P1z) - (2205.0 * dP1x_exp * P2y * P3z) + (315.0 * dP1x_exp * P5z);
    double dterm_1_dy = (7245.0 * P3x * dP2y_exp * P1z) - (630.0 * P5x * dP0y_exp * P1z) - (315.0 * P3x * dP0y_exp * P3z) - (2520.0 * P1x * dP4y_exp * P1z) - (2205.0 * P1x * dP2y_exp * P3z) + (315.0 * P1x * dP0y_exp * P5z);
    double dterm_1_dz = (7245.0 * P3x * P2y * dP1z_exp) - (630.0 * P5x * dP1z_exp) - (315.0 * P3x * dP3z_exp) - (2520.0 * P1x * P4y * dP1z_exp) - (2205.0 * P1x * P2y * dP3z_exp) + (315.0 * P1x * dP5z_exp);

    double dterm_2_dx = (7245.0 * dP2x_exp * P3y * P1z) - (630.0 * dP0x_exp * P5y * P1z) - (315.0 * dP0x_exp * P3y * P3z) - (2520.0 * dP4x_exp * P1y * P1z) - (2205.0 * dP2x_exp * P1y * P3z) + (315.0 * dP0x_exp * P1y * P5z);
    double dterm_2_dy = (7245.0 * P2x * dP3y_exp * P1z) - (630.0 * dP5y_exp * P1z) - (315.0 * dP3y_exp * P3z) - (2520.0 * P4x * dP1y_exp * P1z) - (2205.0 * P2x * dP1y_exp * P3z) + (315.0 * dP1y_exp * P5z);
    double dterm_2_dz = (7245.0 * P2x * P3y * dP1z_exp) - (630.0 * P5y * dP1z_exp) - (315.0 * P3y * dP3z_exp) - (2520.0 * P4x * P1y * dP1z_exp) - (2205.0 * P2x * P1y * dP3z_exp) + (315.0 * P1y * dP5z_exp);

    double dterm_3_dx = (7245.0 * dP3x_exp * P1y * P2z) - (630.0 * dP5x_exp * P1y) - (315.0 * dP3x_exp * P3y) - (2520.0 * dP1x_exp * P1y * P4z) - (2205.0 * dP1x_exp * P3y * P2z) + (315.0 * dP1x_exp * P5y);
    double dterm_3_dy = (7245.0 * P3x * dP1y_exp * P2z) - (630.0 * P5x * dP1y_exp) - (315.0 * P3x * dP3y_exp) - (2520.0 * P1x * dP1y_exp * P4z) - (2205.0 * P1x * dP3y_exp * P2z) + (315.0 * P1x * dP5y_exp);
    double dterm_3_dz = (7245.0 * P3x * P1y * dP2z_exp) - (630.0 * P5x * P1y * dP0z_exp) - (315.0 * P3x * P3y * dP0z_exp) - (2520.0 * P1x * P1y * dP4z_exp) - (2205.0 * P1x * P3y * dP2z_exp) + (315.0 * P1x * P5y * dP0z_exp);

    double dterm_4_dx = (7245.0 * dP2x_exp * P1y * P3z) - (630.0 * dP0x_exp * P1y * P5z) - (315.0 * dP0x_exp * P3y * P3z) - (2520.0 * dP4x_exp * P1y * P1z) - (2205.0 * dP2x_exp * P3y * P1z) + (315.0 * dP0x_exp * P5y * P1z);
    double dterm_4_dy = (7245.0 * P2x * dP1y_exp * P3z) - (630.0 * dP1y_exp * P5z) - (315.0 * dP3y_exp * P3z) - (2520.0 * P4x * dP1y_exp * P1z) - (2205.0 * P2x * dP3y_exp * P1z) + (315.0 * dP5y_exp * P1z);
    double dterm_4_dz = (7245.0 * P2x * P1y * dP3z_exp) - (630.0 * P1y * dP5z_exp) - (315.0 * P3y * dP3z_exp) - (2520.0 * P4x * P1y * dP1z_exp) - (2205.0 * P2x * P3y * dP1z_exp) + (315.0 * P5y * dP1z_exp);

    double dterm_5_dx = (7245.0 * dP1x_exp * P3y * P2z) - (630.0 * dP1x_exp * P5y) - (315.0 * dP3x_exp * P3y) - (2520.0 * dP1x_exp * P1y * P4z) - (2205.0 * dP3x_exp * P1y * P2z) + (315.0 * dP5x_exp * P1y);
    double dterm_5_dy = (7245.0 * P1x * dP3y_exp * P2z) - (630.0 * P1x * dP5y_exp) - (315.0 * P3x * dP3y_exp) - (2520.0 * P1x * dP1y_exp * P4z) - (2205.0 * P3x * dP1y_exp * P2z) + (315.0 * P5x * dP1y_exp);
    double dterm_5_dz = (7245.0 * P1x * P3y * dP2z_exp) - (630.0 * P1x * P5y * dP0z_exp) - (315.0 * P3x * P3y * dP0z_exp) - (2520.0 * P1x * P1y * dP4z_exp) - (2205.0 * P3x * P1y * dP2z_exp) + (315.0 * P5x * P1y * dP0z_exp);

    double dterm_6_dx = (7245.0 * dP1x_exp * P2y * P3z) - (630.0 * dP1x_exp * P5z) - (315.0 * dP3x_exp * P3z) - (2520.0 * dP1x_exp * P4y * P1z) - (2205.0 * dP3x_exp * P2y * P1z) + (315.0 * dP5x_exp * P1z);
    double dterm_6_dy = (7245.0 * P1x * dP2y_exp * P3z) - (630.0 * P1x * dP0y_exp * P5z) - (315.0 * P3x * dP0y_exp * P3z) - (2520.0 * P1x * dP4y_exp * P1z) - (2205.0 * P3x * dP2y_exp * P1z) + (315.0 * P5x * dP0y_exp * P1z);
    double dterm_6_dz = (7245.0 * P1x * P2y * dP3z_exp) - (630.0 * P1x * dP5z_exp) - (315.0 * P3x * dP3z_exp) - (2520.0 * P1x * P4y * dP1z_exp) - (2205.0 * P3x * P2y * dP1z_exp) + (315.0 * P5x * dP1z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

    // dmiu4 dx/dy/dz
    deriv[9] = temp * dterm_4_dx;
    deriv[10] = temp * dterm_4_dy;
    deriv[11] = temp * dterm_4_dz;

    // dmiu5 dx/dy/dz
    deriv[12] = temp * dterm_5_dx;
    deriv[13] = temp * dterm_5_dy;
    deriv[14] = temp * dterm_5_dz;

    // dmiu6 dx/dy/dz
    deriv[15] = temp * dterm_6_dx;
    deriv[16] = temp * dterm_6_dy;
    deriv[17] = temp * dterm_6_dz;

}

void calc_solid_MCSH_6_7(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

    double P2x = P2(lambda, x0, gamma);
    double P2y = P2(lambda, y0, gamma);
    double P2z = P2(lambda, z0, gamma);

    double P4x = P4(lambda, x0, gamma);
    double P4y = P4(lambda, y0, gamma);
    double P4z = P4(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double term_1 = (8100.0 * P2x * P2y * P2z) - (675.0 * P4x * P2y) - (675.0 * P2x * P4y) - (675.0 * P4x * P2z) - (675.0 * P2x * P4z) - (675.0 * P4y * P2z) - (675.0 * P2y * P4z) + (90.0 * P6x) + (90.0 * P6y) + (90.0 * P6z);

    double m = temp * term_1;

    value[0] = m;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP4x_exp = dP4_exp(P4x, C2, lambda, x0, gamma);
    double dP4y_exp = dP4_exp(P4y, C2, lambda, y0, gamma);
    double dP4z_exp = dP4_exp(P4z, C2, lambda, z0, gamma);

    double dP6x_exp = dP6_exp(P6x, C2, lambda, x0, gamma);
    double dP6y_exp = dP6_exp(P6y, C2, lambda, y0, gamma);
    double dP6z_exp = dP6_exp(P6z, C2, lambda, z0, gamma);

    double dterm_1_dx = (8100.0 * dP2x_exp * P2y * P2z) - (675.0 * dP4x_exp * P2y) - (675.0 * dP2x_exp * P4y) - (675.0 * dP4x_exp * P2z) - (675.0 * dP2x_exp * P4z) - (675.0 * dP0x_exp * P4y * P2z) - (675.0 * dP0x_exp * P2y * P4z) + (90.0 * dP6x_exp) + (90.0 * dP0x_exp * P6y) + (90.0 * dP0x_exp * P6z);
    double dterm_1_dy = (8100.0 * P2x * dP2y_exp * P2z) - (675.0 * P4x * dP2y_exp) - (675.0 * P2x * dP4y_exp) - (675.0 * P4x * dP0y_exp * P2z) - (675.0 * P2x * dP0y_exp * P4z) - (675.0 * dP4y_exp * P2z) - (675.0 * dP2y_exp * P4z) + (90.0 * P6x * dP0y_exp) + (90.0 * dP6y_exp) + (90.0 * dP0y_exp * P6z);
    double dterm_1_dz = (8100.0 * P2x * P2y * dP2z_exp) - (675.0 * P4x * P2y * dP0z_exp) - (675.0 * P2x * P4y * dP0z_exp) - (675.0 * P4x * dP2z_exp) - (675.0 * P2x * dP4z_exp) - (675.0 * P4y * dP2z_exp) - (675.0 * P2y * dP4z_exp) + (90.0 * P6x * dP0z_exp) + (90.0 * P6y * dP0z_exp) + (90.0 * dP6z_exp);



    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

}


void calc_solid_MCSH_7_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double term_1 = (5040.0 * P7x) - (52920.0 * P5x * P2y) - (52920.0 * P5x * P2z) + (66150.0 * P3x * P4y) + (132300.0 * P3x * P2y * P2z) + (66150.0 * P3x * P4z) - (11025.0 * P1x * P6y) - (33075.0 * P1x * P4y * P2z) - (33075.0 * P1x * P2y * P4z) - (11025.0 * P1x * P6z);
    double term_2 = (5040.0 * P7y) - (52920.0 * P2x * P5y) - (52920.0 * P5y * P2z) + (66150.0 * P4x * P3y) + (132300.0 * P2x * P3y * P2z) + (66150.0 * P3y * P4z) - (11025.0 * P6x * P1y) - (33075.0 * P4x * P1y * P2z) - (33075.0 * P2x * P1y * P4z) - (11025.0 * P1y * P6z);
    double term_3 = (5040.0 * P7z) - (52920.0 * P2x * P5z) - (52920.0 * P2y * P5z) + (66150.0 * P4x * P3z) + (132300.0 * P2x * P2y * P3z) + (66150.0 * P4y * P3z) - (11025.0 * P6x * P1z) - (33075.0 * P4x * P2y * P1z) - (33075.0 * P2x * P4y * P1z) - (11025.0 * P6y * P1z);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP3x_exp = dP3_exp(P3x, C2, lambda, x0, gamma);
    double dP3y_exp = dP3_exp(P3y, C2, lambda, y0, gamma);
    double dP3z_exp = dP3_exp(P3z, C2, lambda, z0, gamma);

    double dP4x_exp = dP4_exp(P4x, C2, lambda, x0, gamma);
    double dP4y_exp = dP4_exp(P4y, C2, lambda, y0, gamma);
    double dP4z_exp = dP4_exp(P4z, C2, lambda, z0, gamma);

    double dP5x_exp = dP5_exp(P5x, C2, lambda, x0, gamma);
    double dP5y_exp = dP5_exp(P5y, C2, lambda, y0, gamma);
    double dP5z_exp = dP5_exp(P5z, C2, lambda, z0, gamma);

    double dP6x_exp = dP6_exp(P6x, C2, lambda, x0, gamma);
    double dP6y_exp = dP6_exp(P6y, C2, lambda, y0, gamma);
    double dP6z_exp = dP6_exp(P6z, C2, lambda, z0, gamma);

    double dP7x_exp = dP7_exp(P7x, C2, lambda, x0, gamma);
    double dP7y_exp = dP7_exp(P7y, C2, lambda, y0, gamma);
    double dP7z_exp = dP7_exp(P7z, C2, lambda, z0, gamma);

    double dterm_1_dx = (5040.0 * dP7x_exp) - (52920.0 * dP5x_exp * P2y) - (52920.0 * dP5x_exp * P2z) + (66150.0 * dP3x_exp * P4y) + (132300.0 * dP3x_exp * P2y * P2z) + (66150.0 * dP3x_exp * P4z) - (11025.0 * dP1x_exp * P6y) - (33075.0 * dP1x_exp * P4y * P2z) - (33075.0 * dP1x_exp * P2y * P4z) - (11025.0 * dP1x_exp * P6z);
    double dterm_1_dy = (5040.0 * P7x * dP0y_exp) - (52920.0 * P5x * dP2y_exp) - (52920.0 * P5x * dP0y_exp * P2z) + (66150.0 * P3x * dP4y_exp) + (132300.0 * P3x * dP2y_exp * P2z) + (66150.0 * P3x * dP0y_exp * P4z) - (11025.0 * P1x * dP6y_exp) - (33075.0 * P1x * dP4y_exp * P2z) - (33075.0 * P1x * dP2y_exp * P4z) - (11025.0 * P1x * dP0y_exp * P6z);
    double dterm_1_dz = (5040.0 * P7x * dP0z_exp) - (52920.0 * P5x * P2y * dP0z_exp) - (52920.0 * P5x * dP2z_exp) + (66150.0 * P3x * P4y * dP0z_exp) + (132300.0 * P3x * P2y * dP2z_exp) + (66150.0 * P3x * dP4z_exp) - (11025.0 * P1x * P6y * dP0z_exp) - (33075.0 * P1x * P4y * dP2z_exp) - (33075.0 * P1x * P2y * dP4z_exp) - (11025.0 * P1x * dP6z_exp);

    double dterm_2_dx = (5040.0 * dP0x_exp * P7y) - (52920.0 * dP2x_exp * P5y) - (52920.0 * dP0x_exp * P5y * P2z) + (66150.0 * dP4x_exp * P3y) + (132300.0 * dP2x_exp * P3y * P2z) + (66150.0 * dP0x_exp * P3y * P4z) - (11025.0 * dP6x_exp * P1y) - (33075.0 * dP4x_exp * P1y * P2z) - (33075.0 * dP2x_exp * P1y * P4z) - (11025.0 * dP0x_exp * P1y * P6z);
    double dterm_2_dy = (5040.0 * dP7y_exp) - (52920.0 * P2x * dP5y_exp) - (52920.0 * dP5y_exp * P2z) + (66150.0 * P4x * dP3y_exp) + (132300.0 * P2x * dP3y_exp * P2z) + (66150.0 * dP3y_exp * P4z) - (11025.0 * P6x * dP1y_exp) - (33075.0 * P4x * dP1y_exp * P2z) - (33075.0 * P2x * dP1y_exp * P4z) - (11025.0 * dP1y_exp * P6z);
    double dterm_2_dz = (5040.0 * P7y * dP0z_exp) - (52920.0 * P2x * P5y * dP0z_exp) - (52920.0 * P5y * dP2z_exp) + (66150.0 * P4x * P3y * dP0z_exp) + (132300.0 * P2x * P3y * dP2z_exp) + (66150.0 * P3y * dP4z_exp) - (11025.0 * P6x * P1y * dP0z_exp) - (33075.0 * P4x * P1y * dP2z_exp) - (33075.0 * P2x * P1y * dP4z_exp) - (11025.0 * P1y * dP6z_exp);

    double dterm_3_dx = (5040.0 * dP0x_exp * P7z) - (52920.0 * dP2x_exp * P5z) - (52920.0 * dP0x_exp * P2y * P5z) + (66150.0 * dP4x_exp * P3z) + (132300.0 * dP2x_exp * P2y * P3z) + (66150.0 * dP0x_exp * P4y * P3z) - (11025.0 * dP6x_exp * P1z) - (33075.0 * dP4x_exp * P2y * P1z) - (33075.0 * dP2x_exp * P4y * P1z) - (11025.0 * dP0x_exp * P6y * P1z);
    double dterm_3_dy = (5040.0 * dP0y_exp * P7z) - (52920.0 * P2x * dP0y_exp * P5z) - (52920.0 * dP2y_exp * P5z) + (66150.0 * P4x * dP0y_exp * P3z) + (132300.0 * P2x * dP2y_exp * P3z) + (66150.0 * dP4y_exp * P3z) - (11025.0 * P6x * dP0y_exp * P1z) - (33075.0 * P4x * dP2y_exp * P1z) - (33075.0 * P2x * dP4y_exp * P1z) - (11025.0 * dP6y_exp * P1z);
    double dterm_3_dz = (5040.0 * dP7z_exp) - (52920.0 * P2x * dP5z_exp) - (52920.0 * P2y * dP5z_exp) + (66150.0 * P4x * dP3z_exp) + (132300.0 * P2x * P2y * dP3z_exp) + (66150.0 * P4y * dP3z_exp) - (11025.0 * P6x * dP1z_exp) - (33075.0 * P4x * P2y * dP1z_exp) - (33075.0 * P2x * P4y * dP1z_exp) - (11025.0 * P6y * dP1z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

}

void calc_solid_MCSH_7_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double term_1 = (20160.0 * P6x * P1y) - (75600.0 * P4x * P3y) - (75600.0 * P4x * P1y * P2z) + (37800.0 * P2x * P5y) + (75600.0 * P2x * P3y * P2z) + (37800.0 * P2x * P1y * P4z) - (1575.0 * P7y) - (4725.0 * P5y * P2z) - (4725.0 * P3y * P4z) - (1575.0 * P1y * P6z);
    double term_2 = (20160.0 * P1x * P6y) - (75600.0 * P3x * P4y) - (75600.0 * P1x * P4y * P2z) + (37800.0 * P5x * P2y) + (75600.0 * P3x * P2y * P2z) + (37800.0 * P1x * P2y * P4z) - (1575.0 * P7x) - (4725.0 * P5x * P2z) - (4725.0 * P3x * P4z) - (1575.0 * P1x * P6z);
    double term_3 = (20160.0 * P6x * P1z) - (75600.0 * P4x * P3z) - (75600.0 * P4x * P2y * P1z) + (37800.0 * P2x * P5z) + (75600.0 * P2x * P2y * P3z) + (37800.0 * P2x * P4y * P1z) - (1575.0 * P7z) - (4725.0 * P2y * P5z) - (4725.0 * P4y * P3z) - (1575.0 * P6y * P1z);
    double term_4 = (20160.0 * P1x * P6z) - (75600.0 * P3x * P4z) - (75600.0 * P1x * P2y * P4z) + (37800.0 * P5x * P2z) + (75600.0 * P3x * P2y * P2z) + (37800.0 * P1x * P4y * P2z) - (1575.0 * P7x) - (4725.0 * P5x * P2y) - (4725.0 * P3x * P4y) - (1575.0 * P1x * P6y);
    double term_5 = (20160.0 * P6y * P1z) - (75600.0 * P4y * P3z) - (75600.0 * P2x * P4y * P1z) + (37800.0 * P2y * P5z) + (75600.0 * P2x * P2y * P3z) + (37800.0 * P4x * P2y * P1z) - (1575.0 * P7z) - (4725.0 * P2x * P5z) - (4725.0 * P4x * P3z) - (1575.0 * P6x * P1z);
    double term_6 = (20160.0 * P1y * P6z) - (75600.0 * P3y * P4z) - (75600.0 * P2x * P1y * P4z) + (37800.0 * P5y * P2z) + (75600.0 * P2x * P3y * P2z) + (37800.0 * P4x * P1y * P2z) - (1575.0 * P7y) - (4725.0 * P2x * P5y) - (4725.0 * P4x * P3y) - (1575.0 * P6x * P1y);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP3x_exp = dP3_exp(P3x, C2, lambda, x0, gamma);
    double dP3y_exp = dP3_exp(P3y, C2, lambda, y0, gamma);
    double dP3z_exp = dP3_exp(P3z, C2, lambda, z0, gamma);

    double dP4x_exp = dP4_exp(P4x, C2, lambda, x0, gamma);
    double dP4y_exp = dP4_exp(P4y, C2, lambda, y0, gamma);
    double dP4z_exp = dP4_exp(P4z, C2, lambda, z0, gamma);

    double dP5x_exp = dP5_exp(P5x, C2, lambda, x0, gamma);
    double dP5y_exp = dP5_exp(P5y, C2, lambda, y0, gamma);
    double dP5z_exp = dP5_exp(P5z, C2, lambda, z0, gamma);

    double dP6x_exp = dP6_exp(P6x, C2, lambda, x0, gamma);
    double dP6y_exp = dP6_exp(P6y, C2, lambda, y0, gamma);
    double dP6z_exp = dP6_exp(P6z, C2, lambda, z0, gamma);

    double dP7x_exp = dP7_exp(P7x, C2, lambda, x0, gamma);
    double dP7y_exp = dP7_exp(P7y, C2, lambda, y0, gamma);
    double dP7z_exp = dP7_exp(P7z, C2, lambda, z0, gamma);

    double dterm_1_dx = (20160.0 * dP6x_exp * P1y) - (75600.0 * dP4x_exp * P3y) - (75600.0 * dP4x_exp * P1y * P2z) + (37800.0 * dP2x_exp * P5y) + (75600.0 * dP2x_exp * P3y * P2z) + (37800.0 * dP2x_exp * P1y * P4z) - (1575.0 * dP0x_exp * P7y) - (4725.0 * dP0x_exp * P5y * P2z) - (4725.0 * dP0x_exp * P3y * P4z) - (1575.0 * dP0x_exp * P1y * P6z);
    double dterm_1_dy = (20160.0 * P6x * dP1y_exp) - (75600.0 * P4x * dP3y_exp) - (75600.0 * P4x * dP1y_exp * P2z) + (37800.0 * P2x * dP5y_exp) + (75600.0 * P2x * dP3y_exp * P2z) + (37800.0 * P2x * dP1y_exp * P4z) - (1575.0 * dP7y_exp) - (4725.0 * dP5y_exp * P2z) - (4725.0 * dP3y_exp * P4z) - (1575.0 * dP1y_exp * P6z);
    double dterm_1_dz = (20160.0 * P6x * P1y * dP0z_exp) - (75600.0 * P4x * P3y * dP0z_exp) - (75600.0 * P4x * P1y * dP2z_exp) + (37800.0 * P2x * P5y * dP0z_exp) + (75600.0 * P2x * P3y * dP2z_exp) + (37800.0 * P2x * P1y * dP4z_exp) - (1575.0 * P7y * dP0z_exp) - (4725.0 * P5y * dP2z_exp) - (4725.0 * P3y * dP4z_exp) - (1575.0 * P1y * dP6z_exp);

    double dterm_2_dx = (20160.0 * dP1x_exp * P6y) - (75600.0 * dP3x_exp * P4y) - (75600.0 * dP1x_exp * P4y * P2z) + (37800.0 * dP5x_exp * P2y) + (75600.0 * dP3x_exp * P2y * P2z) + (37800.0 * dP1x_exp * P2y * P4z) - (1575.0 * dP7x_exp) - (4725.0 * dP5x_exp * P2z) - (4725.0 * dP3x_exp * P4z) - (1575.0 * dP1x_exp * P6z);
    double dterm_2_dy = (20160.0 * P1x * dP6y_exp) - (75600.0 * P3x * dP4y_exp) - (75600.0 * P1x * dP4y_exp * P2z) + (37800.0 * P5x * dP2y_exp) + (75600.0 * P3x * dP2y_exp * P2z) + (37800.0 * P1x * dP2y_exp * P4z) - (1575.0 * P7x * dP0y_exp) - (4725.0 * P5x * dP0y_exp * P2z) - (4725.0 * P3x * dP0y_exp * P4z) - (1575.0 * P1x * dP0y_exp * P6z);
    double dterm_2_dz = (20160.0 * P1x * P6y * dP0z_exp) - (75600.0 * P3x * P4y * dP0z_exp) - (75600.0 * P1x * P4y * dP2z_exp) + (37800.0 * P5x * P2y * dP0z_exp) + (75600.0 * P3x * P2y * dP2z_exp) + (37800.0 * P1x * P2y * dP4z_exp) - (1575.0 * P7x * dP0z_exp) - (4725.0 * P5x * dP2z_exp) - (4725.0 * P3x * dP4z_exp) - (1575.0 * P1x * dP6z_exp);

    double dterm_3_dx = (20160.0 * dP6x_exp * P1z) - (75600.0 * dP4x_exp * P3z) - (75600.0 * dP4x_exp * P2y * P1z) + (37800.0 * dP2x_exp * P5z) + (75600.0 * dP2x_exp * P2y * P3z) + (37800.0 * dP2x_exp * P4y * P1z) - (1575.0 * dP0x_exp * P7z) - (4725.0 * dP0x_exp * P2y * P5z) - (4725.0 * dP0x_exp * P4y * P3z) - (1575.0 * dP0x_exp * P6y * P1z);
    double dterm_3_dy = (20160.0 * P6x * dP0y_exp * P1z) - (75600.0 * P4x * dP0y_exp * P3z) - (75600.0 * P4x * dP2y_exp * P1z) + (37800.0 * P2x * dP0y_exp * P5z) + (75600.0 * P2x * dP2y_exp * P3z) + (37800.0 * P2x * dP4y_exp * P1z) - (1575.0 * dP0y_exp * P7z) - (4725.0 * dP2y_exp * P5z) - (4725.0 * dP4y_exp * P3z) - (1575.0 * dP6y_exp * P1z);
    double dterm_3_dz = (20160.0 * P6x * dP1z_exp) - (75600.0 * P4x * dP3z_exp) - (75600.0 * P4x * P2y * dP1z_exp) + (37800.0 * P2x * dP5z_exp) + (75600.0 * P2x * P2y * dP3z_exp) + (37800.0 * P2x * P4y * dP1z_exp) - (1575.0 * dP7z_exp) - (4725.0 * P2y * dP5z_exp) - (4725.0 * P4y * dP3z_exp) - (1575.0 * P6y * dP1z_exp);

    double dterm_4_dx = (20160.0 * dP1x_exp * P6z) - (75600.0 * dP3x_exp * P4z) - (75600.0 * dP1x_exp * P2y * P4z) + (37800.0 * dP5x_exp * P2z) + (75600.0 * dP3x_exp * P2y * P2z) + (37800.0 * dP1x_exp * P4y * P2z) - (1575.0 * dP7x_exp) - (4725.0 * dP5x_exp * P2y) - (4725.0 * dP3x_exp * P4y) - (1575.0 * dP1x_exp * P6y);
    double dterm_4_dy = (20160.0 * P1x * dP0y_exp * P6z) - (75600.0 * P3x * dP0y_exp * P4z) - (75600.0 * P1x * dP2y_exp * P4z) + (37800.0 * P5x * dP0y_exp * P2z) + (75600.0 * P3x * dP2y_exp * P2z) + (37800.0 * P1x * dP4y_exp * P2z) - (1575.0 * P7x * dP0y_exp) - (4725.0 * P5x * dP2y_exp) - (4725.0 * P3x * dP4y_exp) - (1575.0 * P1x * dP6y_exp);
    double dterm_4_dz = (20160.0 * P1x * dP6z_exp) - (75600.0 * P3x * dP4z_exp) - (75600.0 * P1x * P2y * dP4z_exp) + (37800.0 * P5x * dP2z_exp) + (75600.0 * P3x * P2y * dP2z_exp) + (37800.0 * P1x * P4y * dP2z_exp) - (1575.0 * P7x * dP0z_exp) - (4725.0 * P5x * P2y * dP0z_exp) - (4725.0 * P3x * P4y * dP0z_exp) - (1575.0 * P1x * P6y * dP0z_exp);

    double dterm_5_dx = (20160.0 * dP0x_exp * P6y * P1z) - (75600.0 * dP0x_exp * P4y * P3z) - (75600.0 * dP2x_exp * P4y * P1z) + (37800.0 * dP0x_exp * P2y * P5z) + (75600.0 * dP2x_exp * P2y * P3z) + (37800.0 * dP4x_exp * P2y * P1z) - (1575.0 * dP0x_exp * P7z) - (4725.0 * dP2x_exp * P5z) - (4725.0 * dP4x_exp * P3z) - (1575.0 * dP6x_exp * P1z);
    double dterm_5_dy = (20160.0 * dP6y_exp * P1z) - (75600.0 * dP4y_exp * P3z) - (75600.0 * P2x * dP4y_exp * P1z) + (37800.0 * dP2y_exp * P5z) + (75600.0 * P2x * dP2y_exp * P3z) + (37800.0 * P4x * dP2y_exp * P1z) - (1575.0 * dP0y_exp * P7z) - (4725.0 * P2x * dP0y_exp * P5z) - (4725.0 * P4x * dP0y_exp * P3z) - (1575.0 * P6x * dP0y_exp * P1z);
    double dterm_5_dz = (20160.0 * P6y * dP1z_exp) - (75600.0 * P4y * dP3z_exp) - (75600.0 * P2x * P4y * dP1z_exp) + (37800.0 * P2y * dP5z_exp) + (75600.0 * P2x * P2y * dP3z_exp) + (37800.0 * P4x * P2y * dP1z_exp) - (1575.0 * dP7z_exp) - (4725.0 * P2x * dP5z_exp) - (4725.0 * P4x * dP3z_exp) - (1575.0 * P6x * dP1z_exp);

    double dterm_6_dx = (20160.0 * dP0x_exp * P1y * P6z) - (75600.0 * dP0x_exp * P3y * P4z) - (75600.0 * dP2x_exp * P1y * P4z) + (37800.0 * dP0x_exp * P5y * P2z) + (75600.0 * dP2x_exp * P3y * P2z) + (37800.0 * dP4x_exp * P1y * P2z) - (1575.0 * dP0x_exp * P7y) - (4725.0 * dP2x_exp * P5y) - (4725.0 * dP4x_exp * P3y) - (1575.0 * dP6x_exp * P1y);
    double dterm_6_dy = (20160.0 * dP1y_exp * P6z) - (75600.0 * dP3y_exp * P4z) - (75600.0 * P2x * dP1y_exp * P4z) + (37800.0 * dP5y_exp * P2z) + (75600.0 * P2x * dP3y_exp * P2z) + (37800.0 * P4x * dP1y_exp * P2z) - (1575.0 * dP7y_exp) - (4725.0 * P2x * dP5y_exp) - (4725.0 * P4x * dP3y_exp) - (1575.0 * P6x * dP1y_exp);
    double dterm_6_dz = (20160.0 * P1y * dP6z_exp) - (75600.0 * P3y * dP4z_exp) - (75600.0 * P2x * P1y * dP4z_exp) + (37800.0 * P5y * dP2z_exp) + (75600.0 * P2x * P3y * dP2z_exp) + (37800.0 * P4x * P1y * dP2z_exp) - (1575.0 * P7y * dP0z_exp) - (4725.0 * P2x * P5y * dP0z_exp) - (4725.0 * P4x * P3y * dP0z_exp) - (1575.0 * P6x * P1y * dP0z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

    // dmiu4 dx/dy/dz
    deriv[9] = temp * dterm_4_dx;
    deriv[10] = temp * dterm_4_dy;
    deriv[11] = temp * dterm_4_dz;

    // dmiu5 dx/dy/dz
    deriv[12] = temp * dterm_5_dx;
    deriv[13] = temp * dterm_5_dy;
    deriv[14] = temp * dterm_5_dz;

    // dmiu6 dx/dy/dz
    deriv[15] = temp * dterm_6_dx;
    deriv[16] = temp * dterm_6_dy;
    deriv[17] = temp * dterm_6_dz;

}

void calc_solid_MCSH_7_3(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double term_1 = (49140.0 * P5x * P2y) - (2520.0 * P7x) + (3780.0 * P5x * P2z) - (70875.0 * P3x * P4y) - (66150.0 * P3x * P2y * P2z) + (4725.0 * P3x * P4z) + (12600.0 * P1x * P6y) + (23625.0 * P1x * P4y * P2z) + (9450.0 * P1x * P2y * P4z) - (1575.0 * P1x * P6z);
    double term_2 = (49140.0 * P2x * P5y) - (2520.0 * P7y) + (3780.0 * P5y * P2z) - (70875.0 * P4x * P3y) - (66150.0 * P2x * P3y * P2z) + (4725.0 * P3y * P4z) + (12600.0 * P6x * P1y) + (23625.0 * P4x * P1y * P2z) + (9450.0 * P2x * P1y * P4z) - (1575.0 * P1y * P6z);
    double term_3 = (49140.0 * P5x * P2z) - (2520.0 * P7x) + (3780.0 * P5x * P2y) - (70875.0 * P3x * P4z) - (66150.0 * P3x * P2y * P2z) + (4725.0 * P3x * P4y) + (12600.0 * P1x * P6z) + (23625.0 * P1x * P2y * P4z) + (9450.0 * P1x * P4y * P2z) - (1575.0 * P1x * P6y);
    double term_4 = (49140.0 * P2x * P5z) - (2520.0 * P7z) + (3780.0 * P2y * P5z) - (70875.0 * P4x * P3z) - (66150.0 * P2x * P2y * P3z) + (4725.0 * P4y * P3z) + (12600.0 * P6x * P1z) + (23625.0 * P4x * P2y * P1z) + (9450.0 * P2x * P4y * P1z) - (1575.0 * P6y * P1z);
    double term_5 = (49140.0 * P5y * P2z) - (2520.0 * P7y) + (3780.0 * P2x * P5y) - (70875.0 * P3y * P4z) - (66150.0 * P2x * P3y * P2z) + (4725.0 * P4x * P3y) + (12600.0 * P1y * P6z) + (23625.0 * P2x * P1y * P4z) + (9450.0 * P4x * P1y * P2z) - (1575.0 * P6x * P1y);
    double term_6 = (49140.0 * P2y * P5z) - (2520.0 * P7z) + (3780.0 * P2x * P5z) - (70875.0 * P4y * P3z) - (66150.0 * P2x * P2y * P3z) + (4725.0 * P4x * P3z) + (12600.0 * P6y * P1z) + (23625.0 * P2x * P4y * P1z) + (9450.0 * P4x * P2y * P1z) - (1575.0 * P6x * P1z);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP3x_exp = dP3_exp(P3x, C2, lambda, x0, gamma);
    double dP3y_exp = dP3_exp(P3y, C2, lambda, y0, gamma);
    double dP3z_exp = dP3_exp(P3z, C2, lambda, z0, gamma);

    double dP4x_exp = dP4_exp(P4x, C2, lambda, x0, gamma);
    double dP4y_exp = dP4_exp(P4y, C2, lambda, y0, gamma);
    double dP4z_exp = dP4_exp(P4z, C2, lambda, z0, gamma);

    double dP5x_exp = dP5_exp(P5x, C2, lambda, x0, gamma);
    double dP5y_exp = dP5_exp(P5y, C2, lambda, y0, gamma);
    double dP5z_exp = dP5_exp(P5z, C2, lambda, z0, gamma);

    double dP6x_exp = dP6_exp(P6x, C2, lambda, x0, gamma);
    double dP6y_exp = dP6_exp(P6y, C2, lambda, y0, gamma);
    double dP6z_exp = dP6_exp(P6z, C2, lambda, z0, gamma);

    double dP7x_exp = dP7_exp(P7x, C2, lambda, x0, gamma);
    double dP7y_exp = dP7_exp(P7y, C2, lambda, y0, gamma);
    double dP7z_exp = dP7_exp(P7z, C2, lambda, z0, gamma);

    double dterm_1_dx = (49140.0 * dP5x_exp * P2y) - (2520.0 * dP7x_exp) + (3780.0 * dP5x_exp * P2z) - (70875.0 * dP3x_exp * P4y) - (66150.0 * dP3x_exp * P2y * P2z) + (4725.0 * dP3x_exp * P4z) + (12600.0 * dP1x_exp * P6y) + (23625.0 * dP1x_exp * P4y * P2z) + (9450.0 * dP1x_exp * P2y * P4z) - (1575.0 * dP1x_exp * P6z);
    double dterm_1_dy = (49140.0 * P5x * dP2y_exp) - (2520.0 * P7x * dP0y_exp) + (3780.0 * P5x * dP0y_exp * P2z) - (70875.0 * P3x * dP4y_exp) - (66150.0 * P3x * dP2y_exp * P2z) + (4725.0 * P3x * dP0y_exp * P4z) + (12600.0 * P1x * dP6y_exp) + (23625.0 * P1x * dP4y_exp * P2z) + (9450.0 * P1x * dP2y_exp * P4z) - (1575.0 * P1x * dP0y_exp * P6z);
    double dterm_1_dz = (49140.0 * P5x * P2y * dP0z_exp) - (2520.0 * P7x * dP0z_exp) + (3780.0 * P5x * dP2z_exp) - (70875.0 * P3x * P4y * dP0z_exp) - (66150.0 * P3x * P2y * dP2z_exp) + (4725.0 * P3x * dP4z_exp) + (12600.0 * P1x * P6y * dP0z_exp) + (23625.0 * P1x * P4y * dP2z_exp) + (9450.0 * P1x * P2y * dP4z_exp) - (1575.0 * P1x * dP6z_exp);

    double dterm_2_dx = (49140.0 * dP2x_exp * P5y) - (2520.0 * dP0x_exp * P7y) + (3780.0 * dP0x_exp * P5y * P2z) - (70875.0 * dP4x_exp * P3y) - (66150.0 * dP2x_exp * P3y * P2z) + (4725.0 * dP0x_exp * P3y * P4z) + (12600.0 * dP6x_exp * P1y) + (23625.0 * dP4x_exp * P1y * P2z) + (9450.0 * dP2x_exp * P1y * P4z) - (1575.0 * dP0x_exp * P1y * P6z);
    double dterm_2_dy = (49140.0 * P2x * dP5y_exp) - (2520.0 * dP7y_exp) + (3780.0 * dP5y_exp * P2z) - (70875.0 * P4x * dP3y_exp) - (66150.0 * P2x * dP3y_exp * P2z) + (4725.0 * dP3y_exp * P4z) + (12600.0 * P6x * dP1y_exp) + (23625.0 * P4x * dP1y_exp * P2z) + (9450.0 * P2x * dP1y_exp * P4z) - (1575.0 * dP1y_exp * P6z);
    double dterm_2_dz = (49140.0 * P2x * P5y * dP0z_exp) - (2520.0 * P7y * dP0z_exp) + (3780.0 * P5y * dP2z_exp) - (70875.0 * P4x * P3y * dP0z_exp) - (66150.0 * P2x * P3y * dP2z_exp) + (4725.0 * P3y * dP4z_exp) + (12600.0 * P6x * P1y * dP0z_exp) + (23625.0 * P4x * P1y * dP2z_exp) + (9450.0 * P2x * P1y * dP4z_exp) - (1575.0 * P1y * dP6z_exp);

    double dterm_3_dx = (49140.0 * dP5x_exp * P2z) - (2520.0 * dP7x_exp) + (3780.0 * dP5x_exp * P2y) - (70875.0 * dP3x_exp * P4z) - (66150.0 * dP3x_exp * P2y * P2z) + (4725.0 * dP3x_exp * P4y) + (12600.0 * dP1x_exp * P6z) + (23625.0 * dP1x_exp * P2y * P4z) + (9450.0 * dP1x_exp * P4y * P2z) - (1575.0 * dP1x_exp * P6y);
    double dterm_3_dy = (49140.0 * P5x * dP0y_exp * P2z) - (2520.0 * P7x * dP0y_exp) + (3780.0 * P5x * dP2y_exp) - (70875.0 * P3x * dP0y_exp * P4z) - (66150.0 * P3x * dP2y_exp * P2z) + (4725.0 * P3x * dP4y_exp) + (12600.0 * P1x * dP0y_exp * P6z) + (23625.0 * P1x * dP2y_exp * P4z) + (9450.0 * P1x * dP4y_exp * P2z) - (1575.0 * P1x * dP6y_exp);
    double dterm_3_dz = (49140.0 * P5x * dP2z_exp) - (2520.0 * P7x * dP0z_exp) + (3780.0 * P5x * P2y * dP0z_exp) - (70875.0 * P3x * dP4z_exp) - (66150.0 * P3x * P2y * dP2z_exp) + (4725.0 * P3x * P4y * dP0z_exp) + (12600.0 * P1x * dP6z_exp) + (23625.0 * P1x * P2y * dP4z_exp) + (9450.0 * P1x * P4y * dP2z_exp) - (1575.0 * P1x * P6y * dP0z_exp);

    double dterm_4_dx = (49140.0 * dP2x_exp * P5z) - (2520.0 * dP0x_exp * P7z) + (3780.0 * dP0x_exp * P2y * P5z) - (70875.0 * dP4x_exp * P3z) - (66150.0 * dP2x_exp * P2y * P3z) + (4725.0 * dP0x_exp * P4y * P3z) + (12600.0 * dP6x_exp * P1z) + (23625.0 * dP4x_exp * P2y * P1z) + (9450.0 * dP2x_exp * P4y * P1z) - (1575.0 * dP0x_exp * P6y * P1z);
    double dterm_4_dy = (49140.0 * P2x * dP0y_exp * P5z) - (2520.0 * dP0y_exp * P7z) + (3780.0 * dP2y_exp * P5z) - (70875.0 * P4x * dP0y_exp * P3z) - (66150.0 * P2x * dP2y_exp * P3z) + (4725.0 * dP4y_exp * P3z) + (12600.0 * P6x * dP0y_exp * P1z) + (23625.0 * P4x * dP2y_exp * P1z) + (9450.0 * P2x * dP4y_exp * P1z) - (1575.0 * dP6y_exp * P1z);
    double dterm_4_dz = (49140.0 * P2x * dP5z_exp) - (2520.0 * dP7z_exp) + (3780.0 * P2y * dP5z_exp) - (70875.0 * P4x * dP3z_exp) - (66150.0 * P2x * P2y * dP3z_exp) + (4725.0 * P4y * dP3z_exp) + (12600.0 * P6x * dP1z_exp) + (23625.0 * P4x * P2y * dP1z_exp) + (9450.0 * P2x * P4y * dP1z_exp) - (1575.0 * P6y * dP1z_exp);

    double dterm_5_dx = (49140.0 * dP0x_exp * P5y * P2z) - (2520.0 * dP0x_exp * P7y) + (3780.0 * dP2x_exp * P5y) - (70875.0 * dP0x_exp * P3y * P4z) - (66150.0 * dP2x_exp * P3y * P2z) + (4725.0 * dP4x_exp * P3y) + (12600.0 * dP0x_exp * P1y * P6z) + (23625.0 * dP2x_exp * P1y * P4z) + (9450.0 * dP4x_exp * P1y * P2z) - (1575.0 * dP6x_exp * P1y);
    double dterm_5_dy = (49140.0 * dP5y_exp * P2z) - (2520.0 * dP7y_exp) + (3780.0 * P2x * dP5y_exp) - (70875.0 * dP3y_exp * P4z) - (66150.0 * P2x * dP3y_exp * P2z) + (4725.0 * P4x * dP3y_exp) + (12600.0 * dP1y_exp * P6z) + (23625.0 * P2x * dP1y_exp * P4z) + (9450.0 * P4x * dP1y_exp * P2z) - (1575.0 * P6x * dP1y_exp);
    double dterm_5_dz = (49140.0 * P5y * dP2z_exp) - (2520.0 * P7y * dP0z_exp) + (3780.0 * P2x * P5y * dP0z_exp) - (70875.0 * P3y * dP4z_exp) - (66150.0 * P2x * P3y * dP2z_exp) + (4725.0 * P4x * P3y * dP0z_exp) + (12600.0 * P1y * dP6z_exp) + (23625.0 * P2x * P1y * dP4z_exp) + (9450.0 * P4x * P1y * dP2z_exp) - (1575.0 * P6x * P1y * dP0z_exp);

    double dterm_6_dx = (49140.0 * dP0x_exp * P2y * P5z) - (2520.0 * dP0x_exp * P7z) + (3780.0 * dP2x_exp * P5z) - (70875.0 * dP0x_exp * P4y * P3z) - (66150.0 * dP2x_exp * P2y * P3z) + (4725.0 * dP4x_exp * P3z) + (12600.0 * dP0x_exp * P6y * P1z) + (23625.0 * dP2x_exp * P4y * P1z) + (9450.0 * dP4x_exp * P2y * P1z) - (1575.0 * dP6x_exp * P1z);
    double dterm_6_dy = (49140.0 * dP2y_exp * P5z) - (2520.0 * dP0y_exp * P7z) + (3780.0 * P2x * dP0y_exp * P5z) - (70875.0 * dP4y_exp * P3z) - (66150.0 * P2x * dP2y_exp * P3z) + (4725.0 * P4x * dP0y_exp * P3z) + (12600.0 * dP6y_exp * P1z) + (23625.0 * P2x * dP4y_exp * P1z) + (9450.0 * P4x * dP2y_exp * P1z) - (1575.0 * P6x * dP0y_exp * P1z);
    double dterm_6_dz = (49140.0 * P2y * dP5z_exp) - (2520.0 * dP7z_exp) + (3780.0 * P2x * dP5z_exp) - (70875.0 * P4y * dP3z_exp) - (66150.0 * P2x * P2y * dP3z_exp) + (4725.0 * P4x * dP3z_exp) + (12600.0 * P6y * dP1z_exp) + (23625.0 * P2x * P4y * dP1z_exp) + (9450.0 * P4x * P2y * dP1z_exp) - (1575.0 * P6x * dP1z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

    // dmiu4 dx/dy/dz
    deriv[9] = temp * dterm_4_dx;
    deriv[10] = temp * dterm_4_dy;
    deriv[11] = temp * dterm_4_dz;

    // dmiu5 dx/dy/dz
    deriv[12] = temp * dterm_5_dx;
    deriv[13] = temp * dterm_5_dy;
    deriv[14] = temp * dterm_5_dz;

    // dmiu6 dx/dy/dz
    deriv[15] = temp * dterm_6_dx;
    deriv[16] = temp * dterm_6_dy;
    deriv[17] = temp * dterm_6_dz;

}

void calc_solid_MCSH_7_4(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

    double P1x = P1(lambda, x0, gamma);
    double P1y = P1(lambda, y0, gamma);
    double P1z = P1(lambda, z0, gamma);

    double P3x = P3(lambda, x0, gamma);
    double P3y = P3(lambda, y0, gamma);
    double P3z = P3(lambda, z0, gamma);

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double term_1 = (45360.0 * P5x * P1y * P1z) - (75600.0 * P3x * P3y * P1z) - (75600.0 * P3x * P1y * P3z) + (14175.0 * P1x * P5y * P1z) + (28350.0 * P1x * P3y * P3z) + (14175.0 * P1x * P1y * P5z);
    double term_2 = (45360.0 * P1x * P5y * P1z) - (75600.0 * P3x * P3y * P1z) - (75600.0 * P1x * P3y * P3z) + (14175.0 * P5x * P1y * P1z) + (28350.0 * P3x * P1y * P3z) + (14175.0 * P1x * P1y * P5z);
    double term_3 = (45360.0 * P1x * P1y * P5z) - (75600.0 * P3x * P1y * P3z) - (75600.0 * P1x * P3y * P3z) + (14175.0 * P5x * P1y * P1z) + (28350.0 * P3x * P3y * P1z) + (14175.0 * P1x * P5y * P1z);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dP3x_exp = dP3_exp(P3x, C2, lambda, x0, gamma);
    double dP3y_exp = dP3_exp(P3y, C2, lambda, y0, gamma);
    double dP3z_exp = dP3_exp(P3z, C2, lambda, z0, gamma);

    double dP5x_exp = dP5_exp(P5x, C2, lambda, x0, gamma);
    double dP5y_exp = dP5_exp(P5y, C2, lambda, y0, gamma);
    double dP5z_exp = dP5_exp(P5z, C2, lambda, z0, gamma);

    double dterm_1_dx = (45360.0 * dP5x_exp * P1y * P1z) - (75600.0 * dP3x_exp * P3y * P1z) - (75600.0 * dP3x_exp * P1y * P3z) + (14175.0 * dP1x_exp * P5y * P1z) + (28350.0 * dP1x_exp * P3y * P3z) + (14175.0 * dP1x_exp * P1y * P5z);
    double dterm_1_dy = (45360.0 * P5x * dP1y_exp * P1z) - (75600.0 * P3x * dP3y_exp * P1z) - (75600.0 * P3x * dP1y_exp * P3z) + (14175.0 * P1x * dP5y_exp * P1z) + (28350.0 * P1x * dP3y_exp * P3z) + (14175.0 * P1x * dP1y_exp * P5z);
    double dterm_1_dz = (45360.0 * P5x * P1y * dP1z_exp) - (75600.0 * P3x * P3y * dP1z_exp) - (75600.0 * P3x * P1y * dP3z_exp) + (14175.0 * P1x * P5y * dP1z_exp) + (28350.0 * P1x * P3y * dP3z_exp) + (14175.0 * P1x * P1y * dP5z_exp);

    double dterm_2_dx = (45360.0 * dP1x_exp * P5y * P1z) - (75600.0 * dP3x_exp * P3y * P1z) - (75600.0 * dP1x_exp * P3y * P3z) + (14175.0 * dP5x_exp * P1y * P1z) + (28350.0 * dP3x_exp * P1y * P3z) + (14175.0 * dP1x_exp * P1y * P5z);
    double dterm_2_dy = (45360.0 * P1x * dP5y_exp * P1z) - (75600.0 * P3x * dP3y_exp * P1z) - (75600.0 * P1x * dP3y_exp * P3z) + (14175.0 * P5x * dP1y_exp * P1z) + (28350.0 * P3x * dP1y_exp * P3z) + (14175.0 * P1x * dP1y_exp * P5z);
    double dterm_2_dz = (45360.0 * P1x * P5y * dP1z_exp) - (75600.0 * P3x * P3y * dP1z_exp) - (75600.0 * P1x * P3y * dP3z_exp) + (14175.0 * P5x * P1y * dP1z_exp) + (28350.0 * P3x * P1y * dP3z_exp) + (14175.0 * P1x * P1y * dP5z_exp);

    double dterm_3_dx = (45360.0 * dP1x_exp * P1y * P5z) - (75600.0 * dP3x_exp * P1y * P3z) - (75600.0 * dP1x_exp * P3y * P3z) + (14175.0 * dP5x_exp * P1y * P1z) + (28350.0 * dP3x_exp * P3y * P1z) + (14175.0 * dP1x_exp * P5y * P1z);
    double dterm_3_dy = (45360.0 * P1x * dP1y_exp * P5z) - (75600.0 * P3x * dP1y_exp * P3z) - (75600.0 * P1x * dP3y_exp * P3z) + (14175.0 * P5x * dP1y_exp * P1z) + (28350.0 * P3x * dP3y_exp * P1z) + (14175.0 * P1x * dP5y_exp * P1z);
    double dterm_3_dz = (45360.0 * P1x * P1y * dP5z_exp) - (75600.0 * P3x * P1y * dP3z_exp) - (75600.0 * P1x * P3y * dP3z_exp) + (14175.0 * P5x * P1y * dP1z_exp) + (28350.0 * P3x * P3y * dP1z_exp) + (14175.0 * P1x * P5y * dP1z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

}

void calc_solid_MCSH_7_5(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double term_1 = (75600.0 * P4x * P3y) - (15120.0 * P6x * P1y) - (42525.0 * P2x * P5y) - (28350.0 * P2x * P3y * P2z) + (14175.0 * P2x * P1y * P4z) + (1890.0 * P7y) + (2835.0 * P5y * P2z) - (945.0 * P1y * P6z);
    double term_2 = (75600.0 * P3x * P4y) - (15120.0 * P1x * P6y) - (42525.0 * P5x * P2y) - (28350.0 * P3x * P2y * P2z) + (14175.0 * P1x * P2y * P4z) + (1890.0 * P7x) + (2835.0 * P5x * P2z) - (945.0 * P1x * P6z);
    double term_3 = (75600.0 * P4x * P3z) - (15120.0 * P6x * P1z) - (42525.0 * P2x * P5z) - (28350.0 * P2x * P2y * P3z) + (14175.0 * P2x * P4y * P1z) + (1890.0 * P7z) + (2835.0 * P2y * P5z) - (945.0 * P6y * P1z);
    double term_4 = (75600.0 * P3x * P4z) - (15120.0 * P1x * P6z) - (42525.0 * P5x * P2z) - (28350.0 * P3x * P2y * P2z) + (14175.0 * P1x * P4y * P2z) + (1890.0 * P7x) + (2835.0 * P5x * P2y) - (945.0 * P1x * P6y);
    double term_5 = (75600.0 * P4y * P3z) - (15120.0 * P6y * P1z) - (42525.0 * P2y * P5z) - (28350.0 * P2x * P2y * P3z) + (14175.0 * P4x * P2y * P1z) + (1890.0 * P7z) + (2835.0 * P2x * P5z) - (945.0 * P6x * P1z);
    double term_6 = (75600.0 * P3y * P4z) - (15120.0 * P1y * P6z) - (42525.0 * P5y * P2z) - (28350.0 * P2x * P3y * P2z) + (14175.0 * P4x * P1y * P2z) + (1890.0 * P7y) + (2835.0 * P2x * P5y) - (945.0 * P6x * P1y);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP3x_exp = dP3_exp(P3x, C2, lambda, x0, gamma);
    double dP3y_exp = dP3_exp(P3y, C2, lambda, y0, gamma);
    double dP3z_exp = dP3_exp(P3z, C2, lambda, z0, gamma);

    double dP4x_exp = dP4_exp(P4x, C2, lambda, x0, gamma);
    double dP4y_exp = dP4_exp(P4y, C2, lambda, y0, gamma);
    double dP4z_exp = dP4_exp(P4z, C2, lambda, z0, gamma);

    double dP5x_exp = dP5_exp(P5x, C2, lambda, x0, gamma);
    double dP5y_exp = dP5_exp(P5y, C2, lambda, y0, gamma);
    double dP5z_exp = dP5_exp(P5z, C2, lambda, z0, gamma);

    double dP6x_exp = dP6_exp(P6x, C2, lambda, x0, gamma);
    double dP6y_exp = dP6_exp(P6y, C2, lambda, y0, gamma);
    double dP6z_exp = dP6_exp(P6z, C2, lambda, z0, gamma);

    double dP7x_exp = dP7_exp(P7x, C2, lambda, x0, gamma);
    double dP7y_exp = dP7_exp(P7y, C2, lambda, y0, gamma);
    double dP7z_exp = dP7_exp(P7z, C2, lambda, z0, gamma);

    double dterm_1_dx = (75600.0 * dP4x_exp * P3y) - (15120.0 * dP6x_exp * P1y) - (42525.0 * dP2x_exp * P5y) - (28350.0 * dP2x_exp * P3y * P2z) + (14175.0 * dP2x_exp * P1y * P4z) + (1890.0 * dP0x_exp * P7y) + (2835.0 * dP0x_exp * P5y * P2z) - (945.0 * dP0x_exp * P1y * P6z);
    double dterm_1_dy = (75600.0 * P4x * dP3y_exp) - (15120.0 * P6x * dP1y_exp) - (42525.0 * P2x * dP5y_exp) - (28350.0 * P2x * dP3y_exp * P2z) + (14175.0 * P2x * dP1y_exp * P4z) + (1890.0 * dP7y_exp) + (2835.0 * dP5y_exp * P2z) - (945.0 * dP1y_exp * P6z);
    double dterm_1_dz = (75600.0 * P4x * P3y * dP0z_exp) - (15120.0 * P6x * P1y * dP0z_exp) - (42525.0 * P2x * P5y * dP0z_exp) - (28350.0 * P2x * P3y * dP2z_exp) + (14175.0 * P2x * P1y * dP4z_exp) + (1890.0 * P7y * dP0z_exp) + (2835.0 * P5y * dP2z_exp) - (945.0 * P1y * dP6z_exp);

    double dterm_2_dx = (75600.0 * dP3x_exp * P4y) - (15120.0 * dP1x_exp * P6y) - (42525.0 * dP5x_exp * P2y) - (28350.0 * dP3x_exp * P2y * P2z) + (14175.0 * dP1x_exp * P2y * P4z) + (1890.0 * dP7x_exp) + (2835.0 * dP5x_exp * P2z) - (945.0 * dP1x_exp * P6z);
    double dterm_2_dy = (75600.0 * P3x * dP4y_exp) - (15120.0 * P1x * dP6y_exp) - (42525.0 * P5x * dP2y_exp) - (28350.0 * P3x * dP2y_exp * P2z) + (14175.0 * P1x * dP2y_exp * P4z) + (1890.0 * P7x * dP0y_exp) + (2835.0 * P5x * dP0y_exp * P2z) - (945.0 * P1x * dP0y_exp * P6z);
    double dterm_2_dz = (75600.0 * P3x * P4y * dP0z_exp) - (15120.0 * P1x * P6y * dP0z_exp) - (42525.0 * P5x * P2y * dP0z_exp) - (28350.0 * P3x * P2y * dP2z_exp) + (14175.0 * P1x * P2y * dP4z_exp) + (1890.0 * P7x * dP0z_exp) + (2835.0 * P5x * dP2z_exp) - (945.0 * P1x * dP6z_exp);

    double dterm_3_dx = (75600.0 * dP4x_exp * P3z) - (15120.0 * dP6x_exp * P1z) - (42525.0 * dP2x_exp * P5z) - (28350.0 * dP2x_exp * P2y * P3z) + (14175.0 * dP2x_exp * P4y * P1z) + (1890.0 * dP0x_exp * P7z) + (2835.0 * dP0x_exp * P2y * P5z) - (945.0 * dP0x_exp * P6y * P1z);
    double dterm_3_dy = (75600.0 * P4x * dP0y_exp * P3z) - (15120.0 * P6x * dP0y_exp * P1z) - (42525.0 * P2x * dP0y_exp * P5z) - (28350.0 * P2x * dP2y_exp * P3z) + (14175.0 * P2x * dP4y_exp * P1z) + (1890.0 * dP0y_exp * P7z) + (2835.0 * dP2y_exp * P5z) - (945.0 * dP6y_exp * P1z);
    double dterm_3_dz = (75600.0 * P4x * dP3z_exp) - (15120.0 * P6x * dP1z_exp) - (42525.0 * P2x * dP5z_exp) - (28350.0 * P2x * P2y * dP3z_exp) + (14175.0 * P2x * P4y * dP1z_exp) + (1890.0 * dP7z_exp) + (2835.0 * P2y * dP5z_exp) - (945.0 * P6y * dP1z_exp);

    double dterm_4_dx = (75600.0 * dP3x_exp * P4z) - (15120.0 * dP1x_exp * P6z) - (42525.0 * dP5x_exp * P2z) - (28350.0 * dP3x_exp * P2y * P2z) + (14175.0 * dP1x_exp * P4y * P2z) + (1890.0 * dP7x_exp) + (2835.0 * dP5x_exp * P2y) - (945.0 * dP1x_exp * P6y);
    double dterm_4_dy = (75600.0 * P3x * dP0y_exp * P4z) - (15120.0 * P1x * dP0y_exp * P6z) - (42525.0 * P5x * dP0y_exp * P2z) - (28350.0 * P3x * dP2y_exp * P2z) + (14175.0 * P1x * dP4y_exp * P2z) + (1890.0 * P7x * dP0y_exp) + (2835.0 * P5x * dP2y_exp) - (945.0 * P1x * dP6y_exp);
    double dterm_4_dz = (75600.0 * P3x * dP4z_exp) - (15120.0 * P1x * dP6z_exp) - (42525.0 * P5x * dP2z_exp) - (28350.0 * P3x * P2y * dP2z_exp) + (14175.0 * P1x * P4y * dP2z_exp) + (1890.0 * P7x * dP0z_exp) + (2835.0 * P5x * P2y * dP0z_exp) - (945.0 * P1x * P6y * dP0z_exp);

    double dterm_5_dx = (75600.0 * dP0x_exp * P4y * P3z) - (15120.0 * dP0x_exp * P6y * P1z) - (42525.0 * dP0x_exp * P2y * P5z) - (28350.0 * dP2x_exp * P2y * P3z) + (14175.0 * dP4x_exp * P2y * P1z) + (1890.0 * dP0x_exp * P7z) + (2835.0 * dP2x_exp * P5z) - (945.0 * dP6x_exp * P1z);
    double dterm_5_dy = (75600.0 * dP4y_exp * P3z) - (15120.0 * dP6y_exp * P1z) - (42525.0 * dP2y_exp * P5z) - (28350.0 * P2x * dP2y_exp * P3z) + (14175.0 * P4x * dP2y_exp * P1z) + (1890.0 * dP0y_exp * P7z) + (2835.0 * P2x * dP0y_exp * P5z) - (945.0 * P6x * dP0y_exp * P1z);
    double dterm_5_dz = (75600.0 * P4y * dP3z_exp) - (15120.0 * P6y * dP1z_exp) - (42525.0 * P2y * dP5z_exp) - (28350.0 * P2x * P2y * dP3z_exp) + (14175.0 * P4x * P2y * dP1z_exp) + (1890.0 * dP7z_exp) + (2835.0 * P2x * dP5z_exp) - (945.0 * P6x * dP1z_exp);

    double dterm_6_dx = (75600.0 * dP0x_exp * P3y * P4z) - (15120.0 * dP0x_exp * P1y * P6z) - (42525.0 * dP0x_exp * P5y * P2z) - (28350.0 * dP2x_exp * P3y * P2z) + (14175.0 * dP4x_exp * P1y * P2z) + (1890.0 * dP0x_exp * P7y) + (2835.0 * dP2x_exp * P5y) - (945.0 * dP6x_exp * P1y);
    double dterm_6_dy = (75600.0 * dP3y_exp * P4z) - (15120.0 * dP1y_exp * P6z) - (42525.0 * dP5y_exp * P2z) - (28350.0 * P2x * dP3y_exp * P2z) + (14175.0 * P4x * dP1y_exp * P2z) + (1890.0 * dP7y_exp) + (2835.0 * P2x * dP5y_exp) - (945.0 * P6x * dP1y_exp);
    double dterm_6_dz = (75600.0 * P3y * dP4z_exp) - (15120.0 * P1y * dP6z_exp) - (42525.0 * P5y * dP2z_exp) - (28350.0 * P2x * P3y * dP2z_exp) + (14175.0 * P4x * P1y * dP2z_exp) + (1890.0 * P7y * dP0z_exp) + (2835.0 * P2x * P5y * dP0z_exp) - (945.0 * P6x * P1y * dP0z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

    // dmiu4 dx/dy/dz
    deriv[9] = temp * dterm_4_dx;
    deriv[10] = temp * dterm_4_dy;
    deriv[11] = temp * dterm_4_dz;

    // dmiu5 dx/dy/dz
    deriv[12] = temp * dterm_5_dx;
    deriv[13] = temp * dterm_5_dy;
    deriv[14] = temp * dterm_5_dz;

    // dmiu6 dx/dy/dz
    deriv[15] = temp * dterm_6_dx;
    deriv[16] = temp * dterm_6_dy;
    deriv[17] = temp * dterm_6_dz;

}

void calc_solid_MCSH_7_6(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double term_1 = (75600.0 * P4x * P2y * P1z) - (5040.0 * P6x * P1z) - (51975.0 * P2x * P4y * P1z) - (47250.0 * P2x * P2y * P3z) + (4725.0 * P2x * P5z) + (2520.0 * P6y * P1z) + (4725.0 * P4y * P3z) + (1890.0 * P2y * P5z) - (315.0 * P7z);
    double term_2 = (75600.0 * P2x * P4y * P1z) - (5040.0 * P6y * P1z) - (51975.0 * P4x * P2y * P1z) - (47250.0 * P2x * P2y * P3z) + (4725.0 * P2y * P5z) + (2520.0 * P6x * P1z) + (4725.0 * P4x * P3z) + (1890.0 * P2x * P5z) - (315.0 * P7z);
    double term_3 = (75600.0 * P4x * P1y * P2z) - (5040.0 * P6x * P1y) - (51975.0 * P2x * P1y * P4z) - (47250.0 * P2x * P3y * P2z) + (4725.0 * P2x * P5y) + (2520.0 * P1y * P6z) + (4725.0 * P3y * P4z) + (1890.0 * P5y * P2z) - (315.0 * P7y);
    double term_4 = (75600.0 * P2x * P1y * P4z) - (5040.0 * P1y * P6z) - (51975.0 * P4x * P1y * P2z) - (47250.0 * P2x * P3y * P2z) + (4725.0 * P5y * P2z) + (2520.0 * P6x * P1y) + (4725.0 * P4x * P3y) + (1890.0 * P2x * P5y) - (315.0 * P7y);
    double term_5 = (75600.0 * P1x * P4y * P2z) - (5040.0 * P1x * P6y) - (51975.0 * P1x * P2y * P4z) - (47250.0 * P3x * P2y * P2z) + (4725.0 * P5x * P2y) + (2520.0 * P1x * P6z) + (4725.0 * P3x * P4z) + (1890.0 * P5x * P2z) - (315.0 * P7x);
    double term_6 = (75600.0 * P1x * P2y * P4z) - (5040.0 * P1x * P6z) - (51975.0 * P1x * P4y * P2z) - (47250.0 * P3x * P2y * P2z) + (4725.0 * P5x * P2z) + (2520.0 * P1x * P6y) + (4725.0 * P3x * P4y) + (1890.0 * P5x * P2y) - (315.0 * P7x);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP3x_exp = dP3_exp(P3x, C2, lambda, x0, gamma);
    double dP3y_exp = dP3_exp(P3y, C2, lambda, y0, gamma);
    double dP3z_exp = dP3_exp(P3z, C2, lambda, z0, gamma);

    double dP4x_exp = dP4_exp(P4x, C2, lambda, x0, gamma);
    double dP4y_exp = dP4_exp(P4y, C2, lambda, y0, gamma);
    double dP4z_exp = dP4_exp(P4z, C2, lambda, z0, gamma);

    double dP5x_exp = dP5_exp(P5x, C2, lambda, x0, gamma);
    double dP5y_exp = dP5_exp(P5y, C2, lambda, y0, gamma);
    double dP5z_exp = dP5_exp(P5z, C2, lambda, z0, gamma);

    double dP6x_exp = dP6_exp(P6x, C2, lambda, x0, gamma);
    double dP6y_exp = dP6_exp(P6y, C2, lambda, y0, gamma);
    double dP6z_exp = dP6_exp(P6z, C2, lambda, z0, gamma);

    double dP7x_exp = dP7_exp(P7x, C2, lambda, x0, gamma);
    double dP7y_exp = dP7_exp(P7y, C2, lambda, y0, gamma);
    double dP7z_exp = dP7_exp(P7z, C2, lambda, z0, gamma);

    double dterm_1_dx = (75600.0 * dP4x_exp * P2y * P1z) - (5040.0 * dP6x_exp * P1z) - (51975.0 * dP2x_exp * P4y * P1z) - (47250.0 * dP2x_exp * P2y * P3z) + (4725.0 * dP2x_exp * P5z) + (2520.0 * dP0x_exp * P6y * P1z) + (4725.0 * dP0x_exp * P4y * P3z) + (1890.0 * dP0x_exp * P2y * P5z) - (315.0 * dP0x_exp * P7z);
    double dterm_1_dy = (75600.0 * P4x * dP2y_exp * P1z) - (5040.0 * P6x * dP0y_exp * P1z) - (51975.0 * P2x * dP4y_exp * P1z) - (47250.0 * P2x * dP2y_exp * P3z) + (4725.0 * P2x * dP0y_exp * P5z) + (2520.0 * dP6y_exp * P1z) + (4725.0 * dP4y_exp * P3z) + (1890.0 * dP2y_exp * P5z) - (315.0 * dP0y_exp * P7z);
    double dterm_1_dz = (75600.0 * P4x * P2y * dP1z_exp) - (5040.0 * P6x * dP1z_exp) - (51975.0 * P2x * P4y * dP1z_exp) - (47250.0 * P2x * P2y * dP3z_exp) + (4725.0 * P2x * dP5z_exp) + (2520.0 * P6y * dP1z_exp) + (4725.0 * P4y * dP3z_exp) + (1890.0 * P2y * dP5z_exp) - (315.0 * dP7z_exp);

    double dterm_2_dx = (75600.0 * dP2x_exp * P4y * P1z) - (5040.0 * dP0x_exp * P6y * P1z) - (51975.0 * dP4x_exp * P2y * P1z) - (47250.0 * dP2x_exp * P2y * P3z) + (4725.0 * dP0x_exp * P2y * P5z) + (2520.0 * dP6x_exp * P1z) + (4725.0 * dP4x_exp * P3z) + (1890.0 * dP2x_exp * P5z) - (315.0 * dP0x_exp * P7z);
    double dterm_2_dy = (75600.0 * P2x * dP4y_exp * P1z) - (5040.0 * dP6y_exp * P1z) - (51975.0 * P4x * dP2y_exp * P1z) - (47250.0 * P2x * dP2y_exp * P3z) + (4725.0 * dP2y_exp * P5z) + (2520.0 * P6x * dP0y_exp * P1z) + (4725.0 * P4x * dP0y_exp * P3z) + (1890.0 * P2x * dP0y_exp * P5z) - (315.0 * dP0y_exp * P7z);
    double dterm_2_dz = (75600.0 * P2x * P4y * dP1z_exp) - (5040.0 * P6y * dP1z_exp) - (51975.0 * P4x * P2y * dP1z_exp) - (47250.0 * P2x * P2y * dP3z_exp) + (4725.0 * P2y * dP5z_exp) + (2520.0 * P6x * dP1z_exp) + (4725.0 * P4x * dP3z_exp) + (1890.0 * P2x * dP5z_exp) - (315.0 * dP7z_exp);

    double dterm_3_dx = (75600.0 * dP4x_exp * P1y * P2z) - (5040.0 * dP6x_exp * P1y) - (51975.0 * dP2x_exp * P1y * P4z) - (47250.0 * dP2x_exp * P3y * P2z) + (4725.0 * dP2x_exp * P5y) + (2520.0 * dP0x_exp * P1y * P6z) + (4725.0 * dP0x_exp * P3y * P4z) + (1890.0 * dP0x_exp * P5y * P2z) - (315.0 * dP0x_exp * P7y);
    double dterm_3_dy = (75600.0 * P4x * dP1y_exp * P2z) - (5040.0 * P6x * dP1y_exp) - (51975.0 * P2x * dP1y_exp * P4z) - (47250.0 * P2x * dP3y_exp * P2z) + (4725.0 * P2x * dP5y_exp) + (2520.0 * dP1y_exp * P6z) + (4725.0 * dP3y_exp * P4z) + (1890.0 * dP5y_exp * P2z) - (315.0 * dP7y_exp);
    double dterm_3_dz = (75600.0 * P4x * P1y * dP2z_exp) - (5040.0 * P6x * P1y * dP0z_exp) - (51975.0 * P2x * P1y * dP4z_exp) - (47250.0 * P2x * P3y * dP2z_exp) + (4725.0 * P2x * P5y * dP0z_exp) + (2520.0 * P1y * dP6z_exp) + (4725.0 * P3y * dP4z_exp) + (1890.0 * P5y * dP2z_exp) - (315.0 * P7y * dP0z_exp);

    double dterm_4_dx = (75600.0 * dP2x_exp * P1y * P4z) - (5040.0 * dP0x_exp * P1y * P6z) - (51975.0 * dP4x_exp * P1y * P2z) - (47250.0 * dP2x_exp * P3y * P2z) + (4725.0 * dP0x_exp * P5y * P2z) + (2520.0 * dP6x_exp * P1y) + (4725.0 * dP4x_exp * P3y) + (1890.0 * dP2x_exp * P5y) - (315.0 * dP0x_exp * P7y);
    double dterm_4_dy = (75600.0 * P2x * dP1y_exp * P4z) - (5040.0 * dP1y_exp * P6z) - (51975.0 * P4x * dP1y_exp * P2z) - (47250.0 * P2x * dP3y_exp * P2z) + (4725.0 * dP5y_exp * P2z) + (2520.0 * P6x * dP1y_exp) + (4725.0 * P4x * dP3y_exp) + (1890.0 * P2x * dP5y_exp) - (315.0 * dP7y_exp);
    double dterm_4_dz = (75600.0 * P2x * P1y * dP4z_exp) - (5040.0 * P1y * dP6z_exp) - (51975.0 * P4x * P1y * dP2z_exp) - (47250.0 * P2x * P3y * dP2z_exp) + (4725.0 * P5y * dP2z_exp) + (2520.0 * P6x * P1y * dP0z_exp) + (4725.0 * P4x * P3y * dP0z_exp) + (1890.0 * P2x * P5y * dP0z_exp) - (315.0 * P7y * dP0z_exp);

    double dterm_5_dx = (75600.0 * dP1x_exp * P4y * P2z) - (5040.0 * dP1x_exp * P6y) - (51975.0 * dP1x_exp * P2y * P4z) - (47250.0 * dP3x_exp * P2y * P2z) + (4725.0 * dP5x_exp * P2y) + (2520.0 * dP1x_exp * P6z) + (4725.0 * dP3x_exp * P4z) + (1890.0 * dP5x_exp * P2z) - (315.0 * dP7x_exp);
    double dterm_5_dy = (75600.0 * P1x * dP4y_exp * P2z) - (5040.0 * P1x * dP6y_exp) - (51975.0 * P1x * dP2y_exp * P4z) - (47250.0 * P3x * dP2y_exp * P2z) + (4725.0 * P5x * dP2y_exp) + (2520.0 * P1x * dP0y_exp * P6z) + (4725.0 * P3x * dP0y_exp * P4z) + (1890.0 * P5x * dP0y_exp * P2z) - (315.0 * P7x * dP0y_exp);
    double dterm_5_dz = (75600.0 * P1x * P4y * dP2z_exp) - (5040.0 * P1x * P6y * dP0z_exp) - (51975.0 * P1x * P2y * dP4z_exp) - (47250.0 * P3x * P2y * dP2z_exp) + (4725.0 * P5x * P2y * dP0z_exp) + (2520.0 * P1x * dP6z_exp) + (4725.0 * P3x * dP4z_exp) + (1890.0 * P5x * dP2z_exp) - (315.0 * P7x * dP0z_exp);

    double dterm_6_dx = (75600.0 * dP1x_exp * P2y * P4z) - (5040.0 * dP1x_exp * P6z) - (51975.0 * dP1x_exp * P4y * P2z) - (47250.0 * dP3x_exp * P2y * P2z) + (4725.0 * dP5x_exp * P2z) + (2520.0 * dP1x_exp * P6y) + (4725.0 * dP3x_exp * P4y) + (1890.0 * dP5x_exp * P2y) - (315.0 * dP7x_exp);
    double dterm_6_dy = (75600.0 * P1x * dP2y_exp * P4z) - (5040.0 * P1x * dP0y_exp * P6z) - (51975.0 * P1x * dP4y_exp * P2z) - (47250.0 * P3x * dP2y_exp * P2z) + (4725.0 * P5x * dP0y_exp * P2z) + (2520.0 * P1x * dP6y_exp) + (4725.0 * P3x * dP4y_exp) + (1890.0 * P5x * dP2y_exp) - (315.0 * P7x * dP0y_exp);
    double dterm_6_dz = (75600.0 * P1x * P2y * dP4z_exp) - (5040.0 * P1x * dP6z_exp) - (51975.0 * P1x * P4y * dP2z_exp) - (47250.0 * P3x * P2y * dP2z_exp) + (4725.0 * P5x * dP2z_exp) + (2520.0 * P1x * P6y * dP0z_exp) + (4725.0 * P3x * P4y * dP0z_exp) + (1890.0 * P5x * P2y * dP0z_exp) - (315.0 * P7x * dP0z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

    // dmiu4 dx/dy/dz
    deriv[9] = temp * dterm_4_dx;
    deriv[10] = temp * dterm_4_dy;
    deriv[11] = temp * dterm_4_dz;

    // dmiu5 dx/dy/dz
    deriv[12] = temp * dterm_5_dx;
    deriv[13] = temp * dterm_5_dy;
    deriv[14] = temp * dterm_5_dz;

    // dmiu6 dx/dy/dz
    deriv[15] = temp * dterm_6_dx;
    deriv[16] = temp * dterm_6_dy;
    deriv[17] = temp * dterm_6_dz;

}

void calc_solid_MCSH_7_7(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

    double P1x = P1(lambda, x0, gamma);
    double P1y = P1(lambda, y0, gamma);
    double P1z = P1(lambda, z0, gamma);

    double P3x = P3(lambda, x0, gamma);
    double P3y = P3(lambda, y0, gamma);
    double P3z = P3(lambda, z0, gamma);

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double term_1 = (89775.0 * P3x * P3y * P1z) - (22680.0 * P5x * P1y * P1z) - (14175.0 * P3x * P1y * P3z) - (22680.0 * P1x * P5y * P1z) - (14175.0 * P1x * P3y * P3z) + (8505.0 * P1x * P1y * P5z);
    double term_2 = (89775.0 * P3x * P1y * P3z) - (22680.0 * P5x * P1y * P1z) - (14175.0 * P3x * P3y * P1z) - (22680.0 * P1x * P1y * P5z) - (14175.0 * P1x * P3y * P3z) + (8505.0 * P1x * P5y * P1z);
    double term_3 = (89775.0 * P1x * P3y * P3z) - (22680.0 * P1x * P5y * P1z) - (14175.0 * P3x * P3y * P1z) - (22680.0 * P1x * P1y * P5z) - (14175.0 * P3x * P1y * P3z) + (8505.0 * P5x * P1y * P1z);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dP3x_exp = dP3_exp(P3x, C2, lambda, x0, gamma);
    double dP3y_exp = dP3_exp(P3y, C2, lambda, y0, gamma);
    double dP3z_exp = dP3_exp(P3z, C2, lambda, z0, gamma);

    double dP5x_exp = dP5_exp(P5x, C2, lambda, x0, gamma);
    double dP5y_exp = dP5_exp(P5y, C2, lambda, y0, gamma);
    double dP5z_exp = dP5_exp(P5z, C2, lambda, z0, gamma);

    double dterm_1_dx = (89775.0 * dP3x_exp * P3y * P1z) - (22680.0 * dP5x_exp * P1y * P1z) - (14175.0 * dP3x_exp * P1y * P3z) - (22680.0 * dP1x_exp * P5y * P1z) - (14175.0 * dP1x_exp * P3y * P3z) + (8505.0 * dP1x_exp * P1y * P5z);
    double dterm_1_dy = (89775.0 * P3x * dP3y_exp * P1z) - (22680.0 * P5x * dP1y_exp * P1z) - (14175.0 * P3x * dP1y_exp * P3z) - (22680.0 * P1x * dP5y_exp * P1z) - (14175.0 * P1x * dP3y_exp * P3z) + (8505.0 * P1x * dP1y_exp * P5z);
    double dterm_1_dz = (89775.0 * P3x * P3y * dP1z_exp) - (22680.0 * P5x * P1y * dP1z_exp) - (14175.0 * P3x * P1y * dP3z_exp) - (22680.0 * P1x * P5y * dP1z_exp) - (14175.0 * P1x * P3y * dP3z_exp) + (8505.0 * P1x * P1y * dP5z_exp);

    double dterm_2_dx = (89775.0 * dP3x_exp * P1y * P3z) - (22680.0 * dP5x_exp * P1y * P1z) - (14175.0 * dP3x_exp * P3y * P1z) - (22680.0 * dP1x_exp * P1y * P5z) - (14175.0 * dP1x_exp * P3y * P3z) + (8505.0 * dP1x_exp * P5y * P1z);
    double dterm_2_dy = (89775.0 * P3x * dP1y_exp * P3z) - (22680.0 * P5x * dP1y_exp * P1z) - (14175.0 * P3x * dP3y_exp * P1z) - (22680.0 * P1x * dP1y_exp * P5z) - (14175.0 * P1x * dP3y_exp * P3z) + (8505.0 * P1x * dP5y_exp * P1z);
    double dterm_2_dz = (89775.0 * P3x * P1y * dP3z_exp) - (22680.0 * P5x * P1y * dP1z_exp) - (14175.0 * P3x * P3y * dP1z_exp) - (22680.0 * P1x * P1y * dP5z_exp) - (14175.0 * P1x * P3y * dP3z_exp) + (8505.0 * P1x * P5y * dP1z_exp);

    double dterm_3_dx = (89775.0 * dP1x_exp * P3y * P3z) - (22680.0 * dP1x_exp * P5y * P1z) - (14175.0 * dP3x_exp * P3y * P1z) - (22680.0 * dP1x_exp * P1y * P5z) - (14175.0 * dP3x_exp * P1y * P3z) + (8505.0 * dP5x_exp * P1y * P1z);
    double dterm_3_dy = (89775.0 * P1x * dP3y_exp * P3z) - (22680.0 * P1x * dP5y_exp * P1z) - (14175.0 * P3x * dP3y_exp * P1z) - (22680.0 * P1x * dP1y_exp * P5z) - (14175.0 * P3x * dP1y_exp * P3z) + (8505.0 * P5x * dP1y_exp * P1z);
    double dterm_3_dz = (89775.0 * P1x * P3y * dP3z_exp) - (22680.0 * P1x * P5y * dP1z_exp) - (14175.0 * P3x * P3y * dP1z_exp) - (22680.0 * P1x * P1y * dP5z_exp) - (14175.0 * P3x * P1y * dP3z_exp) + (8505.0 * P5x * P1y * dP1z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

}

void calc_solid_MCSH_7_8(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double term_1 = (94500.0 * P3x * P2y * P2z) - (6615.0 * P5x * P2y) - (4725.0 * P3x * P4y) - (6615.0 * P5x * P2z) - (4725.0 * P3x * P4z) + (630.0 * P7x) - (23625.0 * P1x * P4y * P2z) - (23625.0 * P1x * P2y * P4z) + (2520.0 * P1x * P6y) + (2520.0 * P1x * P6z);
    double term_2 = (94500.0 * P2x * P3y * P2z) - (6615.0 * P2x * P5y) - (4725.0 * P4x * P3y) - (6615.0 * P5y * P2z) - (4725.0 * P3y * P4z) + (630.0 * P7y) - (23625.0 * P4x * P1y * P2z) - (23625.0 * P2x * P1y * P4z) + (2520.0 * P6x * P1y) + (2520.0 * P1y * P6z);
    double term_3 = (94500.0 * P2x * P2y * P3z) - (6615.0 * P2x * P5z) - (4725.0 * P4x * P3z) - (6615.0 * P2y * P5z) - (4725.0 * P4y * P3z) + (630.0 * P7z) - (23625.0 * P4x * P2y * P1z) - (23625.0 * P2x * P4y * P1z) + (2520.0 * P6x * P1z) + (2520.0 * P6y * P1z);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP3x_exp = dP3_exp(P3x, C2, lambda, x0, gamma);
    double dP3y_exp = dP3_exp(P3y, C2, lambda, y0, gamma);
    double dP3z_exp = dP3_exp(P3z, C2, lambda, z0, gamma);

    double dP4x_exp = dP4_exp(P4x, C2, lambda, x0, gamma);
    double dP4y_exp = dP4_exp(P4y, C2, lambda, y0, gamma);
    double dP4z_exp = dP4_exp(P4z, C2, lambda, z0, gamma);

    double dP5x_exp = dP5_exp(P5x, C2, lambda, x0, gamma);
    double dP5y_exp = dP5_exp(P5y, C2, lambda, y0, gamma);
    double dP5z_exp = dP5_exp(P5z, C2, lambda, z0, gamma);

    double dP6x_exp = dP6_exp(P6x, C2, lambda, x0, gamma);
    double dP6y_exp = dP6_exp(P6y, C2, lambda, y0, gamma);
    double dP6z_exp = dP6_exp(P6z, C2, lambda, z0, gamma);

    double dP7x_exp = dP7_exp(P7x, C2, lambda, x0, gamma);
    double dP7y_exp = dP7_exp(P7y, C2, lambda, y0, gamma);
    double dP7z_exp = dP7_exp(P7z, C2, lambda, z0, gamma);

    double dterm_1_dx = (94500.0 * dP3x_exp * P2y * P2z) - (6615.0 * dP5x_exp * P2y) - (4725.0 * dP3x_exp * P4y) - (6615.0 * dP5x_exp * P2z) - (4725.0 * dP3x_exp * P4z) + (630.0 * dP7x_exp) - (23625.0 * dP1x_exp * P4y * P2z) - (23625.0 * dP1x_exp * P2y * P4z) + (2520.0 * dP1x_exp * P6y) + (2520.0 * dP1x_exp * P6z);
    double dterm_1_dy = (94500.0 * P3x * dP2y_exp * P2z) - (6615.0 * P5x * dP2y_exp) - (4725.0 * P3x * dP4y_exp) - (6615.0 * P5x * dP0y_exp * P2z) - (4725.0 * P3x * dP0y_exp * P4z) + (630.0 * P7x * dP0y_exp) - (23625.0 * P1x * dP4y_exp * P2z) - (23625.0 * P1x * dP2y_exp * P4z) + (2520.0 * P1x * dP6y_exp) + (2520.0 * P1x * dP0y_exp * P6z);
    double dterm_1_dz = (94500.0 * P3x * P2y * dP2z_exp) - (6615.0 * P5x * P2y * dP0z_exp) - (4725.0 * P3x * P4y * dP0z_exp) - (6615.0 * P5x * dP2z_exp) - (4725.0 * P3x * dP4z_exp) + (630.0 * P7x * dP0z_exp) - (23625.0 * P1x * P4y * dP2z_exp) - (23625.0 * P1x * P2y * dP4z_exp) + (2520.0 * P1x * P6y * dP0z_exp) + (2520.0 * P1x * dP6z_exp);

    double dterm_2_dx = (94500.0 * dP2x_exp * P3y * P2z) - (6615.0 * dP2x_exp * P5y) - (4725.0 * dP4x_exp * P3y) - (6615.0 * dP0x_exp * P5y * P2z) - (4725.0 * dP0x_exp * P3y * P4z) + (630.0 * dP0x_exp * P7y) - (23625.0 * dP4x_exp * P1y * P2z) - (23625.0 * dP2x_exp * P1y * P4z) + (2520.0 * dP6x_exp * P1y) + (2520.0 * dP0x_exp * P1y * P6z);
    double dterm_2_dy = (94500.0 * P2x * dP3y_exp * P2z) - (6615.0 * P2x * dP5y_exp) - (4725.0 * P4x * dP3y_exp) - (6615.0 * dP5y_exp * P2z) - (4725.0 * dP3y_exp * P4z) + (630.0 * dP7y_exp) - (23625.0 * P4x * dP1y_exp * P2z) - (23625.0 * P2x * dP1y_exp * P4z) + (2520.0 * P6x * dP1y_exp) + (2520.0 * dP1y_exp * P6z);
    double dterm_2_dz = (94500.0 * P2x * P3y * dP2z_exp) - (6615.0 * P2x * P5y * dP0z_exp) - (4725.0 * P4x * P3y * dP0z_exp) - (6615.0 * P5y * dP2z_exp) - (4725.0 * P3y * dP4z_exp) + (630.0 * P7y * dP0z_exp) - (23625.0 * P4x * P1y * dP2z_exp) - (23625.0 * P2x * P1y * dP4z_exp) + (2520.0 * P6x * P1y * dP0z_exp) + (2520.0 * P1y * dP6z_exp);

    double dterm_3_dx = (94500.0 * dP2x_exp * P2y * P3z) - (6615.0 * dP2x_exp * P5z) - (4725.0 * dP4x_exp * P3z) - (6615.0 * dP0x_exp * P2y * P5z) - (4725.0 * dP0x_exp * P4y * P3z) + (630.0 * dP0x_exp * P7z) - (23625.0 * dP4x_exp * P2y * P1z) - (23625.0 * dP2x_exp * P4y * P1z) + (2520.0 * dP6x_exp * P1z) + (2520.0 * dP0x_exp * P6y * P1z);
    double dterm_3_dy = (94500.0 * P2x * dP2y_exp * P3z) - (6615.0 * P2x * dP0y_exp * P5z) - (4725.0 * P4x * dP0y_exp * P3z) - (6615.0 * dP2y_exp * P5z) - (4725.0 * dP4y_exp * P3z) + (630.0 * dP0y_exp * P7z) - (23625.0 * P4x * dP2y_exp * P1z) - (23625.0 * P2x * dP4y_exp * P1z) + (2520.0 * P6x * dP0y_exp * P1z) + (2520.0 * dP6y_exp * P1z);
    double dterm_3_dz = (94500.0 * P2x * P2y * dP3z_exp) - (6615.0 * P2x * dP5z_exp) - (4725.0 * P4x * dP3z_exp) - (6615.0 * P2y * dP5z_exp) - (4725.0 * P4y * dP3z_exp) + (630.0 * dP7z_exp) - (23625.0 * P4x * P2y * dP1z_exp) - (23625.0 * P2x * P4y * dP1z_exp) + (2520.0 * P6x * dP1z_exp) + (2520.0 * P6y * dP1z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

}


void calc_solid_MCSH_8_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double term_1 = (40320.0 * P8x) - (564480.0 * P6x * P2y) - (564480.0 * P6x * P2z) + (1058400.0 * P4x * P4y) + (2116800.0 * P4x * P2y * P2z) + (1058400.0 * P4x * P4z) - (352800.0 * P2x * P6y) - (1058400.0 * P2x * P4y * P2z) - (1058400.0 * P2x * P2y * P4z) - (352800.0 * P2x * P6z) + (11025.0 * P8y) + (44100.0 * P6y * P2z) + (66150.0 * P4y * P4z) + (44100.0 * P2y * P6z) + (11025.0 * P8z);
    double term_2 = (40320.0 * P8y) - (564480.0 * P2x * P6y) - (564480.0 * P6y * P2z) + (1058400.0 * P4x * P4y) + (2116800.0 * P2x * P4y * P2z) + (1058400.0 * P4y * P4z) - (352800.0 * P6x * P2y) - (1058400.0 * P4x * P2y * P2z) - (1058400.0 * P2x * P2y * P4z) - (352800.0 * P2y * P6z) + (11025.0 * P8x) + (44100.0 * P6x * P2z) + (66150.0 * P4x * P4z) + (44100.0 * P2x * P6z) + (11025.0 * P8z);
    double term_3 = (40320.0 * P8z) - (564480.0 * P2x * P6z) - (564480.0 * P2y * P6z) + (1058400.0 * P4x * P4z) + (2116800.0 * P2x * P2y * P4z) + (1058400.0 * P4y * P4z) - (352800.0 * P6x * P2z) - (1058400.0 * P4x * P2y * P2z) - (1058400.0 * P2x * P4y * P2z) - (352800.0 * P6y * P2z) + (11025.0 * P8x) + (44100.0 * P6x * P2y) + (66150.0 * P4x * P4y) + (44100.0 * P2x * P6y) + (11025.0 * P8y);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP4x_exp = dP4_exp(P4x, C2, lambda, x0, gamma);
    double dP4y_exp = dP4_exp(P4y, C2, lambda, y0, gamma);
    double dP4z_exp = dP4_exp(P4z, C2, lambda, z0, gamma);

    double dP6x_exp = dP6_exp(P6x, C2, lambda, x0, gamma);
    double dP6y_exp = dP6_exp(P6y, C2, lambda, y0, gamma);
    double dP6z_exp = dP6_exp(P6z, C2, lambda, z0, gamma);

    double dP8x_exp = dP8_exp(P8x, C2, lambda, x0, gamma);
    double dP8y_exp = dP8_exp(P8y, C2, lambda, y0, gamma);
    double dP8z_exp = dP8_exp(P8z, C2, lambda, z0, gamma);

    double dterm_1_dx = (40320.0 * dP8x_exp) - (564480.0 * dP6x_exp * P2y) - (564480.0 * dP6x_exp * P2z) + (1058400.0 * dP4x_exp * P4y) + (2116800.0 * dP4x_exp * P2y * P2z) + (1058400.0 * dP4x_exp * P4z) - (352800.0 * dP2x_exp * P6y) - (1058400.0 * dP2x_exp * P4y * P2z) - (1058400.0 * dP2x_exp * P2y * P4z) - (352800.0 * dP2x_exp * P6z) + (11025.0 * dP0x_exp * P8y) + (44100.0 * dP0x_exp * P6y * P2z) + (66150.0 * dP0x_exp * P4y * P4z) + (44100.0 * dP0x_exp * P2y * P6z) + (11025.0 * dP0x_exp * P8z);
    double dterm_1_dy = (40320.0 * P8x * dP0y_exp) - (564480.0 * P6x * dP2y_exp) - (564480.0 * P6x * dP0y_exp * P2z) + (1058400.0 * P4x * dP4y_exp) + (2116800.0 * P4x * dP2y_exp * P2z) + (1058400.0 * P4x * dP0y_exp * P4z) - (352800.0 * P2x * dP6y_exp) - (1058400.0 * P2x * dP4y_exp * P2z) - (1058400.0 * P2x * dP2y_exp * P4z) - (352800.0 * P2x * dP0y_exp * P6z) + (11025.0 * dP8y_exp) + (44100.0 * dP6y_exp * P2z) + (66150.0 * dP4y_exp * P4z) + (44100.0 * dP2y_exp * P6z) + (11025.0 * dP0y_exp * P8z);
    double dterm_1_dz = (40320.0 * P8x * dP0z_exp) - (564480.0 * P6x * P2y * dP0z_exp) - (564480.0 * P6x * dP2z_exp) + (1058400.0 * P4x * P4y * dP0z_exp) + (2116800.0 * P4x * P2y * dP2z_exp) + (1058400.0 * P4x * dP4z_exp) - (352800.0 * P2x * P6y * dP0z_exp) - (1058400.0 * P2x * P4y * dP2z_exp) - (1058400.0 * P2x * P2y * dP4z_exp) - (352800.0 * P2x * dP6z_exp) + (11025.0 * P8y * dP0z_exp) + (44100.0 * P6y * dP2z_exp) + (66150.0 * P4y * dP4z_exp) + (44100.0 * P2y * dP6z_exp) + (11025.0 * dP8z_exp);

    double dterm_2_dx = (40320.0 * dP0x_exp * P8y) - (564480.0 * dP2x_exp * P6y) - (564480.0 * dP0x_exp * P6y * P2z) + (1058400.0 * dP4x_exp * P4y) + (2116800.0 * dP2x_exp * P4y * P2z) + (1058400.0 * dP0x_exp * P4y * P4z) - (352800.0 * dP6x_exp * P2y) - (1058400.0 * dP4x_exp * P2y * P2z) - (1058400.0 * dP2x_exp * P2y * P4z) - (352800.0 * dP0x_exp * P2y * P6z) + (11025.0 * dP8x_exp) + (44100.0 * dP6x_exp * P2z) + (66150.0 * dP4x_exp * P4z) + (44100.0 * dP2x_exp * P6z) + (11025.0 * dP0x_exp * P8z);
    double dterm_2_dy = (40320.0 * dP8y_exp) - (564480.0 * P2x * dP6y_exp) - (564480.0 * dP6y_exp * P2z) + (1058400.0 * P4x * dP4y_exp) + (2116800.0 * P2x * dP4y_exp * P2z) + (1058400.0 * dP4y_exp * P4z) - (352800.0 * P6x * dP2y_exp) - (1058400.0 * P4x * dP2y_exp * P2z) - (1058400.0 * P2x * dP2y_exp * P4z) - (352800.0 * dP2y_exp * P6z) + (11025.0 * P8x * dP0y_exp) + (44100.0 * P6x * dP0y_exp * P2z) + (66150.0 * P4x * dP0y_exp * P4z) + (44100.0 * P2x * dP0y_exp * P6z) + (11025.0 * dP0y_exp * P8z);
    double dterm_2_dz = (40320.0 * P8y * dP0z_exp) - (564480.0 * P2x * P6y * dP0z_exp) - (564480.0 * P6y * dP2z_exp) + (1058400.0 * P4x * P4y * dP0z_exp) + (2116800.0 * P2x * P4y * dP2z_exp) + (1058400.0 * P4y * dP4z_exp) - (352800.0 * P6x * P2y * dP0z_exp) - (1058400.0 * P4x * P2y * dP2z_exp) - (1058400.0 * P2x * P2y * dP4z_exp) - (352800.0 * P2y * dP6z_exp) + (11025.0 * P8x * dP0z_exp) + (44100.0 * P6x * dP2z_exp) + (66150.0 * P4x * dP4z_exp) + (44100.0 * P2x * dP6z_exp) + (11025.0 * dP8z_exp);

    double dterm_3_dx = (40320.0 * dP0x_exp * P8z) - (564480.0 * dP2x_exp * P6z) - (564480.0 * dP0x_exp * P2y * P6z) + (1058400.0 * dP4x_exp * P4z) + (2116800.0 * dP2x_exp * P2y * P4z) + (1058400.0 * dP0x_exp * P4y * P4z) - (352800.0 * dP6x_exp * P2z) - (1058400.0 * dP4x_exp * P2y * P2z) - (1058400.0 * dP2x_exp * P4y * P2z) - (352800.0 * dP0x_exp * P6y * P2z) + (11025.0 * dP8x_exp) + (44100.0 * dP6x_exp * P2y) + (66150.0 * dP4x_exp * P4y) + (44100.0 * dP2x_exp * P6y) + (11025.0 * dP0x_exp * P8y);
    double dterm_3_dy = (40320.0 * dP0y_exp * P8z) - (564480.0 * P2x * dP0y_exp * P6z) - (564480.0 * dP2y_exp * P6z) + (1058400.0 * P4x * dP0y_exp * P4z) + (2116800.0 * P2x * dP2y_exp * P4z) + (1058400.0 * dP4y_exp * P4z) - (352800.0 * P6x * dP0y_exp * P2z) - (1058400.0 * P4x * dP2y_exp * P2z) - (1058400.0 * P2x * dP4y_exp * P2z) - (352800.0 * dP6y_exp * P2z) + (11025.0 * P8x * dP0y_exp) + (44100.0 * P6x * dP2y_exp) + (66150.0 * P4x * dP4y_exp) + (44100.0 * P2x * dP6y_exp) + (11025.0 * dP8y_exp);
    double dterm_3_dz = (40320.0 * dP8z_exp) - (564480.0 * P2x * dP6z_exp) - (564480.0 * P2y * dP6z_exp) + (1058400.0 * P4x * dP4z_exp) + (2116800.0 * P2x * P2y * dP4z_exp) + (1058400.0 * P4y * dP4z_exp) - (352800.0 * P6x * dP2z_exp) - (1058400.0 * P4x * P2y * dP2z_exp) - (1058400.0 * P2x * P4y * dP2z_exp) - (352800.0 * P6y * dP2z_exp) + (11025.0 * P8x * dP0z_exp) + (44100.0 * P6x * P2y * dP0z_exp) + (66150.0 * P4x * P4y * dP0z_exp) + (44100.0 * P2x * P6y * dP0z_exp) + (11025.0 * P8y * dP0z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

}

void calc_solid_MCSH_8_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double term_1 = (181440.0 * P7x * P1y) - (952560.0 * P5x * P3y) - (952560.0 * P5x * P1y * P2z) + (793800.0 * P3x * P5y) + (1587600.0 * P3x * P3y * P2z) + (793800.0 * P3x * P1y * P4z) - (99225.0 * P1x * P7y) - (297675.0 * P1x * P5y * P2z) - (297675.0 * P1x * P3y * P4z) - (99225.0 * P1x * P1y * P6z);
    double term_2 = (181440.0 * P1x * P7y) - (952560.0 * P3x * P5y) - (952560.0 * P1x * P5y * P2z) + (793800.0 * P5x * P3y) + (1587600.0 * P3x * P3y * P2z) + (793800.0 * P1x * P3y * P4z) - (99225.0 * P7x * P1y) - (297675.0 * P5x * P1y * P2z) - (297675.0 * P3x * P1y * P4z) - (99225.0 * P1x * P1y * P6z);
    double term_3 = (181440.0 * P7x * P1z) - (952560.0 * P5x * P3z) - (952560.0 * P5x * P2y * P1z) + (793800.0 * P3x * P5z) + (1587600.0 * P3x * P2y * P3z) + (793800.0 * P3x * P4y * P1z) - (99225.0 * P1x * P7z) - (297675.0 * P1x * P2y * P5z) - (297675.0 * P1x * P4y * P3z) - (99225.0 * P1x * P6y * P1z);
    double term_4 = (181440.0 * P1x * P7z) - (952560.0 * P3x * P5z) - (952560.0 * P1x * P2y * P5z) + (793800.0 * P5x * P3z) + (1587600.0 * P3x * P2y * P3z) + (793800.0 * P1x * P4y * P3z) - (99225.0 * P7x * P1z) - (297675.0 * P5x * P2y * P1z) - (297675.0 * P3x * P4y * P1z) - (99225.0 * P1x * P6y * P1z);
    double term_5 = (181440.0 * P7y * P1z) - (952560.0 * P5y * P3z) - (952560.0 * P2x * P5y * P1z) + (793800.0 * P3y * P5z) + (1587600.0 * P2x * P3y * P3z) + (793800.0 * P4x * P3y * P1z) - (99225.0 * P1y * P7z) - (297675.0 * P2x * P1y * P5z) - (297675.0 * P4x * P1y * P3z) - (99225.0 * P6x * P1y * P1z);
    double term_6 = (181440.0 * P1y * P7z) - (952560.0 * P3y * P5z) - (952560.0 * P2x * P1y * P5z) + (793800.0 * P5y * P3z) + (1587600.0 * P2x * P3y * P3z) + (793800.0 * P4x * P1y * P3z) - (99225.0 * P7y * P1z) - (297675.0 * P2x * P5y * P1z) - (297675.0 * P4x * P3y * P1z) - (99225.0 * P6x * P1y * P1z);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP3x_exp = dP3_exp(P3x, C2, lambda, x0, gamma);
    double dP3y_exp = dP3_exp(P3y, C2, lambda, y0, gamma);
    double dP3z_exp = dP3_exp(P3z, C2, lambda, z0, gamma);

    double dP4x_exp = dP4_exp(P4x, C2, lambda, x0, gamma);
    double dP4y_exp = dP4_exp(P4y, C2, lambda, y0, gamma);
    double dP4z_exp = dP4_exp(P4z, C2, lambda, z0, gamma);

    double dP5x_exp = dP5_exp(P5x, C2, lambda, x0, gamma);
    double dP5y_exp = dP5_exp(P5y, C2, lambda, y0, gamma);
    double dP5z_exp = dP5_exp(P5z, C2, lambda, z0, gamma);

    double dP6x_exp = dP6_exp(P6x, C2, lambda, x0, gamma);
    double dP6y_exp = dP6_exp(P6y, C2, lambda, y0, gamma);
    double dP6z_exp = dP6_exp(P6z, C2, lambda, z0, gamma);

    double dP7x_exp = dP7_exp(P7x, C2, lambda, x0, gamma);
    double dP7y_exp = dP7_exp(P7y, C2, lambda, y0, gamma);
    double dP7z_exp = dP7_exp(P7z, C2, lambda, z0, gamma);

    double dterm_1_dx = (181440.0 * dP7x_exp * P1y) - (952560.0 * dP5x_exp * P3y) - (952560.0 * dP5x_exp * P1y * P2z) + (793800.0 * dP3x_exp * P5y) + (1587600.0 * dP3x_exp * P3y * P2z) + (793800.0 * dP3x_exp * P1y * P4z) - (99225.0 * dP1x_exp * P7y) - (297675.0 * dP1x_exp * P5y * P2z) - (297675.0 * dP1x_exp * P3y * P4z) - (99225.0 * dP1x_exp * P1y * P6z);
    double dterm_1_dy = (181440.0 * P7x * dP1y_exp) - (952560.0 * P5x * dP3y_exp) - (952560.0 * P5x * dP1y_exp * P2z) + (793800.0 * P3x * dP5y_exp) + (1587600.0 * P3x * dP3y_exp * P2z) + (793800.0 * P3x * dP1y_exp * P4z) - (99225.0 * P1x * dP7y_exp) - (297675.0 * P1x * dP5y_exp * P2z) - (297675.0 * P1x * dP3y_exp * P4z) - (99225.0 * P1x * dP1y_exp * P6z);
    double dterm_1_dz = (181440.0 * P7x * P1y * dP0z_exp) - (952560.0 * P5x * P3y * dP0z_exp) - (952560.0 * P5x * P1y * dP2z_exp) + (793800.0 * P3x * P5y * dP0z_exp) + (1587600.0 * P3x * P3y * dP2z_exp) + (793800.0 * P3x * P1y * dP4z_exp) - (99225.0 * P1x * P7y * dP0z_exp) - (297675.0 * P1x * P5y * dP2z_exp) - (297675.0 * P1x * P3y * dP4z_exp) - (99225.0 * P1x * P1y * dP6z_exp);

    double dterm_2_dx = (181440.0 * dP1x_exp * P7y) - (952560.0 * dP3x_exp * P5y) - (952560.0 * dP1x_exp * P5y * P2z) + (793800.0 * dP5x_exp * P3y) + (1587600.0 * dP3x_exp * P3y * P2z) + (793800.0 * dP1x_exp * P3y * P4z) - (99225.0 * dP7x_exp * P1y) - (297675.0 * dP5x_exp * P1y * P2z) - (297675.0 * dP3x_exp * P1y * P4z) - (99225.0 * dP1x_exp * P1y * P6z);
    double dterm_2_dy = (181440.0 * P1x * dP7y_exp) - (952560.0 * P3x * dP5y_exp) - (952560.0 * P1x * dP5y_exp * P2z) + (793800.0 * P5x * dP3y_exp) + (1587600.0 * P3x * dP3y_exp * P2z) + (793800.0 * P1x * dP3y_exp * P4z) - (99225.0 * P7x * dP1y_exp) - (297675.0 * P5x * dP1y_exp * P2z) - (297675.0 * P3x * dP1y_exp * P4z) - (99225.0 * P1x * dP1y_exp * P6z);
    double dterm_2_dz = (181440.0 * P1x * P7y * dP0z_exp) - (952560.0 * P3x * P5y * dP0z_exp) - (952560.0 * P1x * P5y * dP2z_exp) + (793800.0 * P5x * P3y * dP0z_exp) + (1587600.0 * P3x * P3y * dP2z_exp) + (793800.0 * P1x * P3y * dP4z_exp) - (99225.0 * P7x * P1y * dP0z_exp) - (297675.0 * P5x * P1y * dP2z_exp) - (297675.0 * P3x * P1y * dP4z_exp) - (99225.0 * P1x * P1y * dP6z_exp);

    double dterm_3_dx = (181440.0 * dP7x_exp * P1z) - (952560.0 * dP5x_exp * P3z) - (952560.0 * dP5x_exp * P2y * P1z) + (793800.0 * dP3x_exp * P5z) + (1587600.0 * dP3x_exp * P2y * P3z) + (793800.0 * dP3x_exp * P4y * P1z) - (99225.0 * dP1x_exp * P7z) - (297675.0 * dP1x_exp * P2y * P5z) - (297675.0 * dP1x_exp * P4y * P3z) - (99225.0 * dP1x_exp * P6y * P1z);
    double dterm_3_dy = (181440.0 * P7x * dP0y_exp * P1z) - (952560.0 * P5x * dP0y_exp * P3z) - (952560.0 * P5x * dP2y_exp * P1z) + (793800.0 * P3x * dP0y_exp * P5z) + (1587600.0 * P3x * dP2y_exp * P3z) + (793800.0 * P3x * dP4y_exp * P1z) - (99225.0 * P1x * dP0y_exp * P7z) - (297675.0 * P1x * dP2y_exp * P5z) - (297675.0 * P1x * dP4y_exp * P3z) - (99225.0 * P1x * dP6y_exp * P1z);
    double dterm_3_dz = (181440.0 * P7x * dP1z_exp) - (952560.0 * P5x * dP3z_exp) - (952560.0 * P5x * P2y * dP1z_exp) + (793800.0 * P3x * dP5z_exp) + (1587600.0 * P3x * P2y * dP3z_exp) + (793800.0 * P3x * P4y * dP1z_exp) - (99225.0 * P1x * dP7z_exp) - (297675.0 * P1x * P2y * dP5z_exp) - (297675.0 * P1x * P4y * dP3z_exp) - (99225.0 * P1x * P6y * dP1z_exp);

    double dterm_4_dx = (181440.0 * dP1x_exp * P7z) - (952560.0 * dP3x_exp * P5z) - (952560.0 * dP1x_exp * P2y * P5z) + (793800.0 * dP5x_exp * P3z) + (1587600.0 * dP3x_exp * P2y * P3z) + (793800.0 * dP1x_exp * P4y * P3z) - (99225.0 * dP7x_exp * P1z) - (297675.0 * dP5x_exp * P2y * P1z) - (297675.0 * dP3x_exp * P4y * P1z) - (99225.0 * dP1x_exp * P6y * P1z);
    double dterm_4_dy = (181440.0 * P1x * dP0y_exp * P7z) - (952560.0 * P3x * dP0y_exp * P5z) - (952560.0 * P1x * dP2y_exp * P5z) + (793800.0 * P5x * dP0y_exp * P3z) + (1587600.0 * P3x * dP2y_exp * P3z) + (793800.0 * P1x * dP4y_exp * P3z) - (99225.0 * P7x * dP0y_exp * P1z) - (297675.0 * P5x * dP2y_exp * P1z) - (297675.0 * P3x * dP4y_exp * P1z) - (99225.0 * P1x * dP6y_exp * P1z);
    double dterm_4_dz = (181440.0 * P1x * dP7z_exp) - (952560.0 * P3x * dP5z_exp) - (952560.0 * P1x * P2y * dP5z_exp) + (793800.0 * P5x * dP3z_exp) + (1587600.0 * P3x * P2y * dP3z_exp) + (793800.0 * P1x * P4y * dP3z_exp) - (99225.0 * P7x * dP1z_exp) - (297675.0 * P5x * P2y * dP1z_exp) - (297675.0 * P3x * P4y * dP1z_exp) - (99225.0 * P1x * P6y * dP1z_exp);

    double dterm_5_dx = (181440.0 * dP0x_exp * P7y * P1z) - (952560.0 * dP0x_exp * P5y * P3z) - (952560.0 * dP2x_exp * P5y * P1z) + (793800.0 * dP0x_exp * P3y * P5z) + (1587600.0 * dP2x_exp * P3y * P3z) + (793800.0 * dP4x_exp * P3y * P1z) - (99225.0 * dP0x_exp * P1y * P7z) - (297675.0 * dP2x_exp * P1y * P5z) - (297675.0 * dP4x_exp * P1y * P3z) - (99225.0 * dP6x_exp * P1y * P1z);
    double dterm_5_dy = (181440.0 * dP7y_exp * P1z) - (952560.0 * dP5y_exp * P3z) - (952560.0 * P2x * dP5y_exp * P1z) + (793800.0 * dP3y_exp * P5z) + (1587600.0 * P2x * dP3y_exp * P3z) + (793800.0 * P4x * dP3y_exp * P1z) - (99225.0 * dP1y_exp * P7z) - (297675.0 * P2x * dP1y_exp * P5z) - (297675.0 * P4x * dP1y_exp * P3z) - (99225.0 * P6x * dP1y_exp * P1z);
    double dterm_5_dz = (181440.0 * P7y * dP1z_exp) - (952560.0 * P5y * dP3z_exp) - (952560.0 * P2x * P5y * dP1z_exp) + (793800.0 * P3y * dP5z_exp) + (1587600.0 * P2x * P3y * dP3z_exp) + (793800.0 * P4x * P3y * dP1z_exp) - (99225.0 * P1y * dP7z_exp) - (297675.0 * P2x * P1y * dP5z_exp) - (297675.0 * P4x * P1y * dP3z_exp) - (99225.0 * P6x * P1y * dP1z_exp);

    double dterm_6_dx = (181440.0 * dP0x_exp * P1y * P7z) - (952560.0 * dP0x_exp * P3y * P5z) - (952560.0 * dP2x_exp * P1y * P5z) + (793800.0 * dP0x_exp * P5y * P3z) + (1587600.0 * dP2x_exp * P3y * P3z) + (793800.0 * dP4x_exp * P1y * P3z) - (99225.0 * dP0x_exp * P7y * P1z) - (297675.0 * dP2x_exp * P5y * P1z) - (297675.0 * dP4x_exp * P3y * P1z) - (99225.0 * dP6x_exp * P1y * P1z);
    double dterm_6_dy = (181440.0 * dP1y_exp * P7z) - (952560.0 * dP3y_exp * P5z) - (952560.0 * P2x * dP1y_exp * P5z) + (793800.0 * dP5y_exp * P3z) + (1587600.0 * P2x * dP3y_exp * P3z) + (793800.0 * P4x * dP1y_exp * P3z) - (99225.0 * dP7y_exp * P1z) - (297675.0 * P2x * dP5y_exp * P1z) - (297675.0 * P4x * dP3y_exp * P1z) - (99225.0 * P6x * dP1y_exp * P1z);
    double dterm_6_dz = (181440.0 * P1y * dP7z_exp) - (952560.0 * P3y * dP5z_exp) - (952560.0 * P2x * P1y * dP5z_exp) + (793800.0 * P5y * dP3z_exp) + (1587600.0 * P2x * P3y * dP3z_exp) + (793800.0 * P4x * P1y * dP3z_exp) - (99225.0 * P7y * dP1z_exp) - (297675.0 * P2x * P5y * dP1z_exp) - (297675.0 * P4x * P3y * dP1z_exp) - (99225.0 * P6x * P1y * dP1z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

    // dmiu4 dx/dy/dz
    deriv[9] = temp * dterm_4_dx;
    deriv[10] = temp * dterm_4_dy;
    deriv[11] = temp * dterm_4_dz;

    // dmiu5 dx/dy/dz
    deriv[12] = temp * dterm_5_dx;
    deriv[13] = temp * dterm_5_dy;
    deriv[14] = temp * dterm_5_dz;

    // dmiu6 dx/dy/dz
    deriv[15] = temp * dterm_6_dx;
    deriv[16] = temp * dterm_6_dy;
    deriv[17] = temp * dterm_6_dz;

}

void calc_solid_MCSH_8_3(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double term_1 = (509040.0 * P6x * P2y) - (20160.0 * P8x) + (55440.0 * P6x * P2z) - (1096200.0 * P4x * P4y) - (1058400.0 * P4x * P2y * P2z) + (37800.0 * P4x * P4z) + (389025.0 * P2x * P6y) + (741825.0 * P2x * P4y * P2z) + (316575.0 * P2x * P2y * P4z) - (36225.0 * P2x * P6z) - (12600.0 * P8y) - (36225.0 * P6y * P2z) - (33075.0 * P4y * P4z) - (7875.0 * P2y * P6z) + (1575.0 * P8z);
    double term_2 = (509040.0 * P2x * P6y) - (20160.0 * P8y) + (55440.0 * P6y * P2z) - (1096200.0 * P4x * P4y) - (1058400.0 * P2x * P4y * P2z) + (37800.0 * P4y * P4z) + (389025.0 * P6x * P2y) + (741825.0 * P4x * P2y * P2z) + (316575.0 * P2x * P2y * P4z) - (36225.0 * P2y * P6z) - (12600.0 * P8x) - (36225.0 * P6x * P2z) - (33075.0 * P4x * P4z) - (7875.0 * P2x * P6z) + (1575.0 * P8z);
    double term_3 = (509040.0 * P6x * P2z) - (20160.0 * P8x) + (55440.0 * P6x * P2y) - (1096200.0 * P4x * P4z) - (1058400.0 * P4x * P2y * P2z) + (37800.0 * P4x * P4y) + (389025.0 * P2x * P6z) + (741825.0 * P2x * P2y * P4z) + (316575.0 * P2x * P4y * P2z) - (36225.0 * P2x * P6y) - (12600.0 * P8z) - (36225.0 * P2y * P6z) - (33075.0 * P4y * P4z) - (7875.0 * P6y * P2z) + (1575.0 * P8y);
    double term_4 = (509040.0 * P2x * P6z) - (20160.0 * P8z) + (55440.0 * P2y * P6z) - (1096200.0 * P4x * P4z) - (1058400.0 * P2x * P2y * P4z) + (37800.0 * P4y * P4z) + (389025.0 * P6x * P2z) + (741825.0 * P4x * P2y * P2z) + (316575.0 * P2x * P4y * P2z) - (36225.0 * P6y * P2z) - (12600.0 * P8x) - (36225.0 * P6x * P2y) - (33075.0 * P4x * P4y) - (7875.0 * P2x * P6y) + (1575.0 * P8y);
    double term_5 = (509040.0 * P6y * P2z) - (20160.0 * P8y) + (55440.0 * P2x * P6y) - (1096200.0 * P4y * P4z) - (1058400.0 * P2x * P4y * P2z) + (37800.0 * P4x * P4y) + (389025.0 * P2y * P6z) + (741825.0 * P2x * P2y * P4z) + (316575.0 * P4x * P2y * P2z) - (36225.0 * P6x * P2y) - (12600.0 * P8z) - (36225.0 * P2x * P6z) - (33075.0 * P4x * P4z) - (7875.0 * P6x * P2z) + (1575.0 * P8x);
    double term_6 = (509040.0 * P2y * P6z) - (20160.0 * P8z) + (55440.0 * P2x * P6z) - (1096200.0 * P4y * P4z) - (1058400.0 * P2x * P2y * P4z) + (37800.0 * P4x * P4z) + (389025.0 * P6y * P2z) + (741825.0 * P2x * P4y * P2z) + (316575.0 * P4x * P2y * P2z) - (36225.0 * P6x * P2z) - (12600.0 * P8y) - (36225.0 * P2x * P6y) - (33075.0 * P4x * P4y) - (7875.0 * P6x * P2y) + (1575.0 * P8x);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP4x_exp = dP4_exp(P4x, C2, lambda, x0, gamma);
    double dP4y_exp = dP4_exp(P4y, C2, lambda, y0, gamma);
    double dP4z_exp = dP4_exp(P4z, C2, lambda, z0, gamma);

    double dP6x_exp = dP6_exp(P6x, C2, lambda, x0, gamma);
    double dP6y_exp = dP6_exp(P6y, C2, lambda, y0, gamma);
    double dP6z_exp = dP6_exp(P6z, C2, lambda, z0, gamma);

    double dP8x_exp = dP8_exp(P8x, C2, lambda, x0, gamma);
    double dP8y_exp = dP8_exp(P8y, C2, lambda, y0, gamma);
    double dP8z_exp = dP8_exp(P8z, C2, lambda, z0, gamma);

    double dterm_1_dx = (509040.0 * dP6x_exp * P2y) - (20160.0 * dP8x_exp) + (55440.0 * dP6x_exp * P2z) - (1096200.0 * dP4x_exp * P4y) - (1058400.0 * dP4x_exp * P2y * P2z) + (37800.0 * dP4x_exp * P4z) + (389025.0 * dP2x_exp * P6y) + (741825.0 * dP2x_exp * P4y * P2z) + (316575.0 * dP2x_exp * P2y * P4z) - (36225.0 * dP2x_exp * P6z) - (12600.0 * dP0x_exp * P8y) - (36225.0 * dP0x_exp * P6y * P2z) - (33075.0 * dP0x_exp * P4y * P4z) - (7875.0 * dP0x_exp * P2y * P6z) + (1575.0 * dP0x_exp * P8z);
    double dterm_1_dy = (509040.0 * P6x * dP2y_exp) - (20160.0 * P8x * dP0y_exp) + (55440.0 * P6x * dP0y_exp * P2z) - (1096200.0 * P4x * dP4y_exp) - (1058400.0 * P4x * dP2y_exp * P2z) + (37800.0 * P4x * dP0y_exp * P4z) + (389025.0 * P2x * dP6y_exp) + (741825.0 * P2x * dP4y_exp * P2z) + (316575.0 * P2x * dP2y_exp * P4z) - (36225.0 * P2x * dP0y_exp * P6z) - (12600.0 * dP8y_exp) - (36225.0 * dP6y_exp * P2z) - (33075.0 * dP4y_exp * P4z) - (7875.0 * dP2y_exp * P6z) + (1575.0 * dP0y_exp * P8z);
    double dterm_1_dz = (509040.0 * P6x * P2y * dP0z_exp) - (20160.0 * P8x * dP0z_exp) + (55440.0 * P6x * dP2z_exp) - (1096200.0 * P4x * P4y * dP0z_exp) - (1058400.0 * P4x * P2y * dP2z_exp) + (37800.0 * P4x * dP4z_exp) + (389025.0 * P2x * P6y * dP0z_exp) + (741825.0 * P2x * P4y * dP2z_exp) + (316575.0 * P2x * P2y * dP4z_exp) - (36225.0 * P2x * dP6z_exp) - (12600.0 * P8y * dP0z_exp) - (36225.0 * P6y * dP2z_exp) - (33075.0 * P4y * dP4z_exp) - (7875.0 * P2y * dP6z_exp) + (1575.0 * dP8z_exp);

    double dterm_2_dx = (509040.0 * dP2x_exp * P6y) - (20160.0 * dP0x_exp * P8y) + (55440.0 * dP0x_exp * P6y * P2z) - (1096200.0 * dP4x_exp * P4y) - (1058400.0 * dP2x_exp * P4y * P2z) + (37800.0 * dP0x_exp * P4y * P4z) + (389025.0 * dP6x_exp * P2y) + (741825.0 * dP4x_exp * P2y * P2z) + (316575.0 * dP2x_exp * P2y * P4z) - (36225.0 * dP0x_exp * P2y * P6z) - (12600.0 * dP8x_exp) - (36225.0 * dP6x_exp * P2z) - (33075.0 * dP4x_exp * P4z) - (7875.0 * dP2x_exp * P6z) + (1575.0 * dP0x_exp * P8z);
    double dterm_2_dy = (509040.0 * P2x * dP6y_exp) - (20160.0 * dP8y_exp) + (55440.0 * dP6y_exp * P2z) - (1096200.0 * P4x * dP4y_exp) - (1058400.0 * P2x * dP4y_exp * P2z) + (37800.0 * dP4y_exp * P4z) + (389025.0 * P6x * dP2y_exp) + (741825.0 * P4x * dP2y_exp * P2z) + (316575.0 * P2x * dP2y_exp * P4z) - (36225.0 * dP2y_exp * P6z) - (12600.0 * P8x * dP0y_exp) - (36225.0 * P6x * dP0y_exp * P2z) - (33075.0 * P4x * dP0y_exp * P4z) - (7875.0 * P2x * dP0y_exp * P6z) + (1575.0 * dP0y_exp * P8z);
    double dterm_2_dz = (509040.0 * P2x * P6y * dP0z_exp) - (20160.0 * P8y * dP0z_exp) + (55440.0 * P6y * dP2z_exp) - (1096200.0 * P4x * P4y * dP0z_exp) - (1058400.0 * P2x * P4y * dP2z_exp) + (37800.0 * P4y * dP4z_exp) + (389025.0 * P6x * P2y * dP0z_exp) + (741825.0 * P4x * P2y * dP2z_exp) + (316575.0 * P2x * P2y * dP4z_exp) - (36225.0 * P2y * dP6z_exp) - (12600.0 * P8x * dP0z_exp) - (36225.0 * P6x * dP2z_exp) - (33075.0 * P4x * dP4z_exp) - (7875.0 * P2x * dP6z_exp) + (1575.0 * dP8z_exp);

    double dterm_3_dx = (509040.0 * dP6x_exp * P2z) - (20160.0 * dP8x_exp) + (55440.0 * dP6x_exp * P2y) - (1096200.0 * dP4x_exp * P4z) - (1058400.0 * dP4x_exp * P2y * P2z) + (37800.0 * dP4x_exp * P4y) + (389025.0 * dP2x_exp * P6z) + (741825.0 * dP2x_exp * P2y * P4z) + (316575.0 * dP2x_exp * P4y * P2z) - (36225.0 * dP2x_exp * P6y) - (12600.0 * dP0x_exp * P8z) - (36225.0 * dP0x_exp * P2y * P6z) - (33075.0 * dP0x_exp * P4y * P4z) - (7875.0 * dP0x_exp * P6y * P2z) + (1575.0 * dP0x_exp * P8y);
    double dterm_3_dy = (509040.0 * P6x * dP0y_exp * P2z) - (20160.0 * P8x * dP0y_exp) + (55440.0 * P6x * dP2y_exp) - (1096200.0 * P4x * dP0y_exp * P4z) - (1058400.0 * P4x * dP2y_exp * P2z) + (37800.0 * P4x * dP4y_exp) + (389025.0 * P2x * dP0y_exp * P6z) + (741825.0 * P2x * dP2y_exp * P4z) + (316575.0 * P2x * dP4y_exp * P2z) - (36225.0 * P2x * dP6y_exp) - (12600.0 * dP0y_exp * P8z) - (36225.0 * dP2y_exp * P6z) - (33075.0 * dP4y_exp * P4z) - (7875.0 * dP6y_exp * P2z) + (1575.0 * dP8y_exp);
    double dterm_3_dz = (509040.0 * P6x * dP2z_exp) - (20160.0 * P8x * dP0z_exp) + (55440.0 * P6x * P2y * dP0z_exp) - (1096200.0 * P4x * dP4z_exp) - (1058400.0 * P4x * P2y * dP2z_exp) + (37800.0 * P4x * P4y * dP0z_exp) + (389025.0 * P2x * dP6z_exp) + (741825.0 * P2x * P2y * dP4z_exp) + (316575.0 * P2x * P4y * dP2z_exp) - (36225.0 * P2x * P6y * dP0z_exp) - (12600.0 * dP8z_exp) - (36225.0 * P2y * dP6z_exp) - (33075.0 * P4y * dP4z_exp) - (7875.0 * P6y * dP2z_exp) + (1575.0 * P8y * dP0z_exp);

    double dterm_4_dx = (509040.0 * dP2x_exp * P6z) - (20160.0 * dP0x_exp * P8z) + (55440.0 * dP0x_exp * P2y * P6z) - (1096200.0 * dP4x_exp * P4z) - (1058400.0 * dP2x_exp * P2y * P4z) + (37800.0 * dP0x_exp * P4y * P4z) + (389025.0 * dP6x_exp * P2z) + (741825.0 * dP4x_exp * P2y * P2z) + (316575.0 * dP2x_exp * P4y * P2z) - (36225.0 * dP0x_exp * P6y * P2z) - (12600.0 * dP8x_exp) - (36225.0 * dP6x_exp * P2y) - (33075.0 * dP4x_exp * P4y) - (7875.0 * dP2x_exp * P6y) + (1575.0 * dP0x_exp * P8y);
    double dterm_4_dy = (509040.0 * P2x * dP0y_exp * P6z) - (20160.0 * dP0y_exp * P8z) + (55440.0 * dP2y_exp * P6z) - (1096200.0 * P4x * dP0y_exp * P4z) - (1058400.0 * P2x * dP2y_exp * P4z) + (37800.0 * dP4y_exp * P4z) + (389025.0 * P6x * dP0y_exp * P2z) + (741825.0 * P4x * dP2y_exp * P2z) + (316575.0 * P2x * dP4y_exp * P2z) - (36225.0 * dP6y_exp * P2z) - (12600.0 * P8x * dP0y_exp) - (36225.0 * P6x * dP2y_exp) - (33075.0 * P4x * dP4y_exp) - (7875.0 * P2x * dP6y_exp) + (1575.0 * dP8y_exp);
    double dterm_4_dz = (509040.0 * P2x * dP6z_exp) - (20160.0 * dP8z_exp) + (55440.0 * P2y * dP6z_exp) - (1096200.0 * P4x * dP4z_exp) - (1058400.0 * P2x * P2y * dP4z_exp) + (37800.0 * P4y * dP4z_exp) + (389025.0 * P6x * dP2z_exp) + (741825.0 * P4x * P2y * dP2z_exp) + (316575.0 * P2x * P4y * dP2z_exp) - (36225.0 * P6y * dP2z_exp) - (12600.0 * P8x * dP0z_exp) - (36225.0 * P6x * P2y * dP0z_exp) - (33075.0 * P4x * P4y * dP0z_exp) - (7875.0 * P2x * P6y * dP0z_exp) + (1575.0 * P8y * dP0z_exp);

    double dterm_5_dx = (509040.0 * dP0x_exp * P6y * P2z) - (20160.0 * dP0x_exp * P8y) + (55440.0 * dP2x_exp * P6y) - (1096200.0 * dP0x_exp * P4y * P4z) - (1058400.0 * dP2x_exp * P4y * P2z) + (37800.0 * dP4x_exp * P4y) + (389025.0 * dP0x_exp * P2y * P6z) + (741825.0 * dP2x_exp * P2y * P4z) + (316575.0 * dP4x_exp * P2y * P2z) - (36225.0 * dP6x_exp * P2y) - (12600.0 * dP0x_exp * P8z) - (36225.0 * dP2x_exp * P6z) - (33075.0 * dP4x_exp * P4z) - (7875.0 * dP6x_exp * P2z) + (1575.0 * dP8x_exp);
    double dterm_5_dy = (509040.0 * dP6y_exp * P2z) - (20160.0 * dP8y_exp) + (55440.0 * P2x * dP6y_exp) - (1096200.0 * dP4y_exp * P4z) - (1058400.0 * P2x * dP4y_exp * P2z) + (37800.0 * P4x * dP4y_exp) + (389025.0 * dP2y_exp * P6z) + (741825.0 * P2x * dP2y_exp * P4z) + (316575.0 * P4x * dP2y_exp * P2z) - (36225.0 * P6x * dP2y_exp) - (12600.0 * dP0y_exp * P8z) - (36225.0 * P2x * dP0y_exp * P6z) - (33075.0 * P4x * dP0y_exp * P4z) - (7875.0 * P6x * dP0y_exp * P2z) + (1575.0 * P8x * dP0y_exp);
    double dterm_5_dz = (509040.0 * P6y * dP2z_exp) - (20160.0 * P8y * dP0z_exp) + (55440.0 * P2x * P6y * dP0z_exp) - (1096200.0 * P4y * dP4z_exp) - (1058400.0 * P2x * P4y * dP2z_exp) + (37800.0 * P4x * P4y * dP0z_exp) + (389025.0 * P2y * dP6z_exp) + (741825.0 * P2x * P2y * dP4z_exp) + (316575.0 * P4x * P2y * dP2z_exp) - (36225.0 * P6x * P2y * dP0z_exp) - (12600.0 * dP8z_exp) - (36225.0 * P2x * dP6z_exp) - (33075.0 * P4x * dP4z_exp) - (7875.0 * P6x * dP2z_exp) + (1575.0 * P8x * dP0z_exp);

    double dterm_6_dx = (509040.0 * dP0x_exp * P2y * P6z) - (20160.0 * dP0x_exp * P8z) + (55440.0 * dP2x_exp * P6z) - (1096200.0 * dP0x_exp * P4y * P4z) - (1058400.0 * dP2x_exp * P2y * P4z) + (37800.0 * dP4x_exp * P4z) + (389025.0 * dP0x_exp * P6y * P2z) + (741825.0 * dP2x_exp * P4y * P2z) + (316575.0 * dP4x_exp * P2y * P2z) - (36225.0 * dP6x_exp * P2z) - (12600.0 * dP0x_exp * P8y) - (36225.0 * dP2x_exp * P6y) - (33075.0 * dP4x_exp * P4y) - (7875.0 * dP6x_exp * P2y) + (1575.0 * dP8x_exp);
    double dterm_6_dy = (509040.0 * dP2y_exp * P6z) - (20160.0 * dP0y_exp * P8z) + (55440.0 * P2x * dP0y_exp * P6z) - (1096200.0 * dP4y_exp * P4z) - (1058400.0 * P2x * dP2y_exp * P4z) + (37800.0 * P4x * dP0y_exp * P4z) + (389025.0 * dP6y_exp * P2z) + (741825.0 * P2x * dP4y_exp * P2z) + (316575.0 * P4x * dP2y_exp * P2z) - (36225.0 * P6x * dP0y_exp * P2z) - (12600.0 * dP8y_exp) - (36225.0 * P2x * dP6y_exp) - (33075.0 * P4x * dP4y_exp) - (7875.0 * P6x * dP2y_exp) + (1575.0 * P8x * dP0y_exp);
    double dterm_6_dz = (509040.0 * P2y * dP6z_exp) - (20160.0 * dP8z_exp) + (55440.0 * P2x * dP6z_exp) - (1096200.0 * P4y * dP4z_exp) - (1058400.0 * P2x * P2y * dP4z_exp) + (37800.0 * P4x * dP4z_exp) + (389025.0 * P6y * dP2z_exp) + (741825.0 * P2x * P4y * dP2z_exp) + (316575.0 * P4x * P2y * dP2z_exp) - (36225.0 * P6x * dP2z_exp) - (12600.0 * P8y * dP0z_exp) - (36225.0 * P2x * P6y * dP0z_exp) - (33075.0 * P4x * P4y * dP0z_exp) - (7875.0 * P6x * P2y * dP0z_exp) + (1575.0 * P8x * dP0z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

    // dmiu4 dx/dy/dz
    deriv[9] = temp * dterm_4_dx;
    deriv[10] = temp * dterm_4_dy;
    deriv[11] = temp * dterm_4_dz;

    // dmiu5 dx/dy/dz
    deriv[12] = temp * dterm_5_dx;
    deriv[13] = temp * dterm_5_dy;
    deriv[14] = temp * dterm_5_dz;

    // dmiu6 dx/dy/dz
    deriv[15] = temp * dterm_6_dx;
    deriv[16] = temp * dterm_6_dy;
    deriv[17] = temp * dterm_6_dz;

}

void calc_solid_MCSH_8_4(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double term_1 = (453600.0 * P6x * P1y * P1z) - (1134000.0 * P4x * P3y * P1z) - (1134000.0 * P4x * P1y * P3z) + (425250.0 * P2x * P5y * P1z) + (850500.0 * P2x * P3y * P3z) + (425250.0 * P2x * P1y * P5z) - (14175.0 * P7y * P1z) - (42525.0 * P5y * P3z) - (42525.0 * P3y * P5z) - (14175.0 * P1y * P7z);
    double term_2 = (453600.0 * P1x * P6y * P1z) - (1134000.0 * P3x * P4y * P1z) - (1134000.0 * P1x * P4y * P3z) + (425250.0 * P5x * P2y * P1z) + (850500.0 * P3x * P2y * P3z) + (425250.0 * P1x * P2y * P5z) - (14175.0 * P7x * P1z) - (42525.0 * P5x * P3z) - (42525.0 * P3x * P5z) - (14175.0 * P1x * P7z);
    double term_3 = (453600.0 * P1x * P1y * P6z) - (1134000.0 * P3x * P1y * P4z) - (1134000.0 * P1x * P3y * P4z) + (425250.0 * P5x * P1y * P2z) + (850500.0 * P3x * P3y * P2z) + (425250.0 * P1x * P5y * P2z) - (14175.0 * P7x * P1y) - (42525.0 * P5x * P3y) - (42525.0 * P3x * P5y) - (14175.0 * P1x * P7y);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP3x_exp = dP3_exp(P3x, C2, lambda, x0, gamma);
    double dP3y_exp = dP3_exp(P3y, C2, lambda, y0, gamma);
    double dP3z_exp = dP3_exp(P3z, C2, lambda, z0, gamma);

    double dP4x_exp = dP4_exp(P4x, C2, lambda, x0, gamma);
    double dP4y_exp = dP4_exp(P4y, C2, lambda, y0, gamma);
    double dP4z_exp = dP4_exp(P4z, C2, lambda, z0, gamma);

    double dP5x_exp = dP5_exp(P5x, C2, lambda, x0, gamma);
    double dP5y_exp = dP5_exp(P5y, C2, lambda, y0, gamma);
    double dP5z_exp = dP5_exp(P5z, C2, lambda, z0, gamma);

    double dP6x_exp = dP6_exp(P6x, C2, lambda, x0, gamma);
    double dP6y_exp = dP6_exp(P6y, C2, lambda, y0, gamma);
    double dP6z_exp = dP6_exp(P6z, C2, lambda, z0, gamma);

    double dP7x_exp = dP7_exp(P7x, C2, lambda, x0, gamma);
    double dP7y_exp = dP7_exp(P7y, C2, lambda, y0, gamma);
    double dP7z_exp = dP7_exp(P7z, C2, lambda, z0, gamma);

    double dterm_1_dx = (453600.0 * dP6x_exp * P1y * P1z) - (1134000.0 * dP4x_exp * P3y * P1z) - (1134000.0 * dP4x_exp * P1y * P3z) + (425250.0 * dP2x_exp * P5y * P1z) + (850500.0 * dP2x_exp * P3y * P3z) + (425250.0 * dP2x_exp * P1y * P5z) - (14175.0 * dP0x_exp * P7y * P1z) - (42525.0 * dP0x_exp * P5y * P3z) - (42525.0 * dP0x_exp * P3y * P5z) - (14175.0 * dP0x_exp * P1y * P7z);
    double dterm_1_dy = (453600.0 * P6x * dP1y_exp * P1z) - (1134000.0 * P4x * dP3y_exp * P1z) - (1134000.0 * P4x * dP1y_exp * P3z) + (425250.0 * P2x * dP5y_exp * P1z) + (850500.0 * P2x * dP3y_exp * P3z) + (425250.0 * P2x * dP1y_exp * P5z) - (14175.0 * dP7y_exp * P1z) - (42525.0 * dP5y_exp * P3z) - (42525.0 * dP3y_exp * P5z) - (14175.0 * dP1y_exp * P7z);
    double dterm_1_dz = (453600.0 * P6x * P1y * dP1z_exp) - (1134000.0 * P4x * P3y * dP1z_exp) - (1134000.0 * P4x * P1y * dP3z_exp) + (425250.0 * P2x * P5y * dP1z_exp) + (850500.0 * P2x * P3y * dP3z_exp) + (425250.0 * P2x * P1y * dP5z_exp) - (14175.0 * P7y * dP1z_exp) - (42525.0 * P5y * dP3z_exp) - (42525.0 * P3y * dP5z_exp) - (14175.0 * P1y * dP7z_exp);

    double dterm_2_dx = (453600.0 * dP1x_exp * P6y * P1z) - (1134000.0 * dP3x_exp * P4y * P1z) - (1134000.0 * dP1x_exp * P4y * P3z) + (425250.0 * dP5x_exp * P2y * P1z) + (850500.0 * dP3x_exp * P2y * P3z) + (425250.0 * dP1x_exp * P2y * P5z) - (14175.0 * dP7x_exp * P1z) - (42525.0 * dP5x_exp * P3z) - (42525.0 * dP3x_exp * P5z) - (14175.0 * dP1x_exp * P7z);
    double dterm_2_dy = (453600.0 * P1x * dP6y_exp * P1z) - (1134000.0 * P3x * dP4y_exp * P1z) - (1134000.0 * P1x * dP4y_exp * P3z) + (425250.0 * P5x * dP2y_exp * P1z) + (850500.0 * P3x * dP2y_exp * P3z) + (425250.0 * P1x * dP2y_exp * P5z) - (14175.0 * P7x * dP0y_exp * P1z) - (42525.0 * P5x * dP0y_exp * P3z) - (42525.0 * P3x * dP0y_exp * P5z) - (14175.0 * P1x * dP0y_exp * P7z);
    double dterm_2_dz = (453600.0 * P1x * P6y * dP1z_exp) - (1134000.0 * P3x * P4y * dP1z_exp) - (1134000.0 * P1x * P4y * dP3z_exp) + (425250.0 * P5x * P2y * dP1z_exp) + (850500.0 * P3x * P2y * dP3z_exp) + (425250.0 * P1x * P2y * dP5z_exp) - (14175.0 * P7x * dP1z_exp) - (42525.0 * P5x * dP3z_exp) - (42525.0 * P3x * dP5z_exp) - (14175.0 * P1x * dP7z_exp);

    double dterm_3_dx = (453600.0 * dP1x_exp * P1y * P6z) - (1134000.0 * dP3x_exp * P1y * P4z) - (1134000.0 * dP1x_exp * P3y * P4z) + (425250.0 * dP5x_exp * P1y * P2z) + (850500.0 * dP3x_exp * P3y * P2z) + (425250.0 * dP1x_exp * P5y * P2z) - (14175.0 * dP7x_exp * P1y) - (42525.0 * dP5x_exp * P3y) - (42525.0 * dP3x_exp * P5y) - (14175.0 * dP1x_exp * P7y);
    double dterm_3_dy = (453600.0 * P1x * dP1y_exp * P6z) - (1134000.0 * P3x * dP1y_exp * P4z) - (1134000.0 * P1x * dP3y_exp * P4z) + (425250.0 * P5x * dP1y_exp * P2z) + (850500.0 * P3x * dP3y_exp * P2z) + (425250.0 * P1x * dP5y_exp * P2z) - (14175.0 * P7x * dP1y_exp) - (42525.0 * P5x * dP3y_exp) - (42525.0 * P3x * dP5y_exp) - (14175.0 * P1x * dP7y_exp);
    double dterm_3_dz = (453600.0 * P1x * P1y * dP6z_exp) - (1134000.0 * P3x * P1y * dP4z_exp) - (1134000.0 * P1x * P3y * dP4z_exp) + (425250.0 * P5x * P1y * dP2z_exp) + (850500.0 * P3x * P3y * dP2z_exp) + (425250.0 * P1x * P5y * dP2z_exp) - (14175.0 * P7x * P1y * dP0z_exp) - (42525.0 * P5x * P3y * dP0z_exp) - (42525.0 * P3x * P5y * dP0z_exp) - (14175.0 * P1x * P7y * dP0z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

}

void calc_solid_MCSH_8_5(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double term_1 = (922320.0 * P5x * P3y) - (136080.0 * P7x * P1y) + (90720.0 * P5x * P1y * P2z) - (855225.0 * P3x * P5y) - (670950.0 * P3x * P3y * P2z) + (184275.0 * P3x * P1y * P4z) + (113400.0 * P1x * P7y) + (184275.0 * P1x * P5y * P2z) + (28350.0 * P1x * P3y * P4z) - (42525.0 * P1x * P1y * P6z);
    double term_2 = (922320.0 * P3x * P5y) - (136080.0 * P1x * P7y) + (90720.0 * P1x * P5y * P2z) - (855225.0 * P5x * P3y) - (670950.0 * P3x * P3y * P2z) + (184275.0 * P1x * P3y * P4z) + (113400.0 * P7x * P1y) + (184275.0 * P5x * P1y * P2z) + (28350.0 * P3x * P1y * P4z) - (42525.0 * P1x * P1y * P6z);
    double term_3 = (922320.0 * P5x * P3z) - (136080.0 * P7x * P1z) + (90720.0 * P5x * P2y * P1z) - (855225.0 * P3x * P5z) - (670950.0 * P3x * P2y * P3z) + (184275.0 * P3x * P4y * P1z) + (113400.0 * P1x * P7z) + (184275.0 * P1x * P2y * P5z) + (28350.0 * P1x * P4y * P3z) - (42525.0 * P1x * P6y * P1z);
    double term_4 = (922320.0 * P3x * P5z) - (136080.0 * P1x * P7z) + (90720.0 * P1x * P2y * P5z) - (855225.0 * P5x * P3z) - (670950.0 * P3x * P2y * P3z) + (184275.0 * P1x * P4y * P3z) + (113400.0 * P7x * P1z) + (184275.0 * P5x * P2y * P1z) + (28350.0 * P3x * P4y * P1z) - (42525.0 * P1x * P6y * P1z);
    double term_5 = (922320.0 * P5y * P3z) - (136080.0 * P7y * P1z) + (90720.0 * P2x * P5y * P1z) - (855225.0 * P3y * P5z) - (670950.0 * P2x * P3y * P3z) + (184275.0 * P4x * P3y * P1z) + (113400.0 * P1y * P7z) + (184275.0 * P2x * P1y * P5z) + (28350.0 * P4x * P1y * P3z) - (42525.0 * P6x * P1y * P1z);
    double term_6 = (922320.0 * P3y * P5z) - (136080.0 * P1y * P7z) + (90720.0 * P2x * P1y * P5z) - (855225.0 * P5y * P3z) - (670950.0 * P2x * P3y * P3z) + (184275.0 * P4x * P1y * P3z) + (113400.0 * P7y * P1z) + (184275.0 * P2x * P5y * P1z) + (28350.0 * P4x * P3y * P1z) - (42525.0 * P6x * P1y * P1z);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP3x_exp = dP3_exp(P3x, C2, lambda, x0, gamma);
    double dP3y_exp = dP3_exp(P3y, C2, lambda, y0, gamma);
    double dP3z_exp = dP3_exp(P3z, C2, lambda, z0, gamma);

    double dP4x_exp = dP4_exp(P4x, C2, lambda, x0, gamma);
    double dP4y_exp = dP4_exp(P4y, C2, lambda, y0, gamma);
    double dP4z_exp = dP4_exp(P4z, C2, lambda, z0, gamma);

    double dP5x_exp = dP5_exp(P5x, C2, lambda, x0, gamma);
    double dP5y_exp = dP5_exp(P5y, C2, lambda, y0, gamma);
    double dP5z_exp = dP5_exp(P5z, C2, lambda, z0, gamma);

    double dP6x_exp = dP6_exp(P6x, C2, lambda, x0, gamma);
    double dP6y_exp = dP6_exp(P6y, C2, lambda, y0, gamma);
    double dP6z_exp = dP6_exp(P6z, C2, lambda, z0, gamma);

    double dP7x_exp = dP7_exp(P7x, C2, lambda, x0, gamma);
    double dP7y_exp = dP7_exp(P7y, C2, lambda, y0, gamma);
    double dP7z_exp = dP7_exp(P7z, C2, lambda, z0, gamma);

    double dterm_1_dx = (922320.0 * dP5x_exp * P3y) - (136080.0 * dP7x_exp * P1y) + (90720.0 * dP5x_exp * P1y * P2z) - (855225.0 * dP3x_exp * P5y) - (670950.0 * dP3x_exp * P3y * P2z) + (184275.0 * dP3x_exp * P1y * P4z) + (113400.0 * dP1x_exp * P7y) + (184275.0 * dP1x_exp * P5y * P2z) + (28350.0 * dP1x_exp * P3y * P4z) - (42525.0 * dP1x_exp * P1y * P6z);
    double dterm_1_dy = (922320.0 * P5x * dP3y_exp) - (136080.0 * P7x * dP1y_exp) + (90720.0 * P5x * dP1y_exp * P2z) - (855225.0 * P3x * dP5y_exp) - (670950.0 * P3x * dP3y_exp * P2z) + (184275.0 * P3x * dP1y_exp * P4z) + (113400.0 * P1x * dP7y_exp) + (184275.0 * P1x * dP5y_exp * P2z) + (28350.0 * P1x * dP3y_exp * P4z) - (42525.0 * P1x * dP1y_exp * P6z);
    double dterm_1_dz = (922320.0 * P5x * P3y * dP0z_exp) - (136080.0 * P7x * P1y * dP0z_exp) + (90720.0 * P5x * P1y * dP2z_exp) - (855225.0 * P3x * P5y * dP0z_exp) - (670950.0 * P3x * P3y * dP2z_exp) + (184275.0 * P3x * P1y * dP4z_exp) + (113400.0 * P1x * P7y * dP0z_exp) + (184275.0 * P1x * P5y * dP2z_exp) + (28350.0 * P1x * P3y * dP4z_exp) - (42525.0 * P1x * P1y * dP6z_exp);

    double dterm_2_dx = (922320.0 * dP3x_exp * P5y) - (136080.0 * dP1x_exp * P7y) + (90720.0 * dP1x_exp * P5y * P2z) - (855225.0 * dP5x_exp * P3y) - (670950.0 * dP3x_exp * P3y * P2z) + (184275.0 * dP1x_exp * P3y * P4z) + (113400.0 * dP7x_exp * P1y) + (184275.0 * dP5x_exp * P1y * P2z) + (28350.0 * dP3x_exp * P1y * P4z) - (42525.0 * dP1x_exp * P1y * P6z);
    double dterm_2_dy = (922320.0 * P3x * dP5y_exp) - (136080.0 * P1x * dP7y_exp) + (90720.0 * P1x * dP5y_exp * P2z) - (855225.0 * P5x * dP3y_exp) - (670950.0 * P3x * dP3y_exp * P2z) + (184275.0 * P1x * dP3y_exp * P4z) + (113400.0 * P7x * dP1y_exp) + (184275.0 * P5x * dP1y_exp * P2z) + (28350.0 * P3x * dP1y_exp * P4z) - (42525.0 * P1x * dP1y_exp * P6z);
    double dterm_2_dz = (922320.0 * P3x * P5y * dP0z_exp) - (136080.0 * P1x * P7y * dP0z_exp) + (90720.0 * P1x * P5y * dP2z_exp) - (855225.0 * P5x * P3y * dP0z_exp) - (670950.0 * P3x * P3y * dP2z_exp) + (184275.0 * P1x * P3y * dP4z_exp) + (113400.0 * P7x * P1y * dP0z_exp) + (184275.0 * P5x * P1y * dP2z_exp) + (28350.0 * P3x * P1y * dP4z_exp) - (42525.0 * P1x * P1y * dP6z_exp);

    double dterm_3_dx = (922320.0 * dP5x_exp * P3z) - (136080.0 * dP7x_exp * P1z) + (90720.0 * dP5x_exp * P2y * P1z) - (855225.0 * dP3x_exp * P5z) - (670950.0 * dP3x_exp * P2y * P3z) + (184275.0 * dP3x_exp * P4y * P1z) + (113400.0 * dP1x_exp * P7z) + (184275.0 * dP1x_exp * P2y * P5z) + (28350.0 * dP1x_exp * P4y * P3z) - (42525.0 * dP1x_exp * P6y * P1z);
    double dterm_3_dy = (922320.0 * P5x * dP0y_exp * P3z) - (136080.0 * P7x * dP0y_exp * P1z) + (90720.0 * P5x * dP2y_exp * P1z) - (855225.0 * P3x * dP0y_exp * P5z) - (670950.0 * P3x * dP2y_exp * P3z) + (184275.0 * P3x * dP4y_exp * P1z) + (113400.0 * P1x * dP0y_exp * P7z) + (184275.0 * P1x * dP2y_exp * P5z) + (28350.0 * P1x * dP4y_exp * P3z) - (42525.0 * P1x * dP6y_exp * P1z);
    double dterm_3_dz = (922320.0 * P5x * dP3z_exp) - (136080.0 * P7x * dP1z_exp) + (90720.0 * P5x * P2y * dP1z_exp) - (855225.0 * P3x * dP5z_exp) - (670950.0 * P3x * P2y * dP3z_exp) + (184275.0 * P3x * P4y * dP1z_exp) + (113400.0 * P1x * dP7z_exp) + (184275.0 * P1x * P2y * dP5z_exp) + (28350.0 * P1x * P4y * dP3z_exp) - (42525.0 * P1x * P6y * dP1z_exp);

    double dterm_4_dx = (922320.0 * dP3x_exp * P5z) - (136080.0 * dP1x_exp * P7z) + (90720.0 * dP1x_exp * P2y * P5z) - (855225.0 * dP5x_exp * P3z) - (670950.0 * dP3x_exp * P2y * P3z) + (184275.0 * dP1x_exp * P4y * P3z) + (113400.0 * dP7x_exp * P1z) + (184275.0 * dP5x_exp * P2y * P1z) + (28350.0 * dP3x_exp * P4y * P1z) - (42525.0 * dP1x_exp * P6y * P1z);
    double dterm_4_dy = (922320.0 * P3x * dP0y_exp * P5z) - (136080.0 * P1x * dP0y_exp * P7z) + (90720.0 * P1x * dP2y_exp * P5z) - (855225.0 * P5x * dP0y_exp * P3z) - (670950.0 * P3x * dP2y_exp * P3z) + (184275.0 * P1x * dP4y_exp * P3z) + (113400.0 * P7x * dP0y_exp * P1z) + (184275.0 * P5x * dP2y_exp * P1z) + (28350.0 * P3x * dP4y_exp * P1z) - (42525.0 * P1x * dP6y_exp * P1z);
    double dterm_4_dz = (922320.0 * P3x * dP5z_exp) - (136080.0 * P1x * dP7z_exp) + (90720.0 * P1x * P2y * dP5z_exp) - (855225.0 * P5x * dP3z_exp) - (670950.0 * P3x * P2y * dP3z_exp) + (184275.0 * P1x * P4y * dP3z_exp) + (113400.0 * P7x * dP1z_exp) + (184275.0 * P5x * P2y * dP1z_exp) + (28350.0 * P3x * P4y * dP1z_exp) - (42525.0 * P1x * P6y * dP1z_exp);

    double dterm_5_dx = (922320.0 * dP0x_exp * P5y * P3z) - (136080.0 * dP0x_exp * P7y * P1z) + (90720.0 * dP2x_exp * P5y * P1z) - (855225.0 * dP0x_exp * P3y * P5z) - (670950.0 * dP2x_exp * P3y * P3z) + (184275.0 * dP4x_exp * P3y * P1z) + (113400.0 * dP0x_exp * P1y * P7z) + (184275.0 * dP2x_exp * P1y * P5z) + (28350.0 * dP4x_exp * P1y * P3z) - (42525.0 * dP6x_exp * P1y * P1z);
    double dterm_5_dy = (922320.0 * dP5y_exp * P3z) - (136080.0 * dP7y_exp * P1z) + (90720.0 * P2x * dP5y_exp * P1z) - (855225.0 * dP3y_exp * P5z) - (670950.0 * P2x * dP3y_exp * P3z) + (184275.0 * P4x * dP3y_exp * P1z) + (113400.0 * dP1y_exp * P7z) + (184275.0 * P2x * dP1y_exp * P5z) + (28350.0 * P4x * dP1y_exp * P3z) - (42525.0 * P6x * dP1y_exp * P1z);
    double dterm_5_dz = (922320.0 * P5y * dP3z_exp) - (136080.0 * P7y * dP1z_exp) + (90720.0 * P2x * P5y * dP1z_exp) - (855225.0 * P3y * dP5z_exp) - (670950.0 * P2x * P3y * dP3z_exp) + (184275.0 * P4x * P3y * dP1z_exp) + (113400.0 * P1y * dP7z_exp) + (184275.0 * P2x * P1y * dP5z_exp) + (28350.0 * P4x * P1y * dP3z_exp) - (42525.0 * P6x * P1y * dP1z_exp);

    double dterm_6_dx = (922320.0 * dP0x_exp * P3y * P5z) - (136080.0 * dP0x_exp * P1y * P7z) + (90720.0 * dP2x_exp * P1y * P5z) - (855225.0 * dP0x_exp * P5y * P3z) - (670950.0 * dP2x_exp * P3y * P3z) + (184275.0 * dP4x_exp * P1y * P3z) + (113400.0 * dP0x_exp * P7y * P1z) + (184275.0 * dP2x_exp * P5y * P1z) + (28350.0 * dP4x_exp * P3y * P1z) - (42525.0 * dP6x_exp * P1y * P1z);
    double dterm_6_dy = (922320.0 * dP3y_exp * P5z) - (136080.0 * dP1y_exp * P7z) + (90720.0 * P2x * dP1y_exp * P5z) - (855225.0 * dP5y_exp * P3z) - (670950.0 * P2x * dP3y_exp * P3z) + (184275.0 * P4x * dP1y_exp * P3z) + (113400.0 * dP7y_exp * P1z) + (184275.0 * P2x * dP5y_exp * P1z) + (28350.0 * P4x * dP3y_exp * P1z) - (42525.0 * P6x * dP1y_exp * P1z);
    double dterm_6_dz = (922320.0 * P3y * dP5z_exp) - (136080.0 * P1y * dP7z_exp) + (90720.0 * P2x * P1y * dP5z_exp) - (855225.0 * P5y * dP3z_exp) - (670950.0 * P2x * P3y * dP3z_exp) + (184275.0 * P4x * P1y * dP3z_exp) + (113400.0 * P7y * dP1z_exp) + (184275.0 * P2x * P5y * dP1z_exp) + (28350.0 * P4x * P3y * dP1z_exp) - (42525.0 * P6x * P1y * dP1z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

    // dmiu4 dx/dy/dz
    deriv[9] = temp * dterm_4_dx;
    deriv[10] = temp * dterm_4_dy;
    deriv[11] = temp * dterm_4_dz;

    // dmiu5 dx/dy/dz
    deriv[12] = temp * dterm_5_dx;
    deriv[13] = temp * dterm_5_dy;
    deriv[14] = temp * dterm_5_dz;

    // dmiu6 dx/dy/dz
    deriv[15] = temp * dterm_6_dx;
    deriv[16] = temp * dterm_6_dy;
    deriv[17] = temp * dterm_6_dz;

}

void calc_solid_MCSH_8_6(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double term_1 = (861840.0 * P5x * P2y * P1z) - (45360.0 * P7x * P1z) + (30240.0 * P5x * P3z) - (978075.0 * P3x * P4y * P1z) - (916650.0 * P3x * P2y * P3z) + (61425.0 * P3x * P5z) + (141750.0 * P1x * P6y * P1z) + (269325.0 * P1x * P4y * P3z) + (113400.0 * P1x * P2y * P5z) - (14175.0 * P1x * P7z);
    double term_2 = (861840.0 * P2x * P5y * P1z) - (45360.0 * P7y * P1z) + (30240.0 * P5y * P3z) - (978075.0 * P4x * P3y * P1z) - (916650.0 * P2x * P3y * P3z) + (61425.0 * P3y * P5z) + (141750.0 * P6x * P1y * P1z) + (269325.0 * P4x * P1y * P3z) + (113400.0 * P2x * P1y * P5z) - (14175.0 * P1y * P7z);
    double term_3 = (861840.0 * P5x * P1y * P2z) - (45360.0 * P7x * P1y) + (30240.0 * P5x * P3y) - (978075.0 * P3x * P1y * P4z) - (916650.0 * P3x * P3y * P2z) + (61425.0 * P3x * P5y) + (141750.0 * P1x * P1y * P6z) + (269325.0 * P1x * P3y * P4z) + (113400.0 * P1x * P5y * P2z) - (14175.0 * P1x * P7y);
    double term_4 = (861840.0 * P2x * P1y * P5z) - (45360.0 * P1y * P7z) + (30240.0 * P3y * P5z) - (978075.0 * P4x * P1y * P3z) - (916650.0 * P2x * P3y * P3z) + (61425.0 * P5y * P3z) + (141750.0 * P6x * P1y * P1z) + (269325.0 * P4x * P3y * P1z) + (113400.0 * P2x * P5y * P1z) - (14175.0 * P7y * P1z);
    double term_5 = (861840.0 * P1x * P5y * P2z) - (45360.0 * P1x * P7y) + (30240.0 * P3x * P5y) - (978075.0 * P1x * P3y * P4z) - (916650.0 * P3x * P3y * P2z) + (61425.0 * P5x * P3y) + (141750.0 * P1x * P1y * P6z) + (269325.0 * P3x * P1y * P4z) + (113400.0 * P5x * P1y * P2z) - (14175.0 * P7x * P1y);
    double term_6 = (861840.0 * P1x * P2y * P5z) - (45360.0 * P1x * P7z) + (30240.0 * P3x * P5z) - (978075.0 * P1x * P4y * P3z) - (916650.0 * P3x * P2y * P3z) + (61425.0 * P5x * P3z) + (141750.0 * P1x * P6y * P1z) + (269325.0 * P3x * P4y * P1z) + (113400.0 * P5x * P2y * P1z) - (14175.0 * P7x * P1z);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP3x_exp = dP3_exp(P3x, C2, lambda, x0, gamma);
    double dP3y_exp = dP3_exp(P3y, C2, lambda, y0, gamma);
    double dP3z_exp = dP3_exp(P3z, C2, lambda, z0, gamma);

    double dP4x_exp = dP4_exp(P4x, C2, lambda, x0, gamma);
    double dP4y_exp = dP4_exp(P4y, C2, lambda, y0, gamma);
    double dP4z_exp = dP4_exp(P4z, C2, lambda, z0, gamma);

    double dP5x_exp = dP5_exp(P5x, C2, lambda, x0, gamma);
    double dP5y_exp = dP5_exp(P5y, C2, lambda, y0, gamma);
    double dP5z_exp = dP5_exp(P5z, C2, lambda, z0, gamma);

    double dP6x_exp = dP6_exp(P6x, C2, lambda, x0, gamma);
    double dP6y_exp = dP6_exp(P6y, C2, lambda, y0, gamma);
    double dP6z_exp = dP6_exp(P6z, C2, lambda, z0, gamma);

    double dP7x_exp = dP7_exp(P7x, C2, lambda, x0, gamma);
    double dP7y_exp = dP7_exp(P7y, C2, lambda, y0, gamma);
    double dP7z_exp = dP7_exp(P7z, C2, lambda, z0, gamma);

    double dterm_1_dx = (861840.0 * dP5x_exp * P2y * P1z) - (45360.0 * dP7x_exp * P1z) + (30240.0 * dP5x_exp * P3z) - (978075.0 * dP3x_exp * P4y * P1z) - (916650.0 * dP3x_exp * P2y * P3z) + (61425.0 * dP3x_exp * P5z) + (141750.0 * dP1x_exp * P6y * P1z) + (269325.0 * dP1x_exp * P4y * P3z) + (113400.0 * dP1x_exp * P2y * P5z) - (14175.0 * dP1x_exp * P7z);
    double dterm_1_dy = (861840.0 * P5x * dP2y_exp * P1z) - (45360.0 * P7x * dP0y_exp * P1z) + (30240.0 * P5x * dP0y_exp * P3z) - (978075.0 * P3x * dP4y_exp * P1z) - (916650.0 * P3x * dP2y_exp * P3z) + (61425.0 * P3x * dP0y_exp * P5z) + (141750.0 * P1x * dP6y_exp * P1z) + (269325.0 * P1x * dP4y_exp * P3z) + (113400.0 * P1x * dP2y_exp * P5z) - (14175.0 * P1x * dP0y_exp * P7z);
    double dterm_1_dz = (861840.0 * P5x * P2y * dP1z_exp) - (45360.0 * P7x * dP1z_exp) + (30240.0 * P5x * dP3z_exp) - (978075.0 * P3x * P4y * dP1z_exp) - (916650.0 * P3x * P2y * dP3z_exp) + (61425.0 * P3x * dP5z_exp) + (141750.0 * P1x * P6y * dP1z_exp) + (269325.0 * P1x * P4y * dP3z_exp) + (113400.0 * P1x * P2y * dP5z_exp) - (14175.0 * P1x * dP7z_exp);

    double dterm_2_dx = (861840.0 * dP2x_exp * P5y * P1z) - (45360.0 * dP0x_exp * P7y * P1z) + (30240.0 * dP0x_exp * P5y * P3z) - (978075.0 * dP4x_exp * P3y * P1z) - (916650.0 * dP2x_exp * P3y * P3z) + (61425.0 * dP0x_exp * P3y * P5z) + (141750.0 * dP6x_exp * P1y * P1z) + (269325.0 * dP4x_exp * P1y * P3z) + (113400.0 * dP2x_exp * P1y * P5z) - (14175.0 * dP0x_exp * P1y * P7z);
    double dterm_2_dy = (861840.0 * P2x * dP5y_exp * P1z) - (45360.0 * dP7y_exp * P1z) + (30240.0 * dP5y_exp * P3z) - (978075.0 * P4x * dP3y_exp * P1z) - (916650.0 * P2x * dP3y_exp * P3z) + (61425.0 * dP3y_exp * P5z) + (141750.0 * P6x * dP1y_exp * P1z) + (269325.0 * P4x * dP1y_exp * P3z) + (113400.0 * P2x * dP1y_exp * P5z) - (14175.0 * dP1y_exp * P7z);
    double dterm_2_dz = (861840.0 * P2x * P5y * dP1z_exp) - (45360.0 * P7y * dP1z_exp) + (30240.0 * P5y * dP3z_exp) - (978075.0 * P4x * P3y * dP1z_exp) - (916650.0 * P2x * P3y * dP3z_exp) + (61425.0 * P3y * dP5z_exp) + (141750.0 * P6x * P1y * dP1z_exp) + (269325.0 * P4x * P1y * dP3z_exp) + (113400.0 * P2x * P1y * dP5z_exp) - (14175.0 * P1y * dP7z_exp);

    double dterm_3_dx = (861840.0 * dP5x_exp * P1y * P2z) - (45360.0 * dP7x_exp * P1y) + (30240.0 * dP5x_exp * P3y) - (978075.0 * dP3x_exp * P1y * P4z) - (916650.0 * dP3x_exp * P3y * P2z) + (61425.0 * dP3x_exp * P5y) + (141750.0 * dP1x_exp * P1y * P6z) + (269325.0 * dP1x_exp * P3y * P4z) + (113400.0 * dP1x_exp * P5y * P2z) - (14175.0 * dP1x_exp * P7y);
    double dterm_3_dy = (861840.0 * P5x * dP1y_exp * P2z) - (45360.0 * P7x * dP1y_exp) + (30240.0 * P5x * dP3y_exp) - (978075.0 * P3x * dP1y_exp * P4z) - (916650.0 * P3x * dP3y_exp * P2z) + (61425.0 * P3x * dP5y_exp) + (141750.0 * P1x * dP1y_exp * P6z) + (269325.0 * P1x * dP3y_exp * P4z) + (113400.0 * P1x * dP5y_exp * P2z) - (14175.0 * P1x * dP7y_exp);
    double dterm_3_dz = (861840.0 * P5x * P1y * dP2z_exp) - (45360.0 * P7x * P1y * dP0z_exp) + (30240.0 * P5x * P3y * dP0z_exp) - (978075.0 * P3x * P1y * dP4z_exp) - (916650.0 * P3x * P3y * dP2z_exp) + (61425.0 * P3x * P5y * dP0z_exp) + (141750.0 * P1x * P1y * dP6z_exp) + (269325.0 * P1x * P3y * dP4z_exp) + (113400.0 * P1x * P5y * dP2z_exp) - (14175.0 * P1x * P7y * dP0z_exp);

    double dterm_4_dx = (861840.0 * dP2x_exp * P1y * P5z) - (45360.0 * dP0x_exp * P1y * P7z) + (30240.0 * dP0x_exp * P3y * P5z) - (978075.0 * dP4x_exp * P1y * P3z) - (916650.0 * dP2x_exp * P3y * P3z) + (61425.0 * dP0x_exp * P5y * P3z) + (141750.0 * dP6x_exp * P1y * P1z) + (269325.0 * dP4x_exp * P3y * P1z) + (113400.0 * dP2x_exp * P5y * P1z) - (14175.0 * dP0x_exp * P7y * P1z);
    double dterm_4_dy = (861840.0 * P2x * dP1y_exp * P5z) - (45360.0 * dP1y_exp * P7z) + (30240.0 * dP3y_exp * P5z) - (978075.0 * P4x * dP1y_exp * P3z) - (916650.0 * P2x * dP3y_exp * P3z) + (61425.0 * dP5y_exp * P3z) + (141750.0 * P6x * dP1y_exp * P1z) + (269325.0 * P4x * dP3y_exp * P1z) + (113400.0 * P2x * dP5y_exp * P1z) - (14175.0 * dP7y_exp * P1z);
    double dterm_4_dz = (861840.0 * P2x * P1y * dP5z_exp) - (45360.0 * P1y * dP7z_exp) + (30240.0 * P3y * dP5z_exp) - (978075.0 * P4x * P1y * dP3z_exp) - (916650.0 * P2x * P3y * dP3z_exp) + (61425.0 * P5y * dP3z_exp) + (141750.0 * P6x * P1y * dP1z_exp) + (269325.0 * P4x * P3y * dP1z_exp) + (113400.0 * P2x * P5y * dP1z_exp) - (14175.0 * P7y * dP1z_exp);

    double dterm_5_dx = (861840.0 * dP1x_exp * P5y * P2z) - (45360.0 * dP1x_exp * P7y) + (30240.0 * dP3x_exp * P5y) - (978075.0 * dP1x_exp * P3y * P4z) - (916650.0 * dP3x_exp * P3y * P2z) + (61425.0 * dP5x_exp * P3y) + (141750.0 * dP1x_exp * P1y * P6z) + (269325.0 * dP3x_exp * P1y * P4z) + (113400.0 * dP5x_exp * P1y * P2z) - (14175.0 * dP7x_exp * P1y);
    double dterm_5_dy = (861840.0 * P1x * dP5y_exp * P2z) - (45360.0 * P1x * dP7y_exp) + (30240.0 * P3x * dP5y_exp) - (978075.0 * P1x * dP3y_exp * P4z) - (916650.0 * P3x * dP3y_exp * P2z) + (61425.0 * P5x * dP3y_exp) + (141750.0 * P1x * dP1y_exp * P6z) + (269325.0 * P3x * dP1y_exp * P4z) + (113400.0 * P5x * dP1y_exp * P2z) - (14175.0 * P7x * dP1y_exp);
    double dterm_5_dz = (861840.0 * P1x * P5y * dP2z_exp) - (45360.0 * P1x * P7y * dP0z_exp) + (30240.0 * P3x * P5y * dP0z_exp) - (978075.0 * P1x * P3y * dP4z_exp) - (916650.0 * P3x * P3y * dP2z_exp) + (61425.0 * P5x * P3y * dP0z_exp) + (141750.0 * P1x * P1y * dP6z_exp) + (269325.0 * P3x * P1y * dP4z_exp) + (113400.0 * P5x * P1y * dP2z_exp) - (14175.0 * P7x * P1y * dP0z_exp);

    double dterm_6_dx = (861840.0 * dP1x_exp * P2y * P5z) - (45360.0 * dP1x_exp * P7z) + (30240.0 * dP3x_exp * P5z) - (978075.0 * dP1x_exp * P4y * P3z) - (916650.0 * dP3x_exp * P2y * P3z) + (61425.0 * dP5x_exp * P3z) + (141750.0 * dP1x_exp * P6y * P1z) + (269325.0 * dP3x_exp * P4y * P1z) + (113400.0 * dP5x_exp * P2y * P1z) - (14175.0 * dP7x_exp * P1z);
    double dterm_6_dy = (861840.0 * P1x * dP2y_exp * P5z) - (45360.0 * P1x * dP0y_exp * P7z) + (30240.0 * P3x * dP0y_exp * P5z) - (978075.0 * P1x * dP4y_exp * P3z) - (916650.0 * P3x * dP2y_exp * P3z) + (61425.0 * P5x * dP0y_exp * P3z) + (141750.0 * P1x * dP6y_exp * P1z) + (269325.0 * P3x * dP4y_exp * P1z) + (113400.0 * P5x * dP2y_exp * P1z) - (14175.0 * P7x * dP0y_exp * P1z);
    double dterm_6_dz = (861840.0 * P1x * P2y * dP5z_exp) - (45360.0 * P1x * dP7z_exp) + (30240.0 * P3x * dP5z_exp) - (978075.0 * P1x * P4y * dP3z_exp) - (916650.0 * P3x * P2y * dP3z_exp) + (61425.0 * P5x * dP3z_exp) + (141750.0 * P1x * P6y * dP1z_exp) + (269325.0 * P3x * P4y * dP1z_exp) + (113400.0 * P5x * P2y * dP1z_exp) - (14175.0 * P7x * dP1z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

    // dmiu4 dx/dy/dz
    deriv[9] = temp * dterm_4_dx;
    deriv[10] = temp * dterm_4_dy;
    deriv[11] = temp * dterm_4_dz;

    // dmiu5 dx/dy/dz
    deriv[12] = temp * dterm_5_dx;
    deriv[13] = temp * dterm_5_dy;
    deriv[14] = temp * dterm_5_dz;

    // dmiu6 dx/dy/dz
    deriv[15] = temp * dterm_6_dx;
    deriv[16] = temp * dterm_6_dy;
    deriv[17] = temp * dterm_6_dz;

}

void calc_solid_MCSH_8_7(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double term_1 = (1119825.0 * P4x * P4y) - (438480.0 * P6x * P2y) - (141750.0 * P4x * P2y * P2z) + (15120.0 * P8x) + (15120.0 * P6x * P2z) - (14175.0 * P4x * P4z) - (438480.0 * P2x * P6y) - (141750.0 * P2x * P4y * P2z) + (283500.0 * P2x * P2y * P4z) - (13230.0 * P2x * P6z) + (15120.0 * P8y) + (15120.0 * P6y * P2z) - (14175.0 * P4y * P4z) - (13230.0 * P2y * P6z) + (945.0 * P8z);
    double term_2 = (1119825.0 * P4x * P4z) - (438480.0 * P6x * P2z) - (141750.0 * P4x * P2y * P2z) + (15120.0 * P8x) + (15120.0 * P6x * P2y) - (14175.0 * P4x * P4y) - (438480.0 * P2x * P6z) - (141750.0 * P2x * P2y * P4z) + (283500.0 * P2x * P4y * P2z) - (13230.0 * P2x * P6y) + (15120.0 * P8z) + (15120.0 * P2y * P6z) - (14175.0 * P4y * P4z) - (13230.0 * P6y * P2z) + (945.0 * P8y);
    double term_3 = (1119825.0 * P4y * P4z) - (438480.0 * P6y * P2z) - (141750.0 * P2x * P4y * P2z) + (15120.0 * P8y) + (15120.0 * P2x * P6y) - (14175.0 * P4x * P4y) - (438480.0 * P2y * P6z) - (141750.0 * P2x * P2y * P4z) + (283500.0 * P4x * P2y * P2z) - (13230.0 * P6x * P2y) + (15120.0 * P8z) + (15120.0 * P2x * P6z) - (14175.0 * P4x * P4z) - (13230.0 * P6x * P2z) + (945.0 * P8x);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP4x_exp = dP4_exp(P4x, C2, lambda, x0, gamma);
    double dP4y_exp = dP4_exp(P4y, C2, lambda, y0, gamma);
    double dP4z_exp = dP4_exp(P4z, C2, lambda, z0, gamma);

    double dP6x_exp = dP6_exp(P6x, C2, lambda, x0, gamma);
    double dP6y_exp = dP6_exp(P6y, C2, lambda, y0, gamma);
    double dP6z_exp = dP6_exp(P6z, C2, lambda, z0, gamma);

    double dP8x_exp = dP8_exp(P8x, C2, lambda, x0, gamma);
    double dP8y_exp = dP8_exp(P8y, C2, lambda, y0, gamma);
    double dP8z_exp = dP8_exp(P8z, C2, lambda, z0, gamma);

    double dterm_1_dx = (1119825.0 * dP4x_exp * P4y) - (438480.0 * dP6x_exp * P2y) - (141750.0 * dP4x_exp * P2y * P2z) + (15120.0 * dP8x_exp) + (15120.0 * dP6x_exp * P2z) - (14175.0 * dP4x_exp * P4z) - (438480.0 * dP2x_exp * P6y) - (141750.0 * dP2x_exp * P4y * P2z) + (283500.0 * dP2x_exp * P2y * P4z) - (13230.0 * dP2x_exp * P6z) + (15120.0 * dP0x_exp * P8y) + (15120.0 * dP0x_exp * P6y * P2z) - (14175.0 * dP0x_exp * P4y * P4z) - (13230.0 * dP0x_exp * P2y * P6z) + (945.0 * dP0x_exp * P8z);
    double dterm_1_dy = (1119825.0 * P4x * dP4y_exp) - (438480.0 * P6x * dP2y_exp) - (141750.0 * P4x * dP2y_exp * P2z) + (15120.0 * P8x * dP0y_exp) + (15120.0 * P6x * dP0y_exp * P2z) - (14175.0 * P4x * dP0y_exp * P4z) - (438480.0 * P2x * dP6y_exp) - (141750.0 * P2x * dP4y_exp * P2z) + (283500.0 * P2x * dP2y_exp * P4z) - (13230.0 * P2x * dP0y_exp * P6z) + (15120.0 * dP8y_exp) + (15120.0 * dP6y_exp * P2z) - (14175.0 * dP4y_exp * P4z) - (13230.0 * dP2y_exp * P6z) + (945.0 * dP0y_exp * P8z);
    double dterm_1_dz = (1119825.0 * P4x * P4y * dP0z_exp) - (438480.0 * P6x * P2y * dP0z_exp) - (141750.0 * P4x * P2y * dP2z_exp) + (15120.0 * P8x * dP0z_exp) + (15120.0 * P6x * dP2z_exp) - (14175.0 * P4x * dP4z_exp) - (438480.0 * P2x * P6y * dP0z_exp) - (141750.0 * P2x * P4y * dP2z_exp) + (283500.0 * P2x * P2y * dP4z_exp) - (13230.0 * P2x * dP6z_exp) + (15120.0 * P8y * dP0z_exp) + (15120.0 * P6y * dP2z_exp) - (14175.0 * P4y * dP4z_exp) - (13230.0 * P2y * dP6z_exp) + (945.0 * dP8z_exp);

    double dterm_2_dx = (1119825.0 * dP4x_exp * P4z) - (438480.0 * dP6x_exp * P2z) - (141750.0 * dP4x_exp * P2y * P2z) + (15120.0 * dP8x_exp) + (15120.0 * dP6x_exp * P2y) - (14175.0 * dP4x_exp * P4y) - (438480.0 * dP2x_exp * P6z) - (141750.0 * dP2x_exp * P2y * P4z) + (283500.0 * dP2x_exp * P4y * P2z) - (13230.0 * dP2x_exp * P6y) + (15120.0 * dP0x_exp * P8z) + (15120.0 * dP0x_exp * P2y * P6z) - (14175.0 * dP0x_exp * P4y * P4z) - (13230.0 * dP0x_exp * P6y * P2z) + (945.0 * dP0x_exp * P8y);
    double dterm_2_dy = (1119825.0 * P4x * dP0y_exp * P4z) - (438480.0 * P6x * dP0y_exp * P2z) - (141750.0 * P4x * dP2y_exp * P2z) + (15120.0 * P8x * dP0y_exp) + (15120.0 * P6x * dP2y_exp) - (14175.0 * P4x * dP4y_exp) - (438480.0 * P2x * dP0y_exp * P6z) - (141750.0 * P2x * dP2y_exp * P4z) + (283500.0 * P2x * dP4y_exp * P2z) - (13230.0 * P2x * dP6y_exp) + (15120.0 * dP0y_exp * P8z) + (15120.0 * dP2y_exp * P6z) - (14175.0 * dP4y_exp * P4z) - (13230.0 * dP6y_exp * P2z) + (945.0 * dP8y_exp);
    double dterm_2_dz = (1119825.0 * P4x * dP4z_exp) - (438480.0 * P6x * dP2z_exp) - (141750.0 * P4x * P2y * dP2z_exp) + (15120.0 * P8x * dP0z_exp) + (15120.0 * P6x * P2y * dP0z_exp) - (14175.0 * P4x * P4y * dP0z_exp) - (438480.0 * P2x * dP6z_exp) - (141750.0 * P2x * P2y * dP4z_exp) + (283500.0 * P2x * P4y * dP2z_exp) - (13230.0 * P2x * P6y * dP0z_exp) + (15120.0 * dP8z_exp) + (15120.0 * P2y * dP6z_exp) - (14175.0 * P4y * dP4z_exp) - (13230.0 * P6y * dP2z_exp) + (945.0 * P8y * dP0z_exp);

    double dterm_3_dx = (1119825.0 * dP0x_exp * P4y * P4z) - (438480.0 * dP0x_exp * P6y * P2z) - (141750.0 * dP2x_exp * P4y * P2z) + (15120.0 * dP0x_exp * P8y) + (15120.0 * dP2x_exp * P6y) - (14175.0 * dP4x_exp * P4y) - (438480.0 * dP0x_exp * P2y * P6z) - (141750.0 * dP2x_exp * P2y * P4z) + (283500.0 * dP4x_exp * P2y * P2z) - (13230.0 * dP6x_exp * P2y) + (15120.0 * dP0x_exp * P8z) + (15120.0 * dP2x_exp * P6z) - (14175.0 * dP4x_exp * P4z) - (13230.0 * dP6x_exp * P2z) + (945.0 * dP8x_exp);
    double dterm_3_dy = (1119825.0 * dP4y_exp * P4z) - (438480.0 * dP6y_exp * P2z) - (141750.0 * P2x * dP4y_exp * P2z) + (15120.0 * dP8y_exp) + (15120.0 * P2x * dP6y_exp) - (14175.0 * P4x * dP4y_exp) - (438480.0 * dP2y_exp * P6z) - (141750.0 * P2x * dP2y_exp * P4z) + (283500.0 * P4x * dP2y_exp * P2z) - (13230.0 * P6x * dP2y_exp) + (15120.0 * dP0y_exp * P8z) + (15120.0 * P2x * dP0y_exp * P6z) - (14175.0 * P4x * dP0y_exp * P4z) - (13230.0 * P6x * dP0y_exp * P2z) + (945.0 * P8x * dP0y_exp);
    double dterm_3_dz = (1119825.0 * P4y * dP4z_exp) - (438480.0 * P6y * dP2z_exp) - (141750.0 * P2x * P4y * dP2z_exp) + (15120.0 * P8y * dP0z_exp) + (15120.0 * P2x * P6y * dP0z_exp) - (14175.0 * P4x * P4y * dP0z_exp) - (438480.0 * P2y * dP6z_exp) - (141750.0 * P2x * P2y * dP4z_exp) + (283500.0 * P4x * P2y * dP2z_exp) - (13230.0 * P6x * P2y * dP0z_exp) + (15120.0 * dP8z_exp) + (15120.0 * P2x * dP6z_exp) - (14175.0 * P4x * dP4z_exp) - (13230.0 * P6x * dP2z_exp) + (945.0 * P8x * dP0z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

}

void calc_solid_MCSH_8_8(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double term_1 = (1190700.0 * P4x * P3y * P1z) - (226800.0 * P6x * P1y * P1z) - (56700.0 * P4x * P1y * P3z) - (586845.0 * P2x * P5y * P1z) - (425250.0 * P2x * P3y * P3z) + (161595.0 * P2x * P1y * P5z) + (22680.0 * P7y * P1z) + (36855.0 * P5y * P3z) + (5670.0 * P3y * P5z) - (8505.0 * P1y * P7z);
    double term_2 = (1190700.0 * P3x * P4y * P1z) - (226800.0 * P1x * P6y * P1z) - (56700.0 * P1x * P4y * P3z) - (586845.0 * P5x * P2y * P1z) - (425250.0 * P3x * P2y * P3z) + (161595.0 * P1x * P2y * P5z) + (22680.0 * P7x * P1z) + (36855.0 * P5x * P3z) + (5670.0 * P3x * P5z) - (8505.0 * P1x * P7z);
    double term_3 = (1190700.0 * P4x * P1y * P3z) - (226800.0 * P6x * P1y * P1z) - (56700.0 * P4x * P3y * P1z) - (586845.0 * P2x * P1y * P5z) - (425250.0 * P2x * P3y * P3z) + (161595.0 * P2x * P5y * P1z) + (22680.0 * P1y * P7z) + (36855.0 * P3y * P5z) + (5670.0 * P5y * P3z) - (8505.0 * P7y * P1z);
    double term_4 = (1190700.0 * P3x * P1y * P4z) - (226800.0 * P1x * P1y * P6z) - (56700.0 * P1x * P3y * P4z) - (586845.0 * P5x * P1y * P2z) - (425250.0 * P3x * P3y * P2z) + (161595.0 * P1x * P5y * P2z) + (22680.0 * P7x * P1y) + (36855.0 * P5x * P3y) + (5670.0 * P3x * P5y) - (8505.0 * P1x * P7y);
    double term_5 = (1190700.0 * P1x * P4y * P3z) - (226800.0 * P1x * P6y * P1z) - (56700.0 * P3x * P4y * P1z) - (586845.0 * P1x * P2y * P5z) - (425250.0 * P3x * P2y * P3z) + (161595.0 * P5x * P2y * P1z) + (22680.0 * P1x * P7z) + (36855.0 * P3x * P5z) + (5670.0 * P5x * P3z) - (8505.0 * P7x * P1z);
    double term_6 = (1190700.0 * P1x * P3y * P4z) - (226800.0 * P1x * P1y * P6z) - (56700.0 * P3x * P1y * P4z) - (586845.0 * P1x * P5y * P2z) - (425250.0 * P3x * P3y * P2z) + (161595.0 * P5x * P1y * P2z) + (22680.0 * P1x * P7y) + (36855.0 * P3x * P5y) + (5670.0 * P5x * P3y) - (8505.0 * P7x * P1y);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP3x_exp = dP3_exp(P3x, C2, lambda, x0, gamma);
    double dP3y_exp = dP3_exp(P3y, C2, lambda, y0, gamma);
    double dP3z_exp = dP3_exp(P3z, C2, lambda, z0, gamma);

    double dP4x_exp = dP4_exp(P4x, C2, lambda, x0, gamma);
    double dP4y_exp = dP4_exp(P4y, C2, lambda, y0, gamma);
    double dP4z_exp = dP4_exp(P4z, C2, lambda, z0, gamma);

    double dP5x_exp = dP5_exp(P5x, C2, lambda, x0, gamma);
    double dP5y_exp = dP5_exp(P5y, C2, lambda, y0, gamma);
    double dP5z_exp = dP5_exp(P5z, C2, lambda, z0, gamma);

    double dP6x_exp = dP6_exp(P6x, C2, lambda, x0, gamma);
    double dP6y_exp = dP6_exp(P6y, C2, lambda, y0, gamma);
    double dP6z_exp = dP6_exp(P6z, C2, lambda, z0, gamma);

    double dP7x_exp = dP7_exp(P7x, C2, lambda, x0, gamma);
    double dP7y_exp = dP7_exp(P7y, C2, lambda, y0, gamma);
    double dP7z_exp = dP7_exp(P7z, C2, lambda, z0, gamma);

    double dterm_1_dx = (1190700.0 * dP4x_exp * P3y * P1z) - (226800.0 * dP6x_exp * P1y * P1z) - (56700.0 * dP4x_exp * P1y * P3z) - (586845.0 * dP2x_exp * P5y * P1z) - (425250.0 * dP2x_exp * P3y * P3z) + (161595.0 * dP2x_exp * P1y * P5z) + (22680.0 * dP0x_exp * P7y * P1z) + (36855.0 * dP0x_exp * P5y * P3z) + (5670.0 * dP0x_exp * P3y * P5z) - (8505.0 * dP0x_exp * P1y * P7z);
    double dterm_1_dy = (1190700.0 * P4x * dP3y_exp * P1z) - (226800.0 * P6x * dP1y_exp * P1z) - (56700.0 * P4x * dP1y_exp * P3z) - (586845.0 * P2x * dP5y_exp * P1z) - (425250.0 * P2x * dP3y_exp * P3z) + (161595.0 * P2x * dP1y_exp * P5z) + (22680.0 * dP7y_exp * P1z) + (36855.0 * dP5y_exp * P3z) + (5670.0 * dP3y_exp * P5z) - (8505.0 * dP1y_exp * P7z);
    double dterm_1_dz = (1190700.0 * P4x * P3y * dP1z_exp) - (226800.0 * P6x * P1y * dP1z_exp) - (56700.0 * P4x * P1y * dP3z_exp) - (586845.0 * P2x * P5y * dP1z_exp) - (425250.0 * P2x * P3y * dP3z_exp) + (161595.0 * P2x * P1y * dP5z_exp) + (22680.0 * P7y * dP1z_exp) + (36855.0 * P5y * dP3z_exp) + (5670.0 * P3y * dP5z_exp) - (8505.0 * P1y * dP7z_exp);

    double dterm_2_dx = (1190700.0 * dP3x_exp * P4y * P1z) - (226800.0 * dP1x_exp * P6y * P1z) - (56700.0 * dP1x_exp * P4y * P3z) - (586845.0 * dP5x_exp * P2y * P1z) - (425250.0 * dP3x_exp * P2y * P3z) + (161595.0 * dP1x_exp * P2y * P5z) + (22680.0 * dP7x_exp * P1z) + (36855.0 * dP5x_exp * P3z) + (5670.0 * dP3x_exp * P5z) - (8505.0 * dP1x_exp * P7z);
    double dterm_2_dy = (1190700.0 * P3x * dP4y_exp * P1z) - (226800.0 * P1x * dP6y_exp * P1z) - (56700.0 * P1x * dP4y_exp * P3z) - (586845.0 * P5x * dP2y_exp * P1z) - (425250.0 * P3x * dP2y_exp * P3z) + (161595.0 * P1x * dP2y_exp * P5z) + (22680.0 * P7x * dP0y_exp * P1z) + (36855.0 * P5x * dP0y_exp * P3z) + (5670.0 * P3x * dP0y_exp * P5z) - (8505.0 * P1x * dP0y_exp * P7z);
    double dterm_2_dz = (1190700.0 * P3x * P4y * dP1z_exp) - (226800.0 * P1x * P6y * dP1z_exp) - (56700.0 * P1x * P4y * dP3z_exp) - (586845.0 * P5x * P2y * dP1z_exp) - (425250.0 * P3x * P2y * dP3z_exp) + (161595.0 * P1x * P2y * dP5z_exp) + (22680.0 * P7x * dP1z_exp) + (36855.0 * P5x * dP3z_exp) + (5670.0 * P3x * dP5z_exp) - (8505.0 * P1x * dP7z_exp);

    double dterm_3_dx = (1190700.0 * dP4x_exp * P1y * P3z) - (226800.0 * dP6x_exp * P1y * P1z) - (56700.0 * dP4x_exp * P3y * P1z) - (586845.0 * dP2x_exp * P1y * P5z) - (425250.0 * dP2x_exp * P3y * P3z) + (161595.0 * dP2x_exp * P5y * P1z) + (22680.0 * dP0x_exp * P1y * P7z) + (36855.0 * dP0x_exp * P3y * P5z) + (5670.0 * dP0x_exp * P5y * P3z) - (8505.0 * dP0x_exp * P7y * P1z);
    double dterm_3_dy = (1190700.0 * P4x * dP1y_exp * P3z) - (226800.0 * P6x * dP1y_exp * P1z) - (56700.0 * P4x * dP3y_exp * P1z) - (586845.0 * P2x * dP1y_exp * P5z) - (425250.0 * P2x * dP3y_exp * P3z) + (161595.0 * P2x * dP5y_exp * P1z) + (22680.0 * dP1y_exp * P7z) + (36855.0 * dP3y_exp * P5z) + (5670.0 * dP5y_exp * P3z) - (8505.0 * dP7y_exp * P1z);
    double dterm_3_dz = (1190700.0 * P4x * P1y * dP3z_exp) - (226800.0 * P6x * P1y * dP1z_exp) - (56700.0 * P4x * P3y * dP1z_exp) - (586845.0 * P2x * P1y * dP5z_exp) - (425250.0 * P2x * P3y * dP3z_exp) + (161595.0 * P2x * P5y * dP1z_exp) + (22680.0 * P1y * dP7z_exp) + (36855.0 * P3y * dP5z_exp) + (5670.0 * P5y * dP3z_exp) - (8505.0 * P7y * dP1z_exp);

    double dterm_4_dx = (1190700.0 * dP3x_exp * P1y * P4z) - (226800.0 * dP1x_exp * P1y * P6z) - (56700.0 * dP1x_exp * P3y * P4z) - (586845.0 * dP5x_exp * P1y * P2z) - (425250.0 * dP3x_exp * P3y * P2z) + (161595.0 * dP1x_exp * P5y * P2z) + (22680.0 * dP7x_exp * P1y) + (36855.0 * dP5x_exp * P3y) + (5670.0 * dP3x_exp * P5y) - (8505.0 * dP1x_exp * P7y);
    double dterm_4_dy = (1190700.0 * P3x * dP1y_exp * P4z) - (226800.0 * P1x * dP1y_exp * P6z) - (56700.0 * P1x * dP3y_exp * P4z) - (586845.0 * P5x * dP1y_exp * P2z) - (425250.0 * P3x * dP3y_exp * P2z) + (161595.0 * P1x * dP5y_exp * P2z) + (22680.0 * P7x * dP1y_exp) + (36855.0 * P5x * dP3y_exp) + (5670.0 * P3x * dP5y_exp) - (8505.0 * P1x * dP7y_exp);
    double dterm_4_dz = (1190700.0 * P3x * P1y * dP4z_exp) - (226800.0 * P1x * P1y * dP6z_exp) - (56700.0 * P1x * P3y * dP4z_exp) - (586845.0 * P5x * P1y * dP2z_exp) - (425250.0 * P3x * P3y * dP2z_exp) + (161595.0 * P1x * P5y * dP2z_exp) + (22680.0 * P7x * P1y * dP0z_exp) + (36855.0 * P5x * P3y * dP0z_exp) + (5670.0 * P3x * P5y * dP0z_exp) - (8505.0 * P1x * P7y * dP0z_exp);

    double dterm_5_dx = (1190700.0 * dP1x_exp * P4y * P3z) - (226800.0 * dP1x_exp * P6y * P1z) - (56700.0 * dP3x_exp * P4y * P1z) - (586845.0 * dP1x_exp * P2y * P5z) - (425250.0 * dP3x_exp * P2y * P3z) + (161595.0 * dP5x_exp * P2y * P1z) + (22680.0 * dP1x_exp * P7z) + (36855.0 * dP3x_exp * P5z) + (5670.0 * dP5x_exp * P3z) - (8505.0 * dP7x_exp * P1z);
    double dterm_5_dy = (1190700.0 * P1x * dP4y_exp * P3z) - (226800.0 * P1x * dP6y_exp * P1z) - (56700.0 * P3x * dP4y_exp * P1z) - (586845.0 * P1x * dP2y_exp * P5z) - (425250.0 * P3x * dP2y_exp * P3z) + (161595.0 * P5x * dP2y_exp * P1z) + (22680.0 * P1x * dP0y_exp * P7z) + (36855.0 * P3x * dP0y_exp * P5z) + (5670.0 * P5x * dP0y_exp * P3z) - (8505.0 * P7x * dP0y_exp * P1z);
    double dterm_5_dz = (1190700.0 * P1x * P4y * dP3z_exp) - (226800.0 * P1x * P6y * dP1z_exp) - (56700.0 * P3x * P4y * dP1z_exp) - (586845.0 * P1x * P2y * dP5z_exp) - (425250.0 * P3x * P2y * dP3z_exp) + (161595.0 * P5x * P2y * dP1z_exp) + (22680.0 * P1x * dP7z_exp) + (36855.0 * P3x * dP5z_exp) + (5670.0 * P5x * dP3z_exp) - (8505.0 * P7x * dP1z_exp);

    double dterm_6_dx = (1190700.0 * dP1x_exp * P3y * P4z) - (226800.0 * dP1x_exp * P1y * P6z) - (56700.0 * dP3x_exp * P1y * P4z) - (586845.0 * dP1x_exp * P5y * P2z) - (425250.0 * dP3x_exp * P3y * P2z) + (161595.0 * dP5x_exp * P1y * P2z) + (22680.0 * dP1x_exp * P7y) + (36855.0 * dP3x_exp * P5y) + (5670.0 * dP5x_exp * P3y) - (8505.0 * dP7x_exp * P1y);
    double dterm_6_dy = (1190700.0 * P1x * dP3y_exp * P4z) - (226800.0 * P1x * dP1y_exp * P6z) - (56700.0 * P3x * dP1y_exp * P4z) - (586845.0 * P1x * dP5y_exp * P2z) - (425250.0 * P3x * dP3y_exp * P2z) + (161595.0 * P5x * dP1y_exp * P2z) + (22680.0 * P1x * dP7y_exp) + (36855.0 * P3x * dP5y_exp) + (5670.0 * P5x * dP3y_exp) - (8505.0 * P7x * dP1y_exp);
    double dterm_6_dz = (1190700.0 * P1x * P3y * dP4z_exp) - (226800.0 * P1x * P1y * dP6z_exp) - (56700.0 * P3x * P1y * dP4z_exp) - (586845.0 * P1x * P5y * dP2z_exp) - (425250.0 * P3x * P3y * dP2z_exp) + (161595.0 * P5x * P1y * dP2z_exp) + (22680.0 * P1x * P7y * dP0z_exp) + (36855.0 * P3x * P5y * dP0z_exp) + (5670.0 * P5x * P3y * dP0z_exp) - (8505.0 * P7x * P1y * dP0z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

    // dmiu4 dx/dy/dz
    deriv[9] = temp * dterm_4_dx;
    deriv[10] = temp * dterm_4_dy;
    deriv[11] = temp * dterm_4_dz;

    // dmiu5 dx/dy/dz
    deriv[12] = temp * dterm_5_dx;
    deriv[13] = temp * dterm_5_dy;
    deriv[14] = temp * dterm_5_dz;

    // dmiu6 dx/dy/dz
    deriv[15] = temp * dterm_6_dx;
    deriv[16] = temp * dterm_6_dy;
    deriv[17] = temp * dterm_6_dz;

}

void calc_solid_MCSH_8_9(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double term_1 = (1200150.0 * P4x * P2y * P2z) - (70560.0 * P6x * P2y) - (23625.0 * P4x * P4y) - (70560.0 * P6x * P2z) - (23625.0 * P4x * P4z) + (5040.0 * P8x) - (600075.0 * P2x * P4y * P2z) - (600075.0 * P2x * P2y * P4z) + (49455.0 * P2x * P6y) + (49455.0 * P2x * P6z) + (21105.0 * P6y * P2z) + (47250.0 * P4y * P4z) + (21105.0 * P2y * P6z) - (2520.0 * P8y) - (2520.0 * P8z);
    double term_2 = (1200150.0 * P2x * P4y * P2z) - (70560.0 * P2x * P6y) - (23625.0 * P4x * P4y) - (70560.0 * P6y * P2z) - (23625.0 * P4y * P4z) + (5040.0 * P8y) - (600075.0 * P4x * P2y * P2z) - (600075.0 * P2x * P2y * P4z) + (49455.0 * P6x * P2y) + (49455.0 * P2y * P6z) + (21105.0 * P6x * P2z) + (47250.0 * P4x * P4z) + (21105.0 * P2x * P6z) - (2520.0 * P8x) - (2520.0 * P8z);
    double term_3 = (1200150.0 * P2x * P2y * P4z) - (70560.0 * P2x * P6z) - (23625.0 * P4x * P4z) - (70560.0 * P2y * P6z) - (23625.0 * P4y * P4z) + (5040.0 * P8z) - (600075.0 * P4x * P2y * P2z) - (600075.0 * P2x * P4y * P2z) + (49455.0 * P6x * P2z) + (49455.0 * P6y * P2z) + (21105.0 * P6x * P2y) + (47250.0 * P4x * P4y) + (21105.0 * P2x * P6y) - (2520.0 * P8x) - (2520.0 * P8y);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP4x_exp = dP4_exp(P4x, C2, lambda, x0, gamma);
    double dP4y_exp = dP4_exp(P4y, C2, lambda, y0, gamma);
    double dP4z_exp = dP4_exp(P4z, C2, lambda, z0, gamma);

    double dP6x_exp = dP6_exp(P6x, C2, lambda, x0, gamma);
    double dP6y_exp = dP6_exp(P6y, C2, lambda, y0, gamma);
    double dP6z_exp = dP6_exp(P6z, C2, lambda, z0, gamma);

    double dP8x_exp = dP8_exp(P8x, C2, lambda, x0, gamma);
    double dP8y_exp = dP8_exp(P8y, C2, lambda, y0, gamma);
    double dP8z_exp = dP8_exp(P8z, C2, lambda, z0, gamma);

    double dterm_1_dx = (1200150.0 * dP4x_exp * P2y * P2z) - (70560.0 * dP6x_exp * P2y) - (23625.0 * dP4x_exp * P4y) - (70560.0 * dP6x_exp * P2z) - (23625.0 * dP4x_exp * P4z) + (5040.0 * dP8x_exp) - (600075.0 * dP2x_exp * P4y * P2z) - (600075.0 * dP2x_exp * P2y * P4z) + (49455.0 * dP2x_exp * P6y) + (49455.0 * dP2x_exp * P6z) + (21105.0 * dP0x_exp * P6y * P2z) + (47250.0 * dP0x_exp * P4y * P4z) + (21105.0 * dP0x_exp * P2y * P6z) - (2520.0 * dP0x_exp * P8y) - (2520.0 * dP0x_exp * P8z);
    double dterm_1_dy = (1200150.0 * P4x * dP2y_exp * P2z) - (70560.0 * P6x * dP2y_exp) - (23625.0 * P4x * dP4y_exp) - (70560.0 * P6x * dP0y_exp * P2z) - (23625.0 * P4x * dP0y_exp * P4z) + (5040.0 * P8x * dP0y_exp) - (600075.0 * P2x * dP4y_exp * P2z) - (600075.0 * P2x * dP2y_exp * P4z) + (49455.0 * P2x * dP6y_exp) + (49455.0 * P2x * dP0y_exp * P6z) + (21105.0 * dP6y_exp * P2z) + (47250.0 * dP4y_exp * P4z) + (21105.0 * dP2y_exp * P6z) - (2520.0 * dP8y_exp) - (2520.0 * dP0y_exp * P8z);
    double dterm_1_dz = (1200150.0 * P4x * P2y * dP2z_exp) - (70560.0 * P6x * P2y * dP0z_exp) - (23625.0 * P4x * P4y * dP0z_exp) - (70560.0 * P6x * dP2z_exp) - (23625.0 * P4x * dP4z_exp) + (5040.0 * P8x * dP0z_exp) - (600075.0 * P2x * P4y * dP2z_exp) - (600075.0 * P2x * P2y * dP4z_exp) + (49455.0 * P2x * P6y * dP0z_exp) + (49455.0 * P2x * dP6z_exp) + (21105.0 * P6y * dP2z_exp) + (47250.0 * P4y * dP4z_exp) + (21105.0 * P2y * dP6z_exp) - (2520.0 * P8y * dP0z_exp) - (2520.0 * dP8z_exp);

    double dterm_2_dx = (1200150.0 * dP2x_exp * P4y * P2z) - (70560.0 * dP2x_exp * P6y) - (23625.0 * dP4x_exp * P4y) - (70560.0 * dP0x_exp * P6y * P2z) - (23625.0 * dP0x_exp * P4y * P4z) + (5040.0 * dP0x_exp * P8y) - (600075.0 * dP4x_exp * P2y * P2z) - (600075.0 * dP2x_exp * P2y * P4z) + (49455.0 * dP6x_exp * P2y) + (49455.0 * dP0x_exp * P2y * P6z) + (21105.0 * dP6x_exp * P2z) + (47250.0 * dP4x_exp * P4z) + (21105.0 * dP2x_exp * P6z) - (2520.0 * dP8x_exp) - (2520.0 * dP0x_exp * P8z);
    double dterm_2_dy = (1200150.0 * P2x * dP4y_exp * P2z) - (70560.0 * P2x * dP6y_exp) - (23625.0 * P4x * dP4y_exp) - (70560.0 * dP6y_exp * P2z) - (23625.0 * dP4y_exp * P4z) + (5040.0 * dP8y_exp) - (600075.0 * P4x * dP2y_exp * P2z) - (600075.0 * P2x * dP2y_exp * P4z) + (49455.0 * P6x * dP2y_exp) + (49455.0 * dP2y_exp * P6z) + (21105.0 * P6x * dP0y_exp * P2z) + (47250.0 * P4x * dP0y_exp * P4z) + (21105.0 * P2x * dP0y_exp * P6z) - (2520.0 * P8x * dP0y_exp) - (2520.0 * dP0y_exp * P8z);
    double dterm_2_dz = (1200150.0 * P2x * P4y * dP2z_exp) - (70560.0 * P2x * P6y * dP0z_exp) - (23625.0 * P4x * P4y * dP0z_exp) - (70560.0 * P6y * dP2z_exp) - (23625.0 * P4y * dP4z_exp) + (5040.0 * P8y * dP0z_exp) - (600075.0 * P4x * P2y * dP2z_exp) - (600075.0 * P2x * P2y * dP4z_exp) + (49455.0 * P6x * P2y * dP0z_exp) + (49455.0 * P2y * dP6z_exp) + (21105.0 * P6x * dP2z_exp) + (47250.0 * P4x * dP4z_exp) + (21105.0 * P2x * dP6z_exp) - (2520.0 * P8x * dP0z_exp) - (2520.0 * dP8z_exp);

    double dterm_3_dx = (1200150.0 * dP2x_exp * P2y * P4z) - (70560.0 * dP2x_exp * P6z) - (23625.0 * dP4x_exp * P4z) - (70560.0 * dP0x_exp * P2y * P6z) - (23625.0 * dP0x_exp * P4y * P4z) + (5040.0 * dP0x_exp * P8z) - (600075.0 * dP4x_exp * P2y * P2z) - (600075.0 * dP2x_exp * P4y * P2z) + (49455.0 * dP6x_exp * P2z) + (49455.0 * dP0x_exp * P6y * P2z) + (21105.0 * dP6x_exp * P2y) + (47250.0 * dP4x_exp * P4y) + (21105.0 * dP2x_exp * P6y) - (2520.0 * dP8x_exp) - (2520.0 * dP0x_exp * P8y);
    double dterm_3_dy = (1200150.0 * P2x * dP2y_exp * P4z) - (70560.0 * P2x * dP0y_exp * P6z) - (23625.0 * P4x * dP0y_exp * P4z) - (70560.0 * dP2y_exp * P6z) - (23625.0 * dP4y_exp * P4z) + (5040.0 * dP0y_exp * P8z) - (600075.0 * P4x * dP2y_exp * P2z) - (600075.0 * P2x * dP4y_exp * P2z) + (49455.0 * P6x * dP0y_exp * P2z) + (49455.0 * dP6y_exp * P2z) + (21105.0 * P6x * dP2y_exp) + (47250.0 * P4x * dP4y_exp) + (21105.0 * P2x * dP6y_exp) - (2520.0 * P8x * dP0y_exp) - (2520.0 * dP8y_exp);
    double dterm_3_dz = (1200150.0 * P2x * P2y * dP4z_exp) - (70560.0 * P2x * dP6z_exp) - (23625.0 * P4x * dP4z_exp) - (70560.0 * P2y * dP6z_exp) - (23625.0 * P4y * dP4z_exp) + (5040.0 * dP8z_exp) - (600075.0 * P4x * P2y * dP2z_exp) - (600075.0 * P2x * P4y * dP2z_exp) + (49455.0 * P6x * dP2z_exp) + (49455.0 * P6y * dP2z_exp) + (21105.0 * P6x * P2y * dP0z_exp) + (47250.0 * P4x * P4y * dP0z_exp) + (21105.0 * P2x * P6y * dP0z_exp) - (2520.0 * P8x * dP0z_exp) - (2520.0 * P8y * dP0z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

}

void calc_solid_MCSH_8_10(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double term_1 = (1341900.0 * P3x * P3y * P2z) - (67095.0 * P5x * P3y) - (67095.0 * P3x * P5y) - (274995.0 * P5x * P1y * P2z) - (212625.0 * P3x * P1y * P4z) + (22680.0 * P7x * P1y) - (274995.0 * P1x * P5y * P2z) - (212625.0 * P1x * P3y * P4z) + (22680.0 * P1x * P7y) + (85050.0 * P1x * P1y * P6z);
    double term_2 = (1341900.0 * P3x * P2y * P3z) - (67095.0 * P5x * P3z) - (67095.0 * P3x * P5z) - (274995.0 * P5x * P2y * P1z) - (212625.0 * P3x * P4y * P1z) + (22680.0 * P7x * P1z) - (274995.0 * P1x * P2y * P5z) - (212625.0 * P1x * P4y * P3z) + (22680.0 * P1x * P7z) + (85050.0 * P1x * P6y * P1z);
    double term_3 = (1341900.0 * P2x * P3y * P3z) - (67095.0 * P5y * P3z) - (67095.0 * P3y * P5z) - (274995.0 * P2x * P5y * P1z) - (212625.0 * P4x * P3y * P1z) + (22680.0 * P7y * P1z) - (274995.0 * P2x * P1y * P5z) - (212625.0 * P4x * P1y * P3z) + (22680.0 * P1y * P7z) + (85050.0 * P6x * P1y * P1z);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP3x_exp = dP3_exp(P3x, C2, lambda, x0, gamma);
    double dP3y_exp = dP3_exp(P3y, C2, lambda, y0, gamma);
    double dP3z_exp = dP3_exp(P3z, C2, lambda, z0, gamma);

    double dP4x_exp = dP4_exp(P4x, C2, lambda, x0, gamma);
    double dP4y_exp = dP4_exp(P4y, C2, lambda, y0, gamma);
    double dP4z_exp = dP4_exp(P4z, C2, lambda, z0, gamma);

    double dP5x_exp = dP5_exp(P5x, C2, lambda, x0, gamma);
    double dP5y_exp = dP5_exp(P5y, C2, lambda, y0, gamma);
    double dP5z_exp = dP5_exp(P5z, C2, lambda, z0, gamma);

    double dP6x_exp = dP6_exp(P6x, C2, lambda, x0, gamma);
    double dP6y_exp = dP6_exp(P6y, C2, lambda, y0, gamma);
    double dP6z_exp = dP6_exp(P6z, C2, lambda, z0, gamma);

    double dP7x_exp = dP7_exp(P7x, C2, lambda, x0, gamma);
    double dP7y_exp = dP7_exp(P7y, C2, lambda, y0, gamma);
    double dP7z_exp = dP7_exp(P7z, C2, lambda, z0, gamma);

    double dterm_1_dx = (1341900.0 * dP3x_exp * P3y * P2z) - (67095.0 * dP5x_exp * P3y) - (67095.0 * dP3x_exp * P5y) - (274995.0 * dP5x_exp * P1y * P2z) - (212625.0 * dP3x_exp * P1y * P4z) + (22680.0 * dP7x_exp * P1y) - (274995.0 * dP1x_exp * P5y * P2z) - (212625.0 * dP1x_exp * P3y * P4z) + (22680.0 * dP1x_exp * P7y) + (85050.0 * dP1x_exp * P1y * P6z);
    double dterm_1_dy = (1341900.0 * P3x * dP3y_exp * P2z) - (67095.0 * P5x * dP3y_exp) - (67095.0 * P3x * dP5y_exp) - (274995.0 * P5x * dP1y_exp * P2z) - (212625.0 * P3x * dP1y_exp * P4z) + (22680.0 * P7x * dP1y_exp) - (274995.0 * P1x * dP5y_exp * P2z) - (212625.0 * P1x * dP3y_exp * P4z) + (22680.0 * P1x * dP7y_exp) + (85050.0 * P1x * dP1y_exp * P6z);
    double dterm_1_dz = (1341900.0 * P3x * P3y * dP2z_exp) - (67095.0 * P5x * P3y * dP0z_exp) - (67095.0 * P3x * P5y * dP0z_exp) - (274995.0 * P5x * P1y * dP2z_exp) - (212625.0 * P3x * P1y * dP4z_exp) + (22680.0 * P7x * P1y * dP0z_exp) - (274995.0 * P1x * P5y * dP2z_exp) - (212625.0 * P1x * P3y * dP4z_exp) + (22680.0 * P1x * P7y * dP0z_exp) + (85050.0 * P1x * P1y * dP6z_exp);

    double dterm_2_dx = (1341900.0 * dP3x_exp * P2y * P3z) - (67095.0 * dP5x_exp * P3z) - (67095.0 * dP3x_exp * P5z) - (274995.0 * dP5x_exp * P2y * P1z) - (212625.0 * dP3x_exp * P4y * P1z) + (22680.0 * dP7x_exp * P1z) - (274995.0 * dP1x_exp * P2y * P5z) - (212625.0 * dP1x_exp * P4y * P3z) + (22680.0 * dP1x_exp * P7z) + (85050.0 * dP1x_exp * P6y * P1z);
    double dterm_2_dy = (1341900.0 * P3x * dP2y_exp * P3z) - (67095.0 * P5x * dP0y_exp * P3z) - (67095.0 * P3x * dP0y_exp * P5z) - (274995.0 * P5x * dP2y_exp * P1z) - (212625.0 * P3x * dP4y_exp * P1z) + (22680.0 * P7x * dP0y_exp * P1z) - (274995.0 * P1x * dP2y_exp * P5z) - (212625.0 * P1x * dP4y_exp * P3z) + (22680.0 * P1x * dP0y_exp * P7z) + (85050.0 * P1x * dP6y_exp * P1z);
    double dterm_2_dz = (1341900.0 * P3x * P2y * dP3z_exp) - (67095.0 * P5x * dP3z_exp) - (67095.0 * P3x * dP5z_exp) - (274995.0 * P5x * P2y * dP1z_exp) - (212625.0 * P3x * P4y * dP1z_exp) + (22680.0 * P7x * dP1z_exp) - (274995.0 * P1x * P2y * dP5z_exp) - (212625.0 * P1x * P4y * dP3z_exp) + (22680.0 * P1x * dP7z_exp) + (85050.0 * P1x * P6y * dP1z_exp);

    double dterm_3_dx = (1341900.0 * dP2x_exp * P3y * P3z) - (67095.0 * dP0x_exp * P5y * P3z) - (67095.0 * dP0x_exp * P3y * P5z) - (274995.0 * dP2x_exp * P5y * P1z) - (212625.0 * dP4x_exp * P3y * P1z) + (22680.0 * dP0x_exp * P7y * P1z) - (274995.0 * dP2x_exp * P1y * P5z) - (212625.0 * dP4x_exp * P1y * P3z) + (22680.0 * dP0x_exp * P1y * P7z) + (85050.0 * dP6x_exp * P1y * P1z);
    double dterm_3_dy = (1341900.0 * P2x * dP3y_exp * P3z) - (67095.0 * dP5y_exp * P3z) - (67095.0 * dP3y_exp * P5z) - (274995.0 * P2x * dP5y_exp * P1z) - (212625.0 * P4x * dP3y_exp * P1z) + (22680.0 * dP7y_exp * P1z) - (274995.0 * P2x * dP1y_exp * P5z) - (212625.0 * P4x * dP1y_exp * P3z) + (22680.0 * dP1y_exp * P7z) + (85050.0 * P6x * dP1y_exp * P1z);
    double dterm_3_dz = (1341900.0 * P2x * P3y * dP3z_exp) - (67095.0 * P5y * dP3z_exp) - (67095.0 * P3y * dP5z_exp) - (274995.0 * P2x * P5y * dP1z_exp) - (212625.0 * P4x * P3y * dP1z_exp) + (22680.0 * P7y * dP1z_exp) - (274995.0 * P2x * P1y * dP5z_exp) - (212625.0 * P4x * P1y * dP3z_exp) + (22680.0 * P1y * dP7z_exp) + (85050.0 * P6x * P1y * dP1z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

}


void calc_solid_MCSH_9_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double P8x = P8(lambda, x0, gamma);
    double P8y = P8(lambda, y0, gamma);
    double P8z = P8(lambda, z0, gamma);

    double P9x = P9(lambda, x0, gamma);
    double P9y = P9(lambda, y0, gamma);
    double P9z = P9(lambda, z0, gamma);

    double term_1 = (362880.0 * P9x) - (6531840.0 * P7x * P2y) - (6531840.0 * P7x * P2z) + (17146080.0 * P5x * P4y) + (34292160.0 * P5x * P2y * P2z) + (17146080.0 * P5x * P4z) - (9525600.0 * P3x * P6y) - (28576800.0 * P3x * P4y * P2z) - (28576800.0 * P3x * P2y * P4z) - (9525600.0 * P3x * P6z) + (893025.0 * P1x * P8y) + (3572100.0 * P1x * P6y * P2z) + (5358150.0 * P1x * P4y * P4z) + (3572100.0 * P1x * P2y * P6z) + (893025.0 * P1x * P8z);
    double term_2 = (362880.0 * P9y) - (6531840.0 * P2x * P7y) - (6531840.0 * P7y * P2z) + (17146080.0 * P4x * P5y) + (34292160.0 * P2x * P5y * P2z) + (17146080.0 * P5y * P4z) - (9525600.0 * P6x * P3y) - (28576800.0 * P4x * P3y * P2z) - (28576800.0 * P2x * P3y * P4z) - (9525600.0 * P3y * P6z) + (893025.0 * P8x * P1y) + (3572100.0 * P6x * P1y * P2z) + (5358150.0 * P4x * P1y * P4z) + (3572100.0 * P2x * P1y * P6z) + (893025.0 * P1y * P8z);
    double term_3 = (362880.0 * P9z) - (6531840.0 * P2x * P7z) - (6531840.0 * P2y * P7z) + (17146080.0 * P4x * P5z) + (34292160.0 * P2x * P2y * P5z) + (17146080.0 * P4y * P5z) - (9525600.0 * P6x * P3z) - (28576800.0 * P4x * P2y * P3z) - (28576800.0 * P2x * P4y * P3z) - (9525600.0 * P6y * P3z) + (893025.0 * P8x * P1z) + (3572100.0 * P6x * P2y * P1z) + (5358150.0 * P4x * P4y * P1z) + (3572100.0 * P2x * P6y * P1z) + (893025.0 * P8y * P1z);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP3x_exp = dP3_exp(P3x, C2, lambda, x0, gamma);
    double dP3y_exp = dP3_exp(P3y, C2, lambda, y0, gamma);
    double dP3z_exp = dP3_exp(P3z, C2, lambda, z0, gamma);

    double dP4x_exp = dP4_exp(P4x, C2, lambda, x0, gamma);
    double dP4y_exp = dP4_exp(P4y, C2, lambda, y0, gamma);
    double dP4z_exp = dP4_exp(P4z, C2, lambda, z0, gamma);

    double dP5x_exp = dP5_exp(P5x, C2, lambda, x0, gamma);
    double dP5y_exp = dP5_exp(P5y, C2, lambda, y0, gamma);
    double dP5z_exp = dP5_exp(P5z, C2, lambda, z0, gamma);

    double dP6x_exp = dP6_exp(P6x, C2, lambda, x0, gamma);
    double dP6y_exp = dP6_exp(P6y, C2, lambda, y0, gamma);
    double dP6z_exp = dP6_exp(P6z, C2, lambda, z0, gamma);

    double dP7x_exp = dP7_exp(P7x, C2, lambda, x0, gamma);
    double dP7y_exp = dP7_exp(P7y, C2, lambda, y0, gamma);
    double dP7z_exp = dP7_exp(P7z, C2, lambda, z0, gamma);

    double dP8x_exp = dP8_exp(P8x, C2, lambda, x0, gamma);
    double dP8y_exp = dP8_exp(P8y, C2, lambda, y0, gamma);
    double dP8z_exp = dP8_exp(P8z, C2, lambda, z0, gamma);

    double dP9x_exp = dP9_exp(P9x, C2, lambda, x0, gamma);
    double dP9y_exp = dP9_exp(P9y, C2, lambda, y0, gamma);
    double dP9z_exp = dP9_exp(P9z, C2, lambda, z0, gamma);

    double dterm_1_dx = (362880.0 * dP9x_exp) - (6531840.0 * dP7x_exp * P2y) - (6531840.0 * dP7x_exp * P2z) + (17146080.0 * dP5x_exp * P4y) + (34292160.0 * dP5x_exp * P2y * P2z) + (17146080.0 * dP5x_exp * P4z) - (9525600.0 * dP3x_exp * P6y) - (28576800.0 * dP3x_exp * P4y * P2z) - (28576800.0 * dP3x_exp * P2y * P4z) - (9525600.0 * dP3x_exp * P6z) + (893025.0 * dP1x_exp * P8y) + (3572100.0 * dP1x_exp * P6y * P2z) + (5358150.0 * dP1x_exp * P4y * P4z) + (3572100.0 * dP1x_exp * P2y * P6z) + (893025.0 * dP1x_exp * P8z);
    double dterm_1_dy = (362880.0 * P9x * dP0y_exp) - (6531840.0 * P7x * dP2y_exp) - (6531840.0 * P7x * dP0y_exp * P2z) + (17146080.0 * P5x * dP4y_exp) + (34292160.0 * P5x * dP2y_exp * P2z) + (17146080.0 * P5x * dP0y_exp * P4z) - (9525600.0 * P3x * dP6y_exp) - (28576800.0 * P3x * dP4y_exp * P2z) - (28576800.0 * P3x * dP2y_exp * P4z) - (9525600.0 * P3x * dP0y_exp * P6z) + (893025.0 * P1x * dP8y_exp) + (3572100.0 * P1x * dP6y_exp * P2z) + (5358150.0 * P1x * dP4y_exp * P4z) + (3572100.0 * P1x * dP2y_exp * P6z) + (893025.0 * P1x * dP0y_exp * P8z);
    double dterm_1_dz = (362880.0 * P9x * dP0z_exp) - (6531840.0 * P7x * P2y * dP0z_exp) - (6531840.0 * P7x * dP2z_exp) + (17146080.0 * P5x * P4y * dP0z_exp) + (34292160.0 * P5x * P2y * dP2z_exp) + (17146080.0 * P5x * dP4z_exp) - (9525600.0 * P3x * P6y * dP0z_exp) - (28576800.0 * P3x * P4y * dP2z_exp) - (28576800.0 * P3x * P2y * dP4z_exp) - (9525600.0 * P3x * dP6z_exp) + (893025.0 * P1x * P8y * dP0z_exp) + (3572100.0 * P1x * P6y * dP2z_exp) + (5358150.0 * P1x * P4y * dP4z_exp) + (3572100.0 * P1x * P2y * dP6z_exp) + (893025.0 * P1x * dP8z_exp);

    double dterm_2_dx = (362880.0 * dP0x_exp * P9y) - (6531840.0 * dP2x_exp * P7y) - (6531840.0 * dP0x_exp * P7y * P2z) + (17146080.0 * dP4x_exp * P5y) + (34292160.0 * dP2x_exp * P5y * P2z) + (17146080.0 * dP0x_exp * P5y * P4z) - (9525600.0 * dP6x_exp * P3y) - (28576800.0 * dP4x_exp * P3y * P2z) - (28576800.0 * dP2x_exp * P3y * P4z) - (9525600.0 * dP0x_exp * P3y * P6z) + (893025.0 * dP8x_exp * P1y) + (3572100.0 * dP6x_exp * P1y * P2z) + (5358150.0 * dP4x_exp * P1y * P4z) + (3572100.0 * dP2x_exp * P1y * P6z) + (893025.0 * dP0x_exp * P1y * P8z);
    double dterm_2_dy = (362880.0 * dP9y_exp) - (6531840.0 * P2x * dP7y_exp) - (6531840.0 * dP7y_exp * P2z) + (17146080.0 * P4x * dP5y_exp) + (34292160.0 * P2x * dP5y_exp * P2z) + (17146080.0 * dP5y_exp * P4z) - (9525600.0 * P6x * dP3y_exp) - (28576800.0 * P4x * dP3y_exp * P2z) - (28576800.0 * P2x * dP3y_exp * P4z) - (9525600.0 * dP3y_exp * P6z) + (893025.0 * P8x * dP1y_exp) + (3572100.0 * P6x * dP1y_exp * P2z) + (5358150.0 * P4x * dP1y_exp * P4z) + (3572100.0 * P2x * dP1y_exp * P6z) + (893025.0 * dP1y_exp * P8z);
    double dterm_2_dz = (362880.0 * P9y * dP0z_exp) - (6531840.0 * P2x * P7y * dP0z_exp) - (6531840.0 * P7y * dP2z_exp) + (17146080.0 * P4x * P5y * dP0z_exp) + (34292160.0 * P2x * P5y * dP2z_exp) + (17146080.0 * P5y * dP4z_exp) - (9525600.0 * P6x * P3y * dP0z_exp) - (28576800.0 * P4x * P3y * dP2z_exp) - (28576800.0 * P2x * P3y * dP4z_exp) - (9525600.0 * P3y * dP6z_exp) + (893025.0 * P8x * P1y * dP0z_exp) + (3572100.0 * P6x * P1y * dP2z_exp) + (5358150.0 * P4x * P1y * dP4z_exp) + (3572100.0 * P2x * P1y * dP6z_exp) + (893025.0 * P1y * dP8z_exp);

    double dterm_3_dx = (362880.0 * dP0x_exp * P9z) - (6531840.0 * dP2x_exp * P7z) - (6531840.0 * dP0x_exp * P2y * P7z) + (17146080.0 * dP4x_exp * P5z) + (34292160.0 * dP2x_exp * P2y * P5z) + (17146080.0 * dP0x_exp * P4y * P5z) - (9525600.0 * dP6x_exp * P3z) - (28576800.0 * dP4x_exp * P2y * P3z) - (28576800.0 * dP2x_exp * P4y * P3z) - (9525600.0 * dP0x_exp * P6y * P3z) + (893025.0 * dP8x_exp * P1z) + (3572100.0 * dP6x_exp * P2y * P1z) + (5358150.0 * dP4x_exp * P4y * P1z) + (3572100.0 * dP2x_exp * P6y * P1z) + (893025.0 * dP0x_exp * P8y * P1z);
    double dterm_3_dy = (362880.0 * dP0y_exp * P9z) - (6531840.0 * P2x * dP0y_exp * P7z) - (6531840.0 * dP2y_exp * P7z) + (17146080.0 * P4x * dP0y_exp * P5z) + (34292160.0 * P2x * dP2y_exp * P5z) + (17146080.0 * dP4y_exp * P5z) - (9525600.0 * P6x * dP0y_exp * P3z) - (28576800.0 * P4x * dP2y_exp * P3z) - (28576800.0 * P2x * dP4y_exp * P3z) - (9525600.0 * dP6y_exp * P3z) + (893025.0 * P8x * dP0y_exp * P1z) + (3572100.0 * P6x * dP2y_exp * P1z) + (5358150.0 * P4x * dP4y_exp * P1z) + (3572100.0 * P2x * dP6y_exp * P1z) + (893025.0 * dP8y_exp * P1z);
    double dterm_3_dz = (362880.0 * dP9z_exp) - (6531840.0 * P2x * dP7z_exp) - (6531840.0 * P2y * dP7z_exp) + (17146080.0 * P4x * dP5z_exp) + (34292160.0 * P2x * P2y * dP5z_exp) + (17146080.0 * P4y * dP5z_exp) - (9525600.0 * P6x * dP3z_exp) - (28576800.0 * P4x * P2y * dP3z_exp) - (28576800.0 * P2x * P4y * dP3z_exp) - (9525600.0 * P6y * dP3z_exp) + (893025.0 * P8x * dP1z_exp) + (3572100.0 * P6x * P2y * dP1z_exp) + (5358150.0 * P4x * P4y * dP1z_exp) + (3572100.0 * P2x * P6y * dP1z_exp) + (893025.0 * P8y * dP1z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

}

void calc_solid_MCSH_9_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double P8x = P8(lambda, x0, gamma);
    double P8y = P8(lambda, y0, gamma);
    double P8z = P8(lambda, z0, gamma);

    double P9x = P9(lambda, x0, gamma);
    double P9y = P9(lambda, y0, gamma);
    double P9z = P9(lambda, z0, gamma);

    double term_1 = (1814400.0 * P8x * P1y) - (12700800.0 * P6x * P3y) - (12700800.0 * P6x * P1y * P2z) + (15876000.0 * P4x * P5y) + (31752000.0 * P4x * P3y * P2z) + (15876000.0 * P4x * P1y * P4z) - (3969000.0 * P2x * P7y) - (11907000.0 * P2x * P5y * P2z) - (11907000.0 * P2x * P3y * P4z) - (3969000.0 * P2x * P1y * P6z) + (99225.0 * P9y) + (396900.0 * P7y * P2z) + (595350.0 * P5y * P4z) + (396900.0 * P3y * P6z) + (99225.0 * P1y * P8z);
    double term_2 = (1814400.0 * P1x * P8y) - (12700800.0 * P3x * P6y) - (12700800.0 * P1x * P6y * P2z) + (15876000.0 * P5x * P4y) + (31752000.0 * P3x * P4y * P2z) + (15876000.0 * P1x * P4y * P4z) - (3969000.0 * P7x * P2y) - (11907000.0 * P5x * P2y * P2z) - (11907000.0 * P3x * P2y * P4z) - (3969000.0 * P1x * P2y * P6z) + (99225.0 * P9x) + (396900.0 * P7x * P2z) + (595350.0 * P5x * P4z) + (396900.0 * P3x * P6z) + (99225.0 * P1x * P8z);
    double term_3 = (1814400.0 * P8x * P1z) - (12700800.0 * P6x * P3z) - (12700800.0 * P6x * P2y * P1z) + (15876000.0 * P4x * P5z) + (31752000.0 * P4x * P2y * P3z) + (15876000.0 * P4x * P4y * P1z) - (3969000.0 * P2x * P7z) - (11907000.0 * P2x * P2y * P5z) - (11907000.0 * P2x * P4y * P3z) - (3969000.0 * P2x * P6y * P1z) + (99225.0 * P9z) + (396900.0 * P2y * P7z) + (595350.0 * P4y * P5z) + (396900.0 * P6y * P3z) + (99225.0 * P8y * P1z);
    double term_4 = (1814400.0 * P1x * P8z) - (12700800.0 * P3x * P6z) - (12700800.0 * P1x * P2y * P6z) + (15876000.0 * P5x * P4z) + (31752000.0 * P3x * P2y * P4z) + (15876000.0 * P1x * P4y * P4z) - (3969000.0 * P7x * P2z) - (11907000.0 * P5x * P2y * P2z) - (11907000.0 * P3x * P4y * P2z) - (3969000.0 * P1x * P6y * P2z) + (99225.0 * P9x) + (396900.0 * P7x * P2y) + (595350.0 * P5x * P4y) + (396900.0 * P3x * P6y) + (99225.0 * P1x * P8y);
    double term_5 = (1814400.0 * P8y * P1z) - (12700800.0 * P6y * P3z) - (12700800.0 * P2x * P6y * P1z) + (15876000.0 * P4y * P5z) + (31752000.0 * P2x * P4y * P3z) + (15876000.0 * P4x * P4y * P1z) - (3969000.0 * P2y * P7z) - (11907000.0 * P2x * P2y * P5z) - (11907000.0 * P4x * P2y * P3z) - (3969000.0 * P6x * P2y * P1z) + (99225.0 * P9z) + (396900.0 * P2x * P7z) + (595350.0 * P4x * P5z) + (396900.0 * P6x * P3z) + (99225.0 * P8x * P1z);
    double term_6 = (1814400.0 * P1y * P8z) - (12700800.0 * P3y * P6z) - (12700800.0 * P2x * P1y * P6z) + (15876000.0 * P5y * P4z) + (31752000.0 * P2x * P3y * P4z) + (15876000.0 * P4x * P1y * P4z) - (3969000.0 * P7y * P2z) - (11907000.0 * P2x * P5y * P2z) - (11907000.0 * P4x * P3y * P2z) - (3969000.0 * P6x * P1y * P2z) + (99225.0 * P9y) + (396900.0 * P2x * P7y) + (595350.0 * P4x * P5y) + (396900.0 * P6x * P3y) + (99225.0 * P8x * P1y);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP3x_exp = dP3_exp(P3x, C2, lambda, x0, gamma);
    double dP3y_exp = dP3_exp(P3y, C2, lambda, y0, gamma);
    double dP3z_exp = dP3_exp(P3z, C2, lambda, z0, gamma);

    double dP4x_exp = dP4_exp(P4x, C2, lambda, x0, gamma);
    double dP4y_exp = dP4_exp(P4y, C2, lambda, y0, gamma);
    double dP4z_exp = dP4_exp(P4z, C2, lambda, z0, gamma);

    double dP5x_exp = dP5_exp(P5x, C2, lambda, x0, gamma);
    double dP5y_exp = dP5_exp(P5y, C2, lambda, y0, gamma);
    double dP5z_exp = dP5_exp(P5z, C2, lambda, z0, gamma);

    double dP6x_exp = dP6_exp(P6x, C2, lambda, x0, gamma);
    double dP6y_exp = dP6_exp(P6y, C2, lambda, y0, gamma);
    double dP6z_exp = dP6_exp(P6z, C2, lambda, z0, gamma);

    double dP7x_exp = dP7_exp(P7x, C2, lambda, x0, gamma);
    double dP7y_exp = dP7_exp(P7y, C2, lambda, y0, gamma);
    double dP7z_exp = dP7_exp(P7z, C2, lambda, z0, gamma);

    double dP8x_exp = dP8_exp(P8x, C2, lambda, x0, gamma);
    double dP8y_exp = dP8_exp(P8y, C2, lambda, y0, gamma);
    double dP8z_exp = dP8_exp(P8z, C2, lambda, z0, gamma);

    double dP9x_exp = dP9_exp(P9x, C2, lambda, x0, gamma);
    double dP9y_exp = dP9_exp(P9y, C2, lambda, y0, gamma);
    double dP9z_exp = dP9_exp(P9z, C2, lambda, z0, gamma);

    double dterm_1_dx = (1814400.0 * dP8x_exp * P1y) - (12700800.0 * dP6x_exp * P3y) - (12700800.0 * dP6x_exp * P1y * P2z) + (15876000.0 * dP4x_exp * P5y) + (31752000.0 * dP4x_exp * P3y * P2z) + (15876000.0 * dP4x_exp * P1y * P4z) - (3969000.0 * dP2x_exp * P7y) - (11907000.0 * dP2x_exp * P5y * P2z) - (11907000.0 * dP2x_exp * P3y * P4z) - (3969000.0 * dP2x_exp * P1y * P6z) + (99225.0 * dP0x_exp * P9y) + (396900.0 * dP0x_exp * P7y * P2z) + (595350.0 * dP0x_exp * P5y * P4z) + (396900.0 * dP0x_exp * P3y * P6z) + (99225.0 * dP0x_exp * P1y * P8z);
    double dterm_1_dy = (1814400.0 * P8x * dP1y_exp) - (12700800.0 * P6x * dP3y_exp) - (12700800.0 * P6x * dP1y_exp * P2z) + (15876000.0 * P4x * dP5y_exp) + (31752000.0 * P4x * dP3y_exp * P2z) + (15876000.0 * P4x * dP1y_exp * P4z) - (3969000.0 * P2x * dP7y_exp) - (11907000.0 * P2x * dP5y_exp * P2z) - (11907000.0 * P2x * dP3y_exp * P4z) - (3969000.0 * P2x * dP1y_exp * P6z) + (99225.0 * dP9y_exp) + (396900.0 * dP7y_exp * P2z) + (595350.0 * dP5y_exp * P4z) + (396900.0 * dP3y_exp * P6z) + (99225.0 * dP1y_exp * P8z);
    double dterm_1_dz = (1814400.0 * P8x * P1y * dP0z_exp) - (12700800.0 * P6x * P3y * dP0z_exp) - (12700800.0 * P6x * P1y * dP2z_exp) + (15876000.0 * P4x * P5y * dP0z_exp) + (31752000.0 * P4x * P3y * dP2z_exp) + (15876000.0 * P4x * P1y * dP4z_exp) - (3969000.0 * P2x * P7y * dP0z_exp) - (11907000.0 * P2x * P5y * dP2z_exp) - (11907000.0 * P2x * P3y * dP4z_exp) - (3969000.0 * P2x * P1y * dP6z_exp) + (99225.0 * P9y * dP0z_exp) + (396900.0 * P7y * dP2z_exp) + (595350.0 * P5y * dP4z_exp) + (396900.0 * P3y * dP6z_exp) + (99225.0 * P1y * dP8z_exp);

    double dterm_2_dx = (1814400.0 * dP1x_exp * P8y) - (12700800.0 * dP3x_exp * P6y) - (12700800.0 * dP1x_exp * P6y * P2z) + (15876000.0 * dP5x_exp * P4y) + (31752000.0 * dP3x_exp * P4y * P2z) + (15876000.0 * dP1x_exp * P4y * P4z) - (3969000.0 * dP7x_exp * P2y) - (11907000.0 * dP5x_exp * P2y * P2z) - (11907000.0 * dP3x_exp * P2y * P4z) - (3969000.0 * dP1x_exp * P2y * P6z) + (99225.0 * dP9x_exp) + (396900.0 * dP7x_exp * P2z) + (595350.0 * dP5x_exp * P4z) + (396900.0 * dP3x_exp * P6z) + (99225.0 * dP1x_exp * P8z);
    double dterm_2_dy = (1814400.0 * P1x * dP8y_exp) - (12700800.0 * P3x * dP6y_exp) - (12700800.0 * P1x * dP6y_exp * P2z) + (15876000.0 * P5x * dP4y_exp) + (31752000.0 * P3x * dP4y_exp * P2z) + (15876000.0 * P1x * dP4y_exp * P4z) - (3969000.0 * P7x * dP2y_exp) - (11907000.0 * P5x * dP2y_exp * P2z) - (11907000.0 * P3x * dP2y_exp * P4z) - (3969000.0 * P1x * dP2y_exp * P6z) + (99225.0 * P9x * dP0y_exp) + (396900.0 * P7x * dP0y_exp * P2z) + (595350.0 * P5x * dP0y_exp * P4z) + (396900.0 * P3x * dP0y_exp * P6z) + (99225.0 * P1x * dP0y_exp * P8z);
    double dterm_2_dz = (1814400.0 * P1x * P8y * dP0z_exp) - (12700800.0 * P3x * P6y * dP0z_exp) - (12700800.0 * P1x * P6y * dP2z_exp) + (15876000.0 * P5x * P4y * dP0z_exp) + (31752000.0 * P3x * P4y * dP2z_exp) + (15876000.0 * P1x * P4y * dP4z_exp) - (3969000.0 * P7x * P2y * dP0z_exp) - (11907000.0 * P5x * P2y * dP2z_exp) - (11907000.0 * P3x * P2y * dP4z_exp) - (3969000.0 * P1x * P2y * dP6z_exp) + (99225.0 * P9x * dP0z_exp) + (396900.0 * P7x * dP2z_exp) + (595350.0 * P5x * dP4z_exp) + (396900.0 * P3x * dP6z_exp) + (99225.0 * P1x * dP8z_exp);

    double dterm_3_dx = (1814400.0 * dP8x_exp * P1z) - (12700800.0 * dP6x_exp * P3z) - (12700800.0 * dP6x_exp * P2y * P1z) + (15876000.0 * dP4x_exp * P5z) + (31752000.0 * dP4x_exp * P2y * P3z) + (15876000.0 * dP4x_exp * P4y * P1z) - (3969000.0 * dP2x_exp * P7z) - (11907000.0 * dP2x_exp * P2y * P5z) - (11907000.0 * dP2x_exp * P4y * P3z) - (3969000.0 * dP2x_exp * P6y * P1z) + (99225.0 * dP0x_exp * P9z) + (396900.0 * dP0x_exp * P2y * P7z) + (595350.0 * dP0x_exp * P4y * P5z) + (396900.0 * dP0x_exp * P6y * P3z) + (99225.0 * dP0x_exp * P8y * P1z);
    double dterm_3_dy = (1814400.0 * P8x * dP0y_exp * P1z) - (12700800.0 * P6x * dP0y_exp * P3z) - (12700800.0 * P6x * dP2y_exp * P1z) + (15876000.0 * P4x * dP0y_exp * P5z) + (31752000.0 * P4x * dP2y_exp * P3z) + (15876000.0 * P4x * dP4y_exp * P1z) - (3969000.0 * P2x * dP0y_exp * P7z) - (11907000.0 * P2x * dP2y_exp * P5z) - (11907000.0 * P2x * dP4y_exp * P3z) - (3969000.0 * P2x * dP6y_exp * P1z) + (99225.0 * dP0y_exp * P9z) + (396900.0 * dP2y_exp * P7z) + (595350.0 * dP4y_exp * P5z) + (396900.0 * dP6y_exp * P3z) + (99225.0 * dP8y_exp * P1z);
    double dterm_3_dz = (1814400.0 * P8x * dP1z_exp) - (12700800.0 * P6x * dP3z_exp) - (12700800.0 * P6x * P2y * dP1z_exp) + (15876000.0 * P4x * dP5z_exp) + (31752000.0 * P4x * P2y * dP3z_exp) + (15876000.0 * P4x * P4y * dP1z_exp) - (3969000.0 * P2x * dP7z_exp) - (11907000.0 * P2x * P2y * dP5z_exp) - (11907000.0 * P2x * P4y * dP3z_exp) - (3969000.0 * P2x * P6y * dP1z_exp) + (99225.0 * dP9z_exp) + (396900.0 * P2y * dP7z_exp) + (595350.0 * P4y * dP5z_exp) + (396900.0 * P6y * dP3z_exp) + (99225.0 * P8y * dP1z_exp);

    double dterm_4_dx = (1814400.0 * dP1x_exp * P8z) - (12700800.0 * dP3x_exp * P6z) - (12700800.0 * dP1x_exp * P2y * P6z) + (15876000.0 * dP5x_exp * P4z) + (31752000.0 * dP3x_exp * P2y * P4z) + (15876000.0 * dP1x_exp * P4y * P4z) - (3969000.0 * dP7x_exp * P2z) - (11907000.0 * dP5x_exp * P2y * P2z) - (11907000.0 * dP3x_exp * P4y * P2z) - (3969000.0 * dP1x_exp * P6y * P2z) + (99225.0 * dP9x_exp) + (396900.0 * dP7x_exp * P2y) + (595350.0 * dP5x_exp * P4y) + (396900.0 * dP3x_exp * P6y) + (99225.0 * dP1x_exp * P8y);
    double dterm_4_dy = (1814400.0 * P1x * dP0y_exp * P8z) - (12700800.0 * P3x * dP0y_exp * P6z) - (12700800.0 * P1x * dP2y_exp * P6z) + (15876000.0 * P5x * dP0y_exp * P4z) + (31752000.0 * P3x * dP2y_exp * P4z) + (15876000.0 * P1x * dP4y_exp * P4z) - (3969000.0 * P7x * dP0y_exp * P2z) - (11907000.0 * P5x * dP2y_exp * P2z) - (11907000.0 * P3x * dP4y_exp * P2z) - (3969000.0 * P1x * dP6y_exp * P2z) + (99225.0 * P9x * dP0y_exp) + (396900.0 * P7x * dP2y_exp) + (595350.0 * P5x * dP4y_exp) + (396900.0 * P3x * dP6y_exp) + (99225.0 * P1x * dP8y_exp);
    double dterm_4_dz = (1814400.0 * P1x * dP8z_exp) - (12700800.0 * P3x * dP6z_exp) - (12700800.0 * P1x * P2y * dP6z_exp) + (15876000.0 * P5x * dP4z_exp) + (31752000.0 * P3x * P2y * dP4z_exp) + (15876000.0 * P1x * P4y * dP4z_exp) - (3969000.0 * P7x * dP2z_exp) - (11907000.0 * P5x * P2y * dP2z_exp) - (11907000.0 * P3x * P4y * dP2z_exp) - (3969000.0 * P1x * P6y * dP2z_exp) + (99225.0 * P9x * dP0z_exp) + (396900.0 * P7x * P2y * dP0z_exp) + (595350.0 * P5x * P4y * dP0z_exp) + (396900.0 * P3x * P6y * dP0z_exp) + (99225.0 * P1x * P8y * dP0z_exp);

    double dterm_5_dx = (1814400.0 * dP0x_exp * P8y * P1z) - (12700800.0 * dP0x_exp * P6y * P3z) - (12700800.0 * dP2x_exp * P6y * P1z) + (15876000.0 * dP0x_exp * P4y * P5z) + (31752000.0 * dP2x_exp * P4y * P3z) + (15876000.0 * dP4x_exp * P4y * P1z) - (3969000.0 * dP0x_exp * P2y * P7z) - (11907000.0 * dP2x_exp * P2y * P5z) - (11907000.0 * dP4x_exp * P2y * P3z) - (3969000.0 * dP6x_exp * P2y * P1z) + (99225.0 * dP0x_exp * P9z) + (396900.0 * dP2x_exp * P7z) + (595350.0 * dP4x_exp * P5z) + (396900.0 * dP6x_exp * P3z) + (99225.0 * dP8x_exp * P1z);
    double dterm_5_dy = (1814400.0 * dP8y_exp * P1z) - (12700800.0 * dP6y_exp * P3z) - (12700800.0 * P2x * dP6y_exp * P1z) + (15876000.0 * dP4y_exp * P5z) + (31752000.0 * P2x * dP4y_exp * P3z) + (15876000.0 * P4x * dP4y_exp * P1z) - (3969000.0 * dP2y_exp * P7z) - (11907000.0 * P2x * dP2y_exp * P5z) - (11907000.0 * P4x * dP2y_exp * P3z) - (3969000.0 * P6x * dP2y_exp * P1z) + (99225.0 * dP0y_exp * P9z) + (396900.0 * P2x * dP0y_exp * P7z) + (595350.0 * P4x * dP0y_exp * P5z) + (396900.0 * P6x * dP0y_exp * P3z) + (99225.0 * P8x * dP0y_exp * P1z);
    double dterm_5_dz = (1814400.0 * P8y * dP1z_exp) - (12700800.0 * P6y * dP3z_exp) - (12700800.0 * P2x * P6y * dP1z_exp) + (15876000.0 * P4y * dP5z_exp) + (31752000.0 * P2x * P4y * dP3z_exp) + (15876000.0 * P4x * P4y * dP1z_exp) - (3969000.0 * P2y * dP7z_exp) - (11907000.0 * P2x * P2y * dP5z_exp) - (11907000.0 * P4x * P2y * dP3z_exp) - (3969000.0 * P6x * P2y * dP1z_exp) + (99225.0 * dP9z_exp) + (396900.0 * P2x * dP7z_exp) + (595350.0 * P4x * dP5z_exp) + (396900.0 * P6x * dP3z_exp) + (99225.0 * P8x * dP1z_exp);

    double dterm_6_dx = (1814400.0 * dP0x_exp * P1y * P8z) - (12700800.0 * dP0x_exp * P3y * P6z) - (12700800.0 * dP2x_exp * P1y * P6z) + (15876000.0 * dP0x_exp * P5y * P4z) + (31752000.0 * dP2x_exp * P3y * P4z) + (15876000.0 * dP4x_exp * P1y * P4z) - (3969000.0 * dP0x_exp * P7y * P2z) - (11907000.0 * dP2x_exp * P5y * P2z) - (11907000.0 * dP4x_exp * P3y * P2z) - (3969000.0 * dP6x_exp * P1y * P2z) + (99225.0 * dP0x_exp * P9y) + (396900.0 * dP2x_exp * P7y) + (595350.0 * dP4x_exp * P5y) + (396900.0 * dP6x_exp * P3y) + (99225.0 * dP8x_exp * P1y);
    double dterm_6_dy = (1814400.0 * dP1y_exp * P8z) - (12700800.0 * dP3y_exp * P6z) - (12700800.0 * P2x * dP1y_exp * P6z) + (15876000.0 * dP5y_exp * P4z) + (31752000.0 * P2x * dP3y_exp * P4z) + (15876000.0 * P4x * dP1y_exp * P4z) - (3969000.0 * dP7y_exp * P2z) - (11907000.0 * P2x * dP5y_exp * P2z) - (11907000.0 * P4x * dP3y_exp * P2z) - (3969000.0 * P6x * dP1y_exp * P2z) + (99225.0 * dP9y_exp) + (396900.0 * P2x * dP7y_exp) + (595350.0 * P4x * dP5y_exp) + (396900.0 * P6x * dP3y_exp) + (99225.0 * P8x * dP1y_exp);
    double dterm_6_dz = (1814400.0 * P1y * dP8z_exp) - (12700800.0 * P3y * dP6z_exp) - (12700800.0 * P2x * P1y * dP6z_exp) + (15876000.0 * P5y * dP4z_exp) + (31752000.0 * P2x * P3y * dP4z_exp) + (15876000.0 * P4x * P1y * dP4z_exp) - (3969000.0 * P7y * dP2z_exp) - (11907000.0 * P2x * P5y * dP2z_exp) - (11907000.0 * P4x * P3y * dP2z_exp) - (3969000.0 * P6x * P1y * dP2z_exp) + (99225.0 * P9y * dP0z_exp) + (396900.0 * P2x * P7y * dP0z_exp) + (595350.0 * P4x * P5y * dP0z_exp) + (396900.0 * P6x * P3y * dP0z_exp) + (99225.0 * P8x * P1y * dP0z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

    // dmiu4 dx/dy/dz
    deriv[9] = temp * dterm_4_dx;
    deriv[10] = temp * dterm_4_dy;
    deriv[11] = temp * dterm_4_dz;

    // dmiu5 dx/dy/dz
    deriv[12] = temp * dterm_5_dx;
    deriv[13] = temp * dterm_5_dy;
    deriv[14] = temp * dterm_5_dz;

    // dmiu6 dx/dy/dz
    deriv[15] = temp * dterm_6_dx;
    deriv[16] = temp * dterm_6_dy;
    deriv[17] = temp * dterm_6_dz;

}

void calc_solid_MCSH_9_3(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double P8x = P8(lambda, x0, gamma);
    double P8y = P8(lambda, y0, gamma);
    double P8z = P8(lambda, z0, gamma);

    double P9x = P9(lambda, x0, gamma);
    double P9y = P9(lambda, y0, gamma);
    double P9z = P9(lambda, z0, gamma);

    double term_1 = (5760720.0 * P7x * P2y) - (181440.0 * P9x) + (771120.0 * P7x * P2z) - (17304840.0 * P5x * P4y) - (17146080.0 * P5x * P2y * P2z) + (158760.0 * P5x * P4z) + (10220175.0 * P3x * P6y) + (19745775.0 * P3x * P4y * P2z) + (8831025.0 * P3x * P2y * P4z) - (694575.0 * P3x * P6z) - (992250.0 * P1x * P8y) - (2877525.0 * P1x * P6y * P2z) - (2679075.0 * P1x * P4y * P4z) - (694575.0 * P1x * P2y * P6z) + (99225.0 * P1x * P8z);
    double term_2 = (5760720.0 * P2x * P7y) - (181440.0 * P9y) + (771120.0 * P7y * P2z) - (17304840.0 * P4x * P5y) - (17146080.0 * P2x * P5y * P2z) + (158760.0 * P5y * P4z) + (10220175.0 * P6x * P3y) + (19745775.0 * P4x * P3y * P2z) + (8831025.0 * P2x * P3y * P4z) - (694575.0 * P3y * P6z) - (992250.0 * P8x * P1y) - (2877525.0 * P6x * P1y * P2z) - (2679075.0 * P4x * P1y * P4z) - (694575.0 * P2x * P1y * P6z) + (99225.0 * P1y * P8z);
    double term_3 = (5760720.0 * P7x * P2z) - (181440.0 * P9x) + (771120.0 * P7x * P2y) - (17304840.0 * P5x * P4z) - (17146080.0 * P5x * P2y * P2z) + (158760.0 * P5x * P4y) + (10220175.0 * P3x * P6z) + (19745775.0 * P3x * P2y * P4z) + (8831025.0 * P3x * P4y * P2z) - (694575.0 * P3x * P6y) - (992250.0 * P1x * P8z) - (2877525.0 * P1x * P2y * P6z) - (2679075.0 * P1x * P4y * P4z) - (694575.0 * P1x * P6y * P2z) + (99225.0 * P1x * P8y);
    double term_4 = (5760720.0 * P2x * P7z) - (181440.0 * P9z) + (771120.0 * P2y * P7z) - (17304840.0 * P4x * P5z) - (17146080.0 * P2x * P2y * P5z) + (158760.0 * P4y * P5z) + (10220175.0 * P6x * P3z) + (19745775.0 * P4x * P2y * P3z) + (8831025.0 * P2x * P4y * P3z) - (694575.0 * P6y * P3z) - (992250.0 * P8x * P1z) - (2877525.0 * P6x * P2y * P1z) - (2679075.0 * P4x * P4y * P1z) - (694575.0 * P2x * P6y * P1z) + (99225.0 * P8y * P1z);
    double term_5 = (5760720.0 * P7y * P2z) - (181440.0 * P9y) + (771120.0 * P2x * P7y) - (17304840.0 * P5y * P4z) - (17146080.0 * P2x * P5y * P2z) + (158760.0 * P4x * P5y) + (10220175.0 * P3y * P6z) + (19745775.0 * P2x * P3y * P4z) + (8831025.0 * P4x * P3y * P2z) - (694575.0 * P6x * P3y) - (992250.0 * P1y * P8z) - (2877525.0 * P2x * P1y * P6z) - (2679075.0 * P4x * P1y * P4z) - (694575.0 * P6x * P1y * P2z) + (99225.0 * P8x * P1y);
    double term_6 = (5760720.0 * P2y * P7z) - (181440.0 * P9z) + (771120.0 * P2x * P7z) - (17304840.0 * P4y * P5z) - (17146080.0 * P2x * P2y * P5z) + (158760.0 * P4x * P5z) + (10220175.0 * P6y * P3z) + (19745775.0 * P2x * P4y * P3z) + (8831025.0 * P4x * P2y * P3z) - (694575.0 * P6x * P3z) - (992250.0 * P8y * P1z) - (2877525.0 * P2x * P6y * P1z) - (2679075.0 * P4x * P4y * P1z) - (694575.0 * P6x * P2y * P1z) + (99225.0 * P8x * P1z);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP3x_exp = dP3_exp(P3x, C2, lambda, x0, gamma);
    double dP3y_exp = dP3_exp(P3y, C2, lambda, y0, gamma);
    double dP3z_exp = dP3_exp(P3z, C2, lambda, z0, gamma);

    double dP4x_exp = dP4_exp(P4x, C2, lambda, x0, gamma);
    double dP4y_exp = dP4_exp(P4y, C2, lambda, y0, gamma);
    double dP4z_exp = dP4_exp(P4z, C2, lambda, z0, gamma);

    double dP5x_exp = dP5_exp(P5x, C2, lambda, x0, gamma);
    double dP5y_exp = dP5_exp(P5y, C2, lambda, y0, gamma);
    double dP5z_exp = dP5_exp(P5z, C2, lambda, z0, gamma);

    double dP6x_exp = dP6_exp(P6x, C2, lambda, x0, gamma);
    double dP6y_exp = dP6_exp(P6y, C2, lambda, y0, gamma);
    double dP6z_exp = dP6_exp(P6z, C2, lambda, z0, gamma);

    double dP7x_exp = dP7_exp(P7x, C2, lambda, x0, gamma);
    double dP7y_exp = dP7_exp(P7y, C2, lambda, y0, gamma);
    double dP7z_exp = dP7_exp(P7z, C2, lambda, z0, gamma);

    double dP8x_exp = dP8_exp(P8x, C2, lambda, x0, gamma);
    double dP8y_exp = dP8_exp(P8y, C2, lambda, y0, gamma);
    double dP8z_exp = dP8_exp(P8z, C2, lambda, z0, gamma);

    double dP9x_exp = dP9_exp(P9x, C2, lambda, x0, gamma);
    double dP9y_exp = dP9_exp(P9y, C2, lambda, y0, gamma);
    double dP9z_exp = dP9_exp(P9z, C2, lambda, z0, gamma);

    double dterm_1_dx = (5760720.0 * dP7x_exp * P2y) - (181440.0 * dP9x_exp) + (771120.0 * dP7x_exp * P2z) - (17304840.0 * dP5x_exp * P4y) - (17146080.0 * dP5x_exp * P2y * P2z) + (158760.0 * dP5x_exp * P4z) + (10220175.0 * dP3x_exp * P6y) + (19745775.0 * dP3x_exp * P4y * P2z) + (8831025.0 * dP3x_exp * P2y * P4z) - (694575.0 * dP3x_exp * P6z) - (992250.0 * dP1x_exp * P8y) - (2877525.0 * dP1x_exp * P6y * P2z) - (2679075.0 * dP1x_exp * P4y * P4z) - (694575.0 * dP1x_exp * P2y * P6z) + (99225.0 * dP1x_exp * P8z);
    double dterm_1_dy = (5760720.0 * P7x * dP2y_exp) - (181440.0 * P9x * dP0y_exp) + (771120.0 * P7x * dP0y_exp * P2z) - (17304840.0 * P5x * dP4y_exp) - (17146080.0 * P5x * dP2y_exp * P2z) + (158760.0 * P5x * dP0y_exp * P4z) + (10220175.0 * P3x * dP6y_exp) + (19745775.0 * P3x * dP4y_exp * P2z) + (8831025.0 * P3x * dP2y_exp * P4z) - (694575.0 * P3x * dP0y_exp * P6z) - (992250.0 * P1x * dP8y_exp) - (2877525.0 * P1x * dP6y_exp * P2z) - (2679075.0 * P1x * dP4y_exp * P4z) - (694575.0 * P1x * dP2y_exp * P6z) + (99225.0 * P1x * dP0y_exp * P8z);
    double dterm_1_dz = (5760720.0 * P7x * P2y * dP0z_exp) - (181440.0 * P9x * dP0z_exp) + (771120.0 * P7x * dP2z_exp) - (17304840.0 * P5x * P4y * dP0z_exp) - (17146080.0 * P5x * P2y * dP2z_exp) + (158760.0 * P5x * dP4z_exp) + (10220175.0 * P3x * P6y * dP0z_exp) + (19745775.0 * P3x * P4y * dP2z_exp) + (8831025.0 * P3x * P2y * dP4z_exp) - (694575.0 * P3x * dP6z_exp) - (992250.0 * P1x * P8y * dP0z_exp) - (2877525.0 * P1x * P6y * dP2z_exp) - (2679075.0 * P1x * P4y * dP4z_exp) - (694575.0 * P1x * P2y * dP6z_exp) + (99225.0 * P1x * dP8z_exp);

    double dterm_2_dx = (5760720.0 * dP2x_exp * P7y) - (181440.0 * dP0x_exp * P9y) + (771120.0 * dP0x_exp * P7y * P2z) - (17304840.0 * dP4x_exp * P5y) - (17146080.0 * dP2x_exp * P5y * P2z) + (158760.0 * dP0x_exp * P5y * P4z) + (10220175.0 * dP6x_exp * P3y) + (19745775.0 * dP4x_exp * P3y * P2z) + (8831025.0 * dP2x_exp * P3y * P4z) - (694575.0 * dP0x_exp * P3y * P6z) - (992250.0 * dP8x_exp * P1y) - (2877525.0 * dP6x_exp * P1y * P2z) - (2679075.0 * dP4x_exp * P1y * P4z) - (694575.0 * dP2x_exp * P1y * P6z) + (99225.0 * dP0x_exp * P1y * P8z);
    double dterm_2_dy = (5760720.0 * P2x * dP7y_exp) - (181440.0 * dP9y_exp) + (771120.0 * dP7y_exp * P2z) - (17304840.0 * P4x * dP5y_exp) - (17146080.0 * P2x * dP5y_exp * P2z) + (158760.0 * dP5y_exp * P4z) + (10220175.0 * P6x * dP3y_exp) + (19745775.0 * P4x * dP3y_exp * P2z) + (8831025.0 * P2x * dP3y_exp * P4z) - (694575.0 * dP3y_exp * P6z) - (992250.0 * P8x * dP1y_exp) - (2877525.0 * P6x * dP1y_exp * P2z) - (2679075.0 * P4x * dP1y_exp * P4z) - (694575.0 * P2x * dP1y_exp * P6z) + (99225.0 * dP1y_exp * P8z);
    double dterm_2_dz = (5760720.0 * P2x * P7y * dP0z_exp) - (181440.0 * P9y * dP0z_exp) + (771120.0 * P7y * dP2z_exp) - (17304840.0 * P4x * P5y * dP0z_exp) - (17146080.0 * P2x * P5y * dP2z_exp) + (158760.0 * P5y * dP4z_exp) + (10220175.0 * P6x * P3y * dP0z_exp) + (19745775.0 * P4x * P3y * dP2z_exp) + (8831025.0 * P2x * P3y * dP4z_exp) - (694575.0 * P3y * dP6z_exp) - (992250.0 * P8x * P1y * dP0z_exp) - (2877525.0 * P6x * P1y * dP2z_exp) - (2679075.0 * P4x * P1y * dP4z_exp) - (694575.0 * P2x * P1y * dP6z_exp) + (99225.0 * P1y * dP8z_exp);

    double dterm_3_dx = (5760720.0 * dP7x_exp * P2z) - (181440.0 * dP9x_exp) + (771120.0 * dP7x_exp * P2y) - (17304840.0 * dP5x_exp * P4z) - (17146080.0 * dP5x_exp * P2y * P2z) + (158760.0 * dP5x_exp * P4y) + (10220175.0 * dP3x_exp * P6z) + (19745775.0 * dP3x_exp * P2y * P4z) + (8831025.0 * dP3x_exp * P4y * P2z) - (694575.0 * dP3x_exp * P6y) - (992250.0 * dP1x_exp * P8z) - (2877525.0 * dP1x_exp * P2y * P6z) - (2679075.0 * dP1x_exp * P4y * P4z) - (694575.0 * dP1x_exp * P6y * P2z) + (99225.0 * dP1x_exp * P8y);
    double dterm_3_dy = (5760720.0 * P7x * dP0y_exp * P2z) - (181440.0 * P9x * dP0y_exp) + (771120.0 * P7x * dP2y_exp) - (17304840.0 * P5x * dP0y_exp * P4z) - (17146080.0 * P5x * dP2y_exp * P2z) + (158760.0 * P5x * dP4y_exp) + (10220175.0 * P3x * dP0y_exp * P6z) + (19745775.0 * P3x * dP2y_exp * P4z) + (8831025.0 * P3x * dP4y_exp * P2z) - (694575.0 * P3x * dP6y_exp) - (992250.0 * P1x * dP0y_exp * P8z) - (2877525.0 * P1x * dP2y_exp * P6z) - (2679075.0 * P1x * dP4y_exp * P4z) - (694575.0 * P1x * dP6y_exp * P2z) + (99225.0 * P1x * dP8y_exp);
    double dterm_3_dz = (5760720.0 * P7x * dP2z_exp) - (181440.0 * P9x * dP0z_exp) + (771120.0 * P7x * P2y * dP0z_exp) - (17304840.0 * P5x * dP4z_exp) - (17146080.0 * P5x * P2y * dP2z_exp) + (158760.0 * P5x * P4y * dP0z_exp) + (10220175.0 * P3x * dP6z_exp) + (19745775.0 * P3x * P2y * dP4z_exp) + (8831025.0 * P3x * P4y * dP2z_exp) - (694575.0 * P3x * P6y * dP0z_exp) - (992250.0 * P1x * dP8z_exp) - (2877525.0 * P1x * P2y * dP6z_exp) - (2679075.0 * P1x * P4y * dP4z_exp) - (694575.0 * P1x * P6y * dP2z_exp) + (99225.0 * P1x * P8y * dP0z_exp);

    double dterm_4_dx = (5760720.0 * dP2x_exp * P7z) - (181440.0 * dP0x_exp * P9z) + (771120.0 * dP0x_exp * P2y * P7z) - (17304840.0 * dP4x_exp * P5z) - (17146080.0 * dP2x_exp * P2y * P5z) + (158760.0 * dP0x_exp * P4y * P5z) + (10220175.0 * dP6x_exp * P3z) + (19745775.0 * dP4x_exp * P2y * P3z) + (8831025.0 * dP2x_exp * P4y * P3z) - (694575.0 * dP0x_exp * P6y * P3z) - (992250.0 * dP8x_exp * P1z) - (2877525.0 * dP6x_exp * P2y * P1z) - (2679075.0 * dP4x_exp * P4y * P1z) - (694575.0 * dP2x_exp * P6y * P1z) + (99225.0 * dP0x_exp * P8y * P1z);
    double dterm_4_dy = (5760720.0 * P2x * dP0y_exp * P7z) - (181440.0 * dP0y_exp * P9z) + (771120.0 * dP2y_exp * P7z) - (17304840.0 * P4x * dP0y_exp * P5z) - (17146080.0 * P2x * dP2y_exp * P5z) + (158760.0 * dP4y_exp * P5z) + (10220175.0 * P6x * dP0y_exp * P3z) + (19745775.0 * P4x * dP2y_exp * P3z) + (8831025.0 * P2x * dP4y_exp * P3z) - (694575.0 * dP6y_exp * P3z) - (992250.0 * P8x * dP0y_exp * P1z) - (2877525.0 * P6x * dP2y_exp * P1z) - (2679075.0 * P4x * dP4y_exp * P1z) - (694575.0 * P2x * dP6y_exp * P1z) + (99225.0 * dP8y_exp * P1z);
    double dterm_4_dz = (5760720.0 * P2x * dP7z_exp) - (181440.0 * dP9z_exp) + (771120.0 * P2y * dP7z_exp) - (17304840.0 * P4x * dP5z_exp) - (17146080.0 * P2x * P2y * dP5z_exp) + (158760.0 * P4y * dP5z_exp) + (10220175.0 * P6x * dP3z_exp) + (19745775.0 * P4x * P2y * dP3z_exp) + (8831025.0 * P2x * P4y * dP3z_exp) - (694575.0 * P6y * dP3z_exp) - (992250.0 * P8x * dP1z_exp) - (2877525.0 * P6x * P2y * dP1z_exp) - (2679075.0 * P4x * P4y * dP1z_exp) - (694575.0 * P2x * P6y * dP1z_exp) + (99225.0 * P8y * dP1z_exp);

    double dterm_5_dx = (5760720.0 * dP0x_exp * P7y * P2z) - (181440.0 * dP0x_exp * P9y) + (771120.0 * dP2x_exp * P7y) - (17304840.0 * dP0x_exp * P5y * P4z) - (17146080.0 * dP2x_exp * P5y * P2z) + (158760.0 * dP4x_exp * P5y) + (10220175.0 * dP0x_exp * P3y * P6z) + (19745775.0 * dP2x_exp * P3y * P4z) + (8831025.0 * dP4x_exp * P3y * P2z) - (694575.0 * dP6x_exp * P3y) - (992250.0 * dP0x_exp * P1y * P8z) - (2877525.0 * dP2x_exp * P1y * P6z) - (2679075.0 * dP4x_exp * P1y * P4z) - (694575.0 * dP6x_exp * P1y * P2z) + (99225.0 * dP8x_exp * P1y);
    double dterm_5_dy = (5760720.0 * dP7y_exp * P2z) - (181440.0 * dP9y_exp) + (771120.0 * P2x * dP7y_exp) - (17304840.0 * dP5y_exp * P4z) - (17146080.0 * P2x * dP5y_exp * P2z) + (158760.0 * P4x * dP5y_exp) + (10220175.0 * dP3y_exp * P6z) + (19745775.0 * P2x * dP3y_exp * P4z) + (8831025.0 * P4x * dP3y_exp * P2z) - (694575.0 * P6x * dP3y_exp) - (992250.0 * dP1y_exp * P8z) - (2877525.0 * P2x * dP1y_exp * P6z) - (2679075.0 * P4x * dP1y_exp * P4z) - (694575.0 * P6x * dP1y_exp * P2z) + (99225.0 * P8x * dP1y_exp);
    double dterm_5_dz = (5760720.0 * P7y * dP2z_exp) - (181440.0 * P9y * dP0z_exp) + (771120.0 * P2x * P7y * dP0z_exp) - (17304840.0 * P5y * dP4z_exp) - (17146080.0 * P2x * P5y * dP2z_exp) + (158760.0 * P4x * P5y * dP0z_exp) + (10220175.0 * P3y * dP6z_exp) + (19745775.0 * P2x * P3y * dP4z_exp) + (8831025.0 * P4x * P3y * dP2z_exp) - (694575.0 * P6x * P3y * dP0z_exp) - (992250.0 * P1y * dP8z_exp) - (2877525.0 * P2x * P1y * dP6z_exp) - (2679075.0 * P4x * P1y * dP4z_exp) - (694575.0 * P6x * P1y * dP2z_exp) + (99225.0 * P8x * P1y * dP0z_exp);

    double dterm_6_dx = (5760720.0 * dP0x_exp * P2y * P7z) - (181440.0 * dP0x_exp * P9z) + (771120.0 * dP2x_exp * P7z) - (17304840.0 * dP0x_exp * P4y * P5z) - (17146080.0 * dP2x_exp * P2y * P5z) + (158760.0 * dP4x_exp * P5z) + (10220175.0 * dP0x_exp * P6y * P3z) + (19745775.0 * dP2x_exp * P4y * P3z) + (8831025.0 * dP4x_exp * P2y * P3z) - (694575.0 * dP6x_exp * P3z) - (992250.0 * dP0x_exp * P8y * P1z) - (2877525.0 * dP2x_exp * P6y * P1z) - (2679075.0 * dP4x_exp * P4y * P1z) - (694575.0 * dP6x_exp * P2y * P1z) + (99225.0 * dP8x_exp * P1z);
    double dterm_6_dy = (5760720.0 * dP2y_exp * P7z) - (181440.0 * dP0y_exp * P9z) + (771120.0 * P2x * dP0y_exp * P7z) - (17304840.0 * dP4y_exp * P5z) - (17146080.0 * P2x * dP2y_exp * P5z) + (158760.0 * P4x * dP0y_exp * P5z) + (10220175.0 * dP6y_exp * P3z) + (19745775.0 * P2x * dP4y_exp * P3z) + (8831025.0 * P4x * dP2y_exp * P3z) - (694575.0 * P6x * dP0y_exp * P3z) - (992250.0 * dP8y_exp * P1z) - (2877525.0 * P2x * dP6y_exp * P1z) - (2679075.0 * P4x * dP4y_exp * P1z) - (694575.0 * P6x * dP2y_exp * P1z) + (99225.0 * P8x * dP0y_exp * P1z);
    double dterm_6_dz = (5760720.0 * P2y * dP7z_exp) - (181440.0 * dP9z_exp) + (771120.0 * P2x * dP7z_exp) - (17304840.0 * P4y * dP5z_exp) - (17146080.0 * P2x * P2y * dP5z_exp) + (158760.0 * P4x * dP5z_exp) + (10220175.0 * P6y * dP3z_exp) + (19745775.0 * P2x * P4y * dP3z_exp) + (8831025.0 * P4x * P2y * dP3z_exp) - (694575.0 * P6x * dP3z_exp) - (992250.0 * P8y * dP1z_exp) - (2877525.0 * P2x * P6y * dP1z_exp) - (2679075.0 * P4x * P4y * dP1z_exp) - (694575.0 * P6x * P2y * dP1z_exp) + (99225.0 * P8x * dP1z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

    // dmiu4 dx/dy/dz
    deriv[9] = temp * dterm_4_dx;
    deriv[10] = temp * dterm_4_dy;
    deriv[11] = temp * dterm_4_dz;

    // dmiu5 dx/dy/dz
    deriv[12] = temp * dterm_5_dx;
    deriv[13] = temp * dterm_5_dy;
    deriv[14] = temp * dterm_5_dz;

    // dmiu6 dx/dy/dz
    deriv[15] = temp * dterm_6_dx;
    deriv[16] = temp * dterm_6_dy;
    deriv[17] = temp * dterm_6_dz;

}

void calc_solid_MCSH_9_4(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double term_1 = (4989600.0 * P7x * P1y * P1z) - (17463600.0 * P5x * P3y * P1z) - (17463600.0 * P5x * P1y * P3z) + (10914750.0 * P3x * P5y * P1z) + (21829500.0 * P3x * P3y * P3z) + (10914750.0 * P3x * P1y * P5z) - (1091475.0 * P1x * P7y * P1z) - (3274425.0 * P1x * P5y * P3z) - (3274425.0 * P1x * P3y * P5z) - (1091475.0 * P1x * P1y * P7z);
    double term_2 = (4989600.0 * P1x * P7y * P1z) - (17463600.0 * P3x * P5y * P1z) - (17463600.0 * P1x * P5y * P3z) + (10914750.0 * P5x * P3y * P1z) + (21829500.0 * P3x * P3y * P3z) + (10914750.0 * P1x * P3y * P5z) - (1091475.0 * P7x * P1y * P1z) - (3274425.0 * P5x * P1y * P3z) - (3274425.0 * P3x * P1y * P5z) - (1091475.0 * P1x * P1y * P7z);
    double term_3 = (4989600.0 * P1x * P1y * P7z) - (17463600.0 * P3x * P1y * P5z) - (17463600.0 * P1x * P3y * P5z) + (10914750.0 * P5x * P1y * P3z) + (21829500.0 * P3x * P3y * P3z) + (10914750.0 * P1x * P5y * P3z) - (1091475.0 * P7x * P1y * P1z) - (3274425.0 * P5x * P3y * P1z) - (3274425.0 * P3x * P5y * P1z) - (1091475.0 * P1x * P7y * P1z);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dP3x_exp = dP3_exp(P3x, C2, lambda, x0, gamma);
    double dP3y_exp = dP3_exp(P3y, C2, lambda, y0, gamma);
    double dP3z_exp = dP3_exp(P3z, C2, lambda, z0, gamma);

    double dP5x_exp = dP5_exp(P5x, C2, lambda, x0, gamma);
    double dP5y_exp = dP5_exp(P5y, C2, lambda, y0, gamma);
    double dP5z_exp = dP5_exp(P5z, C2, lambda, z0, gamma);

    double dP7x_exp = dP7_exp(P7x, C2, lambda, x0, gamma);
    double dP7y_exp = dP7_exp(P7y, C2, lambda, y0, gamma);
    double dP7z_exp = dP7_exp(P7z, C2, lambda, z0, gamma);

    double dterm_1_dx = (4989600.0 * dP7x_exp * P1y * P1z) - (17463600.0 * dP5x_exp * P3y * P1z) - (17463600.0 * dP5x_exp * P1y * P3z) + (10914750.0 * dP3x_exp * P5y * P1z) + (21829500.0 * dP3x_exp * P3y * P3z) + (10914750.0 * dP3x_exp * P1y * P5z) - (1091475.0 * dP1x_exp * P7y * P1z) - (3274425.0 * dP1x_exp * P5y * P3z) - (3274425.0 * dP1x_exp * P3y * P5z) - (1091475.0 * dP1x_exp * P1y * P7z);
    double dterm_1_dy = (4989600.0 * P7x * dP1y_exp * P1z) - (17463600.0 * P5x * dP3y_exp * P1z) - (17463600.0 * P5x * dP1y_exp * P3z) + (10914750.0 * P3x * dP5y_exp * P1z) + (21829500.0 * P3x * dP3y_exp * P3z) + (10914750.0 * P3x * dP1y_exp * P5z) - (1091475.0 * P1x * dP7y_exp * P1z) - (3274425.0 * P1x * dP5y_exp * P3z) - (3274425.0 * P1x * dP3y_exp * P5z) - (1091475.0 * P1x * dP1y_exp * P7z);
    double dterm_1_dz = (4989600.0 * P7x * P1y * dP1z_exp) - (17463600.0 * P5x * P3y * dP1z_exp) - (17463600.0 * P5x * P1y * dP3z_exp) + (10914750.0 * P3x * P5y * dP1z_exp) + (21829500.0 * P3x * P3y * dP3z_exp) + (10914750.0 * P3x * P1y * dP5z_exp) - (1091475.0 * P1x * P7y * dP1z_exp) - (3274425.0 * P1x * P5y * dP3z_exp) - (3274425.0 * P1x * P3y * dP5z_exp) - (1091475.0 * P1x * P1y * dP7z_exp);

    double dterm_2_dx = (4989600.0 * dP1x_exp * P7y * P1z) - (17463600.0 * dP3x_exp * P5y * P1z) - (17463600.0 * dP1x_exp * P5y * P3z) + (10914750.0 * dP5x_exp * P3y * P1z) + (21829500.0 * dP3x_exp * P3y * P3z) + (10914750.0 * dP1x_exp * P3y * P5z) - (1091475.0 * dP7x_exp * P1y * P1z) - (3274425.0 * dP5x_exp * P1y * P3z) - (3274425.0 * dP3x_exp * P1y * P5z) - (1091475.0 * dP1x_exp * P1y * P7z);
    double dterm_2_dy = (4989600.0 * P1x * dP7y_exp * P1z) - (17463600.0 * P3x * dP5y_exp * P1z) - (17463600.0 * P1x * dP5y_exp * P3z) + (10914750.0 * P5x * dP3y_exp * P1z) + (21829500.0 * P3x * dP3y_exp * P3z) + (10914750.0 * P1x * dP3y_exp * P5z) - (1091475.0 * P7x * dP1y_exp * P1z) - (3274425.0 * P5x * dP1y_exp * P3z) - (3274425.0 * P3x * dP1y_exp * P5z) - (1091475.0 * P1x * dP1y_exp * P7z);
    double dterm_2_dz = (4989600.0 * P1x * P7y * dP1z_exp) - (17463600.0 * P3x * P5y * dP1z_exp) - (17463600.0 * P1x * P5y * dP3z_exp) + (10914750.0 * P5x * P3y * dP1z_exp) + (21829500.0 * P3x * P3y * dP3z_exp) + (10914750.0 * P1x * P3y * dP5z_exp) - (1091475.0 * P7x * P1y * dP1z_exp) - (3274425.0 * P5x * P1y * dP3z_exp) - (3274425.0 * P3x * P1y * dP5z_exp) - (1091475.0 * P1x * P1y * dP7z_exp);

    double dterm_3_dx = (4989600.0 * dP1x_exp * P1y * P7z) - (17463600.0 * dP3x_exp * P1y * P5z) - (17463600.0 * dP1x_exp * P3y * P5z) + (10914750.0 * dP5x_exp * P1y * P3z) + (21829500.0 * dP3x_exp * P3y * P3z) + (10914750.0 * dP1x_exp * P5y * P3z) - (1091475.0 * dP7x_exp * P1y * P1z) - (3274425.0 * dP5x_exp * P3y * P1z) - (3274425.0 * dP3x_exp * P5y * P1z) - (1091475.0 * dP1x_exp * P7y * P1z);
    double dterm_3_dy = (4989600.0 * P1x * dP1y_exp * P7z) - (17463600.0 * P3x * dP1y_exp * P5z) - (17463600.0 * P1x * dP3y_exp * P5z) + (10914750.0 * P5x * dP1y_exp * P3z) + (21829500.0 * P3x * dP3y_exp * P3z) + (10914750.0 * P1x * dP5y_exp * P3z) - (1091475.0 * P7x * dP1y_exp * P1z) - (3274425.0 * P5x * dP3y_exp * P1z) - (3274425.0 * P3x * dP5y_exp * P1z) - (1091475.0 * P1x * dP7y_exp * P1z);
    double dterm_3_dz = (4989600.0 * P1x * P1y * dP7z_exp) - (17463600.0 * P3x * P1y * dP5z_exp) - (17463600.0 * P1x * P3y * dP5z_exp) + (10914750.0 * P5x * P1y * dP3z_exp) + (21829500.0 * P3x * P3y * dP3z_exp) + (10914750.0 * P1x * P5y * dP3z_exp) - (1091475.0 * P7x * P1y * dP1z_exp) - (3274425.0 * P5x * P3y * dP1z_exp) - (3274425.0 * P3x * P5y * dP1z_exp) - (1091475.0 * P1x * P7y * dP1z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

}

void calc_solid_MCSH_9_5(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double P8x = P8(lambda, x0, gamma);
    double P8y = P8(lambda, y0, gamma);
    double P8z = P8(lambda, z0, gamma);

    double P9x = P9(lambda, x0, gamma);
    double P9y = P9(lambda, y0, gamma);
    double P9z = P9(lambda, z0, gamma);

    double term_1 = (12020400.0 * P6x * P3y) - (1360800.0 * P8x * P1y) + (2041200.0 * P6x * P1y * P2z) - (16584750.0 * P4x * P5y) - (14458500.0 * P4x * P3y * P2z) + (2126250.0 * P4x * P1y * P4z) + (4380075.0 * P2x * P7y) + (7526925.0 * P2x * P5y * P2z) + (1913625.0 * P2x * P3y * P4z) - (1233225.0 * P2x * P1y * P6z) - (113400.0 * P9y) - (297675.0 * P7y * P2z) - (212625.0 * P5y * P4z) + (14175.0 * P3y * P6z) + (42525.0 * P1y * P8z);
    double term_2 = (12020400.0 * P3x * P6y) - (1360800.0 * P1x * P8y) + (2041200.0 * P1x * P6y * P2z) - (16584750.0 * P5x * P4y) - (14458500.0 * P3x * P4y * P2z) + (2126250.0 * P1x * P4y * P4z) + (4380075.0 * P7x * P2y) + (7526925.0 * P5x * P2y * P2z) + (1913625.0 * P3x * P2y * P4z) - (1233225.0 * P1x * P2y * P6z) - (113400.0 * P9x) - (297675.0 * P7x * P2z) - (212625.0 * P5x * P4z) + (14175.0 * P3x * P6z) + (42525.0 * P1x * P8z);
    double term_3 = (12020400.0 * P6x * P3z) - (1360800.0 * P8x * P1z) + (2041200.0 * P6x * P2y * P1z) - (16584750.0 * P4x * P5z) - (14458500.0 * P4x * P2y * P3z) + (2126250.0 * P4x * P4y * P1z) + (4380075.0 * P2x * P7z) + (7526925.0 * P2x * P2y * P5z) + (1913625.0 * P2x * P4y * P3z) - (1233225.0 * P2x * P6y * P1z) - (113400.0 * P9z) - (297675.0 * P2y * P7z) - (212625.0 * P4y * P5z) + (14175.0 * P6y * P3z) + (42525.0 * P8y * P1z);
    double term_4 = (12020400.0 * P3x * P6z) - (1360800.0 * P1x * P8z) + (2041200.0 * P1x * P2y * P6z) - (16584750.0 * P5x * P4z) - (14458500.0 * P3x * P2y * P4z) + (2126250.0 * P1x * P4y * P4z) + (4380075.0 * P7x * P2z) + (7526925.0 * P5x * P2y * P2z) + (1913625.0 * P3x * P4y * P2z) - (1233225.0 * P1x * P6y * P2z) - (113400.0 * P9x) - (297675.0 * P7x * P2y) - (212625.0 * P5x * P4y) + (14175.0 * P3x * P6y) + (42525.0 * P1x * P8y);
    double term_5 = (12020400.0 * P6y * P3z) - (1360800.0 * P8y * P1z) + (2041200.0 * P2x * P6y * P1z) - (16584750.0 * P4y * P5z) - (14458500.0 * P2x * P4y * P3z) + (2126250.0 * P4x * P4y * P1z) + (4380075.0 * P2y * P7z) + (7526925.0 * P2x * P2y * P5z) + (1913625.0 * P4x * P2y * P3z) - (1233225.0 * P6x * P2y * P1z) - (113400.0 * P9z) - (297675.0 * P2x * P7z) - (212625.0 * P4x * P5z) + (14175.0 * P6x * P3z) + (42525.0 * P8x * P1z);
    double term_6 = (12020400.0 * P3y * P6z) - (1360800.0 * P1y * P8z) + (2041200.0 * P2x * P1y * P6z) - (16584750.0 * P5y * P4z) - (14458500.0 * P2x * P3y * P4z) + (2126250.0 * P4x * P1y * P4z) + (4380075.0 * P7y * P2z) + (7526925.0 * P2x * P5y * P2z) + (1913625.0 * P4x * P3y * P2z) - (1233225.0 * P6x * P1y * P2z) - (113400.0 * P9y) - (297675.0 * P2x * P7y) - (212625.0 * P4x * P5y) + (14175.0 * P6x * P3y) + (42525.0 * P8x * P1y);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP3x_exp = dP3_exp(P3x, C2, lambda, x0, gamma);
    double dP3y_exp = dP3_exp(P3y, C2, lambda, y0, gamma);
    double dP3z_exp = dP3_exp(P3z, C2, lambda, z0, gamma);

    double dP4x_exp = dP4_exp(P4x, C2, lambda, x0, gamma);
    double dP4y_exp = dP4_exp(P4y, C2, lambda, y0, gamma);
    double dP4z_exp = dP4_exp(P4z, C2, lambda, z0, gamma);

    double dP5x_exp = dP5_exp(P5x, C2, lambda, x0, gamma);
    double dP5y_exp = dP5_exp(P5y, C2, lambda, y0, gamma);
    double dP5z_exp = dP5_exp(P5z, C2, lambda, z0, gamma);

    double dP6x_exp = dP6_exp(P6x, C2, lambda, x0, gamma);
    double dP6y_exp = dP6_exp(P6y, C2, lambda, y0, gamma);
    double dP6z_exp = dP6_exp(P6z, C2, lambda, z0, gamma);

    double dP7x_exp = dP7_exp(P7x, C2, lambda, x0, gamma);
    double dP7y_exp = dP7_exp(P7y, C2, lambda, y0, gamma);
    double dP7z_exp = dP7_exp(P7z, C2, lambda, z0, gamma);

    double dP8x_exp = dP8_exp(P8x, C2, lambda, x0, gamma);
    double dP8y_exp = dP8_exp(P8y, C2, lambda, y0, gamma);
    double dP8z_exp = dP8_exp(P8z, C2, lambda, z0, gamma);

    double dP9x_exp = dP9_exp(P9x, C2, lambda, x0, gamma);
    double dP9y_exp = dP9_exp(P9y, C2, lambda, y0, gamma);
    double dP9z_exp = dP9_exp(P9z, C2, lambda, z0, gamma);

    double dterm_1_dx = (12020400.0 * dP6x_exp * P3y) - (1360800.0 * dP8x_exp * P1y) + (2041200.0 * dP6x_exp * P1y * P2z) - (16584750.0 * dP4x_exp * P5y) - (14458500.0 * dP4x_exp * P3y * P2z) + (2126250.0 * dP4x_exp * P1y * P4z) + (4380075.0 * dP2x_exp * P7y) + (7526925.0 * dP2x_exp * P5y * P2z) + (1913625.0 * dP2x_exp * P3y * P4z) - (1233225.0 * dP2x_exp * P1y * P6z) - (113400.0 * dP0x_exp * P9y) - (297675.0 * dP0x_exp * P7y * P2z) - (212625.0 * dP0x_exp * P5y * P4z) + (14175.0 * dP0x_exp * P3y * P6z) + (42525.0 * dP0x_exp * P1y * P8z);
    double dterm_1_dy = (12020400.0 * P6x * dP3y_exp) - (1360800.0 * P8x * dP1y_exp) + (2041200.0 * P6x * dP1y_exp * P2z) - (16584750.0 * P4x * dP5y_exp) - (14458500.0 * P4x * dP3y_exp * P2z) + (2126250.0 * P4x * dP1y_exp * P4z) + (4380075.0 * P2x * dP7y_exp) + (7526925.0 * P2x * dP5y_exp * P2z) + (1913625.0 * P2x * dP3y_exp * P4z) - (1233225.0 * P2x * dP1y_exp * P6z) - (113400.0 * dP9y_exp) - (297675.0 * dP7y_exp * P2z) - (212625.0 * dP5y_exp * P4z) + (14175.0 * dP3y_exp * P6z) + (42525.0 * dP1y_exp * P8z);
    double dterm_1_dz = (12020400.0 * P6x * P3y * dP0z_exp) - (1360800.0 * P8x * P1y * dP0z_exp) + (2041200.0 * P6x * P1y * dP2z_exp) - (16584750.0 * P4x * P5y * dP0z_exp) - (14458500.0 * P4x * P3y * dP2z_exp) + (2126250.0 * P4x * P1y * dP4z_exp) + (4380075.0 * P2x * P7y * dP0z_exp) + (7526925.0 * P2x * P5y * dP2z_exp) + (1913625.0 * P2x * P3y * dP4z_exp) - (1233225.0 * P2x * P1y * dP6z_exp) - (113400.0 * P9y * dP0z_exp) - (297675.0 * P7y * dP2z_exp) - (212625.0 * P5y * dP4z_exp) + (14175.0 * P3y * dP6z_exp) + (42525.0 * P1y * dP8z_exp);

    double dterm_2_dx = (12020400.0 * dP3x_exp * P6y) - (1360800.0 * dP1x_exp * P8y) + (2041200.0 * dP1x_exp * P6y * P2z) - (16584750.0 * dP5x_exp * P4y) - (14458500.0 * dP3x_exp * P4y * P2z) + (2126250.0 * dP1x_exp * P4y * P4z) + (4380075.0 * dP7x_exp * P2y) + (7526925.0 * dP5x_exp * P2y * P2z) + (1913625.0 * dP3x_exp * P2y * P4z) - (1233225.0 * dP1x_exp * P2y * P6z) - (113400.0 * dP9x_exp) - (297675.0 * dP7x_exp * P2z) - (212625.0 * dP5x_exp * P4z) + (14175.0 * dP3x_exp * P6z) + (42525.0 * dP1x_exp * P8z);
    double dterm_2_dy = (12020400.0 * P3x * dP6y_exp) - (1360800.0 * P1x * dP8y_exp) + (2041200.0 * P1x * dP6y_exp * P2z) - (16584750.0 * P5x * dP4y_exp) - (14458500.0 * P3x * dP4y_exp * P2z) + (2126250.0 * P1x * dP4y_exp * P4z) + (4380075.0 * P7x * dP2y_exp) + (7526925.0 * P5x * dP2y_exp * P2z) + (1913625.0 * P3x * dP2y_exp * P4z) - (1233225.0 * P1x * dP2y_exp * P6z) - (113400.0 * P9x * dP0y_exp) - (297675.0 * P7x * dP0y_exp * P2z) - (212625.0 * P5x * dP0y_exp * P4z) + (14175.0 * P3x * dP0y_exp * P6z) + (42525.0 * P1x * dP0y_exp * P8z);
    double dterm_2_dz = (12020400.0 * P3x * P6y * dP0z_exp) - (1360800.0 * P1x * P8y * dP0z_exp) + (2041200.0 * P1x * P6y * dP2z_exp) - (16584750.0 * P5x * P4y * dP0z_exp) - (14458500.0 * P3x * P4y * dP2z_exp) + (2126250.0 * P1x * P4y * dP4z_exp) + (4380075.0 * P7x * P2y * dP0z_exp) + (7526925.0 * P5x * P2y * dP2z_exp) + (1913625.0 * P3x * P2y * dP4z_exp) - (1233225.0 * P1x * P2y * dP6z_exp) - (113400.0 * P9x * dP0z_exp) - (297675.0 * P7x * dP2z_exp) - (212625.0 * P5x * dP4z_exp) + (14175.0 * P3x * dP6z_exp) + (42525.0 * P1x * dP8z_exp);

    double dterm_3_dx = (12020400.0 * dP6x_exp * P3z) - (1360800.0 * dP8x_exp * P1z) + (2041200.0 * dP6x_exp * P2y * P1z) - (16584750.0 * dP4x_exp * P5z) - (14458500.0 * dP4x_exp * P2y * P3z) + (2126250.0 * dP4x_exp * P4y * P1z) + (4380075.0 * dP2x_exp * P7z) + (7526925.0 * dP2x_exp * P2y * P5z) + (1913625.0 * dP2x_exp * P4y * P3z) - (1233225.0 * dP2x_exp * P6y * P1z) - (113400.0 * dP0x_exp * P9z) - (297675.0 * dP0x_exp * P2y * P7z) - (212625.0 * dP0x_exp * P4y * P5z) + (14175.0 * dP0x_exp * P6y * P3z) + (42525.0 * dP0x_exp * P8y * P1z);
    double dterm_3_dy = (12020400.0 * P6x * dP0y_exp * P3z) - (1360800.0 * P8x * dP0y_exp * P1z) + (2041200.0 * P6x * dP2y_exp * P1z) - (16584750.0 * P4x * dP0y_exp * P5z) - (14458500.0 * P4x * dP2y_exp * P3z) + (2126250.0 * P4x * dP4y_exp * P1z) + (4380075.0 * P2x * dP0y_exp * P7z) + (7526925.0 * P2x * dP2y_exp * P5z) + (1913625.0 * P2x * dP4y_exp * P3z) - (1233225.0 * P2x * dP6y_exp * P1z) - (113400.0 * dP0y_exp * P9z) - (297675.0 * dP2y_exp * P7z) - (212625.0 * dP4y_exp * P5z) + (14175.0 * dP6y_exp * P3z) + (42525.0 * dP8y_exp * P1z);
    double dterm_3_dz = (12020400.0 * P6x * dP3z_exp) - (1360800.0 * P8x * dP1z_exp) + (2041200.0 * P6x * P2y * dP1z_exp) - (16584750.0 * P4x * dP5z_exp) - (14458500.0 * P4x * P2y * dP3z_exp) + (2126250.0 * P4x * P4y * dP1z_exp) + (4380075.0 * P2x * dP7z_exp) + (7526925.0 * P2x * P2y * dP5z_exp) + (1913625.0 * P2x * P4y * dP3z_exp) - (1233225.0 * P2x * P6y * dP1z_exp) - (113400.0 * dP9z_exp) - (297675.0 * P2y * dP7z_exp) - (212625.0 * P4y * dP5z_exp) + (14175.0 * P6y * dP3z_exp) + (42525.0 * P8y * dP1z_exp);

    double dterm_4_dx = (12020400.0 * dP3x_exp * P6z) - (1360800.0 * dP1x_exp * P8z) + (2041200.0 * dP1x_exp * P2y * P6z) - (16584750.0 * dP5x_exp * P4z) - (14458500.0 * dP3x_exp * P2y * P4z) + (2126250.0 * dP1x_exp * P4y * P4z) + (4380075.0 * dP7x_exp * P2z) + (7526925.0 * dP5x_exp * P2y * P2z) + (1913625.0 * dP3x_exp * P4y * P2z) - (1233225.0 * dP1x_exp * P6y * P2z) - (113400.0 * dP9x_exp) - (297675.0 * dP7x_exp * P2y) - (212625.0 * dP5x_exp * P4y) + (14175.0 * dP3x_exp * P6y) + (42525.0 * dP1x_exp * P8y);
    double dterm_4_dy = (12020400.0 * P3x * dP0y_exp * P6z) - (1360800.0 * P1x * dP0y_exp * P8z) + (2041200.0 * P1x * dP2y_exp * P6z) - (16584750.0 * P5x * dP0y_exp * P4z) - (14458500.0 * P3x * dP2y_exp * P4z) + (2126250.0 * P1x * dP4y_exp * P4z) + (4380075.0 * P7x * dP0y_exp * P2z) + (7526925.0 * P5x * dP2y_exp * P2z) + (1913625.0 * P3x * dP4y_exp * P2z) - (1233225.0 * P1x * dP6y_exp * P2z) - (113400.0 * P9x * dP0y_exp) - (297675.0 * P7x * dP2y_exp) - (212625.0 * P5x * dP4y_exp) + (14175.0 * P3x * dP6y_exp) + (42525.0 * P1x * dP8y_exp);
    double dterm_4_dz = (12020400.0 * P3x * dP6z_exp) - (1360800.0 * P1x * dP8z_exp) + (2041200.0 * P1x * P2y * dP6z_exp) - (16584750.0 * P5x * dP4z_exp) - (14458500.0 * P3x * P2y * dP4z_exp) + (2126250.0 * P1x * P4y * dP4z_exp) + (4380075.0 * P7x * dP2z_exp) + (7526925.0 * P5x * P2y * dP2z_exp) + (1913625.0 * P3x * P4y * dP2z_exp) - (1233225.0 * P1x * P6y * dP2z_exp) - (113400.0 * P9x * dP0z_exp) - (297675.0 * P7x * P2y * dP0z_exp) - (212625.0 * P5x * P4y * dP0z_exp) + (14175.0 * P3x * P6y * dP0z_exp) + (42525.0 * P1x * P8y * dP0z_exp);

    double dterm_5_dx = (12020400.0 * dP0x_exp * P6y * P3z) - (1360800.0 * dP0x_exp * P8y * P1z) + (2041200.0 * dP2x_exp * P6y * P1z) - (16584750.0 * dP0x_exp * P4y * P5z) - (14458500.0 * dP2x_exp * P4y * P3z) + (2126250.0 * dP4x_exp * P4y * P1z) + (4380075.0 * dP0x_exp * P2y * P7z) + (7526925.0 * dP2x_exp * P2y * P5z) + (1913625.0 * dP4x_exp * P2y * P3z) - (1233225.0 * dP6x_exp * P2y * P1z) - (113400.0 * dP0x_exp * P9z) - (297675.0 * dP2x_exp * P7z) - (212625.0 * dP4x_exp * P5z) + (14175.0 * dP6x_exp * P3z) + (42525.0 * dP8x_exp * P1z);
    double dterm_5_dy = (12020400.0 * dP6y_exp * P3z) - (1360800.0 * dP8y_exp * P1z) + (2041200.0 * P2x * dP6y_exp * P1z) - (16584750.0 * dP4y_exp * P5z) - (14458500.0 * P2x * dP4y_exp * P3z) + (2126250.0 * P4x * dP4y_exp * P1z) + (4380075.0 * dP2y_exp * P7z) + (7526925.0 * P2x * dP2y_exp * P5z) + (1913625.0 * P4x * dP2y_exp * P3z) - (1233225.0 * P6x * dP2y_exp * P1z) - (113400.0 * dP0y_exp * P9z) - (297675.0 * P2x * dP0y_exp * P7z) - (212625.0 * P4x * dP0y_exp * P5z) + (14175.0 * P6x * dP0y_exp * P3z) + (42525.0 * P8x * dP0y_exp * P1z);
    double dterm_5_dz = (12020400.0 * P6y * dP3z_exp) - (1360800.0 * P8y * dP1z_exp) + (2041200.0 * P2x * P6y * dP1z_exp) - (16584750.0 * P4y * dP5z_exp) - (14458500.0 * P2x * P4y * dP3z_exp) + (2126250.0 * P4x * P4y * dP1z_exp) + (4380075.0 * P2y * dP7z_exp) + (7526925.0 * P2x * P2y * dP5z_exp) + (1913625.0 * P4x * P2y * dP3z_exp) - (1233225.0 * P6x * P2y * dP1z_exp) - (113400.0 * dP9z_exp) - (297675.0 * P2x * dP7z_exp) - (212625.0 * P4x * dP5z_exp) + (14175.0 * P6x * dP3z_exp) + (42525.0 * P8x * dP1z_exp);

    double dterm_6_dx = (12020400.0 * dP0x_exp * P3y * P6z) - (1360800.0 * dP0x_exp * P1y * P8z) + (2041200.0 * dP2x_exp * P1y * P6z) - (16584750.0 * dP0x_exp * P5y * P4z) - (14458500.0 * dP2x_exp * P3y * P4z) + (2126250.0 * dP4x_exp * P1y * P4z) + (4380075.0 * dP0x_exp * P7y * P2z) + (7526925.0 * dP2x_exp * P5y * P2z) + (1913625.0 * dP4x_exp * P3y * P2z) - (1233225.0 * dP6x_exp * P1y * P2z) - (113400.0 * dP0x_exp * P9y) - (297675.0 * dP2x_exp * P7y) - (212625.0 * dP4x_exp * P5y) + (14175.0 * dP6x_exp * P3y) + (42525.0 * dP8x_exp * P1y);
    double dterm_6_dy = (12020400.0 * dP3y_exp * P6z) - (1360800.0 * dP1y_exp * P8z) + (2041200.0 * P2x * dP1y_exp * P6z) - (16584750.0 * dP5y_exp * P4z) - (14458500.0 * P2x * dP3y_exp * P4z) + (2126250.0 * P4x * dP1y_exp * P4z) + (4380075.0 * dP7y_exp * P2z) + (7526925.0 * P2x * dP5y_exp * P2z) + (1913625.0 * P4x * dP3y_exp * P2z) - (1233225.0 * P6x * dP1y_exp * P2z) - (113400.0 * dP9y_exp) - (297675.0 * P2x * dP7y_exp) - (212625.0 * P4x * dP5y_exp) + (14175.0 * P6x * dP3y_exp) + (42525.0 * P8x * dP1y_exp);
    double dterm_6_dz = (12020400.0 * P3y * dP6z_exp) - (1360800.0 * P1y * dP8z_exp) + (2041200.0 * P2x * P1y * dP6z_exp) - (16584750.0 * P5y * dP4z_exp) - (14458500.0 * P2x * P3y * dP4z_exp) + (2126250.0 * P4x * P1y * dP4z_exp) + (4380075.0 * P7y * dP2z_exp) + (7526925.0 * P2x * P5y * dP2z_exp) + (1913625.0 * P4x * P3y * dP2z_exp) - (1233225.0 * P6x * P1y * dP2z_exp) - (113400.0 * P9y * dP0z_exp) - (297675.0 * P2x * P7y * dP0z_exp) - (212625.0 * P4x * P5y * dP0z_exp) + (14175.0 * P6x * P3y * dP0z_exp) + (42525.0 * P8x * P1y * dP0z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

    // dmiu4 dx/dy/dz
    deriv[9] = temp * dterm_4_dx;
    deriv[10] = temp * dterm_4_dy;
    deriv[11] = temp * dterm_4_dz;

    // dmiu5 dx/dy/dz
    deriv[12] = temp * dterm_5_dx;
    deriv[13] = temp * dterm_5_dy;
    deriv[14] = temp * dterm_5_dz;

    // dmiu6 dx/dy/dz
    deriv[15] = temp * dterm_6_dx;
    deriv[16] = temp * dterm_6_dy;
    deriv[17] = temp * dterm_6_dz;

}

void calc_solid_MCSH_9_6(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double P8x = P8(lambda, x0, gamma);
    double P8y = P8(lambda, y0, gamma);
    double P8z = P8(lambda, z0, gamma);

    double P9x = P9(lambda, x0, gamma);
    double P9y = P9(lambda, y0, gamma);
    double P9z = P9(lambda, z0, gamma);

    double term_1 = (10659600.0 * P6x * P2y * P1z) - (453600.0 * P8x * P1z) + (680400.0 * P6x * P3z) - (18002250.0 * P4x * P4y * P1z) - (17293500.0 * P4x * P2y * P3z) + (708750.0 * P4x * P5z) + (5202225.0 * P2x * P6y * P1z) + (9993375.0 * P2x * P4y * P3z) + (4380075.0 * P2x * P2y * P5z) - (411075.0 * P2x * P7z) - (141750.0 * P8y * P1z) - (411075.0 * P6y * P3z) - (382725.0 * P4y * P5z) - (99225.0 * P2y * P7z) + (14175.0 * P9z);
    double term_2 = (10659600.0 * P2x * P6y * P1z) - (453600.0 * P8y * P1z) + (680400.0 * P6y * P3z) - (18002250.0 * P4x * P4y * P1z) - (17293500.0 * P2x * P4y * P3z) + (708750.0 * P4y * P5z) + (5202225.0 * P6x * P2y * P1z) + (9993375.0 * P4x * P2y * P3z) + (4380075.0 * P2x * P2y * P5z) - (411075.0 * P2y * P7z) - (141750.0 * P8x * P1z) - (411075.0 * P6x * P3z) - (382725.0 * P4x * P5z) - (99225.0 * P2x * P7z) + (14175.0 * P9z);
    double term_3 = (10659600.0 * P6x * P1y * P2z) - (453600.0 * P8x * P1y) + (680400.0 * P6x * P3y) - (18002250.0 * P4x * P1y * P4z) - (17293500.0 * P4x * P3y * P2z) + (708750.0 * P4x * P5y) + (5202225.0 * P2x * P1y * P6z) + (9993375.0 * P2x * P3y * P4z) + (4380075.0 * P2x * P5y * P2z) - (411075.0 * P2x * P7y) - (141750.0 * P1y * P8z) - (411075.0 * P3y * P6z) - (382725.0 * P5y * P4z) - (99225.0 * P7y * P2z) + (14175.0 * P9y);
    double term_4 = (10659600.0 * P2x * P1y * P6z) - (453600.0 * P1y * P8z) + (680400.0 * P3y * P6z) - (18002250.0 * P4x * P1y * P4z) - (17293500.0 * P2x * P3y * P4z) + (708750.0 * P5y * P4z) + (5202225.0 * P6x * P1y * P2z) + (9993375.0 * P4x * P3y * P2z) + (4380075.0 * P2x * P5y * P2z) - (411075.0 * P7y * P2z) - (141750.0 * P8x * P1y) - (411075.0 * P6x * P3y) - (382725.0 * P4x * P5y) - (99225.0 * P2x * P7y) + (14175.0 * P9y);
    double term_5 = (10659600.0 * P1x * P6y * P2z) - (453600.0 * P1x * P8y) + (680400.0 * P3x * P6y) - (18002250.0 * P1x * P4y * P4z) - (17293500.0 * P3x * P4y * P2z) + (708750.0 * P5x * P4y) + (5202225.0 * P1x * P2y * P6z) + (9993375.0 * P3x * P2y * P4z) + (4380075.0 * P5x * P2y * P2z) - (411075.0 * P7x * P2y) - (141750.0 * P1x * P8z) - (411075.0 * P3x * P6z) - (382725.0 * P5x * P4z) - (99225.0 * P7x * P2z) + (14175.0 * P9x);
    double term_6 = (10659600.0 * P1x * P2y * P6z) - (453600.0 * P1x * P8z) + (680400.0 * P3x * P6z) - (18002250.0 * P1x * P4y * P4z) - (17293500.0 * P3x * P2y * P4z) + (708750.0 * P5x * P4z) + (5202225.0 * P1x * P6y * P2z) + (9993375.0 * P3x * P4y * P2z) + (4380075.0 * P5x * P2y * P2z) - (411075.0 * P7x * P2z) - (141750.0 * P1x * P8y) - (411075.0 * P3x * P6y) - (382725.0 * P5x * P4y) - (99225.0 * P7x * P2y) + (14175.0 * P9x);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP3x_exp = dP3_exp(P3x, C2, lambda, x0, gamma);
    double dP3y_exp = dP3_exp(P3y, C2, lambda, y0, gamma);
    double dP3z_exp = dP3_exp(P3z, C2, lambda, z0, gamma);

    double dP4x_exp = dP4_exp(P4x, C2, lambda, x0, gamma);
    double dP4y_exp = dP4_exp(P4y, C2, lambda, y0, gamma);
    double dP4z_exp = dP4_exp(P4z, C2, lambda, z0, gamma);

    double dP5x_exp = dP5_exp(P5x, C2, lambda, x0, gamma);
    double dP5y_exp = dP5_exp(P5y, C2, lambda, y0, gamma);
    double dP5z_exp = dP5_exp(P5z, C2, lambda, z0, gamma);

    double dP6x_exp = dP6_exp(P6x, C2, lambda, x0, gamma);
    double dP6y_exp = dP6_exp(P6y, C2, lambda, y0, gamma);
    double dP6z_exp = dP6_exp(P6z, C2, lambda, z0, gamma);

    double dP7x_exp = dP7_exp(P7x, C2, lambda, x0, gamma);
    double dP7y_exp = dP7_exp(P7y, C2, lambda, y0, gamma);
    double dP7z_exp = dP7_exp(P7z, C2, lambda, z0, gamma);

    double dP8x_exp = dP8_exp(P8x, C2, lambda, x0, gamma);
    double dP8y_exp = dP8_exp(P8y, C2, lambda, y0, gamma);
    double dP8z_exp = dP8_exp(P8z, C2, lambda, z0, gamma);

    double dP9x_exp = dP9_exp(P9x, C2, lambda, x0, gamma);
    double dP9y_exp = dP9_exp(P9y, C2, lambda, y0, gamma);
    double dP9z_exp = dP9_exp(P9z, C2, lambda, z0, gamma);

    double dterm_1_dx = (10659600.0 * dP6x_exp * P2y * P1z) - (453600.0 * dP8x_exp * P1z) + (680400.0 * dP6x_exp * P3z) - (18002250.0 * dP4x_exp * P4y * P1z) - (17293500.0 * dP4x_exp * P2y * P3z) + (708750.0 * dP4x_exp * P5z) + (5202225.0 * dP2x_exp * P6y * P1z) + (9993375.0 * dP2x_exp * P4y * P3z) + (4380075.0 * dP2x_exp * P2y * P5z) - (411075.0 * dP2x_exp * P7z) - (141750.0 * dP0x_exp * P8y * P1z) - (411075.0 * dP0x_exp * P6y * P3z) - (382725.0 * dP0x_exp * P4y * P5z) - (99225.0 * dP0x_exp * P2y * P7z) + (14175.0 * dP0x_exp * P9z);
    double dterm_1_dy = (10659600.0 * P6x * dP2y_exp * P1z) - (453600.0 * P8x * dP0y_exp * P1z) + (680400.0 * P6x * dP0y_exp * P3z) - (18002250.0 * P4x * dP4y_exp * P1z) - (17293500.0 * P4x * dP2y_exp * P3z) + (708750.0 * P4x * dP0y_exp * P5z) + (5202225.0 * P2x * dP6y_exp * P1z) + (9993375.0 * P2x * dP4y_exp * P3z) + (4380075.0 * P2x * dP2y_exp * P5z) - (411075.0 * P2x * dP0y_exp * P7z) - (141750.0 * dP8y_exp * P1z) - (411075.0 * dP6y_exp * P3z) - (382725.0 * dP4y_exp * P5z) - (99225.0 * dP2y_exp * P7z) + (14175.0 * dP0y_exp * P9z);
    double dterm_1_dz = (10659600.0 * P6x * P2y * dP1z_exp) - (453600.0 * P8x * dP1z_exp) + (680400.0 * P6x * dP3z_exp) - (18002250.0 * P4x * P4y * dP1z_exp) - (17293500.0 * P4x * P2y * dP3z_exp) + (708750.0 * P4x * dP5z_exp) + (5202225.0 * P2x * P6y * dP1z_exp) + (9993375.0 * P2x * P4y * dP3z_exp) + (4380075.0 * P2x * P2y * dP5z_exp) - (411075.0 * P2x * dP7z_exp) - (141750.0 * P8y * dP1z_exp) - (411075.0 * P6y * dP3z_exp) - (382725.0 * P4y * dP5z_exp) - (99225.0 * P2y * dP7z_exp) + (14175.0 * dP9z_exp);

    double dterm_2_dx = (10659600.0 * dP2x_exp * P6y * P1z) - (453600.0 * dP0x_exp * P8y * P1z) + (680400.0 * dP0x_exp * P6y * P3z) - (18002250.0 * dP4x_exp * P4y * P1z) - (17293500.0 * dP2x_exp * P4y * P3z) + (708750.0 * dP0x_exp * P4y * P5z) + (5202225.0 * dP6x_exp * P2y * P1z) + (9993375.0 * dP4x_exp * P2y * P3z) + (4380075.0 * dP2x_exp * P2y * P5z) - (411075.0 * dP0x_exp * P2y * P7z) - (141750.0 * dP8x_exp * P1z) - (411075.0 * dP6x_exp * P3z) - (382725.0 * dP4x_exp * P5z) - (99225.0 * dP2x_exp * P7z) + (14175.0 * dP0x_exp * P9z);
    double dterm_2_dy = (10659600.0 * P2x * dP6y_exp * P1z) - (453600.0 * dP8y_exp * P1z) + (680400.0 * dP6y_exp * P3z) - (18002250.0 * P4x * dP4y_exp * P1z) - (17293500.0 * P2x * dP4y_exp * P3z) + (708750.0 * dP4y_exp * P5z) + (5202225.0 * P6x * dP2y_exp * P1z) + (9993375.0 * P4x * dP2y_exp * P3z) + (4380075.0 * P2x * dP2y_exp * P5z) - (411075.0 * dP2y_exp * P7z) - (141750.0 * P8x * dP0y_exp * P1z) - (411075.0 * P6x * dP0y_exp * P3z) - (382725.0 * P4x * dP0y_exp * P5z) - (99225.0 * P2x * dP0y_exp * P7z) + (14175.0 * dP0y_exp * P9z);
    double dterm_2_dz = (10659600.0 * P2x * P6y * dP1z_exp) - (453600.0 * P8y * dP1z_exp) + (680400.0 * P6y * dP3z_exp) - (18002250.0 * P4x * P4y * dP1z_exp) - (17293500.0 * P2x * P4y * dP3z_exp) + (708750.0 * P4y * dP5z_exp) + (5202225.0 * P6x * P2y * dP1z_exp) + (9993375.0 * P4x * P2y * dP3z_exp) + (4380075.0 * P2x * P2y * dP5z_exp) - (411075.0 * P2y * dP7z_exp) - (141750.0 * P8x * dP1z_exp) - (411075.0 * P6x * dP3z_exp) - (382725.0 * P4x * dP5z_exp) - (99225.0 * P2x * dP7z_exp) + (14175.0 * dP9z_exp);

    double dterm_3_dx = (10659600.0 * dP6x_exp * P1y * P2z) - (453600.0 * dP8x_exp * P1y) + (680400.0 * dP6x_exp * P3y) - (18002250.0 * dP4x_exp * P1y * P4z) - (17293500.0 * dP4x_exp * P3y * P2z) + (708750.0 * dP4x_exp * P5y) + (5202225.0 * dP2x_exp * P1y * P6z) + (9993375.0 * dP2x_exp * P3y * P4z) + (4380075.0 * dP2x_exp * P5y * P2z) - (411075.0 * dP2x_exp * P7y) - (141750.0 * dP0x_exp * P1y * P8z) - (411075.0 * dP0x_exp * P3y * P6z) - (382725.0 * dP0x_exp * P5y * P4z) - (99225.0 * dP0x_exp * P7y * P2z) + (14175.0 * dP0x_exp * P9y);
    double dterm_3_dy = (10659600.0 * P6x * dP1y_exp * P2z) - (453600.0 * P8x * dP1y_exp) + (680400.0 * P6x * dP3y_exp) - (18002250.0 * P4x * dP1y_exp * P4z) - (17293500.0 * P4x * dP3y_exp * P2z) + (708750.0 * P4x * dP5y_exp) + (5202225.0 * P2x * dP1y_exp * P6z) + (9993375.0 * P2x * dP3y_exp * P4z) + (4380075.0 * P2x * dP5y_exp * P2z) - (411075.0 * P2x * dP7y_exp) - (141750.0 * dP1y_exp * P8z) - (411075.0 * dP3y_exp * P6z) - (382725.0 * dP5y_exp * P4z) - (99225.0 * dP7y_exp * P2z) + (14175.0 * dP9y_exp);
    double dterm_3_dz = (10659600.0 * P6x * P1y * dP2z_exp) - (453600.0 * P8x * P1y * dP0z_exp) + (680400.0 * P6x * P3y * dP0z_exp) - (18002250.0 * P4x * P1y * dP4z_exp) - (17293500.0 * P4x * P3y * dP2z_exp) + (708750.0 * P4x * P5y * dP0z_exp) + (5202225.0 * P2x * P1y * dP6z_exp) + (9993375.0 * P2x * P3y * dP4z_exp) + (4380075.0 * P2x * P5y * dP2z_exp) - (411075.0 * P2x * P7y * dP0z_exp) - (141750.0 * P1y * dP8z_exp) - (411075.0 * P3y * dP6z_exp) - (382725.0 * P5y * dP4z_exp) - (99225.0 * P7y * dP2z_exp) + (14175.0 * P9y * dP0z_exp);

    double dterm_4_dx = (10659600.0 * dP2x_exp * P1y * P6z) - (453600.0 * dP0x_exp * P1y * P8z) + (680400.0 * dP0x_exp * P3y * P6z) - (18002250.0 * dP4x_exp * P1y * P4z) - (17293500.0 * dP2x_exp * P3y * P4z) + (708750.0 * dP0x_exp * P5y * P4z) + (5202225.0 * dP6x_exp * P1y * P2z) + (9993375.0 * dP4x_exp * P3y * P2z) + (4380075.0 * dP2x_exp * P5y * P2z) - (411075.0 * dP0x_exp * P7y * P2z) - (141750.0 * dP8x_exp * P1y) - (411075.0 * dP6x_exp * P3y) - (382725.0 * dP4x_exp * P5y) - (99225.0 * dP2x_exp * P7y) + (14175.0 * dP0x_exp * P9y);
    double dterm_4_dy = (10659600.0 * P2x * dP1y_exp * P6z) - (453600.0 * dP1y_exp * P8z) + (680400.0 * dP3y_exp * P6z) - (18002250.0 * P4x * dP1y_exp * P4z) - (17293500.0 * P2x * dP3y_exp * P4z) + (708750.0 * dP5y_exp * P4z) + (5202225.0 * P6x * dP1y_exp * P2z) + (9993375.0 * P4x * dP3y_exp * P2z) + (4380075.0 * P2x * dP5y_exp * P2z) - (411075.0 * dP7y_exp * P2z) - (141750.0 * P8x * dP1y_exp) - (411075.0 * P6x * dP3y_exp) - (382725.0 * P4x * dP5y_exp) - (99225.0 * P2x * dP7y_exp) + (14175.0 * dP9y_exp);
    double dterm_4_dz = (10659600.0 * P2x * P1y * dP6z_exp) - (453600.0 * P1y * dP8z_exp) + (680400.0 * P3y * dP6z_exp) - (18002250.0 * P4x * P1y * dP4z_exp) - (17293500.0 * P2x * P3y * dP4z_exp) + (708750.0 * P5y * dP4z_exp) + (5202225.0 * P6x * P1y * dP2z_exp) + (9993375.0 * P4x * P3y * dP2z_exp) + (4380075.0 * P2x * P5y * dP2z_exp) - (411075.0 * P7y * dP2z_exp) - (141750.0 * P8x * P1y * dP0z_exp) - (411075.0 * P6x * P3y * dP0z_exp) - (382725.0 * P4x * P5y * dP0z_exp) - (99225.0 * P2x * P7y * dP0z_exp) + (14175.0 * P9y * dP0z_exp);

    double dterm_5_dx = (10659600.0 * dP1x_exp * P6y * P2z) - (453600.0 * dP1x_exp * P8y) + (680400.0 * dP3x_exp * P6y) - (18002250.0 * dP1x_exp * P4y * P4z) - (17293500.0 * dP3x_exp * P4y * P2z) + (708750.0 * dP5x_exp * P4y) + (5202225.0 * dP1x_exp * P2y * P6z) + (9993375.0 * dP3x_exp * P2y * P4z) + (4380075.0 * dP5x_exp * P2y * P2z) - (411075.0 * dP7x_exp * P2y) - (141750.0 * dP1x_exp * P8z) - (411075.0 * dP3x_exp * P6z) - (382725.0 * dP5x_exp * P4z) - (99225.0 * dP7x_exp * P2z) + (14175.0 * dP9x_exp);
    double dterm_5_dy = (10659600.0 * P1x * dP6y_exp * P2z) - (453600.0 * P1x * dP8y_exp) + (680400.0 * P3x * dP6y_exp) - (18002250.0 * P1x * dP4y_exp * P4z) - (17293500.0 * P3x * dP4y_exp * P2z) + (708750.0 * P5x * dP4y_exp) + (5202225.0 * P1x * dP2y_exp * P6z) + (9993375.0 * P3x * dP2y_exp * P4z) + (4380075.0 * P5x * dP2y_exp * P2z) - (411075.0 * P7x * dP2y_exp) - (141750.0 * P1x * dP0y_exp * P8z) - (411075.0 * P3x * dP0y_exp * P6z) - (382725.0 * P5x * dP0y_exp * P4z) - (99225.0 * P7x * dP0y_exp * P2z) + (14175.0 * P9x * dP0y_exp);
    double dterm_5_dz = (10659600.0 * P1x * P6y * dP2z_exp) - (453600.0 * P1x * P8y * dP0z_exp) + (680400.0 * P3x * P6y * dP0z_exp) - (18002250.0 * P1x * P4y * dP4z_exp) - (17293500.0 * P3x * P4y * dP2z_exp) + (708750.0 * P5x * P4y * dP0z_exp) + (5202225.0 * P1x * P2y * dP6z_exp) + (9993375.0 * P3x * P2y * dP4z_exp) + (4380075.0 * P5x * P2y * dP2z_exp) - (411075.0 * P7x * P2y * dP0z_exp) - (141750.0 * P1x * dP8z_exp) - (411075.0 * P3x * dP6z_exp) - (382725.0 * P5x * dP4z_exp) - (99225.0 * P7x * dP2z_exp) + (14175.0 * P9x * dP0z_exp);

    double dterm_6_dx = (10659600.0 * dP1x_exp * P2y * P6z) - (453600.0 * dP1x_exp * P8z) + (680400.0 * dP3x_exp * P6z) - (18002250.0 * dP1x_exp * P4y * P4z) - (17293500.0 * dP3x_exp * P2y * P4z) + (708750.0 * dP5x_exp * P4z) + (5202225.0 * dP1x_exp * P6y * P2z) + (9993375.0 * dP3x_exp * P4y * P2z) + (4380075.0 * dP5x_exp * P2y * P2z) - (411075.0 * dP7x_exp * P2z) - (141750.0 * dP1x_exp * P8y) - (411075.0 * dP3x_exp * P6y) - (382725.0 * dP5x_exp * P4y) - (99225.0 * dP7x_exp * P2y) + (14175.0 * dP9x_exp);
    double dterm_6_dy = (10659600.0 * P1x * dP2y_exp * P6z) - (453600.0 * P1x * dP0y_exp * P8z) + (680400.0 * P3x * dP0y_exp * P6z) - (18002250.0 * P1x * dP4y_exp * P4z) - (17293500.0 * P3x * dP2y_exp * P4z) + (708750.0 * P5x * dP0y_exp * P4z) + (5202225.0 * P1x * dP6y_exp * P2z) + (9993375.0 * P3x * dP4y_exp * P2z) + (4380075.0 * P5x * dP2y_exp * P2z) - (411075.0 * P7x * dP0y_exp * P2z) - (141750.0 * P1x * dP8y_exp) - (411075.0 * P3x * dP6y_exp) - (382725.0 * P5x * dP4y_exp) - (99225.0 * P7x * dP2y_exp) + (14175.0 * P9x * dP0y_exp);
    double dterm_6_dz = (10659600.0 * P1x * P2y * dP6z_exp) - (453600.0 * P1x * dP8z_exp) + (680400.0 * P3x * dP6z_exp) - (18002250.0 * P1x * P4y * dP4z_exp) - (17293500.0 * P3x * P2y * dP4z_exp) + (708750.0 * P5x * dP4z_exp) + (5202225.0 * P1x * P6y * dP2z_exp) + (9993375.0 * P3x * P4y * dP2z_exp) + (4380075.0 * P5x * P2y * dP2z_exp) - (411075.0 * P7x * dP2z_exp) - (141750.0 * P1x * P8y * dP0z_exp) - (411075.0 * P3x * P6y * dP0z_exp) - (382725.0 * P5x * P4y * dP0z_exp) - (99225.0 * P7x * P2y * dP0z_exp) + (14175.0 * P9x * dP0z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

    // dmiu4 dx/dy/dz
    deriv[9] = temp * dterm_4_dx;
    deriv[10] = temp * dterm_4_dy;
    deriv[11] = temp * dterm_4_dz;

    // dmiu5 dx/dy/dz
    deriv[12] = temp * dterm_5_dx;
    deriv[13] = temp * dterm_5_dy;
    deriv[14] = temp * dterm_5_dz;

    // dmiu6 dx/dy/dz
    deriv[15] = temp * dterm_6_dx;
    deriv[16] = temp * dterm_6_dy;
    deriv[17] = temp * dterm_6_dz;

}

void calc_solid_MCSH_9_7(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double P8x = P8(lambda, x0, gamma);
    double P8y = P8(lambda, y0, gamma);
    double P8z = P8(lambda, z0, gamma);

    double P9x = P9(lambda, x0, gamma);
    double P9y = P9(lambda, y0, gamma);
    double P9z = P9(lambda, z0, gamma);

    double term_1 = (17188605.0 * P5x * P4y) - (4944240.0 * P7x * P2y) + (697410.0 * P5x * P2y * P2z) + (136080.0 * P9x) + (45360.0 * P7x * P2z) - (274995.0 * P5x * P4z) - (11056500.0 * P3x * P6y) - (6038550.0 * P3x * P4y * P2z) + (4876200.0 * P3x * P2y * P4z) - (141750.0 * P3x * P6z) + (1134000.0 * P1x * P8y) + (1417500.0 * P1x * P6y * P2z) - (524475.0 * P1x * P4y * P4z) - (765450.0 * P1x * P2y * P6z) + (42525.0 * P1x * P8z);
    double term_2 = (17188605.0 * P4x * P5y) - (4944240.0 * P2x * P7y) + (697410.0 * P2x * P5y * P2z) + (136080.0 * P9y) + (45360.0 * P7y * P2z) - (274995.0 * P5y * P4z) - (11056500.0 * P6x * P3y) - (6038550.0 * P4x * P3y * P2z) + (4876200.0 * P2x * P3y * P4z) - (141750.0 * P3y * P6z) + (1134000.0 * P8x * P1y) + (1417500.0 * P6x * P1y * P2z) - (524475.0 * P4x * P1y * P4z) - (765450.0 * P2x * P1y * P6z) + (42525.0 * P1y * P8z);
    double term_3 = (17188605.0 * P5x * P4z) - (4944240.0 * P7x * P2z) + (697410.0 * P5x * P2y * P2z) + (136080.0 * P9x) + (45360.0 * P7x * P2y) - (274995.0 * P5x * P4y) - (11056500.0 * P3x * P6z) - (6038550.0 * P3x * P2y * P4z) + (4876200.0 * P3x * P4y * P2z) - (141750.0 * P3x * P6y) + (1134000.0 * P1x * P8z) + (1417500.0 * P1x * P2y * P6z) - (524475.0 * P1x * P4y * P4z) - (765450.0 * P1x * P6y * P2z) + (42525.0 * P1x * P8y);
    double term_4 = (17188605.0 * P4x * P5z) - (4944240.0 * P2x * P7z) + (697410.0 * P2x * P2y * P5z) + (136080.0 * P9z) + (45360.0 * P2y * P7z) - (274995.0 * P4y * P5z) - (11056500.0 * P6x * P3z) - (6038550.0 * P4x * P2y * P3z) + (4876200.0 * P2x * P4y * P3z) - (141750.0 * P6y * P3z) + (1134000.0 * P8x * P1z) + (1417500.0 * P6x * P2y * P1z) - (524475.0 * P4x * P4y * P1z) - (765450.0 * P2x * P6y * P1z) + (42525.0 * P8y * P1z);
    double term_5 = (17188605.0 * P5y * P4z) - (4944240.0 * P7y * P2z) + (697410.0 * P2x * P5y * P2z) + (136080.0 * P9y) + (45360.0 * P2x * P7y) - (274995.0 * P4x * P5y) - (11056500.0 * P3y * P6z) - (6038550.0 * P2x * P3y * P4z) + (4876200.0 * P4x * P3y * P2z) - (141750.0 * P6x * P3y) + (1134000.0 * P1y * P8z) + (1417500.0 * P2x * P1y * P6z) - (524475.0 * P4x * P1y * P4z) - (765450.0 * P6x * P1y * P2z) + (42525.0 * P8x * P1y);
    double term_6 = (17188605.0 * P4y * P5z) - (4944240.0 * P2y * P7z) + (697410.0 * P2x * P2y * P5z) + (136080.0 * P9z) + (45360.0 * P2x * P7z) - (274995.0 * P4x * P5z) - (11056500.0 * P6y * P3z) - (6038550.0 * P2x * P4y * P3z) + (4876200.0 * P4x * P2y * P3z) - (141750.0 * P6x * P3z) + (1134000.0 * P8y * P1z) + (1417500.0 * P2x * P6y * P1z) - (524475.0 * P4x * P4y * P1z) - (765450.0 * P6x * P2y * P1z) + (42525.0 * P8x * P1z);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP3x_exp = dP3_exp(P3x, C2, lambda, x0, gamma);
    double dP3y_exp = dP3_exp(P3y, C2, lambda, y0, gamma);
    double dP3z_exp = dP3_exp(P3z, C2, lambda, z0, gamma);

    double dP4x_exp = dP4_exp(P4x, C2, lambda, x0, gamma);
    double dP4y_exp = dP4_exp(P4y, C2, lambda, y0, gamma);
    double dP4z_exp = dP4_exp(P4z, C2, lambda, z0, gamma);

    double dP5x_exp = dP5_exp(P5x, C2, lambda, x0, gamma);
    double dP5y_exp = dP5_exp(P5y, C2, lambda, y0, gamma);
    double dP5z_exp = dP5_exp(P5z, C2, lambda, z0, gamma);

    double dP6x_exp = dP6_exp(P6x, C2, lambda, x0, gamma);
    double dP6y_exp = dP6_exp(P6y, C2, lambda, y0, gamma);
    double dP6z_exp = dP6_exp(P6z, C2, lambda, z0, gamma);

    double dP7x_exp = dP7_exp(P7x, C2, lambda, x0, gamma);
    double dP7y_exp = dP7_exp(P7y, C2, lambda, y0, gamma);
    double dP7z_exp = dP7_exp(P7z, C2, lambda, z0, gamma);

    double dP8x_exp = dP8_exp(P8x, C2, lambda, x0, gamma);
    double dP8y_exp = dP8_exp(P8y, C2, lambda, y0, gamma);
    double dP8z_exp = dP8_exp(P8z, C2, lambda, z0, gamma);

    double dP9x_exp = dP9_exp(P9x, C2, lambda, x0, gamma);
    double dP9y_exp = dP9_exp(P9y, C2, lambda, y0, gamma);
    double dP9z_exp = dP9_exp(P9z, C2, lambda, z0, gamma);

    double dterm_1_dx = (17188605.0 * dP5x_exp * P4y) - (4944240.0 * dP7x_exp * P2y) + (697410.0 * dP5x_exp * P2y * P2z) + (136080.0 * dP9x_exp) + (45360.0 * dP7x_exp * P2z) - (274995.0 * dP5x_exp * P4z) - (11056500.0 * dP3x_exp * P6y) - (6038550.0 * dP3x_exp * P4y * P2z) + (4876200.0 * dP3x_exp * P2y * P4z) - (141750.0 * dP3x_exp * P6z) + (1134000.0 * dP1x_exp * P8y) + (1417500.0 * dP1x_exp * P6y * P2z) - (524475.0 * dP1x_exp * P4y * P4z) - (765450.0 * dP1x_exp * P2y * P6z) + (42525.0 * dP1x_exp * P8z);
    double dterm_1_dy = (17188605.0 * P5x * dP4y_exp) - (4944240.0 * P7x * dP2y_exp) + (697410.0 * P5x * dP2y_exp * P2z) + (136080.0 * P9x * dP0y_exp) + (45360.0 * P7x * dP0y_exp * P2z) - (274995.0 * P5x * dP0y_exp * P4z) - (11056500.0 * P3x * dP6y_exp) - (6038550.0 * P3x * dP4y_exp * P2z) + (4876200.0 * P3x * dP2y_exp * P4z) - (141750.0 * P3x * dP0y_exp * P6z) + (1134000.0 * P1x * dP8y_exp) + (1417500.0 * P1x * dP6y_exp * P2z) - (524475.0 * P1x * dP4y_exp * P4z) - (765450.0 * P1x * dP2y_exp * P6z) + (42525.0 * P1x * dP0y_exp * P8z);
    double dterm_1_dz = (17188605.0 * P5x * P4y * dP0z_exp) - (4944240.0 * P7x * P2y * dP0z_exp) + (697410.0 * P5x * P2y * dP2z_exp) + (136080.0 * P9x * dP0z_exp) + (45360.0 * P7x * dP2z_exp) - (274995.0 * P5x * dP4z_exp) - (11056500.0 * P3x * P6y * dP0z_exp) - (6038550.0 * P3x * P4y * dP2z_exp) + (4876200.0 * P3x * P2y * dP4z_exp) - (141750.0 * P3x * dP6z_exp) + (1134000.0 * P1x * P8y * dP0z_exp) + (1417500.0 * P1x * P6y * dP2z_exp) - (524475.0 * P1x * P4y * dP4z_exp) - (765450.0 * P1x * P2y * dP6z_exp) + (42525.0 * P1x * dP8z_exp);

    double dterm_2_dx = (17188605.0 * dP4x_exp * P5y) - (4944240.0 * dP2x_exp * P7y) + (697410.0 * dP2x_exp * P5y * P2z) + (136080.0 * dP0x_exp * P9y) + (45360.0 * dP0x_exp * P7y * P2z) - (274995.0 * dP0x_exp * P5y * P4z) - (11056500.0 * dP6x_exp * P3y) - (6038550.0 * dP4x_exp * P3y * P2z) + (4876200.0 * dP2x_exp * P3y * P4z) - (141750.0 * dP0x_exp * P3y * P6z) + (1134000.0 * dP8x_exp * P1y) + (1417500.0 * dP6x_exp * P1y * P2z) - (524475.0 * dP4x_exp * P1y * P4z) - (765450.0 * dP2x_exp * P1y * P6z) + (42525.0 * dP0x_exp * P1y * P8z);
    double dterm_2_dy = (17188605.0 * P4x * dP5y_exp) - (4944240.0 * P2x * dP7y_exp) + (697410.0 * P2x * dP5y_exp * P2z) + (136080.0 * dP9y_exp) + (45360.0 * dP7y_exp * P2z) - (274995.0 * dP5y_exp * P4z) - (11056500.0 * P6x * dP3y_exp) - (6038550.0 * P4x * dP3y_exp * P2z) + (4876200.0 * P2x * dP3y_exp * P4z) - (141750.0 * dP3y_exp * P6z) + (1134000.0 * P8x * dP1y_exp) + (1417500.0 * P6x * dP1y_exp * P2z) - (524475.0 * P4x * dP1y_exp * P4z) - (765450.0 * P2x * dP1y_exp * P6z) + (42525.0 * dP1y_exp * P8z);
    double dterm_2_dz = (17188605.0 * P4x * P5y * dP0z_exp) - (4944240.0 * P2x * P7y * dP0z_exp) + (697410.0 * P2x * P5y * dP2z_exp) + (136080.0 * P9y * dP0z_exp) + (45360.0 * P7y * dP2z_exp) - (274995.0 * P5y * dP4z_exp) - (11056500.0 * P6x * P3y * dP0z_exp) - (6038550.0 * P4x * P3y * dP2z_exp) + (4876200.0 * P2x * P3y * dP4z_exp) - (141750.0 * P3y * dP6z_exp) + (1134000.0 * P8x * P1y * dP0z_exp) + (1417500.0 * P6x * P1y * dP2z_exp) - (524475.0 * P4x * P1y * dP4z_exp) - (765450.0 * P2x * P1y * dP6z_exp) + (42525.0 * P1y * dP8z_exp);

    double dterm_3_dx = (17188605.0 * dP5x_exp * P4z) - (4944240.0 * dP7x_exp * P2z) + (697410.0 * dP5x_exp * P2y * P2z) + (136080.0 * dP9x_exp) + (45360.0 * dP7x_exp * P2y) - (274995.0 * dP5x_exp * P4y) - (11056500.0 * dP3x_exp * P6z) - (6038550.0 * dP3x_exp * P2y * P4z) + (4876200.0 * dP3x_exp * P4y * P2z) - (141750.0 * dP3x_exp * P6y) + (1134000.0 * dP1x_exp * P8z) + (1417500.0 * dP1x_exp * P2y * P6z) - (524475.0 * dP1x_exp * P4y * P4z) - (765450.0 * dP1x_exp * P6y * P2z) + (42525.0 * dP1x_exp * P8y);
    double dterm_3_dy = (17188605.0 * P5x * dP0y_exp * P4z) - (4944240.0 * P7x * dP0y_exp * P2z) + (697410.0 * P5x * dP2y_exp * P2z) + (136080.0 * P9x * dP0y_exp) + (45360.0 * P7x * dP2y_exp) - (274995.0 * P5x * dP4y_exp) - (11056500.0 * P3x * dP0y_exp * P6z) - (6038550.0 * P3x * dP2y_exp * P4z) + (4876200.0 * P3x * dP4y_exp * P2z) - (141750.0 * P3x * dP6y_exp) + (1134000.0 * P1x * dP0y_exp * P8z) + (1417500.0 * P1x * dP2y_exp * P6z) - (524475.0 * P1x * dP4y_exp * P4z) - (765450.0 * P1x * dP6y_exp * P2z) + (42525.0 * P1x * dP8y_exp);
    double dterm_3_dz = (17188605.0 * P5x * dP4z_exp) - (4944240.0 * P7x * dP2z_exp) + (697410.0 * P5x * P2y * dP2z_exp) + (136080.0 * P9x * dP0z_exp) + (45360.0 * P7x * P2y * dP0z_exp) - (274995.0 * P5x * P4y * dP0z_exp) - (11056500.0 * P3x * dP6z_exp) - (6038550.0 * P3x * P2y * dP4z_exp) + (4876200.0 * P3x * P4y * dP2z_exp) - (141750.0 * P3x * P6y * dP0z_exp) + (1134000.0 * P1x * dP8z_exp) + (1417500.0 * P1x * P2y * dP6z_exp) - (524475.0 * P1x * P4y * dP4z_exp) - (765450.0 * P1x * P6y * dP2z_exp) + (42525.0 * P1x * P8y * dP0z_exp);

    double dterm_4_dx = (17188605.0 * dP4x_exp * P5z) - (4944240.0 * dP2x_exp * P7z) + (697410.0 * dP2x_exp * P2y * P5z) + (136080.0 * dP0x_exp * P9z) + (45360.0 * dP0x_exp * P2y * P7z) - (274995.0 * dP0x_exp * P4y * P5z) - (11056500.0 * dP6x_exp * P3z) - (6038550.0 * dP4x_exp * P2y * P3z) + (4876200.0 * dP2x_exp * P4y * P3z) - (141750.0 * dP0x_exp * P6y * P3z) + (1134000.0 * dP8x_exp * P1z) + (1417500.0 * dP6x_exp * P2y * P1z) - (524475.0 * dP4x_exp * P4y * P1z) - (765450.0 * dP2x_exp * P6y * P1z) + (42525.0 * dP0x_exp * P8y * P1z);
    double dterm_4_dy = (17188605.0 * P4x * dP0y_exp * P5z) - (4944240.0 * P2x * dP0y_exp * P7z) + (697410.0 * P2x * dP2y_exp * P5z) + (136080.0 * dP0y_exp * P9z) + (45360.0 * dP2y_exp * P7z) - (274995.0 * dP4y_exp * P5z) - (11056500.0 * P6x * dP0y_exp * P3z) - (6038550.0 * P4x * dP2y_exp * P3z) + (4876200.0 * P2x * dP4y_exp * P3z) - (141750.0 * dP6y_exp * P3z) + (1134000.0 * P8x * dP0y_exp * P1z) + (1417500.0 * P6x * dP2y_exp * P1z) - (524475.0 * P4x * dP4y_exp * P1z) - (765450.0 * P2x * dP6y_exp * P1z) + (42525.0 * dP8y_exp * P1z);
    double dterm_4_dz = (17188605.0 * P4x * dP5z_exp) - (4944240.0 * P2x * dP7z_exp) + (697410.0 * P2x * P2y * dP5z_exp) + (136080.0 * dP9z_exp) + (45360.0 * P2y * dP7z_exp) - (274995.0 * P4y * dP5z_exp) - (11056500.0 * P6x * dP3z_exp) - (6038550.0 * P4x * P2y * dP3z_exp) + (4876200.0 * P2x * P4y * dP3z_exp) - (141750.0 * P6y * dP3z_exp) + (1134000.0 * P8x * dP1z_exp) + (1417500.0 * P6x * P2y * dP1z_exp) - (524475.0 * P4x * P4y * dP1z_exp) - (765450.0 * P2x * P6y * dP1z_exp) + (42525.0 * P8y * dP1z_exp);

    double dterm_5_dx = (17188605.0 * dP0x_exp * P5y * P4z) - (4944240.0 * dP0x_exp * P7y * P2z) + (697410.0 * dP2x_exp * P5y * P2z) + (136080.0 * dP0x_exp * P9y) + (45360.0 * dP2x_exp * P7y) - (274995.0 * dP4x_exp * P5y) - (11056500.0 * dP0x_exp * P3y * P6z) - (6038550.0 * dP2x_exp * P3y * P4z) + (4876200.0 * dP4x_exp * P3y * P2z) - (141750.0 * dP6x_exp * P3y) + (1134000.0 * dP0x_exp * P1y * P8z) + (1417500.0 * dP2x_exp * P1y * P6z) - (524475.0 * dP4x_exp * P1y * P4z) - (765450.0 * dP6x_exp * P1y * P2z) + (42525.0 * dP8x_exp * P1y);
    double dterm_5_dy = (17188605.0 * dP5y_exp * P4z) - (4944240.0 * dP7y_exp * P2z) + (697410.0 * P2x * dP5y_exp * P2z) + (136080.0 * dP9y_exp) + (45360.0 * P2x * dP7y_exp) - (274995.0 * P4x * dP5y_exp) - (11056500.0 * dP3y_exp * P6z) - (6038550.0 * P2x * dP3y_exp * P4z) + (4876200.0 * P4x * dP3y_exp * P2z) - (141750.0 * P6x * dP3y_exp) + (1134000.0 * dP1y_exp * P8z) + (1417500.0 * P2x * dP1y_exp * P6z) - (524475.0 * P4x * dP1y_exp * P4z) - (765450.0 * P6x * dP1y_exp * P2z) + (42525.0 * P8x * dP1y_exp);
    double dterm_5_dz = (17188605.0 * P5y * dP4z_exp) - (4944240.0 * P7y * dP2z_exp) + (697410.0 * P2x * P5y * dP2z_exp) + (136080.0 * P9y * dP0z_exp) + (45360.0 * P2x * P7y * dP0z_exp) - (274995.0 * P4x * P5y * dP0z_exp) - (11056500.0 * P3y * dP6z_exp) - (6038550.0 * P2x * P3y * dP4z_exp) + (4876200.0 * P4x * P3y * dP2z_exp) - (141750.0 * P6x * P3y * dP0z_exp) + (1134000.0 * P1y * dP8z_exp) + (1417500.0 * P2x * P1y * dP6z_exp) - (524475.0 * P4x * P1y * dP4z_exp) - (765450.0 * P6x * P1y * dP2z_exp) + (42525.0 * P8x * P1y * dP0z_exp);

    double dterm_6_dx = (17188605.0 * dP0x_exp * P4y * P5z) - (4944240.0 * dP0x_exp * P2y * P7z) + (697410.0 * dP2x_exp * P2y * P5z) + (136080.0 * dP0x_exp * P9z) + (45360.0 * dP2x_exp * P7z) - (274995.0 * dP4x_exp * P5z) - (11056500.0 * dP0x_exp * P6y * P3z) - (6038550.0 * dP2x_exp * P4y * P3z) + (4876200.0 * dP4x_exp * P2y * P3z) - (141750.0 * dP6x_exp * P3z) + (1134000.0 * dP0x_exp * P8y * P1z) + (1417500.0 * dP2x_exp * P6y * P1z) - (524475.0 * dP4x_exp * P4y * P1z) - (765450.0 * dP6x_exp * P2y * P1z) + (42525.0 * dP8x_exp * P1z);
    double dterm_6_dy = (17188605.0 * dP4y_exp * P5z) - (4944240.0 * dP2y_exp * P7z) + (697410.0 * P2x * dP2y_exp * P5z) + (136080.0 * dP0y_exp * P9z) + (45360.0 * P2x * dP0y_exp * P7z) - (274995.0 * P4x * dP0y_exp * P5z) - (11056500.0 * dP6y_exp * P3z) - (6038550.0 * P2x * dP4y_exp * P3z) + (4876200.0 * P4x * dP2y_exp * P3z) - (141750.0 * P6x * dP0y_exp * P3z) + (1134000.0 * dP8y_exp * P1z) + (1417500.0 * P2x * dP6y_exp * P1z) - (524475.0 * P4x * dP4y_exp * P1z) - (765450.0 * P6x * dP2y_exp * P1z) + (42525.0 * P8x * dP0y_exp * P1z);
    double dterm_6_dz = (17188605.0 * P4y * dP5z_exp) - (4944240.0 * P2y * dP7z_exp) + (697410.0 * P2x * P2y * dP5z_exp) + (136080.0 * dP9z_exp) + (45360.0 * P2x * dP7z_exp) - (274995.0 * P4x * dP5z_exp) - (11056500.0 * P6y * dP3z_exp) - (6038550.0 * P2x * P4y * dP3z_exp) + (4876200.0 * P4x * P2y * dP3z_exp) - (141750.0 * P6x * dP3z_exp) + (1134000.0 * P8y * dP1z_exp) + (1417500.0 * P2x * P6y * dP1z_exp) - (524475.0 * P4x * P4y * dP1z_exp) - (765450.0 * P6x * P2y * dP1z_exp) + (42525.0 * P8x * dP1z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

    // dmiu4 dx/dy/dz
    deriv[9] = temp * dterm_4_dx;
    deriv[10] = temp * dterm_4_dy;
    deriv[11] = temp * dterm_4_dz;

    // dmiu5 dx/dy/dz
    deriv[12] = temp * dterm_5_dx;
    deriv[13] = temp * dterm_5_dy;
    deriv[14] = temp * dterm_5_dz;

    // dmiu6 dx/dy/dz
    deriv[15] = temp * dterm_6_dx;
    deriv[16] = temp * dterm_6_dy;
    deriv[17] = temp * dterm_6_dz;

}

void calc_solid_MCSH_9_8(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double term_1 = (16839900.0 * P5x * P3y * P1z) - (2494800.0 * P7x * P1y * P1z) + (623700.0 * P5x * P1y * P3z) - (13565475.0 * P3x * P5y * P1z) - (10914750.0 * P3x * P3y * P3z) + (2650725.0 * P3x * P1y * P5z) + (1559250.0 * P1x * P7y * P1z) + (2650725.0 * P1x * P5y * P3z) + (623700.0 * P1x * P3y * P5z) - (467775.0 * P1x * P1y * P7z);
    double term_2 = (16839900.0 * P3x * P5y * P1z) - (2494800.0 * P1x * P7y * P1z) + (623700.0 * P1x * P5y * P3z) - (13565475.0 * P5x * P3y * P1z) - (10914750.0 * P3x * P3y * P3z) + (2650725.0 * P1x * P3y * P5z) + (1559250.0 * P7x * P1y * P1z) + (2650725.0 * P5x * P1y * P3z) + (623700.0 * P3x * P1y * P5z) - (467775.0 * P1x * P1y * P7z);
    double term_3 = (16839900.0 * P5x * P1y * P3z) - (2494800.0 * P7x * P1y * P1z) + (623700.0 * P5x * P3y * P1z) - (13565475.0 * P3x * P1y * P5z) - (10914750.0 * P3x * P3y * P3z) + (2650725.0 * P3x * P5y * P1z) + (1559250.0 * P1x * P1y * P7z) + (2650725.0 * P1x * P3y * P5z) + (623700.0 * P1x * P5y * P3z) - (467775.0 * P1x * P7y * P1z);
    double term_4 = (16839900.0 * P3x * P1y * P5z) - (2494800.0 * P1x * P1y * P7z) + (623700.0 * P1x * P3y * P5z) - (13565475.0 * P5x * P1y * P3z) - (10914750.0 * P3x * P3y * P3z) + (2650725.0 * P1x * P5y * P3z) + (1559250.0 * P7x * P1y * P1z) + (2650725.0 * P5x * P3y * P1z) + (623700.0 * P3x * P5y * P1z) - (467775.0 * P1x * P7y * P1z);
    double term_5 = (16839900.0 * P1x * P5y * P3z) - (2494800.0 * P1x * P7y * P1z) + (623700.0 * P3x * P5y * P1z) - (13565475.0 * P1x * P3y * P5z) - (10914750.0 * P3x * P3y * P3z) + (2650725.0 * P5x * P3y * P1z) + (1559250.0 * P1x * P1y * P7z) + (2650725.0 * P3x * P1y * P5z) + (623700.0 * P5x * P1y * P3z) - (467775.0 * P7x * P1y * P1z);
    double term_6 = (16839900.0 * P1x * P3y * P5z) - (2494800.0 * P1x * P1y * P7z) + (623700.0 * P3x * P1y * P5z) - (13565475.0 * P1x * P5y * P3z) - (10914750.0 * P3x * P3y * P3z) + (2650725.0 * P5x * P1y * P3z) + (1559250.0 * P1x * P7y * P1z) + (2650725.0 * P3x * P5y * P1z) + (623700.0 * P5x * P3y * P1z) - (467775.0 * P7x * P1y * P1z);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dP3x_exp = dP3_exp(P3x, C2, lambda, x0, gamma);
    double dP3y_exp = dP3_exp(P3y, C2, lambda, y0, gamma);
    double dP3z_exp = dP3_exp(P3z, C2, lambda, z0, gamma);

    double dP5x_exp = dP5_exp(P5x, C2, lambda, x0, gamma);
    double dP5y_exp = dP5_exp(P5y, C2, lambda, y0, gamma);
    double dP5z_exp = dP5_exp(P5z, C2, lambda, z0, gamma);

    double dP7x_exp = dP7_exp(P7x, C2, lambda, x0, gamma);
    double dP7y_exp = dP7_exp(P7y, C2, lambda, y0, gamma);
    double dP7z_exp = dP7_exp(P7z, C2, lambda, z0, gamma);

    double dterm_1_dx = (16839900.0 * dP5x_exp * P3y * P1z) - (2494800.0 * dP7x_exp * P1y * P1z) + (623700.0 * dP5x_exp * P1y * P3z) - (13565475.0 * dP3x_exp * P5y * P1z) - (10914750.0 * dP3x_exp * P3y * P3z) + (2650725.0 * dP3x_exp * P1y * P5z) + (1559250.0 * dP1x_exp * P7y * P1z) + (2650725.0 * dP1x_exp * P5y * P3z) + (623700.0 * dP1x_exp * P3y * P5z) - (467775.0 * dP1x_exp * P1y * P7z);
    double dterm_1_dy = (16839900.0 * P5x * dP3y_exp * P1z) - (2494800.0 * P7x * dP1y_exp * P1z) + (623700.0 * P5x * dP1y_exp * P3z) - (13565475.0 * P3x * dP5y_exp * P1z) - (10914750.0 * P3x * dP3y_exp * P3z) + (2650725.0 * P3x * dP1y_exp * P5z) + (1559250.0 * P1x * dP7y_exp * P1z) + (2650725.0 * P1x * dP5y_exp * P3z) + (623700.0 * P1x * dP3y_exp * P5z) - (467775.0 * P1x * dP1y_exp * P7z);
    double dterm_1_dz = (16839900.0 * P5x * P3y * dP1z_exp) - (2494800.0 * P7x * P1y * dP1z_exp) + (623700.0 * P5x * P1y * dP3z_exp) - (13565475.0 * P3x * P5y * dP1z_exp) - (10914750.0 * P3x * P3y * dP3z_exp) + (2650725.0 * P3x * P1y * dP5z_exp) + (1559250.0 * P1x * P7y * dP1z_exp) + (2650725.0 * P1x * P5y * dP3z_exp) + (623700.0 * P1x * P3y * dP5z_exp) - (467775.0 * P1x * P1y * dP7z_exp);

    double dterm_2_dx = (16839900.0 * dP3x_exp * P5y * P1z) - (2494800.0 * dP1x_exp * P7y * P1z) + (623700.0 * dP1x_exp * P5y * P3z) - (13565475.0 * dP5x_exp * P3y * P1z) - (10914750.0 * dP3x_exp * P3y * P3z) + (2650725.0 * dP1x_exp * P3y * P5z) + (1559250.0 * dP7x_exp * P1y * P1z) + (2650725.0 * dP5x_exp * P1y * P3z) + (623700.0 * dP3x_exp * P1y * P5z) - (467775.0 * dP1x_exp * P1y * P7z);
    double dterm_2_dy = (16839900.0 * P3x * dP5y_exp * P1z) - (2494800.0 * P1x * dP7y_exp * P1z) + (623700.0 * P1x * dP5y_exp * P3z) - (13565475.0 * P5x * dP3y_exp * P1z) - (10914750.0 * P3x * dP3y_exp * P3z) + (2650725.0 * P1x * dP3y_exp * P5z) + (1559250.0 * P7x * dP1y_exp * P1z) + (2650725.0 * P5x * dP1y_exp * P3z) + (623700.0 * P3x * dP1y_exp * P5z) - (467775.0 * P1x * dP1y_exp * P7z);
    double dterm_2_dz = (16839900.0 * P3x * P5y * dP1z_exp) - (2494800.0 * P1x * P7y * dP1z_exp) + (623700.0 * P1x * P5y * dP3z_exp) - (13565475.0 * P5x * P3y * dP1z_exp) - (10914750.0 * P3x * P3y * dP3z_exp) + (2650725.0 * P1x * P3y * dP5z_exp) + (1559250.0 * P7x * P1y * dP1z_exp) + (2650725.0 * P5x * P1y * dP3z_exp) + (623700.0 * P3x * P1y * dP5z_exp) - (467775.0 * P1x * P1y * dP7z_exp);

    double dterm_3_dx = (16839900.0 * dP5x_exp * P1y * P3z) - (2494800.0 * dP7x_exp * P1y * P1z) + (623700.0 * dP5x_exp * P3y * P1z) - (13565475.0 * dP3x_exp * P1y * P5z) - (10914750.0 * dP3x_exp * P3y * P3z) + (2650725.0 * dP3x_exp * P5y * P1z) + (1559250.0 * dP1x_exp * P1y * P7z) + (2650725.0 * dP1x_exp * P3y * P5z) + (623700.0 * dP1x_exp * P5y * P3z) - (467775.0 * dP1x_exp * P7y * P1z);
    double dterm_3_dy = (16839900.0 * P5x * dP1y_exp * P3z) - (2494800.0 * P7x * dP1y_exp * P1z) + (623700.0 * P5x * dP3y_exp * P1z) - (13565475.0 * P3x * dP1y_exp * P5z) - (10914750.0 * P3x * dP3y_exp * P3z) + (2650725.0 * P3x * dP5y_exp * P1z) + (1559250.0 * P1x * dP1y_exp * P7z) + (2650725.0 * P1x * dP3y_exp * P5z) + (623700.0 * P1x * dP5y_exp * P3z) - (467775.0 * P1x * dP7y_exp * P1z);
    double dterm_3_dz = (16839900.0 * P5x * P1y * dP3z_exp) - (2494800.0 * P7x * P1y * dP1z_exp) + (623700.0 * P5x * P3y * dP1z_exp) - (13565475.0 * P3x * P1y * dP5z_exp) - (10914750.0 * P3x * P3y * dP3z_exp) + (2650725.0 * P3x * P5y * dP1z_exp) + (1559250.0 * P1x * P1y * dP7z_exp) + (2650725.0 * P1x * P3y * dP5z_exp) + (623700.0 * P1x * P5y * dP3z_exp) - (467775.0 * P1x * P7y * dP1z_exp);

    double dterm_4_dx = (16839900.0 * dP3x_exp * P1y * P5z) - (2494800.0 * dP1x_exp * P1y * P7z) + (623700.0 * dP1x_exp * P3y * P5z) - (13565475.0 * dP5x_exp * P1y * P3z) - (10914750.0 * dP3x_exp * P3y * P3z) + (2650725.0 * dP1x_exp * P5y * P3z) + (1559250.0 * dP7x_exp * P1y * P1z) + (2650725.0 * dP5x_exp * P3y * P1z) + (623700.0 * dP3x_exp * P5y * P1z) - (467775.0 * dP1x_exp * P7y * P1z);
    double dterm_4_dy = (16839900.0 * P3x * dP1y_exp * P5z) - (2494800.0 * P1x * dP1y_exp * P7z) + (623700.0 * P1x * dP3y_exp * P5z) - (13565475.0 * P5x * dP1y_exp * P3z) - (10914750.0 * P3x * dP3y_exp * P3z) + (2650725.0 * P1x * dP5y_exp * P3z) + (1559250.0 * P7x * dP1y_exp * P1z) + (2650725.0 * P5x * dP3y_exp * P1z) + (623700.0 * P3x * dP5y_exp * P1z) - (467775.0 * P1x * dP7y_exp * P1z);
    double dterm_4_dz = (16839900.0 * P3x * P1y * dP5z_exp) - (2494800.0 * P1x * P1y * dP7z_exp) + (623700.0 * P1x * P3y * dP5z_exp) - (13565475.0 * P5x * P1y * dP3z_exp) - (10914750.0 * P3x * P3y * dP3z_exp) + (2650725.0 * P1x * P5y * dP3z_exp) + (1559250.0 * P7x * P1y * dP1z_exp) + (2650725.0 * P5x * P3y * dP1z_exp) + (623700.0 * P3x * P5y * dP1z_exp) - (467775.0 * P1x * P7y * dP1z_exp);

    double dterm_5_dx = (16839900.0 * dP1x_exp * P5y * P3z) - (2494800.0 * dP1x_exp * P7y * P1z) + (623700.0 * dP3x_exp * P5y * P1z) - (13565475.0 * dP1x_exp * P3y * P5z) - (10914750.0 * dP3x_exp * P3y * P3z) + (2650725.0 * dP5x_exp * P3y * P1z) + (1559250.0 * dP1x_exp * P1y * P7z) + (2650725.0 * dP3x_exp * P1y * P5z) + (623700.0 * dP5x_exp * P1y * P3z) - (467775.0 * dP7x_exp * P1y * P1z);
    double dterm_5_dy = (16839900.0 * P1x * dP5y_exp * P3z) - (2494800.0 * P1x * dP7y_exp * P1z) + (623700.0 * P3x * dP5y_exp * P1z) - (13565475.0 * P1x * dP3y_exp * P5z) - (10914750.0 * P3x * dP3y_exp * P3z) + (2650725.0 * P5x * dP3y_exp * P1z) + (1559250.0 * P1x * dP1y_exp * P7z) + (2650725.0 * P3x * dP1y_exp * P5z) + (623700.0 * P5x * dP1y_exp * P3z) - (467775.0 * P7x * dP1y_exp * P1z);
    double dterm_5_dz = (16839900.0 * P1x * P5y * dP3z_exp) - (2494800.0 * P1x * P7y * dP1z_exp) + (623700.0 * P3x * P5y * dP1z_exp) - (13565475.0 * P1x * P3y * dP5z_exp) - (10914750.0 * P3x * P3y * dP3z_exp) + (2650725.0 * P5x * P3y * dP1z_exp) + (1559250.0 * P1x * P1y * dP7z_exp) + (2650725.0 * P3x * P1y * dP5z_exp) + (623700.0 * P5x * P1y * dP3z_exp) - (467775.0 * P7x * P1y * dP1z_exp);

    double dterm_6_dx = (16839900.0 * dP1x_exp * P3y * P5z) - (2494800.0 * dP1x_exp * P1y * P7z) + (623700.0 * dP3x_exp * P1y * P5z) - (13565475.0 * dP1x_exp * P5y * P3z) - (10914750.0 * dP3x_exp * P3y * P3z) + (2650725.0 * dP5x_exp * P1y * P3z) + (1559250.0 * dP1x_exp * P7y * P1z) + (2650725.0 * dP3x_exp * P5y * P1z) + (623700.0 * dP5x_exp * P3y * P1z) - (467775.0 * dP7x_exp * P1y * P1z);
    double dterm_6_dy = (16839900.0 * P1x * dP3y_exp * P5z) - (2494800.0 * P1x * dP1y_exp * P7z) + (623700.0 * P3x * dP1y_exp * P5z) - (13565475.0 * P1x * dP5y_exp * P3z) - (10914750.0 * P3x * dP3y_exp * P3z) + (2650725.0 * P5x * dP1y_exp * P3z) + (1559250.0 * P1x * dP7y_exp * P1z) + (2650725.0 * P3x * dP5y_exp * P1z) + (623700.0 * P5x * dP3y_exp * P1z) - (467775.0 * P7x * dP1y_exp * P1z);
    double dterm_6_dz = (16839900.0 * P1x * P3y * dP5z_exp) - (2494800.0 * P1x * P1y * dP7z_exp) + (623700.0 * P3x * P1y * dP5z_exp) - (13565475.0 * P1x * P5y * dP3z_exp) - (10914750.0 * P3x * P3y * dP3z_exp) + (2650725.0 * P5x * P1y * dP3z_exp) + (1559250.0 * P1x * P7y * dP1z_exp) + (2650725.0 * P3x * P5y * dP1z_exp) + (623700.0 * P5x * P3y * dP1z_exp) - (467775.0 * P7x * P1y * dP1z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

    // dmiu4 dx/dy/dz
    deriv[9] = temp * dterm_4_dx;
    deriv[10] = temp * dterm_4_dy;
    deriv[11] = temp * dterm_4_dz;

    // dmiu5 dx/dy/dz
    deriv[12] = temp * dterm_5_dx;
    deriv[13] = temp * dterm_5_dy;
    deriv[14] = temp * dterm_5_dz;

    // dmiu6 dx/dy/dz
    deriv[15] = temp * dterm_6_dx;
    deriv[16] = temp * dterm_6_dy;
    deriv[17] = temp * dterm_6_dz;

}

void calc_solid_MCSH_9_9(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double P8x = P8(lambda, x0, gamma);
    double P8y = P8(lambda, y0, gamma);
    double P8z = P8(lambda, z0, gamma);

    double P9x = P9(lambda, x0, gamma);
    double P9y = P9(lambda, y0, gamma);
    double P9z = P9(lambda, z0, gamma);

    double term_1 = (16448670.0 * P5x * P2y * P2z) - (816480.0 * P7x * P2y) + (116235.0 * P5x * P4y) - (816480.0 * P7x * P2z) + (116235.0 * P5x * P4z) + (45360.0 * P9x) - (13707225.0 * P3x * P4y * P2z) - (13707225.0 * P3x * P2y * P4z) + (836325.0 * P3x * P6y) + (836325.0 * P3x * P6z) + (1460025.0 * P1x * P6y * P2z) + (3203550.0 * P1x * P4y * P4z) + (1460025.0 * P1x * P2y * P6z) - (141750.0 * P1x * P8y) - (141750.0 * P1x * P8z);
    double term_2 = (16448670.0 * P2x * P5y * P2z) - (816480.0 * P2x * P7y) + (116235.0 * P4x * P5y) - (816480.0 * P7y * P2z) + (116235.0 * P5y * P4z) + (45360.0 * P9y) - (13707225.0 * P4x * P3y * P2z) - (13707225.0 * P2x * P3y * P4z) + (836325.0 * P6x * P3y) + (836325.0 * P3y * P6z) + (1460025.0 * P6x * P1y * P2z) + (3203550.0 * P4x * P1y * P4z) + (1460025.0 * P2x * P1y * P6z) - (141750.0 * P8x * P1y) - (141750.0 * P1y * P8z);
    double term_3 = (16448670.0 * P2x * P2y * P5z) - (816480.0 * P2x * P7z) + (116235.0 * P4x * P5z) - (816480.0 * P2y * P7z) + (116235.0 * P4y * P5z) + (45360.0 * P9z) - (13707225.0 * P4x * P2y * P3z) - (13707225.0 * P2x * P4y * P3z) + (836325.0 * P6x * P3z) + (836325.0 * P6y * P3z) + (1460025.0 * P6x * P2y * P1z) + (3203550.0 * P4x * P4y * P1z) + (1460025.0 * P2x * P6y * P1z) - (141750.0 * P8x * P1z) - (141750.0 * P8y * P1z);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP3x_exp = dP3_exp(P3x, C2, lambda, x0, gamma);
    double dP3y_exp = dP3_exp(P3y, C2, lambda, y0, gamma);
    double dP3z_exp = dP3_exp(P3z, C2, lambda, z0, gamma);

    double dP4x_exp = dP4_exp(P4x, C2, lambda, x0, gamma);
    double dP4y_exp = dP4_exp(P4y, C2, lambda, y0, gamma);
    double dP4z_exp = dP4_exp(P4z, C2, lambda, z0, gamma);

    double dP5x_exp = dP5_exp(P5x, C2, lambda, x0, gamma);
    double dP5y_exp = dP5_exp(P5y, C2, lambda, y0, gamma);
    double dP5z_exp = dP5_exp(P5z, C2, lambda, z0, gamma);

    double dP6x_exp = dP6_exp(P6x, C2, lambda, x0, gamma);
    double dP6y_exp = dP6_exp(P6y, C2, lambda, y0, gamma);
    double dP6z_exp = dP6_exp(P6z, C2, lambda, z0, gamma);

    double dP7x_exp = dP7_exp(P7x, C2, lambda, x0, gamma);
    double dP7y_exp = dP7_exp(P7y, C2, lambda, y0, gamma);
    double dP7z_exp = dP7_exp(P7z, C2, lambda, z0, gamma);

    double dP8x_exp = dP8_exp(P8x, C2, lambda, x0, gamma);
    double dP8y_exp = dP8_exp(P8y, C2, lambda, y0, gamma);
    double dP8z_exp = dP8_exp(P8z, C2, lambda, z0, gamma);

    double dP9x_exp = dP9_exp(P9x, C2, lambda, x0, gamma);
    double dP9y_exp = dP9_exp(P9y, C2, lambda, y0, gamma);
    double dP9z_exp = dP9_exp(P9z, C2, lambda, z0, gamma);

    double dterm_1_dx = (16448670.0 * dP5x_exp * P2y * P2z) - (816480.0 * dP7x_exp * P2y) + (116235.0 * dP5x_exp * P4y) - (816480.0 * dP7x_exp * P2z) + (116235.0 * dP5x_exp * P4z) + (45360.0 * dP9x_exp) - (13707225.0 * dP3x_exp * P4y * P2z) - (13707225.0 * dP3x_exp * P2y * P4z) + (836325.0 * dP3x_exp * P6y) + (836325.0 * dP3x_exp * P6z) + (1460025.0 * dP1x_exp * P6y * P2z) + (3203550.0 * dP1x_exp * P4y * P4z) + (1460025.0 * dP1x_exp * P2y * P6z) - (141750.0 * dP1x_exp * P8y) - (141750.0 * dP1x_exp * P8z);
    double dterm_1_dy = (16448670.0 * P5x * dP2y_exp * P2z) - (816480.0 * P7x * dP2y_exp) + (116235.0 * P5x * dP4y_exp) - (816480.0 * P7x * dP0y_exp * P2z) + (116235.0 * P5x * dP0y_exp * P4z) + (45360.0 * P9x * dP0y_exp) - (13707225.0 * P3x * dP4y_exp * P2z) - (13707225.0 * P3x * dP2y_exp * P4z) + (836325.0 * P3x * dP6y_exp) + (836325.0 * P3x * dP0y_exp * P6z) + (1460025.0 * P1x * dP6y_exp * P2z) + (3203550.0 * P1x * dP4y_exp * P4z) + (1460025.0 * P1x * dP2y_exp * P6z) - (141750.0 * P1x * dP8y_exp) - (141750.0 * P1x * dP0y_exp * P8z);
    double dterm_1_dz = (16448670.0 * P5x * P2y * dP2z_exp) - (816480.0 * P7x * P2y * dP0z_exp) + (116235.0 * P5x * P4y * dP0z_exp) - (816480.0 * P7x * dP2z_exp) + (116235.0 * P5x * dP4z_exp) + (45360.0 * P9x * dP0z_exp) - (13707225.0 * P3x * P4y * dP2z_exp) - (13707225.0 * P3x * P2y * dP4z_exp) + (836325.0 * P3x * P6y * dP0z_exp) + (836325.0 * P3x * dP6z_exp) + (1460025.0 * P1x * P6y * dP2z_exp) + (3203550.0 * P1x * P4y * dP4z_exp) + (1460025.0 * P1x * P2y * dP6z_exp) - (141750.0 * P1x * P8y * dP0z_exp) - (141750.0 * P1x * dP8z_exp);

    double dterm_2_dx = (16448670.0 * dP2x_exp * P5y * P2z) - (816480.0 * dP2x_exp * P7y) + (116235.0 * dP4x_exp * P5y) - (816480.0 * dP0x_exp * P7y * P2z) + (116235.0 * dP0x_exp * P5y * P4z) + (45360.0 * dP0x_exp * P9y) - (13707225.0 * dP4x_exp * P3y * P2z) - (13707225.0 * dP2x_exp * P3y * P4z) + (836325.0 * dP6x_exp * P3y) + (836325.0 * dP0x_exp * P3y * P6z) + (1460025.0 * dP6x_exp * P1y * P2z) + (3203550.0 * dP4x_exp * P1y * P4z) + (1460025.0 * dP2x_exp * P1y * P6z) - (141750.0 * dP8x_exp * P1y) - (141750.0 * dP0x_exp * P1y * P8z);
    double dterm_2_dy = (16448670.0 * P2x * dP5y_exp * P2z) - (816480.0 * P2x * dP7y_exp) + (116235.0 * P4x * dP5y_exp) - (816480.0 * dP7y_exp * P2z) + (116235.0 * dP5y_exp * P4z) + (45360.0 * dP9y_exp) - (13707225.0 * P4x * dP3y_exp * P2z) - (13707225.0 * P2x * dP3y_exp * P4z) + (836325.0 * P6x * dP3y_exp) + (836325.0 * dP3y_exp * P6z) + (1460025.0 * P6x * dP1y_exp * P2z) + (3203550.0 * P4x * dP1y_exp * P4z) + (1460025.0 * P2x * dP1y_exp * P6z) - (141750.0 * P8x * dP1y_exp) - (141750.0 * dP1y_exp * P8z);
    double dterm_2_dz = (16448670.0 * P2x * P5y * dP2z_exp) - (816480.0 * P2x * P7y * dP0z_exp) + (116235.0 * P4x * P5y * dP0z_exp) - (816480.0 * P7y * dP2z_exp) + (116235.0 * P5y * dP4z_exp) + (45360.0 * P9y * dP0z_exp) - (13707225.0 * P4x * P3y * dP2z_exp) - (13707225.0 * P2x * P3y * dP4z_exp) + (836325.0 * P6x * P3y * dP0z_exp) + (836325.0 * P3y * dP6z_exp) + (1460025.0 * P6x * P1y * dP2z_exp) + (3203550.0 * P4x * P1y * dP4z_exp) + (1460025.0 * P2x * P1y * dP6z_exp) - (141750.0 * P8x * P1y * dP0z_exp) - (141750.0 * P1y * dP8z_exp);

    double dterm_3_dx = (16448670.0 * dP2x_exp * P2y * P5z) - (816480.0 * dP2x_exp * P7z) + (116235.0 * dP4x_exp * P5z) - (816480.0 * dP0x_exp * P2y * P7z) + (116235.0 * dP0x_exp * P4y * P5z) + (45360.0 * dP0x_exp * P9z) - (13707225.0 * dP4x_exp * P2y * P3z) - (13707225.0 * dP2x_exp * P4y * P3z) + (836325.0 * dP6x_exp * P3z) + (836325.0 * dP0x_exp * P6y * P3z) + (1460025.0 * dP6x_exp * P2y * P1z) + (3203550.0 * dP4x_exp * P4y * P1z) + (1460025.0 * dP2x_exp * P6y * P1z) - (141750.0 * dP8x_exp * P1z) - (141750.0 * dP0x_exp * P8y * P1z);
    double dterm_3_dy = (16448670.0 * P2x * dP2y_exp * P5z) - (816480.0 * P2x * dP0y_exp * P7z) + (116235.0 * P4x * dP0y_exp * P5z) - (816480.0 * dP2y_exp * P7z) + (116235.0 * dP4y_exp * P5z) + (45360.0 * dP0y_exp * P9z) - (13707225.0 * P4x * dP2y_exp * P3z) - (13707225.0 * P2x * dP4y_exp * P3z) + (836325.0 * P6x * dP0y_exp * P3z) + (836325.0 * dP6y_exp * P3z) + (1460025.0 * P6x * dP2y_exp * P1z) + (3203550.0 * P4x * dP4y_exp * P1z) + (1460025.0 * P2x * dP6y_exp * P1z) - (141750.0 * P8x * dP0y_exp * P1z) - (141750.0 * dP8y_exp * P1z);
    double dterm_3_dz = (16448670.0 * P2x * P2y * dP5z_exp) - (816480.0 * P2x * dP7z_exp) + (116235.0 * P4x * dP5z_exp) - (816480.0 * P2y * dP7z_exp) + (116235.0 * P4y * dP5z_exp) + (45360.0 * dP9z_exp) - (13707225.0 * P4x * P2y * dP3z_exp) - (13707225.0 * P2x * P4y * dP3z_exp) + (836325.0 * P6x * dP3z_exp) + (836325.0 * P6y * dP3z_exp) + (1460025.0 * P6x * P2y * dP1z_exp) + (3203550.0 * P4x * P4y * dP1z_exp) + (1460025.0 * P2x * P6y * dP1z_exp) - (141750.0 * P8x * dP1z_exp) - (141750.0 * P8y * dP1z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

}

void calc_solid_MCSH_9_10(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double P8x = P8(lambda, x0, gamma);
    double P8y = P8(lambda, y0, gamma);
    double P8z = P8(lambda, z0, gamma);

    double P9x = P9(lambda, x0, gamma);
    double P9y = P9(lambda, y0, gamma);
    double P9z = P9(lambda, z0, gamma);

    double term_1 = (19604025.0 * P4x * P4y * P1z) - (7200900.0 * P6x * P2y * P1z) - (3203550.0 * P4x * P2y * P3z) + (226800.0 * P8x * P1z) + (283500.0 * P6x * P3z) - (104895.0 * P4x * P5z) - (7200900.0 * P2x * P6y * P1z) - (3203550.0 * P2x * P4y * P3z) + (3844260.0 * P2x * P2y * P5z) - (153090.0 * P2x * P7z) + (226800.0 * P8y * P1z) + (283500.0 * P6y * P3z) - (104895.0 * P4y * P5z) - (153090.0 * P2y * P7z) + (8505.0 * P9z);
    double term_2 = (19604025.0 * P4x * P1y * P4z) - (7200900.0 * P6x * P1y * P2z) - (3203550.0 * P4x * P3y * P2z) + (226800.0 * P8x * P1y) + (283500.0 * P6x * P3y) - (104895.0 * P4x * P5y) - (7200900.0 * P2x * P1y * P6z) - (3203550.0 * P2x * P3y * P4z) + (3844260.0 * P2x * P5y * P2z) - (153090.0 * P2x * P7y) + (226800.0 * P1y * P8z) + (283500.0 * P3y * P6z) - (104895.0 * P5y * P4z) - (153090.0 * P7y * P2z) + (8505.0 * P9y);
    double term_3 = (19604025.0 * P1x * P4y * P4z) - (7200900.0 * P1x * P6y * P2z) - (3203550.0 * P3x * P4y * P2z) + (226800.0 * P1x * P8y) + (283500.0 * P3x * P6y) - (104895.0 * P5x * P4y) - (7200900.0 * P1x * P2y * P6z) - (3203550.0 * P3x * P2y * P4z) + (3844260.0 * P5x * P2y * P2z) - (153090.0 * P7x * P2y) + (226800.0 * P1x * P8z) + (283500.0 * P3x * P6z) - (104895.0 * P5x * P4z) - (153090.0 * P7x * P2z) + (8505.0 * P9x);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP3x_exp = dP3_exp(P3x, C2, lambda, x0, gamma);
    double dP3y_exp = dP3_exp(P3y, C2, lambda, y0, gamma);
    double dP3z_exp = dP3_exp(P3z, C2, lambda, z0, gamma);

    double dP4x_exp = dP4_exp(P4x, C2, lambda, x0, gamma);
    double dP4y_exp = dP4_exp(P4y, C2, lambda, y0, gamma);
    double dP4z_exp = dP4_exp(P4z, C2, lambda, z0, gamma);

    double dP5x_exp = dP5_exp(P5x, C2, lambda, x0, gamma);
    double dP5y_exp = dP5_exp(P5y, C2, lambda, y0, gamma);
    double dP5z_exp = dP5_exp(P5z, C2, lambda, z0, gamma);

    double dP6x_exp = dP6_exp(P6x, C2, lambda, x0, gamma);
    double dP6y_exp = dP6_exp(P6y, C2, lambda, y0, gamma);
    double dP6z_exp = dP6_exp(P6z, C2, lambda, z0, gamma);

    double dP7x_exp = dP7_exp(P7x, C2, lambda, x0, gamma);
    double dP7y_exp = dP7_exp(P7y, C2, lambda, y0, gamma);
    double dP7z_exp = dP7_exp(P7z, C2, lambda, z0, gamma);

    double dP8x_exp = dP8_exp(P8x, C2, lambda, x0, gamma);
    double dP8y_exp = dP8_exp(P8y, C2, lambda, y0, gamma);
    double dP8z_exp = dP8_exp(P8z, C2, lambda, z0, gamma);

    double dP9x_exp = dP9_exp(P9x, C2, lambda, x0, gamma);
    double dP9y_exp = dP9_exp(P9y, C2, lambda, y0, gamma);
    double dP9z_exp = dP9_exp(P9z, C2, lambda, z0, gamma);

    double dterm_1_dx = (19604025.0 * dP4x_exp * P4y * P1z) - (7200900.0 * dP6x_exp * P2y * P1z) - (3203550.0 * dP4x_exp * P2y * P3z) + (226800.0 * dP8x_exp * P1z) + (283500.0 * dP6x_exp * P3z) - (104895.0 * dP4x_exp * P5z) - (7200900.0 * dP2x_exp * P6y * P1z) - (3203550.0 * dP2x_exp * P4y * P3z) + (3844260.0 * dP2x_exp * P2y * P5z) - (153090.0 * dP2x_exp * P7z) + (226800.0 * dP0x_exp * P8y * P1z) + (283500.0 * dP0x_exp * P6y * P3z) - (104895.0 * dP0x_exp * P4y * P5z) - (153090.0 * dP0x_exp * P2y * P7z) + (8505.0 * dP0x_exp * P9z);
    double dterm_1_dy = (19604025.0 * P4x * dP4y_exp * P1z) - (7200900.0 * P6x * dP2y_exp * P1z) - (3203550.0 * P4x * dP2y_exp * P3z) + (226800.0 * P8x * dP0y_exp * P1z) + (283500.0 * P6x * dP0y_exp * P3z) - (104895.0 * P4x * dP0y_exp * P5z) - (7200900.0 * P2x * dP6y_exp * P1z) - (3203550.0 * P2x * dP4y_exp * P3z) + (3844260.0 * P2x * dP2y_exp * P5z) - (153090.0 * P2x * dP0y_exp * P7z) + (226800.0 * dP8y_exp * P1z) + (283500.0 * dP6y_exp * P3z) - (104895.0 * dP4y_exp * P5z) - (153090.0 * dP2y_exp * P7z) + (8505.0 * dP0y_exp * P9z);
    double dterm_1_dz = (19604025.0 * P4x * P4y * dP1z_exp) - (7200900.0 * P6x * P2y * dP1z_exp) - (3203550.0 * P4x * P2y * dP3z_exp) + (226800.0 * P8x * dP1z_exp) + (283500.0 * P6x * dP3z_exp) - (104895.0 * P4x * dP5z_exp) - (7200900.0 * P2x * P6y * dP1z_exp) - (3203550.0 * P2x * P4y * dP3z_exp) + (3844260.0 * P2x * P2y * dP5z_exp) - (153090.0 * P2x * dP7z_exp) + (226800.0 * P8y * dP1z_exp) + (283500.0 * P6y * dP3z_exp) - (104895.0 * P4y * dP5z_exp) - (153090.0 * P2y * dP7z_exp) + (8505.0 * dP9z_exp);

    double dterm_2_dx = (19604025.0 * dP4x_exp * P1y * P4z) - (7200900.0 * dP6x_exp * P1y * P2z) - (3203550.0 * dP4x_exp * P3y * P2z) + (226800.0 * dP8x_exp * P1y) + (283500.0 * dP6x_exp * P3y) - (104895.0 * dP4x_exp * P5y) - (7200900.0 * dP2x_exp * P1y * P6z) - (3203550.0 * dP2x_exp * P3y * P4z) + (3844260.0 * dP2x_exp * P5y * P2z) - (153090.0 * dP2x_exp * P7y) + (226800.0 * dP0x_exp * P1y * P8z) + (283500.0 * dP0x_exp * P3y * P6z) - (104895.0 * dP0x_exp * P5y * P4z) - (153090.0 * dP0x_exp * P7y * P2z) + (8505.0 * dP0x_exp * P9y);
    double dterm_2_dy = (19604025.0 * P4x * dP1y_exp * P4z) - (7200900.0 * P6x * dP1y_exp * P2z) - (3203550.0 * P4x * dP3y_exp * P2z) + (226800.0 * P8x * dP1y_exp) + (283500.0 * P6x * dP3y_exp) - (104895.0 * P4x * dP5y_exp) - (7200900.0 * P2x * dP1y_exp * P6z) - (3203550.0 * P2x * dP3y_exp * P4z) + (3844260.0 * P2x * dP5y_exp * P2z) - (153090.0 * P2x * dP7y_exp) + (226800.0 * dP1y_exp * P8z) + (283500.0 * dP3y_exp * P6z) - (104895.0 * dP5y_exp * P4z) - (153090.0 * dP7y_exp * P2z) + (8505.0 * dP9y_exp);
    double dterm_2_dz = (19604025.0 * P4x * P1y * dP4z_exp) - (7200900.0 * P6x * P1y * dP2z_exp) - (3203550.0 * P4x * P3y * dP2z_exp) + (226800.0 * P8x * P1y * dP0z_exp) + (283500.0 * P6x * P3y * dP0z_exp) - (104895.0 * P4x * P5y * dP0z_exp) - (7200900.0 * P2x * P1y * dP6z_exp) - (3203550.0 * P2x * P3y * dP4z_exp) + (3844260.0 * P2x * P5y * dP2z_exp) - (153090.0 * P2x * P7y * dP0z_exp) + (226800.0 * P1y * dP8z_exp) + (283500.0 * P3y * dP6z_exp) - (104895.0 * P5y * dP4z_exp) - (153090.0 * P7y * dP2z_exp) + (8505.0 * P9y * dP0z_exp);

    double dterm_3_dx = (19604025.0 * dP1x_exp * P4y * P4z) - (7200900.0 * dP1x_exp * P6y * P2z) - (3203550.0 * dP3x_exp * P4y * P2z) + (226800.0 * dP1x_exp * P8y) + (283500.0 * dP3x_exp * P6y) - (104895.0 * dP5x_exp * P4y) - (7200900.0 * dP1x_exp * P2y * P6z) - (3203550.0 * dP3x_exp * P2y * P4z) + (3844260.0 * dP5x_exp * P2y * P2z) - (153090.0 * dP7x_exp * P2y) + (226800.0 * dP1x_exp * P8z) + (283500.0 * dP3x_exp * P6z) - (104895.0 * dP5x_exp * P4z) - (153090.0 * dP7x_exp * P2z) + (8505.0 * dP9x_exp);
    double dterm_3_dy = (19604025.0 * P1x * dP4y_exp * P4z) - (7200900.0 * P1x * dP6y_exp * P2z) - (3203550.0 * P3x * dP4y_exp * P2z) + (226800.0 * P1x * dP8y_exp) + (283500.0 * P3x * dP6y_exp) - (104895.0 * P5x * dP4y_exp) - (7200900.0 * P1x * dP2y_exp * P6z) - (3203550.0 * P3x * dP2y_exp * P4z) + (3844260.0 * P5x * dP2y_exp * P2z) - (153090.0 * P7x * dP2y_exp) + (226800.0 * P1x * dP0y_exp * P8z) + (283500.0 * P3x * dP0y_exp * P6z) - (104895.0 * P5x * dP0y_exp * P4z) - (153090.0 * P7x * dP0y_exp * P2z) + (8505.0 * P9x * dP0y_exp);
    double dterm_3_dz = (19604025.0 * P1x * P4y * dP4z_exp) - (7200900.0 * P1x * P6y * dP2z_exp) - (3203550.0 * P3x * P4y * dP2z_exp) + (226800.0 * P1x * P8y * dP0z_exp) + (283500.0 * P3x * P6y * dP0z_exp) - (104895.0 * P5x * P4y * dP0z_exp) - (7200900.0 * P1x * P2y * dP6z_exp) - (3203550.0 * P3x * P2y * dP4z_exp) + (3844260.0 * P5x * P2y * dP2z_exp) - (153090.0 * P7x * P2y * dP0z_exp) + (226800.0 * P1x * dP8z_exp) + (283500.0 * P3x * dP6z_exp) - (104895.0 * P5x * dP4z_exp) - (153090.0 * P7x * dP2z_exp) + (8505.0 * P9x * dP0z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

}

void calc_solid_MCSH_9_11(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double P5x = P5(lambda, x0, gamma);
    double P5y = P5(lambda, y0, gamma);
    double P5z = P5(lambda, z0, gamma);

    double P6x = P6(lambda, x0, gamma);
    double P6y = P6(lambda, y0, gamma);
    double P6z = P6(lambda, z0, gamma);

    double P7x = P7(lambda, x0, gamma);
    double P7y = P7(lambda, y0, gamma);
    double P7z = P7(lambda, z0, gamma);

    double P8x = P8(lambda, x0, gamma);
    double P8y = P8(lambda, y0, gamma);
    double P8z = P8(lambda, z0, gamma);

    double P9x = P9(lambda, x0, gamma);
    double P9y = P9(lambda, y0, gamma);
    double P9z = P9(lambda, z0, gamma);

    double term_1 = (20497050.0 * P4x * P3y * P2z) - (963900.0 * P6x * P3y) - (603855.0 * P4x * P5y) - (3458700.0 * P6x * P1y * P2z) - (1601775.0 * P4x * P1y * P4z) + (226800.0 * P8x * P1y) - (8224335.0 * P2x * P5y * P2z) - (6789825.0 * P2x * P3y * P4z) + (564165.0 * P2x * P7y) + (1998675.0 * P2x * P1y * P6z) + (252315.0 * P7y * P2z) + (487620.0 * P5y * P4z) + (127575.0 * P3y * P6z) - (22680.0 * P9y) - (85050.0 * P1y * P8z);
    double term_2 = (20497050.0 * P3x * P4y * P2z) - (963900.0 * P3x * P6y) - (603855.0 * P5x * P4y) - (3458700.0 * P1x * P6y * P2z) - (1601775.0 * P1x * P4y * P4z) + (226800.0 * P1x * P8y) - (8224335.0 * P5x * P2y * P2z) - (6789825.0 * P3x * P2y * P4z) + (564165.0 * P7x * P2y) + (1998675.0 * P1x * P2y * P6z) + (252315.0 * P7x * P2z) + (487620.0 * P5x * P4z) + (127575.0 * P3x * P6z) - (22680.0 * P9x) - (85050.0 * P1x * P8z);
    double term_3 = (20497050.0 * P4x * P2y * P3z) - (963900.0 * P6x * P3z) - (603855.0 * P4x * P5z) - (3458700.0 * P6x * P2y * P1z) - (1601775.0 * P4x * P4y * P1z) + (226800.0 * P8x * P1z) - (8224335.0 * P2x * P2y * P5z) - (6789825.0 * P2x * P4y * P3z) + (564165.0 * P2x * P7z) + (1998675.0 * P2x * P6y * P1z) + (252315.0 * P2y * P7z) + (487620.0 * P4y * P5z) + (127575.0 * P6y * P3z) - (22680.0 * P9z) - (85050.0 * P8y * P1z);
    double term_4 = (20497050.0 * P3x * P2y * P4z) - (963900.0 * P3x * P6z) - (603855.0 * P5x * P4z) - (3458700.0 * P1x * P2y * P6z) - (1601775.0 * P1x * P4y * P4z) + (226800.0 * P1x * P8z) - (8224335.0 * P5x * P2y * P2z) - (6789825.0 * P3x * P4y * P2z) + (564165.0 * P7x * P2z) + (1998675.0 * P1x * P6y * P2z) + (252315.0 * P7x * P2y) + (487620.0 * P5x * P4y) + (127575.0 * P3x * P6y) - (22680.0 * P9x) - (85050.0 * P1x * P8y);
    double term_5 = (20497050.0 * P2x * P4y * P3z) - (963900.0 * P6y * P3z) - (603855.0 * P4y * P5z) - (3458700.0 * P2x * P6y * P1z) - (1601775.0 * P4x * P4y * P1z) + (226800.0 * P8y * P1z) - (8224335.0 * P2x * P2y * P5z) - (6789825.0 * P4x * P2y * P3z) + (564165.0 * P2y * P7z) + (1998675.0 * P6x * P2y * P1z) + (252315.0 * P2x * P7z) + (487620.0 * P4x * P5z) + (127575.0 * P6x * P3z) - (22680.0 * P9z) - (85050.0 * P8x * P1z);
    double term_6 = (20497050.0 * P2x * P3y * P4z) - (963900.0 * P3y * P6z) - (603855.0 * P5y * P4z) - (3458700.0 * P2x * P1y * P6z) - (1601775.0 * P4x * P1y * P4z) + (226800.0 * P1y * P8z) - (8224335.0 * P2x * P5y * P2z) - (6789825.0 * P4x * P3y * P2z) + (564165.0 * P7y * P2z) + (1998675.0 * P6x * P1y * P2z) + (252315.0 * P2x * P7y) + (487620.0 * P4x * P5y) + (127575.0 * P6x * P3y) - (22680.0 * P9y) - (85050.0 * P8x * P1y);

    double miu_1 = temp * term_1;
    double miu_2 = temp * term_2;
    double miu_3 = temp * term_3;
    double miu_4 = temp * term_4;
    double miu_5 = temp * term_5;
    double miu_6 = temp * term_6;

    value[0] = miu_1;
    value[1] = miu_2;
    value[2] = miu_3;
    value[3] = miu_4;
    value[4] = miu_5;
    value[5] = miu_6;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dP2x_exp = dP2_exp(P2x, C2, lambda, x0, gamma);
    double dP2y_exp = dP2_exp(P2y, C2, lambda, y0, gamma);
    double dP2z_exp = dP2_exp(P2z, C2, lambda, z0, gamma);

    double dP3x_exp = dP3_exp(P3x, C2, lambda, x0, gamma);
    double dP3y_exp = dP3_exp(P3y, C2, lambda, y0, gamma);
    double dP3z_exp = dP3_exp(P3z, C2, lambda, z0, gamma);

    double dP4x_exp = dP4_exp(P4x, C2, lambda, x0, gamma);
    double dP4y_exp = dP4_exp(P4y, C2, lambda, y0, gamma);
    double dP4z_exp = dP4_exp(P4z, C2, lambda, z0, gamma);

    double dP5x_exp = dP5_exp(P5x, C2, lambda, x0, gamma);
    double dP5y_exp = dP5_exp(P5y, C2, lambda, y0, gamma);
    double dP5z_exp = dP5_exp(P5z, C2, lambda, z0, gamma);

    double dP6x_exp = dP6_exp(P6x, C2, lambda, x0, gamma);
    double dP6y_exp = dP6_exp(P6y, C2, lambda, y0, gamma);
    double dP6z_exp = dP6_exp(P6z, C2, lambda, z0, gamma);

    double dP7x_exp = dP7_exp(P7x, C2, lambda, x0, gamma);
    double dP7y_exp = dP7_exp(P7y, C2, lambda, y0, gamma);
    double dP7z_exp = dP7_exp(P7z, C2, lambda, z0, gamma);

    double dP8x_exp = dP8_exp(P8x, C2, lambda, x0, gamma);
    double dP8y_exp = dP8_exp(P8y, C2, lambda, y0, gamma);
    double dP8z_exp = dP8_exp(P8z, C2, lambda, z0, gamma);

    double dP9x_exp = dP9_exp(P9x, C2, lambda, x0, gamma);
    double dP9y_exp = dP9_exp(P9y, C2, lambda, y0, gamma);
    double dP9z_exp = dP9_exp(P9z, C2, lambda, z0, gamma);

    double dterm_1_dx = (20497050.0 * dP4x_exp * P3y * P2z) - (963900.0 * dP6x_exp * P3y) - (603855.0 * dP4x_exp * P5y) - (3458700.0 * dP6x_exp * P1y * P2z) - (1601775.0 * dP4x_exp * P1y * P4z) + (226800.0 * dP8x_exp * P1y) - (8224335.0 * dP2x_exp * P5y * P2z) - (6789825.0 * dP2x_exp * P3y * P4z) + (564165.0 * dP2x_exp * P7y) + (1998675.0 * dP2x_exp * P1y * P6z) + (252315.0 * dP0x_exp * P7y * P2z) + (487620.0 * dP0x_exp * P5y * P4z) + (127575.0 * dP0x_exp * P3y * P6z) - (22680.0 * dP0x_exp * P9y) - (85050.0 * dP0x_exp * P1y * P8z);
    double dterm_1_dy = (20497050.0 * P4x * dP3y_exp * P2z) - (963900.0 * P6x * dP3y_exp) - (603855.0 * P4x * dP5y_exp) - (3458700.0 * P6x * dP1y_exp * P2z) - (1601775.0 * P4x * dP1y_exp * P4z) + (226800.0 * P8x * dP1y_exp) - (8224335.0 * P2x * dP5y_exp * P2z) - (6789825.0 * P2x * dP3y_exp * P4z) + (564165.0 * P2x * dP7y_exp) + (1998675.0 * P2x * dP1y_exp * P6z) + (252315.0 * dP7y_exp * P2z) + (487620.0 * dP5y_exp * P4z) + (127575.0 * dP3y_exp * P6z) - (22680.0 * dP9y_exp) - (85050.0 * dP1y_exp * P8z);
    double dterm_1_dz = (20497050.0 * P4x * P3y * dP2z_exp) - (963900.0 * P6x * P3y * dP0z_exp) - (603855.0 * P4x * P5y * dP0z_exp) - (3458700.0 * P6x * P1y * dP2z_exp) - (1601775.0 * P4x * P1y * dP4z_exp) + (226800.0 * P8x * P1y * dP0z_exp) - (8224335.0 * P2x * P5y * dP2z_exp) - (6789825.0 * P2x * P3y * dP4z_exp) + (564165.0 * P2x * P7y * dP0z_exp) + (1998675.0 * P2x * P1y * dP6z_exp) + (252315.0 * P7y * dP2z_exp) + (487620.0 * P5y * dP4z_exp) + (127575.0 * P3y * dP6z_exp) - (22680.0 * P9y * dP0z_exp) - (85050.0 * P1y * dP8z_exp);

    double dterm_2_dx = (20497050.0 * dP3x_exp * P4y * P2z) - (963900.0 * dP3x_exp * P6y) - (603855.0 * dP5x_exp * P4y) - (3458700.0 * dP1x_exp * P6y * P2z) - (1601775.0 * dP1x_exp * P4y * P4z) + (226800.0 * dP1x_exp * P8y) - (8224335.0 * dP5x_exp * P2y * P2z) - (6789825.0 * dP3x_exp * P2y * P4z) + (564165.0 * dP7x_exp * P2y) + (1998675.0 * dP1x_exp * P2y * P6z) + (252315.0 * dP7x_exp * P2z) + (487620.0 * dP5x_exp * P4z) + (127575.0 * dP3x_exp * P6z) - (22680.0 * dP9x_exp) - (85050.0 * dP1x_exp * P8z);
    double dterm_2_dy = (20497050.0 * P3x * dP4y_exp * P2z) - (963900.0 * P3x * dP6y_exp) - (603855.0 * P5x * dP4y_exp) - (3458700.0 * P1x * dP6y_exp * P2z) - (1601775.0 * P1x * dP4y_exp * P4z) + (226800.0 * P1x * dP8y_exp) - (8224335.0 * P5x * dP2y_exp * P2z) - (6789825.0 * P3x * dP2y_exp * P4z) + (564165.0 * P7x * dP2y_exp) + (1998675.0 * P1x * dP2y_exp * P6z) + (252315.0 * P7x * dP0y_exp * P2z) + (487620.0 * P5x * dP0y_exp * P4z) + (127575.0 * P3x * dP0y_exp * P6z) - (22680.0 * P9x * dP0y_exp) - (85050.0 * P1x * dP0y_exp * P8z);
    double dterm_2_dz = (20497050.0 * P3x * P4y * dP2z_exp) - (963900.0 * P3x * P6y * dP0z_exp) - (603855.0 * P5x * P4y * dP0z_exp) - (3458700.0 * P1x * P6y * dP2z_exp) - (1601775.0 * P1x * P4y * dP4z_exp) + (226800.0 * P1x * P8y * dP0z_exp) - (8224335.0 * P5x * P2y * dP2z_exp) - (6789825.0 * P3x * P2y * dP4z_exp) + (564165.0 * P7x * P2y * dP0z_exp) + (1998675.0 * P1x * P2y * dP6z_exp) + (252315.0 * P7x * dP2z_exp) + (487620.0 * P5x * dP4z_exp) + (127575.0 * P3x * dP6z_exp) - (22680.0 * P9x * dP0z_exp) - (85050.0 * P1x * dP8z_exp);

    double dterm_3_dx = (20497050.0 * dP4x_exp * P2y * P3z) - (963900.0 * dP6x_exp * P3z) - (603855.0 * dP4x_exp * P5z) - (3458700.0 * dP6x_exp * P2y * P1z) - (1601775.0 * dP4x_exp * P4y * P1z) + (226800.0 * dP8x_exp * P1z) - (8224335.0 * dP2x_exp * P2y * P5z) - (6789825.0 * dP2x_exp * P4y * P3z) + (564165.0 * dP2x_exp * P7z) + (1998675.0 * dP2x_exp * P6y * P1z) + (252315.0 * dP0x_exp * P2y * P7z) + (487620.0 * dP0x_exp * P4y * P5z) + (127575.0 * dP0x_exp * P6y * P3z) - (22680.0 * dP0x_exp * P9z) - (85050.0 * dP0x_exp * P8y * P1z);
    double dterm_3_dy = (20497050.0 * P4x * dP2y_exp * P3z) - (963900.0 * P6x * dP0y_exp * P3z) - (603855.0 * P4x * dP0y_exp * P5z) - (3458700.0 * P6x * dP2y_exp * P1z) - (1601775.0 * P4x * dP4y_exp * P1z) + (226800.0 * P8x * dP0y_exp * P1z) - (8224335.0 * P2x * dP2y_exp * P5z) - (6789825.0 * P2x * dP4y_exp * P3z) + (564165.0 * P2x * dP0y_exp * P7z) + (1998675.0 * P2x * dP6y_exp * P1z) + (252315.0 * dP2y_exp * P7z) + (487620.0 * dP4y_exp * P5z) + (127575.0 * dP6y_exp * P3z) - (22680.0 * dP0y_exp * P9z) - (85050.0 * dP8y_exp * P1z);
    double dterm_3_dz = (20497050.0 * P4x * P2y * dP3z_exp) - (963900.0 * P6x * dP3z_exp) - (603855.0 * P4x * dP5z_exp) - (3458700.0 * P6x * P2y * dP1z_exp) - (1601775.0 * P4x * P4y * dP1z_exp) + (226800.0 * P8x * dP1z_exp) - (8224335.0 * P2x * P2y * dP5z_exp) - (6789825.0 * P2x * P4y * dP3z_exp) + (564165.0 * P2x * dP7z_exp) + (1998675.0 * P2x * P6y * dP1z_exp) + (252315.0 * P2y * dP7z_exp) + (487620.0 * P4y * dP5z_exp) + (127575.0 * P6y * dP3z_exp) - (22680.0 * dP9z_exp) - (85050.0 * P8y * dP1z_exp);

    double dterm_4_dx = (20497050.0 * dP3x_exp * P2y * P4z) - (963900.0 * dP3x_exp * P6z) - (603855.0 * dP5x_exp * P4z) - (3458700.0 * dP1x_exp * P2y * P6z) - (1601775.0 * dP1x_exp * P4y * P4z) + (226800.0 * dP1x_exp * P8z) - (8224335.0 * dP5x_exp * P2y * P2z) - (6789825.0 * dP3x_exp * P4y * P2z) + (564165.0 * dP7x_exp * P2z) + (1998675.0 * dP1x_exp * P6y * P2z) + (252315.0 * dP7x_exp * P2y) + (487620.0 * dP5x_exp * P4y) + (127575.0 * dP3x_exp * P6y) - (22680.0 * dP9x_exp) - (85050.0 * dP1x_exp * P8y);
    double dterm_4_dy = (20497050.0 * P3x * dP2y_exp * P4z) - (963900.0 * P3x * dP0y_exp * P6z) - (603855.0 * P5x * dP0y_exp * P4z) - (3458700.0 * P1x * dP2y_exp * P6z) - (1601775.0 * P1x * dP4y_exp * P4z) + (226800.0 * P1x * dP0y_exp * P8z) - (8224335.0 * P5x * dP2y_exp * P2z) - (6789825.0 * P3x * dP4y_exp * P2z) + (564165.0 * P7x * dP0y_exp * P2z) + (1998675.0 * P1x * dP6y_exp * P2z) + (252315.0 * P7x * dP2y_exp) + (487620.0 * P5x * dP4y_exp) + (127575.0 * P3x * dP6y_exp) - (22680.0 * P9x * dP0y_exp) - (85050.0 * P1x * dP8y_exp);
    double dterm_4_dz = (20497050.0 * P3x * P2y * dP4z_exp) - (963900.0 * P3x * dP6z_exp) - (603855.0 * P5x * dP4z_exp) - (3458700.0 * P1x * P2y * dP6z_exp) - (1601775.0 * P1x * P4y * dP4z_exp) + (226800.0 * P1x * dP8z_exp) - (8224335.0 * P5x * P2y * dP2z_exp) - (6789825.0 * P3x * P4y * dP2z_exp) + (564165.0 * P7x * dP2z_exp) + (1998675.0 * P1x * P6y * dP2z_exp) + (252315.0 * P7x * P2y * dP0z_exp) + (487620.0 * P5x * P4y * dP0z_exp) + (127575.0 * P3x * P6y * dP0z_exp) - (22680.0 * P9x * dP0z_exp) - (85050.0 * P1x * P8y * dP0z_exp);

    double dterm_5_dx = (20497050.0 * dP2x_exp * P4y * P3z) - (963900.0 * dP0x_exp * P6y * P3z) - (603855.0 * dP0x_exp * P4y * P5z) - (3458700.0 * dP2x_exp * P6y * P1z) - (1601775.0 * dP4x_exp * P4y * P1z) + (226800.0 * dP0x_exp * P8y * P1z) - (8224335.0 * dP2x_exp * P2y * P5z) - (6789825.0 * dP4x_exp * P2y * P3z) + (564165.0 * dP0x_exp * P2y * P7z) + (1998675.0 * dP6x_exp * P2y * P1z) + (252315.0 * dP2x_exp * P7z) + (487620.0 * dP4x_exp * P5z) + (127575.0 * dP6x_exp * P3z) - (22680.0 * dP0x_exp * P9z) - (85050.0 * dP8x_exp * P1z);
    double dterm_5_dy = (20497050.0 * P2x * dP4y_exp * P3z) - (963900.0 * dP6y_exp * P3z) - (603855.0 * dP4y_exp * P5z) - (3458700.0 * P2x * dP6y_exp * P1z) - (1601775.0 * P4x * dP4y_exp * P1z) + (226800.0 * dP8y_exp * P1z) - (8224335.0 * P2x * dP2y_exp * P5z) - (6789825.0 * P4x * dP2y_exp * P3z) + (564165.0 * dP2y_exp * P7z) + (1998675.0 * P6x * dP2y_exp * P1z) + (252315.0 * P2x * dP0y_exp * P7z) + (487620.0 * P4x * dP0y_exp * P5z) + (127575.0 * P6x * dP0y_exp * P3z) - (22680.0 * dP0y_exp * P9z) - (85050.0 * P8x * dP0y_exp * P1z);
    double dterm_5_dz = (20497050.0 * P2x * P4y * dP3z_exp) - (963900.0 * P6y * dP3z_exp) - (603855.0 * P4y * dP5z_exp) - (3458700.0 * P2x * P6y * dP1z_exp) - (1601775.0 * P4x * P4y * dP1z_exp) + (226800.0 * P8y * dP1z_exp) - (8224335.0 * P2x * P2y * dP5z_exp) - (6789825.0 * P4x * P2y * dP3z_exp) + (564165.0 * P2y * dP7z_exp) + (1998675.0 * P6x * P2y * dP1z_exp) + (252315.0 * P2x * dP7z_exp) + (487620.0 * P4x * dP5z_exp) + (127575.0 * P6x * dP3z_exp) - (22680.0 * dP9z_exp) - (85050.0 * P8x * dP1z_exp);

    double dterm_6_dx = (20497050.0 * dP2x_exp * P3y * P4z) - (963900.0 * dP0x_exp * P3y * P6z) - (603855.0 * dP0x_exp * P5y * P4z) - (3458700.0 * dP2x_exp * P1y * P6z) - (1601775.0 * dP4x_exp * P1y * P4z) + (226800.0 * dP0x_exp * P1y * P8z) - (8224335.0 * dP2x_exp * P5y * P2z) - (6789825.0 * dP4x_exp * P3y * P2z) + (564165.0 * dP0x_exp * P7y * P2z) + (1998675.0 * dP6x_exp * P1y * P2z) + (252315.0 * dP2x_exp * P7y) + (487620.0 * dP4x_exp * P5y) + (127575.0 * dP6x_exp * P3y) - (22680.0 * dP0x_exp * P9y) - (85050.0 * dP8x_exp * P1y);
    double dterm_6_dy = (20497050.0 * P2x * dP3y_exp * P4z) - (963900.0 * dP3y_exp * P6z) - (603855.0 * dP5y_exp * P4z) - (3458700.0 * P2x * dP1y_exp * P6z) - (1601775.0 * P4x * dP1y_exp * P4z) + (226800.0 * dP1y_exp * P8z) - (8224335.0 * P2x * dP5y_exp * P2z) - (6789825.0 * P4x * dP3y_exp * P2z) + (564165.0 * dP7y_exp * P2z) + (1998675.0 * P6x * dP1y_exp * P2z) + (252315.0 * P2x * dP7y_exp) + (487620.0 * P4x * dP5y_exp) + (127575.0 * P6x * dP3y_exp) - (22680.0 * dP9y_exp) - (85050.0 * P8x * dP1y_exp);
    double dterm_6_dz = (20497050.0 * P2x * P3y * dP4z_exp) - (963900.0 * P3y * dP6z_exp) - (603855.0 * P5y * dP4z_exp) - (3458700.0 * P2x * P1y * dP6z_exp) - (1601775.0 * P4x * P1y * dP4z_exp) + (226800.0 * P1y * dP8z_exp) - (8224335.0 * P2x * P5y * dP2z_exp) - (6789825.0 * P4x * P3y * dP2z_exp) + (564165.0 * P7y * dP2z_exp) + (1998675.0 * P6x * P1y * dP2z_exp) + (252315.0 * P2x * P7y * dP0z_exp) + (487620.0 * P4x * P5y * dP0z_exp) + (127575.0 * P6x * P3y * dP0z_exp) - (22680.0 * P9y * dP0z_exp) - (85050.0 * P8x * P1y * dP0z_exp);



    // dmiu1 dx/dy/dz
    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

    // dmiu2 dx/dy/dz
    deriv[3] = temp * dterm_2_dx;
    deriv[4] = temp * dterm_2_dy;
    deriv[5] = temp * dterm_2_dz;

    // dmiu3 dx/dy/dz
    deriv[6] = temp * dterm_3_dx;
    deriv[7] = temp * dterm_3_dy;
    deriv[8] = temp * dterm_3_dz;

    // dmiu4 dx/dy/dz
    deriv[9] = temp * dterm_4_dx;
    deriv[10] = temp * dterm_4_dy;
    deriv[11] = temp * dterm_4_dz;

    // dmiu5 dx/dy/dz
    deriv[12] = temp * dterm_5_dx;
    deriv[13] = temp * dterm_5_dy;
    deriv[14] = temp * dterm_5_dz;

    // dmiu6 dx/dy/dz
    deriv[15] = temp * dterm_6_dx;
    deriv[16] = temp * dterm_6_dy;
    deriv[17] = temp * dterm_6_dz;

}

void calc_solid_MCSH_9_12(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double *value, double *deriv)
{
    double C1 = calc_C1(A,B,alpha,beta);
    double C2 = calc_C2(alpha, beta);
    double temp = C1 * exp( C2 * r0_sqr);

    double lambda = calc_lambda(alpha, beta);
    double gamma = calc_gamma(alpha, beta);

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

    double term_1 = (21829500.0 * P3x * P3y * P3z) - (3274425.0 * P5x * P3y * P1z) - (3274425.0 * P3x * P5y * P1z) - (3274425.0 * P5x * P1y * P3z) - (3274425.0 * P3x * P1y * P5z) + (935550.0 * P7x * P1y * P1z) - (3274425.0 * P1x * P5y * P3z) - (3274425.0 * P1x * P3y * P5z) + (935550.0 * P1x * P7y * P1z) + (935550.0 * P1x * P1y * P7z);

    double m = temp * term_1;

    value[0] = m;


    double dP0x_exp = 2.0 * C2 * x0;
    double dP0y_exp = 2.0 * C2 * y0;
    double dP0z_exp = 2.0 * C2 * z0;
    
    double dP1x_exp = dP1_exp(P1x, C2, lambda, x0, gamma);
    double dP1y_exp = dP1_exp(P1y, C2, lambda, y0, gamma);
    double dP1z_exp = dP1_exp(P1z, C2, lambda, z0, gamma);

    double dP3x_exp = dP3_exp(P3x, C2, lambda, x0, gamma);
    double dP3y_exp = dP3_exp(P3y, C2, lambda, y0, gamma);
    double dP3z_exp = dP3_exp(P3z, C2, lambda, z0, gamma);

    double dP5x_exp = dP5_exp(P5x, C2, lambda, x0, gamma);
    double dP5y_exp = dP5_exp(P5y, C2, lambda, y0, gamma);
    double dP5z_exp = dP5_exp(P5z, C2, lambda, z0, gamma);

    double dP7x_exp = dP7_exp(P7x, C2, lambda, x0, gamma);
    double dP7y_exp = dP7_exp(P7y, C2, lambda, y0, gamma);
    double dP7z_exp = dP7_exp(P7z, C2, lambda, z0, gamma);

    double dterm_1_dx = (21829500.0 * dP3x_exp * P3y * P3z) - (3274425.0 * dP5x_exp * P3y * P1z) - (3274425.0 * dP3x_exp * P5y * P1z) - (3274425.0 * dP5x_exp * P1y * P3z) - (3274425.0 * dP3x_exp * P1y * P5z) + (935550.0 * dP7x_exp * P1y * P1z) - (3274425.0 * dP1x_exp * P5y * P3z) - (3274425.0 * dP1x_exp * P3y * P5z) + (935550.0 * dP1x_exp * P7y * P1z) + (935550.0 * dP1x_exp * P1y * P7z);
    double dterm_1_dy = (21829500.0 * P3x * dP3y_exp * P3z) - (3274425.0 * P5x * dP3y_exp * P1z) - (3274425.0 * P3x * dP5y_exp * P1z) - (3274425.0 * P5x * dP1y_exp * P3z) - (3274425.0 * P3x * dP1y_exp * P5z) + (935550.0 * P7x * dP1y_exp * P1z) - (3274425.0 * P1x * dP5y_exp * P3z) - (3274425.0 * P1x * dP3y_exp * P5z) + (935550.0 * P1x * dP7y_exp * P1z) + (935550.0 * P1x * dP1y_exp * P7z);
    double dterm_1_dz = (21829500.0 * P3x * P3y * dP3z_exp) - (3274425.0 * P5x * P3y * dP1z_exp) - (3274425.0 * P3x * P5y * dP1z_exp) - (3274425.0 * P5x * P1y * dP3z_exp) - (3274425.0 * P3x * P1y * dP5z_exp) + (935550.0 * P7x * P1y * dP1z_exp) - (3274425.0 * P1x * P5y * dP3z_exp) - (3274425.0 * P1x * P3y * dP5z_exp) + (935550.0 * P1x * P7y * dP1z_exp) + (935550.0 * P1x * P1y * dP7z_exp);



    deriv[0] = temp * dterm_1_dx;
    deriv[1] = temp * dterm_1_dy;
    deriv[2] = temp * dterm_1_dz;

}