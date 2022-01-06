#include <math.h>
#include <stdio.h>
#include "helper.h"
#include <iostream>


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

double dP1(double lambda, double x0, double gamma){
    return lambda;
}

double dP2(double lambda, double x0, double gamma){
    return 2.0 * lambda * lambda * x0;
}

double dP3(double lambda, double x0, double gamma){
    double lambda_x0 = lambda * x0;
    double lambda_x0_2 = lambda_x0 * lambda_x0;
    return lambda * (3.0 * lambda_x0_2 + (3.0 / (2.0 * gamma)));
}

double dP4(double lambda, double x0, double gamma){
    double lambda_x0 = lambda * x0;
    double lambda_x0_3 = lambda_x0 * lambda_x0 * lambda_x0;
    return lambda * (4.0 * lambda_x0_3 + ((6.0 * lambda_x0) / gamma));
}

double dP5(double lambda, double x0, double gamma){
    double lambda_x0 = lambda * x0;
    double lambda_x0_2 = lambda_x0 * lambda_x0;
    double lambda_x0_4 = lambda_x0_2 * lambda_x0_2;
    return lambda * (5.0 * lambda_x0_4 + ((15.0 * lambda_x0_2) / gamma) + (15.0 / (4.0 * gamma * gamma)));
}

double dP6(double lambda, double x0, double gamma){
    double lambda_x0 = lambda * x0;
    double lambda_x0_2 = lambda_x0 * lambda_x0;
    double lambda_x0_3 = lambda_x0 * lambda_x0_2;
    double lambda_x0_5 = lambda_x0_3 * lambda_x0_2;
    return lambda * (6.0 * lambda_x0_5 + ((30.0 * lambda_x0_3) / gamma) + (45.0 * lambda_x0 / (2.0 * gamma * gamma)));
}

double dP7(double lambda, double x0, double gamma){
    double lambda_x0 = lambda * x0;
    double lambda_x0_2 = lambda_x0 * lambda_x0;
    double lambda_x0_4 = lambda_x0_2 * lambda_x0_2;
    double lambda_x0_6 = lambda_x0_4 * lambda_x0_2;
    return lambda * (7.0 * lambda_x0_6 + ((105.0 * lambda_x0_4) / (2.0 * gamma)) + ((315.0 * lambda_x0_2) / (4.0 * gamma * gamma)) + (105.0 / (8.0 * gamma * gamma * gamma)));
}

double dP8(double lambda, double x0, double gamma){
    double lambda_x0 = lambda * x0;
    double lambda_x0_2 = lambda_x0 * lambda_x0;
    double lambda_x0_3 = lambda_x0 * lambda_x0_2;
    double lambda_x0_5 = lambda_x0_3 * lambda_x0_2;
    double lambda_x0_7 = lambda_x0_5 * lambda_x0_2;
    return lambda * (8.0 * lambda_x0_7 + ((84.0 * lambda_x0_5) / gamma) + (210.0 * lambda_x0_3 / (gamma * gamma)) + ((105.0 * lambda_x0) / (gamma * gamma * gamma)));
}

double dP9(double lambda, double x0, double gamma){
    double lambda_x0 = lambda * x0;
    double lambda_x0_2 = lambda_x0 * lambda_x0;
    double lambda_x0_4 = lambda_x0_2 * lambda_x0_2;
    double lambda_x0_6 = lambda_x0_4 * lambda_x0_2;
    double lambda_x0_8 = lambda_x0_6 * lambda_x0_2;
    return lambda * (9.0 * lambda_x0_8 + ((126.0 * lambda_x0_6) / gamma) + ((945.0 * lambda_x0_4) / (2.0 * gamma * gamma)) + ((945.0 * lambda_x0_2) / (2.0 * gamma * gamma * gamma)) + (945.0 / (16.0 * gamma * gamma * gamma * gamma)));
}
 

double dP1_exp(double P, double C2, double lambda, double x0, double gamma){
    double temp = 2.0 * C2 * x0;
    return P * temp + dP1(lambda, x0, gamma);
}

double dP2_exp(double P, double C2, double lambda, double x0, double gamma){
    double temp = 2.0 * C2 * x0;
    return P * temp + dP2(lambda, x0, gamma);
}

double dP3_exp(double P, double C2, double lambda, double x0, double gamma){
    double temp = 2.0 * C2 * x0;
    return P * temp + dP3(lambda, x0, gamma);
}

double dP4_exp(double P, double C2, double lambda, double x0, double gamma){
    double temp = 2.0 * C2 * x0;
    return P * temp + dP4(lambda, x0, gamma);
}

double dP5_exp(double P, double C2, double lambda, double x0, double gamma){
    double temp = 2.0 * C2 * x0;
    return P * temp + dP5(lambda, x0, gamma);
}

double dP6_exp(double P, double C2, double lambda, double x0, double gamma){
    double temp = 2.0 * C2 * x0;
    return P * temp + dP6(lambda, x0, gamma);
}

double dP7_exp(double P, double C2, double lambda, double x0, double gamma){
    double temp = 2.0 * C2 * x0;
    return P * temp + dP7(lambda, x0, gamma);
}

double dP8_exp(double P, double C2, double lambda, double x0, double gamma){
    double temp = 2.0 * C2 * x0;
    return P * temp + dP8(lambda, x0, gamma);
}

double dP9_exp(double P, double C2, double lambda, double x0, double gamma){
    double temp = 2.0 * C2 * x0;
    return P * temp + dP9(lambda, x0, gamma);
}


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
    } else if (mcsh_order == 9) {
        if (group_num == 1) {
            return 1.0; // 9 0 0 -> 362880.0 / (362880.0 * 1.0 * 1.0)
        } else if (group_num == 2){
            return 9.0; // 8 1 0 -> 362880.0 / (40320.0 * 1.0 * 1.0)
        } else if (group_num == 3){
            return 36.0; // 7 2 0 -> 362880.0 / (5040.0 * 2.0 * 1.0)
        } else if (group_num == 4){
            return 72.0; // 7 1 1 -> 362880.0 / (5040.0 * 1.0 * 1.0)
        } else if (group_num == 5){
            return 84.0; // 6 3 0 -> 362880.0 / (720.0 * 6.0 * 1.0)
        } else if (group_num == 6){
            return 252.0; // 6 2 1 -> 362880.0 / (720.0 * 2.0 * 1.0)
        } else if (group_num == 7){
            return 126.0; // 5 4 0 -> 362880.0 / (120.0 * 24.0 * 1.0)
        } else if (group_num == 8){
            return 504.0; // 5 3 1 -> 362880.0 / (120.0 * 6.0 * 1.0)
        } else if (group_num == 9){
            return 756.0; // 5 2 2 -> 362880.0 / (120.0 * 2.0 * 2.0)
        } else if (group_num == 10){
            return 630.0; // 4 4 1 -> 362880.0 / (24.0 * 24.0 * 1.0)
        } else if (group_num == 11){
            return 1260.0; // 4 3 2 -> 362880.0 / (24.0 * 6.0 * 2.0)
        } else if (group_num == 12){
            return 1680.0; // 3 3 3 -> 362880.0 / (6.0 * 6.0 * 6.0)
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
