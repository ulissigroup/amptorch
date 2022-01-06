#include <math.h>





double calc_C1(double A, double B, double alpha, double beta);

double calc_C2(double alpha, double beta);

double calc_lambda(double alpha, double beta);

double calc_gamma(double alpha, double beta);

double P1(double lambda, double x0, double gamma);
double P2(double lambda, double x0, double gamma);
double P3(double lambda, double x0, double gamma);
double P4(double lambda, double x0, double gamma);
double P5(double lambda, double x0, double gamma);
double P6(double lambda, double x0, double gamma);
double P7(double lambda, double x0, double gamma);
double P8(double lambda, double x0, double gamma);
double P9(double lambda, double x0, double gamma);

double dP1(double lambda, double x0, double gamma);
double dP2(double lambda, double x0, double gamma);
double dP3(double lambda, double x0, double gamma);
double dP4(double lambda, double x0, double gamma);
double dP5(double lambda, double x0, double gamma);
double dP6(double lambda, double x0, double gamma);
double dP7(double lambda, double x0, double gamma);
double dP8(double lambda, double x0, double gamma);
double dP9(double lambda, double x0, double gamma);

double dP1_exp(double P, double C2, double lambda, double x0, double gamma);
double dP2_exp(double P, double C2, double lambda, double x0, double gamma);
double dP3_exp(double P, double C2, double lambda, double x0, double gamma);
double dP4_exp(double P, double C2, double lambda, double x0, double gamma);
double dP5_exp(double P, double C2, double lambda, double x0, double gamma);
double dP6_exp(double P, double C2, double lambda, double x0, double gamma);
double dP7_exp(double P, double C2, double lambda, double x0, double gamma);
double dP8_exp(double P, double C2, double lambda, double x0, double gamma);
double dP9_exp(double P, double C2, double lambda, double x0, double gamma);

double get_group_coefficients(int mcsh_order, int group_num);
int get_num_groups(int mcsh_order);
int get_mcsh_type(int mcsh_order, int group_num);