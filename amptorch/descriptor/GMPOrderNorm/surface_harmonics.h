#include "helper.h"


typedef void (*GMPFunction) (double, double, double, double, double, double, double, double, double, double *, double *);
GMPFunction get_mcsh_function(int mcsh_order, int group_num);

void calc_MCSH_0_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);

void calc_MCSH_1_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);

void calc_MCSH_2_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_2_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);

void calc_MCSH_3_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_3_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_3_3(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);

void calc_MCSH_4_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_4_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_4_3(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_4_4(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);

void calc_MCSH_5_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_5_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_5_3(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_5_4(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_5_5(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);

void calc_MCSH_6_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_6_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_6_3(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_6_4(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_6_5(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_6_6(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_6_7(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);

void calc_MCSH_7_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_7_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_7_3(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_7_4(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_7_5(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_7_6(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_7_7(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_7_8(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);

void calc_MCSH_8_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_8_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_8_3(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_8_4(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_8_5(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_8_6(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_8_7(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_8_8(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_8_9(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_8_10(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);

void calc_MCSH_9_1(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_9_2(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_9_3(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_9_4(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_9_5(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_9_6(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_9_7(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_9_8(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_9_9(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_9_10(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_9_11(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);
void calc_MCSH_9_12(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value, double *deriv);


typedef void (*GMPFunctionNoderiv) (double, double, double, double, double, double, double, double, double, double *);
GMPFunctionNoderiv get_mcsh_function_noderiv(int mcsh_order, int group_num);

void calc_MCSH_0_1_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);

void calc_MCSH_1_1_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);

void calc_MCSH_2_1_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_2_2_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);

void calc_MCSH_3_1_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_3_2_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_3_3_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);

void calc_MCSH_4_1_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_4_2_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_4_3_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_4_4_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);

void calc_MCSH_5_1_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_5_2_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_5_3_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_5_4_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_5_5_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);

void calc_MCSH_6_1_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_6_2_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_6_3_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_6_4_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_6_5_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_6_6_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_6_7_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);

void calc_MCSH_7_1_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_7_2_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_7_3_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_7_4_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_7_5_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_7_6_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_7_7_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_7_8_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);

void calc_MCSH_8_1_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_8_2_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_8_3_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_8_4_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_8_5_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_8_6_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_8_7_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_8_8_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_8_9_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_8_10_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);

void calc_MCSH_9_1_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_9_2_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_9_3_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_9_4_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_9_5_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_9_6_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_9_7_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_9_8_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_9_9_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_9_10_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_9_11_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);
void calc_MCSH_9_12_noderiv(double x0, double y0, double z0, double r0_sqr, double A, double B, double alpha, double beta, double inv_rs, double *value);

