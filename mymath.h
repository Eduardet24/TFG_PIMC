/*mymath.h*/

#ifndef __MYMATH_H_
#define __MYMATH_H_

#include "main.h"

DOUBLE K0opt(DOUBLE r);
DOUBLE K1opt(DOUBLE r);
DOUBLE K2opt(DOUBLE r);
#define K0(r) K0opt(4./((r)*(r)))
#define K1(r) K1opt(4./((r)*(r)))
#define K2(r) K2opt(4./((r)*(r)))

DOUBLE BesselI0(DOUBLE x);
DOUBLE BesselI1(DOUBLE x);
DOUBLE BesselIN(int N, DOUBLE x);
DOUBLE BesselJ0(DOUBLE x);
DOUBLE BesselJ1(DOUBLE x);
DOUBLE BesselY0(DOUBLE x);
DOUBLE BesselY1(DOUBLE x);
DOUBLE BesselK0(DOUBLE x);
DOUBLE BesselK1(DOUBLE x);
#define BesselK2(r) K2opt(4./((r)*(r)))
double digamma(double x);
double erf_my(double x);
double erfc_my(double x);

double E1(double x);
double      Exponential_Integral_Ei(double x);
long double xExponential_Integral_Ei(long double x);
static long double Continued_Fraction_Ei(long double x);
static long double Power_Series_Ei(long double x);
static long double Argument_Addition_Series_Ei(long double x);

DOUBLE Erfcx(DOUBLE x);

int compare_double(const void *a, const void *b);

#endif

