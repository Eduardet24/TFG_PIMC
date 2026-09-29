/*spline.h*/

#ifndef __SPLINE_H_
#define __SPLINE_H_

#include "main.h"

void SaveWaveFunctionSpline(struct Grid *G, char *file_out);
void SplineConstruct(struct Grid *G, DOUBLE yp1, DOUBLE ypn);
void InterpolateSpline(struct Grid *G, DOUBLE r, DOUBLE *y, DOUBLE *yp, DOUBLE *ypp);
DOUBLE InterpolateSplineU(struct Grid *G, DOUBLE r);
DOUBLE InterpolateSplineF(struct Grid *G, DOUBLE r);
DOUBLE InterpolateSplineFp(struct Grid *G, DOUBLE r);
DOUBLE InterpolateSplineE(struct Grid *G, DOUBLE r);
double InterpolateKurbakovU(double r);
double InterpolateKurbakovFp(double r);
double InterpolateKurbakovE(double r);

void LinearInterpolationConstruct(struct Grid *G);

void ConstructGridKurbakovNumeric(struct Grid *G);
void SavingSplineArrays(void);

#endif
