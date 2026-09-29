/*optimiz.h*/

#ifndef _OPTIMIZ_H_
#define _OPTIMIZ_H_

#include "main.h"

void Optimization(void);
DOUBLE OptimizationGradient(double *par);
void OptimizationAnnealing(double *par);
void OptimizationEnergy(double Etarget, double *Emean, double *sigma);
int SaveCoordinatesOptimization(void);
int LoadCoordinatesOptimization(void);

extern DOUBLE opt_Etarget;
extern DOUBLE opt_sigma;

#endif

