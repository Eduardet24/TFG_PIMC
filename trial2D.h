/*trial2D.h*/

#ifndef __TRIAL2D_H_
#define __TRIAL2D_H_

#include "main.h"

struct Grid2D {
  int sizex1, sizex2;        // number of elements in the grid 
  DOUBLE stepx1, stepx2;     // grid spacing dx
  DOUBLE I_stepx1, I_stepx2; // 1/dx
  DOUBLE minx1, minx2;       // min
  DOUBLE maxx1, maxx2;       // max
  DOUBLE maxx12, maxx22;     // max2
  DOUBLE *x1, *x2;           // axis
  DOUBLE **fp;               // f'/f  
  DOUBLE **lnf;              // ln(f) 
  DOUBLE **E;                // (f"+2f'/r)/f - (f'/f)^2
};

extern DOUBLE Uwell; // depth of the well

void LoadTrialWaveFunction2D(struct Grid2D *G);
DOUBLE Interpolate2DU(struct Grid2D *Grid2D, DOUBLE x1, DOUBLE x2);

#endif
