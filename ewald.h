/*ewald.h*/

#ifndef _EWALD_H
#define _EWALD_H

#include "compatab.h"

DOUBLE EwaldSumWalker(struct Walker *Walker);
int EwaldCheck(void);
void EwaldInit(void);
DOUBLE EwaldXi(void);
DOUBLE EwaldXi2(void);
DOUBLE EwaldSumAllWalkers(void);
DOUBLE InteractionSelfEnergy(void);
DOUBLE EwaldPsi(DOUBLE x, DOUBLE y, DOUBLE z);
DOUBLE EwaldOptimize(int Nx, int Nk);
void EwaldOptimizeInitExact(void);
void EwaldOptimizeSigma(DOUBLE *E, DOUBLE *sigma);
void EwaldInitCoef(void);

extern int Ewald_Nx;
extern int Ewald_Nk;
extern DOUBLE Ewald_Xi;
extern DOUBLE Ewald_L,Ewald_Lk;

//DOUBLE EwaldSumPair(DOUBLE x, DOUBLE y, DOUBLE z);
#define EwaldSumPair(x,y,z) (EwaldPsi(x*Ewald_L, y*Ewald_L, z*Ewald_L) + Ewald_Xi/(double) (N-1))*Ewald_Lk
//INLINE DOUBLE EwaldSumPair(DOUBLE x, DOUBLE y, DOUBLE z) { return (EwaldPsi(x*Ewald_L, y*Ewald_L, z*Ewald_L) + Ewald_Xi/(double)(N-1))*Ewald_L3;}

int SmoothCutoffInit(void);
double SmoothCutoffSum(double x, double y, double z);

#endif
  