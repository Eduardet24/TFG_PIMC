/*crystal.h*/

#ifndef _CRYSTAL_H_
#define _CRYSTAL_H_

#include "main.h"

#ifndef DOUBLE
#  define DOUBLE double
#endif

extern DOUBLE crystal_side_x;
extern DOUBLE crystal_side_y;

#ifdef ONE_BODY_TRIAL_TERMS

struct Cryst {
  DOUBLE *x; // coordinates of the cells
  DOUBLE *y; 
  DOUBLE *z; 
  int size; // number of grid cells
  // Solid phase
#ifdef CRYSTAL_WIDTH_ARRAY
  DOUBLE *Rx; // Gaussian width in (x,y,z) directions
  DOUBLE *Ry; // 
  DOUBLE *Rz; // 
  DOUBLE *weight; // 
#else
  DOUBLE Rx; // Gaussian width in (x,y,z) directions
  DOUBLE Ry; // 
  DOUBLE Rz; // 
  DOUBLE weight; // for symmetrized w.f.
#endif
  DOUBLE Rx_av; // (sum_{i} Rx[i]) / NCrystal
  DOUBLE Ry_av; // 
  DOUBLE Rz_av; // 
  // Impurities
  DOUBLE b; // scattering length
  DOUBLE k; // momentum (fixed by Rpar)
  DOUBLE k2;
  // Non-orthogonal basis
  DOUBLE e1x;
  DOUBLE e1y;
  DOUBLE e1abs_inv2; // basis: e1 = (e1x,e1y) 
  DOUBLE e2x;
  DOUBLE e2y;
  DOUBLE e2abs_inv2; //        e2 = (e2x,e2y)
};

#define CRYSTAL_Zi Crystal.z[i]

extern struct Cryst Crystal;

void ConstructGridCrystal(struct Cryst *Crystal, int size, DOUBLE n);
void ConstructGridImpurity(struct Cryst *Crystal, const DOUBLE R, int size);
void ConstructGridCrystalOneBody(void);

#ifdef ONE_BODY_IMPURITY
DOUBLE ImpurityInterpolateU(DOUBLE x);
DOUBLE ImpurityInterpolateFp(DOUBLE x);
DOUBLE ImpurityInterpolateE(DOUBLE x);
#endif

#ifdef CHAINS_PINNED_TO_CENTER_OF_MASS
extern struct Walker *walker;
int ChainsPinnedCenterMass(void);
void CalculateAllChainsCenterOfMass(DOUBLE** R, int w);
void CalculateAllChainsCenterOfMassWalker(struct Walker *walker);
void CalculateOneChainCenterOfMass(DOUBLE x, DOUBLE y, DOUBLE z, int i, int w);
#endif

#endif

#endif
