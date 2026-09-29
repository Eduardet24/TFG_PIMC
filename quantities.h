/*quantities.h*/

#ifndef _QUANTITIES_H_
#define _QUANTITIES_H_

#include "main.h"

#ifdef INTERACTION_LENNARD_JONES

// T_lj = h^2/2m A^2
// U_lj = 4 eps_lj / T_lj 
// H2
/*
#  define eps_lj 10.7 //= 42.8/4
#  define sigma_lj 2.97
#  define T_lj 24.0668
#  define U_lj 1.77839 */

//4He
//#  define eps_lj 10.22
//#  define sigma_lj 2.556
//#  define T_lj 12.115922801724114 // (2*6.0599278)
//#  define U_lj 3.374072340093048 // (4*eps_lj/T_lj)

// LJ units
#  define eps_lj 1.
#  define sigma_lj 1.
#  define T_lj 1.
//#  define U_lj 4.
#  define U_lj 1.

//#  define energy_unit T_lj
#endif

DOUBLE Energy(DOUBLE *Epot, DOUBLE *Ekin, DOUBLE *Eff, DOUBLE *Edamping, DOUBLE *Eint, DOUBLE *Eext);
DOUBLE WalkerEnergy(struct Walker *walker, DOUBLE *Epot, DOUBLE *Ekin, DOUBLE *EFF, DOUBLE *Edamping, DOUBLE *Eint, DOUBLE *Eext);
DOUBLE WalkerEnergy0(struct Walker *W);
DOUBLE Vext(DOUBLE x, DOUBLE y, DOUBLE z);
DOUBLE VextTilde(DOUBLE x, DOUBLE y, DOUBLE z); // [V, [T, V]] , needed for quartic order
DOUBLE InteractionEnergy(DOUBLE r);
DOUBLE InteractionEnergy_ij(DOUBLE r, int spin1, int spin2); // taking in account spin
DOUBLE InteractionEnergyWithImages_ij(DOUBLE x, DOUBLE y, DOUBLE z);
DOUBLE InteractionEnergyTilde(DOUBLE r);  // [V, [T, V]] , needed for quartic order
DOUBLE InteractionEnergySumOverImages(DOUBLE x, DOUBLE y, DOUBLE z);
DOUBLE InteractionDamping(DOUBLE r);
DOUBLE WalkerPotentialEnergy(struct Walker *W);
DOUBLE WalkerPotentialEnergyTilde(struct Walker *W);
void MeasureHessianMatrix(void);

int PermutationParity(DOUBLE *r);

int MeasureOBDM(int w);
int MeasureOBDMpure(void);
int MeasureOBDMShortRange(void);
void MeasureOBDMMatrix(void);
int MeasureTBDM(void);
int MeasureTBDMpure(void);
void MeasureSuperfluidDensity(void);
void MeasurePairDistribution(void);
void MeasureRadialDistribution(void);
void MeasureRadialDistributionWalker(struct Walker *W);
void EmptyDistributionArrays(void);
void MeasureStaticStructureFactor(int w);
void MeasureOrderParameter(void);
void MeasureLindemannRatio(void);
void InitializeFormFactor(void);
void MeasureFormFactor(int w);
int MeasureOBDMorbital(void);
void MeasureWaveFunctionProjection(void);
void MeasureHyperRadius(void);
void MeasureEffectivePotential(void);
void MeasurePurePotentialEnergy(void);
DOUBLE MeasureKineticEnergyNumerically(struct Walker *walker);
DOUBLE SpinPolarization(void);
int MeasureOBDMfromU(void);
int MeasureOBDMfromUWalker(int w);
int MeasureTBDMfromU(void);
int MeasureTBDMfromUWalker(int w);

struct sOBDM {
  int size;
  DOUBLE step;
  DOUBLE min;
  DOUBLE max;
  DOUBLE *f;
  DOUBLE *N;
  int Nksize; // parameters of momentum distribution are stored here
  DOUBLE kmin; 
  DOUBLE kmax; 
  DOUBLE *k;  // grid for nk
  DOUBLE *Nk; // Nk
  int times_measured; // needed for Nk
  int Nktimes_measured; // needed for Nk
  DOUBLE CF; // condensate fraction calculated from projection to an external orbital
  int CFN;
};

struct sMATRIX {
  int size;
  DOUBLE step;
  DOUBLE step_x;
  DOUBLE step_y;
  DOUBLE min;
  DOUBLE max;
  DOUBLE maxx;
  DOUBLE maxy;
  DOUBLE **f;
  int **N;
  DOUBLE Norma;
  int times_measured;
};

struct sSk {
  int size; // length of input array
  unsigned long int *N;        // times measured
  DOUBLE *k;     // k^2 = kx^2+ky^2+kz^2
  DOUBLE *kx;
  DOUBLE *ky;
  DOUBLE *kz;
  int *index;    // number of non-degenerate values of momenta
  DOUBLE *f;     // same spin
  DOUBLE L;      // maximal distance
  DOUBLE *cos, *sin; // <cos(kr)>; <sin(kr)>, important for optical lattice
  int *degeneracy; // 1-no degeneracy, set d>1 for the first degerate momentum and d=0 to other 
  unsigned long int times_measured;
};

struct Distribution {
  int *N;
  DOUBLE *f; // ready to be plotted !
  //double times_measured;
  int times_measured;
  int size;
  DOUBLE step;
  DOUBLE min;
  DOUBLE max;
  DOUBLE r2;
  DOUBLE r2recent;
  DOUBLE r4;
  DOUBLE width;
  int g3;
#ifdef SPINFULL
   int **PDSpin;
#endif 
};

struct sSD {
  DOUBLE *CM2;
  int *N;
  int size;
  int spacing;
  int times_measured;
};

struct sOrderParameter {
  DOUBLE cos;
  DOUBLE sin;
  int times_measured;
  DOUBLE cos_pure;
  DOUBLE sin_pure;
  int times_measured_pure;
};

struct sLindemannRatio {
  int N;    // times measured
  DOUBLE F; // quantity measured (mixed)
  int *Npure; // length grid_pure_block
  DOUBLE *Fpure;
};

struct sPureEstimator {
  unsigned long int size;
  DOUBLE step;
  DOUBLE min;
  DOUBLE max;
  DOUBLE *f;
  DOUBLE *x;
  int times_measured;
};

#endif
