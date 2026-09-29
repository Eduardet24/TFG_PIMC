/*trial.h*/

#ifndef __TRIAL_H_
#define __TRIAL_H_

#include "main.h"

extern DOUBLE trial_Vo, trial_Kappa, trial_Kappa2, trial_k, trial_k2;
extern DOUBLE A1SW, A2SW, A3SW, A4SW;

struct Grid {
  long int size; //number of elements in the grid
  DOUBLE step;
  DOUBLE I_step;
  DOUBLE min;
  DOUBLE max;
  DOUBLE max2;
  DOUBLE *x; 
  DOUBLE *fp;  // f'/f
  DOUBLE *lnf; // ln(f)
  DOUBLE *E;   // (f"+(D-1)f'/r)/f - (f'/f)^2
  DOUBLE *f;   // f
  DOUBLE *fpp; // f''
//#ifdef INTERACTION_WITH_DAMPING
  DOUBLE *ReV, *ImV; // Real and imaginary part of interaction potential
  DOUBLE k; // sqrt(|E|)
//#endif
};

void LoadTrialWaveFunction(struct Grid *G, char *file_wf, int normalize_to_one);
void LoadInteractionPotential(struct Grid *G, char *file_wf);
void SampleTrialWaveFunction(struct Grid *G,  DOUBLE min, DOUBLE max, long int size, DOUBLE (*U)(DOUBLE x), DOUBLE (*Fp)(DOUBLE x), DOUBLE (*E)(DOUBLE x));
void CheckTrialWaveFunctionInterpolation(struct Grid *G, DOUBLE (*U)(DOUBLE x), DOUBLE (*Fp)(DOUBLE x), DOUBLE (*E)(DOUBLE x));
void CheckTrialWaveFunction(struct Grid *G);

DOUBLE ConstructSoftSphereHeight(DOUBLE Ro, DOUBLE a);
DOUBLE ConstructHardSphereSquareWellDepth(DOUBLE a, DOUBLE Ro, DOUBLE Apar);
void ConstructGridLim(struct Grid *G, const DOUBLE xmin, const DOUBLE xmax, const long int size);
void ConstructGridPower(struct Grid *G, const DOUBLE xmin, const DOUBLE xmax, const long int size);
void ConstructGridHS(struct Grid *G);
void ConstructGridHSSimple(struct Grid *G);
void ConstructGridPseudopotentialSimple(struct Grid *G);
void ConstructGridHS1D(struct Grid *G, const DOUBLE R, const DOUBLE xi, const DOUBLE D, const int size, const DOUBLE Rmin, const DOUBLE Rmax);
void ConstructGridFree(struct Grid *G, const DOUBLE R, const DOUBLE Bpar, const DOUBLE a, const int size, const DOUBLE Rmin, const DOUBLE Rmax);
void ConstructGridLieb(struct Grid *G, const int size);
void ConstructGridSS(struct Grid *G, const DOUBLE a, const DOUBLE R,const DOUBLE Rpar, const DOUBLE Bpar, const int size, const DOUBLE Rmin, const DOUBLE Rmax);
void ConstructGridSchiffVerlet(struct Grid *G, const DOUBLE Apar, const DOUBLE Bpar,const DOUBLE n);
void ConstructGridShortRange3D(struct Grid *G, const DOUBLE a, const DOUBLE Rpar);
void ConstructGridSutherland(struct Grid *G);
void ConstructGridCalogero(struct Grid *G);
void ConstructGridYukawa2D(struct Grid *G, const DOUBLE Apar, const DOUBLE Bpar, const DOUBLE n);
void ConstructGridYukawa3D(struct Grid *G, const DOUBLE Apar, const DOUBLE Bpar, const DOUBLE Cpar);

void ConstructGridPowerKo(struct Grid *G);
void ConstructGridDipoleKo(struct Grid *G);
void ConstructGridQuadrupoleKo(struct Grid *G);
void ConstructGridPhonons(struct Grid *G);
void ConstructGridDipolePhonon(struct Grid *G);
void ConstructGridDecayPhonon(struct Grid *G, const DOUBLE Apar, const DOUBLE Bpar);
void ConstructGridCalogero2D(struct Grid *G, const DOUBLE Apar);
void ConstructGridCoulomb1DPhonon(struct Grid *G);
void ConstructGridDipole1DPhonon(struct Grid *G);
void ConstructGridCoulomb1DTrap(struct Grid *G);
void ConstructGridCoulomb2DTrap(struct Grid *G);
void ConstructGridExpPhonon(struct Grid *G);
void ConstructGridHS2DFiniteEnergy(struct Grid *G);
void ConstructGridDipolePhononFiniteD(struct Grid *G);
void ConstructGridHardRods3D(struct Grid *G);
void ConstructSqWellDepth(void);
void ConstructGridSquareWellZeroConst(void);
void ConstructGridBCSSquareWellZeroEnergy(void);
void ConstructGridSSlarge(struct Grid *G);
void ConstructGridSSExp(struct Grid *G);
void ConstructGridSquareWellBoundStatePhonons(struct Grid *G);
void ConstructGridSquareWellFreeStatePhonons(struct Grid *G);
void ConstructGridSquareWellBoundStateTrap(struct Grid *G);
void ConstructGridDipoleKo11(struct Grid *G);
void ConstructGridLiebLuttingerPhonon(struct Grid *G, DOUBLE K);
void ConstructGridPhononLuttingerAuto(struct Grid *G);
void ConstructGridPhononLuttingerAuto11(struct Grid *G);
void ConstructGridPhononLuttingerAuto12(struct Grid *G);
void ConstructGridLiebExponentialDecay(struct Grid *G);
void ConstructGridPhononLuttinger(struct Grid *G);
void ConstructGridZeroRangeUnitaryPhonons(struct Grid *G);
void ConstructGridCSM3DPlasmon(struct Grid *G, DOUBLE Apar);
void ConstructVextSqWellDepth(void);
void ConstructGridSquareWellTrap(void);
void ConstructGridSquareWellTrapMcGuire1D(void);
void ConstructGridCoulombPhonon(struct Grid *G);
void ConstructGrid3BodyZeroEnergyPhonons(void);
void ConstructGridMcGuirePhonon(void);
void ConstructGridSoftSpherePhonon(struct Grid *G);
void ConstructGridSoftSpherePhononAuto(struct Grid *G);
void ConstructGridZeroRange3D(struct Grid *G);

DOUBLE InterpolateU(struct Grid *Grid, DOUBLE r, int i, int j);
DOUBLE InterpolateFp(struct Grid *Grid, DOUBLE r, int i, int j);
DOUBLE InterpolateE(struct Grid *Grid, DOUBLE r, int i, int j);
DOUBLE InterpolateF(struct Grid *Grid, DOUBLE r);

DOUBLE InterpolateExactU(DOUBLE r);
DOUBLE InterpolateExactFp(DOUBLE r);
DOUBLE InterpolateExactE(DOUBLE r);

DOUBLE InterpolateExactU11(DOUBLE r);
DOUBLE InterpolateExactFp11(DOUBLE r);
DOUBLE InterpolateExactE11(DOUBLE r);

DOUBLE InterpolateGridU(struct Grid *Grid, DOUBLE r);
DOUBLE InterpolateGridFp(struct Grid *Grid, DOUBLE r);
DOUBLE InterpolateGridE(struct Grid *Grid, DOUBLE r);

DOUBLE InterpolateGridMatchedToScatteringSolutionU(struct Grid *Grid, DOUBLE r);
DOUBLE InterpolateGridMatchedToScatteringSolutionFp(struct Grid *Grid, DOUBLE r);
DOUBLE InterpolateGridMatchedToScatteringSolutionE(struct Grid *Grid, DOUBLE r);

DOUBLE InterpolateEG(DOUBLE r);
DOUBLE IntegrateSimpson(DOUBLE xmin, DOUBLE xmax, DOUBLE (*f)(DOUBLE x), long int Npoints);

/*DOUBLE InterpolateDerivU(struct Grid *Grid, DOUBLE r);
void InterpolateDerivUError(struct Grid *Grid, DOUBLE *Erel, DOUBLE *Eabs, DOUBLE *rrel, DOUBLE *rabs);
DOUBLE InterpolateDerivFp(struct Grid *Grid, DOUBLE r);
void InterpolateDerivFpError(struct Grid *Grid, DOUBLE *Erel, DOUBLE *Eabs, DOUBLE *rrel, DOUBLE *rabs);*/

extern DOUBLE sqE, sqE2, expBpar, alpha_1, Xi, Atrial, Btrial, alpha2E, sqEkappa2;

void ConstructGridSpline2Body(struct Grid *G, DOUBLE (*InteractionEnergy)(DOUBLE x));
void ConstructGridSpline2Body_ij(struct Grid *G, int spin1, int spin2);
void ConstructGridSpline2BodyPhonons(struct Grid *G);
DOUBLE SolveScatteringProblem(struct Grid *G, DOUBLE (*InteractionEnergy)(DOUBLE x), DOUBLE E, int flag_boundary_condition_right, int flag_boundary_condition_left, int *function_was_negative, int flag_lnf_instead_of_f);
DOUBLE SolveScatteringProblem_ij(struct Grid *G, int spin1, int spin2, DOUBLE E, int flag_boundary_condition_right, int flag_boundary_condition_left, int *function_was_negative, int flag_lnf_instead_of_f);
void FindSWBindingEnergy(void);

DOUBLE u3(DOUBLE R); // three-body terms
DOUBLE Fp3(DOUBLE R);
DOUBLE E3(DOUBLE R);

#ifdef SPINFULL
// functions are defined in "trial.c"
#else // spinless
#ifdef INTERPOLATE_LINEAR_JASTROW_WF
#ifdef JASTROW_RIGHT_BOUNDARY_BOUND_STATE
#  define InterpolateU(G,r,i,j) InterpolateGridMatchedToScatteringSolutionU(G,r)
#  define InterpolateFp(G,r,i,j) InterpolateGridMatchedToScatteringSolutionFp(G,r)
#  define InterpolateE(G,r,i,j) InterpolateGridMatchedToScatteringSolutionE(G,r)
#else
#  define InterpolateU(G,r,i,j) InterpolateGridU(G,r)
#  define InterpolateFp(G,r,i,j) InterpolateGridFp(G,r)
#  define InterpolateE(G,r,i,j) InterpolateGridE(G,r)
#endif
#else // no linear interp.
#ifdef INTERPOLATE_SPLINE_JASTROW_WF
#ifdef JASTROW_RIGHT_BOUNDARY_PHONONS 
// defined in trial.c
#else
#    define InterpolateU(G,r,i,j)  InterpolateSplineU(G,r)
#    define InterpolateFp(G,r,i,j) InterpolateSplineFp(G,r)
#    define InterpolateE(G,r,i,j)  InterpolateSplineE(G,r)
#endif
#    define orbitalF orbitalExactF
#    define orbitalFp orbitalExactFp
#    define orbitalEloc orbitalExactEloc
#else //no linear, no spline interp.
#ifdef TRIAL_KURBAKOV_NUMERIC
#      define InterpolateU(G,r,i,j)  InterpolateKurbakovU(r)
#      define InterpolateFp(G,r,i,j) InterpolateKurbakovFp(r)
#      define InterpolateE(G,r,i,j)  InterpolateKurbakovE(r)
#else
#ifndef SYMMETRIZE_TRIAL_WF
#ifndef POWER_TRIAL_WF
#        define InterpolateU(G,r,i,j)  InterpolateExactU(r)
#        define InterpolateFp(G,r,i,j) InterpolateExactFp(r)
#        define InterpolateE(G,r,i,j)  InterpolateExactE(r)
#endif // POWER_TRIAL_WF
#endif // SYMMETRIZE_TRIAL_WF
#endif // TRIAL_KURBAKOV_NUMERIC
#endif // INTERPOLATE_SPLINE_JASTROW_WF
#endif // INTERPOLATE_LINEAR_JASTROW_WF
#endif // SPINFULL


#endif
