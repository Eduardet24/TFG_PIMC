/*pigs.h*/
#ifndef _PIGS_H_
#define _PIGS_H_

#define VextUp(x, y, z, i) Vext(x, y, z)
#define VextDn(x, y, z, i) Vext(x, y, z)

// one/two components
void PathIntegralMoveAll(void);
DOUBLE PathIntegralPseudopotentialSdetailed(DOUBLE *Skincheck, DOUBLE *Sdeltacheck, DOUBLE *Swfcheck);
DOUBLE PathIntegralPseudopotentialS(void);
double PropagatorPseudopotenial1D_ij(double m, double a, double dt, double xi0, double xj0, double xi, double xj);
double PropagatorLaplacian(double m, double dt, double dr);
void PathIntegralPseudopotentialStaging(void);

#ifdef ONE_COMPONENT_CODE
  extern DOUBLE mA;
#endif

#endif
