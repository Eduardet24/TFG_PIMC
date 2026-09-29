/*trial.c*/
/* best trial wavefunction Psi_T = Prod_{i<j} f(r_ij)
   U = ln f
   Fp = f' / f
3D Eloc = [-(f"+2f'/r)/f] + (f'/f)^2
2D Eloc = [-(f" + f'/r)/f] + (f'/f)^2
1D Eloc = [-f"/f] + (f'/f)^2

  3body
1D E3(R) = -f''(R)/f(R) - 1/R f'(R)/f(R) + (f'(R)/f(R))^2

1)        {0                   (r<a)
   f(r) = {A Sin(sqE*(r-D))/r  (a<r<R)
          {1-B*Exp(-r/alpha)   (R<r)
   R <-> Rpar,  xi <-> Bpar, u <-> x

2)

3) f(r) = A Sin(sqE*(r-a))    (r<R)
   f(r) = 1-B*Exp(-r/alpha)   (r>R)

SS)        [ A Sin(sqE*r)/r
          {[ A sh(sqE*r)/r,        (r<Atrial)
   f(r) = { B Sin(sqE'*r+delta)/r, (Atrial<r<Rpar)
          { 1- C Exp(-r/alpha),    (Rpar<r)
*/
#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include "trial.h"
#include "main.h"
#include "utils.h"
#include "memory.h"
#include "spline.h"
#include "compatab.h"
#include "mymath.h"
#include MATHINCLUDE

#define DOUBLE double

DOUBLE sqE, sqE2, expBpar, alpha_1, Xi, Atrial, Btrial, Ctrial, Dtrial, Etrial, Ktrial, alpha2E;
//SS  
DOUBLE sqEkappa, sqEkappa2, trial_delta;
DOUBLE trial_Vo, trial_Kappa, trial_Kappa2, trial_k, trial_k2;
//DOUBLE A1SW, A2SW, A3SW, A4SW;

int high_energy;

DOUBLE A1SW, A2SW, A3SW, A4SW, A5SW;
DOUBLE A1SD, A2SD, A3SD, A4SD;
DOUBLE alpha11_1, alpha12_1;
DOUBLE Atrial12, Btrial12, Ctrial12;
DOUBLE Atrial11, Btrial11, Ctrial11;
DOUBLE Atrial22, Btrial22, Ctrial22;
DOUBLE sqE11, sqE12, sqE211, sqE212;

#ifdef SPINFULL
#ifdef TRIAL_TWO_BODY_SCATTERING_NUMERICAL_SPINFULL
DOUBLE InterpolateU(struct Grid *Grid, DOUBLE r, int spin1, int spin2) {
#ifdef JASTROW_RIGHT_BOUNDARY_BOUND_STATE
  int bound_state_asymptotics = OFF;
#endif

#ifdef SAME_SPIN_CONSTANT_JASTROW
  if(spin1 == spin2) return 0.;
#endif

#ifdef JASTROW_RIGHT_BOUNDARY_BOUND_STATE
#ifdef JASTROW_RIGHT_BOUNDARY_SPINFULL_SAME_SPIN_BOUND_STATE
  if(spin1 == spin2 && r>GridSpin[spin1][spin2].max-10*GridSpin[spin1][spin2].step) bound_state_asymptotics = ON;
#endif
#ifdef JASTROW_RIGHT_BOUNDARY_SPINFULL_DIFFERENT_SPIN_BOUND_STATE
  if(spin1 != spin2 && r>GridSpin[spin1][spin2].max-10*GridSpin[spin1][spin2].step) bound_state_asymptotics = ON;
#endif
  if(bound_state_asymptotics) {
#ifdef TRIAL_1D // f(r) = exp(-kr)
    return -GridSpin[spin1][spin2].k*(r-GridSpin[spin1][spin2].max);
#endif
#ifdef TRIAL_2D // f(r) = Ko(kr)
    return log(K0(GridSpin[spin1][spin2].k*r)) - log(K0(GridSpin[spin1][spin2].k*GridSpin[spin1][spin2].max));
#endif
#ifdef TRIAL_3D // f(r) = exp(-kr)/r
    return -GridSpin[spin1][spin2].k*(r-GridSpin[spin1][spin2].max) - log(r/GridSpin[spin1][spin2].max);
#endif
  }
#endif

  return InterpolateSplineU(&GridSpin[spin1][spin2], r);
}

DOUBLE InterpolateFp(struct Grid *Grid, DOUBLE r, int spin1, int spin2) {
#ifdef JASTROW_RIGHT_BOUNDARY_BOUND_STATE
  int bound_state_asymptotics = OFF;
#endif

#ifdef SAME_SPIN_CONSTANT_JASTROW
  if(spin1 == spin2) return 0;
#endif

#ifdef JASTROW_RIGHT_BOUNDARY_BOUND_STATE
#ifdef JASTROW_RIGHT_BOUNDARY_SPINFULL_SAME_SPIN_BOUND_STATE
  if(spin1 == spin2 && r>GridSpin[spin1][spin2].max-10*GridSpin[spin1][spin2].step) bound_state_asymptotics = ON;
#endif
#ifdef JASTROW_RIGHT_BOUNDARY_SPINFULL_DIFFERENT_SPIN_BOUND_STATE
  if(spin1 != spin2 && r>GridSpin[spin1][spin2].max-10*GridSpin[spin1][spin2].step) bound_state_asymptotics = ON;
#endif
  if(bound_state_asymptotics) {
#ifdef TRIAL_1D // f(r) = exp(-kr)
    return -GridSpin[spin1][spin2].k;
#endif
#ifdef TRIAL_2D // f(r) = Ko(kr)
    return -GridSpin[spin1][spin2].k*K1(GridSpin[spin1][spin2].k * r)/K0(GridSpin[spin1][spin2].k * r);
#endif
#ifdef TRIAL_3D // f(r) = exp(-kr)/r
    return -GridSpin[spin1][spin2].k - 1./r;
#endif
  }
#endif

  return InterpolateSplineFp(&GridSpin[spin1][spin2], r);
}

DOUBLE InterpolateE(struct Grid *Grid, DOUBLE r, int spin1, int spin2) {
#ifdef JASTROW_RIGHT_BOUNDARY_BOUND_STATE
  int bound_state_asymptotics = OFF;
#endif

#ifdef SAME_SPIN_CONSTANT_JASTROW
  if(spin1 == spin2) return 0.;
#endif

#ifdef JASTROW_RIGHT_BOUNDARY_BOUND_STATE
#ifdef JASTROW_RIGHT_BOUNDARY_SPINFULL_SAME_SPIN_BOUND_STATE
  if(spin1 == spin2 && r>GridSpin[spin1][spin2].max-10*GridSpin[spin1][spin2].step) bound_state_asymptotics = ON;
#endif
#ifdef JASTROW_RIGHT_BOUNDARY_SPINFULL_DIFFERENT_SPIN_BOUND_STATE
  if(spin1 != spin2 && r>GridSpin[spin1][spin2].max-10*GridSpin[spin1][spin2].step) bound_state_asymptotics = ON;
#endif
  if(bound_state_asymptotics) {
#ifdef TRIAL_1D // f(r) = exp(-kr)
    return 0.;
#endif
#ifdef TRIAL_2D // f(r) = Ko(kr)
    DOUBLE Fp;
    Fp = K1(GridSpin[spin1][spin2].k * r)/K0(GridSpin[spin1][spin2].k * r);
    return GridSpin[spin1][spin2].k*GridSpin[spin1][spin2].k*(-1. + Fp*Fp);
#endif
#ifdef TRIAL_3D // f(r) = exp(-kr)/r
    return (1.+2.*GridSpin[spin1][spin2].k*r)/(r*r);
#endif
  }
#endif

  return InterpolateSplineE(&GridSpin[spin1][spin2], r);
}
#else // SPINFULL but not TRIAL_TWO_BODY_SCATTERING_NUMERICAL_SPINFULL
DOUBLE InterpolateU(struct Grid *Grid, DOUBLE r, int spin1, int spin2) {
  if(spin1 == spin2)
    return InterpolateExactU11(r);
  else {
#ifdef JASTROW_RIGHT_BOUNDARY_BOUND_STATE
    int bound_state_asymptotics = OFF;
#ifdef JASTROW_RIGHT_BOUNDARY_SPINFULL_SAME_SPIN_BOUND_STATE
    if(spin1 == spin2 && r>Grid->max-10*Grid->step) bound_state_asymptotics = ON;
#endif
#ifdef JASTROW_RIGHT_BOUNDARY_SPINFULL_DIFFERENT_SPIN_BOUND_STATE
    if(spin1 != spin2 && r>Grid->max-10*Grid->step) bound_state_asymptotics = ON;
#endif
    if(bound_state_asymptotics) {
#ifdef TRIAL_1D // f(r) = exp(-kr)
      return -Grid->k*(r-Grid->max);
#endif
#ifdef TRIAL_2D // f(r) = Ko(kr)
      return log(K0(Grid->k*r)) - log(K0(Grid->k*Grid->max));
#endif
#ifdef TRIAL_3D // f(r) = exp(-kr)/r
      return -Grid->k*(r-Grid->max) - log(r/Grid->max);
#endif
  }
#endif

#ifdef INTERPOLATE_SPLINE_JASTROW_WF
    return InterpolateSplineU(Grid,r);
#else
    return InterpolateExactU(r);
#endif
  }
}

DOUBLE InterpolateFp(struct Grid *Grid, DOUBLE r, int spin1, int spin2) {
  if(spin1 == spin2)
    return InterpolateExactFp11(r);
  else {
#ifdef JASTROW_RIGHT_BOUNDARY_BOUND_STATE
  int bound_state_asymptotics = OFF;
#ifdef JASTROW_RIGHT_BOUNDARY_SPINFULL_SAME_SPIN_BOUND_STATE
  if(spin1 == spin2 && r>Grid->max-10*Grid->step) bound_state_asymptotics = ON;
#endif
#ifdef JASTROW_RIGHT_BOUNDARY_SPINFULL_DIFFERENT_SPIN_BOUND_STATE
  if(spin1 != spin2 && r>Grid->max-10*Grid->step) bound_state_asymptotics = ON;
#endif
  if(bound_state_asymptotics) {
#ifdef TRIAL_1D // f(r) = exp(-kr)
    return -Grid->k;
#endif
#ifdef TRIAL_2D // f(r) = Ko(kr)
    return -Grid->k*K1(Grid->k * r)/K0(Grid->k * r);
#endif
#ifdef TRIAL_3D // f(r) = exp(-kr)/r
    return -Grid->k - 1./r;
#endif
  }
#endif

#ifdef INTERPOLATE_SPLINE_JASTROW_WF
    return InterpolateSplineFp(Grid,r);
#else
    return InterpolateExactFp(r);
#endif
  }
}

DOUBLE InterpolateE(struct Grid *Grid, DOUBLE r, int spin1, int spin2) {
  if(spin1 == spin2)
    return InterpolateExactE11(r);
  else {
#ifdef JASTROW_RIGHT_BOUNDARY_BOUND_STATE
  int bound_state_asymptotics = OFF;
#ifdef JASTROW_RIGHT_BOUNDARY_SPINFULL_SAME_SPIN_BOUND_STATE
  if(spin1 == spin2 && r>Grid->max-10*Grid->step) bound_state_asymptotics = ON;
#endif
#ifdef JASTROW_RIGHT_BOUNDARY_SPINFULL_DIFFERENT_SPIN_BOUND_STATE
  if(spin1 != spin2 && r>Grid->max-10*Grid->step) bound_state_asymptotics = ON;
#endif
  if(bound_state_asymptotics) {
#ifdef TRIAL_1D // f(r) = exp(-kr)
    return 0.;
#endif
#ifdef TRIAL_2D // f(r) = Ko(kr)
    DOUBLE Fp;
    Fp = K1(Grid->k * r)/K0(Grid->k * r);
    return Grid->k*Grid->k*(-1. + Fp*Fp);
#endif
#ifdef TRIAL_3D // f(r) = exp(-kr)/r
    return (1.+2.*Grid->k*r)/(r*r);
#endif
  }
#endif

#ifdef INTERPOLATE_SPLINE_JASTROW_WF
    return InterpolateSplineE(Grid,r);
#else
    return InterpolateExactE(r);
#endif
  }
}
#endif
#else // i.e. spinless

/************************** W F Symmetrization *****************************/
#ifdef SYMMETRIZE_TRIAL_WF
DOUBLE InterpolateU(struct Grid *Grid, DOUBLE r, int i, int j) {
  // u(r) + u(L-r)
  return ((r>Lhalfwf)?(0):(InterpolateExactU(r)+InterpolateExactU(Lwf-r)-2.*InterpolateExactU(Lhalfwf)));
}
DOUBLE InterpolateFp(struct Grid *Grid, DOUBLE r, int i, int j) {
  // fp(r) - fp(L-r)
  return ((r>Lhalfwf)?(0):(InterpolateExactFp(r)-InterpolateExactFp(Lwf-r)));
}
DOUBLE InterpolateE(struct Grid *Grid, DOUBLE r, int i, int j) {
#ifdef TRIAL_1D
  return (r>Lhalfwf)?(0):(InterpolateExactE(r) + InterpolateExactE(Lwf-r));
#endif
#ifdef TRIAL_2D
  return (r>Lhalfwf)?(0):(InterpolateExactE(r) + InterpolateExactE(Lwf-r) + 1./r*InterpolateExactFp(Lwf-r)*Lwf/(Lwf-r));
#endif
#ifdef TRIAL_3D
  return (r>Lhalfwf)?(0):(InterpolateExactE(r) + InterpolateExactE(Lwf-r) + 2./r*InterpolateExactFp(Lwf-r)*Lwf/(Lwf-r));
#endif
}
#endif// SYMMETRIZE_TRIAL_WF

#ifdef POWER_TRIAL_WF
DOUBLE InterpolateU(struct Grid *Grid, DOUBLE r, int i, int j) {
  // alpha_Jastrow u(r)
  return alpha_Jastrow*InterpolateExactU(r);
}
DOUBLE InterpolateFp(struct Grid *Grid, DOUBLE r, int i, int j) {
  return alpha_Jastrow*InterpolateExactFp(r);
}
DOUBLE InterpolateE(struct Grid *Grid, DOUBLE r, int i, int j) {
  return alpha_Jastrow*InterpolateExactE(r);
}
#endif

#ifdef JASTROW_RIGHT_BOUNDARY_PHONONS
#ifdef INTERPOLATE_SPLINE_JASTROW_WF
DOUBLE InterpolateU(struct Grid *Grid, DOUBLE r, int spin1, int spin2) {
  if(r>Grid->max-10*Grid->step) { 
#ifdef TRIAL_3D
    return - Apar/(r*r) - Apar/((L-r)*(L-r)) + 2.*Apar/Lhalf2;
#else
    return - Apar/r - Apar/(L-r) + 2.*Apar/Lhalf;
#endif
  }
  else {
    return InterpolateSplineU(Grid, r);
  }
}

DOUBLE InterpolateFp(struct Grid *Grid, DOUBLE r, int spin1, int spin2) {
  if(r>Grid->max-10*Grid->step) { 
#ifdef TRIAL_3D
    return 2.*Apar*(1./(r*r*r)-1./((L-r)*(L-r)*(L-r)));
#else
    return Apar*(1./(r*r)-1./((L-r)*(L-r)));
#endif
  }
  else {
    return InterpolateSplineFp(Grid, r);
  }
}

DOUBLE InterpolateE(struct Grid *Grid, DOUBLE r, int spin1, int spin2) {
  if(r>Grid->max-10*Grid->step) { 
#ifdef TRIAL_3D
    return (2.*Apar*(L*L*L*L - 4.*L*L*L*r + 6.*L*L*r*r - 2.*L*r*r*r + 2.*r*r*r*r))/((L-r)*(L-r)*(L-r)*(L-r)*r*r*r*r);
#else
    return (Apar*L*(L*L - 3*L*r + 4*r*r))/((L - r)*(L - r)*(L - r)*r*r*r);//??
#endif
  }
  else {
    return InterpolateSplineE(Grid, r);
  }
}
#endif
#endif
#endif// not SPINFULL

/************************** Check Trial Wave Function **********************/
void CheckTrialWaveFunction(struct Grid *G) {
  int Ngrid = 1000;
  DOUBLE tolerance_fp = 1e-3;
  DOUBLE tolerance_Eloc = 1e-2;
  DOUBLE r, dr, epsilon;
  DOUBLE fp, Fp_check, Fp_exact;
  DOUBLE fpp, Eloc_check, Eloc_exact;
  int i;

  if(G->step == 0) {
    G->step = 1e-5;
    G->I_step = 1./G->step;
  }

  Message("\nTesting grids...\n");
#ifndef DONT_CHECK_CONTINUITY_OF_WF_AT_L_2
  Message("continuity of w.f. at Lwf/2:\n    f(Lwf/2) = 1 ... ");
  if(fabs(Exp(InterpolateU(G, Lhalfwf, 0, 0))-1.)>1e-2)
    Error("f(L/2-0) = %" LF ", f(Lwf/2) = 1\n", Exp(InterpolateU(G, Lhalfwf-1e-10, 0, 0)));
  else
    Message("done\n");

  Message("    f'(Lwf/2) = 0 ... ");
  if(fabs(InterpolateFp(G, Lhalfwf-1e-10, 0, 0))>1e-2)
    Error("f'(Lwf/2-0) = %" LF ", f'(Lwf/2) = 0\n", InterpolateFp(G, Lhalfwf-1e-10, 0, 0));
  else
    Message("done\n");
#endif

  Message("  continuity of w.f. at Rpar:\n    f(R) ... ");
#ifndef DONT_CHECK_CONTINUITY_OF_WF
  if(fabs(Exp(InterpolateU(G, Rpar-1e-10, 0, 0)) -Exp(InterpolateU(G, Rpar+1e-10, 0, 0)))>1e-2)
    Error("f(R-0) = %" LF ", f(R+0) = %" LF "\n", Exp(InterpolateU(G, Rpar-1e-10, 0, 0)),
      Exp(InterpolateU(G, Rpar+1e-10, 0, 0)));
  else
    Message("done\n    f'(R) ... ");
  if(fabs(InterpolateFp(G, Rpar-1e-10, 0, 0)-InterpolateFp(G, Rpar+1e-10, 0, 0))>1e-2) 
    Error("f'(R-0) = %" LF ", f'(R+0) = %" LF "\n", InterpolateFp(G, Rpar-1e-10, 0, 0), InterpolateFp(G, Rpar+1e-10, 0, 0));
  else
    Message("done\n");
#ifndef DONT_CHECK_CONTINUITY_OF_ELOC
  Message("    Eloc(R) ... ");
  if(fabs(InterpolateE(G, Rpar-1e-10, 0, 0)-InterpolateE(G, Rpar+1e-10, 0, 0))>1e-2*InterpolateE(G, Rpar, 0, 0))
    Error("Eloc(R-0) = %" LF ", Eloc(R+0) = %" LF "\n", InterpolateE(G, Rpar-1e-10, 0, 0), InterpolateE(G, Rpar+1e-10, 0, 0));
  else
     Message("done\n");
#else
  Message("\n    Continuity of Eloc(Rpar) is not checked.\n");
#endif
#else
  Message("\n    Continuity of Jastrow term at Rpar is not checked.\n");
#endif

  dr = (Lhalfwf-G->min)/(DOUBLE)Ngrid;
  //if(dr>1e-3) dr = 1e-3;
  // dr = 1e-3
  Message("    Testing first derivative (second order difference method), %" LG " <r<%" LG " ...\n", G->min+dr, G->min+dr*(Ngrid+1));
  //epsilon = G->step;
  //if(epsilon>1e-3) epsilon = 1e-3;
  epsilon = 1e-4;
  if(Lhalfwf>1e5) epsilon = 1e-1;

  for(i=0; i<Ngrid; i++) {
#ifdef HARD_SPHERE
    r = a+(Lhalfwf-a)/(DOUBLE)Ngrid*(DOUBLE)(i+1);
#else
    r = G->min+dr*(i+1);
#endif
    // f' = (f(i+1)-f(i-1)) / (2 dx)
    Fp_check = (InterpolateU(G, r+epsilon, 0, 0) - InterpolateU(G, r-epsilon, 0, 0)) / (2.*epsilon);
    Fp_exact = InterpolateFp(G, r, 0, 0);
    if(fabs(Fp_check/Fp_exact-1.)>tolerance_fp) {
      Warning(" r=%.3" LG " (r/L=%.3" LG "): f'(r)=%.2" LG " (exact) %.2" LG " (test), Error %2.1" LG "%%\n", r, r/Lwf, Fp_exact, Fp_check, fabs(Fp_check/Fp_exact-1.)*100);
    }
  }

  Message("\n    Testing local energy, %" LG " <r<%" LG " ...\n", G->min+dr, G->min+dr*(Ngrid+1));
  for(i=0; i<Ngrid; i++) {
#ifdef HARD_SPHERE
    r = a+(Lhalfwf-a)/(DOUBLE)Ngrid*(DOUBLE)(i+1);
#else
    r = G->min+dr*(i+1);
#endif
    fp  = (InterpolateU(G, r+epsilon, 0, 0) - InterpolateU(G, r-epsilon, 0, 0)) / (2.*epsilon);
    fpp = (InterpolateU(G, r+epsilon, 0, 0) - 2.*InterpolateU(G, r, 0, 0) + InterpolateU(G, r-epsilon, 0, 0)) / (epsilon*epsilon) + fp*fp;
#ifdef TRIAL_1D // 1D Eloc = [-f"/f] + (f'/f)^2
    Eloc_check = -fpp + fp*fp;
#endif
#ifdef TRIAL_2D // 2D Eloc = [-(f" + f'/r)/f] + (f'/f)^2
#ifdef INTERACTION_WITH_DAMPING
    Eloc_check = -(fpp + fp/r) + InteractionEnergy(r);
#else
    Eloc_check = -(fpp + fp/r) + fp*fp;
#endif
#endif
#ifdef TRIAL_3D // 3D Eloc = [-(f"+2f'/r)/f] + (f'/f)^2
    Eloc_check = -(fpp + 2*fp/r) + fp*fp;
#endif
    Eloc_exact = InterpolateE(G, r, 0, 0);

    if(fabs(Eloc_check/Eloc_exact-1.)>tolerance_Eloc && Eloc_check!= 0. && Eloc_exact!= 0) {
      Warning("    r=%.3" LE " (r/L=%.3" LF "): E(r)=%.2" LE " (exact) %.2" LE " (test), Error %2.1g%%\n", r, r/Lwf, Eloc_exact, Eloc_check, fabs(Eloc_check/Eloc_exact-1.)*100);
    }
  }

#ifdef INTERACTION_WITH_DAMPING
  for(i=0; i<G->size; i++) G->E[i] += G->fp[i]*G->fp[i] - G->ReV[i]; // add Force squared, subtract potential
#endif

  Message("\ndone\n\n");
}

/************************** Load Trial Wave Function **********************/
void LoadTrialWaveFunction(struct Grid *G, char *file_wf, int normalize_to_one) {
#ifdef INTERPOLATE_TYPE_SPLINE
  FILE *in;
  int size;
  DOUBLE dummy, dx;
  int i;

  Message("Loading trial wavefunction from file %s... \n", file_wf);

  in = fopen(file_wf, "r");

  if(in == NULL) {
    perror("\nError: can't load trial wave function from file,\nexiting ...");
    Error("\nError: can't load " INPATH "%s file,\nexiting ...", file_wf);
  }

  fscanf(in, "N= %i\n", &size);
  Message("  Spline size %i\n", size);
  AllocateWFGrid(G, size);
  G->size = size;

#ifdef JASTROW_RIGHT_BOUNDARY_BOUND_STATE // take the energy from scat. solution and use analytic wf. for no potential
  fscanf(in, "Escat= %" LF "\n", &G->k);
  G->k = sqrt(fabs(G->k));
  gaussian_alpha = gaussian_beta = 0;
  Warning("  Assuming interaction potential is added to local energy, setting gaussian_alpha = gaussian_beta = 0;\n");
#endif

  fscanf(in, "x f lnf fp Eloc E\n");
  for(i=0; i<size; i++) fscanf(in, "%" LF " %" LF " %" LF  "%" LF " %" LF " %" LF  "\n", &G->x[i], &G->f[i], &G->lnf[i], &G->fp[i], &G->E[i], &dummy);
  fclose(in);

  /*if(G->f[0] == 0.) {
    Warning("  Adjusting value of w.f. in zero to 1e-15 for spline stability\n");
    for(i=0; i<size; i++) G->f[i] += 0;
    for(i=0; i<size; i++) G->f[i] /= G->f[size-1];
  }*/

#ifdef INTERPOLATE_LOG
  Message("  Spline will fit: u(r), i.e. logarithm of the w.f. term using ln f(r) column\n");
  //for(i=0; i<size; i++) G->f[i] = log(G->f[i]);
  for(i=0; i<size; i++) {
    G->f[i] = G->lnf[i];
    if(G->lnf[i] == 1) {
      Warning(" strange value of ln f(r), check inwf.in!\n");
    }
  }

#ifndef BC_ABSENT
  if(G->f[size-1] != 0.) { // make sure that f(L/2) = 1
    Warning("  last element is different from 1. Adjusting the w.f. ...\n");
    for(i=0; i<size; i++) G->f[i] -= G->f[size-1];
  }
#endif
#else
  Message("  Spline will fit: f(r), i.e. the w.f. term\n");

#ifndef BC_ABSENT
  if(normalize_to_one && G->f[size-1] != 1.) { // make sure that f(L/2) = 1
    Warning("  last element is different from 1. Adjusting the w.f. ...\n");
    for(i=0; i<size; i++) G->f[i] /= G->f[size-1];
  }
#endif
#endif

  G->min = G->x[0];
  G->max = G->x[size-1];
  G->step = (G->max - G->min) / (DOUBLE) (size-1);
  G->I_step = 1./G->step;
  G->max2 = G->max*G->max;

  // check step in uniform grid
  dx = G->x[1]-G->x[0];
  for(i=1; i<size; i++)
    if(fabs((G->x[i]-G->x[i-1])/dx-1.)>1e-4)
      Error("  The input spline grid must be uniform. Instead element %i has different spacing, (%" LE " , %" LE " ) error %" LE " \n", i, G->x[i], G->x[i-1], fabs((G->x[i]-G->x[i-1])/dx-1));

  if(fabs(dx/G->step-1.)>1e-8)
      Error("  Grid spacing error in spline interpolation\n");

  Message("done\n");

#else // linear interpolation
  FILE *in;
  int size;
  DOUBLE n_check, dummy1, dummy2;
  int N_check;
  int i;

  Message("Loading trial wavefunction ... ");
  in = fopen("wfpar.in", "r");
  if(in == NULL) {
    perror("\nError: can't load trial wave function from file,\nexiting ...");
    Error("\nError: can't load " INPATH "%s file,\nexiting ...", "wfpar.in");
  }
  fscanf(in, "#n= %" LF "\nN= %i\nsize= %i\n", &n_check, &N_check, &size);
  Message(" done\n");
#ifdef JASTROW_RIGHT_BOUNDARY_BOUND_STATE // take the energy from scat. solution and use analytic wf. for no potential
  fscanf(in, "Escat= %" LF "\n", &G->k);
  G->k = sqrt(fabs(G->k));
  gaussian_alpha = gaussian_beta = 0;
  Warning("  Assuming interaction potential is added to local energy, setting gaussian_alpha = gaussian_beta = 0;\n");
#endif
  fclose(in);

  if(n_check != n) Error("  Cannot load the trial wave function as it corresponds to different density %" LF " (and not %" LF ")", n_check, n); 
  if(N_check != N) Error("  Cannot load the trial wave function as it corresponds to number of particles %i (and not %i)", N_check, N); 
  Message("  Allocating grid of size %i\n", size);
  AllocateWFGrid(G, size);
  G->size = size;

  in = fopen("wf.in", "r");
  if(in == NULL) {
    perror("\nError: can't load trial wave function from file,\nexiting ...");
    Error("\nError: can't load " INPATH "%s file,\nexiting ...", "wf.in");
  }
  //fscanf(in, "x  f  lnf  fp_f  Eloc\n");
#ifdef INTERACTION_WITH_DAMPING
  Message("\n    Columns:\n   1 - r\n   2 - ln f\n   3 - f'/f\n   4 - Re V(r)\n   5 - Im V(r)\n   6 - Eloc (no F^2)\n");
  for(i=0; i<size; i++) {
    fscanf(in, "%" LF " %" LF " %" LF  "%" LF " %" LF " %" LF  "\n", &G->x[i], &G->lnf[i], &G->fp[i], &G->ReV[i], &G->ImV[i], &G->E[i]);
    //G->E[i] += G->fp[i]*G->fp[i] - G->ReV[i]; // add Force squared, subtract potential
    if(G->x[i]<Lhalfwf) {
      dummy1 = G->lnf[i]; // remember ln f(L/2);
      N_check = i;
     }
  }

  if(fabs(G->x[size-1]/Lhalfwf -1)>0.01) Warning("    Last grid coordinate (%" LG " ) is different from L/2 (%" LG " )!\n", G->x[size-1], Lhalfwf);
  if(N_check < size) {
    Warning("    Adjusting grid size from %i to %i\n", size, N_check);
    size = N_check;
  }

#else // "standart" file
  Message("  Expecting 6 columns\n");
  for(i=0; i<size; i++) {
    fscanf(in, "%" LF " %" LF " %" LF  "%" LF " %" LF " %" LF  "\n", &G->x[i], &dummy1, &G->lnf[i], &G->fp[i], &dummy2, &G->E[i]);
    if(G->x[i]<1e-5) Warning("  strange x = %lf in wf.in, line %i\n", G->x[i], i+1);
  };
  if(normalize_to_one && fabs(G->x[size-1]/Lhalfwf -1)>0.01) Warning("    Last grid coordinate (%" LG " ) is different from L/2 (%" LG " )!\n", G->x[size-1], Lhalfwf);
#endif
  fclose(in);

  G->min = G->x[0];
  //G->max = Lhalfwf;
  G->max = G->x[size-1];
  if(G->max<1e-5)
    Error("  maximal distance in loaded wf.in is 0. Check number of lines.\n");
  G->step =  (G->max - G->min) / (DOUBLE) (size-1);
  G->I_step = (DOUBLE) (size-1) / (G->max - G->min);
  G->max2 = G->max*G->max;

  Message("done\n");
#endif
}


/************************** Load Interaction Potential ********************/
void LoadInteractionPotential(struct Grid *G, char *file_wf) {
  FILE *in;
  int size;
  DOUBLE dummy, dx;
  int i;

  Message("Loading interaction potential from file %s... \n", file_wf);

  in = fopen(file_wf, "r");

  if(in == NULL) {
    perror("\nError: can't load trial wave function from file,\nexiting ...");
    Error("\nError: can't load " INPATH "%s file,\nexiting ...", file_wf);
  }

  fscanf(in, "N= %i\n", &size);
  Message("  Spline size %i\n", size);
  AllocateWFGrid(G, size);
  G->size = size;

  //fscanf(in, "x f lnf fp Eloc E\n");
  for(i=0; i<size; i++) fscanf(in, "%" LF " %" LF "\n", &G->x[i], &G->f[i]);
  fclose(in);

  G->min = G->x[0];
  G->max = G->x[size-1];
  G->step = (G->max - G->min) / (DOUBLE) (size-1);
  G->I_step = 1./G->step;
  G->max2 = G->max*G->max;

  // check step in uniform grid
  dx = G->x[1]-G->x[0];
  for(i=1; i<size; i++)
    if(fabs((G->x[i]-G->x[i-1])/dx-1.)>1e-4)
      Error("  The input spline grid must be uniform. Instead element %i has different spacing, (%" LE " , %" LE " ) error %" LE " \n", i, G->x[i], G->x[i-1], fabs((G->x[i]-G->x[i-1])/dx-1));

  if(fabs(dx/G->step-1.)>1e-8)
      Error("  Grid spacing error in spline interpolation\n");

  Message("done\n");
}

/************************** Sample Trial Wave Function **********************/
void SampleTrialWaveFunction(struct Grid *G, DOUBLE min, DOUBLE max, long int size, DOUBLE (*U)(DOUBLE x), DOUBLE (*Fp)(DOUBLE x), DOUBLE (*E)(DOUBLE x)) {
  int i;
  DOUBLE x;

  Message("Sampling trial wavefunction ... ");

  AllocateWFGrid(G, size);

  G->size = size;
  G->step = (max - min) / (DOUBLE) size;
  G->I_step = (DOUBLE) size / (max - min);
  G->max = max;
  G->min = min + G->step;
  G->max2 = G->max*G->max;

  for(i=0; i<size; i++) {
    x = G->min+(DOUBLE)(i)*G->step;
    G->x[i] = x;
    G->lnf[i] = U(x);
#ifdef INTERPOLATE_LOG
    G->f[i] = G->lnf[i];
#else
    G->f[i] = exp(G->lnf[i]);
#endif
    G->fp[i] = Fp(x);
    G->E[i] = E(x);
  }
}

/************************** Sample Trial Wave Function **********************/
void CheckTrialWaveFunctionInterpolation(struct Grid *G, DOUBLE (*U)(DOUBLE x), DOUBLE (*Fp)(DOUBLE x), DOUBLE (*E)(DOUBLE x)) {
  int i;
  DOUBLE x;
  DOUBLE error_max1, error_max2, error_max3;
  DOUBLE x1, x2, x3;
  DOUBLE exact, interpolated, error;

  error_max1 = error_max2 = error_max3 = 0.;
  x1 = x2 = x3 = 0;
  for(i=0; i<G->size; i++) {
    x = (DOUBLE)(i+0.5)*G->max/(DOUBLE) G->size;

    exact = U(x);
    interpolated = InterpolateU(G, x, 0, 0);
    if(exact!=0) {
      error = fabs((interpolated-exact)/exact);
      if(error>error_max1) {
        error_max1 = error;
        x1 = x;
      }
    }

    exact = Fp(x);
    interpolated = InterpolateFp(G, x, 0, 0);
    if(exact!=0) {
      error = fabs((interpolated-exact)/exact);
      if(error>error_max2) {
        error_max2 = error;
        x2 = x;
      }
    }

    exact = E(x);
    interpolated = InterpolateE(G, x, 0, 0);
    if(exact!=0) {
      error = fabs((interpolated-exact)/exact);
      if(error>error_max3) {
        error_max3 = error;
        x3 = x;
      }
    }
  }
  Message("\n  Maximal error of w.f. interpolation:");
  Message("\n   F(r) - %" LG "%% at r=%" LG " =%" LG " L (%" LG " instead of %" LG ")", error_max1*100, x1, x1/L, InterpolateU(G, x1, 0, 0), U(x1));
  Message("\n   Fp(r)- %" LG "%% at r=%" LG " =%" LG " L (%" LG " instead of %" LG ")", error_max2*100, x2, x2/L, InterpolateFp(G, x2, 0, 0), Fp(x2));
  Message("\n   E(r) - %" LG "%% at r=%" LG " =%" LG " L (%" LG " instead of %" LG ")", error_max3*100, x3, x3/L, InterpolateE(G, x3, 0, 0), E(x3));
  Message("\ndone\n");
}

// linear interpolation
DOUBLE InterpolateGridU(struct Grid *Grid, DOUBLE r) {
  int i;
  DOUBLE dx;

  if(r<Grid->min) return Grid->lnf[0];

  i = (int) ((r-Grid->min)/Grid->step);
  if(i >= Grid->size) return 0.;

  dx = r - Grid->x[i];

  return Grid->lnf[i] + (Grid->lnf[i+1]-Grid->lnf[i]) * Grid->I_step * dx;
}

DOUBLE InterpolateGridFp(struct Grid *Grid, DOUBLE r) {
  int i;
  DOUBLE dx;

  if(r<=Grid->min) return Grid->fp[0];

  i = (int) ((r-Grid->min)/Grid->step);
  if(i >= Grid->size) return Grid->fp[Grid->size-1];

  dx = r - Grid->x[i];

  return Grid->fp[i] + (Grid->fp[i+1] - Grid->fp[i]) * Grid->I_step * dx;
}

DOUBLE InterpolateGridE(struct Grid *Grid, DOUBLE r) {
  int i;
  DOUBLE dx;

  if(r<=Grid->min) return Grid->E[0];

  i = (int) ((r-Grid->min)/Grid->step);
  if(i>= Grid->size) return 0;

  dx = r - Grid->x[i];
  return Grid->E[i]+(Grid->E[i+1] - Grid->E[i])*Grid->I_step * dx;
}

/************************** Construct Grid Lim ****************************/
void ConstructGridLim(struct Grid *G, const DOUBLE xmin, 
  const DOUBLE xmax, const long int size) {
  DOUBLE dx;

  AllocateWFGrid(G, size);

  G->min = xmin;
  G->max = xmax;
  dx = (G->max - G->min) / (DOUBLE) size;
  G->step = dx;
  G->I_step = (DOUBLE) size / (G->max - G->min);
  G->max2 = G->max*G->max;

  // f(x) = 1. - a/x 
  // f'/f = a / (x*(x-a)); 
  // E = (a / (x(x-a)))^2      [3D]
  // E = a / (x (x-a)^2)       [2D]
  // E = a (2x-a) / (x(x-a))^2 [1D]
}

/************************** Construct Grid Lim ****************************/
void ConstructGridPower(struct Grid *G, const DOUBLE xmin, 
  const DOUBLE xmax, const long int size) {
  DOUBLE dx;

  AllocateWFGrid(G, size);

  G->min = xmin;
  G->max = xmax;
  dx = (G->max - G->min) / (DOUBLE) size;
  G->step = dx;
  G->I_step = (DOUBLE) size / (G->max - G->min);
  G->max2 = G->max*G->max;

  // f(x) = 1. - a/x^n 
  // f'/f = na / (x*(x^n-a)); 
  // E = 1/x^2  na/(x^n-a) [n+2-D + na/(x^n-a)], D=1,2,3
}

/************************** Construct Grid HS *****************************/
DOUBLE A1HS, A2HS, A3HS, kHS, EHS;
// k in units of 2pi/L
void ConstructGridHS(struct Grid *G) {
  int size = 3;
  DOUBLE dx;
  DOUBLE x;

  G->min = a;
  G->max = Lhalf;
  dx = (G->max - G->min) / (DOUBLE) (size-1);
  G->step = dx;
  G->I_step = (DOUBLE) size / (G->max - G->min);
  G->max2 = G->max*G->max;

  Message("  Constructing Hard sphere Jastrow term (TRIAL_HS)\n");
  Message("  (2-body solution + symmetrized phonons) \n");
  Message("  parameters: 0<Rpar<1 \n");
  kHS = PI/(Lhalf-a); // maximal value of k in units of r_0
  kHS /= 2.*PI/L;
  Message("              0<Bpar<1/(1-2a/L) = %lf\n", kHS);
  if(Bpar<=0 || Bpar>kHS)Error("  Wrong value of Bpar = %lf\n", Bpar);
  if(!(boundary == NO_BOUNDARY_CONDITIONS || boundary == ONE_BOUNDARY_CONDITION)) { // else in units of [a]
  if(Rpar<=0 || Rpar>1) Error("  Wrong value of Rpar = %lf\n", Rpar);
    Warning("  Automatically converting matching distance: Rpar = %f [L/2] => Rpar = %f[r_0]\n", Rpar, Rpar*Lhalf);
    Rpar *= Lhalf;
  }
  if(Rpar<a) Error("  too small Rpar = %lf, minimal value is %lf [L/2]\n", Rpar, 1./Lhalf);

  Warning("  Automatically converting 2-body scattering momentum: Bpar = %lf [2pi/L] => Bpar = %lf [r_0]\n", Bpar, Bpar*2.*PI/L);
  Bpar *= 2.*PI/L;

  kHS = Bpar;
  EHS = kHS*kHS;

  if(a<1e-4) {
    Warning("  Assuming a=0\n");
    Rpar = 0.;
    A1HS = A2HS = A3HS = 0.;
  }
  else {
    // f'(Rpar)
    x = Rpar;
    A3HS = (kHS / Tg(kHS*(x-a)) - 1./x)/(2.*(1./(x*x*x)-1./((L-x)*(L-x)*(L-x))));

    // f(L/2)
    x = Lhalf;
    A2HS = A3HS/(x*x) + A3HS/((L-x)*(L-x));

    // f(Rpar)
    x = Rpar;
    A1HS = - log(Sin(kHS*(x-a))/x) + A2HS - A3HS/(x*x) - A3HS/((L-x)*(L-x));
  }
}

/* executed only if TRIAL_HS is defined */
#ifdef TRIAL_HS
DOUBLE InterpolateExactU(DOUBLE x) {
  DOUBLE f;

#ifdef SECURE
  if(x<a) {
    Warning(" call of InterpolateExactU with argument %lf < %lf!\n", x, a);
    return -100;
  }
  else 
#endif
    if(x<Rpar) {
#ifdef SECURE
    if(Sin(kHS*(x-a))<0) Warning("  InterpolateU: log of negative argument!\n");
#endif
    return A1HS + log(Sin(kHS*(x-a))/x);
  }
  else {
    return A2HS - A3HS/(x*x) - A3HS/((L-x)*(L-x));
  }
}

DOUBLE InterpolateExactFp(DOUBLE x) {

#ifdef SECURE
  if(x<a) {
    Warning("InterpolateExactFp : (r<a) call\n");
    return 100;
  }
  else
#endif
  if(x<Rpar) {
    return kHS / Tg(kHS*(x-a)) - 1./x;
  }
  else {
    return 2.*A3HS*(1./(x*x*x)-1./((L-x)*(L-x)*(L-x)));
  }
}

DOUBLE InterpolateExactE(DOUBLE x) {
  DOUBLE g;

#ifdef SECURE
  if(x<a) {
    Warning("Interpolate E : (r<a) call\n");
    return 0;
  }
  else
#endif
  if(x<Rpar) {
    g = kHS/ Tg(kHS*(x-a))-1./x;
    return EHS + g*g;
  }
  else {
   return (2.*A3HS*(L*L*L*L - 4.*L*L*L*x + 6.*L*L*x*x - 2.*L*x*x*x + 2.*x*x*x*x))/((L-x)*(L-x)*(L-x)*(L-x)*x*x*x*x);
  }
}
#endif

// k in units of 2pi/L
void ConstructGridPseudopotentialPhonon(struct Grid *G) {
  int size = 3;
  DOUBLE dx;
  DOUBLE x;

  if(a<0)
    G->min=0;
  else
    G->min = a;
  G->max = Lhalf;
  dx = (G->max - G->min) / (DOUBLE) (size-1);
  G->step = dx;
  G->I_step = (DOUBLE) size / (G->max - G->min);
  G->max2 = G->max*G->max;

  Message("  Constructing pseudopotential Jastrow term (TRIAL_PSEUDOPOTENTIAL_PHONON)\n");
  Message("  (2-body solution + symmetrized phonons) \n");
  Message("  parameters: 0<Rpar<1 \n");
  kHS = PI/(Lhalf-a); // maximal value of k in units of r_0
  kHS /= 2.*PI/L;
  Message("              0<Bpar<1/(1-2a/L) = %lf\n", kHS);
  if(Bpar<=0 || Bpar>kHS)Error("  Wrong value of Bpar = %lf\n", Bpar);
  if(!(boundary == NO_BOUNDARY_CONDITIONS || boundary == ONE_BOUNDARY_CONDITION)) { // else in units of [a]
  if(Rpar<=0 || Rpar>1) Error("  Wrong value of Rpar = %lf\n", Rpar);
    Warning("  Automatically converting matching distance: Rpar = %f [L/2] => Rpar = %f[r_0]\n", Rpar, Rpar*Lhalf);
    Rpar *= Lhalf;
  }
  Warning("  Automatically converting 2-body scattering momentum: Bpar = %lf [2pi/L] => Bpar = %lf [r_0]\n", Bpar, Bpar*2.*PI/L);
  Bpar *= 2.*PI/L;

  kHS = Bpar;
  EHS = kHS*kHS;

  if(fabs(a)<1e-4) {
    Warning("  Assuming a=0\n");
    Rpar = 0.;
    A1HS = A2HS = A3HS = 0.;
  }
  else {
    // f'(Rpar)
    x = Rpar;
    A3HS = (kHS / Tg(kHS*(x-a)) - 1./x)/(2.*(1./(x*x*x)-1./((L-x)*(L-x)*(L-x))));

    // f(L/2)
    x = Lhalf;
    A2HS = A3HS/(x*x) + A3HS/((L-x)*(L-x));

    // f(Rpar)
    x = Rpar;
    A1HS = - log(Sin(kHS*(x-a))/x) + A2HS - A3HS/(x*x) - A3HS/((L-x)*(L-x));
  }
}

#ifdef TRIAL_PSEUDOPOTENTIAL_PHONON
DOUBLE InterpolateExactU(DOUBLE x) {
  DOUBLE f;

  if(x<Rpar) {
    return A1HS + log(Sin(kHS*(x-a))/x);
  }
  else {
    return A2HS - A3HS/(x*x) - A3HS/((L-x)*(L-x));
  }
}

DOUBLE InterpolateExactFp(DOUBLE x) {
  if(x<Rpar) {
    return kHS / Tg(kHS*(x-a)) - 1./x;
  }
  else {
    return 2.*A3HS*(1./(x*x*x)-1./((L-x)*(L-x)*(L-x)));
  }
}

DOUBLE InterpolateExactE(DOUBLE x) {
  DOUBLE g;

  if(x<Rpar) {
    g = kHS/ Tg(kHS*(x-a))-1./x;
    return EHS + g*g;
  }
  else {
   return (2.*A3HS*(L*L*L*L - 4.*L*L*L*x + 6.*L*L*x*x - 2.*L*x*x*x + 2.*x*x*x*x))/((L-x)*(L-x)*(L-x)*(L-x)*x*x*x*x);
  }
}
#endif

/************************** Construct Grid HS Simple **********************/
void ConstructGridHSSimple(struct Grid *G) {
  DOUBLE x, xmin, xmax;
  DOUBLE y, ymin, ymax;
  DOUBLE V;
  DOUBLE precision = 1e-12;
  int search = ON;

  AllocateWFGrid(G, 3);
  // tg k(r-a) / rk = 1
  V = 1;

  Message(" Jastrow HS simple: f(r) = sin(k*(r-a))/r\n");

  xmin = 1e-5;
  xmax = 3.14/2./(Lhalf-a);

  // f'(L/2) = sqE / tg(sqE*(x-a)) - 1./x = 0
  // tg(sqE*(x-a))/(sqE*x) = 1

  ymin = tg(xmin*(Lhalf-a)) / (xmin*Lhalf) - V;
  ymax = tg(xmax*(Lhalf-a)) / (xmax*Lhalf) - V;

 if(ymin*ymax > 0) Error("Can't construct the grid : No solution found");

  while(search) {
    x = (xmin+xmax)/2.;
    y = tg(x*(Lhalf-a)) / (x*Lhalf) - V;

    if(y*ymin < 0) {
      xmax = x;
      ymax = y;
    } else {
      xmin = x;
      ymin = y;
    }

    if(fabs(y)<precision) {
      x = (xmax*ymin-xmin*ymax) / (ymin-ymax);
      y = tg(x*(Lhalf-a)) / (x*Lhalf) - V;
      search = OFF;
    }
  }
  sqE = x;

  Message("Solution is k = %f, error %e\n", sqE, y);

  sqE2 = sqE*sqE;
  A1HS = Lhalf/sin(sqE*(Lhalf-a));
}

#ifdef TRIAL_HS_SIMPLE
DOUBLE InterpolateExactU(DOUBLE x) {

#ifdef SECURE
  if(x>Lhalf) 
    Warning("InterpolateExactU : (r>L/2) call\n");
#endif

  if(x<a) {
#ifdef SECURE
    Warning("InterpolateExactU : (r<a) call\n");
#endif
    return -100;
  }
  else {
    return log(A1HS/x * sin(sqE*(x-a)));
  }
}

DOUBLE InterpolateExactFp(DOUBLE x) {
#ifdef SECURE
  if(x>Lhalf) Warning("InterpolateExactFp : (r>L/2) call\n");
#endif

  if(x<a) {
#ifdef SECURE
    Warning("InterpolateExactFp : (r<a) call\n");
#endif
    return 0.;
  }
  else
    return sqE / tg(sqE*(x-a)) - 1./x;
}

DOUBLE InterpolateExactE(DOUBLE x) {
  DOUBLE g;
#ifdef SECURE
  if(x>Lhalf) Warning("InterpolateExactE12 : (r>L/2) call\n");
#endif

  if(x<a) {
#ifdef SECURE
    Warning("InterpolateExactE12 : (r<a) call\n");
#endif
    return 0.;
  }
  else {
#ifdef TRIAL_3D
    g = sqE/ tg(sqE*(x-a))-1./x;
    return sqE2 + g*g;
#endif
#ifdef TRIAL_2D
    g = sqE/sin(sqE*(x-a));
    g *= g;
    g -= sqE/tg(sqE*(x-a))/x;
    return g;
#endif
#ifdef TRIAL_1D
#  error "wavefunction cannot be used in 1D"
#endif
  }
}
#endif

/************************** Construct Grid HS SW Simple **********************/
void ConstructGridHSSWSimple(struct Grid *G) {
  DOUBLE x, xmin, xmax;
  DOUBLE y, ymin, ymax;
  DOUBLE V;
  DOUBLE precision = 1e-12;
  int search = ON;

  AllocateWFGrid(G, 3);
  //Message(" Jastrow HS simple: f(r) = sin(k*(r-Apar))/r\n");

  /*//xmin = 1e-5;
  //xmax = 3.14/2./(Lhalf-Apar);
  xmin = 1.49*3.1415/(Lhalf-Apar);
  xmax = 0.51*3.1415/(Lhalf-Apar);

  // f'(L/2) = sqE / tg(sqE*(x-a)) - 1./x = 0
  // tg(sqE*(x-a))/(sqE*x) = 1
  // tg(sqE*(L/2-Apar))/(sqE*Lhalf) = 1
  ymin = tg(xmin*(Lhalf-Apar)) / (xmin*Lhalf) - 1.;
  ymax = tg(xmax*(Lhalf-Apar)) / (xmax*Lhalf) - 1.;

  if(ymin*ymax > 0) Error("Can't construct the grid : No solution found");

  while(search) {
    x = (xmin+xmax)/2.;
    y = tg(x*(Lhalf-Apar)) / (x*Lhalf) - 1.;

    if(y*ymin < 0) {
      xmax = x;
      ymax = y;
    } else {
      xmin = x;
      ymin = y;
    }

    if(fabs(y)<precision) {
      x = (xmax*ymin-xmin*ymax) / (ymin-ymax);
      y = tg(x*(Lhalf-Apar)) / (x*Lhalf) - 1.;
      search = OFF;
    }
  }
  sqE = x;*/

  sqE = sqrt(Vo);

  Atrial = Ro / sin(sqE*(Ro-a));

  //Message("Solution is k = %f, error %e\n", sqE, y);

  //sqE2 = sqE*sqE;
  //sqEkappa2 = sqE2 + Vo;
  //sqEkappa = sqrt(sqEkappa2);
  // f(L/2)
  //A1HS = Lhalf/sin(sqE*(Lhalf-Apar));
  // f(Ro)
  //Atrial  = A1HS * sin(sqE*(Ro-Apar)) / sin(sqEkappa*(Ro-a));
}

#ifdef TRIAL_HS_SW_SIMPLE
DOUBLE InterpolateExactU(DOUBLE r) {
  if(r<a) {
    return -100;
  }
  else if(r<Ro) {
    //return log(Atrial/r * sin(sqEkappa*(r-a)));
    return log(Atrial/r * sin(sqE*(r-a)));
  }
  else {
    return 0.;
    //return log(A1HS/r * sin(sqE*(r-Apar)));
  }
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  if(r<a) {
    return 0.;
  }
  else if(r<Ro) {
    //return sqEkappa / tg(sqEkappa*(r-a)) - 1./r;
    return sqE / tg(sqE*(r-a)) - 1./r;
  }
  else
    return 0.;// sqE / tg(sqE*(r-Apar)) - 1./r;
}

DOUBLE InterpolateExactE(DOUBLE r) {
  DOUBLE Fp;

  if(r<a) {
    return 0.;
  }
  else if(r<Ro) {
    //Fp = sqE/ tg(sqEkappa*(r-a))-1./r;
    //return sqEkappa2 + Fp*Fp;
    Fp = sqE/ tg(sqE*(r-a))-1./r;
    return sqE2 + Fp*Fp;
  }
  else {
    //Fp = sqE/ tg(sqE*(r-Apar))-1./r;
    //return sqE2 + Fp*Fp;
    return 0.;
  }
}
#endif


/************************** Construct Grid HS Simple **********************/
void ConstructGridPseudopotentialSimple(struct Grid *G) {
  DOUBLE x, xmin, xmax;
  DOUBLE y, ymin, ymax;
  DOUBLE V;
  DOUBLE precision = 1e-12;
  int search = ON;

  AllocateWFGrid(G, 3);
  // tg k(r-a) / rk = 1
  V = 1;

  Message(" Jastrow pseudopotential simple, a>0: f(r) = sin(k*(r-a))/r\n");
  if(a<0) Error(" w.f. cannot be used for negatve a");

  xmin = 1e-5;
  xmax = 3.14/2./(Lhalf-a);

  // f'(L/2) = sqE / tg(sqE*(x-a)) - 1./x = 0
  // tg(sqE*(x-a))/(sqE*x) = 1

  ymin = tg(xmin*(Lhalf-a)) / (xmin*Lhalf) - V;
  ymax = tg(xmax*(Lhalf-a)) / (xmax*Lhalf) - V;

 if(ymin*ymax > 0) Error("Cannot construct the grid : No solution found for positive energy\n");

  while(search) {
    x = (xmin+xmax)/2.;
    y = tg(x*(Lhalf-a)) / (x*Lhalf) - V;

    if(y*ymin < 0) {
      xmax = x;
      ymax = y;
    } else {
      xmin = x;
      ymin = y;
    }

    if(fabs(y)<precision) {
      x = (xmax*ymin-xmin*ymax) / (ymin-ymax);
      y = tg(x*(Lhalf-a)) / (x*Lhalf) - V;
      search = OFF;
    }
  }
  sqE = x;

  Message("Solution is k = %f, error %e\n", sqE, y);

  sqE2 = sqE*sqE;
  A1SW = Lhalf/sin(sqE*(Lhalf-a));
}
#ifdef TRIAL_PSEUDOPOTENTIAL_SIMPLE
DOUBLE InterpolateExactU(DOUBLE x) {
  return log(A1SW/x * fabs(sin(sqE*(x-a))));
}

DOUBLE InterpolateExactFp(DOUBLE x) {
  return sqE / tg(sqE*(x-a)) - 1./x;
}

DOUBLE InterpolateExactE(DOUBLE x) {
  DOUBLE g;

#ifdef TRIAL_3D
    g = sqE/ tg(sqE*(x-a))-1./x;
    return sqE2 + g*g;
#endif
#ifdef TRIAL_2D
    g = sqE/sin(sqE*(x-a));
    g *= g;
    g -= sqE/tg(sqE*(x-a))/x;
    return g;
#endif
}
#endif

/************************** Construct Grid SS ****************************
R is the size of a soft sphere
Rpar is the matching distance */
void ConstructGridSS(struct Grid *G, const DOUBLE a, const DOUBLE R,const DOUBLE Rpar, const DOUBLE Bpar, const int size, const DOUBLE Rmin, const DOUBLE Rmax) {
#ifdef TRIAL_3D
  DOUBLE x, xmin, xmax;
  DOUBLE y, ymin, ymax;
  DOUBLE V;
  DOUBLE precision = 1e-12;
  DOUBLE k, kappa, K, K2;
  DOUBLE delta, delta_bar;
  int search = ON;
  int repeat = ON;

  Message("Jastrow term: solution of two-body scattering problem\n");
  Message("  parameters: sphere size Apar = %" LE "\n", R);
  Message("  parameters: a = %lf scattering length\n", a);
  Message("  size of the soft sphere %" LE "\n  scattering length %" LE " \n", R, a);
  if(R<a) Error("Size of the soft sphere %" LE " is smaller than the scattering length %" LE " \n", R, a);
#ifdef BC_ABSENT // for PBC
  Message("  Wave function is constructed assuming external trap\n");
#else
  Message("  Wave function is constructed assuming PBC\n");
  Message("  parameters: Rpar = %lf = %lf [L/2] matching distance\n", Rpar, Rpar/L);
  Message("  parameters: Bpar = %lf, decay exp(-r*Bpar / Rpar)\n", Bpar);
  if(R>Rpar) Error("Size of the soft sphere %" LE " is larger than the matching distance %" LE " \n", R, Rpar);
  if(Rpar<a) Error("Matching distance %" LE " is smaller than the scattering length %" LE " \n", Rpar, a);
  if(Rpar>Lhalf) Error("Matching distance %" LE " is larger than L/2=%" LE " \n", Rpar, Lhalf);
#endif

  AllocateWFGrid(G, size);

  // define kappa
  // Solve th x / x = 1 - a/R, E<V_0
  //       tg x / x = 1 - a/R, E>V_0
  //       x = kappa R
  V = 1. - a/R;
  xmin = 1e-5;
  xmax = 100;
  ymin = tanh(xmin)/xmin - V;
  ymax = tanh(xmax)/xmax - V;
  if(ymin*ymax > 0) Error("R = %" LF " Can't construct the grid : No solution found", R);
  while(search) {
    x = (xmin+xmax)/2.;
    y = tanh(x)/x - V;
    if(y*ymin < 0) {
      xmax = x;
      ymax = y;
    } else {
      xmin = x;
      ymin = y;
    }
    if(fabs(y)<precision) {
      x = (xmax*ymin-xmin*ymax) / (ymin-ymax);
      y = tanh(x)/x - V;
      search = OFF;
    }
  }

  kappa = x / R;
  Message("  Solution is kappa = %" LF ", error %" LE " \n", kappa, y);

#ifdef BC_ABSENT // trap
  sqEkappa = kappa;
  k = 0.;
  K2 = k*k - kappa*kappa;
  K = Sqrt(fabs(K2));
  x = Apar;
  Atrial = (1 - a/x) / (sinh(sqEkappa*x) / x);
#else // PBC
  search = ON;
  // initialize delta_bar
  delta_bar = -a/Rpar;

  // loop to determine delta
  // Solve 1 + delta_bar = 1/x * arctg(x(y-2)/(x^2+y-2))
  Message("  determining delta...\n");
  while(repeat) {
    V = 1. + delta_bar;
    xmin = 1e-5;
    xmax = 2;
    ymin = Arctg(xmin*(Bpar-2.)/(xmin*xmin+Bpar-2.))/xmin - V;
    ymax = Arctg(xmax*(Bpar-2.)/(xmax*xmax+Bpar-2.))/xmax - V;
    if(ymin*ymax > 0) Error("Can't construct the grid : No solution found (2)");
    while(search) {
      x = (xmin+xmax)/2.;
      y = Arctg(x*(Bpar-2.)/(x*x+Bpar-2.))/x - V;
      if(y*ymin < 0) {
        xmax = x;
        ymax = y;
      } else {
        xmin = x;
        ymin = y;
      }
      if(fabs(y)<precision) {
        x = (xmax*ymin-xmin*ymax) / (ymin-ymax);
        y = Arctg(x*(Bpar-2.)/(x*x+Bpar-2.))/x - V;
        search = OFF;
      }
    }
    // determine momenta
    k = x / Rpar;
    K2 = k*k - kappa*kappa;
    K = Sqrt(fabs(K2));

    if(K2>0)
      high_energy = ON;
    else
      high_energy = OFF;

    // extract value of delta from momenta (scattering problem)
    if(high_energy) {
      // delta = arctg(k/K tg(KR))-kR
      delta = Arctg(k/K*Tg(K*R))-k*R;
      //Message("  low energy case, delta=%" LE " \n", delta);
    }
    else {
      // delta = arctg(k/K th(KR))-kR
      delta = Arctg(k/K*tanh(K*R))-k*R;
      //Message("  large energy case, delta=%" LE " \n", delta);
    }

    // check if delta has converged
    if(fabs((delta_bar-delta/x)/delta_bar)<precision) {
      repeat = OFF;
    }
    else {
      delta_bar = 0.5*(delta_bar + delta/x);
    }
  }
  Message("  delta=%" LF " kR=%" LF "\n", delta, k*R);

  Btrial = Rpar / Sin(x+delta) * Bpar*(Bpar-2.) / (x*x+Bpar*Bpar-2.*Bpar);
  Ctrial = Exp(Bpar) * x*x / (x*x+Bpar*Bpar-2*Bpar);
  alpha_1 = Bpar / Rpar;
  if(high_energy) {
    Atrial = Btrial * Sin(k*R+delta) / Sin(k*R);
    Message("  High energy\n");
  }
  else {
    Atrial = Btrial * Sin(k*R+delta) / sinh(K*R);
    Message("  Low energy\n");
  }

  trial_delta = delta;
#endif

  G->min = Rmin;
  G->max = Rmax;
  G->step = (G->max - G->min) / (DOUBLE) (size-1);
  G->I_step = (DOUBLE) size / (G->max - G->min);
  G->max2 = G->max*G->max;

  sqE = k;
  sqE2 = sqE*sqE;
  sqEkappa = K;
  sqEkappa2 = K2; // can be negative

  Warning("  Setting height of the soft sphere potential to Vo = %lf\n", sqE2-K2);
  Vo = sqE2-K2;
  Message("  A = %" LF " B = %" LF " C = %" LF "\n", Atrial, Btrial, Ctrial);
#else // 2D 
  DOUBLE x, xmin, xmax;
  DOUBLE y, ymin, ymax;
  DOUBLE V;
  DOUBLE precision = 1e-12;
  DOUBLE kappa;
  int search = ON;

  if(Rpar>Lhalf) Error("Size of the soft sphere %" LE "  is larger than L/2 = %" LE " \n", R, Lhalf);
  if(Rpar<a) Error("Size of the soft sphere %" LE "  is smaller than the scattering length %" LE " \n", R, a);
  if(Rpar<a) Error("Matching distance %" LE "  is smaller than the scattering length %" LE " \n", Rpar, a);

  AllocateWFGrid(G, size);

  Warning("Size of the soft sphere R = %" LE " \n\n", Rpar);

  // define kappa in units of R
  // Solve I0(x) / (x I1(x)) = ln (R/a)  
  //       
  //       x = kappa R
  V = Log(Rpar/a);
  xmin = 1e-5;
  xmax = 100;
  ymin = BesselI0(xmin)/(xmin*BesselI1(xmin)) - V;
  ymax = BesselI0(xmax)/(xmax*BesselI1(xmax)) - V;
  if(ymin*ymax > 0) Error("R = %" LF " Can't construct the grid : No solution found", Rpar);
  while(search) {
    x = (xmin+xmax)/2.;
    y = BesselI0(x)/(x*BesselI1(x)) - V;
    if(y*ymin < 0) {
      xmax = x;
      ymax = y;
    } else {
      xmin = x;
      ymin = y;
    }
    if(fabs(y)<precision) {
      x = (xmax*ymin-xmin*ymax) / (ymin-ymax);
      y = BesselI0(x)/(x*BesselI1(x)) - V;
      search = OFF;
    }
  }
  kappa = x / Rpar;
  Message("  Solution is kappa = %" LF ", error %" LE " \n", kappa, y);
  Message("  Check the value of the s-wave scattering length a/R = %" LE "\n", Exp(-BesselI0(kappa*Rpar)/BesselI1(kappa*Rpar)/(kappa*Rpar)));

  // Atrial Io(sqEkappa r), r<R
  // Btrial  + Ctrial * log(r), r>R

  Atrial = 1.;
  Ctrial = Atrial*x*BesselI1(x);
  Btrial = Atrial*BesselI0(x) - Ctrial*Log(Rpar);

  Atrial /= Btrial + Ctrial*Log(Lhalf);
  Btrial *= Atrial;
  Ctrial *= Atrial;

  G->min = Rmin;
  G->max = Rmax;
  G->step = (G->max - G->min) / (DOUBLE) (size-1);
  G->I_step = (DOUBLE) size / (G->max - G->min);
  G->max2 = G->max*G->max;

  sqEkappa = kappa;
  sqEkappa2 = kappa*kappa;

  Warning("  A = %" LF " B = %" LF " C = %" LF "\n", Atrial, Btrial, Ctrial);
#endif
}

/* executed only if TRIAL_SS is defined */
#ifdef TRIAL_SS
#ifdef TRIAL_3D

DOUBLE InterpolateExactU(DOUBLE x) {
#ifdef BC_ABSENT // 3D trap
  if(x<Apar)
    return Log(Atrial*sinh(sqEkappa*x) / x);
  else
    return Log(1 - a/x);
#else // 3D PBC
  DOUBLE f;

  if(x<Apar) {
    if(high_energy)
      f = Atrial * Sin(sqEkappa*x) / x;
    else
      f = Atrial * sinh(sqEkappa*x) / x;
  }
  else if(x<Rpar) {
    f = Btrial * Sin(sqE*x+trial_delta) / x;
  }
  else {
    f = 1. - Ctrial*Exp(-alpha_1*x);
  }
#ifndef __WATCOMC__
    if(f<1e-10) return -100.;
#endif
  return Log(f);
#endif
}

DOUBLE InterpolateExactFp(DOUBLE x) {
#ifdef BC_ABSENT // 3D trap
  if(x<Apar)
    return (sqEkappa / tanh(sqEkappa*x) - 1./x);
  else
    return a/(x*(x-a));
#else // 3D PBC
  if(x<Apar) {
    if(high_energy) return (sqEkappa / Tg(sqEkappa*x) - 1./x);
    else return (sqEkappa / tanh(sqEkappa*x) - 1./x);
  }
  else if(x<Rpar) {
    return sqE  / Tg(sqE*x+trial_delta) - 1./x;
  }
  else {
    return Ctrial*alpha_1 / (Exp(x*alpha_1) - Ctrial);
  }
#endif
}

DOUBLE InterpolateExactE(DOUBLE x) {
#ifdef BC_ABSENT // 3D trap
  DOUBLE g;
  if(x<Apar) {
    g = sqEkappa / tanh(sqEkappa*x) - 1./x;
    return sqEkappa2 + g*g; // note that sqEkappa2 < 0
  }
  else {
    g = a/(x*(x-a));
    return g*g;
  }
#else // 3D PBC
  DOUBLE g, fexp;
  if(x<Apar) {
    if(high_energy) 
      g = sqEkappa / Tg(sqEkappa*x) - 1./x;
    else 
      g = sqEkappa / tanh(sqEkappa*x) - 1./x;
    return sqEkappa2 + g*g; // no SS potential summed
    //return sqE2 + g*g; // with SS potential summed
  }
  else if(x<Rpar) {
    g =  sqE  / Tg(sqE*x+trial_delta) - 1./x;
    return sqE2 + g*g;
  }
  else {
#ifdef __WATCOMC__
    if(x*alpha_1<700)
      fexp = 1. / (Exp(x*alpha_1)-Ctrial);
    else
      fexp = 1. / (Exp(700)-Ctrial);
#else
    fexp = 1. / (Exp(x*alpha_1)-Ctrial);
#endif
    return Ctrial*alpha_1*((1+Ctrial*fexp)*alpha_1-2./x)* fexp;
  }
#endif
}
#else // 2D

DOUBLE InterpolateExactU(DOUBLE x) {
  DOUBLE f;

#ifndef SYMMETRIZE_TRIAL_WF
  if(x>Lhalf) return 0.;
#endif

  if(x<Rpar) {
    return Log(Atrial*BesselI0(sqEkappa*x));
  }
  else {
    return Log(Btrial + Ctrial*Log(x));
  }
}

DOUBLE InterpolateExactFp(DOUBLE x) {

#ifndef SYMMETRIZE_TRIAL_WF
  if(x>Lhalf) return 0.;
#endif

  if(x<Rpar) {
    return sqEkappa*BesselI1(sqEkappa*x)/BesselI0(sqEkappa*x);
  }
  else {
    return Ctrial / (x*(Btrial + Ctrial*Log(x)));
  }
}

DOUBLE InterpolateExactE(DOUBLE x) {
  DOUBLE fp;
    
#ifndef SYMMETRIZE_TRIAL_WF
  if(x>Lhalf) return 0.;
#endif

  if(x<Rpar) {
    fp = sqEkappa*BesselI1(sqEkappa*x)/BesselI0(sqEkappa*x);
    return fp*fp - sqEkappa2;
  }
  else {
    fp = Ctrial / (x*(Btrial + Ctrial*Log(x)));
    return fp*fp;
  }
}
#endif
#endif //end 2D

/************************** Construct Grid SS ****************************
R is the size of a soft sphere
Rpar is the matching distance */
void ConstructGridSSlarge(struct Grid *G) {
  DOUBLE x;
  int size = 3;

  if(Rpar>Apar) Error("Size of the soft sphere Apar = %" LE " must be larger than matching distance Rpar = %" LE " \n", Apar, Rpar);
  if(Rpar>Lhalf) Error("Matching distance %" LE " is larger than L/2=%" LE " \n", Rpar, Lhalf);

  AllocateWFGrid(G, size);

  sqEkappa2 = Vo;
  sqEkappa = sqrt(Vo);
  Message("  potential height Vo = %lf\n", Vo);
  Message("  potential size   Apar = %lf\n", Apar);
  Message("  s-wave scattering length  a/R = %lf\n", 1-tanh(sqEkappa*Apar)/(sqEkappa*Apar));

  x = Rpar;
  Ctrial = (sqEkappa / tanh(sqEkappa*x) - 1./x)/(2.*(1./(x*x*x)-1./((Lwf-x)*(Lwf-x)*(Lwf-x))));

  x = Lhalfwf;
  Btrial = Ctrial/(x*x) + Ctrial/((L-x)*(L-x));

  x = Rpar;
  Atrial  = exp(Btrial - Ctrial/(x*x) - Ctrial/((Lwf-x)*(Lwf-x)))/ (sinh(sqEkappa*x) / x);

  G->min = 0;
  G->max = Lhalfwf;
  G->step = (G->max - G->min) / (DOUBLE) (size-1);
  G->I_step = (DOUBLE) size / (G->max - G->min);
  G->max2 = G->max*G->max;
}

#ifdef TRIAL_SS_LARGE
DOUBLE InterpolateExactU(DOUBLE x) {
  if(x<Rpar) {
    return log(Atrial * sinh(sqEkappa*x) / x);
  }
  else {
    return Btrial - Ctrial/(x*x) - Ctrial/((Lwf-x)*(Lwf-x));
  }
}

DOUBLE InterpolateExactFp(DOUBLE x) {
  if(x<Rpar) {
    return (sqEkappa / tanh(sqEkappa*x) - 1./x);
  }
  else {
    return 2.*Ctrial*(1./(x*x*x)-1./((Lwf-x)*(Lwf-x)*(Lwf-x)));
  }
}

DOUBLE InterpolateExactE(DOUBLE x) {
  DOUBLE g;
  if(x<Rpar) {
    g = sqEkappa / tanh(sqEkappa*x) - 1./x;
    return -sqEkappa2 + g*g;
  }
  else {
   return (2.*Ctrial*(Lwf*Lwf*Lwf*Lwf - 4.*Lwf*Lwf*Lwf*x + 6.*Lwf*Lwf*x*x - 2.*Lwf*x*x*x + 2.*x*x*x*x))/((Lwf-x)*(Lwf-x)*(Lwf-x)*(Lwf-x)*x*x*x*x);
  }
}
#endif

/************************** Construct Grid SS ****************************
R is the size of a soft sphere
Rpar is the matching distance */
void ConstructGridSSExp(struct Grid *G) {
  DOUBLE x;
  int size = 3;

  if(Rpar>Lhalf) Error("Matching distance %" LE " is larger than L/2=%" LE " \n", Rpar, Lhalf);

  AllocateWFGrid(G, size);

  sqEkappa2 = Vo;
  sqEkappa = sqrt(Vo);
  Message("  potential height Vo = %lf\n", Vo);
  Message("  potential size   Apar = %lf\n", Apar);
  Message("  s-wave scattering length  a/R = %lf\n", 1-tanh(sqEkappa*Apar)/(sqEkappa*Apar));
  Message("  Jastrow term exp(-Bpar r^2) with Bpar = %lf\n", Bpar);

  x = Rpar;
  Ctrial = Bpar*x/((1./(x*x*x)-1./((Lwf-x)*(Lwf-x)*(Lwf-x))));

  x = Lhalfwf;
  Btrial = Ctrial/(x*x) + Ctrial/((L-x)*(L-x));

  x = Rpar;
  Atrial = -Bpar*x*x + Btrial - Ctrial/(x*x) - Ctrial/((Lwf-x)*(Lwf-x));

  G->min = 0;
  G->max = Lhalfwf;
  G->step = (G->max - G->min) / (DOUBLE) (size-1);
  G->I_step = (DOUBLE) size / (G->max - G->min);
  G->max2 = G->max*G->max;
}

#ifdef TRIAL_SS_EXP
DOUBLE InterpolateExactU(DOUBLE x) {
  if(x<Rpar) {
    return Atrial + Bpar*x*x;
  }
  else {
    return Btrial - Ctrial/(x*x) - Ctrial/((Lwf-x)*(Lwf-x));
  }
}

DOUBLE InterpolateExactFp(DOUBLE x) {
  if(x<Rpar) {
    return 2.*Bpar*x;
  }
  else {
    return 2.*Ctrial*(1./(x*x*x)-1./((Lwf-x)*(Lwf-x)*(Lwf-x)));
  }
}

DOUBLE InterpolateExactE(DOUBLE x) {
  if(x<Rpar) {
    return -6.*Bpar;
  }
  else {
   return (2.*Ctrial*(Lwf*Lwf*Lwf*Lwf - 4.*Lwf*Lwf*Lwf*x + 6.*Lwf*Lwf*x*x - 2.*Lwf*x*x*x + 2.*x*x*x*x))/((Lwf-x)*(Lwf-x)*(Lwf-x)*(Lwf-x)*x*x*x*x);
  }
}
#endif

/************************** Construct Grid 1D ******************************/
void ConstructGridHS1D(struct Grid *G, const DOUBLE R, const DOUBLE Bpar, const DOUBLE a, const int size, const DOUBLE Rmin, const DOUBLE Rmax) {
  DOUBLE x;
  DOUBLE xmin, xmax;
  DOUBLE y, ymin, ymax;
  DOUBLE precision = 1e-10;
  int search = ON;
  int i;
  DOUBLE dx;
  DOUBLE alpha;

  if(a >= R) Error("The size of the hard sphere is larger than the matching distance\n");

  AllocateWFGrid(G, size);
  i = (int) G->size;

  alpha = R/Bpar;
  alpha_1 = Bpar/R;

  /* x = Sqrt(E)
   x tg x(R-a) = 1/alpha */

  xmin = 1e-200;
  xmax = (3.14159265358979-1e-200)/(2.*(R-a));
  ymin = xmin * Tg(xmin*(R-a)) - alpha_1;
  ymax = xmax * Tg(xmax*(R-a)) - alpha_1;

  if(ymin*ymax > 0) Error("Can't construct the grid : No solution found");

  while(search) {
    x = (xmin+xmax)/2.;
    y = x * Tg(x*(R-a)) - alpha_1;

    if(y*ymin < 0) {
      xmax = x;
      ymax = y;
    } else {
      xmin = x;
      ymin = y;
    }

    if(fabs(y/alpha_1)<precision) {
      x = (xmax*ymin-xmin*ymax) / (ymin-ymax);
      y = x * Tg(x*(R-a)) - alpha_1;
      search = OFF;
    }
  }

  sqE = x;
  sqE2 = sqE*sqE;
  alpha2E = alpha*alpha*sqE2;

  G->min = Rmin;
  G->max = Rmax;
  dx = (G->max - G->min) / (DOUBLE) (size-1);
  G->step = dx;
  G->I_step = (DOUBLE) size / (G->max - G->min);
  G->max2 = G->max*G->max;

  Atrial = Sin(sqE*(R-a));
  Btrial = alpha2E / (alpha2E+1) * Exp(R*alpha_1);
  if(Btrial > 1e300)  Error("Divergence in Btrial");
}

/************************** Construct Grid 1D ******************************/
void ConstructGridFree(struct Grid *G, const DOUBLE R, const DOUBLE Bpar, const DOUBLE a, const int size, const DOUBLE Rmin, const DOUBLE Rmax) {
  int i;
  DOUBLE dx;

  AllocateWFGrid(G, size);
  i = (int) G->size;

  G->min = Rmin;
  G->max = Rmax;
  dx = (G->max - G->min) / (DOUBLE) (size-1);
  G->step = dx;
  G->I_step = (DOUBLE) size / (G->max - G->min);
  G->max2 = G->max*G->max;

  sqE = Sqrt(Bpar);
  sqE2 = Bpar;
}

/* executed only if TRIAL_TONKS is defined */
#ifdef TRIAL_TONKS
DOUBLE InterpolateExactU(DOUBLE x) {
#ifdef HARD_SPHERE
  if(x<a) return -100;
#endif
  return Log(fabs(Sin(sqE*(x-a))));
}

DOUBLE InterpolateExactFp(DOUBLE x) {
  return sqE / Tg(sqE*(x-a));
}

/* Eloc1D = -[f"/f - (f'/f)^2]*/
DOUBLE InterpolateExactE(DOUBLE x) {
  DOUBLE g;

  g = sqE/Tg(sqE*(x-a));
  //g = 0.; //Real local energy
  return sqE2 + g*g;
}
#endif

/* executed only if TRIAL_TONKS is defined */
#ifdef TRIAL_TONKS_CUTOFF
DOUBLE InterpolateExactU(DOUBLE x) {
  if(x<Rpar)
    return Log(fabs(Sin(sqE*(x-a))));
  else
    return 0;
}

DOUBLE InterpolateExactFp(DOUBLE x) {
  if(x<Rpar)
    return sqE / Tg(sqE*(x-a));
  else
    return 0;
}

/* Eloc1D = -[f"/f - (f'/f)^2]*/
DOUBLE InterpolateExactE(DOUBLE x) {
  DOUBLE g;

  if(x<Rpar) {
    g = sqE/Tg(sqE*(x-a));
    return sqE2 + g*g;
  }
  else
    return 0;
}
#endif

/* executed only if TRIAL_TONKS_SIMPLE is defined */
#ifdef TRIAL_TONKS_SIMPLE
DOUBLE InterpolateExactU(DOUBLE x) {
  return Log(fabs(Sin(sqE*x)));
}

DOUBLE InterpolateExactFp(DOUBLE x) {
  return sqE / Tg(sqE*x);
}

/* Eloc1D = -[f"/f - (f'/f)^2]*/
DOUBLE InterpolateExactE(DOUBLE x) {
  DOUBLE g;

  g = sqE/Tg(sqE*x);
  return sqE2 + g*g;
}
#endif

/************************** Construct Grid Lieb ****************************/
void ConstructGridLieb(struct Grid *G, const int size) {
  DOUBLE x;
  DOUBLE xmin, xmax;
  DOUBLE y, ymin, ymax;
  DOUBLE precision = 1e-10;
  int search = ON;
  int i;
  DOUBLE dx;
  DOUBLE alpha;

  AllocateWFGrid(G, size);
  i = (int) G->size;

  Message("  Jastrow: Lieb scattering state (E>0), matching distance Rpar [L/2], excitation branch Epar = 1,2,... \n");

  if(a>0) {
    Warning("  Wavefunction will have node at distance |a|\n");
  }
  a2 = a*a;

  if(Rpar>1) Error("  adjust Rpar to (0,1) range\n");
  Warning("  changing Rpar  = %lf [L/2] ->  %lf\n", Rpar, Rpar*Lhalf);
  Rpar *= Lhalf;
  Epar *= 2.*PI;
  if(Rpar<1e-5) {
    Warning("  Setting Rpar -> L/2 =  %lf\n", Lhalf);
    Rpar = Lhalf;
  }

  // x tg xR = -1/a1D
  alpha = -1./a;
  if(a>0) { // node
    //xmin = (PI+1e-10)/(2.*Rpar);
    //xmax = (2*PI-1e-10)/(2.*Rpar);
    xmin = (Epar + PI + 1e-10)/(2.*Rpar);
    xmax = (Epar + 2.*PI - 1e-10)/(2.*Rpar);
  }
  else { // no node
    //xmin = 1e-200;
    //xmax = (PI-1e-10)/(2.*Rpar);
    xmin = (Epar + 1e-10)/(2.*Rpar);
    xmax = (Epar + PI - 1e-10)/(2.*Rpar);
  }
  // kR = arctg ka + pi/2 + pi*n
  // kR - pi/2 - pi*n = artctg ka
  // tg(kR - pi/2 - pi*n) = ka
  // -ctg(kR) = ka
  // tg(kR) = -1/ka
  // k tg kR = -1/a
  ymin = xmin * Tg(xmin*Rpar) - alpha;
  ymax = xmax * Tg(xmax*Rpar) - alpha;

  if(ymin*ymax > 0) Error("Can't construct the grid : No solution found");

  while(search) {
    x = (xmin+xmax)/2.;
    y = x * Tg(x*Rpar) - alpha;

    if(y*ymin < 0) {
      xmax = x;
      ymax = y;
    } else {
      xmin = x;
      ymin = y;
    }

    if(fabs(y/alpha)<precision) {
      x = (xmax*ymin-xmin*ymax) / (ymin-ymax);
      y = x * Tg(x*Rpar) - alpha;
      search = OFF;
    }
  }

  sqE = x;
  Message("  k = %" LG "\n", sqE);
  sqE2 = sqE*sqE;

  G->min = 0;
  G->max = 1e100;
  dx = (G->max - G->min) / (DOUBLE) (size-1);
  G->step = dx;
  G->I_step = (DOUBLE) size / (G->max - G->min);
  G->max2 = G->max*G->max;
}

/* executed only if TRIAL_LIEB is defined */
#ifdef TRIAL_LIEB
DOUBLE InterpolateExactU(DOUBLE x) {
#ifdef HARD_SPHERE
  if(x<a) return -100;
#endif
  if(x<Rpar)
    return Log(fabs(Cos(sqE*(x-Rpar))));
  else 
    return 0;
}

DOUBLE InterpolateExactFp(DOUBLE x) {
  if(x<Rpar)
    return -sqE * Tg(sqE*(x-Rpar));
  else 
    return 0;
}

DOUBLE InterpolateExactE(DOUBLE x) {
  DOUBLE g;

  if(x<Rpar) {
    g = Tg(sqE*(x-Rpar));
    return sqE2*(1+g*g);
  }
  else {
    return 0;
  }
}
#endif

/************************** Construct Grid Phonons *************************/
void ConstructGridPhonons(struct Grid *G) {
  DOUBLE k, kmin, kmax;
  DOUBLE y, ymin, ymax;
  DOUBLE precision = 1e-12;
  int search = ON;
  DOUBLE C;
  DOUBLE R;

  if(Rpar<0) Error("Negative Rpar = %lf, adjust to Rpar>0\n", Rpar);
  if(Rpar>Lhalf) Error("Rpar = %lf, adjust to 0<Rpar<L/2=%lf\n", Rpar, Lhalf);

  Message("  Jastrow term: ConstructGridPhonons trial wavefunction\n");
  if(Rpar<1) {
    Warning("  Rpar %" LG "  [Lwf/2] -> %" LG "  [R]\n", Rpar, Rpar*Lhalfwf);
    Rpar *= Lhalfwf;
  }
  //Warning("  Rpar=%" LG "  is used as it is, without rescaling\n", Rpar);
  R = Rpar;

  AllocateWFGrid(G, 3);

  C = -L/(2*PI)*Sin(2*PI*R/L);
  k = kmin = 1e-3;

  Message("Fixing kR (1)... ");
  ymin = (a*k*Cos(k*R)+Sin(k*R))*(a*k*Sin(k*R)-Cos(k*R))/(k*(a*a*k*k+1));
  while(search) {
    k = kmax = kmin + 0.1*kmin;
    y = (a*k*Cos(k*R)+Sin(k*R))*(a*k*Sin(k*R)-Cos(k*R))/(k*(a*a*k*k+1));

    if((ymin-C)*(y-C)<0) {
      search = OFF;
    }
    else {
      kmin = k;
      ymin = y;
    }
  }

  search = ON;
  Message("done\n  Fixing kR (2)... ");

  // Solve (a*k*Cos(k*R)+Sin(k*R))*(a*k*Sin(k*R)-Cos(k*R))/(k*(a*a*k*k+1)) = -L/(2*PI)*Sin(2*PI*R/L)
  // i.e. (2.60)
  while(search) {
    k = (kmin+kmax)/2.;
    y = (a*k*Cos(k*R)+Sin(k*R))*(a*k*Sin(k*R)-Cos(k*R))/(k*(a*a*k*k+1));

    if((y-C)*(ymin-C) < 0) {
      kmax = k;
      ymax = y;
    } else {
      kmin = k;
      ymin = y;
    }

    if(fabs((y-C)/C)<precision) {
      k = (kmax*(ymin-C)-kmin*(ymax-C)) / (ymin-ymax);
      y = (a*k*Cos(k*R)+Sin(k*R))*(a*k*Sin(k*R)-Cos(k*R))/(k*(a*a*k*k+1));
      search = OFF;
    }
  }
  Message("done\n");

  Btrial = Arctg(1/(k*a))/k;
  alpha_1 = 1. + Tg(PI*R/L)*(k*L/PI/Tg(k*(R-Btrial))+Tg(PI*R/L));
  Atrial = Sin(PI*R/L);
  Atrial = pow(Sin(PI*R/L), alpha_1) / Cos(k*(R-Btrial));

  if(a<0) {
    Warning("Negative scattering length ! %" LG "  -> %" LG "\n", a, -a);
    a = -a;
    Warning("Wavefunction will have node at distance |a|\n  ");
  }

  Message("  k = %" LE " , error %" LE " \n", k, fabs((y-C)/C));
  Message("    A = %" LE " \n", Atrial);
  Message("    B = %" LE " \n", Btrial);
  Message("    alpha = %" LE " \n", alpha_1);
  Message("    Equivalent Luttiner parameter K = 1/alpha = %" LE " \n", 1./alpha_1);

  sqE = k;
  sqE2 = sqE*sqE;

  G->min = 0;
  G->max = L/2;
  G->max2 = G->max*G->max;
}

/* executed only if TRIAL_PHONON is defined */
#ifdef TRIAL_PHONON
// ideal bosons K -> \infty, alpha_1 = 1/K = 0
// ideal fermions K = 1, alpha_1 = 1/K = 1
DOUBLE InterpolateExactU(DOUBLE x) {
  if(x<Rpar) 
    return Log(fabs(Atrial*Cos(sqE*(x-Btrial))));
  else 
    return Log(pow(Sin(PI*x/L), alpha_1));
}

DOUBLE InterpolateExactFp(DOUBLE x) {
#ifdef HARD_SPHERE
#ifdef SECURE
      if(x<a) Error("InterpolateExactFp : (r<a) call\n");
#endif
#endif
  if(x<Rpar)
    return -sqE * Tg(sqE*(x-Btrial));
  else
    return alpha_1*PI/(Lwf*Tg(PI*x/L));
}

/* executed only if TRIAL_PHONON is defined */
DOUBLE InterpolateExactE(DOUBLE x) {
  DOUBLE c;

  if(x<Rpar) {
    c = Tg(sqE*(x-Btrial));
    return sqE2*(1.+c*c);
    //return sqE2; // Local energy
  }
  else {
    c = 1. / Tg(PI*x/L);
    return alpha_1*(PI*PI/(L*L))*(1.+c*c);
    //return -alpha_1*(PI*PI/(L*L))*((alpha_1-1)*c*c-1); // Local energy
  }
}
#endif

/************************** Construct Grid Phonon Luttinger *************************/
void ConstructGridPhononLuttinger(struct Grid *G) {
  DOUBLE k, kmin, kmax;
  DOUBLE y, ymin, ymax;
  DOUBLE precision = 1e-12;
  int search = ON;
  DOUBLE C;
  DOUBLE R;
  int iter = 1;

  if(Rpar<0) Error("Negative Rpar = %lf, adjust to Rpar>0\n", Rpar);
  if(Rpar>Lhalf) Error("Rpar = %lf, adjust to 0<Rpar<L/2=%lf\n", Rpar, Lhalf);

  Message("Jastrow term: ConstructGridPhononLuttinger trial wavefunction\n");
  Message("  Rpar=%" LG "  is used as it is, without rescaling\n", Rpar);
  R = Rpar;

  if(a>0) {
    Warning("  assuming repulsive case, with a<0\n");
    a = -a;
  }

  AllocateWFGrid(G, 3);

  alpha_1 = 1./Kpar;

  C = -PI*alpha_1/L/tan(PI*R/L);
  k = kmin = 1e-8;

  //y = k*tan(k*R + atan(1./(k*a))) -> 1/(a-R) for k -> 0
  if(1./(a-R) > C) {
    Warning("  no solution with small k is possible\n");
    alpha_1  = 1./(a-R)/(-PI/L/tan(PI*R/L));
    alpha_1 *= 0.999;
    Kpar = 1./alpha_1;
    Warning("  automatically choosing minimal Kpar possible, Kpar = %lf\n", Kpar);
    C = -PI*alpha_1/L/tan(PI*R/L);
  }

  Message("  Fixing kR (1)... ");
  ymin = k*tan(k*R + atan(1./(k*a))) - C;
  while(search) {
    k = kmax = kmin*1.01;
    y = k*tan(k*R + atan(1./(k*a))) - C;
    if(ymin*y<0) {
      ymax = y;
      search = OFF;
    }
    else {
      kmin = k;
      ymin = y;
    }
  }
  search = ON;
  Message("done\n  Fixing kR (2)... ");

  // Solve k*tan(k*R + atan(1./(k*a))) = -PI*alpha_1/L/tan(PI*R/L);
  // i.e. f'(R)/f(R)
  while(search) {
    k = (kmin+kmax)/2.;
    y = k*tan(k*R + atan(1./(k*a))) - C;

    if(y*ymin < 0) {
      kmax = k;
      ymax = y;
    } else {
      kmin = k;
      ymin = y;
    }

    if(fabs(y)<precision) {
      k = (kmax*ymin-kmin*ymax) / (ymin-ymax);
      y = k*tan(k*R + atan(1./(k*a)));
      search = OFF;
    }

    if(iter++ == 1000) {
      Warning("  try Rpar closer to L/2\n");
    }
    /*if(iter == 10000) {
      Warning("  not able to converge, proceeding anyway\n");
      search = OFF;
    }*/
  }
  Message("done\n");

  Btrial = - Arctg(1/(k*a))/k;
  Atrial = Sin(PI*R/L);
  Atrial = pow(Sin(PI*R/L), alpha_1) / Cos(k*(R-Btrial));

  Message("    k = %" LE " , error %" LE " \n", k, fabs((y-C)/C));
  Message("    A = %" LE " \n", Atrial);
  Message("    B = %" LE " \n", Btrial);
  Message("    Equivalent Luttiner parameter K = 1/alpha = %" LE " \n", 1./alpha_1);

  sqE = k;
  sqE2 = sqE*sqE;

  G->min = 0;
  G->max = L/2;
  G->max2 = G->max*G->max;
}

/* executed only if TRIAL_PHONON is defined */
#ifdef TRIAL_PHONON_LUTTINGER
// ideal bosons K -> \infty, alpha_1 = 1/K = 0
// ideal fermions K = 1, alpha_1 = 1/K = 1
DOUBLE InterpolateExactU(DOUBLE x) {
  if(x<Rpar) 
    return Log(fabs(Atrial*Cos(sqE*(x-Btrial))));
  else 
    return Log(pow(Sin(PI*x/L), alpha_1));
}

DOUBLE InterpolateExactFp(DOUBLE x) {
#ifdef HARD_SPHERE
#ifdef SECURE
      if(x<a) Error("InterpolateExactFp : (r<a) call\n");
#endif
#endif
  if(x<Rpar)
    return -sqE * Tg(sqE*(x-Btrial));
  else
    return alpha_1*PI/(Lwf*Tg(PI*x/L));
}

/* executed only if TRIAL_PHONON is defined */
DOUBLE InterpolateExactE(DOUBLE x) {
  DOUBLE c;

  if(x<Rpar) {
    c = Tg(sqE*(x-Btrial));
    return sqE2*(1.+c*c);
    //return sqE2; // Local energy
  }
  else {
    c = 1 / Tg(PI*x/Lwf);
    return alpha_1*(PI*PI/(Lwf*Lwf))*(1+c*c);
    //return -alpha_1*(PI*PI/(L*L))*((alpha_1-1)*c*c-1); // Local energy
  }
}
#endif

/* executed only if TRIAL_PHONON is defined */
#ifdef TRIAL_PHONON_LUTTINGER_LATTICE
DOUBLE InterpolateExactU(DOUBLE x) {
  if(x<Rpar)
    return Log(fabs(Atrial*Cos(sqE*(x-Btrial))) * (1 + Apar12*(tanh(Bpar12*(x - Rpar12)) + tanh(Bpar12*(-x - Rpar12)))));
  else 
    return Log(pow(Sin(PI*x/L), alpha_1) * (1 + Apar12*(tanh(Bpar12*(x - Rpar12)) + tanh(Bpar12*(-x - Rpar12)))));
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  if(r<Rpar)
    return -(sqE*tan((-Btrial + r)*sqE)) + Apar12*Bpar12*(1./(cosh(Bpar12*(r - Rpar12)*cosh(Bpar12*(r - Rpar12)))) - 1./(cosh(Bpar12*(r + Rpar12)*cosh(Bpar12*(r + Rpar12)))))/(1 + Apar12*tanh(Bpar12*(r - Rpar12)) - Apar12*tanh(Bpar12*(r + Rpar12)));
  else
    return (alpha_1*PI/tan((PI*r)/L))/L + Apar12*Bpar12*(1./(cosh(Bpar12*(r - Rpar12)*cosh(Bpar12*(r - Rpar12)))) - 1./(cosh(Bpar12*(r + Rpar12)*cosh(Bpar12*(r + Rpar12)))))/(1 + Apar12*tanh(Bpar12*(r - Rpar12)) - Apar12*tanh(Bpar12*(r + Rpar12)));
}

/* executed only if TRIAL_PHONON is defined */
DOUBLE InterpolateExactE(DOUBLE r) {
  if(r<Rpar) {
    return (sqE*tan((-Btrial+r)*sqE)-(Apar12*Bpar12*(1./(cosh(Bpar12*(r-Rpar12))*cosh(Bpar12*(r-Rpar12)))-1./(cosh(Bpar12*(r+Rpar12))*cosh(Bpar12*(r+Rpar12)))))/(1+Apar12*tanh(Bpar12*(r-Rpar12))-Apar12*tanh(Bpar12*(r+Rpar12))))*
    (sqE*tan((-Btrial+r)*sqE)-(Apar12*Bpar12*(1./(cosh(Bpar12*(r-Rpar12))*cosh(Bpar12*(r-Rpar12)))-1./(cosh(Bpar12*(r+Rpar12))*cosh(Bpar12*(r+Rpar12)))))/(1+Apar12*tanh(Bpar12*(r-Rpar12))-Apar12*tanh(Bpar12*(r+Rpar12)))) + 
    (2*Apar12*Bpar12*1./(cosh(Bpar12*(r - Rpar12))*cosh(Bpar12*(r - Rpar12)))*(sqE*tan((-Btrial + r)*sqE) + Bpar12*tanh(Bpar12*(r - Rpar12))) + sqE2*(1 + Apar12*tanh(Bpar12*(r - Rpar12)) - Apar12*tanh(Bpar12*(r + Rpar12))) - 
    2*Apar12*Bpar12*1./(cosh(Bpar12*(r + Rpar12))*cosh(Bpar12*(r + Rpar12)))*(sqE*tan((-Btrial + r)*sqE) + Bpar12*tanh(Bpar12*(r + Rpar12))))/(1 + Apar12*tanh(Bpar12*(r - Rpar12)) - Apar12*tanh(Bpar12*(r + Rpar12)));
  }
  else {
    return ((alpha_1*PI/tan((PI*r)/L))/L + (Apar12*Bpar12*(1./(cosh(Bpar12*(r - Rpar12))*cosh(Bpar12*(r - Rpar12))) - 1./(cosh(Bpar12*(r + Rpar12))*cosh(Bpar12*(r + Rpar12)))))/(1 + Apar12*tanh(Bpar12*(r - Rpar12)) - Apar12*tanh(Bpar12*(r + Rpar12))))*
  ((alpha_1*PI/tan((PI*r)/L))/L + (Apar12*Bpar12*(1./(cosh(Bpar12*(r - Rpar12))*cosh(Bpar12*(r - Rpar12))) - 1./(cosh(Bpar12*(r + Rpar12))*cosh(Bpar12*(r + Rpar12)))))/(1 + Apar12*tanh(Bpar12*(r - Rpar12)) - Apar12*tanh(Bpar12*(r + Rpar12)))) - 
    ((alpha_1*PI*PI*(-2 + alpha_1 + alpha_1*Cos((2*PI*r)/L))/sin((PI*r)/L)/sin((PI*r)/L)/cosh(Bpar12*(r - Rpar12))/cosh(Bpar12*(r + Rpar12))*(cosh(2*Bpar12*r) + cosh(2*Bpar12*Rpar12) - 2*Apar12*sinh(2*Bpar12*Rpar12)))/4. + 
       2*Apar12*Bpar12*L*1./(cosh(Bpar12*(r - Rpar12))*cosh(Bpar12*(r - Rpar12)))*(alpha_1*PI/tan((PI*r)/L) - Bpar12*L*tanh(Bpar12*(r - Rpar12))) + 
    2*Apar12*Bpar12*L*1./(cosh(Bpar12*(r + Rpar12))*cosh(Bpar12*(r + Rpar12)))*(-(alpha_1*PI/tan((PI*r)/L)) + Bpar12*L*tanh(Bpar12*(r + Rpar12))))/(L*L*(1 + Apar12*tanh(Bpar12*(r - Rpar12)) - Apar12*tanh(Bpar12*(r + Rpar12))));
  }
}
#endif


/* executed only if TRIAL_PHONON is defined */
#ifdef TRIAL_PHONON_LUTTINGER_PIECEWISE
DOUBLE InterpolateExactU(DOUBLE x) {
  DOUBLE f,x1,x2,A1,A2;

  if(x<Rpar) {
    if(x<Rpar12) {
      f = Apar12;
    } else if (x < 1.-Rpar12) {
      A1 = Apar12;
      A2 = Bpar12;
      x1 = Rpar12;
      x2 = 1.-Rpar12;
      f = ((A1 + A2)/2. + (A1 - A2)/2.*cos(PI*(x - x1)/(x2 - x1)));
    } else if(x < 1.+Rpar12) {
      f = Bpar12;
    } else if (x < 2.-Rpar12) {
      A1 = Bpar12;
      A2 = Cpar12;
      x1 = 1. + Rpar12;
      x2 = 2. - Rpar12;
      f = ((A1 + A2)/2. + (A1 - A2)/2.*cos(PI*(x - x1)/(x2 - x1)));
    } else {
      f = Cpar12;
    }
    return Log(f * fabs(Atrial*Cos(sqE*(x-Btrial))));
  }
  else {
    return Log(pow(Sin(PI*x/L), alpha_1));
  }
}

DOUBLE InterpolateExactFp(DOUBLE x) {
  DOUBLE fp,x1,x2,A1,A2;

  if(x<Rpar) {
    if(x<Rpar12) {
      fp = 0.;
    } else if (x < 1.-Rpar12) {
      A1 = Apar12;
      A2 = Bpar12;
      x1 = Rpar12;
      x2 = 1.-Rpar12;
      fp = ((A1 - A2)*PI*sin((PI*(x1-x))/(x1 - x2)))/((x1 - x2)*(A1 + A2 + (A1 - A2)*cos((PI*(x - x1))/(x1 - x2))));
    } else if(x < 1.+Rpar12) {
      fp = 0.;
    } else if (x < 2.-Rpar12) {
      A1 = Bpar12;
      A2 = Cpar12;
      x1 = 1. + Rpar12;
      x2 = 2. - Rpar12;
      fp = ((A1 - A2)*PI*sin((PI*(x1-x))/(x1 - x2)))/((x1 - x2)*(A1 + A2 + (A1 - A2)*cos((PI*(x - x1))/(x1 - x2))));
    } else {
      fp = 0.;
    }
    return  fp - sqE * Tg(sqE*(x-Btrial));
  }
  else {
    return alpha_1*PI/(Lwf*Tg(PI*x/L));
  }
}

DOUBLE InterpolateExactE(DOUBLE x) {
  DOUBLE c;
  DOUBLE x1,x2,A1,A2;
  DOUBLE f, Fp, fpp;

  if(x<Rpar) {
    if(x<Rpar12) {
      c = Tg(sqE*(x-Btrial));
      return sqE2*(1.+c*c);
    } else if (x < 1.-Rpar12) {
      A1 = Apar12;
      A2 = Bpar12;
      x1 = Rpar12;
      x2 = 1.-Rpar12;

      f = ((A1 + A2)/2. + (A1 - A2)/2.*cos(PI*(x - x1)/(x2 - x1))) * Atrial*Cos(sqE*(x-Btrial));
      Fp = ((A1 - A2)*PI*sin((PI*(x1-x))/(x1 - x2)))/((x1 - x2)*(A1 + A2 + (A1 - A2)*cos((PI*(x - x1))/(x1 - x2)))) - sqE * Tg(sqE*(x-Btrial));
      fpp = 0.5*Atrial*(cos((-Btrial + x)*sqE)*(-(A1 + A2)*sqE2 - ((A1 - A2)*(PI*PI + sqE2*(x1 - x2)*(x1 - x2))*cos((PI*(x - x1))/(x1 - x2)))/ ( (x1 - x2)*(x1 - x2)  )) 
        + (2.*(-A1 + A2)*PI*sqE*sin((-Btrial + x)*sqE)*sin((PI*(x - x1))/(-x1 + x2)))/(x1 - x2));
       return  -fpp / f + Fp*Fp;
    } else if(x < 1.+Rpar12) {
      c = Tg(sqE*(x-Btrial));
      return sqE2*(1.+c*c);
    } else if (x < 2.-Rpar12) {
      A1 = Bpar12;
      A2 = Cpar12;
      x1 = 1. + Rpar12;
      x2 = 2. - Rpar12;

      f = ((A1 + A2)/2. + (A1 - A2)/2.*cos(PI*(x - x1)/(x2 - x1))) * Atrial*Cos(sqE*(x-Btrial));
      Fp = ((A1 - A2)*PI*sin((PI*(x1-x))/(x1 - x2)))/((x1 - x2)*(A1 + A2 + (A1 - A2)*cos((PI*(x - x1))/(x1 - x2)))) - sqE * Tg(sqE*(x-Btrial));
      fpp = 0.5*Atrial*(cos((-Btrial + x)*sqE)*(-(A1 + A2)*sqE2 - ((A1 - A2)*(PI*PI + sqE2*(x1 - x2)*(x1 - x2))*cos((PI*(x - x1))/(x1 - x2)))/ ( (x1 - x2)*(x1 - x2)  )) 
        + (2.*(-A1 + A2)*PI*sqE*sin((-Btrial + x)*sqE)*sin((PI*(x - x1))/(-x1 + x2)))/(x1 - x2));
       return  -fpp / f + Fp*Fp;
    } else {
      c = Tg(sqE*(x-Btrial));
      return sqE2*(1.+c*c);
    }
  }
  else {
    c = 1 / Tg(PI*x/Lwf);
    return alpha_1*(PI*PI/(Lwf*Lwf))*(1+c*c);
  }
}
#endif

/************************** ConstructGridLiebLuttingerPhonon *************************/
void ConstructGridLiebLuttingerPhonon(struct Grid *G, DOUBLE K) {
  DOUBLE x, xmin, xmax;
  DOUBLE y, ymin, ymax;
  DOUBLE precision = 1e-12;
  int search = ON;
  DOUBLE C;
  DOUBLE R;

  Message("  Jastrow term: ConstructGridLiebLuttingerPhonon trial wavefunction\n");
  Warning("  K = %" LG "\n", K);

  if(fabs(a)<1e-3) Error("  Cannot use a=0, set a to a finite value!\n");
  if(a>0) Error("  Expected negative a1D. Set a<0. (Remove warning for sTG gas)\n");

  AllocateWFGrid(G, 3);

  /*// Solve x - K*Lwf/(2.*a)*tan(x) = pi/2
  C = PI/2.;
  //x = xmin = 1e-3;
  x = xmin = -PI/2.001;

  Message("fixing Atrial ... \n");
  ymin = xmin - K*Lwf/(2.*a)*tan(xmin) - C;
  x = xmax = PI/2.001;
  y = x - K*Lwf/(2.*a)*tan(x) - C;
  if(ymin*y>0) Error("Cannot initialize\n");

  while(search) {
    x = (xmin+xmax)/2.;
    y = x - K*Lwf/(2.*a)*tan(x) - C;

    if(y*ymin<0) {
      xmax = x;
      ymax = y;
    } else {
      xmin = x;
      ymin = y;
    }

    if(fabs(y)<precision) {
      x = (xmax*ymin-xmin*ymax) / (ymin-ymax);
      y = x - K*Lwf/(2.*a)*tan(x) - C;
      search = OFF;
    }
  }
  Message("done\n");

  Message("    Atrial = %" LE " , error %" LE " \n", x, y);
  Atrial = x;
  Btrial = -K/a*tan(Atrial);
  Ktrial = K;*/

  if(Apar<0) Error("  Adjust Apar so that 0<Apar<1");
  if(Apar>1) Error("  Adjust Apar so that 0<Apar<1");

  Atrial = asin(pow(Apar, K));
  Btrial = -K/a*tan(Atrial);
  Ctrial = 2.*(pi-2.*Atrial-Btrial*L)/(L*L);
  Ktrial = K;

  Message("    A = %" LE " \n", Atrial);
  Message("    B = %" LE " \n", Btrial);
  Message("    C = %" LE " \n", Ctrial);
  Message("    K = %" LE " \n", Ktrial);
}

#ifdef TRIAL_LIEB_LUTTINGER_PHONON 
// sin[pi/L(r-Atrial) + Btrial r^2]^(1/Rpar), where Rpar is equivalent to the Luttinger parameter K
// ideal bosons K -> \infty, alpha_1 = 1/K = 0
// ideal fermions K = 1, alpha_1 = 1/K = 1
DOUBLE InterpolateExactU(DOUBLE x) {
  //return Log(fabs(Sin(Atrial + Btrial*x)))/Ktrial;
  return Log(fabs(Sin(Atrial + Btrial*x + Ctrial*x*x)))/Ktrial;
}

DOUBLE InterpolateExactFp(DOUBLE x) {
  //return Btrial/(Ktrial * Tg(Atrial + Btrial*x));
  return (Btrial+2.*Ctrial*x)/(Ktrial * Tg(Atrial + Btrial*x + Ctrial*x*x));
}

DOUBLE InterpolateExactE(DOUBLE x) {
// 1D Eloc = [-f"/f] + (f'/f)^2
  DOUBLE c;

  //c = Btrial / Sin(Atrial + Btrial*x);
  //return c*c/Ktrial;

  c = 1. / Sin(Atrial + Btrial*x + Ctrial*x*x);
  return c*c/Ktrial* ((Btrial+2.*Ctrial*x)*(Btrial+2.*Ctrial*x) - Ctrial*Sin(2.*(Atrial + Btrial*x + Ctrial*x*x)));
}
#endif

/************************** Construct Grid Phonon Luttinger *************************/
void ConstructGridPhononLuttingerAuto(struct Grid *G) {
  DOUBLE k, kmin, kmax;
  DOUBLE y, ymin, ymax;
  DOUBLE precision = 1e-12;
  int search = ON;
  DOUBLE R, Rmin, Rmax;
  int iter = 1;

  Message("\n  Jastrow term: ConstructGridPhononLuttingerAuto trial wavefunction\n");

  if(a>0) {
    Warning(" assuming repulsive case, with a = %lf < 0\n", -a);
    a = -a;
  }

  AllocateWFGrid(G, 3);

  Message("  Kpar= %"LE " \n", Kpar);
  alpha_1 = 1./Kpar;

  // search for Rpar
  R = Rmin = 1e-8;
  // (2.59) - k ^ 2 = alpha_1 * (PI / L)*(PI / L)* ((alpha_1 - 1.)*cotan(PI*R/L)*cotan(PI*R/L)-1.)
  k = kmin = sqrt(-alpha_1 * (PI / L)*(PI / L)* ((alpha_1 - 1.)/tan(PI*R/L)/tan(PI*R/L)-1.));
  // (2.58) -k*tan(k*R + atan(1./(k*a))) = alpha_1*PI/L*cotan(PI*R/L)
  ymin = -k*tan(k*R + atan(1. / (k*a))) - alpha_1*PI / L/tan(PI*R / L);

  R = Rmax = Lhalf;
  k = kmax = sqrt(-alpha_1 * (PI / L)*(PI / L)* ((alpha_1 - 1.)/tan(PI*R/L)/tan(PI*R/L)-1.));
  ymax = -k*tan(k*R + atan(1. / (k*a))) - alpha_1*PI / L/tan(PI*R / L);

  if(fabs(ymax) < precision) {
    Message("  no iterations are needed\n");
    R = Rmax;
    search = OFF;
  }
  else {
    if(ymin*ymax > 0) Error("  cannot construct w.f. ()");
  }

  // Solve k*tan(k*R + atan(1./(k*a))) = -PI*alpha_1/L/tan(PI*R/L);
  // i.e. f'(R)/f(R)
  while(search) {
    R = (Rmin+Rmax)/2.;
    k = sqrt(-alpha_1 * (PI / L)*(PI / L)* ((alpha_1 - 1.)/tan(PI*R/L)/tan(PI*R/L)-1.));
    y = -k*tan(k*R + atan(1. / (k*a))) - alpha_1*PI / L/tan(PI*R / L);

    if(y*ymin < 0) {
      Rmax = R;
      ymax = y;
    } else {
      Rmin = R;
      ymin = y;
    }

    if(fabs(y)<precision) {
      R = (Rmax*ymin-Rmin*ymax) / (ymin-ymax);
      search = OFF;
    }

    if(iter++ == 1000) {
      Warning("  maximal number of iterations exceeded");
    }
  }
  k = sqrt(-alpha_1 * (PI / L)*(PI / L)* ((alpha_1 - 1.)/tan(PI*R/L)/tan(PI*R/L)-1.));
  y = -k*tan(k*R + atan(1. / (k*a))) - alpha_1*PI / L/tan(PI*R / L);
  Message("done\n");

  Btrial = - Arctg(1/(k*a))/k;
  Atrial = Sin(PI*R/L);
  Atrial = pow(Sin(PI*R/L), alpha_1) / Cos(k*(R-Btrial));
  Rpar = R;
  sqE = k;
  sqE2 = sqE*sqE;

  Message("    k = %"LE " , error %"LE "\n", k, y);
  Message("    R = %"LE " = %lf L/2\n", R, R/Lhalf);
  Message("    A = %"LE " \n", Atrial);
  Message("    B = %"LE " \n", Btrial);
  Message("    Equivalent Luttiner parameter K = 1/alpha = %"LE " \n", 1./alpha_1);
  Message("    [Eloc(R-0) -Eloc(R+0)] / Eloc(R) = %"LE "\n", (InterpolateE(G,R*0.99999,0,0)-InterpolateE(G,R*1.00001,0,0))/InterpolateE(G,R,0,0));

  G->min = 0;
  G->max = L/2;
  G->max2 = G->max*G->max;
}

/* executed only if TRIAL_PHONON is defined */
#ifdef TRIAL_PHONON_LUTTINGER_AUTO
// ideal bosons K -> \infty, alpha_1 = 1/K = 0
// ideal fermions K = 1, alpha_1 = 1/K = 1
DOUBLE InterpolateExactU(DOUBLE x) {
#ifdef BC_ABSENT
  if(x>Lhalf) return 0;
#endif
  if(x<Rpar)
    return Log(fabs(Atrial*Cos(sqE*(x-Btrial))));
  else 
    return Log(pow(Sin(PI*x/L), alpha_1));
}

DOUBLE InterpolateExactFp(DOUBLE x) {
#ifdef HARD_SPHERE
#ifdef SECURE
      if(x<a) Error("InterpolateExactFp : (r<a) call\n");
#endif
#endif
#ifdef BC_ABSENT
  if(x>Lhalf) return 0;
#endif
  if(x<Rpar)
    return -sqE * Tg(sqE*(x-Btrial));
  else
    return alpha_1*PI/(Lwf*Tg(PI*x/L));
}

/* executed only if TRIAL_PHONON is defined */
DOUBLE InterpolateExactE(DOUBLE x) {
  DOUBLE c;
#ifdef BC_ABSENT
  if(x>Lhalf) return 0;
#endif
  if(x<Rpar) {
    c = Tg(sqE*(x-Btrial));
    return sqE2*(1.+c*c);
    //return sqE2; // Local energy
  }
  else {
    c = 1 / Tg(PI*x/Lwf);
    return alpha_1*(PI*PI/(Lwf*Lwf))*(1+c*c);
    //return -alpha_1*(PI*PI/(L*L))*((alpha_1-1)*c*c-1); // Local energy
  }
}
#endif

/************************** Construct Grid Soft Sphere Phonon ***********************/
void ConstructGridSoftSpherePhonon(struct Grid *G) {
  DOUBLE x;
  DOUBLE k, kmin, kmax;
  DOUBLE y, ymin, ymax;
  DOUBLE precision = 1e-12;
  int search = ON;
  DOUBLE R;
  int iter = 1;

  Message("\n  Jastrow term: ConstructGridSoftSpherePhonon trial wavefunction\n");
  AllocateWFGrid(G, 3);

  Message("  Soft sphere diameter R = %f, height Vo = %f \n", Ro, Vo);
  Message("  Kpar= %"LE " \n", Kpar);
  alpha_1 = 1./Kpar;
  Message("  Converting units Rpar %lf [L/2]-> %lf", Rpar, Rpar*Lhalf);
  Rpar *= Lhalf;
  if(Rpar<Ro) Error("  Rpar = %lf < Ro = %lf, increase Rpar", Rpar, Ro);

  // sqrt(Vo-k*k) * tanh(sqrt(Vo-k*k)*Ro) = k / tan(k*Ro + trial_delta)
  // k / tan(k*Rpar + trial_delta) = alpha_1*PI/(Lwf*Tg(PI*Rpar/L));
  //tan(k*Ro + trial_delta) = k / (sqrt(Vo-k*k) * tanh(sqrt(Vo-k*k)*Ro))
  //tan(k*Rpar + trial_delta) =  k / (alpha_1*PI)*(Lwf*Tg(PI*Rpar/L));
  //k*Ro + trial_delta = atan(k / (sqrt(Vo-k*k) * tanh(sqrt(Vo-k*k)*Ro)))
  //k*Rpar + trial_delta = atan(k / (alpha_1*PI)*(Lwf*Tg(PI*Rpar/L)));
  k = kmin = 1e-8;
  ymin = k*(Ro - Rpar) - atan(k / (sqrt(Vo-k*k) * tanh(sqrt(Vo-k*k)*Ro))) + atan(k / (alpha_1*PI)*(Lwf*Tg(PI*Rpar/L)));
  k = kmax = sqrt(Vo*0.9999);
  ymax = k*(Ro - Rpar) - atan(k / (sqrt(Vo-k*k) * tanh(sqrt(Vo-k*k)*Ro))) + atan(k / (alpha_1*PI)*(Lwf*Tg(PI*Rpar/L)));

  if(ymin*ymax>0) Error("  cannot solve, wrong initial momentum\n");
  while(search) {
    k = (kmin+kmax)/2.;
    y = k*(Ro - Rpar) - atan(k / (sqrt(Vo-k*k) * tanh(sqrt(Vo-k*k)*Ro))) + atan(k / (alpha_1*PI)*(Lwf*Tg(PI*Rpar/L)));

    if(y*ymin < 0.) {
      kmax = k;
      ymax = y;
    } else {
      kmin = k;
      ymin = y;
    }

    if(fabs(y)<precision) {
      k = (kmax*ymin-kmin*ymax) / (ymin-ymax);
      search = OFF;
    }

    if(iter++ == 1000) {
      Warning("  maximal number of iterations exceeded");
    }
  }
  Message("done\n");
  sqE = k;
  sqE2 = sqE*sqE;
  sqEkappa = sqrt(Vo-k*k);
  sqEkappa2 = sqEkappa*sqEkappa;
  trial_delta = atan(k / (sqrt(Vo-k*k) * tanh(sqrt(Vo-k*k)*Ro))) - k*Ro;
  Message("  delta = %lf\n", trial_delta);
  trial_delta = atan(k / (alpha_1*PI)*(Lwf*Tg(PI*Rpar/L))) - k*Rpar;
  Message("  delta = %lf (consistency check)\n", trial_delta);
  Message("  F'(Ro): %lf\n", sqrt(Vo-k*k) * tanh(sqrt(Vo-k*k)*Ro) - k / tan(k*Ro + trial_delta));
  Message("  F'(Rpar): %lf\n", k / tan(k*Rpar + trial_delta) - alpha_1*PI/(Lwf*Tg(PI*Rpar/L)));

  //f(Rpar)
  x = Rpar;
  Btrial = pow(Sin(PI*x/L), alpha_1) / sin(sqE*x + trial_delta);

  //f(Ro)
  x = Ro;
  Atrial = Btrial*sin(sqE*x + trial_delta) / cosh(sqEkappa*x);

  Message("    k = %"LE " \n", k);
  Message("    kappa = sqrt(Vo-k*k) = %"LE " \n", sqEkappa);
  Message("    delta = %"LE " \n", trial_delta);
  Message("    Rpar = %lf = %lf Ro = %lf L/2\n", Rpar, Rpar/Ro, Rpar/Lhalf);
  Message("    A = %"LE " \n", Atrial);
  Message("    B = %"LE " \n", Btrial);
  Message("    Equivalent Luttiner parameter K = 1/alpha = %"LE " \n", 1./alpha_1);
  Message("    [Eloc(Ro-0) - Eloc(Ro+0)] / Eloc(Ro) = %"LE "\n", (InterpolateE(G,Ro*0.99999,0,0)-InterpolateE(G,Ro*1.00001,0,0))/InterpolateE(G,Ro,0,0));

  G->min = 0;
  G->max = L/2;
  G->max2 = G->max*G->max;
}

#ifdef TRIAL_SS_PHONON_LUTINGER
DOUBLE InterpolateExactU(DOUBLE x) {
  if(x<Ro)
    return Log(Atrial*cosh(sqEkappa*x));
  else if(x<Rpar)
      return Log(Btrial*sin(sqE*x + trial_delta));
  else
    return Log(pow(Sin(PI*x/L), alpha_1));
}

DOUBLE InterpolateExactFp(DOUBLE x) {
  if(x<Ro)
    return sqEkappa * tanh(sqEkappa*x);
  else if(x<Rpar)
    return sqE / tan(sqE*x + trial_delta);
  else
    return alpha_1*PI/(Lwf*Tg(PI*x/L));
}

DOUBLE InterpolateExactE(DOUBLE x) { //1D Eloc = [-f"/f] + (f'/f)^2
  DOUBLE c;
  if(x<Ro) {
    c = tanh(sqEkappa*x);
    return sqEkappa2*(-1. + c*c);
  }
  else if(x<Rpar) {
    c = 1. / tan(sqE*x + trial_delta);
    return sqE2*(1. + c*c);
  }
  else {
    c = 1 / Tg(PI*x/Lwf);
    return alpha_1*(PI*PI/(Lwf*Lwf))*(1+c*c);
  }
}
#endif

/************************** Construct Grid Soft Sphere Phonon ***********************/
void ConstructGridSoftSpherePhononAuto(struct Grid *G) {
  DOUBLE x;
  DOUBLE k, kmin, kmax;
  DOUBLE y, ymin, ymax;
  DOUBLE precision = 1e-12;
  int search = ON;
  DOUBLE R;
  int iter = 1;

  Message("\n  Jastrow term: ConstructGridSoftSpherePhononAuto trial wavefunction\n");
  AllocateWFGrid(G, 3);

  sqEkappa2 = Vo;
  R = Ro;
  Message("  Soft sphere diameter R = %f, height Vo = %f \n", Ro, Vo);

  k = kmin = 1e-8;
  ymin = Vo + k*tanh(k*R)*(k*tanh(k*R) - 2.*PI/L/sin(2.*pi*R/L)) - k*k;
  k = kmax = 1e8;
  ymax = Vo + k*tanh(k*R)*(k*tanh(k*R) - 2.*PI/L/sin(2.*pi*R/L)) - k*k;

  if(ymin*ymax>0) Error("  cannot solve, wrong initial momentum\n");

  // Solve k*tan(k*R) = -PI*alpha_1/L/tan(PI*R/L);
  // i.e. f'(R)/f(R)
  while(search) {
    k = (kmin+kmax)/2.;
    y = Vo + k*tanh(k*R)*(k*tanh(k*R) - 2.*PI/L/sin(2.*pi*R/L)) - k*k;

    if(y*ymin < 0) {
      kmax = k;
      ymax = y;
    } else {
      kmin = k;
      ymin = y;
    }

    if(fabs(y)<precision) {
      k = (kmax*ymin-kmin*ymax) / (ymin-ymax);
      search = OFF;
    }

    if(iter++ == 1000) {
      Warning("  maximal number of iterations exceeded");
    }
  }
  Message("done\n");
  sqE = k;
  sqE2 = k*k;

  alpha_1 = k*L/PI*tan(pi*Ro/L)*tanh(k*Ro);

  //f(Ro)
  x = Ro;
  Atrial = pow(Sin(PI*x/L), alpha_1) / cosh(sqE*x);

  Message("    k = %"LE " \n", sqE);
  Message("    R = %"LE " = %lf L/2\n", Ro, Ro/Lhalf);
  Message("    A = %"LE " \n", Atrial);
  Message("    Equivalent Luttiner parameter K = 1/alpha = %"LE " \n", 1./alpha_1);
  Message("    [Eloc(R-0) -Eloc(R+0)] / Eloc(R) = %"LE "\n", (InterpolateE(G,Ro*0.99999,0,0)-InterpolateE(G,Ro*1.00001,0,0))/InterpolateE(G,Ro,0,0));

  G->min = 0;
  G->max = L/2;
  G->max2 = G->max*G->max;
}

#ifdef TRIAL_SS_PHONON_LUTINGER_AUTO
DOUBLE InterpolateExactU(DOUBLE x) {
  if(x<Ro)
    return Log(Atrial*cosh(sqE*x));
  else
    return Log(pow(Sin(PI*x/L), alpha_1));
}

DOUBLE InterpolateExactFp(DOUBLE x) {
  if(x<Ro)
    return sqE * tanh(sqE*x);
  else
    return alpha_1*PI/(Lwf*Tg(PI*x/L));
}

DOUBLE InterpolateExactE(DOUBLE x) {
  DOUBLE c;
  if(x<Ro) {
    c = tanh(sqE*x);
    return sqE2*(-1.+c*c);
  }
  else {
    c = 1 / Tg(PI*x/Lwf);
    return alpha_1*(PI*PI/(Lwf*Lwf))*(1+c*c);
  }
}
#endif

#ifdef TRIAL_GAPPED
DOUBLE InterpolateExactU(DOUBLE x) {
  //int i;
  //DOUBLE u=0.;
  //for(i=1; i<=Nimp; i++) u -= L*L*D/(2.*pi*pi*N*i*i)*cos(2.*pi*i*x/L);
  //return u;
  return 0.5*D*x*(1.-x/L);
}

DOUBLE InterpolateExactFp(DOUBLE x) {
  // Fp = u'
  return 0.5*D*(1.-2.*x/L);
}

DOUBLE InterpolateExactE(DOUBLE x) {
  // Eloc = -u''
  return D/L;
}
#endif

/************************** ConstructGridLiebLuttingerPhonon *************************/
void ConstructGridLiebExponentialDecay(struct Grid *G) {
  DOUBLE x;

  Message("  Jastrow term: ConstructGridLiebExponentialDecay trial wavefunction\n");

  if(a>0) Error("  Expected negative a1D. Set a<0. (Remove warning for sTG gas)\n");
  Warning("  Rpar %" LG "  [Lwf/2] -> %" LG "  [R]\n", Rpar, Rpar*Lhalfwf);
  Rpar *= Lhalfwf;

  AllocateWFGrid(G, 3);

  x = Rpar;
  // f(L/2)
  // initialize with Ctrial = 1
  //Btrial = exp(Apar*(L + Rpar))/(exp(2.*Apar*Rpar)*(1. + a*Apar - Apar*Rpar) + exp(Apar*L)*(1. - a*Apar + Apar*Rpar));
  //Ctrial = 1. + 2.*Btrial*exp(-Apar*L/2.);
  Ctrial = 1./(1. + 1./(-cosh((Apar*(L - 2*Rpar))/2.) + Apar*(a - Rpar)*sinh((Apar*(L - 2*Rpar))/2.)));

  // f'(Rpar)
  Btrial = Ctrial*exp(Apar*(L + Rpar))/(exp(2.*Apar*Rpar)*(1. + a*Apar - Apar*Rpar) + exp(Apar*L)*(1. - a*Apar + Apar*Rpar));

  // f(Rpar)
  Atrial = (Ctrial - Btrial*(exp(-Apar*x)+exp(-Apar*(L-x))))/(x-a);

  Message("    A = %" LE " \n", Atrial);
  Message("    B = %" LE " \n", Btrial);
}

#ifdef TRIAL_LIEB_EXPONENTIAL_DECAY
// sin[pi/L(r-Atrial) + Btrial r^2]^(1/Rpar), where Rpar is equivalent to the Luttinger parameter K
// ideal bosons K -> \infty, alpha_1 = 1/K = 0
// ideal fermions K = 1, alpha_1 = 1/K = 1
DOUBLE InterpolateExactU(DOUBLE x) {
  if(x<Rpar)
    return log(Atrial*(x-a));
  else
    return log(Ctrial - Btrial*(exp(-Apar*x)+exp(-Apar*(L-x))));
}

DOUBLE InterpolateExactFp(DOUBLE x) {
  if(x<Rpar)
    return 1./(x-a);
  else
    return Apar*Btrial*(exp(-Apar*x)-exp(-Apar*(L-x))) / (Ctrial - Btrial*(exp(-Apar*x)+exp(-Apar*(L-x))));
}

DOUBLE InterpolateExactE(DOUBLE x) {
// 1D Eloc = [-f"/f] + (f'/f)^2
  DOUBLE c;
  if(x<Rpar)
    return 1./((x-a)*(x-a));
  else {
    c = Apar*Btrial*(exp(-Apar*x)-exp(-Apar*(L-x))) / (Ctrial - Btrial*(exp(-Apar*x)+exp(-Apar*(L-x)))); // Fp
    return Apar*Apar*Btrial*(exp(-Apar*x)+exp(-Apar*(L-x))) / (Ctrial - Btrial*(exp(-Apar*x)+exp(-Apar*(L-x)))) + c*c;
  }
}
#endif

/************************** Construct Grid Schiff Verlet **************************/
void ConstructGridSchiffVerlet(struct Grid *G, const DOUBLE Apar, const DOUBLE Bpar, const DOUBLE n) {

  Message("  Schiff Verlett w.f. will be used\n");
  AllocateWFGrid(G, 3);

#ifdef TRIAL_1D // 1D
  Atrial = pow(Apar/n, Bpar);
  Ctrial = Atrial * pow(0.5*Lwf, -Bpar);
#endif

#ifdef TRIAL_2D
  // cusp condition in 2D
  // Cpar = 0;
  // Atrial = 2;
#ifdef TRIAL_SCHIFF_VERLET_TRAP
  Message("  2D: no zero derivative at r.h.s.\n");
  Atrial = pow(Apar, Bpar);
#endif
#ifdef TRIAL_SCHIFF_VERLET_PHONON
  Message("  2D + phonons (TRIAL_SCHIFF_VERLET_PHONON)\n");
  Atrial = pow(Apar*Lwf, Bpar);
  //Cpar = Cpar/L;
  Cpar = Atrial*Bpar*pow(Lhalfwf,-Bpar-1);
#endif
#ifdef TRIAL_SCHIFF_VERLET_TRAP_PHONON2D
  Message("  2D + phonons (TRIAL_SCHIFF_VERLET_TRAP_PHONON2D)\n");
  Atrial = pow(Apar, Bpar);
  Ctrial = Atrial*pow(Lhalf, -Bpar) + Cpar*(2./Lhalf);
#endif
#endif

#ifdef TRIAL_3D
  Warning("Rpar -> Rpar L/2 (%lf -> %lf)\n", Rpar, Rpar*Lhalfwf);
  Rpar *= Lhalfwf;
#ifdef TRIAL_SCHIFF_VERLET_PHONON
  Message("  phonon asymptotics in the w.f. will be used\n");
  //Atrial = 0.5*pow(Apar, Bpar);
  Ctrial = Apar*Bpar*(Lwf-Rpar)*(Lwf-Rpar)*(Lwf-Rpar)*pow(Rpar,2.-Bpar)/(2.*(Lwf-2.*Rpar)*(Lwf*Lwf-Lwf*Rpar+Rpar*Rpar));
  Btrial = 8.*Ctrial/(Lwf*Lwf);

  Atrial = Btrial + Apar*pow(Rpar,-Bpar) - Ctrial*(1./(Rpar*Rpar)+1./((Lwf-Rpar)*(Lwf-Rpar)));
#else
  Atrial = Apar*pow(Lhalf,-Bpar);
#endif
#endif

  Message("  Schiff Verlet trial wavefunction\n");

  G->min = 0;
  G->max = Lhalf;
  G->max2 = G->max*G->max;
}

/* executed only if TRIAL_SCHIFF_VERLET is defined */
#ifdef TRIAL_SCHIFF_VERLET
DOUBLE InterpolateExactU(DOUBLE r) {
  //if(fabs(z)<0.3*sigma_lj) z = 0.3*sigma_lj; // Pederiva truncates the two-body term !
#ifdef TRIAL_1D // 1D
  return -Atrial*pow(fabs(r), -Bpar) + Ctrial;
#endif

#ifdef TRIAL_2D // 2D
#ifdef TRIAL_SCHIFF_VERLET_TRAP
  return -Atrial*pow(r, -Bpar);
#endif
#ifdef TRIAL_SCHIFF_VERLET_PHONON
  static int first_time = ON;
  static DOUBLE shift;

  if(first_time) {
    shift = Atrial/pow(Lhalfwf,Bpar)-Log(0.5*(Exp(-Cpar*Lhalfwf)+Exp(-Cpar*(-Lwf+Lhalfwf))));
  }
  if(r>Lhalfwf) return 0.;
  return shift-Atrial/pow(r,Bpar)+Log(0.5*(Exp(-Cpar*r)+Exp(-Cpar*(-Lwf+r))));
#endif
#ifdef TRIAL_SCHIFF_VERLET_TRAP_PHONON2D
  return -Atrial*pow(r, -Bpar) - Cpar*(1./r+1./(L-r)) + Ctrial;
#endif
#endif

#ifdef TRIAL_3D
#ifdef TRIAL_SCHIFF_VERLET_PHONON
  if(r<Rpar)
    return Atrial - Apar*pow(r,-Bpar);
  else if(r<Lhalfwf)
    return Btrial-Ctrial*(1./(r*r)+1./((Lwf-r)*(Lwf-r)));
  else
    return 0.;
#else
   return Atrial - Apar*pow(r,-Bpar);
#endif
#endif
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  //if(fabs(z)<0.3*sigma_lj) z = 0.3*sigma_lj;
#ifdef TRIAL_1D
  return Bpar*Atrial*pow(fabs(r), -Bpar-1);
#endif

#ifdef TRIAL_2D
#ifdef TRIAL_SCHIFF_VERLET_TRAP
  return Bpar*Atrial*pow(r, -Bpar-1);
#endif
#ifdef TRIAL_SCHIFF_VERLET_PHONON
  if(r>Lhalf) return 0.;
  return Atrial*Bpar*pow(r,-Bpar-1)-Cpar;
#endif
#ifdef TRIAL_SCHIFF_VERLET_TRAP_PHONON2D
  return Atrial*Bpar*pow(r, -Bpar-1.) + Cpar*(1./(r*r)-1./((L-r)*(L-r)));
#endif
#endif

#ifdef TRIAL_3D
#ifdef TRIAL_SCHIFF_VERLET_PHONON
  if(r<Rpar)
    return Apar*Bpar*pow(r,-1.-Bpar);
  else if(r<Lhalfwf)
    return 2.*Ctrial*(1./(r*r*r)-1./((Lwf-r)*(Lwf-r)*(Lwf-r)));
  else
    return 0.;
#else
  return Apar*Bpar*pow(r,-1.-Bpar);
#endif
#endif
}

DOUBLE InterpolateExactE(DOUBLE r) {
//1D Eloc = [-f"/f] + (f'/f)^2
#ifdef TRIAL_1D
  //if(fabs(z)<0.3*sigma_lj) z = 0.3*sigma_lj; //!
  return Bpar*(Bpar+1)*Atrial*pow(fabs(r), -Bpar-2);
#endif

//2D Eloc = [-(f"/f + f'/r)/f] + (f'/f)^2
#ifdef TRIAL_2D
#ifdef TRIAL_SCHIFF_VERLET_TRAP
  return Bpar*Bpar*Atrial*pow(r, -Bpar-2);
#endif
#ifdef TRIAL_SCHIFF_VERLET_PHONON
  if(r>Lhalf) return 0;
  return (Atrial*Bpar*Bpar+Cpar*pow(r,Bpar+1))*pow(r,-Bpar-2);
#endif
#ifdef TRIAL_SCHIFF_VERLET_TRAP_PHONON2D
  return Atrial*Bpar*Bpar*pow(r, -Bpar-2.) + Cpar*L*(L*L-3.*L*r+4.*r*r)/((L-r)*(L-r)*(L-r)*r*r*r);
#endif
#endif

//3D Eloc = [-(f"/f + 2f'/r)/f] + (f'/f)^2
#ifdef TRIAL_3D
#ifdef TRIAL_SCHIFF_VERLET_PHONON
  if(r<Rpar)
    return Apar*Bpar*(Bpar-1.)*pow(r,-2.-Bpar);
  else if(r<Lhalfwf)
    return 2.*Ctrial*(Lwf*Lwf*Lwf*Lwf-4.*Lwf*Lwf*Lwf*r+6.*Lwf*Lwf*r*r-2.*Lwf*r*r*r+2.*r*r*r*r)/((Lwf-r)*(Lwf-r)*(Lwf-r)*(Lwf-r)*r*r*r*r);
  else
    return 0.;
#else
    return Apar*Bpar*(Bpar-1.)*pow(r,-2.-Bpar);
#endif
#endif
}
#endif

/* executed only if TRIAL_DIPOLES_PHONON is defined */
#ifdef TRIAL_DIPOLES_PHONON
DOUBLE InterpolateExactU(DOUBLE r) {
  if(r<Rpar)
    return -2./Sqrt(r) + Log(Atrial);
  else if(r>Lhalfwf) 
    return 0;
  return -Ctrial*(1./r+1./(Lwf-r))+Log(Btrial);
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  if(r<Rpar)
    return 1./(r*Sqrt(r));
  else if(r>Lhalfwf)
   return 0;
  return Ctrial*(1./(r*r)-1./((Lwf-r)*(Lwf-r)));
}

DOUBLE InterpolateExactE(DOUBLE r) {
  if(r<Rpar)
    return 0.5/(r*r*Sqrt(r));
  else
    if(r>Lhalfwf) r = Lhalfwf;
  return Ctrial*Lwf*(Lwf*Lwf-3.*Lwf*r+4.*r*r)/(r*r*r*(Lwf-r)*(Lwf-r)*(Lwf-r));
}
#endif

/************************** Construct Grid Exp Phonon ***********************/
void ConstructGridExpPhonon(struct Grid *G) {

  AllocateWFGrid(G, 3);

#ifndef TRIAL_2D
  Error("  ConstructGridExpPhonon is defined only in 2D!\n");
#endif

  Message("  Jastrow trial w.f.: ConstructGridExpPhonon\n");
  if(Rpar>1 || Rpar<0) Error("  Rpar out of range (set 0<Rpar<1)!");
  Warning("  Changing Rpar -> Rpar Lhalfwf (%" LG "  -> %" LG " )\n", Rpar, Rpar*Lhalfwf);
  Rpar *= Lhalfwf;
  Apar = pow(Apar*Lwf, Bpar);

  Ctrial = (Apar*Bpar*(Lwf-Rpar)*(L-Rpar)*pow(Rpar,1-Bpar))/(Lwf*(Lwf-2*Rpar));
  Btrial = 4.*Ctrial/Lwf;
  Atrial = Btrial + Apar/pow(Rpar,Bpar) + (Ctrial*Lwf)/(-(Lwf*Rpar) + Rpar*Rpar);

  G->min = 0;
  G->max = Lhalf;
  G->max2 = G->max*G->max;
}

#ifdef TRIAL_EXP_PHONON
DOUBLE InterpolateExactU(DOUBLE r) {
  if(r<Rpar)
    return Atrial-Apar*pow(r,-Bpar);
  else if(r<Lhalfwf)
    return Btrial-Ctrial*(1/r + 1/(Lwf-r));
  else
    return 0.;
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  if(r<Rpar)
    return Apar*Bpar*pow(r,-Bpar-1);
  else if(r<Lhalfwf)
    return Ctrial*(1./(r*r)-1./((Lwf-r)*(Lwf-r)));
  else
    return 0.;
}

DOUBLE InterpolateExactE(DOUBLE r) {
//2D Eloc = [-(f"/f + f'/r)/f] + (f'/f)^2
  if(r<Rpar)
    return Apar*Bpar*Bpar*pow(r,-Bpar-2);
  else if(r<Lhalfwf)
    return Ctrial*Lwf*(Lwf*Lwf-3.*Lwf*r+4.*r*r)/((Lwf-r)*(Lwf-r)*(Lwf-r)*r*r*r);
  else
    return 0.;
}
#endif

/************************** Construct Grid Dipole **************************/
void ConstructGridDipoleKo(struct Grid *G) {
  AllocateWFGrid(G, 3);

  Message("2D Dipole Ko(r) trial wavefunction\n");
  Message("  Rpar %" LG "  [Lwf/2] -> %" LG "  [R]\n", Rpar, Rpar*Lhalfwf);
  Rpar *= Lhalfwf;

  Atrial = 1./K0opt(Rpar);
  Btrial = 1./K1opt(Rpar);
  Ctrial = 1./K2opt(Rpar);
  alpha_1 = Atrial/Btrial*(Lwf-2*Rpar)*(Lwf-Rpar)/(Lwf*Lwf*Sqrt(Rpar));
  Xi = K1opt(Rpar)*(Lwf-Rpar)*(Lwf-Rpar)*Sqrt(Rpar)/(Lwf*Lwf*(Lwf-2*Rpar)*K0opt(Rpar));

  G->min = 0;
  G->max = Lhalfwf;
  G->max2 = G->max*G->max;
}

/* executed only if TRIAL_DIPOLE_Ko is defined */
#ifdef TRIAL_DIPOLE_Ko
DOUBLE InterpolateExactU(DOUBLE r) {
  DOUBLE K0;
  if(r<Rpar) {
    K0 = K0opt(r);
    return Log(Atrial*K0) - alpha_1;
  }

  if(r>Lhalfwf) return 0; //r = Lhalfwf;
  return -Xi*(Lwf-2*r)*(Lwf-2*r)/((Lwf-r)*r);
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  if(r<Rpar) return K1opt(r)/(K0opt(r)*r*Sqrt(r));
  if(r>Lhalfwf) return 0.; //r = Lhalfwf;
  return Xi*(Lwf-2*r)*Lwf*Lwf/((Lwf-r)*(Lwf-r)*r*r);
}

DOUBLE InterpolateExactE(DOUBLE r) {
  DOUBLE t;
  if(r<Rpar) {
    t = K1opt(r)/K0opt(r);
    return (t*t-1)/(r*r*r);
  }
  if(r>Lhalf) return 0.; //r = Lhalf;
  return Xi*Lwf*Lwf*(Lwf*Lwf-3.*Lwf*r+4.*r*r)/((Lwf-r)*(Lwf-r)*(Lwf-r)*r*r*r);
}
#endif

/************************** Construct Grid Dipole **************************/
void ConstructGridDipoleKo11(struct Grid *G) {
  DOUBLE L = Lwf;

  AllocateWFGrid(G, 3);

  Message("  2D Dipole Ko11(r) trial wavefunction\n");
  Message("  Rpar %" LG "  [Lwf/2] -> %" LG "  [R]\n", Rpar, Rpar*Lhalf);
  Rpar *= Lhalf;

  Atrial = 1./K0opt(Rpar);
  Ctrial = K1opt(Rpar)*(L-Rpar)*(L-Rpar)*Sqrt(Rpar)/(L*L*(L-2*Rpar)*K0opt(Rpar));
  Btrial = K1opt(Rpar)/K0opt(Rpar)*(L-2*Rpar)*(L-Rpar)/(L*L*Sqrt(Rpar));

  G->min = 0;
  G->max = Lhalf;
  G->max2 = G->max*G->max;
}

/* executed only if TRIAL_DIPOLE_Ko is defined */
#ifdef TRIAL_DIPOLE_Ko11
DOUBLE InterpolateExactU11(DOUBLE r) {
  DOUBLE K0;
  if(r<Rpar) {
    K0 = K0opt(r);
    return Log(Atrial*K0) - Btrial;
  }

  if(r>Lhalfwf) return 0; //r = Lhalfwf;
  return -Ctrial*(Lwf-2*r)*(Lwf-2*r)/((Lwf-r)*r);
}

DOUBLE InterpolateExactFp11(DOUBLE r) {
  if(r<Rpar) return K1opt(r)/(K0opt(r)*r*Sqrt(r));
  if(r>Lhalfwf) return 0.; //r = Lhalfwf;
  return Ctrial*(Lwf-2*r)*Lwf*Lwf/((Lwf-r)*(Lwf-r)*r*r);
}

DOUBLE InterpolateExactE11(DOUBLE r) {
  DOUBLE t;
  if(r<Rpar) {
    t = K1opt(r)/K0opt(r);
    return (t*t-1)/(r*r*r);
  }
  if(r>Lhalf) return 0.; //r = Lhalf;
  return Ctrial*Lwf*Lwf*(Lwf*Lwf-3.*Lwf*r+4.*r*r)/((Lwf-r)*(Lwf-r)*(Lwf-r)*r*r*r);
}
#endif

#ifdef TRIAL_DIPOLE_Ko_TRAP
DOUBLE InterpolateExactU(DOUBLE r) {
  return Log(K0opt(r/D));
}

DOUBLE InterpolateExactFp(DOUBLE r) {
// Fp = f'(r)/f(r)
  return K1opt(r/D)*Sqrt(D/r)/(K0opt(r/D)*r);
}

DOUBLE InterpolateExactE(DOUBLE r) {
// Eloc(r) = -f"(r)/f(r) - (D-1)/r f'(r)/f(r) + (f'(r)/f(r))^2
  DOUBLE t;
  t = K1opt(r/D)/K0opt(r/D);
  return D*(t*t-1)/(r*r*r);
}
#endif

#ifdef TRIAL_DIPOLE_Ko11_TRAP
DOUBLE InterpolateExactU11(DOUBLE r) {
  return Log(K0opt(r/D));
}

DOUBLE InterpolateExactFp11(DOUBLE r) {
// Fp = f'(r)/f(r)
  return K1opt(r/D)*Sqrt(D/r)/(K0opt(r/D)*r);
}

DOUBLE InterpolateExactE11(DOUBLE r) {
// Eloc(r) = -f"(r)/f(r) - (D-1)/r f'(r)/f(r) + (f'(r)/f(r))^2
  DOUBLE t;
  t = K1opt(r/D)/K0opt(r/D);
  return D*(t*t-1)/(r*r*r);
}
#endif

/************************** Construct Grid Dipole **************************/
void ConstructGridDipolePhonon(struct Grid *G) {

  AllocateWFGrid(G, 3);

  Message("  2D Dipole phonon trial wavefunction\n");
  Warning("  Rpar %" LG "  [Lwf/2] -> %" LG "  [R]\n", Rpar, Rpar*Lhalfwf);
  Rpar *= Lhalfwf;

#ifdef TRIAL_1D
  Error(" Dipole-phonon w.f. is implemented in 2D only!");
#endif
#ifdef TRIAL_3D
  Error(" Dipole-phonon w.f. is implemented in 2D only!");
#endif

  Ctrial = Sqrt(Rpar)*(Lwf-Rpar)*(Lwf-Rpar)/(Lwf-2*Rpar)/Lwf;
  Btrial = Exp(2*Ctrial/Lhalfwf);
  Atrial = Btrial*Exp(2./Sqrt(Rpar)-Ctrial*(1/(Rpar)+1./(Lwf-Rpar)));

  G->min = 0;
  G->max = Lwf/2;
  G->max2 = G->max*G->max;
}

/************************** Construct Grid Dipole **************************/
void ConstructGridDipolePhononFiniteD(struct Grid *G) {

  AllocateWFGrid(G, 3);

  Message("  2D Dipole phonon finte D trial wavefunction\n");
  Message("    f(r) = C1 Ko[2 (D/2)^1/2], r<Rpar\n");
  Message("         = C2 exp{-C3/r-C3/(L-r)}, r>Rpar\n");
//  Message("    K0opt(r/D), r<Rpar\n");
//  Message("    exp{-Apar*(1./r+1./(Lwf-r))}, r>Rpar\n");
  Warning("    Rpar %" LG "  [Lwf/2] -> %" LG "  [R]\n", Rpar, Rpar*Lhalfwf);
  Rpar *= Lhalfwf;

#ifdef TRIAL_1D
  Error(" Dipole-phonon w.f. is implemented in 2D only!");
#endif
#ifdef TRIAL_3D
  Error(" Dipole-phonon w.f. is implemented in 2D only!");
#endif

  Ctrial = K1opt(Rpar/D)*Sqrt(D/Rpar)/(K0opt(Rpar/D)*Rpar)/ (1./(Rpar*Rpar)-1./((Lwf-Rpar)*(Lwf-Rpar)));
  Btrial = 4*Ctrial/Lwf ;
  Atrial = Btrial-Ctrial*(1/(Rpar)+1./(Lwf-Rpar))-Log(K0opt(Rpar/D));

  Message("    C1 = %.15" LG "\n", exp(Atrial));
  Message("    C2 = %.15" LG "\n", exp(Btrial));
  Message("    C3 = %.15" LG "\n", Ctrial);

  G->min = 0;
  G->max = Lwf/2;
  G->max2 = G->max*G->max;
}

// executed only if TRIAL_DIPOLE_PHONON_FINITE_D is defined
#ifdef TRIAL_DIPOLE_PHONON_FINITE_D
DOUBLE InterpolateExactU(DOUBLE r) {
  if(r<Rpar)
    return Atrial+Log(K0opt(r/D));
  else
    if(r>Lhalfwf) return 0;

  return -Ctrial*(1./r+1./(Lwf-r))+Btrial;
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  if(r<Rpar)
    return K1opt(r/D)*Sqrt(D/r)/(K0opt(r/D)*r);
  else
    if(r>Lhalfwf) return 0;

  return Ctrial*(1./(r*r)-1./((Lwf-r)*(Lwf-r)));
}

DOUBLE InterpolateExactE(DOUBLE r) {
  DOUBLE t;
  if(r<Rpar) {
    t = K1opt(r/D)/K0opt(r/D);
    return D*(t*t-1)/(r*r*r);
  }
  else if(r>Lhalfwf) 
    return 0.;
  return Ctrial*Lwf*(Lwf*Lwf-3.*Lwf*r+4.*r*r)/(r*r*r*(Lwf-r)*(Lwf-r)*(Lwf-r));
}
#endif

/************************** Construct Grid HS *****************************/
void ConstructGridShortRange3D(struct Grid *G, const DOUBLE a, const DOUBLE Rpar) {
  DOUBLE x, xmin, xmax;
  DOUBLE y, ymin, ymax;
  DOUBLE V;
  DOUBLE precision = 1e-12;
  int search = ON;
  int i;
  DOUBLE dx;
  int size = 3;

  AllocateWFGrid(G, size);
  i = (int) G->size;

  V = a/Rpar;

  // x = kR
  xmin = 1e-10;
  xmax = 2.79838;

  ymin = Tg(xmin-Arctg(xmin))/xmin - V;
  ymax = Tg(xmax-Arctg(xmax))/xmax - V;

  if(ymin*ymax > 0) Error("Can't construct the grid : No solution found");

  while(search) {
    x = (xmin+xmax)/2.;
    y = Tg(x-Arctg(x))/x - V;

    if(y*ymin < 0) {
      xmax = x;
      ymax = y;
    } else {
      xmin = x;
      ymin = y;
    }

    if(fabs(y)<precision) {
      x = (xmax*ymin-xmin*ymax) / (ymin-ymax);
      y = Tg(x-Arctg(x))/x - V;
      search = OFF;
    }
  }
  Message("Solution is k = %" LF ", error %" LE " \n", x/Rpar, y);

  sqE = x/Rpar;
  sqE2 = sqE*sqE;
  Btrial = Arctg(x) - x;
  Atrial = Lwf/2./Sin(sqE*Lwf/2.+Btrial);

  if(fabs(Tg(Btrial+x)-x)>1e-10) Error(" Boundary condition violated.");
  if(fabs(Tg(Btrial) + sqE*a)>1e-10) Error(" Matching condition violated.");

  G->min = 0;
  G->max = Lwf;
  dx = (G->max - G->min) / (DOUBLE) (size-1);
  G->step = dx;
  G->I_step = (DOUBLE) size / (G->max - G->min);
  G->max2 = G->max*G->max;
}

/************************ Interpolate U ************************************/
/* U = ln(f) */

/* executed only if TRIAL_LIM is defined */
#ifdef TRIAL_LIM
DOUBLE InterpolateExactU(DOUBLE r) {

  if(r>a) {
    return Log(1. - a / r);
  }
  else {
    return -100.;
  }
}
#endif

/* executed only if TRIAL_POWER is defined */
#ifdef TRIAL_POWER
DOUBLE InterpolateExactU(DOUBLE r) {

  if(r>a) {
    return Log(1. - a *pow(r, -Apar))-Log(1. - a *pow(Lhalfz  , -Apar));
  }
  else {
    return -100.;
  }
}
#endif

/* executed only if TRIAL_HS1D is defined */
#ifdef TRIAL_HS1D
DOUBLE InterpolateExactU(DOUBLE x) {
  if(x<a) {
    return -100;
  }
  else 
  if(x<Rpar) {
    return Log(fabs(Atrial * Sin(sqE*(x-a))));
  }
  else {
#ifndef __WATCOMC__
  if (1. - Btrial*Exp(-alpha_1*x)<1e-10) {
    Warning("log");
    return 0;
  }
#endif
    return Log(1. - Btrial*Exp(-alpha_1*x));
  }
}
#endif

#ifdef TRIAL_SR3D
DOUBLE InterpolateExactU(DOUBLE r) {
  if(r<Rpar)
    return Log(fabs(Atrial/r *Sin(sqE*r+Btrial)));
  else
    return 0.;
}
#endif

/************************ InterpolateExactFp ************************************/
/* fp = f'/f */

/* executed only if TRIAL_LIM is defined */
#ifdef TRIAL_LIM
DOUBLE InterpolateExactFp(DOUBLE r) {
  if(r>a)
    return a / (r*(r-a));
  else
    return 0.;
}
#endif

/* executed only if TRIAL_POWER is defined */
#ifdef TRIAL_POWER
DOUBLE InterpolateExactFp(DOUBLE r) {
  if(r>a)
    return Apar*a / (r*(pow(r,Apar)-a));
  else
    return 0.;
}
#endif

/* executed only if TRIAL_HS1D is defined */
#ifdef TRIAL_HS1D
DOUBLE InterpolateExactFp(DOUBLE x) {

#ifdef SECURE
  if(x<a) Error("InterpolateExactFp : (r<a) call\n");
#endif

  if(x<Rpar) {
    return sqE / Tg(sqE*(x-a));
  }
  else {
#ifdef __WATCOMC__
    if(x*alpha_1<700)
      return Btrial*alpha_1 / (Exp(x*alpha_1) - Btrial);
    else
      return Btrial*alpha_1 / (Exp(700) - Btrial);
#endif
    return Btrial*alpha_1 / (Exp(x*alpha_1) - Btrial);
  }
}
#endif

#ifdef TRIAL_SR3D
DOUBLE InterpolateExactFp(DOUBLE r) {
  if(r<Rpar)
    return sqE / Tg(sqE*r + Btrial) - 1./r;
  else
    return 0.;
}
#endif

/************************ InterpolateE ************************************/
/* E =  -(f"/f + f'/f*(2/x - f'/f)); */

DOUBLE InterpolateEG(DOUBLE r) {
#ifdef TRIAL_3D
  return 4.*PI*r*r*InteractionEnergy(r);
#endif

#ifdef TRIAL_2D
  return 2.*PI*r*InteractionEnergy(r);
#endif

#ifdef TRIAL_1D
  return 2.*InteractionEnergy(r);
#endif
}

/* executed only if TRIAL_LIM is defined */
#ifdef TRIAL_LIM
DOUBLE InterpolateExactE(DOUBLE r) {
  DOUBLE e;

  if(r>a) {
#ifdef TRIAL_3D // 3D case: Eloc = (a / (x(x-a)))^2
    e = a/(r*(r-a));
    return (e*e);
#endif

#ifdef TRIAL_2D // 2D case: Eloc = a / (x (x-a)^2)
    e = r-a;
    return a/(r*e*e);
#endif

#ifdef TRIAL_1D // 1D case: Eloc = a (2x-a) / (x(x-a))^2
    e = r*(r-a);
    return a*(2*r-a)/(e*e);
#endif

  }
  else { // i.e. r <= a
    return 0.;
  }
}
#endif

/* executed only if TRIAL_POWER is defined */
#ifdef TRIAL_POWER
DOUBLE InterpolateExactE(DOUBLE r) {
// E = 1/x^2  na/(x^n-a) [n+2-D + na/(x^n-a)], D=1,2,3
  DOUBLE e;

  if(r>a) {
    e = Apar*a/(pow(r,Apar)-a);
#ifdef TRIAL_3D // 3D case:
    return e*(Apar-1.+e)/(r*r);
#endif
#ifdef TRIAL_2D // 2D case:
    return e*(Apar+e)/(r*r);
#endif
#ifdef TRIAL_1D // 1D case:
    return e*(Apar+1.+e)/(r*r);
#endif
  }
  else { // i.e. r <= a
    return 0.;
  }
}
#endif

/* executed only if TRIAL_HS1D is defined */
#ifdef TRIAL_HS1D
DOUBLE InterpolateExactE(DOUBLE x) {
/* Eloc1D = -[f"/f - (f'/f)^2]*/
  DOUBLE g, fexp;

#ifdef SECURE
  if(x<a) Error("Interpolate E : (r<a) call\n");
#endif

  if(x<Rpar) {
    g = sqE/Tg(sqE*(x-a));
    //g = 0.; //Real local energy
    return sqE2 + g*g;
  }
  else {
#ifdef __WATCOMC__
    if(x*alpha_1<700)
      fexp = 1. / (Exp(x*alpha_1)-Btrial);
    else
      fexp = 1. / (Exp(700)-Btrial);
#else
    fexp = 1. / (Exp(x*alpha_1)-Btrial);
#endif
    //return Btrial*alpha_1*alpha_1*Exp(-x*alpha_1)/((1.-Btrial*Exp(-x*alpha_1)));//Real local energy
    return Btrial*alpha_1*alpha_1*(1.+Btrial*fexp)* fexp;
  }
}
#endif

/* Eloc1D = -[f"/f - (f'/f)^2]*/
/* executed only if TRIAL_SR3D is defined */
#ifdef TRIAL_SR3D
DOUBLE InterpolateExactE(DOUBLE r) {
  DOUBLE g;

  if(r<Rpar) {
    g = sqE/Tg(sqE*r + Btrial) - 1./r;
    return sqE2 + g*g;
  }
  else

    return 0.;
}
#endif

void ConstructGridSutherland(struct Grid *G) {
  AllocateWFGrid(G, 3);

  //alpha_1 = 0.5*(1+Sqrt(1+4*D));
  Warning("  Input parameter D means lambda of the CSM\n");
  alpha_1 = D;
  D = D*(D-1.);
  Message("  Sutherland parameter lambda = %.15" LE "\n", alpha_1);

  sqE = PI/Lwf;
  sqE2 = sqE*sqE;
  Message("  Energy per particle of Calogero-Sutherland gas %.15" LE "\n", alpha_1*alpha_1*PI*PI*n*n/6);

  G->min = 0;
  G->max = Lwf/2;
  G->max2 = G->max*G->max;
}

#ifdef TRIAL_SUTHERLAND
DOUBLE InterpolateExactU(DOUBLE x) {
  if(x<Lhalfwf)
    return Log(pow(Sin(sqE*x), alpha_1));
  else
    return 0;
}

DOUBLE InterpolateExactFp(DOUBLE x) {
  if(x<Lhalfwf)
    return alpha_1*sqE/(Tg(sqE*x));
  else
    return 0;
}

DOUBLE InterpolateExactE(DOUBLE x) {
  DOUBLE c;

  c = 1. / Sin(sqE*x);
  if(x<Lhalfwf)
    return alpha_1*sqE2*c*c;
  else
    return 0;
}
#endif

/*#ifdef TRIAL_SUTHERLAND
DOUBLE InterpolateExactU(DOUBLE x) {
  return Log(pow(Sin(sqE*x), alpha_1));
}

DOUBLE InterpolateExactFp(DOUBLE x) {
  return alpha_1*sqE/(Tg(sqE*x));
}

DOUBLE InterpolateExactE(DOUBLE x) {
  DOUBLE c;

  c = 1 / Tg(sqE*x);
  return alpha_1*sqE2*(1+c*c);
}
#endif*/


void ConstructGridCalogero(struct Grid *G) {
  AllocateWFGrid(G, 3);

  Warning("  Input parameter D means lambda of the CSM\n");
  Message("  Calogero parameter lambda = %.15" LE "\n", D);
  Message("  Energy per particle of Calogero gas %.15" LE "\n", 0.5*(1+(N-1)*D));

  G->min = 0;
  G->max = 1e10;
  G->max2 = G->max*G->max;
}

#ifdef TRIAL_CALOGERO
DOUBLE InterpolateExactU(DOUBLE x) {
  return D*Log(x);
}

DOUBLE InterpolateExactFp(DOUBLE x) {
  return D/x;
}

DOUBLE InterpolateExactE(DOUBLE x) {
#ifdef TRIAL_1D //1D Eloc = [-f"/f] + (f'/f)^2
  return D/(x*x);
#endif

#ifdef TRIAL_2D //2D Eloc = [-(f" + f'/r)/f] + (f'/f)^2
  return 0;
#endif

#ifdef TRIAL_3D //3D Eloc = [-(f"+2f'/r)/f] + (f'/f)^2
  return -D/(x*x);
#endif

}
#endif

/************************** Construct Grid Yukawa 2D **************************/
void ConstructGridYukawa2D(struct Grid *G, const DOUBLE Apar, const DOUBLE Bpar, const DOUBLE n) {

  AllocateWFGrid(G, 3);

#ifndef TRIAL_2D
  Error("  Trial wavefunction can be used only in 2D!\n");
#endif
  Atrial = pow(Apar*Lwf, Bpar);
  Btrial = 1./Exp(-Atrial/pow(Lhalfwf,Bpar));
  Cpar = 0;
  //Cpar = Atrial*Bpar*pow(Lhalf,-Bpar-1);
  //Btrial = 1;
 
  Message("  Yukawa trial wavefunction\n");

  G->min = 0;
  G->max = Lwf/2;
  G->max2 = G->max*G->max;
}

#ifdef TRIAL_YUKAWA_2D
DOUBLE InterpolateExactU(DOUBLE r) {
  if(r>Lhalfwf) r = Lhalfwf;
  if(r<0.5) r = 0.5;
  return Log(Btrial)-Atrial/pow(r,Bpar)+Log(0.5*(Exp(-Cpar*r)+Exp(-Cpar*(-Lwf+r))));
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  if(r>Lhalfwf) return 0.;
  if(r<0.5) return 0;
  return Atrial*Bpar*Btrial*pow(r,-Bpar-1)-Cpar;
}

DOUBLE InterpolateExactE(DOUBLE r) {
  if(r>Lhalfwf) return 0;
  if(r<0.5) return 0;
  return (Atrial*Bpar*Bpar+Cpar*pow(r,Bpar+1))*pow(r,-Bpar-2)*Btrial;
}
#endif

/************************** Construct Grid Yukawa 2D **************************
void ConstructGridYukawa3D(struct Grid *G, const DOUBLE Apar, const DOUBLE Bpar, const DOUBLE Cpar) {
// has a cusp at L/2
  AllocateWFGrid(G, 3);

#ifndef TRIAL_3D
  Error("  Trial wavefunction can be used only in 3D!\n");
#endif
  Message("  Yukawa 3D trial wavefunction\n");

  Atrial = Apar * Exp(-Bpar*Lhalf)/Lhalf*(1-Exp(-Lhalf/Cpar));

  G->min = 0;
  G->max = L/2;
  G->max2 = G->max*G->max;
}

#ifdef TRIAL_YUKAWA_3D
DOUBLE InterpolateExactU(DOUBLE r) {
  if(r>Lhalf) return 0;
  return Atrial - Apar * Exp(-Bpar*r)/r*(1-Exp(-r/Cpar));
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  if(r>Lhalf) return 0.;
  return -Apar*Exp(-Bpar*r)*(-Bpar*r*Cpar+Bpar*r*Cpar*Exp(-r/Cpar)-Cpar+Cpar*Exp(-r/Cpar)+Exp(-r/Cpar)*r)/(r*r)/Cpar;
}

DOUBLE InterpolateExactE(DOUBLE r) {
  if(r>Lhalf) return 0;
  return -Apar*Exp(-Bpar*r)/r*(-Bpar*Bpar*Cpar*Cpar+Exp(-r/Cpar)+2.0*Bpar*Cpar*Exp(-r/Cpar)+Bpar*Bpar*Cpar*Cpar*Exp(-r/Cpar))/(Cpar*Cpar);
}
#endif*/


/************************** Construct Grid Yukawa 2D **************************/
void ConstructGridYukawa3D(struct Grid *G, const DOUBLE Apar, const DOUBLE Bpar, const DOUBLE Cpar) {
// smoothed version
  AllocateWFGrid(G, 3);

#ifndef TRIAL_3D
  Error("  Trial wavefunction can be used only in 3D!\n");
#endif
  Message("  Yukawa 3D trial wavefunction (smooth)\n");

  Btrial = 2.; // parameter of smoothing
  G->min = 0;
  G->max = Lwf/2;
  G->max2 = G->max*G->max;
}

#ifdef TRIAL_YUKAWA_3D
DOUBLE InterpolateExactU(DOUBLE r) {
  static int first_time= ON;
  static DOUBLE offset;

  if(first_time == ON) {
    offset = Apar * Exp(-Bpar*Lhalfwf)/Lhalfwf*(1-Exp(-Lhalfwf/Cpar));
    first_time = OFF;
  }
#ifndef SYMMETRIZE_TRIAL_WF
  if(r>Lhalfwf) return 0;
#endif
  if(r == 0.) return 0.;
  return offset - Apar * Exp(-Bpar*r)/r*(1-Exp(-r/Cpar));
}

DOUBLE InterpolateExactFp(DOUBLE r) {
#ifndef SYMMETRIZE_TRIAL_WF
  if(r>Lhalfwf) return 0;
#endif
  return (Apar*(-r + Cpar*(-1 + Exp(r/Cpar))*(1 + Bpar*r)))/(Cpar*Exp((Bpar + 1/Cpar)*r)*(r*r));
}

DOUBLE InterpolateExactE(DOUBLE r) {
#ifndef SYMMETRIZE_TRIAL_WF
  if(r>Lhalfwf) return 0;
#endif
  return (Apar*(-1 + Bpar*Cpar*(-2 + Bpar*Cpar*(-1 + Exp(r/Cpar)))))/(Cpar*Cpar*Exp((Bpar + 1/Cpar)*r)*r);
}
#endif

#ifdef TRIAL_YUKAWA_3D_SYMM
DOUBLE InterpolateExactU(DOUBLE r) {
  if(r>Lhalfwf) return 0;
  return - Apar * Exp(-Bpar*r)/r*(1-Exp(-r/Cpar))*(1-Exp(-pow((Lhalfwf-r),Btrial)));
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  DOUBLE s1,s2,s3,s4,s5,s6;

  if(r>Lhalfwf) return 0.;
      s2 = Apar;
      s4 = Exp(-Bpar*r);
      s6 = Bpar*r*Cpar*Lhalfwf-Bpar*r*r*Cpar-Bpar*r*Cpar*Exp(-pow(Lhalfwf-r,1.0*
Btrial))*Lhalfwf+Bpar*r*r*Cpar*Exp(-pow(Lhalfwf-r,1.0*Btrial))-Bpar*r*Cpar*Exp(-r/
Cpar)*Lhalfwf+Bpar*r*r*Cpar*Exp(-r/Cpar)+Bpar*r*Cpar*Exp(-r/Cpar)*Exp(-pow(Lhalfwf-
r,1.0*Btrial))*Lhalfwf-Bpar*r*r*Cpar*Exp(-r/Cpar)*Exp(-pow(Lhalfwf-r,1.0*Btrial))+
Cpar*Lhalfwf-r*Cpar-Cpar*Exp(-pow(Lhalfwf-r,1.0*Btrial))*Lhalfwf;
      s5 = s6+Cpar*Exp(-pow(Lhalfwf-r,1.0*Btrial))*r-Cpar*Exp(-r/Cpar)*Lhalfwf+Cpar
*Exp(-r/Cpar)*r+Cpar*Exp(-r/Cpar)*Exp(-pow(Lhalfwf-r,1.0*Btrial))*Lhalfwf-Cpar*Exp(
-r/Cpar)*Exp(-pow(Lhalfwf-r,1.0*Btrial))*r-Exp(-r/Cpar)*r*Lhalfwf+Exp(-r/Cpar)*r*r+
Exp(-r/Cpar)*r*Exp(-pow(Lhalfwf-r,1.0*Btrial))*Lhalfwf-Exp(-r/Cpar)*r*r*Exp(-pow(
Lhalfwf-r,1.0*Btrial))+pow(Lhalfwf-r,1.0*Btrial)*Btrial*Exp(-pow(Lhalfwf-r,1.0*Btrial
))*r*Cpar-pow(Lhalfwf-r,1.0*Btrial)*Btrial*Exp(-pow(Lhalfwf-r,1.0*Btrial))*r*Cpar*
Exp(-r/Cpar);
      s3 = s4*s5;
      s1 = s2*s3;
      s2 = 1/(r*r)/Cpar/(Lhalfwf-r);
  return s1*s2;
}

DOUBLE InterpolateExactE(DOUBLE r) {
  DOUBLE s1,s2,s3,s4,s5;
  if(r>Lhalfwf) return 0;

 s1 = Apar*Exp(-Bpar*r)/r;      s5 = Exp(-r/Cpar)*Exp(-pow(Lhalfwf-r,1.0*Btrial))*Lhalfwf*Lhalfwf+Bpar*Bpar*r*r
*Cpar*Cpar+Exp(-r/Cpar)*r*r*Exp(-pow(Lhalfwf-r,1.0*Btrial))+2.0*Exp(-r/Cpar)*r*
Lhalfwf-Exp(-r/Cpar)*r*r-2.0*Bpar*r*r*Cpar*Exp(-r/Cpar)-2.0*Exp(-r/Cpar)*r*Exp(-
pow(Lhalfwf-r,1.0*Btrial))*Lhalfwf+4.0*Bpar*r*Cpar*Exp(-r/Cpar)*Lhalfwf-4.0*Bpar*r*
Cpar*Exp(-r/Cpar)*Exp(-pow(Lhalfwf-r,1.0*Btrial))*Lhalfwf;
      s4 = s5+2.0*Bpar*r*r*Cpar*Exp(-r/Cpar)*Exp(-pow(Lhalfwf-r,1.0*Btrial))+2.0*
pow(Lhalfwf-r,1.0*Btrial)*Btrial*Exp(-pow(Lhalfwf-r,1.0*Btrial))*r*Cpar*Exp(-r/Cpar
)+Bpar*Bpar*Cpar*Cpar*Lhalfwf*Lhalfwf-2.0*Bpar*Bpar*r*Cpar*Cpar*Lhalfwf-Exp(-r/Cpar)*
Lhalfwf*Lhalfwf-pow(Lhalfwf-r,1.0*Btrial)*Btrial*Exp(-pow(Lhalfwf-r,1.0*Btrial))*Cpar*
Cpar-Bpar*Bpar*r*r*Cpar*Cpar*Exp(-r/Cpar)-Bpar*Bpar*r*r*Cpar*Cpar*Exp(-pow(
Lhalfwf-r,1.0*Btrial))-Bpar*Bpar*Cpar*Cpar*Exp(-r/Cpar)*Lhalfwf*Lhalfwf;
      s5 = s4-Bpar*Bpar*Cpar*Cpar*Exp(-pow(Lhalfwf-r,1.0*Btrial))*Lhalfwf*Lhalfwf-pow
(pow(Lhalfwf-r,1.0*Btrial),2.0)*Btrial*Btrial*Exp(-pow(Lhalfwf-r,1.0*Btrial))*Cpar*
Cpar+pow(Lhalfwf-r,1.0*Btrial)*Btrial*Btrial*Exp(-pow(Lhalfwf-r,1.0*Btrial))*Cpar*
Cpar-2.0*Bpar*Cpar*Exp(-r/Cpar)*Lhalfwf*Lhalfwf-2.0*Cpar*Exp(-r/Cpar)*pow(Lhalfwf-r,
1.0*Btrial)*Btrial*Exp(-pow(Lhalfwf-r,1.0*Btrial))*Lhalfwf+pow(Lhalfwf-r,1.0*Btrial)*
Btrial*Exp(-pow(Lhalfwf-r,1.0*Btrial))*Cpar*Cpar*Exp(-r/Cpar)+2.0*Bpar*r*Cpar*
Cpar*Exp(-r/Cpar)*pow(Lhalfwf-r,1.0*Btrial)*Btrial*Exp(-pow(Lhalfwf-r,1.0*Btrial))+
2.0*Bpar*Cpar*Exp(-r/Cpar)*Exp(-pow(Lhalfwf-r,1.0*Btrial))*Lhalfwf*Lhalfwf;
      s3 = s5-2.0*Bpar*Cpar*Cpar*Exp(-r/Cpar)*pow(Lhalfwf-r,1.0*Btrial)*Btrial*
Exp(-pow(Lhalfwf-r,1.0*Btrial))*Lhalfwf+Bpar*Bpar*r*r*Cpar*Cpar*Exp(-r/Cpar)*Exp(-
pow(Lhalfwf-r,1.0*Btrial))-2.0*Bpar*Bpar*r*Cpar*Cpar*Exp(-r/Cpar)*Exp(-pow(Lhalfwf-
r,1.0*Btrial))*Lhalfwf+2.0*Bpar*Bpar*r*Cpar*Cpar*Exp(-r/Cpar)*Lhalfwf+Bpar*Bpar*
Cpar*Cpar*Exp(-r/Cpar)*Exp(-pow(Lhalfwf-r,1.0*Btrial))*Lhalfwf*Lhalfwf+2.0*Bpar*Bpar*
r*Cpar*Cpar*Exp(-pow(Lhalfwf-r,1.0*Btrial))*Lhalfwf-2.0*Bpar*r*Cpar*Cpar*pow(Lhalfwf-
r,1.0*Btrial)*Btrial*Exp(-pow(Lhalfwf-r,1.0*Btrial))+2.0*Bpar*Cpar*Cpar*pow(Lhalfwf
-r,1.0*Btrial)*Btrial*Exp(-pow(Lhalfwf-r,1.0*Btrial))*Lhalfwf-pow(Lhalfwf-r,1.0*
Btrial)*Btrial*Btrial*Exp(-pow(Lhalfwf-r,1.0*Btrial))*Cpar*Cpar*Exp(-r/Cpar)+pow(
pow(Lhalfwf-r,1.0*Btrial),2.0)*Btrial*Btrial*Exp(-pow(Lhalfwf-r,1.0*Btrial))*Cpar*
Cpar*Exp(-r/Cpar);
      s4 = 1/(Cpar*Cpar)/pow(Lhalfwf-r,2.0);
      s2 = s3*s4;

  return s1*s2;
}
#endif

#ifdef TRIAL_HS2D_ZERO_ENERGY
// zero energy HD: f(r) = const ln r
DOUBLE InterpolateExactU(DOUBLE r) {
  if(r<1) return -100;
  if(r>Lhalf) return 0;
  return Log(Log(r)/Log(Lhalf));
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  if(r<1) return 0.;
  if(r<1 || r>Lhalf) return 0.;
  return 1./(r*Log(r));
}

DOUBLE InterpolateExactE(DOUBLE r) {
  DOUBLE Fp;
  if(r<1 || r>Lhalf) return 0.;
  Fp = 1./(r*Log(r));
  return Fp*Fp;
}
#endif

/************************** Construct Grid HS 2D Finite Energy *************************/
#ifdef TRIAL_HS2D_FINITE_ENERGY
// zero energy HD: f(r) = const ln r
void ConstructGridHS2DFiniteEnergy(struct Grid *G) {
  DOUBLE k, kmin, kmax;
  DOUBLE y, ymin, ymax;
  DOUBLE precision = 1e-12;
  int search = ON;
  DOUBLE C=0.;
  DOUBLE R;

  Message("  Jastrow term: ConstructGridHS2DFiniteEnergy trial wavefunction\n");
  if(Rpar == 0.) {
    Warning("  Rpar -> Lwf/2\n", Lhalfwf);
    Rpar = Lhalfwf;
  }
  R = Rpar;

  AllocateWFGrid(G, 3);

  k = kmin = 1e-6/R;
  Atrial =  BesselY0(k*a)/(BesselJ0(k*R)*BesselY0(k*a)-BesselJ0(k*a)*BesselY0(k*R));
  Btrial = -BesselJ0(k*a)/(BesselJ0(k*R)*BesselY0(k*a)-BesselJ0(k*a)*BesselY0(k*R));

  // Solve J1(kR) Y0(ka) - J0(ka) Y0(kR) = 0
  Message("Fixing kR (1)... ");
  ymin = BesselJ1(k*R)*BesselY0(k*a) - BesselY1(k*R)*BesselJ0(k*a);
  while(search) {
    k = kmax = kmin + 0.1*kmin;
    Atrial =  BesselY0(k*a)/(BesselJ0(k*R)*BesselY0(k*a)-BesselJ0(k*a)*BesselY0(k*R));
    Btrial = -BesselJ0(k*a)/(BesselJ0(k*R)*BesselY0(k*a)-BesselJ0(k*a)*BesselY0(k*R));
    y = BesselJ1(k*R)*BesselY0(k*a) - BesselY1(k*R)*BesselJ0(k*a);

    if((ymin-C)*(y-C)<0) {
      search = OFF;
    }
    else {
      kmin = k;
      ymin = y;
    }
  }

  search = ON;
  Message("done\n  Fixing kR (2)... ");
  while(search) {
    k = (kmin+kmax)/2.;
    Atrial =  BesselY0(k*a)/(BesselJ0(k*R)*BesselY0(k*a)-BesselJ0(k*a)*BesselY0(k*R));
    Btrial = -BesselJ0(k*a)/(BesselJ0(k*R)*BesselY0(k*a)-BesselJ0(k*a)*BesselY0(k*R));
    y = BesselJ1(k*R)*BesselY0(k*a) - BesselY1(k*R)*BesselJ0(k*a);

    if((y-C)*(ymin-C) < 0) {
      kmax = k;
      ymax = y;
    } else {
      kmin = k;
      ymin = y;
    }

    if(fabs(y)<precision) {
      k = (kmax*(ymin-C)-kmin*(ymax-C)) / (ymin-ymax);
      Atrial =  BesselY0(k*a)/(BesselJ0(k*R)*BesselY0(k*a)-BesselJ0(k*a)*BesselY0(k*R));
      Btrial = -BesselJ0(k*a)/(BesselJ0(k*R)*BesselY0(k*a)-BesselJ0(k*a)*BesselY0(k*R));
      y = BesselJ1(k*R)*BesselY0(k*a) - BesselY1(k*R)*BesselJ0(k*a);
      search = OFF;
    }
  }
  Message("done\n");

  Atrial =  BesselY0(k*a)/(BesselJ0(k*R)*BesselY0(k*a)-BesselJ0(k*a)*BesselY0(k*R));
  Btrial = -BesselJ0(k*a)/(BesselJ0(k*R)*BesselY0(k*a)-BesselJ0(k*a)*BesselY0(k*R));

  Message("  k = %" LE " , error %" LE " \n", k, fabs(y));
  Message("    A = %" LE " \n", Atrial);
  Message("    B = %" LE " \n", Btrial);

  sqE = k;
  sqE2 = sqE*sqE;

  G->min = 0;
  G->max = L/2;
  G->max2 = G->max*G->max;
}

DOUBLE InterpolateExactU(DOUBLE r) {
  if(r>Rpar)
    return 0.;
  else if(r>a) {
    return Log(Atrial*BesselJ0(sqE*r)+Btrial*BesselY0(sqE*r));
  }
  else
    return -100;
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  if(r>Rpar || r<a)
    return 0.;
  else
    return -sqE*(Atrial*BesselJ1(sqE*r)+Btrial*BesselY1(sqE*r))/(Atrial*BesselJ0(sqE*r)+Btrial*BesselY0(sqE*r));
}

DOUBLE InterpolateExactE(DOUBLE r) {
  DOUBLE Fp;
  if(r<a)
    return -100;
  else
    if( r>Lhalf) return 0.;

  Fp = -sqE*(Atrial*BesselJ1(sqE*r)+Btrial*BesselY1(sqE*r))/(Atrial*BesselJ0(sqE*r)+Btrial*BesselY0(sqE*r));;
  return sqE*sqE+Fp*Fp;
}
#endif

#ifdef  TRIAL_HS2D_XING
DOUBLE InterpolateExactU(DOUBLE r) {
  if(r<1.) return -100./r;
  if(r>Lhalfwf) return 0;
  return Log(tanh(Apar*(pow(r,Bpar)-1))) - Log(tanh(Apar*(pow(Lhalfwf,Bpar)-1)));
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  if(r<1) return 0;
  if(r>Lhalfwf) return 0;
  return 2*Apar*Bpar*pow(r,Bpar-1)/sinh(2*Apar*(pow(r,Bpar)-1));
}

DOUBLE InterpolateExactE(DOUBLE r) {
  DOUBLE arg;

  if(r<1) return 0.;
  if(r>Lhalfwf) return 0;

  arg = 2*Apar*(pow(r,Bpar)-1);
  return 2*Apar*Bpar*Bpar*pow(r,Bpar-2)/sinh(arg)*(-1+2*Apar*pow(r,Bpar)/tanh(arg));
}
#endif

#ifdef  TRIAL_HS2D
// finite energy soluton 0<Apar<1, Apar is close to 1
DOUBLE InterpolateExactU(DOUBLE r) {
  static int first_time =ON;
  static DOUBLE shift;

  if(first_time ==ON) {
    if(Apar*(Lhalfwf-1)<1./0.04) // ln(1-e^{-1/x}) <10^-10, x<0.04
      shift = 0.5*(1+Apar)*Lhalfwf - Apar*Lhalfwf - Log(1-Exp(Apar*(1-Lhalfwf)));
    else 
      shift = 0.5*(1+Apar)*Lhalfwf - Apar*Lhalfwf;
    first_time = OFF;
  }
  if(r<1) return -100;
  if(r>Lhalfwf) return 0;
  if(Apar*(r-1)<25)
    //Log(Exp(-0.5*(1+Apar)*r)*(Exp(Apar*r)-Exp(Apar)))+shift;
    return -0.5*(1+Apar)*r+Apar*r+Log(1-Exp(Apar*(1.-r)))+shift;
  else
    return -0.5*(1+Apar)*r+Apar*r+shift;
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  if(r<1) return 0;
  if(r>Lhalfwf) return 0;
  return -0.5*(1+Apar/tanh(0.5*Apar*(1-r)));
}

DOUBLE InterpolateExactE(DOUBLE r) {
  if(r<1) return 10*E;
  if(r>Lhalfwf) return 0;
  return ((1 + Apar)*Exp(2*Apar) - (-1 + Apar)*Exp(2*Apar*r) + 2*Exp(Apar + Apar*r)*(-1 + Apar*Apar*r))/(2.*pow(Exp(Apar) - Exp(Apar*r),2)*r);
}
#endif

#ifdef TRIAL_2D_ZERO_RANGE_BOUND_STATE
DOUBLE InterpolateExactU(DOUBLE r) {
  static int first_time = ON;

  if(first_time) {
    trial_k = 2./(exp(0.57721566490153286060651209008240243104215933593992)*a);
    trial_k2 = trial_k*trial_k;
    first_time = OFF;
  }

  return Log(BesselK0(trial_k*r));
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  return -trial_k*BesselK1(trial_k*r) / BesselK0(trial_k*r);
}

DOUBLE InterpolateExactE(DOUBLE r) {
  DOUBLE Fp;
  Fp = -trial_k*BesselK1(trial_k*r) / BesselK0(trial_k*r);
  return -trial_k2 + Fp*Fp;
}
#endif

#ifdef TRIAL_2D_ZERO_RANGE_K0_SCREENING
DOUBLE InterpolateExactU(DOUBLE r) {
  static int first_time = ON;

  if(first_time) {
    trial_k = 2./(exp(0.57721566490153286060651209008240243104215933593992)*a);
    first_time = OFF;
  }

  return Log(fabs(BesselK0(trial_k*r))) - r*r/(b*(Apar+r));
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  return -trial_k*BesselK1(trial_k*r)/BesselK0(trial_k*r) - (r*(2.*Apar+r)/(b*(Apar+r)*(Apar+r)));
}

DOUBLE InterpolateExactE(DOUBLE r) {
  DOUBLE K0,K1,K2,k;
  k = trial_k;
  K0 = BesselK0(trial_k*r);
  K1 = BesselK1(trial_k*r);
  //K2 = BesselK2(trial_k*r);
  return (2*Apar*Apar + Apar*(Apar + r) + (Apar + r)*(Apar + r) - k*k*b*(Apar + r)*(Apar + r)*(Apar + r))/(b*(Apar + r)*(Apar + r)*(Apar + r)) + (k*k*K1*K1)/(K0*K0);
}
#endif

#ifdef TRIAL_SPHERE_ZERO_RANGE_BOUND_STATE
//sphere: Eloc = +f''/f
DOUBLE InterpolateExactU(DOUBLE r) {
  static int first_time = ON;

  if(first_time) {
    trial_k = 2./(exp(0.57721566490153286060651209008240243104215933593992)*a);
    first_time = OFF;
  }

  return Log(BesselK0(trial_k*r));
}

DOUBLE InterpolateExactFp(DOUBLE r) { // f'/f
  return -trial_k*BesselK1(trial_k*r) / BesselK0(trial_k*r);
}

DOUBLE InterpolateExactE(DOUBLE r) { //sphere: Eloc =  f''/f
  return trial_k*(trial_k + BesselK1(trial_k*r) / (r*BesselK0(trial_k*r)));
}
#endif

#ifdef TRIAL_SPHERE_ZERO_RANGE_FREE_STATE // (sphere) abs(ln(r/a))
//sphere: Eloc = +f''/f
DOUBLE InterpolateExactU(DOUBLE r) {
  return Log(a/r);
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  return -1./ (r*Log(a/r));
}

DOUBLE InterpolateExactE(DOUBLE r) {
  DOUBLE Fp;
  Fp = 1./ (r*Log(a/r));
  return Fp*Fp;
}
#endif

#ifdef TRIAL_SPHERE_ZERO_RANGE_LOG_MATCH // Rpar: ln(r/a), r<Rpar; ln(Rpar/a) * exp((r/Rpar-1) / ln(r/a)), Eloc(r) = - f"(r)/f(r)
//sphere: Eloc = +f''/f
DOUBLE InterpolateExactU(DOUBLE r) {
  if(r>a) Error(" w.f. assumes that r<a\n");
  if(r<Rpar)
    return Log(Log(a/r));
  else
    return Log(Log(a/Rpar)) + (r/Rpar-1.) / Log(Rpar/a);
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  DOUBLE logf;
  logf = Log(r/a);
  if(r<Rpar)
    return 1./(r*logf);
  else
    return 1./(Rpar*Log(Rpar/a));
}

DOUBLE InterpolateExactE(DOUBLE r) {
  DOUBLE logf;
  if(r<Rpar) {
    logf = r*Log(r/a);
    return -1. / (logf*logf);
  }
  else {
    logf = Rpar*Log(Rpar/a);
    return 1. / (logf*logf);
  }
}

#endif

#ifdef TRIAL_SPHERE_ZERO_RANGE_K0_SCREENING
DOUBLE InterpolateExactU(DOUBLE r) {
  static int first_time = ON;

  if(first_time) {
    trial_k = 2./(exp(0.57721566490153286060651209008240243104215933593992)*a);
    first_time = OFF;
  }

  return Log(fabs(BesselK0(trial_k*r))) - r*r/(b*(Apar+r));
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  return -trial_k*BesselK1(trial_k*r)/BesselK0(trial_k*r) - (r*(2.*Apar+r)/(b*(Apar+r)*(Apar+r)));
}

DOUBLE InterpolateExactE(DOUBLE r) {
  DOUBLE K0,K1,K2,k;
  k = trial_k;
  K0 = BesselK0(trial_k*r);
  K1 = BesselK1(trial_k*r);
  K2 = BesselK2(trial_k*r);

  return ((Apar*Apar*Apar*Apar*b*b*k*k + 4*Apar*(2 + b*b*k*k)*r*r*r + (2 + b*b*k*k)*r*r*r*r + 4*Apar*Apar*Apar*b*(-1 + b*k*k*r) + 2*Apar*Apar*r*(4*r + b*(-2 + 3*b*k*k*r)))*K0 + b*k*(Apar + r)*(Apar + r)*(4*r*(2*Apar + r)*K1 + b*k*(Apar + r)*(Apar + r)*K2))/(2.*b*b*(Apar + r)*(Apar + r)*(Apar + r)*(Apar + r)*K0);
}
#endif

#ifdef TRIAL_SPHERE_ZERO_RANGE_LOG_SCREENING // (sphere) ln(a/r) * exp(-r/b), Eloc(r) = +f"(r)/f(r)
DOUBLE InterpolateExactU(DOUBLE r) {
  return Log(Log(a/r)) - r/b;
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  return -1./(r*Log(a/r)) - 1./b;
}

DOUBLE InterpolateExactE(DOUBLE r) {
  return 1./(b*b) + (b+2.*r)/(b*r*r*Log(a/r));
}
#endif

/************************** Construct Grid Yukawa 2D **************************/
void ConstructGridCalogero2D(struct Grid *G, const DOUBLE Apar) {

  AllocateWFGrid(G, 3);

#ifndef TRIAL_2D
  Error("  Trial wavefunction can be used only in 2D!\n");
#endif

  Message("  Calogero 2D wf. A r^2, r<Rpar L/2; Exp(B-Apar*(1./r+1./(L-r))) r>Rpar*Lhalf)\n");
  Message("  Rpar %" LG "  [Lwf/2] -> %" LG "  [R]\n", Rpar, Rpar*Lhalfwf);
  Rpar *= Lhalfwf;

  Cpar = 2*(Lwf*Lwf*Rpar-2*Rpar*Rpar*Lwf+Rpar*Rpar*Rpar)/(Lwf*(Lwf-2*Rpar));
  Atrial = Exp(4*Cpar/Lwf - Cpar*(1./Rpar+1./(L-Rpar)))/(Rpar*Rpar);
  Ctrial = 4.*Cpar / Lwf;

  G->min = 0;
  G->max = Lwf/2;
  G->max2 = G->max*G->max;
}

#ifdef TRIAL_CALOGERO_2D
DOUBLE InterpolateExactU(DOUBLE r) {

  if(r>Lhalfwf) return 0;

  if(r<Rpar) {
    return Log(Atrial*r*r);
  }
  else {
    return Ctrial-Cpar*(1/r+1/(Lwf-r));
  }
}

DOUBLE InterpolateExactFp(DOUBLE r) {

  if(r>Lhalfwf) return 0;

  if(r<Rpar) {
    return 2./r;
  }
  else {
    return Cpar*(1/(r*r)-1/((Lwf-r)*(Lwf-r)));
  }
}

DOUBLE InterpolateExactE(DOUBLE r) {

  if(r>Lhalfwf) return 0;

  if(r<Rpar) {
    return 0;
  }
  else {
    return Cpar*Lwf*(Lwf*Lwf-3.*Lwf*r+4.*r*r)/((Lwf-r)*(Lwf-r)*(Lwf-r)*r*r*r);
  }
}
#endif

/************************** Construct Grid Decay 2D ***********************/
#ifdef TRIAL_DECAY
DOUBLE InterpolateExactU(DOUBLE r) {
  if(r>Lhalfwf) return 0;
  return -Apar/(1+Bpar*r*r) + Apar/(1+Bpar*Lhalfwf*Lhalfwf);
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  DOUBLE d;
  if(r>Lhalfwf) return 0.;
  d = Bpar*r*r + 1;
  return 2*Apar*Bpar*r/(d*d);
}

DOUBLE InterpolateExactE(DOUBLE r) {
  DOUBLE d;
  if(r>Lhalfwf) return 0.;
  d = Bpar*r*r + 1;
  return 4*Apar*Bpar*(Bpar*r*r-1.)/(d*d*d);
}
#endif

void ConstructGridDecayPhonon(struct Grid *G, const DOUBLE Apar, const DOUBLE Bpar) {
  DOUBLE r;

  AllocateWFGrid(G, 3);

  Message("  Jatrow term: GridDecayPhonon\n");
  Warning("  Rpar %" LG "  [Lwf/2] -> %" LG "  [R]\n", Rpar, Rpar*Lhalfwf);
  Rpar *= Lhalfwf;

  r = Rpar;
  Cpar = 2*Apar*Bpar*r/((Bpar*r*r + 1)*(Bpar*r*r + 1))/(1./(r*r)-1./((Lwf-r)*(Lwf-r)));
  Message("  Construct Grid decay 2D. Cpar = %" LG "\n", Cpar);
  Btrial = 4*Cpar/Lwf;
  //Btrial = Cpar = 0;
  Atrial=Apar/(1+Bpar*r*r)+Btrial-Cpar*(1./r+1./(Lwf-r));

  G->max = Lhalfwf;
  G->max2 = G->max*G->max;
}

#ifdef TRIAL_DECAY_PHONON
DOUBLE InterpolateExactU(DOUBLE r) {
  if(r>Lhalfwf) 
    return 0;
  else if(r<Rpar)
    return Atrial-Apar/(1+Bpar*r*r);
  else
    return Btrial-Cpar*(1./r+1./(Lwf-r));
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  DOUBLE d;
  if(r>Lhalfwf) 
    return 0.;
  else if(r<Rpar) {
    d = Bpar*r*r + 1;
    return 2*Apar*Bpar*r/(d*d);
  }
  else {
    return Cpar*(1./(r*r)-1./((Lwf-r)*(Lwf-r)));
  }
}

DOUBLE InterpolateExactE(DOUBLE r) {
  DOUBLE d;
  if(r>Lhalfwf)
    return 0.;
  else if(r<Rpar) {
    d = Bpar*r*r + 1;
    return 4*Apar*Bpar*(Bpar*r*r-1.)/(d*d*d);
  }
  else
    return Cpar*Lwf*(Lwf*Lwf-3.*Lwf*r+4.*r*r)/(r*r*r*(Lwf-r)*(Lwf-r)*(Lwf-r));
}
#endif

/************************** Construct Grid Coulomb 1D Phonon **************/
void ConstructGridCoulomb1DPhonon(struct Grid *G) {
  AllocateWFGrid(G, 3);

  Message("  Jastrow term: Coulomb 1D phonon.\n");

  // here Rpar is units of [a]
  Warning("  Changing Rpar -> Rpar Lhalfwf (%" LG "  -> %" LG " )\n", Rpar, Rpar*Lhalfwf);
  Rpar *= Lhalfwf; // do not change Rpar earlier!

  // here Rpar is in units of Lhalf
  Btrial = (Lwf*(Sqrt(Rpar)*BesselI0(2.*Sqrt(Rpar))+BesselI1(2.*Sqrt(Rpar))+Sqrt(Rpar)*BesselIN(2,2.*Sqrt(Rpar)))*Tg((PI*Rpar)/Lwf))/(2.*PI*Rpar*BesselI1(2.*Sqrt(Rpar)));
  Atrial = pow(Sin(PI*Rpar/Lwf),Btrial)/(Sqrt(Rpar)*BesselI1(2.*Sqrt(Rpar)));
  Message("  Atrial = %" LG " , Btrial = %" LG "\n", Atrial, Btrial);

  if(Rpar>Lhalfwf) Error("  Rpar must be smaller than L/2!\n");

  G->max = Lhalfwf;
  G->max2 = G->max*G->max;
}

#ifdef TRIAL_COULOMB_1D_PHONON
DOUBLE InterpolateExactU(DOUBLE r) {
  if(r>Lhalfwf) 
    return 0;
  else if(r<Rpar)
    return Log(Atrial*Sqrt(r)*BesselI1(2.*Sqrt(r)));
  else
    return Log(Sin(PI/Lwf*r))*Btrial;
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  DOUBLE I0,I1,I2,arg;
  if(r>Lhalfwf) 
    return 0.;
  else if(r<Rpar) {
    arg = 2.*Sqrt(r);
    I0 = BesselI0(arg);
    I1 = BesselI1(arg);
    I2 = BesselIN(2,arg);
    return (Sqrt(r)*(I0+I2)+I1)/(2.*r*I1);
  }
  else {
    return Btrial*PI/(Lwf*Tg(PI/Lwf*r));
  }
}

DOUBLE InterpolateExactE(DOUBLE r) {
  //1D Eloc = [-f"/f] + (f'/f)^2
  DOUBLE c;
  if(r>Lhalfwf)
    return 0.;
  else if(r<Rpar) {
    c = InterpolateExactFp(r);
    return -1/r + c*c;
  }
  else {
    c = 1 / Tg(PI*r/Lwf);
    return Btrial*PI*PI/(Lwf*Lwf)*(1+c*c);
  }
}
#endif

/************************** Construct Grid Coulomb 1D Phonon **************/
void ConstructGridDipole1DPhonon(struct Grid *G) {
  AllocateWFGrid(G, 3);

  Message("  Jastrow term: Dipole 1D phonon.\n");

  Warning("  Changing Rpar -> Rpar Lhalfwf (%" LG "  -> %" LG " )\n", Rpar, Rpar*Lhalfwf);
  Rpar *= Lhalfwf;
  if(Rpar>Lhalfwf) Error("  Rpar must be smaller than L/2!\n");

  Btrial = (Lwf*(BesselK0(2/Sqrt(Rpar)) + Sqrt(Rpar)*BesselK1(2/Sqrt(Rpar)) + BesselK2(2/Sqrt(Rpar)))*Tg((PI*Rpar)/Lwf))/
   (2.*PI*pow(Rpar,1.5)*BesselK1(2/Sqrt(Rpar)));

  Atrial = pow(Sin(PI*Rpar/Lwf),Btrial)/(Sqrt(Rpar)*BesselK1(2./Sqrt(Rpar)));
  Message("  Atrial = %" LG " , Btrial = %" LG "\n", Atrial, Btrial);

  G->max = Lhalfwf;
  G->max2 = G->max*G->max;
}

#ifdef TRIAL_DIPOLE_1D_PHONON
DOUBLE InterpolateExactU(DOUBLE r) {
  if(r>Lhalfwf) 
    return 0;
  if(r<0.01)
    return -0.12078223763524543 - 2./Sqrt(r) + 0.1875*Sqrt(r) - 0.046875*r + 0.0205078125*r*Sqrt(r) - 0.013183593749999998*r*r + Log(Atrial) + 0.75*Log(r);
  else if(r<Rpar)
    return Log(Atrial*Sqrt(r)*BesselK1(2./Sqrt(r)));
  else
    return Log(Sin(PI/Lwf*r))*Btrial;
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  DOUBLE K0,K1,K2,arg;
  if(r>Lhalfwf) 
    return 0.;
  else if(r<Rpar) {
    arg = 2./Sqrt(r);
    K0 = BesselK0(arg);
    K1 = BesselK1(arg);
    K2 = BesselK2(arg);
    return ((K0+K2)/Sqrt(r)+K1)/(2.*r*K1);
  }
  else {
    return Btrial*PI/(Lwf*Tg(PI/Lwf*r));
  }
}

DOUBLE InterpolateExactE(DOUBLE r) {
  //1D Eloc = [-f"/f] + (f'/f)^2
  DOUBLE c;
  if(r>Lhalfwf)
    return 0.;
  else if(r<Rpar) {
    c = InterpolateExactFp(r);
    return -1/(r*r*r) + c*c;
  }
  else {
    c = 1 / Tg(PI*r/Lwf);
    return Btrial*PI*PI/(Lwf*Lwf)*(1+c*c);
  }
}
#endif

#ifdef TRIAL_DIPOLE_1D_TRAP
DOUBLE InterpolateExactU(DOUBLE r) {
  r /= D;
  if(r<0.01)
    return -0.12078223763524543 - 2./Sqrt(r) + 0.1875*Sqrt(r) - 0.046875*r + 0.0205078125*r*Sqrt(r) - 0.013183593749999998*r*r + 0.75*Log(r);
  else
    return Log(Sqrt(r)*BesselK1(2./Sqrt(r)));
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  DOUBLE K0,K1,K2,arg;

  arg = 2.*Sqrt(D/r);
  K0 = BesselK0(arg);
  K1 = BesselK1(arg);
  K2 = BesselK2(arg);
  return (1.+(K0+K2)*Sqrt(D/r)/K1)/(2.*r);
}

DOUBLE InterpolateExactE(DOUBLE r) {
  //1D Eloc = [-f"/f] + (f'/f)^2
  DOUBLE c;

  c = InterpolateExactFp(r);
  return -D/(r*r*r) + c*c;
}
#endif


#ifdef TRIAL_ABSENT12
DOUBLE InterpolateExactU(DOUBLE r) {
  return 0;
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  return 0;
}

DOUBLE InterpolateExactE(DOUBLE r) {
  return 0;
}
#endif

#ifdef TRIAL_ABSENT11
DOUBLE InterpolateExactU11(DOUBLE r) {
  return 0;
}

DOUBLE InterpolateExactFp11(DOUBLE r) {
  return 0;
}

DOUBLE InterpolateExactE11(DOUBLE r) {
  return 0;
}
#endif

#ifdef TRIAL_ABSENT22
DOUBLE InterpolateExactU22(DOUBLE r) {
  return 0;
}

DOUBLE InterpolateExactFp22(DOUBLE r) {
  return 0;
}

DOUBLE InterpolateExactE22(DOUBLE r) {
  return 0;
}
#endif

/************************** Construct Grid Coulomb 1D Trap ****************/
void ConstructGridCoulomb1DTrap(struct Grid *G) {
  AllocateWFGrid(G, 3);

  Message("  2-body Jastrow term: Coulomb 1D trap\n");

  Btrial = ((Sqrt(Lwf*Rpar)*BesselI0(Sqrt(2*Lwf*Rpar)) + Sqrt(2)*BesselI1(Sqrt(2*Lwf*Rpar)) + Sqrt(Lwf*Rpar)*BesselIN(2,Sqrt(2*Lwf*Rpar)))*Tg((PI*Rpar)/2.))/(Sqrt(2)*PI*Rpar*BesselI1(Sqrt(2*Lwf*Rpar)));
  Atrial = Sqrt(2)*pow(Sin(PI*Rpar/2.),Btrial)/(Sqrt(Lwf*Rpar)*BesselI1(Sqrt(2*Lwf*Rpar)));
  Message("  Atrial = %" LG " , Btrial = %" LG "\n", Atrial, Btrial);

  //Warning("  Changing Rpar -> Rpar Lhalfwf (%" LG "  -> %" LG " )\n", Rpar, Rpar*Lhalfwf);
  //Rpar *= Lhalfwf; // do not change Rpar earlier!

  G->max = N*10;
  G->max2 = G->max*G->max;
}

#ifdef TRIAL_COULOMB_1D_TRAP
DOUBLE InterpolateExactU(DOUBLE r) {
  return Log(Sqrt(r)*BesselI1(2.*Sqrt(D*r)));
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  DOUBLE I0,I1,I2,arg;
  arg = 2.*Sqrt(D*r);
  I0 = BesselI0(arg);
  I1 = BesselI1(arg);
  I2 = BesselIN(2,arg);
  return (Sqrt(D*r)*(I0+I2)+I1)/(2.*r*I1);
}

DOUBLE InterpolateExactE(DOUBLE r) {
  //1D Eloc = [-f"/f] + (f'/f)^2
  DOUBLE fp;
  fp = InterpolateExactFp(r);
  return -D/r + fp*fp;
}
#endif

/************************** Construct Grid Coulomb 2D Trap ****************/
void ConstructGridCoulomb2DTrap(struct Grid *G) {
  AllocateWFGrid(G, 3);

  Message("  2-body Jastrow 11 term: Coulomb 2D trap: f2(r) = I_0[2(Dr)^1/2]\n");

  G->max = N*10;
  G->max2 = G->max*G->max;
}

#ifdef TRIAL_COULOMB_2D_TRAP
DOUBLE InterpolateExactU11(DOUBLE r) {
// ln f
  return Log(BesselI0(2*Sqrt(D*r)));
}

DOUBLE InterpolateExactFp11(DOUBLE r) {
// f' / f
  DOUBLE arg;
  arg = 2.*Sqrt(D*r);
  return Sqrt(D/r)*BesselI1(arg)/BesselI0(arg);
}

DOUBLE InterpolateExactE11(DOUBLE r) {
  //2D Eloc = -[f"+f'/r]/f + (f'/f)^2
  DOUBLE c;
  c = InterpolateExactFp11(r);
  return -D/r + c*c;
}
#endif

/************************** Construct Grid Coulomb 2D Trap ****************/
void ConstructGridHardRods3D(struct Grid *G) {
  AllocateWFGrid(G, 3);

  Message("Hard-rod type (1D,2D,3D) trial w.f.: sin[pi (r-1)/(L-2)]^Apar\n");
  Message("   Apar = %lf\n", Apar);

  Atrial = PI/(L-2.);

  G->max = Lhalf;
  G->max2 = G->max*G->max;
}

#ifdef TRIAL_HARD_RODS_3D
DOUBLE InterpolateExactU(DOUBLE r) {
  // ln f
  if(r<a)
    return -100;
  else
    return Apar*Log(sin(Atrial*(r-a)));
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  // f' / f
  if(r<a)
    return 100;
  else
    return Apar*Atrial/tan(Atrial*(r-a));
}

DOUBLE InterpolateExactE(DOUBLE r) {
  DOUBLE sin_1;
  //Eloc = -[f"+(D-1)(f'/r]/f + (f'/f)^2
  if(r<a) {
    return 0;
  }
  else {
    sin_1 = 1./sin(Atrial*(r-a));
    return (Apar*Atrial*((DIMENSION-1.)/tan(Atrial*(a-r)) + Atrial*r*sin_1*sin_1))/r;
  }
}
#endif

/***************************** Simpson Integration ***************************/
DOUBLE IntegrateSimpson(DOUBLE xmin, DOUBLE xmax, DOUBLE (*f)(DOUBLE x), long int Npoints) {
// 1/3, 4/3, 2/3, ..., 2/3, 4/3, 1/3
  DOUBLE dx, x;
  DOUBLE Int = 0.;
  long int i;
  long int N;

#ifdef HARD_SPHERE
  if(xmin<a) {
    Warning(" Integrate Simpson xmin<a call");
    xmin = a*1.01;
  }
#endif

  N = (Npoints/2)*2+1; // number of points must be odd
  dx = (xmax-xmin)/ (DOUBLE) (N-1);
  Int = f(xmin);
  Int += 4 * f(xmin+dx);

  for(i=2; i<N-2; i += 2) {
    x = xmin + dx*(DOUBLE) i;
    Int += 2. * f(x);
    Int += 4. * f(x+dx);
  }

  Int += f(x+2.*dx);

  return Int * dx / 3.;
}

/************************** Construct Grid Spline 2 body solution  **************************/
// construct 2body solution for use with spline or linear interpolation
void ConstructGridSpline2Body(struct Grid *G, DOUBLE (*InteractionEnergy)(DOUBLE x)) {
  // choose boundary condition of the scattering problem
  // 0: zero boundary condition f(0) = 0;
  // 1: f'(0) = 0; 
  int flag_left_boundary_condition = 0; // 0 for f(0) = 0; 1 for f'(0) = 0; 2 for f'(0)/f(0) = -1/a
  int flag_right_boundary_condition = 0; // 0 for f'(R) = 0; 2 for scattering solution asymptotic
  int size = 1000;
  DOUBLE E=0.; // scattering energy
  DOUBLE xmin, xmax, ymin, ymax, x, y;
  int search = ON;
  DOUBLE precision = 1e-14; //1e-8
  int Nitermax = 1000000; // maximal allowed number of iterations
  DOUBLE energy_step = 1.1; // 1.01 to find ymin*ymax<0
  int iter=0;
  int i;
  int function_was_negative = OFF;
  int flag_lnf_instead_of_f = OFF; // ON - has some problems in 1D solve Scroedinger equation for ln f instead of f

  Message("  Jastrow term will be constructed as 2-body scattering solution with for interaction potential\n");
  Message("  Spline size grid_trial=%i\n", grid_trial);
  size = grid_trial;
  AllocateWFGrid(G, size);

  if(Epotcutoff<0) {// 1e3 maximal allowed potential energy for V(xmin)
    Epotcutoff = 1e3;
  }
  Message("  maximal allowed potential energy for V(xmin) = Epotcutoff = %e\n", Epotcutoff);

  G->size = size;

 // choose boundary condition of the scattering problem in main.h
#ifdef JASTROW_LEFT_BOUNDARY_ZERO
  // 0: zero boundary condition f(0) = 0;
  flag_left_boundary_condition = 0;
  Message("  zero boundary condition f(0) = 0 will be used\n");
#endif
#ifdef JASTROW_LEFT_BOUNDARY_ZERO_DERIVATIVE
  // 1: f'(0) = 0; 
  flag_left_boundary_condition = 1;
  Message("  zero derivative condition f'(0) = 0 will be used\n");
#endif
#ifdef JASTROW_LEFT_BOUNDARY_PSEUDOPOTENTIAL
  // 2: f'(0)/f(0) = -1/a; 
  flag_left_boundary_condition = 2;
  Message("  pseudopotential derivative condition f'(0)/f(0) = -1/a will be used\n");
#endif

#ifdef TRIAL_TWO_BODY_SCATTERING_NUMERICAL
#ifdef SPINFULL
#ifdef SPINFULL_TUNNELING
  Message("  maximal distance for the spline Rpar12 = %lf\n", Rpar12);
  G->max = Rpar12;
#else
#ifdef BC_ABSENT
  Message("  maximal distance for the spline Rpar12 = %lf\n", Rpar12);
  G->max = Rpar12;
#else // no tunneling, with PBC
  if(Rpar12>1) Error("Check failed: 0<Rpar12=%lf<1 (no SPINFULL_TUNNELING)", Rpar12);
  Message("  maximal distance for the spline Rpar12 x L/2 = %lf (Rpar12 = %lf)\n", Rpar12*Lhalf, Rpar12);
  G->max = Lhalf*Rpar12;
#endif
#endif // not TUNELLING
#else // not SPINFULL
  Message("  maximal distance for the spline Rpar x L/2 = %lf (Rpar = %lf)\n", Rpar*Lhalf, Rpar);
  G->max = Lhalf*Rpar;
  if(G->max>Lhalf) Warning("  Spline maximal distance %lf is larger then the box size L/2 = %lf\n", G->max, Lhalf);
  //G->max = Lhalf;
#endif
#endif

  if(G->max<1e-8) Error("  Cannot solve with G->max = 0\n");

  if(DIMENSION>1)
    G->min = G->max/(DOUBLE)size;
  else
    G->min = 0;

#ifdef HARD_SPHERE
    G->min = 1.;
#endif

#ifdef INTERACTION_Aziz
    Warning("Minimal distance for the spline is x=1 (limited by the hard-core of the potential)\n");
    G->min = 1.;
#endif

#ifdef INTERACTION_HYDROGEN
    Warning("Minimal distance for the spline is x=1 (limited by the hard-core of the potential)\n");
    G->min = 2;
#endif

#ifdef BC_ABSENT // size of the box is not a good number
  Message("  maximal distance for the spline Rpar12 = %lf\n", Rpar12);
  G->max = Rpar12;
  Warning("  changing L/2 = Rpar12 = %lf\n", Rpar12);
  Lhalf = Rpar12;
  L = 2.*Lhalf;
  Lhalf2 = Lhalf*Lhalf;
  Lwf = L;
  Lhalfwf = L/2;
  Lhalfwf2 = Lhalfwf*Lhalfwf;
#endif

  while(InteractionEnergy(G->min)>Epotcutoff) {
    if(G->min<1e-6) G->min = 1e-6;
    G->min *= 1.1;
    Warning("  Increasing G->min to %lf\n", G->min);
  }

  G->max2 = G->max*G->max;
  G->step = (G->max-G->min) / (DOUBLE) (G->size-1);
  for(i=0; i<size; i++) G->x[i] = G->min + (DOUBLE)(i)*G->step;

  Message("  Determining the scattering energy by shooting iterations\n");
#ifdef JASTROW_RIGHT_BOUNDARY_PHONONS
    flag_right_boundary_condition = 3; // solution, corresponding to phonons
    Warning("  Right boundary condition will be that of phonons\n");
#endif
  //xmin = 1e-5*E;
  if(flag_lnf_instead_of_f == OFF) { // solve for f(r)
    xmin = 1e-10;
    //xmin = 1e-6;
  } // solve for u(r)
  else{
    if(DIMENSION == 1)
      E = PI*PI/L/L; // characteristic energy
    else if(DIMENSION == 2)
      E = n; // characteristic energy
    else if(DIMENSION == 3)
      E = 1e-30; // characteristic energy

    xmin = 0.1*E;
  }
  ymin = SolveScatteringProblem(G, InteractionEnergy, xmin, flag_left_boundary_condition, flag_right_boundary_condition, &function_was_negative, flag_lnf_instead_of_f);

  if(isnan(ymin)) { // xmin is too small
    search = ON;
    Message("  Initial guess for xmin gives NaN, trying to adjust xmin...\n");
    while(search) {
      xmin /= 10;
      ymin = SolveScatteringProblem(G, InteractionEnergy, xmin, flag_left_boundary_condition, flag_right_boundary_condition, &function_was_negative, flag_lnf_instead_of_f);
      if(!isnan(ymin)) search = OFF;
    }
    if(iter++>1000)
      Error("Cannot find solution, adjust parameters of the shooting\n");
  }

  xmax = xmin;
  search = ON;
  while(search) {
    xmax *= energy_step;
    ymax = SolveScatteringProblem(G, InteractionEnergy, xmax, flag_left_boundary_condition, flag_right_boundary_condition, &function_was_negative, flag_lnf_instead_of_f);
    if(ymin*ymax<0) {
      xmin = xmax/energy_step;
      ymin = SolveScatteringProblem(G, InteractionEnergy, xmin, flag_left_boundary_condition, flag_right_boundary_condition, &function_was_negative, flag_lnf_instead_of_f);
      search = OFF;

    }
    if(iter++>10000)
      Error("Cannot find solution, adjust parameters of the shooting (1)\n");
  }

  if(G->f[2]<-0.1) {
    Warning("  function is still negative, trying negative energy solution...\n");
    function_was_negative = ON;
  }

  //Message(" xmin = %lf xmax = %lf ymin = %lf ymax = %lf\n", xmin, xmax, ymin, ymax);
  //ymin = SolveScatteringProblem(G, InteractionEnergy, xmin, flag_left_boundary_condition, flag_right_boundary_condition, &function_was_negative, flag_lnf_instead_of_f);
  //SaveGrid("ymin.dat", G);
  //ymax = SolveScatteringProblem(G, InteractionEnergy, xmax, flag_left_boundary_condition, flag_right_boundary_condition, &function_was_negative, flag_lnf_instead_of_f);
  //SaveGrid("ymax.dat", G);

#ifdef SPINFULL_TUNNELING
  Warning("  Looking only for negative energy solution!\n");
  function_was_negative = ON; //
#endif
  if(function_was_negative) {
    Warning("  Positive energy solution has nodes, trying to find negative energy solution\n");

#ifdef JASTROW_RIGHT_BOUNDARY_BOUND_STATE
    flag_right_boundary_condition = 2; // solution, corresponding to a bound state
    Warning("  Right boundary condition will be that of a bound state in free space\n");
#endif

#ifdef JASTROW_RIGHT_BOUNDARY_CONSTANT
    flag_right_boundary_condition = 0; // zero derivative
#endif

    // negative energy solution
    //xmin = -1000; // -100 @ n = 1e-2
    //energy_step = 0.99;
    xmin = -1e-50; // -100 @ n = 1e-2
    energy_step = 1.1;

    // loop again
    iter = 1;
    ymin = SolveScatteringProblem(G, InteractionEnergy, xmin, flag_left_boundary_condition, flag_right_boundary_condition, &function_was_negative, flag_lnf_instead_of_f);

    if(isnan(ymin)) { // xmin is too small
      search = ON;
      Message("  Initial guess for xmin gives NaN, trying to adjust xmin...\n");
      while(search) {
        xmin /= 10;
        ymin = SolveScatteringProblem(G, InteractionEnergy, xmin, flag_left_boundary_condition, flag_right_boundary_condition, &function_was_negative, flag_lnf_instead_of_f);
        if(!isnan(ymin) && function_was_negative == 0) search = OFF; // search for a function without nodes
      }
      if(iter++>1000)
        Error("Cannot find solution, adjust parameters of the shooting (2)\n");
    }
    // end of loop

    if(function_was_negative)
      Message("  negative energy solution will be constructed by shooting iterations...\n");
    else
      Warning("  negative energy solution still has nodes\n");
  }

  xmax = xmin;
  search = ON;
  while(search) {
    xmax *= energy_step;
    ymax = SolveScatteringProblem(G, InteractionEnergy, xmax, flag_left_boundary_condition, flag_right_boundary_condition, &function_was_negative, flag_lnf_instead_of_f);
    if(ymin*ymax<0) {
      xmin = xmax/energy_step;
      ymin = SolveScatteringProblem(G, InteractionEnergy, xmin, flag_left_boundary_condition, flag_right_boundary_condition, &function_was_negative, flag_lnf_instead_of_f);
      search = OFF;

    }
    if(iter++>10000)
      Error("Cannot find solution, adjust parameters of the shooting (3)\n");
  }

  search = ON;
  iter = 0;
  while(search) {
    // sanity check
    if(ymin*ymax>0)
      Message("  problems with shooting: xmin=%lf ymin=%lf, xmax=%lf ymax=%lf\n", xmin, ymin, xmax, ymax);

    x = (xmin+xmax)/2.;
    y = SolveScatteringProblem(G, InteractionEnergy, x, flag_left_boundary_condition, flag_right_boundary_condition, &function_was_negative, flag_lnf_instead_of_f);
    if(y*ymin < 0) {
      xmax = x;
      ymax = y;
    } else {
      xmin = x;
      ymin = y;
    }
    if(fabs(y)<precision) {
      x = (xmax*ymin-xmin*ymax) / (ymin-ymax);
      y = SolveScatteringProblem(G, InteractionEnergy, x, flag_left_boundary_condition, flag_right_boundary_condition, &function_was_negative, flag_lnf_instead_of_f);
      search = OFF;
    }
    if(iter++>Nitermax) {
      Warning("  maximal allowed number of iterations (%i) is exceeded. Have to stop.\n", Nitermax);
      search = OFF;
    }
    if(xmax == (xmin+xmax)/2 || xmin == (xmin+xmax)/2) {
      Warning("  xmax = xmin within numerical accuracy. Have to stop.\n");
      Message("  xmin = %.15e, ymin = %.15e\n", xmin, ymin);
      Message("  xmax = %.15e, ymax = %.15e\n", xmax, ymax);
      search = OFF;
    }
  }
  E = x;

  Message("  Scattering energy Escat= %.15" LE "\n", E);
  if(E>0)
    Message("  or k = sqrt(2E) = %.15" LE " = %.15" LE "(pi/L)\n", sqrt(2.*E), sqrt(2.*E) / (PI/L));
  else
    Message("  or k = sqrt(-2E) = %.15" LE " = %.15" LE "(pi/L)\n", sqrt(-2.*E), sqrt(-2.*E) / (PI/L));

  G->k = sqrt(fabs(E));

  Message("  Achieved precision is %" LE " \n", y);
#ifdef TRIAL_1D
  //flag_lnf_instead_of_f = OFF; // in 1D solve for ln f
  y = SolveScatteringProblem(G, InteractionEnergy, x, flag_left_boundary_condition, 0, &function_was_negative, flag_lnf_instead_of_f);
  Message("  Reiterate for ln f: achieved precision is %" LE " \n", y);
#endif

#ifdef JASTROW_DIMENSIONALITY_FROM_DPAR // restore true dimensionality
  Dpar = DIMENSION;
#endif

  Message("  Prepairing a spline... \n");
  G->size = grid_trial;
  G->min = G->x[0];
  G->max = G->x[grid_trial-1];
  G->step = (G->max - G->min) / (DOUBLE) (grid_trial-1);
  G->I_step = 1./G->step;
  G->max2 = G->max*G->max;

#ifdef INTERPOLATE_LOG
  Message("  Spline will fit: u(r), i.e. logarithm of f(r)\n");
  for(i=0; i<size; i++) {
    if(G->f[i]<0) {
      Warning(" spline element %i is negative, inverting\n", i);
      G->f[i] = fabs(G->f[i]);
    }
    G->f[i] = Log(G->f[i]);
  }

#ifdef JASTROW_GOES_TO_ONE_AT_LARGEST_DISTANCE
  if(G->f[size-1] != 0.) { // make sure that f(L/2) = 1
    Warning("  last element is different from 1. Adjusting the w.f. ...\n");
    for(i=0; i<size; i++) G->f[i] -= G->f[size-1];
  }
#endif
#else
  Message("  Spline will fit: f(r) and not ln f(r) \n");

#ifdef JASTROW_GOES_TO_ONE_AT_LARGEST_DISTANCE
  if(G->f[size-1] != 1.) { // make sure that f(L/2) = 1
    Warning("  last element is different from 1. Adjusting the w.f. ...\n");
    for(i=0; i<size; i++) G->f[i] /= G->f[size-1];
  }
#else
  Warning("  Jastow term is allowed to be different from one at the largest distance\n");
#endif

#endif

  // Check spline for NaN
  for(i=0; i<size; i++) {
    if(isnan(G->f[1])) {
      Warning("Spline error, NaN found at position %i\n", i);
    }
  }

  Message("done\n");
}

/************************** Construct Grid Spline 2 body solution  **************************/
void ConstructGridSpline2Body_ij(struct Grid *G, int spin1, int spin2) {
  // choose boundary condition of the scattering problem
  // 0: zero boundary condition f(0) = 0;
  // 1: f'(0) = 0; 
  int flag_left_boundary_condition = 0; // 0 or 1
  int flag_right_boundary_condition = 0; // 0 for f'(R) = 0 and 1 for b.s.
  int size = 1000;
  DOUBLE E=0.; // scattering energy
  DOUBLE xmin, xmax, ymin, ymax, x, y;
  int search = ON;
  DOUBLE precision = 1e-14; //1e-8
  int Nitermax = 1000000; // maximal allowed number of iterations
  DOUBLE energy_step = 1.01; // to find ymin*ymax<0
  int iter=0;
  int i;
  int function_was_negative = OFF;
  int flag_lnf_instead_of_f = OFF; // ON - has some problems in 1D solve Scroedinger equation for ln f instead of f
  int try_first_positive_energy_solution = ON;

#ifdef SPINFULL
  Message("\n  Jastrow (%i - %i) term will be constructed as 2-body scattering solution with for interaction potential\n", spin1, spin2);
#else
  Message("\n  Jastrow term will be constructed as 2-body scattering solution with for interaction potential\n");
#endif
  Message("  Spline size grid_trial=%i\n", grid_trial);
  size = grid_trial;
  AllocateWFGrid(G, size);
  G->size = size;

//#ifdef INTERPOLATE_LOG // unstable
//  flag_lnf_instead_of_f = ON;
//  Message("  Schroedinger equation will be solved for logarithm u(r) = ln f(r)\n");
//#endif

  if(Epotcutoff<0) {// 1e3 maximal allowed potential energy for V(xmin)
    Epotcutoff = 1e3;
  }
  Message("  maximal allowed potential energy for V(xmin) = Epotcutoff = %e\n", Epotcutoff);

#ifdef SAME_SPIN_CONSTANT_JASTROW
  if(spin1 == spin2) {
    Warning("  spin1 = spin2, no interaction, 2body problem will not be solved, exiting\n");
    return;
  }
#endif

  // choose boundary condition of the scattering problem in main.h
#ifdef JASTROW_LEFT_BOUNDARY_ZERO
  // 0: zero boundary condition f(0) = 0;
  if(spin1 == spin2) flag_left_boundary_condition = 0;
  //flag_left_boundary_condition = 0;
  Message("  zero boundary condition f(0) = 0 will be used\n");
#else
  // 1: f'(0) = 0;
  flag_left_boundary_condition = 1;
  Message("  zero derivative condition f'(0) = 0 will be used\n");
#endif

#ifdef JASTROW_RIGHT_BOUNDARY_SPINFULL_DIFFERENT_SPIN_BOUND_STATE
  if(spin1!=spin2) {
    try_first_positive_energy_solution = OFF;
    flag_right_boundary_condition = 2; // solution, corresponding to a bound state
  }
#endif

#ifdef JASTROW_RIGHT_BOUNDARY_SPINFULL_SAME_SPIN_BOUND_STATE
  if(spin1==spin2) {
    try_first_positive_energy_solution = OFF;
    flag_right_boundary_condition = 2; // solution, corresponding to a bound state
  }
#endif

#ifdef JASTROW_RIGHT_BOUNDARY_PHONONS
  flag_right_boundary_condition = 3; // solution, corresponding to phonons
  Warning("  Right boundary condition will be that of phonons\n");
#endif

#ifdef JASTROW_RIGHT_BOUNDARY_CONSTANT
  flag_right_boundary_condition = 0; // zero derivative
#endif

#ifdef BC_ABSENT // size of the box is not a good number
  if(spin1 == spin2) {
    //flag_left_boundary_condition = 0;
    //Message("  zero boundary condition f(0) = 0 will be used for spin1 = spin2\n");
    G->max = Rpar;
    Message("  max spline distance Rpar = %f\n", G->max);
  }
  else {
    //flag_left_boundary_condition = 1;
    //Message("  zero derivative condition f'(0) = 0 will be used for spin1 different from spin2\n");
    G->max = Rpar12;
    Message("  max spline distance Rpar12 = %f\n", G->max);
  }
  if(Rpar>Rpar12)
    Lhalf = Rpar;
  else
    Lhalf = Rpar12;

  L = 2.*Lhalf;
  Lhalf2 = Lhalf*Lhalf;
  Lwf = L;
  Lhalfwf = L/2;
  Lhalfwf2 = Lhalfwf*Lhalfwf;
#else
  if(spin1 == spin2) {
    if(Rpar<1e-8) Rpar = 1; // if zero
    G->max = Rpar*Lhalf;
    Message("  max spline distance Rpar = %f (converted from L/2)\n", G->max);
  }
  else {
    if(Rpar <1e-8) Rpar12 = 1;
    G->max = Rpar12*Lhalf;
    Message("  max spline distance Rpar12 = %f\n", G->max);
  }
  if(Rpar>Rpar12) // chose the largest Lhalf
    Lhalfwf = Rpar*Lhalf;
  else
    Lhalfwf = Rpar12*Lhalf;

  Lwf = 2.*Lhalfwf;
  //Lhalf2 = Lhalf*Lhalf;
  //Lwf = L;
  //Lhalfwf = L/2;
  Lhalfwf2 = Lhalfwf*Lhalfwf;
#endif

  /*if (spin1 == spin2) {
    flag_left_boundary_condition = 0;
    Message("  zero boundary condition f(0) = 0 will be used for spin1 = spin2\n");
  }
  else {
    flag_left_boundary_condition = 1;
    Message("  zero derivative condition f'(0) = 0 will be used for spin1 different from spin2\n");
  }*/

  if(G->max<1e-8) Error("  Cannot solve with G->max = 0\n");

  if(DIMENSION>1)
    G->min = G->max/(DOUBLE)size;
  else
    G->min = 0;

#ifdef HARD_SPHERE
    G->min = 1.;
#endif

#ifdef INTERACTION_Aziz
    Warning("Minimal distance for the spline is x=1 (limited by the hard-core of the potential)\n");
    G->min = 1.;
#endif

  while(InteractionEnergy_ij(G->min, spin1, spin2)>Epotcutoff) {
    G->min *= 1.1;
    Warning("  Increasing G->min to %lf\n", G->min);
  }

  G->max2 = G->max*G->max;
  G->step = (G->max-G->min) / (DOUBLE) (G->size-1);
  for(i=0; i<size; i++) G->x[i] = G->min + (DOUBLE)(i)*G->step;

  if(try_first_positive_energy_solution) {
    Message("  Determining the scattering energy by shooting iterations\n");
    //xmin = 1e-5*E;
    if(flag_lnf_instead_of_f == OFF) { // solve for f(r)
      //xmin = 1e-50;
      xmin = 1e-6;
      energy_step = 1.1;
    } // solve for u(r)
    else{
      if(DIMENSION == 1)
        E = PI*PI/L/L; // characteristic energy
      else if(DIMENSION == 2)
        E = n; // characteristic energy
      else if(DIMENSION == 3)
        E = 1e-30; // characteristic energy

      xmin = 0.1*E;
      energy_step = 1.1;
    }
    ymin = SolveScatteringProblem_ij(G, spin1, spin2, xmin, flag_left_boundary_condition, flag_right_boundary_condition, &function_was_negative, flag_lnf_instead_of_f);

    if(isnan(ymin)) { // xmin is too small
      search = ON;
      Message("  Initial guess for xmin gives NaN, trying to adjust xmin...\n");
      while(search) {
        xmin /= 10;
        ymin = SolveScatteringProblem_ij(G, spin1, spin2, xmin, flag_left_boundary_condition, flag_right_boundary_condition, &function_was_negative, flag_lnf_instead_of_f);
        if(!isnan(ymin)) search = OFF;
      }
      if(iter++>1000)
        Error("Cannot find solution, adjust parameters of the shooting (4)\n");
    }

    xmax = xmin;
    search = ON;
    while(search) {
      xmax *= energy_step;
      ymax = SolveScatteringProblem_ij(G, spin1, spin2, xmax, flag_left_boundary_condition, flag_right_boundary_condition, &function_was_negative, flag_lnf_instead_of_f);
      if(ymin*ymax<0) {
        xmin = xmax/energy_step;
        ymin = SolveScatteringProblem_ij(G, spin1, spin2, xmin, flag_left_boundary_condition, flag_right_boundary_condition, &function_was_negative, flag_lnf_instead_of_f);
        search = OFF;
      }
      if(iter++>10000)
        Error("Cannot find solution, adjust parameters of the shooting\n");
    }

    if(G->f[2]<-0.1) {
      Warning("  function is still negative, trying negative energy solution...\n");
      function_was_negative = ON;
    }
  }

#ifdef SPINFULL_TUNNELING
  Warning("  Looking only for negative energy solution!\n");
  function_was_negative = ON; //
#endif

  if(function_was_negative || try_first_positive_energy_solution == OFF) {
    Warning("  Positive energy solution has nodes, trying to find negative energy solution\n");
    // negative energy solution
    //xmin = -1000; // -100 @ n = 1e-2
    //energy_step = 0.99;
    xmin = -1e-50; // -100 @ n = 1e-2
    energy_step = 1.1;

    // loop again
    iter = 1;
    ymin = SolveScatteringProblem_ij(G, spin1, spin2, xmin, flag_left_boundary_condition, flag_right_boundary_condition, &function_was_negative, flag_lnf_instead_of_f);

    if(isnan(ymin)) { // xmin is too small
      search = ON;
      Message("  Initial guess for xmin gives NaN, trying to adjust xmin...\n");
      while(search) {
        xmin /= 10;
        ymin = SolveScatteringProblem_ij(G, spin1, spin2, xmin, flag_left_boundary_condition, flag_right_boundary_condition, &function_was_negative, flag_lnf_instead_of_f);
        if(!isnan(ymin) && function_was_negative == 0) 
          search = OFF; // search for a function without nodes
      }
      if(iter++>1000)
        Error("Cannot find solution, adjust parameters of the shooting\n");
    }
    // end of loop

#ifdef INTERPOLATE_LOG // if not u = ln f solution, report sign
    Message("  searching solution for the logarithm u = ln f\n");
#else
    if(function_was_negative)
      Message("  negative energy solution will be constructed by shooting iterations...\n");
    else
      Warning("  negative energy solution still has nodes\n");
#endif
  }

  xmax = xmin;
  search = ON;
  while(search) {
    xmax *= energy_step;
    ymax = SolveScatteringProblem_ij(G, spin1, spin2, xmax, flag_left_boundary_condition, flag_right_boundary_condition, &function_was_negative, flag_lnf_instead_of_f);
    if(ymin*ymax<0) {
      xmin = xmax/energy_step;
      ymin = SolveScatteringProblem_ij(G, spin1, spin2, xmin, flag_left_boundary_condition, flag_right_boundary_condition, &function_was_negative, flag_lnf_instead_of_f);
      search = OFF;
    }
    if(iter++>10000)
      Error("Cannot find solution, adjust parameters of the shooting\n");
  }

  search = ON;
  iter = 0;
  while(search) {
    // sanity check
    if(ymin*ymax>0)
      Message("  problems with shooting: xmin=%lf ymin=%lf, xmax=%lf ymax=%lf\n", xmin, ymin, xmax, ymax);

    x = (xmin+xmax)/2.;
    y = SolveScatteringProblem_ij(G, spin1, spin2, x, flag_left_boundary_condition, flag_right_boundary_condition, &function_was_negative, flag_lnf_instead_of_f);
    if(y*ymin < 0) {
      xmax = x;
      ymax = y;
    } else {
      xmin = x;
      ymin = y;
    }
    if(fabs(y)<precision) {
      x = (xmax*ymin-xmin*ymax) / (ymin-ymax);
      y = SolveScatteringProblem_ij(G, spin1, spin2, x, flag_left_boundary_condition, flag_right_boundary_condition, &function_was_negative, flag_lnf_instead_of_f);
      search = OFF;
    }
    if(iter++>Nitermax) {
      Warning("  maximal allowed number of iterations (%i) is exceeded. Have to stop.\n", Nitermax);
      search = OFF;
    }
    if(xmax == (xmin+xmax)/2 || xmin == (xmin+xmax)/2) {
      Warning("  xmax = xmin within numerical accuracy. Have to stop.\n");
      Message("  xmin = %e, ymin = %e\n", xmin, ymin);
      Message("  xmax = %e, ymax = %e\n", xmax, ymax);
      search = OFF;
    }
  }
  E = x;

  Message("  Scattering energy %.15" LE "\n", E);
  if(E>0)
    Message("  or k = sqrt(2E) = %.15" LE " = %.15" LE "(pi/L)\n", sqrt(2.*E), sqrt(2.*E) / (PI/L));
  else
    Message("  or k = sqrt(-2E) = %.15" LE " = %.15" LE "(pi/L)\n", sqrt(-2.*E), sqrt(-2.*E) / (PI/L));

  G->k = sqrt(fabs(E));

  Message("  Achieved precision is %" LE " \n", y);
#ifdef TRIAL_1D
  //flag_lnf_instead_of_f = OFF; // in 1D solve for ln f
  y = SolveScatteringProblem_ij(G, spin1, spin2, x, flag_left_boundary_condition, flag_right_boundary_condition, &function_was_negative, flag_lnf_instead_of_f);
  Message("  Reiterate for ln f: achieved precision is %" LE " \n", y);
#endif

  Message("  Prepairing a grid for interporation... \n");
  G->size = grid_trial;
  G->min = G->x[0];
  G->max = G->x[grid_trial-1];
  G->step = (G->max - G->min) / (DOUBLE) (grid_trial-1);
  G->I_step = 1./G->step;
  G->max2 = G->max*G->max;

#ifdef INTERPOLATE_LOG
  Message("  Spline will fit: u(r), i.e. logarithm of f(r)\n");
  for(i=0; i<size; i++) {
    if(G->f[i]<0) {
      Warning(" spline element %i is negative, inverting\n", i);
      G->f[i] = fabs(G->f[i]);
    }
    G->f[i] = Log(G->f[i]);
  }

#ifdef JASTROW_GOES_TO_ONE_AT_LARGEST_DISTANCE
  if(G->f[size-1] != 0.) { // make sure that f(L/2) = 1
    Warning("  last element is different from 1. Adjusting the w.f. ...\n");
    for(i=0; i<size; i++) G->f[i] -= G->f[size-1];
  }
#endif

#else
  Message("  Spline will fit: f(r) and not ln f(r) \n");

#ifdef JASTROW_GOES_TO_ONE_AT_LARGEST_DISTANCE
  if(G->f[size-1] != 1.) { // make sure that f(L/2) = 1
    Warning("  last element is different from 1. Adjusting the w.f. ...\n");
    for(i=0; i<size; i++) G->f[i] /= G->f[size-1];
  }
#else
  Warning("  Jastow term is allowed to be different from one at the largest distance\n");
#endif
#endif

  // Check spline for NaN
  for(i=0; i<size; i++) {
    if(isnan(G->f[1])) {
      Warning("Grid error, NaN found at position %i\n", i);
    }
  }

  Message("done\n");
}

/************************** Construct Grid Spline 2 body solution  **************************/
void ConstructGridSpline2BodyPhonons(struct Grid *G) {
  // choose boundary condition of the scattering problem
  // 0: zero boundary condition f(0) = 0;
  // 1: f'(0) = 0;
  int flag_left_boundary_condition = 0; // 0 or 1
  int size = 1000;
  DOUBLE E=0.; // scattering energy
  DOUBLE xmin, xmax, ymin, ymax, x, y;
  int search = ON;
  DOUBLE precision = 1e-14; //1e-8
  int Nitermax = 1000000; // maximal allowed number of iterations
  DOUBLE energy_step = 1.01; // to find ymin*ymax<0
  DOUBLE Epotcutoff = 1e3; // maximal allowed potential energy for V(xmin)
  int iter=0;
  int i;
  int function_was_negative = OFF;
  int flag_lnf_instead_of_f = OFF; // ON - has some problems in 1D solve Scroedinger equation for ln f instead of f
  int size_total_with_phonons; // Grid->size will be used for shooting

  Message("  Jastrow term will be constructed as 2-body scattering solution with for interaction potential +  phonons\n");
  Message("  Spline size grid_trial=%i\n", grid_trial);

  size = grid_trial;
  size_total_with_phonons = size;
  AllocateWFGrid(G, size);

  if(DIMENSION>1)
    G->min = G->max/(DOUBLE)size;
  else
    G->min = 0;

#ifdef HARD_SPHERE
    G->min = 1.;
#endif

#ifdef INTERACTION_Aziz
    Warning("Minimal distance for the spline is x=1 (limited by the hard-core of the potential)\n");
    G->min = 1.;
#endif

  while(InteractionEnergy(G->min)>Epotcutoff) {
    G->min *= 1.1;
    Warning("  Increasing G->min to %lf\n", G->min);
  }
  Message("  Matching distance Rpar12 = %lf\n", Rpar12);
  Message("  Exponent factor Apar12 = %lf\n", Apar12);

  G->max = Lhalf;
  G->max2 = G->max*G->max;
  G->step = (G->max-G->min) / (DOUBLE) (G->size-1);
  for(i=0; i<size; i++) G->x[i] = G->min + (DOUBLE)(i)*G->step;

  while(InteractionEnergy(G->min)>Epotcutoff) {
    G->min *= 1.1;
    Warning("  Increasing G->min to %lf\n", G->min);
  }

  Message("  Initializing with phonons...\n");
#ifdef TRIAL_1D
  Message("  f(r) = |sin(PI x/L)|^Apar12, i.e. Apar = 1/KL where KL is the Luttinger parameter\n");
#endif
#ifdef TRIAL_2D
  Message("  f(r) = exp(-Apar12/x - Apar12/(L-x)) / exp(-2.*Apar12/Lhalf)\n");
#endif
#ifdef TRIAL_3D
  Message("  f(r) = exp(-Apar12/x/x - Apar12/(L-x)/(L-x)) / exp(-2.*Apar12/Lhalf/Lhalf)\n");
#endif

  i = size_total_with_phonons-1;
  while(G->x[i]>Rpar12) {
    x = G->x[i];
#ifdef TRIAL_1D
    G->f[i] = pow(sin(PI*x/L), Apar12);
#endif
#ifdef TRIAL_2D
    G->f[i] = exp(-Apar12/x - Apar12/(L-x)) / exp(-2.*Apar12/Lhalf);
#endif
#ifdef TRIAL_3D
    G->f[i] = exp(-Apar12/x/x - Apar12/(L-x)/(L-x)) / exp(-2.*Apar12/Lhalf/Lhalf);
#endif
    i--;
  }

  G->size = i+3; // now elements [size-1] and [size-2] are initialized
  Message("  Phonon regime, index %i < i < %i\n", G->size, size_total_with_phonons);
  Message("  f(R) = %lf\n", G->f[G->size-1]);
  Message("  f'(R) = %lf\n", (G->f[G->size-1]-G->f[G->size-2])/G->step);

  Message("  Determining the scattering energy by shooting iterations\n");
  //xmin = 1e-5*E;
  if(flag_lnf_instead_of_f == OFF) { // solve for f(r)
    //xmin = 1e-50;
    xmin = 1e-6;
    energy_step = 1.1;
  } // solve for u(r)
  else{
    if(DIMENSION == 1)
      E = PI*PI/L/L; // characteristic energy
    else if(DIMENSION == 2)
      E = n; // characteristic energy
    else if(DIMENSION == 3)
      E = 1e-30; // characteristic energy

    xmin = 0.1*E;
    energy_step = 1.1;
  }
  ymin = SolveScatteringProblem(G, InteractionEnergy, xmin, flag_left_boundary_condition, 1, &function_was_negative, flag_lnf_instead_of_f);

  if(isnan(ymin)) { // xmin is too small
    search = ON;
    Message("  Initial guess for xmin gives NaN, trying to adjust xmin...\n");
    while(search) {
      xmin /= 10;
      ymin = SolveScatteringProblem(G, InteractionEnergy, xmin, flag_left_boundary_condition, 1, &function_was_negative, flag_lnf_instead_of_f);
      if(!isnan(ymin)) search = OFF;
    }
    if(iter++>1000)
      Error("Cannot find solution, adjust parameters of the shooting\n");
  }

  xmax = xmin;
  search = ON;
  while(search) {
    xmax *= energy_step;
    ymax = SolveScatteringProblem(G, InteractionEnergy, xmax, flag_left_boundary_condition, 1, &function_was_negative, flag_lnf_instead_of_f);
    if(ymin*ymax<0) {
      xmin = xmax/energy_step;
      ymin = SolveScatteringProblem(G, InteractionEnergy, xmin, flag_left_boundary_condition, 1, &function_was_negative, flag_lnf_instead_of_f);
      search = OFF;

    }
    if(iter++>10000)
      Error("Cannot find solution, adjust parameters of the shooting\n");
  }

  if(G->f[2]<-0.1) {
    Warning(" function is still negative, trying negative energy solution...\n");
    function_was_negative = ON;
  }

#ifdef SPINFULL_TUNNELING
  Warning("  Looking only for negative energy solution!\n");
  function_was_negative = ON; //
#endif
  if(function_was_negative) {
    Warning("  Positive energy solution has nodes, trying to find negative energy solution\n");
    // negative energy solution
    //xmin = -1000; // -100 @ n = 1e-2
    //energy_step = 0.99;
    xmin = -1e-50; // -100 @ n = 1e-2
    energy_step = 1.1;

    // loop again
    iter = 1;
    ymin = SolveScatteringProblem(G, InteractionEnergy, xmin, flag_left_boundary_condition, 1, &function_was_negative, flag_lnf_instead_of_f);

    if(isnan(ymin)) { // xmin is too small
      search = ON;
      Message("  Initial guess for xmin gives NaN, trying to adjust xmin...\n");
      while(search) {
        xmin /= 10;
        ymin = SolveScatteringProblem(G, InteractionEnergy, xmin, flag_left_boundary_condition, 1, &function_was_negative, flag_lnf_instead_of_f);
        if(!isnan(ymin) && function_was_negative == 0) search = OFF; // search for a function without nodes
      }
      if(iter++>1000)
        Error("Cannot find solution, adjust parameters of the shooting\n");
    }
    // end of loop

    if(function_was_negative)
      Message("  negative energy solution will be constructed by shooting iterations...\n");
    else
      Warning("  negative energy solution still has nodes\n");
  }

  xmax = xmin;
  search = ON;
  while(search) {
    xmax *= energy_step;
    ymax = SolveScatteringProblem(G, InteractionEnergy, xmax, flag_left_boundary_condition, 1, &function_was_negative, flag_lnf_instead_of_f);
    if(ymin*ymax<0) {
      xmin = xmax/energy_step;
      ymin = SolveScatteringProblem(G, InteractionEnergy, xmin, flag_left_boundary_condition, 1, &function_was_negative, flag_lnf_instead_of_f);
      search = OFF;

    }
    if(iter++>10000)
      Error("Cannot find solution, adjust parameters of the shooting\n");
  }

  search = ON;
  iter = 0;
  while(search) {
    // sanity check
    if(ymin*ymax>0)
      Message("  problems with shooting: xmin=%lf ymin=%lf, xmax=%lf ymax=%lf\n", xmin, ymin, xmax, ymax);

    x = (xmin+xmax)/2.;
    y = SolveScatteringProblem(G, InteractionEnergy, x, flag_left_boundary_condition, 1, &function_was_negative, flag_lnf_instead_of_f);
    if(y*ymin < 0) {
      xmax = x;
      ymax = y;
    } else {
      xmin = x;
      ymin = y;
    }
    if(fabs(y)<precision) {
      x = (xmax*ymin-xmin*ymax) / (ymin-ymax);
      y = SolveScatteringProblem(G, InteractionEnergy, x, flag_left_boundary_condition, 1, &function_was_negative, flag_lnf_instead_of_f);
      search = OFF;
    }
    if(iter++>Nitermax) {
      Warning("  maximal allowed number of iterations (%i) is exceeded. Have to stop.\n", Nitermax);
      search = OFF;
    }
    if(xmax == (xmin+xmax)/2 || xmin == (xmin+xmax)/2) {
      Warning("  xmax = xmin within numerical accuracy. Have to stop.\n");
      Message("  xmin = %e, ymin = %e\n", xmin, ymin);
      Message("  xmax = %e, ymax = %e\n", xmax, ymax);
      search = OFF;
    }
  }
  E = x;

  Message("  Scattering energy %.15" LE "\n", E);
  if(E>0)
    Message("  or k = sqrt(2E) = %.15" LE " in units of pi/L\n", sqrt(2.*E) / (PI/L));
  else
    Message("  or k = sqrt(-2E) = %.15" LE " in units of pi/L\n", sqrt(-2.*E) / (PI/L));

  Message("  Achieved precision is %" LE " \n", y);
#ifdef TRIAL_1D
  //flag_lnf_instead_of_f = OFF; // in 1D solve for ln f
  y = SolveScatteringProblem(G, InteractionEnergy, x, flag_left_boundary_condition, 1, &function_was_negative, flag_lnf_instead_of_f);
  Message("  Reiterate for ln f: achieved precision is %" LE " \n", y);
#endif

  Message("  Prepairing a spline... \n");
  G->size = grid_trial;
  G->min = G->x[0];
  G->max = G->x[grid_trial-1];
  G->step = (G->max - G->min) / (DOUBLE) (grid_trial-1);
  G->I_step = 1./G->step;
  G->max2 = G->max*G->max;

#ifdef INTERPOLATE_LOG
  Message("  Spline will fit: u(r), i.e. logarithm of f(r)\n");
  for(i=0; i<size; i++) {
    if(G->f[i]<0) {
      Warning(" spline element %i is negative, inverting\n", i);
      G->f[i] = fabs(G->f[i]);
    }
    G->f[i] = Log(G->f[i]);
  }

  if(G->f[size-1] != 0.) { // make sure that f(L/2) = 1
    Warning("  last element is different from 1. Adjusting the w.f. ...\n");
    for(i=0; i<size; i++) G->f[i] -= G->f[size-1];
  }
#else
  Message("  Spline will fit: f(r) and not ln f(r) \n");

  if(G->f[size-1] != 1.) { // make sure that f(L/2) = 1
    Warning("  last element is different from 1. Adjusting the w.f. ...\n");
    for(i=0; i<size; i++) G->f[i] /= G->f[size-1];
  }
#endif

  // Check spline for NaN
  for(i=0; i<size; i++) {
    if(isnan(G->f[1])) {
      Warning("Spline error, NaN found at position %i\n", i);
    }
  }

  Message("done\n");
}

/************************** Construct Grid Spline 2 body solution  **************************/
// see also SolveScatteringProblem_ij
// flag_boundary_condition_left:
// 0 - zero boundary condition, f(0) = 0
// 1 - zero derivatice, f'(0) = 0
//
// flag_boundary_condition_right:
// 0 - set f(R) = 1; f'(R) = 0
// 1 - f[size-2] and f[size-1] are already initialized
// 2 - bound-state in free space
//     1D f(x) = 
//     2D f(x) = K_1(sqrt(-E) r); f'/f = -k K_1(sqrt(-E) r) / K_0(sqrt(-E) r)
//     3D f(x) = exp(sqrt(-E) r) / r;   f'/f = -k - 1/r
DOUBLE SolveScatteringProblem(struct Grid *G, DOUBLE(*InteractionEnergy)(DOUBLE x), DOUBLE E, int flag_left_boundary_condition, int flag_right_boundary_condition, int *function_was_negative, int flag_lnf_instead_of_f) {
  DOUBLE r,dr,dr2;
  int i;
  //int flag_lnf_instead_of_f = ON; // 1D : OFF - f, ON - u = ln f
  static int first_time = ON;
  DOUBLE det;
  DOUBLE invert_last_point_factor = 1.; // 1 or -1
  DOUBLE Vint; // save read
  DOUBLE D = DIMENSION; // effective dimensionality
#ifdef JASTROW_DIMENSIONALITY_FROM_DPAR
  D = Dpar;
#endif

  *function_was_negative = OFF;

#ifdef TRIAL_1D
  if(first_time == ON) 
    if(flag_lnf_instead_of_f == ON)
      Message("  Scattering problem is solved for log of two-body w.f. u(r),\n  i.e. -u'' - (D-1)/r u' + [V(r)-E] u  - (u')^2 = 0\n");
    else
      Message("  Scattering problem is solved for the two-body w.f. f(r),\n  i.e. -f'' - (D-1)/r f' + [V(r)-E] f = 0\n");
  first_time = OFF;
#endif

  if(flag_lnf_instead_of_f == OFF) { // Not logarithmic (i.e. linear scales): standard Schroedinger equation
    // -f'' - (D-1)/r f' + [V(r)-E] f = 0
    // - (f[i+1]+f[i-1]-2 f[i]) / dr^2 - (f[i+2]-f[i])/(2dr) + (V[i]-E) f[i] = 0
    // f'(L/2) = 0
    // f[N] = f[N-2]
    dr = G->step;
    dr2 = dr*dr;
    i = G->size-1;
    r = G->x[i];
    if(flag_right_boundary_condition == 0) { // zero derivative
      G->f[i] = 1.;
      Vint = InteractionEnergy(r);
      G->f[i-1] = G->f[i]*(1.+0.5*(Vint-E)*dr2);
    } else if(flag_right_boundary_condition == 2) { // solution, corresponding to a bound state
#ifdef TRIAL_2D
      G->f[i] = K0(-fabs(sqrt(fabs(E)))*G->x[i]);
      G->f[i-1] = K0(-fabs(sqrt(fabs(E)))*G->x[i-1]);
#endif
    }
    else if(flag_right_boundary_condition == 3) { // solution, corresponding to phonons
#ifdef JASTROW_RIGHT_BOUNDARY_PHONONS // initialize accoring to the phononic wave function
#ifdef TRIAL_3D
      G->f[i] = exp(- Apar/(G->x[i]*G->x[i]) - Apar/((L-G->x[i])*(L-G->x[i])) + 2.*Apar/Lhalf2);
      G->f[i-1] = exp(- Apar/(G->x[i-1]*G->x[i-1]) - Apar/((L-G->x[i-1])*(L-G->x[i-1])) + 2.*Apar/Lhalf2);
#else
      G->f[i] = exp(- Apar/G->x[i] - Apar/(L-G->x[i]) + 2.*Apar/Lhalf);
      G->f[i-1] = exp(- Apar/G->x[i-1] - Apar/(L-G->x[i-1]) + 2.*Apar/Lhalf);
#endif
#endif
    }
    // otherwize f[i] and f[i-1] are already defined

    for(i=G->size-2; i>=1; i--) {
      r = G->x[i];
      Vint = InteractionEnergy(r);
      G->f[i-1] = (G->f[i]*(2.+(Vint-E)*dr2)*r + G->f[i+1]*(-r-(D-1.)*dr*0.5))/(r-(D-1.)*dr*0.5);
      if(isinf(G->f[i-1]))
        Warning("  inf encountered at point %i at distance %lf, Vint(r) = %lf\n", i, r, Vint);
      if(isnan(G->f[i-1]))
        Warning("  NaN encountered at point %i at distance %lf, Vint(r) = %lf\n", i, r, Vint);
      if(G->f[i-1]<0)
        *function_was_negative = ON;
      //G->f[i-1] = (G->f[i]*(2.*r*r+(InteractionEnergy(r)*r*r-E*r*r)*dr2) + G->f[i+1]*r*(-r-(D-1.)*dr*0.5))/(r*r-(D-1.)*r*dr*0.5);
    }
  }
  else { // equation for u(r) = ln f(r)
    // -u'' - (D-1)/r u' + [V(r)-E] u  - (u')^2 = 0
    // - (u[i+1]-2u[i]+u[i-1])/dt2 - (D-2)/r (u[i+1]-u[i-1])/(2dt) + [V(i)-E] u[i] - [(u[i+1]-u[i-1]) / (2 dt)]^2 = 0
    // u'(L/2) = 0
    // u[N] = u[N-2]
    dr = G->step;
    dr2 = dr*dr;
    i = G->size-1;
    r = G->x[i];

    if(flag_right_boundary_condition == 0) {
      // zero derivative for f' = exp(u)' = f u'
      //G->f[i] = 0;
      //G->f[i-1] = exp(G->f[i])*(1.+0.5*(InteractionEnergy(r)-E)*dr2);
      //G->f[i-1] = log(G->f[i-1]);
      // zero derivative for u
      G->f[i] = 1; // 0 gives a trivial and non-normalizable solution u = 0
      G->f[i-1] = G->f[i]*(1.+0.5*(InteractionEnergy(r)-E)*dr2);
    }
    else if(flag_right_boundary_condition == 2) { // solution, corresponding to a bound state
#ifdef TRIAL_2D
      G->f[i] = log(K0(-fabs(sqrt(fabs(E)))*G->x[i]));
      G->f[i-1] = log(K0(-fabs(sqrt(fabs(E)))*G->x[i-1]));
#endif
    }
    else if(flag_right_boundary_condition == 3) { // solution, corresponding to phonons
      G->f[i] = - Apar/(G->x[i]*G->x[i]) - Apar/((L-G->x[i])*(L-G->x[i])) + 2.*Apar/Lhalf2;
      G->f[i-1] = - Apar/(G->x[i-1]*G->x[i-1]) - Apar/((L-G->x[i-1])*(L-G->x[i-1])) + 2.*Apar/Lhalf2;

#ifdef JASTROW_RIGHT_BOUNDARY_PHONONS // initialize accoring to the phononic wave function
#ifdef TRIAL_3D
      G->f[i] = - Apar/(G->x[i]*G->x[i]) - Apar/((L-G->x[i])*(L-G->x[i])) + 2.*Apar/Lhalf2;
      G->f[i-1] = - Apar/(G->x[i-1]*G->x[i-1]) - Apar/((L-G->x[i-1])*(L-G->x[i-1])) + 2.*Apar/Lhalf2;
#else
      G->f[i] = - Apar/G->x[i] - Apar/(L-G->x[i]) + 2.*Apar/Lhalf;
      G->f[i-1] = - Apar/G->x[i-1] - Apar/(L-G->x[i-1]) + 2.*Apar/Lhalf;
#endif
#endif

    } // otherwize f[i] and f[i-1] are already defined

    for(i=G->size-2; i>=1; i--) {
      r = G->x[i];
      det = 1.-2.*(G->f[i+1]-G->f[i]) + G->f[i]*dr*dr*(InteractionEnergy(r)-E);
      if(det<0) {
        Warning("negative Jastrow, i=%i\n", i);
        det = - det;
        *function_was_negative = ON;
        invert_last_point_factor = -1.;
        return -r;
      }
      G->f[i-1] = G->f[i+1] + 2.*(-1. + sqrt(det));
    }
    for(i=0; i<G->size; i++) G->f[i] -= G->f[G->size-1]; // constant shift to ensure that f(L/2) = 1

    for(i=0; i<G->size; i++) G->f[i] = exp(G->f[i]); // convert back to standard notation
  }

  if(flag_left_boundary_condition == 0) // zero boundary condition
    return G->f[0]*invert_last_point_factor;
  else if(flag_left_boundary_condition == 1) // zero derivative
    return (G->f[1] - G->f[0])*invert_last_point_factor;
  else // f'(0) / f(0) = - 1/a, i.e. a f'(0)/f(0) + 1 = 0
    return a*(G->f[1] - G->f[0])/(G->step*0.5*(G->f[1] + G->f[0])) + 1.;
}

/************************** Construct Grid Spline 2 body solution  **************************/
// flag_boundary_condition_left:
// 0 - zero boundary condition, f(0) = 0
// 1 - zero derivatice, f'(0) = 0
//
// flag_boundary_condition_right:
// 0 - set f(R) = 1; f'(R) = 0
// 1 - f[size-2] and f[size-1] are already initialized
// 2 - bound-state in free space
//     1D f(x) = 
//     2D f(x) = K_1(sqrt(-E) r); f'/f = -k K_1(sqrt(-E) r) / K_0(sqrt(-E) r)
//     3D f(x) = exp(sqrt(-E) r) / r;   f'/f = -k - 1/r
DOUBLE SolveScatteringProblem_ij(struct Grid *G, int spin1, int spin2, DOUBLE E, int flag_left_boundary_condition, int flag_right_boundary_condition, int *function_was_negative, int flag_lnf_instead_of_f) {
  DOUBLE r,dr,dr2;
  int i;
  //int flag_lnf_instead_of_f = ON; // 1D : OFF - f, ON - u = ln f
  static int first_time = OFF;
  DOUBLE det;
  DOUBLE invert_last_point_factor = 1.; // 1 or -1
  DOUBLE D = DIMENSION; // effective dimensionality
  DOUBLE r0 = 1.; //effective repulsive interacction for spin1=spin2 
#ifdef JASTROW_DIMENSIONALITY_FROM_DPAR
  D = Dpar;
#endif
#ifdef JASTROW_EFFECTIVE_REPULSIVE_INTERACTION_FROM_R0PAR
  if (spin1 == spin2) {
	  r0 = R0par; //effective repulsive interacction for spin1=spin2 
	//  Message("  Effective repulsive interacction for spin1=spin2 from R0par");
  }
#endif

  *function_was_negative = OFF;

#ifdef TRIAL_1D
  if(first_time == ON) 
    if(flag_lnf_instead_of_f == ON)
      Message("  Scattering problem is solved for log of two-body wf.,\n  i.e. -u'' - (D-1)/r u' + [V(r)-E] u  - (u')^2 = 0");
    else
      Message("  Scattering problem is solved for the two-body wf. f,\n  i.e. -f'' - (D-1)/r f' + [V(r)-E] f = 0");
  first_time = OFF;
#endif
  if(flag_lnf_instead_of_f == OFF) { // Standard Schroedinger equation
    // -f'' - (D-1)/r f' + [V(r)-E] f = 0
    // - (f[i+1]+f[i-1]-2 f[i]) / dr^2 - (f[i+2]-f[i])/(2dr) + (V[i]-E) f[i] = 0
    // f'(L/2) = 0
    // f[N] = f[N-2]
    dr = G->step;
    dr2 = dr*dr;
    i = G->size-1;
    r = G->x[i];
    if(flag_right_boundary_condition == 0) { // zero derivative
      G->f[i] = 1.;
      //G->f[i-1] = G->f[i]*(1.+0.5*(InteractionEnergy_ij(r, spin1, spin2)-E)*dr2);
      G->f[i-1] = G->f[i]*(1.+0.5*(r0*InteractionEnergy_ij(r, spin1, spin2)-E)*dr2);
    } else if(flag_right_boundary_condition == 2) { // solution, corresponding to a bound state
#ifdef TRIAL_2D
      G->f[i] = K0(-fabs(sqrt(fabs(E)))*G->x[i]);
      G->f[i-1] = K0(-fabs(sqrt(fabs(E)))*G->x[i-1]);
#endif
    }
    else if(flag_right_boundary_condition == 3) { // solution, corresponding to phonons
#ifdef JASTROW_RIGHT_BOUNDARY_PHONONS // initialize accoring to the phononic wave function
#ifdef TRIAL_3D
      G->f[i] = exp(- Apar/(G->x[i]*G->x[i]) - Apar/((L-G->x[i])*(L-G->x[i])) + 2.*Apar/Lhalf2);
      G->f[i-1] = exp(- Apar/(G->x[i-1]*G->x[i-1]) - Apar/((L-G->x[i-1])*(L-G->x[i-1])) + 2.*Apar/Lhalf2);
#else
      G->f[i] = exp(- Apar/G->x[i] - Apar/(L-G->x[i]) + 2.*Apar/Lhalf);
      G->f[i-1] = exp(- Apar/G->x[i-1] - Apar/(L-G->x[i-1]) + 2.*Apar/Lhalf);
#endif
#endif
    } // otherwize f[i] and f[i-1] are already defined

    // otherwize f[i] and f[i-1] are already defined

    for(i=G->size-2; i>=1; i--) {
      r = G->x[i];
      //G->f[i-1] = (G->f[i]*(2.+(InteractionEnergy_ij(r, spin1, spin2)-E)*dr2)*r + G->f[i+1]*(-r-(D-1.)*dr*0.5))/(r-(D-1.)*dr*0.5);
      G->f[i-1] = (G->f[i]*(2.+(r0*InteractionEnergy_ij(r, spin1, spin2)-E)*dr2)*r + G->f[i+1]*(-r-(D-1.)*dr*0.5))/(r-(D-1.)*dr*0.5);
      if(G->f[i-1]<0) *function_was_negative = ON;
      //G->f[i-1] = (G->f[i]*(2.*r*r+(InteractionEnergy(r)*r*r-E*r*r)*dr2) + G->f[i+1]*r*(-r-(D-1.)*dr*0.5))/(r*r-(D-1.)*r*dr*0.5);
    }
  }
  else { // equation for u(r) = ln f(r)
    // -u'' - (D-1)/r u' + [V(r)-E] u  - (u')^2 = 0
    // - (u[i+1]-2u[i]+u[i-1])/dt2 - (D-2)/r (u[i+1]-u[i-1])/(2dt) + [V(i)-E] u[i] - [(u[i+1]-u[i-1]) / (2 dt)]^2 = 0
    // u'(L/2) = 0
    // u[N] = u[N-2]
    dr = G->step;
    dr2 = dr*dr;
    i = G->size-1;
    r = G->x[i];

    if(flag_right_boundary_condition == OFF) {
      // zero derivative for f' = exp(u)' = f u'
      //G->f[i] = 0;
      //G->f[i-1] = exp(G->f[i])*(1.+0.5*(InteractionEnergy(r)-E)*dr2);
      //G->f[i-1] = log(G->f[i-1]);
      // zero derivative for u
      G->f[i] = 1; // 0 gives a trivial and non-normalizable solution u = 0
      //G->f[i-1] = G->f[i]*(1.+0.5*(InteractionEnergy_ij(r, spin1, spin2)-E)*dr2);
      G->f[i-1] = G->f[i]*(1.+0.5*(r0*InteractionEnergy_ij(r, spin1, spin2)-E)*dr2);
    }
    else if(flag_right_boundary_condition == 2) { // solution, corresponding to a bound state
#ifdef TRIAL_2D
      G->f[i] = log(fabs(K0(-fabs(sqrt(fabs(E)))*G->x[i])));
      G->f[i-1] = log(fabs(K0(-fabs(sqrt(fabs(E)))*G->x[i-1])));
#endif
    }
    else if(flag_right_boundary_condition == 3) { // solution, corresponding to phonons
#ifdef JASTROW_RIGHT_BOUNDARY_PHONONS // initialize accoring to the phononic wave function
#ifdef TRIAL_3D
      G->f[i]   = - Apar/(G->x[i]*G->x[i]) - Apar/((L-G->x[i])*(L-G->x[i])) + 2.*Apar/Lhalf2;
      G->f[i-1] = - Apar/(G->x[i-1]*G->x[i-1]) - Apar/((L-G->x[i-1])*(L-G->x[i-1])) + 2.*Apar/Lhalf2;
#else
      G->f[i]   = - Apar/G->x[i] - Apar/(L-G->x[i]) + 2.*Apar/Lhalf;
      G->f[i-1] = - Apar/G->x[i-1] - Apar/(L-G->x[i-1]) + 2.*Apar/Lhalf;
#endif
#endif
    } // otherwize f[i] and f[i-1] are already defined

    for(i=G->size-2; i>=1; i--) {
      r = G->x[i];
      //det = 1. + 2.*(G->f[i]-G->f[i+1]) + G->f[i]*dr*dr*(InteractionEnergy_ij(r, spin1, spin2)-E);
      //det = 1. + 2.*(G->f[i]-G->f[i+1]) + dr*dr*(InteractionEnergy_ij(r, spin1, spin2)-E);
      det = 1. + 2.*(G->f[i]-G->f[i+1]) + dr*dr*(r0*InteractionEnergy_ij(r, spin1, spin2)-E);
      if(det<0) {
        Warning("Negative determinant in quadratic equation for log(f(x)) , i=%i\n", i);
        det = - det;
        *function_was_negative = ON;
        invert_last_point_factor = -1.;
        return -r;
      }
      G->f[i-1] = G->f[i+1] - 2.*(1. - sqrt(det)); // increasing solution
      //G->f[i-1] = G->f[i+1] - 2.*(1. + sqrt(det)); // second, non physical solution
    }

#ifdef JASTROW_GOES_TO_ONE_AT_LARGEST_DISTANCE
    for(i=0; i<G->size; i++) G->f[i] -= G->f[G->size-1]; // constant shift to ensure that f(L/2) = 1
#endif

    for(i=0; i<G->size; i++) G->f[i] = exp(G->f[i]); // convert back to standard notation
  }

  if(flag_left_boundary_condition == 0) // zero boundary condition
    return G->f[0]*invert_last_point_factor;
  else // zero derivative
    return (G->f[1] - G->f[0])*invert_last_point_factor;
}

#ifdef TRIAL_TONKS_TRAP
#ifndef HARD_SPHERE
//U = ln f
DOUBLE InterpolateExactU(DOUBLE x) {
  return Log(fabs(x-a));
}

//Fp = f' / f
DOUBLE InterpolateExactFp(DOUBLE x) {
  return 1./(x-a);
}

DOUBLE InterpolateExactE(DOUBLE x) {
#ifdef TRIAL_1D // 1D Eloc = [-f"/f] + (f'/f)^2 = (f'/f)^2
  return 1./((x-a)*(x-a));
#endif
#ifdef TRIAL_2D // 2D Eloc = [-(f" + f'/r)/f] + (f'/f)^2
  return a/((x-a)*(x-a)*x);
#endif
#ifdef TRIAL_3D // 3D Eloc = [-(f"+2f'/r)/f] + (f'/f)^2
  return (2.*a-x)/((x-a)*(x-a)*x);
#endif
}
#else
//U = ln f
DOUBLE InterpolateExactU(DOUBLE x) {
  if(x<a) return -100;
  return Log(fabs(x-a));
}

//Fp = f' / f
DOUBLE InterpolateExactFp(DOUBLE x) {
  if(x<a) return 10;
  return 1./(x-a);
}

DOUBLE InterpolateExactE(DOUBLE x) {
  if(x<a) return 0;
#ifdef TRIAL_1D // 1D Eloc = [-f"/f] + (f'/f)^2 = (f'/f)^2
  return 1./((x-a)*(x-a));
#endif
#ifdef TRIAL_2D // 2D Eloc = [-(f" + f'/r)/f] + (f'/f)^2
  return a/((x-a)*(x-a)*x);
#endif
#ifdef TRIAL_3D // 3D Eloc = [-(f"+2f'/r)/f] + (f'/f)^2
  return (2.*a-x)/((x-a)*(x-a)*x);
#endif
}
#endif // HARD_SPHERE
#endif // TRIAL_TONKS_TRAP

#ifdef TRIAL_CH2_3D_TRAP
DOUBLE InterpolateExactU(DOUBLE r) {
  return Log(tanh(r)/r);
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  return 2./sinh(2.*r)-1./r;
}

// Eloc = f" +2/r f' = u"/r
DOUBLE InterpolateExactE(DOUBLE r) {
  return (1.+4.*r*(r/tanh(2.*r)-1)/sinh(2.*r))/(r*r);
}
#endif

#ifdef TRIAL_CH2_1D_TRAP_FS
DOUBLE InterpolateExactU(DOUBLE r) {
  return Log(tanh(r));
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  return 2./sinh(2.*r);
}

// Eloc = f" = u"/r
DOUBLE InterpolateExactE(DOUBLE r) {
  return 4./(tanh(2.*r)*sinh(2.*r));
}
#endif

#ifdef TRIAL_CH2_1D_TRAP_BS
DOUBLE InterpolateExactU(DOUBLE r) {
  return Log(1/cosh(r));
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  return -tanh(r);
}

// Eloc = f" = u"/r
DOUBLE InterpolateExactE(DOUBLE r) {
  DOUBLE ch;
  ch = cosh(r);
  return 1./(ch*ch);
}
#endif

#ifdef TRIAL_MCGUIRE
//U = ln f
DOUBLE InterpolateExactU(DOUBLE r) {
#ifdef TRIAL_1D
 return -r/a;
#endif
#ifdef TRIAL_3D
 return -r/a - log(r);
#endif
}

//Fp = f' / f
DOUBLE InterpolateExactFp(DOUBLE r) {
#ifdef TRIAL_1D
  return -1./a;
#endif
#ifdef TRIAL_3D
  return -1./r - 1./a;
#endif
}

/* Eloc1D = -[f"/f - (f'/f)^2]*/
DOUBLE InterpolateExactE(DOUBLE r) {
#ifdef TRIAL_1D
  return 0.;
#endif
#ifdef TRIAL_3D
  return (1./r + 2./a)/r;
#endif
}
#endif

/************************** Tonks-Girardeau trapped 11 22 *****************/
#ifdef TRIAL_TONKS_TRAP11
#ifndef HARD_SPHERE11
//U = ln f
DOUBLE InterpolateExactU11(DOUBLE x) {
  return Log(fabs(x-b));
}

//Fp = f' / f
DOUBLE InterpolateExactFp11(DOUBLE x) {
  return 1./(x-b);
}

DOUBLE InterpolateExactE11(DOUBLE x) {
#ifdef TRIAL_1D // 1D Eloc = [-f"/f] + (f'/f)^2 = (f'/f)^2
  return 1./((x-b)*(x-b));
#endif
#ifdef TRIAL_2D // 2D Eloc = [-(f" + f'/r)/f] + (f'/f)^2
  return b/((x-b)*(x-b)*x);
#endif
#ifdef TRIAL_3D // 3D Eloc = [-(f"+2f'/r)/f] + (f'/f)^2
  return (2.*b-x)/((x-b)*(x-b)*x);
#endif
}
#else
//U = ln f
DOUBLE InterpolateExactU11(DOUBLE x) {
  if(x<b) return -100;
  return Log(fabs(x-b));
}

//Fp = f' / f
DOUBLE InterpolateExactFp11(DOUBLE x) {
  if(x<b) return 10;
  return 1./(x-b);
}

DOUBLE InterpolateExactE11(DOUBLE x) {
  if(x<b) return 0;
#ifdef TRIAL_1D // 1D Eloc = [-f"/f] + (f'/f)^2 = (f'/f)^2
  return 1./((x-b)*(x-b));
#endif
#ifdef TRIAL_2D // 2D Eloc = [-(f" + f'/r)/f] + (f'/f)^2
  return b/((x-b)*(x-b)*x);
#endif
#ifdef TRIAL_3D // 3D Eloc = [-(f"+2f'/r)/f] + (f'/f)^2
  return (2.*b-x)/((x-b)*(x-b)*x);
#endif
}
#endif // HARD_SPHERE
#endif // TRIAL_TONKS_TRAP1122


/************************** Construct SqWell depth ************************/
void ConstructSqWellDepth(void) {
#ifdef TRIAL_3D
  DOUBLE V,y,ymin,ymax,x,xmin,xmax;
  DOUBLE precision = 1e-8;
  int search = ON;
  DOUBLE R;
  DOUBLE kappa;

  R = RoSW; // SW size

  //calculate strength of the potential
  // define kappa
  // Solve tg x / x = 1 - a/R, E>V_0
  //       x = kappa R

  if(fabs(a)<1e9) {
    V = 1. - a/R;
    if(a>0) {
      xmin = (pi+precision)/2.;
      xmax = pi;
    }
    else {
      xmin = precision;
      xmax = (pi-precision)/2.;
    }
    ymin = tg(xmin)/xmin - V;
    ymax = tg(xmax)/xmax - V;
    if(ymin*ymax > 0) Error("R = %lf Can't construct the grid : No solution found", R);
    while(search) {
      x = (xmin+xmax)/2.;
      y = tg(x)/x - V;
      if(y*ymin < 0) {
        xmax = x;
        ymax = y;
      } else {
        xmin = x;
        ymin = y;
      }
      if(fabs(y)<precision) {
        x = (xmax*ymin-xmin*ymax) / (ymin-ymax);
        y = tg(x)/x - V;
        search = OFF;
      }
    }
    Message("  Solution is kappa R  = %f = pi/2 + %e, error %e\n", x, x-pi/2, y);
    kappa = x / R;
    Message("  Potential depth %e\n", kappa*kappa);
  }
  else {
    Warning("  Unitarity limit\n");
    kappa = pi/(2.*R) + precision;
  }

  Message("  Check a/R = %e\n", 1-tg(kappa*R)/(kappa*R));
  Vo = kappa*kappa;
  trial_Vo = Vo;
  Message("  Potential depth Vo = %e\n", Vo);
  if(a<0) {
    Warning("Negative scattering length\n");
    a = -a;
  }
#else // 2D system
  DOUBLE V,y,ymin,ymax,x,xmin,xmax;
  DOUBLE precision = 1e-8;
  int search = ON;
  DOUBLE R;
  DOUBLE kappa;

  R = RoSW;
  //calculate strength of the potential
  // define kappa
  // Solve  exp{ Jo(x) / (x J1(x))} = a/R
  //       x = kappa R

  if(fabs(a)<1e9) {
    V = a/R;
    xmin = 3.8317059702075107; // a = 0
    xmax = 0.4;
    ymin = exp(BesselJ0(xmin)/(xmin*BesselJ1(xmin))) - V;
    ymax = exp(BesselJ0(xmax)/(xmax*BesselJ1(xmax))) - V;
    if(ymin*ymax > 0) Error("R = %f Can't construct the grid : No solution found", R);
    while(search) {
      x = (xmin+xmax)/2.;
        y = exp(BesselJ0(x)/(x*BesselJ1(x))) - V;
      if(y*ymin < 0) {
        xmax = x;
        ymax = y;
      } else {
        xmin = x;
        ymin = y;
      }
      if(fabs(y)<precision) {
        x = (xmax*ymin-xmin*ymax) / (ymin-ymax);
        y = exp(BesselJ0(x)/(x*BesselJ1(x))) - V;
        search = OFF;
      }
      //Message("x = %g y = %g\n", x, y);
    }
    Message("  Solution is kappa R  = %f, error %e\n", x, y);
    kappa = x / R;
    Message("  Potential depth %e (%e in EF)\n", kappa*kappa, kappa*kappa*energy_unit);
  }
  else {
    Message("  Unitarity limit\n");
    //kappa = pi/2. + precision;
    kappa = 3.8317059702075107 + precision;
  }

  Message("  Check a/R = %e\n",  exp(BesselJ0(kappa)/(kappa*BesselJ1(kappa))));
  Vo = kappa*kappa;
  trial_Vo = Vo;
#endif
}

/************************** Construct SqWell depth ************************/
void ConstructVextSqWellDepth(void) {
  DOUBLE V,y,ymin,ymax,x,xmin,xmax;
  DOUBLE precision = 1e-8;
  int search = ON;
  DOUBLE R;
  DOUBLE kappa;
  DOUBLE a;

  a = b;
  Message("Constructing square well of diameter RoSW=%lf and s-wave scattering length b=%lf\n", RoSW, a);
  R = RoSW;

  //calculate strength of the potential
  // Solve  tg x / x = 1 - a/R, i.e. zero energy scattering solution
  //       x = kappa R
  //       kappa = sqrt(2 mu Vo ) / hbar
  if(fabs(a)<1e9) {
    V = 1. - a/R;
    if(a>0) {
      xmin = (pi+precision)/2.;
      xmax = pi;
    }
    else {
      xmin = precision;
      xmax = (pi-precision)/2.;
    }
    ymin = tg(xmin)/xmin - V;
    ymax = tg(xmax)/xmax - V;
    if(ymin*ymax > 0) 
      Error("R = %"LF" Can't construct square well", R);
    while(search) {
      x = (xmin+xmax)/2.;
      y = tg(x)/x - V;
      if(y*ymin < 0) {
        xmax = x;
        ymax = y;
      } else {
        xmin = x;
        ymin = y;
      }
      if(fabs(y)<precision) {
        x = (xmax*ymin-xmin*ymax) / (ymin-ymax);
        y = tg(x)/x - V;
        search = OFF;
      }
    }
    Message("  Solution is kappa R  = %"LF" = pi/2 + %"LE", error %"LE"\n", x, x-pi/2, y);
    kappa = x / R;
    Message("  Potential depth %"LE"\n", kappa*kappa);
  }
  else {
    Warning("  Unitarity limit\n");
    kappa = pi/(2.*R) + precision;
  }

  Message("  Check a = %"LE", goal value %"LE"\n", R*(1.-tg(kappa*R)/(kappa*R)), a);

  //Message(", mu = %"LE"\n", mu);
  // kappa = sqrt(2 mu Vo ) / hbar
  // Vo = \hbar^2 kappa^2 / (2 mu)
  trial_Vo = kappa*kappa/(2.*m_mu);
  trial_Vo = Vo;
  Message("  Potential depth Vo= %"LE"\n", trial_Vo);
}

/************************** Construct Grid Fermion ************************/
DOUBLE trial_Kappa, trial_Kappa2;
DOUBLE A1SW,A3SW,A4SW;
void ConstructGridSquareWellZeroConst(void) {
  static int first_time = ON;
  DOUBLE r;

  trial_Kappa = sqrt(Vo);
  trial_Kappa2 = Vo;

  if(first_time) {
    Message("  BCS square well zero energy + const\n");
    Warning("  Scattering length must be negative in order to use this w.f.\n");
    Warning("  Changing Rpar %e -> %e\n", Rpar, Rpar*Lhalf);
    first_time = OFF;
    Rpar *= Lhalf;
  }

  //A3SW = 1./(1.+a/Rpar*(1.-1./(Rpar*Bpar)*(exp(Bpar*(L-2*Rpar))+1.)/(exp(Bpar*(L-2*Rpar))-1.)));
  //A1SW = A3SW*(1.+a)/sin(trial_Kappa);
  //A4SW = A3SW*a/(Rpar*Rpar*Bpar)/(exp(-Bpar*Rpar)-exp(Bpar*(Rpar-L)));

  r = Lhalf;
  A4SW =  -1./(exp(-Bpar*r)+exp(Bpar*(r-L)));
  r = Rpar;
  A3SW = (1. + A4SW*(exp(-Bpar*r)+exp(Bpar*(r-L))))/(1.+a/r);
  r = 1.;
  A1SW = A3SW*(1.+a/r)/(sin(trial_Kappa*r)/r);
}

#ifdef TRIAL_SQ_WELL_ZERO_CONST
DOUBLE InterpolateExactU(DOUBLE r) {
  if(r<1.)
    return log(A1SW*sin(trial_Kappa*r)/r);
  else if(r<Rpar)
    return log(A3SW*(1.+a/r));
  else if(r<Lhalf)
    return log(1. + A4SW*(exp(-Bpar*r)+exp(Bpar*(r-L))));
  else
    return 0;
}

// Fp = f'/f
DOUBLE InterpolateExactFp(DOUBLE r) {
  if(r<1.)
    return trial_Kappa/tan(trial_Kappa*r) - 1./r;
  else if(r<Rpar)
    return -a/(r*r)/(1.+a/r);
  else if(r<Lhalf)
    return A4SW*Bpar*(-exp(-Bpar*r)+exp(Bpar*(r-L)))/(1. + A4SW*(exp(-Bpar*r)+exp(Bpar*(r-L))));
  else
    return 0.;
}

// Eloc = -(f" +2/r f')/f+(f'/f)^2
DOUBLE InterpolateExactE(DOUBLE r) {
  DOUBLE fp;

  fp = InterpolateExactFp(r);
  if(r<1.)
    return trial_Kappa2 + fp*fp;
  else if(r<Rpar)
    return fp*fp;
  else if(r<Lhalf)
    return -((Bpar*A4SW*(exp(Bpar*L)*(-2 + Bpar*r) + exp(2*Bpar*r)*(2 + Bpar*r)))/((exp(Bpar*(L + r)) + A4SW*(exp(Bpar*L) + exp(2*Bpar*r)))*r)) + fp*fp;
  else
    return 0;
}
#endif

/************************** Construct Grid SW BS phonons ******************/
void ConstructGridSquareWellBoundStatePhonons(struct Grid *G) {
  DOUBLE kappa, k;
  DOUBLE V,y,ymin,ymax,x,xmin,xmax;
  DOUBLE precision = 1e-10;
  int search = ON;
  DOUBLE detuning = 1e-10;

  kappa = sqrt(Vo);

  V = -1;

  Message(" Jastrow term: SquareWellBoundStatePhonons\n");

  //xmin = detuning*pi/2.*0.1;
  //xmax = detuning*pi/2.*10.;
  xmin = 1e-10;
  xmax = 1e0;

  // solve: sqrt(kappa^2-k^2) + k tg(sqrt(kappa^2-k^2)) = 0
  ymin = sqrt(-xmin*xmin+kappa*kappa)/tg(sqrt(-xmin*xmin+kappa*kappa))/xmin - V;
  ymax = sqrt(-xmax*xmax+kappa*kappa)/tg(sqrt(-xmax*xmax+kappa*kappa))/xmax - V;

  if(ymin*ymax > 0) Error("Can't construct the grid : No solution found");

  search = ON;
  while(search) {
    x = (xmin+xmax)/2.;
    y = sqrt(-x*x+kappa*kappa)/tg(sqrt(-x*x+kappa*kappa))/x - V;

    if(y*ymin < 0) {
      xmax = x;
      ymax = y;
    } else {
      xmin = x;
      ymin = y;
    }

    if(fabs(y)<precision) {
      x = (xmax*ymin-xmin*ymax) / (ymin-ymax);
      y = sqrt(-x*x+kappa*kappa)/tg(sqrt(-x*x+kappa*kappa))/x - V;
      search = OFF;
    }
  }
  k = x;
  Message("  Solution k = %e, error %e (k_approx = %e)\n", k, y, detuning*pi/2.);

  sqE = k;
  sqE2 = k*k;
  Message("  Scattering energy %e, approx %e\n", sqE2, (detuning*pi/2.)*(detuning*pi/2.));
  trial_Kappa2 = kappa*kappa - k*k;
  trial_Kappa = sqrt(trial_Kappa2);

  // matching at Rpar
  Warning("  Matching distance Rpar = %lf [L/2] -> %lf [R]\n", Rpar, Rpar*Lhalf);
  if(Rpar>1 || Rpar <0) Error(" Rpar out of (0,1) range!\n");
  Rpar *= Lhalf;
  A5SW = -(1 + k*Rpar)/(2.*Rpar*(1/(Rpar*Rpar*Rpar) - 1/((L-Rpar)*(L-Rpar)*(L-Rpar))));

  // A5SW -> Ctrial
  A2SW = exp(2.*A5SW/(Lhalf*Lhalf));
  A3SW = Rpar*exp(k*Rpar)*A2SW*exp(-A5SW/(Rpar*Rpar)-A5SW/((L-Rpar)*(L-Rpar)));
  A1SW = A3SW*exp(-sqE)/sin(trial_Kappa);

  // take log of coef for U(r)
  A3SW = log(A3SW);
  A2SW = log(A2SW);

  G->min = 0;
  G->max = Lhalf;
  G->max2 = G->max*G->max;
}

#ifdef TRIAL_SQUARE_WELL_BS_PHONONS
DOUBLE InterpolateExactU(DOUBLE r) {
  if(r<1)
    return log(fabs(A1SW*sin(trial_Kappa*r)/r));
  else if(r<Rpar)
    return A3SW-sqE*r-log(r);
  else
    return A2SW - A5SW*(1./(r*r)+1./((Lwf-r)*(Lwf-r)));
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  if(r<1)
    return trial_Kappa/tg(trial_Kappa*r) - 1./r;
  else if(r<Rpar)
    return - sqE - 1./r;
  else
    return 2.*A5SW*(1./(r*r*r)-1./((Lwf-r)*(Lwf-r)*(Lwf-r)));
}

DOUBLE InterpolateExactE(DOUBLE r) {
  DOUBLE fp;

  if(r<1) {
    fp = trial_Kappa/tg(trial_Kappa*r) - 1./r;
    return trial_Kappa2 + fp*fp;
  }
  else if(r<Rpar){
    fp = -sqE - 1./r;
    return -sqE2 + fp*fp;
  }
  else
    return 2.*A5SW*(Lwf*Lwf*Lwf*Lwf - 4.*Lwf*Lwf*Lwf*r + 6.*Lwf*Lwf*r*r - 2.*Lwf*r*r*r + 2.*r*r*r*r)/(r*r*r*r*(Lwf-r)*(Lwf-r)*(Lwf-r)*(Lwf-r));
}
#endif


/************************** Construct Grid Square Well Trap ***************/
//  note that a>0 is actually -a
void ConstructGridSquareWellTrap(void) {
  static int first_time = ON;
  DOUBLE r;
  DOUBLE mu;
  DOUBLE precision = 1e-8;

  // reduced mass
  mu = m_mu;

  // kappa = sqrt(2 mu Vo ) / hbar
  trial_Kappa = sqrt(2.*mu*trial_Vo);
  trial_Kappa2 = 2.*mu*trial_Vo;

  if(first_time) {
    Warning("  Scattering length must be negative in order to use this w.f.\n");
    Warning("  Sin(K r)/r [r<RoSW], 1+a/r\n");
    first_time = OFF;
  }

  if(a < 0) {
    Warning("  Automatically changing scattering length to positive value\n");
    a = -a;
  }

  r = RoSW;
  A3SW = 1.;
  A1SW = A3SW*(r+a)/Sin(trial_Kappa*r); // f(r)
  r = (trial_Kappa / tan(trial_Kappa*r) - 1. / r) / (-a / (r*r) / (1. + a / r)); // testing continuity of the first derivative
  if(fabs(r-1)>0.1) Error("ConstructGridSquareWellZeroConst: first derivative is discontinuous");
}

#ifdef TRIAL_SQUARE_WELL_FS_TRAP
DOUBLE InterpolateExactU(DOUBLE r) {
  if(r<RoSW)
    return Log(A1SW*Sin(trial_Kappa*r)/r);
  else
    return Log(A3SW*(1.+a/r));
}

// Fp = f'/f
DOUBLE InterpolateExactFp(DOUBLE r) {
  if(r<RoSW)
    return trial_Kappa/tan(trial_Kappa*r) - 1./r;
  else
    return -a/(r*r)/(1.+a/r);
}

// Eloc = -(f" +2/r f')/f+(f'/f)^2
DOUBLE InterpolateExactE(DOUBLE r) {
  DOUBLE fp;

  fp = InterpolateExactFp(r);
  if(r<RoSW)
    return trial_Kappa2 + fp*fp;
  else
    return fp*fp;
}
#endif

#ifdef TRIAL_PP_FS_TRAP_MCGUIRE1D
DOUBLE InterpolateExactU(DOUBLE r) {
  return log(1. - a/r) - r*r/(b*(Apar + r));
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  return a/(r*(-a + r)) + (-1. + Apar*Apar/((Apar + r)*(Apar + r)))/b;
}

DOUBLE InterpolateExactE(DOUBLE r) {
  return 3./((a-r)*(a-r)) + 1./(r*r) + 2./(a*r) - 2.*r/(a*(a - r)*(a - r)) + (2.*(3.*Apar*Apar + 3*Apar*r + r*r))/(b*(Apar + r)*(Apar + r)*(Apar + r));
}
#endif

/************************** Construct Grid Square Well Trap ***************/
//  note that a>0 is actually -a
void ConstructGridSquareWellTrapMcGuire1D(void) {
  static int first_time = ON;
  DOUBLE r;
  DOUBLE mu;
  DOUBLE precision = 1e-8;

  // reduced mass
  mu = m_mu;

  // kappa = sqrt(2 mu Vo ) / hbar
  trial_Kappa = sqrt(2.*mu*trial_Vo);
  trial_Kappa2 = 2.*mu*trial_Vo;

  if(first_time) {
    Warning("  Scattering length must be negative in order to use this w.f.\n");
    Warning("  Sin(K r)/r [r<RoSW], 1+a/r\n");
    first_time = OFF;
  }

  if(a < 0) {
    Warning("  Automatically changing scattering length to positive value\n");
    a = -a;
  }

  r = RoSW;
  Atrial = (-(Apar*Apar*b) + 2*Apar*r*(-b + r) + r*r*(-b + r) + b*r*(Apar + r)*(Apar + r)/tan(r*trial_Kappa)*trial_Kappa)/(r*(2*Apar + r) + b*(Apar + r)*(Apar + r)/tan(r*trial_Kappa)*trial_Kappa);
  A1SW = ((1. - Atrial/r) - r*r/(b*(Apar + r))) / Sin(trial_Kappa*r)/r;
}

#ifdef TRIAL_SQUARE_WELL_FS_TRAP_MCGUIRE1D
DOUBLE InterpolateExactU(DOUBLE r) {
  if(r<RoSW)
    return log(A1SW*Sin(trial_Kappa*r)/r);
  else
    return log(1. - Atrial/r) - r*r/(b*(Apar + r));
}

// Fp = f'/f
DOUBLE InterpolateExactFp(DOUBLE r) {
  if(r<RoSW)
    return trial_Kappa/tan(trial_Kappa*r) - 1./r;
  else
    return Atrial/(r*(-Atrial + r)) + (-1. + Apar*Apar/((Apar + r)*(Apar + r)))/b;
}

// Eloc = -(f" +2/r f')/f+(f'/f)^2
DOUBLE InterpolateExactE(DOUBLE r) {
  DOUBLE fp;

  if(r<RoSW) {
    fp = InterpolateExactFp(r);
    return trial_Kappa2 + fp*fp;
  }
  else
    return 3./((Atrial-r)*(Atrial-r)) + 1./(r*r) + 2./(Atrial*r) - 2.*r/(Atrial*(Atrial - r)*(Atrial - r)) + (2.*(3.*Apar*Apar + 3*Apar*r + r*r))/(b*(Apar + r)*(Apar + r)*(Apar + r));
}
#endif


/************************** Construct Grid SW BS trap ******************/
void ConstructGridSquareWellBoundStateTrap(struct Grid *G) {
  DOUBLE kappa, k;
  DOUBLE V,y,ymin,ymax,x,xmin,xmax;
  DOUBLE precision = 1e-10;
  int search = ON;
  DOUBLE detuning = 1e-10;
  DOUBLE r;
  DOUBLE m_mu = 0.5;

  kappa = sqrt(2.*m_mu*trial_Vo);
  V = -1;

  Message(" Jastrow term: Square Well Bound State Trap\n");

  //xmin = detuning*pi/2.*0.1;
  //xmax = detuning*pi/2.*10.;
  xmin = 1e-10;
  xmax = 1e0;

  // solve: sqrt(kappa^2-k^2) + k tg(sqrt(kappa^2-k^2)) = 0
  ymin = sqrt(-xmin*xmin+kappa*kappa)/tg(RoSW*sqrt(-xmin*xmin+kappa*kappa))/xmin - V;
  ymax = sqrt(-xmax*xmax+kappa*kappa)/tg(RoSW*sqrt(-xmax*xmax+kappa*kappa))/xmax - V;

  if(ymin*ymax > 0) Error("Can't construct the grid : No solution found");

  search = ON;
  while(search) {
    x = (xmin+xmax)/2.;
    y = sqrt(-x*x+kappa*kappa)/tg(RoSW*sqrt(-x*x+kappa*kappa))/x - V;

    if(y*ymin < 0) {
      xmax = x;
      ymax = y;
    } else {
      xmin = x;
      ymin = y;
    }

    if(fabs(y)<precision) {
      x = (xmax*ymin-xmin*ymax) / (ymin-ymax);
      y = sqrt(-x*x+kappa*kappa)/tg(RoSW*sqrt(-x*x+kappa*kappa))/x - V;
      search = OFF;
    }
  }
  k = x;
  Message("  Solution k = %e, error %e (k_approx = %e)\n", k, y, detuning*pi/2.);

  sqE = k;
  sqE2 = k*k;
  Message("  Scattering energy %e, approx %e\n", sqE2, (detuning*pi/2.)*(detuning*pi/2.));
  trial_Kappa2 = kappa*kappa - k*k;
  trial_Kappa = sqrt(trial_Kappa2);

  // f(RoSW) matching
  r = RoSW;
  A3SW = log(fabs(sin(trial_Kappa*r)/r)) - (-sqE*r-log(r));

  G->min = 0;
  G->max = Lhalf;
  G->max2 = G->max*G->max;
}

#ifdef TRIAL_SQUARE_WELL_BS_TRAP
DOUBLE InterpolateExactU(DOUBLE r) {
  if(r<RoSW)
    return log(fabs(sin(trial_Kappa*r)/r));
  else
    return A3SW-sqE*r-log(r);
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  if(r<RoSW)
    return trial_Kappa/tg(trial_Kappa*r) - 1./r;
  else
    return - sqE - 1./r;
}

DOUBLE InterpolateExactE(DOUBLE r) {
  DOUBLE fp;

  if(r<RoSW) {
    fp = trial_Kappa/tg(trial_Kappa*r) - 1./r;
    return trial_Kappa2 + fp*fp;
  }
  else {
    fp = -sqE - 1./r;
    return -sqE2 + fp*fp;
  }
}
#endif

/************************** Construct Grid Fermion ************************/
void ConstructGridSquareWellFreeState(struct Grid *G) {
  DOUBLE kappa, k;
  DOUBLE V,y,ymin,ymax,x,xmin,xmax;
  DOUBLE precision = 1e-10;
  DOUBLE trial_delta_check;
  int search = ON;
  //DOUBLE Lhalf = 15;

  // matching conditions
  V = Lhalf - 1;

#ifdef BC_ABSENT
  Warning("  f'(Rpar) = 0, Rpar = %lf\n", Rpar);
  Lhalf = Rpar;
#endif

  //kappa = sqrt(trial_Vo);
  kappa = sqrt(2.*m_mu*trial_Vo);
  xmin = 1e-8; //1e-4
  xmax = 100; // 4
  ymin = (arctg(xmin*Lhalf) - arctg(xmin/sqrt(xmin*xmin+kappa*kappa)*tg(sqrt(xmin*xmin+kappa*kappa))))/xmin - V;
  ymax = (arctg(xmax*Lhalf) - arctg(xmax/sqrt(xmax*xmax+kappa*kappa)*tg(sqrt(xmax*xmax+kappa*kappa))))/xmax - V;

  if(ymin*ymax > 0) Error("ConstructGridSquareWellFreeState: cannot construct the grid : No solution found");
  while(search) {
    x = (xmin+xmax)/2.;
    y = (arctg(x*Lhalf) - arctg(x/sqrt(x*x+kappa*kappa)*tg(sqrt(x*x+kappa*kappa))))/x - V;

    if(y*ymin < 0) {
      xmax = x;
      ymax = y;
    } else {
      xmin = x;
      ymin = y;
    }

    if(fabs(y)<precision) {
      x = (xmax*ymin-xmin*ymax) / (ymin-ymax);
      y = (arctg(x*Lhalf) - arctg(x/sqrt(x*x+kappa*kappa)*tg(sqrt(x*x+kappa*kappa))))/x - V;
      search = OFF;
    }
  }
  k = x;

  Message("  Solution k = %"LE", error %"LE"\n", k, y);

  sqE = k;
  sqE2 = k*k;
  trial_Vo = kappa*kappa/(2.*m_mu);
  trial_Kappa2 = k*k + kappa*kappa;
  trial_Kappa = sqrt(trial_Kappa2);

  trial_delta = arctg(k*Lhalf) - k*Lhalf;
  trial_delta_check = arctg(k/trial_Kappa*tg(trial_Kappa)) - k;
  if(fabs((trial_delta_check-trial_delta)/trial_delta)>1e-8) Error("  Invalid solution\n");

  A3SW = Lhalf/Sin(k*Lhalf+trial_delta);
  A1SW = A3SW*Sin(k*RoSW+trial_delta)/Sin(trial_Kappa*RoSW);

  //Xi = A1SW*Sin(trial_Kappa);  Xi = A3SW*Sin(k+trial_delta);
  Message(" Matching to one -> %"LF"\n", A3SW*Sin(k*Lhalf+trial_delta)/Lhalf);

  G->min = 0;
  G->max = Lhalf;
  G->max2 = G->max*G->max;
}

#ifdef TRIAL_SQUARE_WELL_FS
DOUBLE InterpolateExactU(DOUBLE r) {
  if(r<RoSW)
    return Log(fabs(A1SW*Sin(trial_Kappa*r)/r));
  else if(r<Lhalf)
    return Log(fabs(A3SW*Sin(sqE*r+trial_delta)/r));
  else
    return 0;
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  if(r<RoSW)
    return trial_Kappa/tg(trial_Kappa*r) - 1./r;
  else if(r<Lhalf)
    return sqE/tg(sqE*r+trial_delta) - 1./r;
  else
    return 0;
}

DOUBLE InterpolateExactE(DOUBLE r) {
  DOUBLE fp=0;

  if(r<RoSW) {
    fp = trial_Kappa/tg(trial_Kappa*r) - 1./r;
    return trial_Kappa2 + fp*fp;
  }
  else if(r<Lhalf) {
    fp = sqE/tg(sqE*r+trial_delta) - 1./r;
    return sqE2 + fp*fp;
  }
  else
    return 0;
}
#endif

/************************** Construct Grid SW BS phonons ******************/
void ConstructGridZeroRangeUnitaryPhonons(struct Grid *G) {
  DOUBLE r;

  Message(" wf 12: Zero range unitrary\n");

#ifdef BC_ABSENT
  Message("  trapped geometry, f2(r) = 1/r\n");
#else
  // matching at Rpar
  Warning("  Matching distance Rpar = %lf [L/2] -> %lf [R]\n", Rpar, Rpar*Lhalf);
  if(Rpar>1 || Rpar <0) Error(" Rpar out of (0,1) range!\n");
  Rpar *= Lhalf;

  // f'(Rpar)
  r = Rpar;
  A5SW = - 1./r / (2.*(1./(r*r*r)-1./((Lwf-r)*(Lwf-r)*(Lwf-r))));

  // f(L/2)
  r = Lhalfwf;
  A2SW = A5SW*(1./(r*r)+1./((Lwf-r)*(Lwf-r)));

  // f(Rpar)
  r = Rpar;
  A1SW =log(r) + A2SW - A5SW*(1./(r*r)+1./((Lwf-r)*(Lwf-r)));
#endif

  G->min = 0;
  G->max = Lhalf;
  G->max2 = G->max*G->max;
}

#ifdef TRIAL_ZERO_RANGE_UNITARY
DOUBLE InterpolateExactU(DOUBLE r) {
#ifdef BC_ABSENT
  return -log(r);
#else
  if(r<Rpar)
    return A1SW-log(r);
  else
    return A2SW - A5SW*(1./(r*r)+1./((Lwf-r)*(Lwf-r)));
#endif
}

DOUBLE InterpolateExactFp(DOUBLE r) {
#ifdef BC_ABSENT
  return - 1./r;
#else
  if(r<Rpar)
    return - 1./r;
  else
    return 2.*A5SW*(1./(r*r*r)-1./((Lwf-r)*(Lwf-r)*(Lwf-r)));
#endif
}

DOUBLE InterpolateExactE(DOUBLE r) {
#ifdef BC_ABSENT
  return 1./(r*r);
#else
  DOUBLE fp;

  if(r<Rpar) {
    return 1./(r*r);
  }
  else
    return 2.*A5SW*(Lwf*Lwf*Lwf*Lwf - 4.*Lwf*Lwf*Lwf*r + 6.*Lwf*Lwf*r*r - 2.*Lwf*r*r*r + 2.*r*r*r*r)/(r*r*r*r*(Lwf-r)*(Lwf-r)*(Lwf-r)*(Lwf-r));
#endif
}
#endif

#ifdef TRIAL_ZERO_RANGE_TRAP
// f(r) = (r-a)/r = 1 - a/r
DOUBLE InterpolateExactU(DOUBLE r) {
#ifdef HARD_SPHERE
  if(r<a) return -100;
#endif
  return log(fabs(1.-a/r));
}

DOUBLE InterpolateExactFp(DOUBLE r) {
#ifdef HARD_SPHERE
  if(r<a) return 100;
#endif
  return a/(r*(r-a));
}

DOUBLE InterpolateExactE(DOUBLE r) {
// [-(f"+2f'/r)/f] + (f'/f)^2
#ifdef HARD_SPHERE
  if(r<a) return 100;
#endif
  return a/(r*(r-a))*a/(r*(r-a));
}
#endif

/************************** Construct Grid SW BS phonons ******************/
void ConstructGridZeroRange3D(struct Grid *G) {
  DOUBLE r;
  DOUBLE y,ymin,ymax,x,xmin,xmax;
  DOUBLE precision = 1e-10;
  DOUBLE factor = 1.1;
  int search = ON;

  Message(" wf: Zero range 3D\n");

  r = Lhalf;
  // Bethe Peierls
  //trial_delta = - atan(Ktrial*a);
  xmin = 1e-6;
  ymin = -1./r + xmin/tan(xmin*r- atan(xmin*a));
  xmax = xmin;
  ymax = -1./r + xmax/tan(xmax*r- atan(xmax*a));
  while(ymin*ymax > 0.) {
    xmax = xmax*factor;
    ymax = -1./r + xmax/tan(xmax*r- atan(xmax*a));
  }

  if(ymin*ymax > 0) Error("ConstructGridZeroRange3D: cannot construct the grid : No solution found");
  while(search) {
    x = (xmin+xmax)/2.;
    y = -1./r + x/tan(x*r- atan(x*a));

    if(y*ymin < 0) {
      xmax = x;
      ymax = y;
    } else {
      xmin = x;
      ymin = y;
    }

    if(fabs(y)<precision) {
      x = (xmax*ymin-xmin*ymax) / (ymin-ymax);
      y = -1./r + x/tan(x*r- atan(x*a));
      search = OFF;
    }
  }
  Ktrial = x;
  trial_delta = - atan(Ktrial*a);

  Atrial = 1./ (sin(Ktrial*r+trial_delta)/r);
  Message("  Solution k = %"LE", error %"LE", delta = %lf, A = %lf\n", Ktrial, y, trial_delta, Atrial);

  G->min = 0;
  G->max = Lhalf;
  G->max2 = G->max*G->max;
}

#ifdef TRIAL_ZERO_RANGE_3D
DOUBLE InterpolateExactU(DOUBLE r) {
  if(r<Lhalf)
    return log(fabs(Atrial*sin(Ktrial*r+trial_delta)/r));
  else
    return 0;
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  if(r<Lhalf)
    return -1./r + Ktrial/tan(Ktrial*r+trial_delta);
  else
    return 0;
}

DOUBLE InterpolateExactE(DOUBLE r) {
  DOUBLE Fp;
  if(r<Lhalf) {
    Fp = -1./r + Ktrial/tan(Ktrial*r+trial_delta);
    return Ktrial*Ktrial + Fp*Fp;
  }
  else
    return 0;
}
#endif

/************************** Construct Grid SW FS phonons ******************/
void ConstructGridSquareWellFreeStatePhonons(struct Grid *G) {
  DOUBLE kappa, k, r;

  Message(" Jastrow term: SquareWellFreeStatePhonons\n");
  kappa = sqrt(Vo);

  Message(" Scattering energy parameter k = Apar/(2*pi*L)\n");
  k = Apar/(2.*PI*Lwf);
  if(Apar<1e-8) Error("Set Apar to a finite value!\n");

  sqE = k;
  sqE2 = k*k;
  trial_Kappa2 = kappa*kappa + k*k;
  trial_Kappa = sqrt(trial_Kappa2);

  // matching at Rpar
  if(Rpar>1 || Rpar <1e-8) Error(" Rpar out of (0,1) range!\n");
  Warning("  Matching distance Rpar = %lf [L/2] -> %lf [R]\n", Rpar, Rpar*Lhalf);
  Rpar *= Lhalf;

  // f'(1)
  r = 1.;
  trial_delta = atan(sqE*tg(trial_Kappa*r)/trial_Kappa) - sqE*r;
  if(trial_delta + sqE*Rpar<0) trial_delta += PI; // choose proper phase

  // f'(Rpar)
  r = Rpar;
  A5SW = (sqE/tg(trial_delta + sqE*r) - 1./r)/(2.*(1./(r*r*r)-1./((Lwf-r)*(Lwf-r)*(Lwf-r))));

  // f(L/2)
  r = Lhalfwf;
  A2SW = A5SW*(1./(r*r)+1./((Lwf-r)*(Lwf-r)));

  // f(Rpar)
  r = Rpar;
  A3SW = A2SW - A5SW*(1./(r*r)+1./((Lwf-r)*(Lwf-r))) - log(Sin(trial_delta + sqE*r)/r);

  // f(1)
  r = 1;
  A1SW = A3SW + log(Sin(trial_delta + sqE*r)/r) - log(Sin(trial_Kappa*r)/r);

  G->min = 0;
  G->max = Lhalf;
  G->max2 = G->max*G->max;
}

#ifdef TRIAL_SQUARE_WELL_FS_PHONONS
DOUBLE InterpolateExactU(DOUBLE r) {
  if(r<1)
    return A1SW+log(Sin(trial_Kappa*r)/r);
  else if(r<Rpar)
    return A3SW + log(Sin(trial_delta + sqE*r)/r);
  else
    return A2SW - A5SW*(1./(r*r)+1./((Lwf-r)*(Lwf-r)));
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  if(r<1)
    return trial_Kappa/tg(trial_Kappa*r) - 1./r;
  else if(r<Rpar)
    return sqE/tg(trial_delta + sqE*r) - 1./r;
  else
    return 2.*A5SW*(1./(r*r*r)-1./((Lwf-r)*(Lwf-r)*(Lwf-r)));
}

DOUBLE InterpolateExactE(DOUBLE r) {
  DOUBLE fp;

  if(r<1) {
    fp = trial_Kappa/tg(trial_Kappa*r) - 1./r;
    return trial_Kappa2 + fp*fp;
  }
  else if(r<Rpar){
    fp = sqE/tg(trial_delta + sqE*r) - 1./r;
    return sqE2 + fp*fp;
  }
  else
    return 2.*A5SW*(Lwf*Lwf*Lwf*L - 4.*Lwf*Lwf*Lwf*r + 6.*Lwf*Lwf*r*r - 2.*Lwf*r*r*r + 2.*r*r*r*r)/(r*r*r*r*(Lwf-r)*(Lwf-r)*(Lwf-r)*(Lwf-r));
}
#endif

/************************** Construct Grid Quadrupole **********************/
void ConstructGridQuadrupoleKo(struct Grid *G) {
  AllocateWFGrid(G, 3);

  Message("2D Quadrupole Ko(r) trial wavefunction\n");
  Message("  Rpar %" LG "  [Lwf/2] -> %" LG "  [R]\n", Rpar, Rpar*Lhalfwf);
  Rpar *= Lhalfwf;
//Lwf*=sqrt(N_Rc/N);Lhalfwf*=sqrt(N_Rc/N);Lhalfwf2*=N_Rc/N;
  Atrial = 1./BesselK0(2/3./Rpar/Sqrt(Rpar));
  Btrial = 1./BesselK1(2/3./Rpar/Sqrt(Rpar));
  alpha_1 = Atrial/Btrial*(Lwf-2*Rpar)*(Lwf-Rpar)/(Lwf*Lwf*Rpar*Sqrt(Rpar));
  Xi = BesselK1(2/3./Rpar/Sqrt(Rpar))*(Lwf-Rpar)*(Lwf-Rpar)/Sqrt(Rpar)/(Lwf*Lwf*(Lwf-2*Rpar)*BesselK0(2/3./Rpar/Sqrt(Rpar)));

  G->min = 0;
  G->max = Lhalfwf;
  G->max2 = G->max*G->max;
}

// executed only if TRIAL_QUADRUPOLE_Ko is defined
#ifdef TRIAL_QUADRUPOLE_Ko
DOUBLE InterpolateExactU(DOUBLE r) {
  DOUBLE K0, K1;
  if(r<Rpar) {
    K0 = BesselK0(2/(3.*r*Sqrt(r)));
    return Log(Atrial*K0) - alpha_1;
  }

  if(r>Lhalfwf) return 0; //r = Lhalfwf;

  return -Xi*(Lwf-2*r)*(Lwf-2*r)/((Lwf-r)*r);
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  DOUBLE K0, K1;
  if(r<Rpar) {
    K0 = BesselK0(2/(3.*r*Sqrt(r)));
    K1 = BesselK1(2/(3.*r*Sqrt(r)));
    return K1/(K0*Sqrt(r)*r*r);
  }
  if(r>Lhalfwf) return 0.; //r = Lhalfwf;

  return Xi*(Lwf-2*r)*Lwf*Lwf/((Lwf-r)*(Lwf-r)*r*r);
}

DOUBLE InterpolateExactE(DOUBLE r) {
  DOUBLE K0, K1, t;
  if(r<Rpar) {
    K0 = BesselK0(2/(3.*r*Sqrt(r)));
    K1 = BesselK1(2/(3.*r*Sqrt(r)));
    t = K1/K0;
    return (t*t-1)/(r*r*r*r*r);
  }
  if(r>Lhalf) return 0.; //r = Lhalf;

  return Xi*Lwf*Lwf*(Lwf*Lwf-3.*Lwf*r+4.*r*r)/((Lwf-r)*(Lwf-r)*(Lwf-r)*r*r*r);
}
#endif

/************************** Construct Grid Power **********************/
void ConstructGridPowerKo(struct Grid *G) {
  double beta,gamma,g;
  AllocateWFGrid(G, 3);

  Message("2D Power Ko(r) trial wavefunction\n");
  Message("  Rpar %" LG "  [Lwf/2] -> %" LG "  [R]\n", Rpar, Rpar*Lhalfwf);
  Rpar *= Lhalfwf;
  beta=Apar;
  gamma=(beta-2)*.5;
  g=1/gamma/exp(gamma*log(Rpar));
  Atrial = 1./BesselK0(g);
  Btrial = 1./BesselK1(g);
  alpha_1 = Atrial/Btrial*(Lwf-2*Rpar)*(Lwf-Rpar)/(Lwf*Lwf)/exp(gamma*log(Rpar));
  Xi = Atrial/Btrial*Rpar*(Lwf-Rpar)*(Lwf-Rpar)/(Lwf*Lwf*(Lwf-2*Rpar))/exp(gamma*log(Rpar));

  G->min = 0;
  G->max = Lhalfwf;
  G->max2 = G->max*G->max;
}

// executed only if TRIAL_POWER_Ko is defined
#ifdef TRIAL_POWER_Ko
DOUBLE InterpolateExactU(DOUBLE r) {
  DOUBLE beta, gamma, g, gmax, K0;
  if(r<Rpar) {
    beta=Apar; gamma=(beta-2)*.5; gmax=690;
    g=1/(gamma*exp(gamma*log(r))); if(g>gmax)g=gmax;
    K0 = BesselK0(g);
    return Log(Atrial*K0) - alpha_1;
  }

  if(r>Lhalfwf) return 0;

  return -Xi*(Lwf-2*r)*(Lwf-2*r)/((Lwf-r)*r);
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  DOUBLE beta, gamma, g, gmax, K0, K1;
  if(r<Rpar) {
    beta=Apar; gamma=(beta-2)*.5; gmax=690;
    g=1/(gamma*exp(gamma*log(r))); if(g>gmax)g=gmax;
    K0 = BesselK0(g); K1 = BesselK1(g);
    return K1/(K0*r)*g*gamma;
  }
  if(r>Lhalfwf) return 0.; //r = Lhalfwf;

  return Xi*(Lwf-2*r)*Lwf*Lwf/((Lwf-r)*(Lwf-r)*r*r);
}

DOUBLE InterpolateExactE(DOUBLE r) {
  DOUBLE beta, gamma, g, gmax, K0, K1, t;
  if(r<Rpar) {
    beta=Apar; gamma=(beta-2)*.5; gmax=690;
    g=1/(gamma*exp(gamma*log(r))); if(g>gmax)g=gmax;
    K0 = BesselK0(g); K1 = BesselK1(g); t = K1/K0;
    return (t*t-1)*exp(-beta*log(r));
  }
  if(r>Lhalf) return 0.; //r = Lhalf;

  return Xi*Lwf*Lwf*(Lwf*Lwf-3*Lwf*r+4*r*r)/((Lwf-r)*(Lwf-r)*(Lwf-r)*r*r*r);
}
#endif

/************************** Find SW Binding Energy ************************/
void FindSWBindingEnergy(void) {
  DOUBLE kappa, k;
  DOUBLE V,y,ymin,ymax,x,xmin,xmax;
  DOUBLE precision = 1e-8;
  int search = ON;
  DOUBLE detuning = 1e-8;
  FILE *out;
  static int first_time = ON;
  DOUBLE energy_shift;
  DOUBLE Ro; // range of the potential

  V = -1;
  Ro = RoSW;
  kappa = sqrt(Vo);

  // x^2 = Eb, - binding energy
  //xmin = detuning*pi/2.*0.1;
  //xmax = detuning*pi/2.*10.;
  xmin = 1e-10;
  xmax = kappa-1e-10;
  ymin = sqrt(-xmin*xmin+kappa*kappa)/tg(Ro*sqrt(-xmin*xmin+kappa*kappa))/xmin - V;
  ymax = sqrt(-xmax*xmax+kappa*kappa)/tg(Ro*sqrt(-xmax*xmax+kappa*kappa))/xmax - V;

  if(ymin*ymax>0) {
    Warning("Find SW Binding Energy: no valid solution found, assuming a SW without a bound state\n");
    energy_shift = 0.;
    return;
  }

  search = ON;
  while(search) {
    x = (xmin+xmax)/2.;
    y = sqrt(-x*x+kappa*kappa)/tg(Ro*sqrt(-x*x+kappa*kappa))/x - V;

    if(y*ymin < 0) {
      xmax = x;
      ymax = y;
    } else {
      xmin = x;
      ymin = y;
    }

    if(fabs(y)<precision) {
      x = (xmax*ymin-xmin*ymax) / (ymin-ymax);
      y = sqrt(-x*x+kappa*kappa)/tg(Ro*sqrt(-x*x+kappa*kappa))/x - V;
      search = OFF;
    }
  }
  k = x;
  if(first_time) Message("  Solution k = %e, error %e (k_approx = %e)\n", k, y, detuning*pi/2.);

  if(first_time) {
    out = fopen("Eb.dat", "w");
    fprintf(out, "%.15e", energy_unit*k*k);
    fclose(out);
  }

  energy_shift = k*k;
  if(first_time) Warning("  Binding energy %e should be subtracted from the total energy !\n", energy_shift);
  first_time = OFF;
}

#ifdef THREE_BODY_HYPERRADIUS_QUADRATIC  // f_3(R): 0, R<R3; 1 - (R - Apar)^2 / (R3 - Apar)^2 , R3<R < Rpar; 1, R>Apar
DOUBLE u3(DOUBLE R) {
  if(R<R3)
    return -100.;
  else if(R<Apar)
    return log(fabs(1. - (R - Apar)*(R - Apar) / ((R3 - Apar)*(R3 - Apar))));
  else
    return 0.;
}
#endif

#ifdef THREE_BODY_ZERO_ENERGY_SOLUTION_3D // f_3(R): 0, R<R3; 1 - (R3/R)^4
DOUBLE u3(DOUBLE R) {
  DOUBLE x;

  x = R3/R;// inverse distance

  if(R<R3)
    return -100.;
  else
    return log(1.-x*x*x*x);
}

DOUBLE Fp3(DOUBLE R) { // f'(R)/f(R) = 4/[r((r/R)^4-1)]
  DOUBLE x;

  if(R<R3)
    return 1.;
  else {
    x = R/R3; // distance
    return 4./(R*(x*x*x*x-1.));
  }
}

DOUBLE E3(DOUBLE R) { // Eloc(R) = f''(R)/f(R) + (f'(R)/f(R))^2 - 5/F f'(R)/f(R)
  DOUBLE x;

  if(R<R3)
    return 0.;
  else {
    x = R/R3; // distance
    x = 4./(R*(x*x*x*x-1.)); // Force
    return x*x;
  }
}
#endif

#ifdef THREE_BODY_ZERO_ENERGY_SOLUTION_2D // f_3(R): 0, R<R3; 1-(R3/R)^2
DOUBLE u3(DOUBLE R) {
  DOUBLE x;

  if(R<R3)
    return -100.;
  else {
    x = R3/R; // inverse distance
    return log(1-x*x);
  }
}

DOUBLE Fp3(DOUBLE R) { // f'(R)/f(R) = 2/(R(x*x-1))
  DOUBLE x;

  if(R<R3)
    return 1.;
  else {
    x = R/R3; // distance
    return 2./(R*(x*x-1.));
  }
}

DOUBLE E3(DOUBLE R) { // Eloc(R) = f''(R)/f(R) + (f'(R)/f(R))^2 + 1/R f'(R)/f(R)
  DOUBLE x;

  if(R<R3)
    return 0;
  else {
    x = R/R3; // distance
    x = 2./(R*(x*x-1.)); // drift force
    return x*x;
  }
}
#endif

#ifdef THREE_BODY_ZERO_ENERGY_SOLUTION_1D // f_3(R): 0, R<R3; log(R/R3)
DOUBLE u3(DOUBLE R) {
  if(R<R3)
    return -100.;
  else
    return log(log(R/R3));
}

DOUBLE Fp3(DOUBLE R) { // f'(R)/f(R) = 1/(R*log(R/R3))
  if(R<R3)
    return 1.;
  else
    return 1./(R*log(R/R3));
}

DOUBLE E3(DOUBLE R) { // Eloc(R) = -f''(R)/f(R) - 1/R f'(R)/f(R) + (f'(R)/f(R))^2
  DOUBLE lnf;

  if(R<R3)
    return 0;
  else {
    lnf = Log(R/R3);
    return 1./(R*R*lnf*lnf);
  }
}
#endif

#ifdef THREE_BODY_ZERO_ENERGY_SOLUTION_1D_SYMMETRIZED // f_3(R): 0, R<R3; ln(R/R3) + ln((L-R)/R3)
DOUBLE u3(DOUBLE R) {
#ifdef HARD_SPHERE_HYPERRADIUS
  if(R<R3) return -100.;
#endif
  if(R<Lhalf)
    return log(fabs(log(R*(L-R)/(R3*(L-R3))))) - log(log(Lhalf*Lhalf/(R3*(L-R3))));
  else
    return 0.;
}

DOUBLE Fp3(DOUBLE R) { // f'(R)/f(R) = 1/(R*log(R/R3))
  DOUBLE x,l;

#ifdef HARD_SPHERE_HYPERRADIUS
  if(R<R3) return 100.;
#endif

  x = R/R3;
  l = (L-R)/R3;

  if(R<Lhalf)
    return (l-x)/(l*R*log(x*l/(x+l-1.)));
  else
    return 0;
}

DOUBLE E3(DOUBLE R) { // Eloc(R) = f''(R)/f(R) + (f'(R)/f(R))^2 + 1/R f'(R)/f(R)
  DOUBLE x,l, lnxl;

#ifdef HARD_SPHERE_HYPERRADIUS
  if(R<R3) return 0.;
#endif

  x = R/R3;
  l = (L-R)/R3;
  lnxl = log(x*l/(x+l-1.));

  if(R<Lhalf)
    return ((l-x)*(l-x) + x*(l+x)*lnxl) / (l*R*lnxl*l*R*lnxl);
  else
    return 0;
}
#endif

#ifdef THREE_BODY_BOUND_STATE_SOLUTION_1D
DOUBLE u3(DOUBLE R) {
  return log(BesselK0(R/R3));
}

DOUBLE Fp3(DOUBLE R) { // f'(R)/f(R)
  return -BesselK1(R / R3) / (R3* BesselK0(R / R3));
}

DOUBLE E3(DOUBLE R) { // Eloc(R) = f''(R)/f(R) + (f'(R)/f(R))^2 + 1/R f'(R)/f(R)
  DOUBLE Fp;
  Fp = -BesselK1(R / R3) / (R3* BesselK0(R / R3));
  return -1. / (R3*R3) + Fp*Fp;
}
#endif


#ifdef THREE_BODY_ZERO_ENERGY_SOLUTION_3D_SYMMETRIZED
DOUBLE u3(DOUBLE R) {
  DOUBLE x, l, l2;

  x = R3/R; // inverse distance
  l = R3/(L-R);
  l2 = R3/Lhalf;

  if(R<R3)
    return -100.;
  else if(R<Lhalf)
    return log(1.-x*x*x*x - l*l*l*l) - log(1.-2.*l2*l2*l2*l2);
  else
    return 0;
}

DOUBLE Fp3(DOUBLE R) { // f'(R)/f(R) = 4*(l^5-x^5) / (l R x^4 + l^5 R (-1 + x^4))
  DOUBLE x,l;

  x = R/R3; // distance
  l = (L-R)/R3;

  if(R<R3)
    return 1.;
  else if(R<Lhalf)
    return 4.*(l*l*l*l*l-x*x*x*x*x)/(-l*R3*x*x*x*x*x + l*l*l*l*l*R3*x*(-1. + x*x*x*x));
  else 
    return 0;
}

DOUBLE E3(DOUBLE R) { // Eloc(R) = f''(R)/f(R) + (f'(R)/f(R))^2 - 5/R f'(R)/f(R)
  DOUBLE x,x2,x4,x5,x6,x9,x10;
  DOUBLE l,l2,l4,l5,l6,l10;
  DOUBLE R2;
  DOUBLE c;

  c = DIMENSION*2.-1.;
  //c = 3.*(DIMENSION-1.)-1.;

  x = R/R3; // distance
  l = (L-R)/R3;

  x2 = x*x;
  x4 = x2*x2;
  x5 = x4*x;
  x6 = x4*x2;
  x9 = x4*x5;
  x10 = x5*x5;

  l2 = l*l;
  l4 = l2*l2;
  l5 = l4*l;
  l6 = l4*l2;
  l10 = l5*l5;

  R2 = R3*R3;

  if(R<R3)
    return 1.;
  else if(R<Lhalf)
    //return (4.*(4.*l10 - 5.*l*x9 - x10 + 5.*l4*x6*(-1. + x4) +l5*x5*(-13. + 5.*x4)))
    //        / (l2*R2*x2*(x4 - l4*(-1. + x4))*(x4 - l4*(-1. + x4)));
  return (4.*((c-1.)*l10 - (c-5.)*l6*(l4-1.)*x4 - (c+8.)*l5*x5 - 5.*l4*x6 + c*l*(l4-1.)*x9 + (5.*l4-1)*x10))/(l2*R*R*(x4 - l4*(x4-1.))*(x4 - l4*(x4-1.)));
  else
    return 0;
}
#endif

/************************** Construct Grid 3Body Zero Energy Phonons *****************/
// zero energy HD: f(r) = const ln r
void ConstructGrid3BodyZeroEnergyPhonons(void) {
  DOUBLE r;

  Message("  3-body term: THREE_BODY_ZERO_ENERGY_SOLUTION_1D_PHONONS trial wavefunction\n");
  Message("     f(r) = A ln(r/R3), r<Rpar; B Exp(- C/r + D/r^2); Rpar - var. parameter\n");

  Warning("  Changing Rpar to units of L/2, Rpar = %lf -> %lf\n", Rpar, Rpar*Lhalfwf);
  Rpar *= Lhalfwf;

  if(Rpar == 0.) {
    Warning("  Rpar -> Lwf/2\n", Lhalfwf);
    Rpar = Lhalfwf;
  }
  if(R3<=0) Error("  parameter 'R3' must be positive, R3 = %lf", R3);
  if(Rpar < R3) Error("  Rpar = %lf < a = %lf, must be otherwise", Rpar, R3);
  if(Rpar > Lhalfwf) Error("  Rpar = %lf > L/2 = %lf, must be otherwise", Rpar, Lhalfwf);

  Ctrial = -Rpar*Rpar / ((Lhalfwf-Rpar)*Log(Rpar/R3));
  r = Lhalfwf;
  Btrial = Ctrial/r*(1. - 0.5*Lhalfwf/r);

  r = Rpar;
  Atrial = Exp(Btrial - Ctrial/r*(1. - 0.5*Lhalfwf/r)) / Log(r/R3);
  Message("    A = %"LE " \n", Atrial);
  Message("    B = %"LE " \n", Btrial);
  Message("    C = %"LE " \n", Ctrial);
}

#ifdef THREE_BODY_ZERO_ENERGY_SOLUTION_1D_PHONONS
DOUBLE u3(DOUBLE r) {
  if(r<R3) {
    return -100;
  }
  else if(r<Rpar) {
    return Log(Atrial*Log(r/R3));
  }
  else if(r<Lhalfwf) {
    return Btrial - Ctrial/r*(1. - 0.5*Lhalfwf/r);
  }
  else {
    return 0;
  }
}

DOUBLE Fp3(DOUBLE r) { // f'(R)/f(R) = 1/(R*log(R/R3))
  if(r<R3) {
    return 100;
  }
  else if(r<Rpar) {
    return 1./(r*Log(r/R3));
  }
  else if(r<Lhalfwf) {
    return Ctrial*(r-Lhalfwf)/(r*r*r);
  }
  else {
    return 0;
  }
}

DOUBLE E3(DOUBLE r) { // E3(R) = -f''(R)/f(R) - 1/R f'(R)/f(R) + (f'(R)/f(R))^2
  DOUBLE logf1;

  if(r<R3) {
    return 100;
  }
  else if(r<Rpar) {
    logf1 = 1./(r*Log(r/R3));
    return logf1*logf1;
  }
  else if(r<Lhalfwf) {
    return Ctrial*(r-Lwf)/(r*r*r*r);
  }
  else {
    return 0;
  }
}
#endif


DOUBLE InterpolateGridMatchedToScatteringSolutionU(struct Grid *Grid, DOUBLE r) {
  if(r>Grid->max-10*Grid->step) // 2D bound state
    return log(K0(Grid->k * r));
  else
#ifdef INTERPOLATE_SPLINE_JASTROW_WF
    return InterpolateSplineU(Grid, r);
#else
    return InterpolateGridU(Grid, r);
#endif
}

DOUBLE InterpolateGridMatchedToScatteringSolutionFp(struct Grid *Grid, DOUBLE r) {
  if(r>Grid->max-10*Grid->step) // 2D bound state
    return -Grid->k*K1(Grid->k * r)/K0(Grid->k * r);
  else
#ifdef INTERPOLATE_SPLINE_JASTROW_WF
    return InterpolateSplineFp(Grid, r);
#else
    return InterpolateGridFp(Grid, r);
#endif
}

DOUBLE InterpolateGridMatchedToScatteringSolutionE(struct Grid *Grid, DOUBLE r) {
  DOUBLE Fp;
  if(r>Grid->max-10*Grid->step) { // 2D bound state
    Fp = K1(Grid->k * r)/K0(Grid->k * r);
    return Grid->k*Grid->k*(-1. + Fp*Fp);
  }
  else
#ifdef INTERPOLATE_SPLINE_JASTROW_WF
    return InterpolateSplineE(Grid, r);
#else
    return InterpolateGridE(Grid, r);
#endif
}

/************************** Construct Grid CSM Plasmon ************************
void ConstructGridCSM3DPlasmon(struct Grid *G, DOUBLE Apar) {
  DOUBLE r;

  AllocateWFGrid(G, 3);

#ifndef TRIAL_3D
  Error("  Trial wavefunction can be used only in 3D!\n");
#endif
  Message("  TRIAL_CSM_3D_PLASMON  Rpar: sqrt(r); exp(-const/r^3 - const / r^2)\n");

  if(Rpar > Lhalf) Error(" Rpar = %lf must be smaller or equal to L/2 = %lf\n", Rpar, Lhalf);

  //Ktrial = 0.5*(sqrt(1. + 4.*D)-1.);
  Ktrial = lambda;
  Message("  CSM lambda is taken from lambda parameter\n");
  Message("  CSM lambda = %lf\n", Ktrial);

  Ctrial = (2.*Ktrial*Rpar*Rpar*Rpar)/(3.*(-Lhalf*sqrt(Lhalf) + Rpar*sqrt(Rpar)));
  Dtrial = -(Ktrial*Lhalf*sqrt(Lhalf)*Rpar*Rpar*Rpar) / (3.*(-Lhalf*sqrt(Lhalf) + Rpar*sqrt(Rpar)));
  
  // f(L/2) = 1
  r = Lhalf;
  Btrial = Ctrial/sqrt(r*r*r) + Dtrial/(r*r*r);

  // f(Rpar)
  r = Rpar;
  Atrial = Btrial - Ctrial/sqrt(r*r*r) - Dtrial/(r*r*r) - Ktrial*log(r);

  G->min = 0;
  G->max = Lwf/2;
  G->max2 = G->max*G->max;
}

#ifdef TRIAL_CSM_3D_PLASMON
DOUBLE InterpolateExactU(DOUBLE r) {
  if(r < Rpar)
    return Atrial + Ktrial*log(r);
  else 
    return Btrial - Ctrial/sqrt(r*r*r) - Dtrial/(r*r*r);
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  if (r < Rpar)
    return Ktrial / r;
  else 
    return (1.5*Ctrial + 3.*Dtrial/(r*sqrt(r))) / (r*r*sqrt(r));
}

DOUBLE InterpolateExactE(DOUBLE r) {
  if(r < Rpar) 
    return -Ktrial/(r*r);
  else 
    return (0.75*Ctrial + 6.*Dtrial/(r*sqrt(r))) / (r*r*r*sqrt(r));
}
#endif*/ 

/************************** Construct Grid CSM Plasmon ************************/
void ConstructGridCSM3DPlasmon(struct Grid *G, DOUBLE Apar) {
  DOUBLE r;

  AllocateWFGrid(G, 3);

#ifndef TRIAL_3D
  Error("  Trial wavefunction can be used only in 3D!\n");
#endif
  Message("  TRIAL_CSM_3D_PLASMON  Rpar, lambda: r^lambda(1 + const r^2); exp(-const/r^3 - const / r^2)\n");

  if(Rpar > Lhalf) Error(" Rpar = %lf must be smaller or equal to L/2 = %lf\n", Rpar, Lhalf);

  //Ktrial = 0.5*(sqrt(1. + 4.*D)-1.);
  Message("  CSM lambda = 0.5*(sqrt(1. + 4.*D)-1.) = %lf\n", 0.5*(sqrt(1. + 4.*D)-1.));

  Ktrial = lambda;
  Message("  CSM lambda is taken from lambda parameter\n");

  Ctrial = (16*(2+Ktrial)*pow(Rpar,3)*(2*pow(Lhalf,1.5)+pow(Rpar,1.5))*(-4*(5+3*Ktrial)*pow(Lhalf,3)+4*pow(Lhalf,1.5)*pow(Rpar,1.5)+pow(Rpar,2)*(7*Rpar+3*Ktrial*Rpar+Sqrt((pow(2*pow(Lhalf,1.5)+pow(Rpar,1.5),2)*(4*(25+12*Ktrial)*pow(Lhalf,3)-4*(35+18*Ktrial)*pow(Lhalf,1.5)*pow(Rpar,1.5)+(49+24*Ktrial)*pow(Rpar,3)))/pow(Rpar,4)))))/(3.*pow(4*pow(Lhalf,3)+4*pow(Lhalf,1.5)*pow(Rpar,1.5)+pow(Rpar,2)*(Rpar+Sqrt((pow(2*pow(Lhalf,1.5)+pow(Rpar,1.5),2)*(4*(25+12*Ktrial)*pow(Lhalf,3)-4*(35+18*Ktrial)*pow(Lhalf,1.5)*pow(Rpar,1.5)+(49+24*Ktrial)*pow(Rpar,3)))/pow(Rpar,4))),2)); 
  Dtrial = (8 * (2 + Ktrial)*pow(Lhalf, 1.5)*pow(Rpar, 3)*(2 * pow(Lhalf, 1.5) + pow(Rpar, 1.5))*(4 * (5 + 3 * Ktrial)*pow(Lhalf, 3) - 4 * pow(Lhalf, 1.5)*pow(Rpar, 1.5) - pow(Rpar, 2)*(7 * Rpar + 3 * Ktrial*Rpar + Sqrt((pow(2 * pow(Lhalf, 1.5) + pow(Rpar, 1.5), 2)*(4 * (25 + 12 * Ktrial)*pow(Lhalf, 3) - 4 * (35 + 18 * Ktrial)*pow(Lhalf, 1.5)*pow(Rpar, 1.5) + (49 + 24 * Ktrial)*pow(Rpar, 3))) / pow(Rpar, 4))))) / (3.*pow(4 * pow(Lhalf, 3) + 4 * pow(Lhalf, 1.5)*pow(Rpar, 1.5) + pow(Rpar, 2)*(Rpar + Sqrt((pow(2 * pow(Lhalf, 1.5) + pow(Rpar, 1.5), 2)*(4 * (25 + 12 * Ktrial)*pow(Lhalf, 3) - 4 * (35 + 18 * Ktrial)*pow(Lhalf, 1.5)*pow(Rpar, 1.5) + (49 + 24 * Ktrial)*pow(Rpar, 3))) / pow(Rpar, 4))), 2));
  Etrial = -(-4 * (5 + 3 * Ktrial)*pow(Lhalf, 3) + 4 * pow(Lhalf, 1.5)*pow(Rpar, 1.5) + pow(Rpar, 2)*(7 * Rpar + 3 * Ktrial*Rpar + Sqrt((pow(2 * pow(Lhalf, 1.5) + pow(Rpar, 1.5), 2)*(4 * (25 + 12 * Ktrial)*pow(Lhalf, 3) - 4 * (35 + 18 * Ktrial)*pow(Lhalf, 1.5)*pow(Rpar, 1.5) + (49 + 24 * Ktrial)*pow(Rpar, 3))) / pow(Rpar, 4)))) / (3.*(2 + Ktrial)*pow(Rpar, 2)*(-4 * pow(Lhalf, 3) + pow(Rpar, 3)));

  // f(L/2) = 1
  r = Lhalf;
  Btrial = Ctrial/sqrt(r*r*r) + Dtrial/(r*r*r);

  // f(Rpar)
  r = Rpar;
  Atrial = Btrial - Ctrial/sqrt(r*r*r) - Dtrial/(r*r*r) - Ktrial*log(r) - log(1. + Etrial*r*r);

  G->min = 0;
  G->max = Lwf/2;
  G->max2 = G->max*G->max;
}

#ifdef TRIAL_CSM_3D_PLASMON
DOUBLE InterpolateExactU(DOUBLE r) {
  if(r < Rpar)
    return Atrial + Ktrial*log(r) + log(1. + Etrial*r*r);
  else 
    return Btrial - Ctrial/sqrt(r*r*r) - Dtrial/(r*r*r);
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  if (r < Rpar)
    return (Ktrial + Etrial*(2.+Ktrial)*r*r) / (r*(1. + Etrial*r*r));
  else 
    return (1.5*Ctrial + 3.*Dtrial/(r*sqrt(r))) / (r*r*sqrt(r));
}

DOUBLE InterpolateExactE(DOUBLE r) {
  if(r < Rpar) 
    return -Ktrial/(r*r) - 2.*Etrial*(3.+Etrial*r*r)/((1. + Etrial*r*r)*(1. + Etrial*r*r));
  else 
    return (0.75*Ctrial + 6.*Dtrial/(r*sqrt(r))) / (r*r*r*sqrt(r));
}
#endif

/************************** Construct Grid Phonon Luttinger *************************/
void ConstructGridPhononLuttingerAuto11(struct Grid *G) {
  DOUBLE k, kmin, kmax;
  DOUBLE y, ymin, ymax;
  DOUBLE precision = 1e-12;
  int search = ON;
  DOUBLE R, Rmin, Rmax;
  DOUBLE a;
  int iter = 1;

  a = aA;
  Message("\n  Jastrow term: ConstructGridPhononLuttingerAuto11 trial wavefunction\n");

  if(a>0) {
    Warning(" assuming repulsive case, with aA = %lf < 0\n", -a);
    a = -a;
  }

  AllocateWFGrid(G, 3);

  Message("  Kpar11= %"LE " \n", Kpar11);
  alpha11_1 = 1. / Kpar11;

  // search for Rpar
  R = Rmin = 1e-8;
  // (2.59) - k ^ 2 = alpha11_1 * (PI / L)*(PI / L)* ((alpha11_1 - 1.)*cotan(PI*R/L)*cotan(PI*R/L)-1.)
  k = kmin = sqrt(-alpha11_1 * (PI / L)*(PI / L)* ((alpha11_1 - 1.) / tan(PI*R / L) / tan(PI*R / L) - 1.));
  // (2.58) -k*tan(k*R + atan(1./(k*a))) = alpha11_1*PI/L*cotan(PI*R/L)
  ymin = -k * tan(k*R + atan(1. / (k*a))) - alpha11_1 * PI / L / tan(PI*R / L);

  R = Rmax = Lhalf;
  k = kmax = sqrt(-alpha11_1 * (PI / L)*(PI / L)* ((alpha11_1 - 1.) / tan(PI*R / L) / tan(PI*R / L) - 1.));
  ymax = -k * tan(k*R + atan(1. / (k*a))) - alpha11_1 * PI / L / tan(PI*R / L);

  if (fabs(ymax) < precision) {
    Message("  no iterations are needed\n");
    R = Rmax;
    search = OFF;
  }
  else {
    if (ymin*ymax > 0) Error("  cannot construct w.f. (11)");
  }

  // Solve k*tan(k*R + atan(1./(k*a))) = -PI*alpha11_1/L/tan(PI*R/L);
  // i.e. f'(R)/f(R)
  while (search) {
    R = (Rmin + Rmax) / 2.;
    k = sqrt(-alpha11_1 * (PI / L)*(PI / L)* ((alpha11_1 - 1.) / tan(PI*R / L) / tan(PI*R / L) - 1.));
    y = -k * tan(k*R + atan(1. / (k*a))) - alpha11_1 * PI / L / tan(PI*R / L);

    if (y*ymin < 0) {
      Rmax = R;
      ymax = y;
    }
    else {
      Rmin = R;
      ymin = y;
    }

    if (fabs(y) < precision) {
      R = (Rmax*ymin - Rmin * ymax) / (ymin - ymax);
      search = OFF;
    }

    if (iter++ == 1000) {
      Warning("  maximal number of iterations exceeded");
    }
  }
  k = sqrt(-alpha11_1 * (PI / L)*(PI / L)* ((alpha11_1 - 1.) / tan(PI*R / L) / tan(PI*R / L) - 1.));
  y = -k * tan(k*R + atan(1. / (k*a))) - alpha11_1 * PI / L / tan(PI*R / L);
  Message("done\n");

  Btrial11 = -Arctg(1 / (k*a)) / k;
  Atrial11 = Sin(PI*R / L);
  Atrial11 = pow(Sin(PI*R / L), alpha11_1) / Cos(k*(R - Btrial11));
  Rpar11 = R;
  sqE11 = k;
  sqE211 = sqE11 * sqE11;

  Message("    k = %"LE " , error %"LE "\n", k, y);
  Message("    R = %"LE " = %lf L/2\n", R, R / Lhalf);
  Message("    A = %"LE " \n", Atrial11);
  Message("    B = %"LE " \n", Btrial11);
  Message("    Equivalent Luttiner parameter K = 1/alpha = %"LE " \n", 1. / alpha11_1);
  //Message("    [Eloc(R-0) -Eloc(R+0)] / Eloc(R) = %"LE "\n", (InterpolateExactE11(G, R*0.99999) - InterpolateExactE11(G, R*1.00001)) / InterpolateExactE11(G, R));

  G->min = 0;
  G->max = L / 2;
  G->max2 = G->max*G->max;
}

#ifdef TRIAL_PHONON_LUTTINGER_AUTO11_SPINFULL
// ideal bosons K -> \infty, alpha11_1 = 1/K = 0
// ideal fermions K = 1, alpha11_1 = 1/K = 1
DOUBLE InterpolateExactU11(DOUBLE x) {
#ifdef BC_ABSENT
  if(x>Lhalf) return 0;
#endif
  if(x<Rpar11)
    return Log(fabs(Atrial11*Cos(sqE11*(x - Btrial11))));
  else
    return Log(pow(Sin(PI*x / L), alpha11_1));
}

DOUBLE InterpolateExactFp11(DOUBLE x) {
#ifdef HARD_SPHERE
#ifdef SECURE
  if(x<aA) Error("InterpolateExactFp11 : (r<aA) call\n");
#endif
#endif
#ifdef BC_ABSENT
  if(x>Lhalf) return 0;
#endif
  if(x<Rpar11)
    return -sqE11 * Tg(sqE11*(x - Btrial11));
  else
    return alpha11_1 * PI / (Lwf*Tg(PI*x / L));
}

DOUBLE InterpolateExactE11(DOUBLE x) {
  DOUBLE c;
#ifdef BC_ABSENT
  if(x>Lhalf) return 0;
#endif
  if(x<Rpar11) {
    c = Tg(sqE11*(x - Btrial11));
    return sqE211 * (1. + c*c);
  }
  else {
    c = 1 / Tg(PI*x / Lwf);
    return alpha11_1 * (PI*PI / (Lwf*Lwf))*(1 + c*c);
  }
}
#endif

/************************** Construct Grid Phonon Luttinger *************************/
void ConstructGridPhononLuttingerAuto12(struct Grid *G) {
  DOUBLE k, kmin, kmax;
  DOUBLE y, ymin, ymax;
  DOUBLE precision = 1e-12;
  int search = ON;
  DOUBLE R, Rmin, Rmax;
  int iter = 1;

  Message("\n  Jastrow term: ConstructGridPhononLuttingerAuto12 trial wavefunction\n");

  if(a > 0) {
    Warning(" assuming repulsive case, with a = %lf < 0\n", -a);
    a = -a;
  }

  AllocateWFGrid(G, 3);

  Message("  Kpar12= %"LE " \n", Kpar12);
  alpha12_1 = 1. / Kpar12;

  // search for Rpar
  R = Rmin = 1e-8;
  // (2.59) - k ^ 2 = alpha12_1 * (PI / L)*(PI / L)* ((alpha12_1 - 1.)*cotan(PI*R/L)*cotan(PI*R/L)-1.)
  k = kmin = sqrt(-alpha12_1 * (PI / L)*(PI / L)* ((alpha12_1 - 1.) / tan(PI*R / L) / tan(PI*R / L) - 1.));
  // (2.58) -k*tan(k*R + atan(1./(k*a))) = alpha12_1*PI/L*cotan(PI*R/L)
  ymin = -k * tan(k*R + atan(1. / (k*a))) - alpha12_1 * PI / L / tan(PI*R / L);

  R = Rmax = Lhalf;
  k = kmax = sqrt(-alpha12_1 * (PI / L)*(PI / L)* ((alpha12_1 - 1.) / tan(PI*R / L) / tan(PI*R / L) - 1.));
  ymax = -k * tan(k*R + atan(1. / (k*a))) - alpha12_1 * PI / L / tan(PI*R / L);

  if (ymin*ymax > 0) Error("  cannot construct w.f. (12)");

  // Solve k*tan(k*R + atan(1./(k*a))) = -PI*alpha12_1/L/tan(PI*R/L);
  // i.e. f'(R)/f(R)
  while (search) {
    R = (Rmin + Rmax) / 2.;
    k = sqrt(-alpha12_1 * (PI / L)*(PI / L)* ((alpha12_1 - 1.) / tan(PI*R / L) / tan(PI*R / L) - 1.));
    y = -k * tan(k*R + atan(1. / (k*a))) - alpha12_1 * PI / L / tan(PI*R / L);

    if (y*ymin < 0) {
      Rmax = R;
      ymax = y;
    }
    else {
      Rmin = R;
      ymin = y;
    }

    if (fabs(y) < precision) {
      R = (Rmax*ymin - Rmin * ymax) / (ymin - ymax);
      k = sqrt(-alpha12_1 * (PI / L)*(PI / L)* ((alpha12_1 - 1.) / tan(PI*R / L) / tan(PI*R / L) - 1.));
      y = -k * tan(k*R + atan(1. / (k*a))) - alpha12_1 * PI / L / tan(PI*R / L);
      search = OFF;
    }

    if (iter++ == 1000) {
      Warning("  number of iterations exceeded");
    }
  }
  Message("done\n");

  Btrial12 = -Arctg(1 / (k*a)) / k;
  Atrial12 = Sin(PI*R / L);
  Atrial12 = pow(Sin(PI*R / L), alpha12_1) / Cos(k*(R - Btrial12));
  Rpar12 = R;
  sqE12 = k;
  sqE212 = sqE12 * sqE12;

  Message("    k = %"LE " , error %"LE "\n", k, y);
  Message("    R = %"LE " = %lf L/2\n", R, R / Lhalf);
  Message("    A = %"LE " \n", Atrial12);
  Message("    B = %"LE " \n", Btrial12);
  Message("    Equivalent Luttiner parameter K = 1/alpha = %"LE " \n", 1. / alpha12_1);
  //Message("    [Eloc(R-0) -Eloc(R+0)] / Eloc(R) = %"LE "\n", (InterpolateExactE12(G, R*0.99999) - InterpolateExactE12(G, R*1.00001)) / InterpolateExactE12(G, R));

  G->min = 0;
  G->max = L / 2;
  G->max2 = G->max*G->max;
}

/* executed only if TRIAL_PHONON is defined */
#ifdef TRIAL_PHONON_LUTTINGER_AUTO12_SPINFULL
// ideal bosons K -> \infty, alpha12_1 = 1/K = 0
// ideal fermions K = 1, alpha12_1 = 1/K = 1
//DOUBLE InterpolateExactU12(DOUBLE x) {
DOUBLE InterpolateExactU(DOUBLE x) {
  if (x < Rpar12)
    return Log(fabs(Atrial12*Cos(sqE12*(x - Btrial12))));
  else
    return Log(pow(Sin(PI*x / L), alpha12_1));
}

//DOUBLE InterpolateExactFp12(DOUBLE x) {
DOUBLE InterpolateExactFp(DOUBLE x) {
#ifdef HARD_SPHERE
#ifdef SECURE
  if (x < a) Error("InterpolateExactFp12 : (r<a) call\n");
#endif
#endif
  if (x < Rpar12)
    return -sqE12 * Tg(sqE12*(x - Btrial12));
  else
    return alpha12_1 * PI / (Lwf*Tg(PI*x / L));
}

//DOUBLE InterpolateExactE12(DOUBLE x) {
DOUBLE InterpolateExactE(DOUBLE x) {
  DOUBLE c;

  if (x < Rpar12) {
    c = Tg(sqE12*(x - Btrial12));
    return sqE212 * (1. + c * c);
  }
  else {
    c = 1. / Tg(PI*x / Lwf);
    return alpha12_1 * (PI*PI / (Lwf*Lwf))*(1. + c * c);
    //return -alpha12_1*(PI*PI/(L*L))*((alpha12_1-1)*c*c-1); // Local energy
  }
}
#endif

#ifdef TRIAL_ZERO_RANGE_K0_TRAP
DOUBLE InterpolateExactU(DOUBLE r) {
  return Log(fabs(BesselK0(Atrial*r)));
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  return -Atrial*BesselK1(Atrial*r)/BesselK0(Atrial*r);
}

DOUBLE InterpolateExactE(DOUBLE r) {
  DOUBLE fp;
  fp = InterpolateExactFp(r);
  return -Atrial*Atrial + fp*fp;
}
#endif

#ifdef TRIAL_ZERO_RANGE_LOG_TRAP
DOUBLE InterpolateExactU(DOUBLE r) {
  //if(r<a)
  return Log(fabs(Log(r/a)));
  //else
  //  return -100*r;
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  return 1./(r*Log(r/a));
}

DOUBLE InterpolateExactE(DOUBLE r) {
  DOUBLE logf1;
  logf1 = Log(r/a);
  return 1. / (r*r*logf1*logf1);
}
#endif

#ifdef TRIAL_ZERO_RANGE_LOG_EXP_TAIL_TRAP
DOUBLE InterpolateExactU(DOUBLE r) {
  return Log(fabs(Log(a/r))) - r*r/(b*(Apar+r));
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  return -1./(r*Log(a/r)) - (r*(2.*Apar+r)/(b*(Apar+r)*(Apar+r)));
}

DOUBLE InterpolateExactE(DOUBLE r) {
  DOUBLE logf1;
  logf1 = Log(a/r);
  return 1. / (r*r*logf1*logf1) + (4.*Apar*Apar + 3.*Apar*r + r*r) / (b*(Apar+r)*(Apar+r)*(Apar+r));
}
#endif

#ifdef TRIAL_ZERO_RANGE_K0_EXP_TAIL_TRAP
DOUBLE InterpolateExactU(DOUBLE r) {
  return Log(fabs(BesselK0(Atrial*r))) - r*r/(b*(Apar+r));
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  return -Atrial*BesselK1(Atrial*r)/BesselK0(Atrial*r) - (r*(2.*Apar+r)/(b*(Apar+r)*(Apar+r)));
}

DOUBLE InterpolateExactE(DOUBLE r) {
  return (2*Apar*Apar + Apar*(Apar + r) + (Apar + r)*(Apar + r) - Atrial*Atrial*b*(Apar + r)*(Apar + r)*(Apar + r))/(b*(Apar + r)*(Apar + r)*(Apar + r)) + (Atrial*Atrial*BesselK1(Atrial*r)*BesselK1(Atrial*r))/(BesselK0(Atrial*r)*BesselK0(Atrial*r));
}
#endif

#ifdef TRIAL_ZERO_RANGE_LOG_NO_NODE_TRAP
DOUBLE InterpolateExactU(DOUBLE r) {
  return log((log(1. + a/r))*exp(-r/Rpar)*(1. + r/Rpar));
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  return -(r/(Rpar*(r + Rpar))) - a/(r*(a + r)*log((a + r)/r));
}

DOUBLE InterpolateExactE(DOUBLE r) {
  return (r + 2.*Rpar)/(Rpar*(r + Rpar)*(r + Rpar)) +(a*(a - r*log((a + r)/r)))/(r*r*(a + r)*(a + r)*log((a + r)/r)*log((a + r)/r));
}
#endif

#ifdef TRIAL_ZERO_RANGE_LOG_MATCH_TRAP
DOUBLE InterpolateExactU(DOUBLE r) {
  if(r<Rpar)
    return Log(Log(a/r));
  else
    return Log(Log(a/r)) - (r/Rpar-1.) / Log(r/a);
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  DOUBLE logf;
  logf = Log(r/a);
  if(r<Rpar)
    return 1./(r*logf);
  else
    return (r + Rpar + (Rpar-r)/logf)/(r*Rpar*logf);
}

DOUBLE InterpolateExactE(DOUBLE r) {
  DOUBLE logf;
  logf = Log(r/a);
  if(r<Rpar)
    return 1. / (r*r*logf*logf);
  else
    return (2.*(Rpar-r) + logf*(2.*r + Rpar - r*logf))/(r*r*Rpar*logf*logf*logf);
}
#endif

/************************** Construct Grid Coulomb 3D Phonon **************/
void ConstructGridCoulombPhonon(struct Grid *G) {
  DOUBLE arg,r,I0,I1,I2;
  AllocateWFGrid(G, 3);

  Message("  Jastrow term 11,22: 3D repulsive Coulomb + phonons\n");
  Warning("  Changing Rpar %lf -> %lf [L/2] units\n", Rpar, Rpar*Lhalf);
  Rpar *= Lhalf;

  // matching f'(Rpar)/f(Rpar)
  //Btrial = -(-Rpar*Rpar*(Rpar-L)*(Rpar-L)*(Rpar-L)*(Sqrt(Rpar)*BesselI0(2*Sqrt(Rpar)) - BesselI1(2*Sqrt(Rpar)) + Sqrt(Rpar)*BesselIN(2,2*Sqrt(Rpar))))/(4.*(L*L*L - 3.*L*L*Rpar + 3.*L*Rpar*Rpar - 2.*Rpar*Rpar*Rpar)*BesselI1(2*Sqrt(Rpar)));
  r = Rpar;
  arg = 2.*Sqrt(r);
  I0 = BesselI0(arg);
  I1 = BesselI1(arg);
  I2 = BesselIN(2,arg);
  Btrial = (Sqrt(r)*(I0+I2)-I1)/(2.*r*I1)/(2.*(1./(r*r*r)-1./((L-r)*(L-r)*(L-r))));
  // matching f(Rpar)
  Atrial = Btrial*(2./Lhalf2-1./(Rpar*Rpar)-1./((L-Rpar)*(L-Rpar))) - Log(BesselI1(2*Sqrt(Rpar)))+0.5*Log(Rpar);

  Message("  Atrial = %"LG" , Btrial = %"LG"\n", Atrial, Btrial);

  if(Rpar>Lhalf) Error("  Rpar must be smaller than L/2!\n");

  G->max = Lhalf;
  G->max2 = G->max*G->max;
}

#ifdef TRIAL_COULOMB_3D_PHONON
DOUBLE InterpolateExactU(DOUBLE r) {
  if(r>Lhalf)
    return 0;
  else if(r<Rpar)
    return Atrial+Log(BesselI1(2*Sqrt(r)))-0.5*Log(r);
  else
    return Btrial*(2./Lhalf2-1./(r*r)-1./((L-r)*(L-r)));
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  DOUBLE I0,I1,I2,arg;
  if(r>Lhalf) 
    return 0.;
  else if(r<Rpar) {
    arg = 2.*Sqrt(r);
    I0 = BesselI0(arg);
    I1 = BesselI1(arg);
    I2 = BesselIN(2,arg);
    return (Sqrt(r)*(I0+I2)-I1)/(2.*r*I1);
  }
  else {
    return 2.*Btrial*(1./(r*r*r)-1./((L-r)*(L-r)*(L-r)));
  }
}

DOUBLE InterpolateExactE(DOUBLE r) {
  //3D Eloc = [-f"/f - 2/r f'/f] + (f'/f)^2
  DOUBLE c;
  if(r>Lhalf)
    return 0.;
  else if(r<Rpar) {
    c = InterpolateExactFp(r);
    return -1/r + c*c;
  }
  else {
    return 2.*Btrial*(L*L*L*L - 4.*L*L*L*r + 6.*L*L*r*r - 2.*L*r*r*r + 2.*r*r*r*r)/(r*r*r*r*(L-r)*(L-r)*(L-r)*(L-r));
  }
}
#endif

 /************************** Construct Grid Mc Guire PBC *******************/
void ConstructGridMcGuirePhonon(void) {
  DOUBLE x;

  Message("  ConstructGridMcGuirePhonon12: finding adjusted value of b\n");
  Message("  f12: ((1D) Exp(x/b); (x<Rpar12), sin^a(pi*x/L) otherwise\n");

  if(Rpar12 > Lhalf) {
    Warning("  Rpar12 = %"LG" > L/2 = %"LG", setting Rpar12 = L/2\n", Rpar12, Lhalf);
    Rpar12 = Lhalf;
  }

  //Atrial12 = -(L*tg(PI*Rpar12/L))/b/PI;
  x = Rpar12;
  Atrial12 = -1./b/(PI/(L*tg(PI*x/L)));
  Btrial12 = Atrial12*Log(Sin(PI*x/L)) + x/b;

  Message("done\n");
}

#ifdef TRIAL_MCGUIRE_PHONON // (1D) Exp(x/b); (x<Rpar12), sin^a(pi*x/L) otherwise
DOUBLE InterpolateExactU(DOUBLE x) {
  if(x<Rpar12)
    return -x/b + Btrial12;
  else
    //return Log(pow(Sin(PI*x/L), Atrial12));
    return Atrial12*Log(Sin(PI*x/L));
}

DOUBLE InterpolateExactFp(DOUBLE x) {
  if(x<Rpar12)
    return -1./b;
  else
    return Atrial12*PI/(L*tg(PI*x/L));
}

// Eloc1D = -[f"/f - (f'/f)^2]
DOUBLE InterpolateExactE(DOUBLE x) {
  DOUBLE c;

  if(x<Rpar12) {
    return 0;
  }
  else {
    c = 1. / tg(PI*x/L);
    return Atrial12*(PI*PI/(L*L))*(1.+c*c);
  }
}
#endif

#ifdef TRIAL_COS
DOUBLE InterpolateExactU(DOUBLE x) {
  return Log(fabs(Cos(sqE*x)));
}

DOUBLE InterpolateExactFp(DOUBLE x) {
  return -sqE * tan(sqE*x);
}

DOUBLE InterpolateExactE(DOUBLE x) {
  DOUBLE c;
  c = Tg(sqE*x);
  return sqE*sqE*(1.+c*c);
}
#endif

#ifdef TRIAL_SIN_EXCITED_BOSONS
int InterpolateExactSign(DOUBLE x) {
  if(Sin(sqE*x+Bpar)>0.)
    return 1;
  else
    return -1;
}

DOUBLE InterpolateExactU(DOUBLE x) {
  return (1./Kpar)*Log(fabs(Sin(sqE*x+Bpar)));
}

DOUBLE InterpolateExactFp(DOUBLE x) {
  return (1./Kpar)*sqE / Tg(sqE*x+Bpar);
}

DOUBLE InterpolateExactE(DOUBLE x) {
  DOUBLE c;
  c = 1./Tg(sqE*x+Bpar);
  return (1./Kpar)*sqE*sqE*(1.+c*c);
}
#endif

#ifdef TRIAL_COS_OFFSET
DOUBLE InterpolateExactU(DOUBLE x) {
  return Log(fabs(Apar+Bpar*Cos(sqE*x)));
}

DOUBLE InterpolateExactFp(DOUBLE x) {
  return -sqE*Bpar*Sin(sqE*x) /(Apar+Bpar*Cos(sqE*x));
}

DOUBLE InterpolateExactE(DOUBLE x) {
  DOUBLE f, fp, fpp;
  f = Apar+Bpar*Cos(sqE*x);
  fp = -sqE*Bpar*Sin(sqE*x);
  fpp = -sqE*sqE*Bpar*Cos(sqE*x);
  return -fpp/f + (fp / f)*(fp / f);
}
#endif

/************************** Construct Grid Square well 2D Zero Energy *************************/
// zero energy HD: f(r) = const ln r
void ConstructGridSoftDisks2DZeroEnergyPhonons(struct Grid *G) {
  DOUBLE r;

  Message("  Jastrow term: TRIAL_SOFT_DISKS_2D_ZERO_ENERGY_PHONONS trial wavefunction\n");
  Message("     f(r) = Io(kappa r), r<Ro; A ln(r/a), r<Rpar; B Exp(-C/r + D/r^2); Rpar - var. parameter\n");
  Warning("  Changing Rpar to units of L/2, Rpar = %lf -> %lf\n", Rpar, Rpar*Lhalfwf);
  Rpar *= Lhalfwf;

  if(Rpar == 0.) {
    Warning("  Rpar -> Lwf/2\n", Lhalfwf);
    Rpar = Lhalfwf;
  }
  if(a<=0) Error("  parameter 'a' must be positive, a = %lf", a);
  if(Rpar > Lhalfwf) Error("  Rpar = %lf > L/2 = %lf, must be otherwise", Rpar, Lhalfwf);
  if(Rpar < Ro) Error("  Rpar = %lf < Ro = %lf, must be otherwise", Rpar, Ro);

  AllocateWFGrid(G, 3);

  A4SD = -Rpar*Rpar / ((Lhalfwf-Rpar)*Log(Rpar/a));
  r = Lhalfwf;
  A3SD = A4SD/r*(1. - 0.5*Lhalfwf/r);

  r = Rpar;
  A2SD = Exp(A3SD - A4SD/r*(1. - 0.5*Lhalfwf/r)) / Log(r/a);

  if(fabs(trial_Vo)<1e-6) Warning("Vo=0\n");

  trial_Kappa = sqrt(2.*m_mu*trial_Vo);
  trial_Kappa2 = 2.*m_mu*trial_Vo;
  r = Ro;
  A1SD = A2SD*Log(r/a)/BesselI0(trial_Kappa*r);

  Message("    A1 = %"LE " \n", A1SD);
  Message("    A2 = %"LE " \n", A2SD);
  Message("    A3 = %"LE " \n", A3SD);
  Message("    A4 = %"LE " \n", A4SD);

  G->min = 0;
  G->max = Lhalfwf;
  G->max2 = G->max*G->max;
}

#ifdef TRIAL_SOFT_DISKS_2D_ZERO_ENERGY_PHONONS
DOUBLE InterpolateExactU(DOUBLE r) {
  if(r<Ro) {
    return Log(fabs(A1SD*BesselI0(trial_Kappa*r)));
  }
  else if(r<Rpar) {
    return Log(fabs(A2SD*Log(r/a)));
  }
  else if(r<Lhalfwf) {
    return A3SD - A4SD/r*(1. - 0.5*Lhalfwf/r);
  }
  else {
    return 0;
  }
}

DOUBLE InterpolateExactFp(DOUBLE r) {
  if(r<Ro) {
    return trial_Kappa*BesselI1(trial_Kappa*r)/BesselI0(trial_Kappa*r);
  }
  else if(r<Rpar) {
    return 1./(r*Log(r/a));
  }
  else if(r<Lhalfwf) {
    return A4SD*(r-Lhalfwf)/(r*r*r);
  }
  else {
    return 0;
  }
}

DOUBLE InterpolateExactE(DOUBLE r) {
// 2D Eloc = [-(f"/f + f'/r)/f] + (f'/f)^2
  DOUBLE logf1,fp;

  if(r<Ro) {
    fp = InterpolateExactFp(r);
    return -trial_Kappa2 + fp*fp;
  }
  else if(r<Rpar) {
    logf1 = Log(r/a);
    return 1. / (r*r*logf1*logf1);
  }
  else if(r<Lhalfwf) {
    return A4SD*(r-Lwf)/(r*r*r*r);
  }
  else {
    return 0;
  }
}
#endif

/************************** Construct Soft Sphere Height ************************/
DOUBLE ConstructSoftSphereHeight(DOUBLE Ro, DOUBLE a) {
  DOUBLE trial_Vo;

#ifdef TRIAL_3D
  DOUBLE x, xmin, xmax;
  DOUBLE y, ymin, ymax;
  DOUBLE V;
  DOUBLE precision = 1e-12;
  DOUBLE k, kappa, K, K2;
  DOUBLE delta, delta_bar;
  int search = ON;
  int repeat = ON;
  DOUBLE reff,R;

  Message("Jastrow term: solution of two-body scattering problem\n");
  Message("  parameters: sphere size Ro = %"LE "\n", Ro);
  Message("  parameters: a = %lf scattering length\n", a);
  if(Ro<a) Error("Size of the soft sphere %"LE " is smaller than the scattering length %"LE " \n", Ro, a);

  // define kappa
  // Solve th x / x = 1 - a/R, E<V_0
  //       tg x / x = 1 - a/R, E>V_0
  //       x = kappa R
  V = 1. - a/Ro;
  xmin = 1e-5;
  xmax = 100;
  ymin = tanh(xmin)/xmin - V;
  ymax = tanh(xmax)/xmax - V;
  if(ymin*ymax > 0) Error("Ro = %%"LF" Can't construct the grid : No solution found", Ro);
  while(search) {
    x = (xmin+xmax)/2.;
    y = tanh(x)/x - V;
    if(y*ymin < 0) {
      xmax = x;
      ymax = y;
    } else {
      xmin = x;
      ymin = y;
    }
    if(fabs(y)<precision) {
      x = (xmax*ymin-xmin*ymax) / (ymin-ymax);
      y = tanh(x)/x - V;
      search = OFF;
    }
  }

  kappa = x / Ro;
  trial_Vo = kappa*kappa;
  Message("  Soft-sphere height Vo = %lf\n", trial_Vo);

  // effective range
  K = kappa;
  R = Ro;
  reff = (2.*((R*(3.*a*a - 3.*a*R + R*R))/3. - ((a - R)*(a - R)*(-2.*K*R + sinh(2.*K*R)))/(4.*K*sinh(K*R)*sinh(K*R))))/(a*a);

  Message("  effective range %lf\n", reff);
#endif
#ifdef TRIAL_2D // 2D system
  DOUBLE V,y,ymin,ymax,x,xmin,xmax;
  DOUBLE precision = 1e-8;
  int search = ON;
  DOUBLE kappa;

  //calculate strength of the potential
  // define kappa
  // x = kappa R
  Message("\n  Constructing soft disk potential with Ro = %"LF" and aA = %"LE"\n", Ro, a);

  // soft disk:   Solve Jo(x) / (x J1(x)) = ln(a/R)
  // square well: Solve Io(x) / (x I1(x)) = ln(R/a)
  V = Log(Ro/a);

  xmin = 100./V; // Io(x) / (x I1(x)) = 1/x, x -> inf
  xmax = 0.01*sqrt(2./V); //  Io(x) / (x I1(x)) = 2/x^2, x -> 0
  ymin = BesselI0(xmin)/(xmin*BesselI1(xmin)) - V;
  ymax = BesselI0(xmax)/(xmax*BesselI1(xmax)) - V;
  if(ymin*ymax > 0) Error("Ro = %"LF" Can't construct the grid : No solution found", Ro);
  while(search) {
    x = (xmin+xmax)/2.;
    y = BesselI0(x)/(x*BesselI1(x)) - V;
    if(y*ymin < 0) {
      xmax = x;
      ymax = y;
    } else {
      xmin = x;
      ymin = y;
    }
    if(fabs(y)<precision) {
      x = (xmax*ymin-xmin*ymax) / (ymin-ymax);
      y = BesselI0(x)/(x*BesselI1(x)) - V;
      search = OFF;
    }
    //Message("x = %"LG" y = %"LG"\n", x, y);
  }
  Message("  Solution is kappa R = %"LF", error %"LE"\n", x, y);
  kappa = x / Ro;
  trial_Vo = kappa*kappa;
  Message("  Soft-disk height Vo = %lf\n", trial_Vo);

  Message("  Check a = %"LE"\n", Ro*Exp(-BesselI0(kappa*Ro)/(kappa*Ro*BesselI1(kappa*Ro))));

  trial_Vo = kappa*kappa;
#endif

#ifdef TRIAL_1D // 1D system
  DOUBLE x, xmin, xmax;
  DOUBLE y, ymin, ymax;
  DOUBLE V;
  DOUBLE precision = 1e-12;
  DOUBLE k, kappa, K, K2;
  DOUBLE delta, delta_bar;
  int search = ON;
  int repeat = ON;
  DOUBLE reff,R,a_check;

  Message("Jastrow term: solution of two-body scattering problem\n");
  Message("  parameters: sphere rod diameter Ro = %"LE "\n", Ro);
  Message("  parameters: a = %lf scattering length\n", a);
  if(Ro<a) Error("Size of the soft sphere %"LE " is smaller than the scattering length %"LE " \n", Ro, a);

  // define kappa
  // Solve th x / x = 1 - a/R, E<V_0 (3D)
  //       tg x / x = 1 - a/R, E>V_0 (3D)
  //       1 / (x th x) = 1 - a/R (1D)
  //       x = kappa R
  V = 1. - a/Ro;
  xmin = 1e-5;
  xmax = 100;
  ymin = 1./(xmin*tanh(xmin)) - V;
  ymax = 1./(xmax*tanh(xmax)) - V;
  if(ymin*ymax > 0) Error("Ro = %%"LF" Can't construct the grid : No solution found", Ro);
  while(search) {
    x = (xmin+xmax)/2.;
    y = 1./(x*tanh(x)) - V;
    if(y*ymin < 0) {
      xmax = x;
      ymax = y;
    } else {
      xmin = x;
      ymin = y;
    }
    if(fabs(y)<precision) {
      x = (xmax*ymin-xmin*ymax) / (ymin-ymax);
      y = 1./(x*tanh(x)) - V;
      search = OFF;
    }
  }

  kappa = x / Ro;
  trial_Vo = kappa*kappa;
  Message("  Soft-sphere height Vo = %lf\n", trial_Vo);

  R = Ro;
  // scattering length (consistency check)
  a_check = R*(1.-1./(kappa*R*tanh(kappa*R)));
  Message("  Obtained scattering length %e\n", a_check);
  // effective range
  if(fabs(a)>1e-3) {
    reff = R - R*R*R/(3.*a*a) + 1./(kappa*kappa*a);
    Message("  effective range %lf\n", reff);
  }
  else { // a = 0,
    reff = - R*R*R/3.;
    Message("  effective range r*a*a = %lf\n", reff);
  }
#endif

  Vo = trial_Vo;
  return trial_Vo;
}

DOUBLE ConstructHardSphereSquareWellDepth(DOUBLE R1, DOUBLE R2, DOUBLE a3D) {
  DOUBLE trial_Vo;

  DOUBLE x, xmin, xmax;
  DOUBLE y, ymin, ymax;
  DOUBLE V;
  DOUBLE precision = 1e-12;
  DOUBLE k, kappa, K, K2;
  DOUBLE delta, delta_bar;
  int search = ON;
  int repeat = ON;

  Message("Jastrow term: solution of two-body scattering problem\n");
  Message("  parameters: hard sphere size a = %"LE "\n", R1);
  Message("  parameters: square well size Ro = %"LE "\n", R2);
  Message("  parameters: Apar = %lf scattering length\n", a3D);

  // define k
  // Solve: tg [x(R2-R1)] / x = R2 - a3D
  V = R2-a3D;
  xmin = 1e-5;
  xmax = 0.5*3.14/(R2-R1);
  ymin = tan(xmin*(R2-R1))/xmin - V;
  ymax = tan(xmax*(R2-R1))/xmax - V;
  if(ymin*ymax > 0) Error("Ro = %%"LF" Can't construct the grid : No solution found", Ro);
  while(search) {
    x = (xmin+xmax)/2.;
    y = tan(x*(R2-R1))/x - V;
    if(y*ymin < 0) {
      xmax = x;
      ymax = y;
    } else {
      xmin = x;
      ymin = y;
    }
    if(fabs(y)<precision) {
      x = (xmax*ymin-xmin*ymax) / (ymin-ymax);
      y = tan(x*(R2-R1))/x - V;
      search = OFF;
    }
  }

  trial_Vo = x*x;
  Message("  Square well depth Vo = %lf\n", trial_Vo);

  Vo = trial_Vo;
  return trial_Vo;
}

/************************** Construct Grid Square well 2D Zero Energy *************************/
// zero energy HD: f(r) = const ln r
void ConstructSS_Simple(struct Grid *G) {
  DOUBLE x, xmin, xmax;
  DOUBLE y, ymin, ymax;
  DOUBLE V;
  DOUBLE precision = 1e-12;
  DOUBLE k, kappa, K, K2;
  DOUBLE delta, delta_bar;
  int search = ON;
  int repeat = ON;

  Message("  Jastrow term: TRIAL_SS_1D trial wavefunction\n");
  Message("    exp(k*x) + Apar*exp(-k*x) x<Ro, sin(pi*x/L)^(1/Kpar) x>Ro, var. parameters Apar, Kpar\n");

  AllocateWFGrid(G, 3);

  sqE = PI/L;
  sqE2 = sqE*sqE;
  alpha_1 = 1./Kpar;
  Message("    power : 1/Kpar = %"LE " \n", alpha_1);

  // match f'(x) at Ro
  // x*tanh(x*Ro) == alpha_1*sqE/Tg(sqE*x);
  V = alpha_1*sqE/(Tg(sqE*Ro));
  xmin = 1e-5;
  xmax = 1e5;
  x = xmin;
  ymin = x*tanh(x*Ro) - V;
  x = xmax;
  ymax = x*tanh(x*Ro) - V;
  if(ymin*ymax > 0) Error("Can't construct the grid : No solution found (ymin =%lf, ymax = %lf)", ymin, ymax);
  while(search) {
    x = (xmin+xmax)/2.;
    y = x*tanh(x*Ro) - V;
    if(y*ymin < 0) {
      xmax = x;
      ymax = y;
    } else {
      xmin = x;
      ymin = y;
    }
    if(fabs(y)<precision) {
      x = (xmax*ymin-xmin*ymax) / (ymin-ymax);
      y = x*tanh(x*Ro) - V;
      search = OFF;
    }
  }

  Btrial = x;
  Message("    characteristic momentum Btrial = %lf\n", Btrial);

  x = Ro;
  Atrial = Log(pow(Sin(sqE*x), alpha_1)) - Log(cosh(Btrial*x));

  G->min = 0;
  G->max = Lhalfwf;
  G->max2 = G->max*G->max;
}

#ifdef TRIAL_SS_1D
DOUBLE InterpolateExactU(DOUBLE x) {
  if(x<Ro)
    return Atrial  + Log(cosh(Btrial*x));
  else if(x<Lhalf)
    return Log(pow(Sin(sqE*x), alpha_1));
  else
    return 0.;
}

DOUBLE InterpolateExactFp(DOUBLE x) {
  if(x<Ro)
    return Btrial*tanh(Btrial*x);
  else if(x<Lhalf)
    return alpha_1*sqE/Tg(sqE*x);
  else
    return 0.;
  }

DOUBLE InterpolateExactE(DOUBLE x) {
  // Eloc = -f"/f + (f'/f)^2
  DOUBLE c;

  if(x<Ro) {
    c = 1. / cosh(Btrial*x);
    return -Btrial*Btrial*c*c;
  }
  else if(x<Lhalf){
    c = 1. / Sin(sqE*x);
    return alpha_1*sqE2*c*c;
  }
  else
    return 0.;
}

/*DOUBLE InterpolateExactU(DOUBLE x) {
#ifdef HARD_SPHERE
  if(x<a) return -100;
#endif
  return Log(fabs(Sin(sqE*(x-a))));
}

DOUBLE InterpolateExactFp(DOUBLE x) {
  return sqE / Tg(sqE*(x-a));
}

// Eloc1D = -[f"/f - (f'/f)^2]
DOUBLE InterpolateExactE(DOUBLE x) {
  DOUBLE g;

  g = sqE/Tg(sqE*(x-a));
  //g = 0.; //Real local energy
  return sqE2 + g*g;
}*/

#endif