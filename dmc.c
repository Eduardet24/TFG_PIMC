/*DMC.c*/
#include <stdio.h>
#include <limits.h>
#include "main.h"
#include "dmc.h"
#include "randnorm.h"
#include "utils.h"
#include "vmc.h"
#include "trial.h"
#include "rw.h"
#include "display.h"
#include "crystal.h"
#include "spline.h"
#include "speckles.h"
#include "compatab.h"
#include "memory.h"
#include "mymath.h"
#include MATHINCLUDE

#ifdef HARD_SPHERE
#  define DMC_OVERLAP_RETHROW
#  define Overlap_check_vector Overlapping
#endif
#ifdef HARD_SPHERE_HYPERRADIUS
#  define DMC_OVERLAP_RETHROW
#  define Overlap_check_vector CheckVectorHyperradiusOverlapping
#endif

/******************** Coordinates Distance ***********************************/
double CoordinatesDistanceRpRpp(int w) {
  DOUBLE dr[3],r2=0.;
  int i;

  for(i=0; i<N; i++) {
    CaseX(dr[0] = W[w].Rpp[i][0] - W[w].Rp[i][0]);
    CaseY(dr[1] = W[w].Rpp[i][1] - W[w].Rp[i][1]);
    CaseZ(dr[2] = W[w].Rpp[i][2] - W[w].Rp[i][2]);
    r2 += FindNearestImage(&dr[0], &dr[1], &dr[2]);
  }
  return r2;
}

/**************************** DMC Quadratic Move *****************************/
void DMCQuadraticMove(int w) {
  DOUBLE dt_sqrt = Sqrt(dt);
#ifndef MEMORY_CONTIGUOUS
  int i,j;
#endif
#ifdef SPINFULL_TUNNELING
  DOUBLE flip_probability, x;
  int i_spin;
#endif

#ifdef DMC_OVERLAP_RETHROW
 int go = ON;
 int copy_initial_coordinates = ON; // in case of repeating Gaussian jump, no need to copy Rppp to R

  while(go == ON) {
    tried++;
#endif

#ifdef MEASURE_CORRELATION_FUNCTIONS_IN_THE_MIDDLE_OF_DRIFT
#ifdef DMC_OVERLAP_RETHROW
  if(copy_initial_coordinates) {
    copy_initial_coordinates = OFF;
#endif
  if(W[w].status_end_of_step_initialized) { // restore positions from the end of time step, change the code in virtual walkers!
    CopyVectorToWalker(&W[w], W[w].R_dmc_end_of_step, W[w].spin); // copy the position in the middle of the step
  }
  else { // when not initialized
    W[w].status_end_of_step_initialized = 1; // next time will be initialized
  }
#ifdef DMC_OVERLAP_RETHROW
  }
#endif
#endif

  // R = R(i-1) + chi, f(chi) = const Exp(-R^2/2dt)
  GaussianJump(w, dt_sqrt); // sqrt(dt)

#ifdef SECURE
#ifdef HARD_SPHERE
  if(Overlapping(W[w].R)) Error("Gaussian jump : impossible situation: (1)");
#endif
#endif
#ifdef SCALABLE
  UpdateScalableTable(&W[w]);
#endif

  // R' = R + F(R) dt/2
  DriftJump(0.5*dt, w);

#ifdef DMC_OVERLAP_RETHROW
  if(Overlap_check_vector(W[w].Rp)) {
    overlaped++;
#ifdef HARD_SPHERE_KILL_OVERLAPPED
    go = OFF; W[w].status = DEAD;
#endif
    if(verbosity) Message("DMC move: particles overlap (Rp)\n");
  }
  else {
#endif

  // R" = R + (F(R)+F(R') dt/4
  DriftJump1(0.25*dt, w);

#ifdef SECURE // R' = R + F(R) dt
  DriftJump(0.5*dt, w); // linear order move
  if(CoordinatesDistanceRpRpp(w)>1e-7) { // coordinates must not be very different doing linear or quadratic move
    SaveCoordinates("coord_bad.in");
    Error("walker %i (%.15e %.15e %.15e) has inconsistent trajectory, see coord_bad.in", w, W[w].x[0], W[w].y[0], W[w].z[0]);
  }
#endif

#ifdef DMC_OVERLAP_RETHROW
  if(Overlap_check_vector(W[w].Rpp)) {
    overlaped++;
#ifdef HARD_SPHERE_KILL_OVERLAPPED
    go = OFF; W[w].status = DEAD;
#endif
    if(verbosity) Message("DMC move: particles overlap (Rpp)\n");
  }
  else {
#endif

  // Energy in the middle of drift
  ReduceVectorToTheBox(W[w].Rpp); // energy calculation might need that all particles are in the box
  CopyVectorToWalker(&Wp[w], W[w].Rpp, W[w].spin);
  Wp[w].E = WalkerEnergy0(&Wp[w]);

  if(Wp[w].E>Eck) {
    //if(verbosity) 
    Warning("DMC WARNING: large energy %" LF "\n", Wp[w].E);
    if(branchng_present) {
      overlaped++;
      W[w].status = DEAD;
#ifdef DMC_OVERLAP_RETHROW
      go = OFF;
#endif
    }
  }
  else if(Wp[w].E<-Eck) {
    //if(verbosity) 
    Warning("DMC WARNING: negative energy %" LF "\n", Wp[w].E);
    if(branchng_present) {
      overlaped++;
      W[w].status = DEAD;
#ifdef DMC_OVERLAP_RETHROW
      go = OFF;
#endif
    }
  }
  else {
    W[w].Eold = W[w].E;
    W[w].E = Wp[w].E;
    W[w].EFF = Wp[w].EFF;
    W[w].Epot = Wp[w].Epot;
    W[w].Ekin = Wp[w].Ekin;
#ifdef INTERACTION_WITH_DAMPING
    W[w].Edamping = Wp[w].Edamping;
#endif

    // R"' = R + F(R") dt
    DriftJump2(dt, w);

#ifdef DMC_OVERLAP_RETHROW
  if(!Overlap_check_vector(W[w].Rppp)) {
#endif
    ReduceVectorToTheBox(W[w].Rppp);
    CopyVector(W[w].R_dmc_end_of_step, W[w].Rppp);

#ifdef MEASURE_CORRELATION_FUNCTIONS_AT_THE_END_OF_DRIFT
    // standard algorithm copies position at the end of the step, so that corr. functions are measured at the end of the step
    CopyVectorToWalker(&W[w], W[w].Rppp);
#endif
#ifdef MEASURE_CORRELATION_FUNCTIONS_IN_THE_MIDDLE_OF_DRIFT
    //ReduceVectorToTheBox(W[w].Rpp); // done earlier in the calculation of the energy
    CopyVectorToWalker(&W[w], W[w].Rpp, W[w].spin); // copy the position in the middle of the step, so that corr. functions will be measured in the middle of the step
#endif

#ifdef SPINFULL_TUNNELING
  U_psiT_sigma_inv(W[w]); // calculate psiT_sigma_inv

  for(i_spin=0; i_spin<N; i_spin++) {
    x = exp(W[w].psiT_sigma_inv[i_spin]); // psi_T(R, -sigma_i) / psi_T(R, sigma_i)
    flip_probability = x*x/(1.+x*x)*(1. - exp(-t_tunneling*(x+1./x)*dt));
    if(Random()<flip_probability) { // succesfull flip, change spin
      W[w].spin[i_spin] = !W[w].spin[i_spin];
    } 
  }
#endif

#ifdef HARD_SPHERE_HYPERRADIUS
    if(CheckWalkerHyperradiusOverlapping(W[w])) {
      overlaped++;
      W[w].status = DEAD;
    }
#endif

#ifdef DMC_OVERLAP_RETHROW
      go = OFF;
    }
    else {
      overlaped++;
#ifdef HARD_SPHERE_KILL_OVERLAPPED
      go = OFF; W[w].status = DEAD;
#endif
      if(verbosity) Message("DMC move: particles overlap (Rppp)\n");
    }
    }}}
#endif
  }
}

/**************************** DMC Simple Move *********************************/
void DMCLinearMoveDrift(int w) {
  int go = ON;
  int i,j;

  while(go == ON) {
    tried++;
    // R = R(i-1) + chi, f(chi) = Exp(-R^2/2dt)
    // Exp(-(x-mu)^2/(2 sigma^2))
    GaussianJump(w, Sqrt(dt));

    // R' = R + F(R) dt
    DriftJump(dt, w);

#ifdef HARD_SPHERE
    if(Overlapping(W[w].Rp)) {
      overlaped++;
    }
    else
#endif
    {
      CopyVectorToWalker(&W[w], W[w].Rp, W[w].spin);
      ReduceWalkerToTheBox(&W[w]);
#ifdef CHAINS_PINNED_TO_CENTER_OF_MASS
      CalculateAllChainsCenterOfMassWalker(&W[w]);
#endif
      W[w].Eold = W[w].E;
      W[w].E = WalkerEnergy0(&W[w]);
      go = OFF;
    }
  }
}

/**************************** DMC Simple Move *********************************/
void DMCLinearMoveGaussian(int w) {
  int go = ON;
  int i,j;

  while(go == ON) {
    tried++;
    // R = R(i-1) + chi, f(chi) = Exp(-R^2/2dt)
    // Exp(-(x-mu)^2/(2 sigma^2))
    GaussianJump(w, Sqrt(dt));

    CopyVectorToWalker(&W[w], W[w].R, W[w].spin);
    ReduceWalkerToTheBox(&W[w]);
    W[w].Eold = W[w].E;
    W[w].E = WalkerEnergy0(&W[w]);

    // R' = R + F(R) dt
    DriftJump(dt, w);

#ifdef HARD_SPHERE
    if(Overlapping(W[w].Rp)) {
      overlaped++;
    }
    else
#endif
    {
      CopyVectorToWalker(&W[w], W[w].Rp, W[w].spin);
      ReduceWalkerToTheBox(&W[w]);
      go = OFF;
    }
  }
}

/**************************** DMC Linear Metropolis Move ******************/
// Rp = R + dR
// dR = F(R)dt + chi (used also in SD calculation)
void DMCLinearMetropolisMove(int w) {
  DOUBLE dexp = 0.;
  DOUBLE xi,dx,dy,dz;
  int move_accepted;
#ifndef MEMORY_CONTIGUOUS
  int i,j;
#endif
  //DOUBLE WE, WEFF, WEpot, WEkin; // store energy calculated in the middle of drift

#ifdef MEASURE_CORRELATION_FUNCTIONS_IN_THE_MIDDLE_OF_DRIFT
  if(!W[w].status_end_of_step_initialized) // first iteration needs to be initialized
#endif
  {
    CopyWalkerToVector(W[w].R, W[w]);
    DriftForce(W[w].F, W[w].R, w);
    for(i=0; i<N; i++) { // one step backwards
      W[w].dR_drift_old[i][0] = dt*W[w].F[i][0]; // initialize with an approximate value
      W[w].dR_drift_old[i][1] = dt*W[w].F[i][1];
      W[w].dR_drift_old[i][2] = dt*W[w].F[i][2];
    }
#ifdef MEASURE_CORRELATION_FUNCTIONS_IN_THE_MIDDLE_OF_DRIFT
    W[w].status_end_of_step_initialized = ON;
#endif
  }

  // R = R(i-1) + chi, f(chi) = const Exp(-R^2/2dt)
  for(i=0; i<N; i++) {
    RandomNormal3(&dx, &dy, &dz, 0., sqrt(dt)); //(x,y,z, mu, sigma)
    CaseX(W[w].R[i][0] = W[w].x[i] + W[w].dR_drift_old[i][0] + dx); // R = R(i-1) + F(i-1)dt + chi, step starts after Gaussian
    CaseY(W[w].R[i][1] = W[w].y[i] + W[w].dR_drift_old[i][1] + dy);
    CaseZ(W[w].R[i][2] = W[w].z[i] + W[w].dR_drift_old[i][2] + dz);

    ReduceToTheBox(W[w].R[i]);

    CaseX(W[w].dR[i][0] = dx); // (1/2): dR' = chi
    CaseY(W[w].dR[i][1] = dy);
    CaseZ(W[w].dR[i][2] = dz);

    CaseX(dexp += dx*dx); // dR: direct transition probability, P(R->R') = exp(-chi^2/2dt)
    CaseY(dexp += dy*dy);
    CaseZ(dexp += dz*dz);
  }

  // R' = R + F(R) dt
  DriftForce(W[w].F, W[w].R, w);
  for(i=0; i<N; i++) {
    CaseX(dx = dt*W[w].F[i][0]);
    CaseY(dy = dt*W[w].F[i][1]);
    CaseZ(dz = dt*W[w].F[i][2]);

    CaseX(W[w].dR[i][0] += dx); // (2/2): dR' = chi + F(R')dt
    CaseY(W[w].dR[i][1] += dy);
    CaseZ(W[w].dR[i][2] += dz);

    CaseX(W[w].dR_drift[i][0] = dx); // needed for the next step
    CaseY(W[w].dR_drift[i][1] = dy);
    CaseZ(W[w].dR_drift[i][2] = dz);

    CaseX(W[w].Rp[i][0] = W[w].R[i][0] + dx); // R' = R + F(R)dt
    CaseY(W[w].Rp[i][1] = W[w].R[i][1] + dy);
    CaseZ(W[w].Rp[i][2] = W[w].R[i][2] + dz);
#ifdef PARTICLES_NEVER_LEAVE_THE_BOX
    ReduceToTheBox(W[w].Rp[i]);
#endif
  }

  // Metropolis move : w = 2[U(R')-U(R)] + [dR^2-dR'^2]/(2dt)
  CopyVectorToWalker(&Wp[w], W[w].R, W[w].spin); // calculate new weight in R' (after Gaussian)
  ReduceWalkerToTheBox(&Wp[w]);
  Wp[w].U = U(Wp[w]);

  for(i=0; i<N; i++) {
    // with saved previous displacement:
    CaseX(dx = W[w].dR[i][0] + W[w].dR_drift_old[i][0]); // (3/3): dR' = chi + F(R')dt + F(R)dt
    CaseY(dy = W[w].dR[i][1] + W[w].dR_drift_old[i][1]);
    CaseZ(dz = W[w].dR[i][2] + W[w].dR_drift_old[i][2]);

    CaseX(dexp -= dx*dx); // dR': inverse transition probability P(R'->R)
    CaseY(dexp -= dy*dy);
    CaseZ(dexp -= dz*dz);
  }
  dexp = 2.*(Wp[w].U - W[w].U) + dexp/(2.*dt); // detailed balance condition

  // The Metropolis code
  if(dexp > 0.) {
    move_accepted = ON;
  }
  else {
    xi = Random();
    if(xi<Exp(dexp)) {
      move_accepted = ON;
    }
    else {
      rejected++;
      return;
    }
  }

  if(move_accepted) {
    accepted++;
    CopyVectorToWalker(&W[w], W[w].R, W[w].spin); // move is accepted, update walker coordinates (after Gaussian)
    ReduceWalkerToTheBox(&W[w]);
    CopyVector(W[w].dR_drift_old, W[w].dR_drift); // store the deterministic shift due to the force
    W[w].U = Wp[w].U;
    W[w].Eold = W[w].E;
    W[w].E = WalkerEnergy(&W[w], &W[w].Epot, &W[w].Ekin, &W[w].EFF, &W[w].Edamping, &W[w].Eint, &W[w].Eext);
    if(measure_SD) { // update SD with dR = chi + F(R)dt
      for(i=0; i<N; i++) {
        CaseX(W[w].rreal[i][0] += W[w].dR[i][0]);
        CaseY(W[w].rreal[i][1] += W[w].dR[i][1]);
        CaseZ(W[w].rreal[i][2] += W[w].dR[i][2]);
      }
    }
  }
}

/**************************** DMC Quadratic Move *****************************/
// provides quadratic proparation in time
// but not symmetric (linear) in the point where measurements are done
void DMCQuadraticMetropolisMove(int w) {
  DOUBLE dexp=0., xi, Up;
  DOUBLE dx,dy,dz;
  int move_accepted;
#ifndef MEMORY_CONTIGUOUS
  int i,j;
#endif
  //DOUBLE WE, WEFF, WEpot, WEkin; // store energy calculated in the middle of drift

#ifdef MEASURE_CORRELATION_FUNCTIONS_IN_THE_MIDDLE_OF_DRIFT
  if(!W[w].status_end_of_step_initialized) // first iteration needs to be initialized
#endif
  {
    CopyWalkerToVector(W[w].R, W[w]);
    DriftForce(W[w].F, W[w].R, w);
    for(i=0; i<N; i++) { // one step backwards
      W[w].dR_drift_old[i][0] = dt*W[w].F[i][0]; // initialize with an approximate value
      W[w].dR_drift_old[i][1] = dt*W[w].F[i][1];
      W[w].dR_drift_old[i][2] = dt*W[w].F[i][2];
    }
#ifdef MEASURE_CORRELATION_FUNCTIONS_IN_THE_MIDDLE_OF_DRIFT
    W[w].status_end_of_step_initialized = ON;
#endif
  }

  // R = R(i-1) + chi, f(chi) = const Exp(-R^2/2dt)
  for(i=0; i<N; i++) {
    RandomNormal3(&dx, &dy, &dz, 0., sqrt(dt)); //(x,y,z, mu, sigma)
    CaseX(W[w].R[i][0] = W[w].x[i] + W[w].dR_drift_old[i][0] + dx); // R = R(i-1) + F(i-1)dt + chi, step starts after Gaussian
    CaseY(W[w].R[i][1] = W[w].y[i] + W[w].dR_drift_old[i][1] + dy);
    CaseZ(W[w].R[i][2] = W[w].z[i] + W[w].dR_drift_old[i][2] + dz);

    ReduceToTheBox(W[w].R[i]); // U(R) is measured here, so put it to the box

    CaseX(W[w].dR[i][0] = dx); // (1/2): dR' = chi
    CaseY(W[w].dR[i][1] = dy);
    CaseZ(W[w].dR[i][2] = dz);

    CaseX(dexp += dx*dx); // dR: direct transition probability, P(R->R') = exp(-chi^2/2dt)
    CaseY(dexp += dy*dy);
    CaseZ(dexp += dz*dz);
  }
  // end of Gaussian jump

  // R' = R + F(R) dt/2
  DriftForce(W[w].F, W[w].R, w);
  for(i=0; i<N; i++) {
    CaseX(W[w].Rp[i][0] = W[w].R[i][0] + 0.5*dt*W[w].F[i][0]); // 0.5 dt
    CaseY(W[w].Rp[i][1] = W[w].R[i][1] + 0.5*dt*W[w].F[i][1]);
    CaseZ(W[w].Rp[i][2] = W[w].R[i][2] + 0.5*dt*W[w].F[i][2]);
#ifdef PARTICLES_NEVER_LEAVE_THE_BOX
    ReduceToTheBox(W[w].Rp[i]);
#endif
  }

  // R'' = R + (F(R)+F(R') dt/4
  DriftForce(W[w].Fp, W[w].Rp, w);
  for(i=0; i<N; i++) {
    CaseX(W[w].Rpp[i][0] = W[w].R[i][0] + 0.25*dt*(W[w].F[i][0]+W[w].Fp[i][0]));
    CaseY(W[w].Rpp[i][1] = W[w].R[i][1] + 0.25*dt*(W[w].F[i][1]+W[w].Fp[i][1]));
    CaseZ(W[w].Rpp[i][2] = W[w].R[i][2] + 0.25*dt*(W[w].F[i][2]+W[w].Fp[i][2]));
#ifdef PARTICLES_NEVER_LEAVE_THE_BOX
    ReduceToTheBox(W[w].Rpp[i]);
#endif
  }

  // R''' = R + F(R'') dt
  DriftForce(W[w].Fpp, W[w].Rpp, w);
  for(i=0; i<N; i++) {
    dx = dt*W[w].Fpp[i][0]; // dt
    dy = dt*W[w].Fpp[i][1];
    dz = dt*W[w].Fpp[i][2];

    CaseX(W[w].dR_drift[i][0] = dx); // needed for the next step
    CaseY(W[w].dR_drift[i][1] = dy);
    CaseZ(W[w].dR_drift[i][2] = dz);

    CaseX(W[w].dR[i][0] += dx); // (2/2): dR' = chi + F(R'')dt
    CaseY(W[w].dR[i][1] += dy);
    CaseZ(W[w].dR[i][2] += dz);

    CaseX(W[w].Rppp[i][0] = W[w].R[i][0] + dx); // R''' = R + F(R'') dt
    CaseY(W[w].Rppp[i][1] = W[w].R[i][1] + dy);
    CaseZ(W[w].Rppp[i][2] = W[w].R[i][2] + dz);
#ifdef PARTICLES_NEVER_LEAVE_THE_BOX
    ReduceToTheBox(W[w].Rppp[i]);
#endif
  }

  // Metropolis move R
  // R(end) = R(i) + F(R") dt + chi
  // R' = R(end) + F(end) dt + chi'
  // w = 2[U(R')-U(R)] + 1/(2dt) [(R'-R-F(R)dt)^2-(R-R'-F(R')dt)^2]
  // w = 2[U(R')-U(R)] + [dR^2-dR'^2]/(2dt)
  // NB: can be used only if R' was not reduced to the box!
  CopyVectorToWalker(&Wp[w], W[w].R, W[w].spin); // calculate new weight, linear algorithm
  //CopyVectorToWalker(&Wp[w], W[w].Rpp); // calculate new weight in R''
  ReduceWalkerToTheBox(&Wp[w]);//not really neccessary
  Wp[w].U = U(Wp[w]);
  for(i=0; i<N; i++) {
    // with saved previous displacement:
    CaseX(dx = W[w].dR[i][0] + W[w].dR_drift_old[i][0]); // (3/3): dR' = chi + F(R')dt + F(R)dt
    CaseY(dy = W[w].dR[i][1] + W[w].dR_drift_old[i][1]);
    CaseZ(dz = W[w].dR[i][2] + W[w].dR_drift_old[i][2]);

    CaseX(dexp -= dx*dx); // dR': inverse transition probability P(R'->R)
    CaseY(dexp -= dy*dy);
    CaseZ(dexp -= dz*dz);
  }
  dexp = 2.*(Wp[w].U - W[w].U) + dexp/(2.*dt); // detailed balance condition

  // The Metropolis code
  if(dexp > 0.) {
    move_accepted = ON;
  }
  else {
    xi = Random();
    if(xi<Exp(dexp)) {
      move_accepted = ON;
    }
    else {
      rejected++;
      return;
    }
  }

  if(move_accepted) {
    accepted++;

    //CopyVectorToWalker(&W[w], W[w].Rpp); // move is accepted, update walker coordinates (measurements are done in Rpp)
    CopyVectorToWalker(&W[w], W[w].R, W[w].spin); // linear algorithm
    //ReduceWalkerToTheBox(&W[w]); // not necessary
    CopyVector(W[w].dR_drift_old, W[w].dR_drift); // store the deterministic shift due to the force
    //CopyVector(W[w].R_dmc_end_of_step, W[w].Rppp); // store the end of step
    W[w].U = Wp[w].U;
    W[w].Eold = W[w].E;
    W[w].E = WalkerEnergy(&W[w], &W[w].Epot, &W[w].Ekin, &W[w].EFF, &W[w].Edamping, &W[w].Eint, &W[w].Eext);
    if(measure_SD) { // update SD with dR = chi + F(R)dt
      for(i=0; i<N; i++) {
        CaseX(W[w].rreal[i][0] += W[w].dR[i][0]);
        CaseY(W[w].rreal[i][1] += W[w].dR[i][1]);
        CaseZ(W[w].rreal[i][2] += W[w].dR[i][2]);
      }
    }
  }
#ifdef MEASURE_CORRELATION_FUNCTIONS_IN_THE_MIDDLE_OF_DRIFT
  W[w].status_end_of_step_initialized = ON;
#endif
}

/**************************** DMC Quartic Move *********************************/
// note that in quadratic algorithm the time step is [2 dt]
// see Eqs. (33-36) in Siu A. Chin and C. R. Chen Gradient symplectic algorithms for solving the Schrodinger equation with time-dependent potentials J. Chem. Phys. 117, 1409 (2002), 
//
void DMCQuarticMove(int w) {
#ifndef MEMORY_CONTIGUOUS
  int i,j;
#endif
  DOUBLE t1, t2, t3, v0, v1, v2, v3, u0; // time steps
  DOUBLE log_weight = 0.; // logarithm of the branching weight

  // primitive, only two steps:
  //GaussianJumpWalker(w, Sqrt(dt));
  //log_weight += dt*WalkerPotentialEnergy(&W[w]);

  // t1 is a free parameter :
  // algorithm 4A:  t1 = 1/2
  // algorithm 4D:  t1 = 1/3 
  t1 = 0.01; // t1 = 0.28 - 6th order solution for harmonic oscillator

  t3 = t1;
  t2 = 1. - 2.*t1;
  v3 = v0 = (6.*t1*(t1-1.)+1.)/(12.*(t1-1.)*t1);
  v2 = v1 = 0.5 - v0;
  u0 = (1./(6.*t1*(1.-t1)*(1.-t1))-1.)/48.;

  // step 1
  log_weight += v0*dt*WalkerPotentialEnergy(&W[w]);

  // step 2
  GaussianJumpWalker(w, Sqrt(dt*t1));

  // step 3
  //log_weight += v1*dt*(WalkerPotentialEnergy(&W[w]) + u0/v2*dt*dt*WalkerPotentialEnergyTilde(&W[w]));
  log_weight += v1*dt*WalkerPotentialEnergy(&W[w]);
  log_weight += v1*dt*u0/v2*dt*dt*WalkerPotentialEnergyTilde(&W[w]);

  // step 4
  GaussianJumpWalker(w, Sqrt(dt*t2));

  // step 5
  log_weight += v2*dt*(WalkerPotentialEnergy(&W[w]) + u0/v2*dt*dt*WalkerPotentialEnergyTilde(&W[w]));

  // step 6
  GaussianJumpWalker(w, Sqrt(dt*t3));

  // step 7
  log_weight += v3*dt*WalkerPotentialEnergy(&W[w]);

  // invert the sign and subtract the energy shift
  log_weight = dt*Eo-log_weight;
  W[w].weight = Exp(log_weight);

  W[w].E = WalkerEnergy0(&W[w]);
}

/**************************** DMC Simple Move *********************************/
void DMCLinearMoveOne(int w) {
  int go = ON;
  int i,j;

#define Ns 10 // number of single-particle moves
  DOUBLE b[Ns]; // weight of a walker with one particle dispaced
  DOUBLE Ej[Ns]; // energy of a walker with one particle dispaced
  DOUBLE ei; // total weight of particle i
  DOUBLE walker_weight = 1.; // the new weight of the walker, adjusted in the branching
  DOUBLE xi; // random variable
  DOUBLE bsum;
  DOUBLE dx,dy,dz; // Gaussian displacement
  DOUBLE sigma; // Gaussian width
  DOUBLE x[Ns],y[Ns],z[Ns]; // store position of the moved particle

  DOUBLE Eold, Enew;
  DOUBLE dt_local, dt_local2;

#ifdef DMC_LINEAR_ONE_SYMMETRIC_MOVE
  dt_local = 0.5*dt;
  //dt_local = dt;
  dt_local2 = dt - dt_local;
#else
  dt_local = dt;
  //dt_local = dt / (DOUBLE) N;
#endif

  sigma = Sqrt(dt_local); // Gaussian width

  CopyWalker(&Wp[w], &W[w]); // primed walker with have a displaced particle
  W[w].E = WalkerEnergy0(&W[w]); // recalculate energy
  for(i=0; i<N; i++) { // loop over all particles
    Eold = Vext(W[w].x[i], W[w].y[i], W[w].z[i]);
    // W[w].E = WalkerEnergy0(&W[w]); // recalculate energy
    ei = 0.;
    for(j=0; j<Ns; j++) { // loop over single-particle moves
      RandomNormal3(&dx, &dy, &dz, 0., sigma); //(x,y,z, mu, sigma)

      CaseX(x[j] = W[w].x[i] + dx); 
      CaseY(y[j] = W[w].y[i] + dy); 
      CaseZ(z[j] = W[w].z[i] + dz); 

      //ReduceToTheBox(&x[j]);
      // R' = R + F(R) dt
      //DriftJump(dt, w);

      CaseX(Wp[w].x[i] = x[j]);
      CaseY(Wp[w].y[i] = y[j]);
      CaseZ(Wp[w].z[i] = z[j]);

      //Wp[w].E = WalkerEnergy0(&Wp[w]); // WalkerEnergy with the displaced particle
      //b[j] = exp(-dt*(0.5*(Wp[w].E + W[w].E) - Eo)); // here Eo is the energy per particle

      //Enew = W[w].E - Eold + Vext(x[j], y[j], z[j]);
      //b[j] = exp(-dt_local*(0.5*(Enew + W[w].E) - Eo));
      Enew = Vext(x[j], y[j], z[j]);
      b[j] = exp(-dt_local*(0.5*(Enew + Eold) - Eo/(DOUBLE) N));
      //Ej[j] = Enew;
      Ej[j] = W[w].E - Eold + Enew;

      ei += b[j];
    }
    ei /= (DOUBLE) Ns;

    // decide which move is accepted
    xi = ei*Random(); // random variable from [0; ei];
    bsum = b[0];
    j=0;
    while(xi>bsum) {
      j++;
      bsum += b[j];
    }
    // now j is the index of the accepted move;
    CaseX(W[w].x[i] = x[j]); 
    CaseY(W[w].y[i] = y[j]); 
    CaseZ(W[w].z[i] = z[j]); 
    CaseX(Wp[w].x[i] = x[j]); // update the primed walker
    CaseY(Wp[w].y[i] = y[j]);
    CaseZ(Wp[w].z[i] = z[j]);
    W[w].E = Ej[j];

    walker_weight *= ei; // store the total weight
  }

#ifdef DMC_LINEAR_ONE_SYMMETRIC_MOVE
  //// N -> 1
  //CopyWalker(&Wp[w], &W[w]); // primed walker with have a displaced particle
  sigma = Sqrt(dt_local2); // Gaussian width
  //for(i=0; i<N; i++) { // loop over all particles
  for(i=N-1; i>=0; i--) { // loop over all particles
    W[w].E = WalkerEnergy0(&W[w]); // recalculate energy
    Eold = Vext(W[w].x[i], W[w].y[i], W[w].z[i]);

    ei = 0.;
    for(j=0; j<Ns; j++) { // loop over single-particle moves
      RandomNormal3(&dx, &dy, &dz, 0., sigma); //(x,y,z, mu, sigma)

      CaseX(x[j] = W[w].x[i] + dx); 
      CaseY(y[j] = W[w].y[i] + dy); 
      CaseZ(z[j] = W[w].z[i] + dz); 

      //ReduceToTheBox(&x[j]);
      // R' = R + F(R) dt
      //DriftJump(dt, w);

      CaseX(Wp[w].x[i] = x[j]);
      CaseY(Wp[w].y[i] = y[j]);
      CaseZ(Wp[w].z[i] = z[j]);

      //Wp[w].E = WalkerEnergy0(&Wp[w]); // WalkerEnergy with the displaced particle
      //b[j] = exp(-dt*(0.5*(Wp[w].E + W[w].E) - Eo)); // here Eo is the energy per particle

      //Enew = W[w].E - Eold + Vext(x[j], y[j], z[j]);
      //b[j] = exp(-dt_local*(0.5*(Enew + W[w].E) - Eo));
      Enew = Vext(x[j], y[j], z[j]);
      b[j] = exp(-dt_local*(0.5*(Enew + Eold) - Eo/(DOUBLE) N));
      //Ej[j] = Enew;
      Ej[j] = W[w].E - Eold + Enew;

      ei += b[j];
    }
    ei /= (DOUBLE) Ns;

    // decide which move is accepted
    xi = ei*Random(); // random variable from [0; ei];
    bsum = b[0];
    j=0;
    while(xi>bsum) {
      j++;
      bsum += b[j];
    }
    // now j is the index of the accepted move;
    CaseX(W[w].x[i] = x[j]); 
    CaseY(W[w].y[i] = y[j]); 
    CaseZ(W[w].z[i] = z[j]); 
    CaseX(Wp[w].x[i] = x[j]); // update the primed walker
    CaseY(Wp[w].y[i] = y[j]);
    CaseZ(Wp[w].z[i] = z[j]);
    W[w].E = Ej[j];

    walker_weight *= ei; // store the total weight
    //// N -> 1 ends
  }
#endif

#ifdef DMC_LINEAR_ONE_SYMMETRIC_MOVE
  //walker_weight /= (DOUBLE) (Ns*2*N);
#else
  //walker_weight /= (DOUBLE) (Ns*N);
#endif

  W[w].weight = walker_weight;
  //W[w].weight = pow(walker_weight, 1./(DOUBLE) N);

  //W[w].weight = 1.;
}

/**************************** DMC Pseudoptential Move *********************************/
// Normalized Green's function, Eq. (9)
// arguments are arrays Nx3
double GreensFunctionNormalized(DOUBLE **Ro, DOUBLE **R, double tau) {
  int i,j;
  DOUBLE dR2 = 0.;
  DOUBLE ch_tau, th_tau;

  ch_tau = cosh(tau);
  th_tau = tanh(tau);

  for(i=0; i<N; i++) for(j=0; j<3; j++) dR2 += (R[i][j] - Ro[i][j]/ch_tau) * (R[i][j] - Ro[i][j]/ch_tau);
  return pow(2.*PI*th_tau, -3./2.*(DOUBLE)N)*exp(-0.5*dR2/th_tau);
}

// Eq. (7)
double GreensFunction(DOUBLE **Ro, DOUBLE **R, double tau) {
  int i,j;
  DOUBLE dR2 = 0., Ro2 = 0.;
  DOUBLE sh_tau, ch_tau, th_tau;

  sh_tau = sinh(tau);
  ch_tau = cosh(tau);
  th_tau = tanh(tau);

  for(i=0; i<N; i++) for(j=0; j<3; j++) dR2 += (R[i][j] - Ro[i][j]/ch_tau) * (R[i][j] - Ro[i][j]/ch_tau);
  for(i=0; i<N; i++) for(j=0; j<3; j++) Ro2 += Ro[i][j] * Ro[i][j];

  return pow(2.*PI*sh_tau, -3./2.*(DOUBLE)N)*exp(-0.5*Ro2*th_tau - 0.5*dR2/th_tau);
}

/**************************** DMC Pseudooptential Move *********************************/
#define CENTER_OF_MASSES_DOES_NOT_MOVE
void DMCPseudopotentialMove(int dummy) {
// 3 fermions in PBC
  int i,j,w,wp;
  DOUBLE dr[3],r2,rij_abs, rij0_abs;
  DOUBLE weight_positive, weight_negative; // correction to the branching weight
  DOUBLE dx,dy,dz;
  DOUBLE th_tau, ch_tau;
  DOUBLE weight_sum, weight_average;
  DOUBLE xi;
  DOUBLE R_CM[3];
  int lazy_particle; // particle which does not make Gaussian move

  if(N!=2 && N!=3) Error("DMCPseudopotentialMoveonly for 2 or 3 particles\n");

  // generate Gaussian distribution shift
  // R = R(i-1) + chi, f(chi) = const Exp(-R^2/2dt)
  th_tau = tanh(dt);
  ch_tau = cosh(dt);

  for(w=0; w<Nwalkers; w++) { // move all walkers
    // R is the original (Ro) vector
    W[w].R[0][0] = W[w].x[0];
    W[w].R[0][1] = W[w].y[0];
    W[w].R[0][2] = W[w].z[0];
    W[w].R[1][0] = W[w].x[1];
    W[w].R[1][1] = W[w].y[1];
    W[w].R[1][2] = W[w].z[1];
    if(N == 3) {
      W[w].R[2][0] = W[w].x[2];
      W[w].R[2][1] = W[w].y[2];
      W[w].R[2][2] = W[w].z[2];
    }

    // calculate initial weight, Eq. (10)
    W[w].r2 = 0.;
    for(i=0; i<N; i++) for(j=0; j<3; j++) W[w].r2 += W[w].R[i][j] * W[w].R[i][j];
    W[w].weight = pow(ch_tau, -3./2.*(DOUBLE)N)*exp(-0.5*W[w].r2*th_tau); // weight according to Eq. (10)

#ifdef CENTER_OF_MASSES_DOES_NOT_MOVE
    // Rp is Ro plus Gaussian move, Eq. (9)
    lazy_particle = rand() % N; // one of N particles is not moved by a Gaussian
    for(i=0; i<N; i++) { // move all particles but the last one
    //  if(i != lazy_particle) {
        RandomNormal3(&dr[0], &dr[1], &dr[2], 0., Sqrt(th_tau)); // (x,y,z, mu, sigma)
        for(j=0; j<3; j++) W[w].Rp[i][j] = W[w].R[i][j] / ch_tau + dr[j];
    //  }
    }

    R_CM[0] = R_CM[1] = R_CM[2] = 0.;
    for(i=0; i<N; i++) for(j=0; j<3; j++) R_CM[j] += W[w].Rp[i][j];
    for(j=0; j<3; j++) R_CM[j] /= (DOUBLE) N; // now R_CM is the center of mass position
    for(i=0; i<N; i++) for(j=0; j<3; j++) W[w].Rp[i][j] -= R_CM[j]; // but CM back to (0,0,0)
#else // standard way, move all particles
    // Rp is Ro plus Gaussian move, Eq. (9)
    for(i=0; i<N; i++) {
      RandomNormal3(&dr[0], &dr[1], &dr[2], 0., Sqrt(th_tau)); // (x,y,z, mu, sigma)
      for(j=0; j<3; j++) W[w].Rp[i][j] = W[w].R[i][j] / ch_tau + dr[j];
    }
#endif
  }

  // Rpp: create pool {Rpp} in a random way
  // first calculate sum of the weights
  W[0].U = W[0].weight;
  for(w=1; w<Nwalkers; w++) W[w].U = W[w-1].U + W[w].weight;
  weight_sum = W[Nwalkers-1].U;
  weight_average = weight_sum / (DOUBLE) Nwalkers;

  // generate random point on a line [0, weight_sum]
  for(w=0; w<Nwalkers; w++) { // move all walkers
    xi = Random()*weight_sum;
    //simple     wp = w; // direct copy, no change
    for(wp=0; W[wp].U<xi; wp++);// find index of the corresponding walker
    // copy walker's coordinates
    for(i=0; i<N; i++) for(j=0; j<3; j++) W[w].Rpp[i][j] = W[wp].Rp[i][j]; // copy moved coordinates
  }

  // Calculate adjusted weight according to Eq. (11)
  for(w=0; w<Nwalkers; w++) {
    weight_positive = weight_negative = 0.;
    for(wp=0; wp<Nwalkers; wp++) {
      // direct
      weight_positive += GreensFunction(W[wp].R, W[w].Rpp, dt);

      // permutations
      if(N == 2) {
        // negative permutations
        // (2,1)
        for(j=0; j<3; j++) W[wp].Rppp[1][j] = W[wp].R[0][j];
        for(j=0; j<3; j++) W[wp].Rppp[0][j] = W[wp].R[1][j];
        weight_negative += GreensFunction(W[wp].Rppp, W[w].Rpp, dt);
      }
      else { // N = 3
        // (2,3,1)
        for(j=0; j<3; j++) W[wp].Rppp[1][j] = W[wp].R[0][j];
        for(j=0; j<3; j++) W[wp].Rppp[2][j] = W[wp].R[1][j];
        for(j=0; j<3; j++) W[wp].Rppp[0][j] = W[wp].R[2][j];
        weight_positive += GreensFunction(W[wp].Rppp, W[w].Rpp, dt);

        // (3,1,2)
        for(j=0; j<3; j++) W[wp].Rppp[2][j] = W[wp].R[0][j];
        for(j=0; j<3; j++) W[wp].Rppp[0][j] = W[wp].R[1][j];
        for(j=0; j<3; j++) W[wp].Rppp[1][j] = W[wp].R[2][j];
        weight_positive += GreensFunction(W[wp].Rppp, W[w].Rpp, dt);

        // negative permutations
        // (1,3,2)
        for(j=0; j<3; j++) W[wp].Rppp[0][j] = W[wp].R[0][j];
        for(j=0; j<3; j++) W[wp].Rppp[2][j] = W[wp].R[1][j];
        for(j=0; j<3; j++) W[wp].Rppp[1][j] = W[wp].R[2][j];
        weight_negative += GreensFunction(W[wp].Rppp, W[w].Rpp, dt);

        // (3,2,1)
        for(j=0; j<3; j++) W[wp].Rppp[2][j] = W[wp].R[0][j];
        for(j=0; j<3; j++) W[wp].Rppp[1][j] = W[wp].R[1][j];
        for(j=0; j<3; j++) W[wp].Rppp[0][j] = W[wp].R[2][j];
        weight_negative += GreensFunction(W[wp].Rppp, W[w].Rpp, dt);

        // (2,1,3)
        for(j=0; j<3; j++) W[wp].Rppp[1][j] = W[wp].R[0][j];
        for(j=0; j<3; j++) W[wp].Rppp[0][j] = W[wp].R[1][j];
        for(j=0; j<3; j++) W[wp].Rppp[2][j] = W[wp].R[2][j];
        weight_negative += GreensFunction(W[wp].Rppp, W[w].Rpp, dt);
      }
    }
    // Calculate the new weight according to Eq. (11)
    //W[w].weight = weight_average; // bosons
    W[w].weight = weight_average*(weight_positive - weight_negative) / weight_positive; // smart fermions
    //simple fermions        W[w].weight = 1./(DOUBLE) Nwalkers*(weight_positive - weight_negative) / GreensFunctionNormalized(W[w].R, W[w].Rpp, dt);
    W[w].weight *= ch_tau*sqrt(ch_tau); // 2a, get rid of the CM motion

    if(W[w].weight<0) {
      //Warning("  negative weight! Killing walker\n");
      W[w].weight = 0.;
    }
  }
}

/**************************** DMC Pseudooptential Move *********************************/
void DMCPseudopotentialMove2FermionsHOSmart(int dummy) {
// 2 fermions in PBC
  int i,j,w,wp;
  DOUBLE dr[3],r2,rij_abs, rij0_abs;
  DOUBLE weight_positive, weight_negative; // correction to the branching weight
  DOUBLE dx,dy,dz;
  DOUBLE th_tau, ch_tau;
  DOUBLE weight_sum, weight_average;
  DOUBLE xi;

  if(N!=2) Error("DMCPseudopotentialMoveonly for 2 fermions\n");

  // generate Gaussian distribution shift
  // R = R(i-1) + chi, f(chi) = const Exp(-R^2/2dt)
  th_tau = tanh(dt);
  ch_tau = cosh(dt);

  for(w=0; w<Nwalkers; w++) { // move all walkers
    // R is the original (Ro) vector
    W[w].R[0][0] = W[w].x[0];
    W[w].R[0][1] = W[w].y[0];
    W[w].R[0][2] = W[w].z[0];
    W[w].R[1][0] = W[w].x[1];
    W[w].R[1][1] = W[w].y[1];
    W[w].R[1][2] = W[w].z[1];

    // Rppp are permuted coordinates
    W[w].Rppp[0][0] = W[w].x[1];
    W[w].Rppp[0][1] = W[w].y[1];
    W[w].Rppp[0][2] = W[w].z[1];
    W[w].Rppp[1][0] = W[w].x[0];
    W[w].Rppp[1][1] = W[w].y[0];
    W[w].Rppp[1][2] = W[w].z[0];

    // calculate initial weight, Eq. (10)
    W[w].r2 = 0.;
    for(i=0; i<N; i++) for(j=0; j<3; j++) W[w].r2 += W[w].R[i][j] * W[w].R[i][j];
    W[w].weight = pow(ch_tau, -3./2.*(DOUBLE)N)*exp(-0.5*W[w].r2*th_tau); // weight according to Eq. (10)

    // Rp is Ro plus Gaussian move, Eq. (9)
    RandomNormal3(&dr[0], &dr[1], &dr[2], 0., Sqrt(th_tau)); // (x,y,z, mu, sigma)
    for(i=0; i<3; i++) W[w].Rp[0][i] = W[w].R[0][i] / ch_tau + dr[i];

    RandomNormal3(&dr[0], &dr[1], &dr[2], 0., Sqrt(th_tau)); // (x,y,z, mu, sigma)
    for(i=0; i<3; i++) W[w].Rp[1][i] = W[w].R[1][i] / ch_tau + dr[i];
  }

  // Rpp: create pool {Rpp} in a random way
  // first calculate sum of the weights
  W[0].U = W[0].weight;
  for(w=1; w<Nwalkers; w++) W[w].U = W[w-1].U + W[w].weight;
  weight_sum = W[Nwalkers-1].U;
  weight_average = weight_sum / (DOUBLE) Nwalkers;

  // generate random point on a line [0, weight_sum]
  for(w=0; w<Nwalkers; w++) { // move all walkers
    xi = Random()*weight_sum;
    //simple wp = w; // direct copy, no change
    for(wp=0; W[wp].U<xi; wp++);// find index of the corresponding walker
    // copy walker's coordinates
    for(i=0; i<N; i++) for(j=0; j<3; j++) W[w].Rpp[i][j] = W[wp].Rp[i][j]; // copy moved coordinates
  }

  // Calculate adjusted weight accroding to Eq. (11)
  for(w=0; w<Nwalkers; w++) {
    weight_positive = weight_negative = 0.;
    for(wp=0; wp<Nwalkers; wp++) {
      // direct
      weight_positive += GreensFunction(W[wp].R, W[w].Rpp, dt);
      // permuted
      weight_negative += GreensFunction(W[wp].Rppp, W[w].Rpp, dt);
    }
    // Calculate the new weight according to Eq. (11)
    W[w].weight = weight_average*(weight_positive - weight_negative) / weight_positive;
    //simple     W[w].weight = 1./(DOUBLE) Nwalkers*(weight_positive - weight_negative) / GreensFunctionNormalized(W[w].R, W[w].Rpp, dt);

    if(W[w].weight<0) {
      //Warning("  negative weight! Killing walker\n");
      W[w].weight = 0.;
    }
  }
}

/**************************** DMC Pseudooptential Move *********************************/
void DMCPseudopotentialMove2FermionsHOSimple(int w) {
// 2 fermions in PBC
  int i,j,wp;
  DOUBLE dr[3],r2,rij_abs, rij0_abs;
  DOUBLE weight = 0.; // correction to the branching weight
  DOUBLE erf_arg, mult_arg, exp_arg;
  DOUBLE a1D;
  DOUBLE dx,dy,dz;
  DOUBLE dR2 = 0., R2diag = 0.;
  DOUBLE th_tau, ch_tau;

  if(N!=2) Error("DMCPseudopotentialMoveonly for 2 fermions\n");

  // generate Gaussian distribution shift
  // R = R(i-1) + chi, f(chi) = const Exp(-R^2/2dt)
  th_tau = tanh(dt);
  ch_tau = cosh(dt);

  RandomNormal3(&dx, &dy, &dz, 0., Sqrt(th_tau)); // (x,y,z, mu, sigma)
  R2diag = FindNearestImage(&dx, &dy, &dz);
  W[w].R[0][0] = W[w].x[0] / ch_tau + dx;
  W[w].R[0][1] = W[w].y[0] / ch_tau + dy;
  W[w].R[0][2] = W[w].z[0] / ch_tau + dz;

  RandomNormal3(&dx, &dy, &dz, 0., Sqrt(th_tau)); // (x,y,z, mu, sigma)
  R2diag += FindNearestImage(&dx, &dy, &dz);
  W[w].R[1][0] = W[w].x[1] / ch_tau + dx;
  W[w].R[1][1] = W[w].y[1] / ch_tau + dy;
  W[w].R[1][2] = W[w].z[1] / ch_tau + dz;

  // Eq. (5)
  for(wp=0; wp<Nwalkers; wp++) {
    // direct
    W[wp].Rp[0][0] = W[wp].x[0];
    W[wp].Rp[0][1] = W[wp].y[0];
    W[wp].Rp[0][2] = W[wp].z[0];
    W[wp].Rp[1][0] = W[wp].x[1];
    W[wp].Rp[1][1] = W[wp].y[1];
    W[wp].Rp[1][2] = W[wp].z[1];
    weight += GreensFunction(W[wp].Rp, W[w].R, dt);
    // permuted
    W[wp].Rp[0][0] = W[wp].x[1];
    W[wp].Rp[0][1] = W[wp].y[1];
    W[wp].Rp[0][2] = W[wp].z[1];
    W[wp].Rp[1][0] = W[wp].x[0];
    W[wp].Rp[1][1] = W[wp].y[0];
    W[wp].Rp[1][2] = W[wp].z[0];
    weight -= GreensFunction(W[wp].Rp, W[w].R, dt);
  }
  weight /= (DOUBLE) Nwalkers;

  weight /= pow(2.*PI*th_tau, -3./2.*(DOUBLE)N)*exp(-0.5*R2diag/th_tau); // convert the initial Gaussian distribution to uniform

  if(weight<0) {
    //Warning("  negative weight! Killing walker\n");
    weight = 0.;
  }

  //Message("%e\n", weight);
  W[w].weight = weight;
  //W[w].E = Eo;
}


/**************************** DMC Pseudooptential Move *********************************/
void DMCPseudopotentialMove2fermionsPBC(int w) {
// 2 fermions in PBC
  int i,j,wp;
  DOUBLE dr[3],r2,rij_abs, rij0_abs;
  DOUBLE weight = 0.; // correction to the branching weight
  DOUBLE erf_arg, mult_arg, exp_arg;
  DOUBLE a1D;
  DOUBLE dx,dy,dz;
  DOUBLE R2diag = 0.;
  DOUBLE R2permuted  = 0.;
  DOUBLE R2_12, R2_21, R2_11, R2_22;
  DOUBLE x1p,x2p,y1p,y2p,z1p,z2p;

  if(N!=2) Error("DMCPseudopotentialMoveonly for 2 fermions\n");

  // generate Gaussian distribution shift
  // R = R(i-1) + chi, f(chi) = const Exp(-R^2/2dt)
  RandomNormal3(&dx, &dy, &dz, 0., Sqrt(dt)); //(x,y,z, mu, sigma)
  R2diag = FindNearestImage(&dx, &dy, &dz);
  x1p = W[w].x[0] + dx;
  y1p = W[w].y[0] + dy;
  z1p = W[w].z[0] + dz;

  RandomNormal3(&dx, &dy, &dz, 0., Sqrt(dt)); //(x,y,z, mu, sigma)
  R2diag = R2diag + FindNearestImage(&dx, &dy, &dz);
  x2p = W[w].x[1] + dx;
  y2p = W[w].y[1] + dy;
  z2p = W[w].z[1] + dz;

  /*x1p = L*Random();
  y1p = L*Random();
  z1p = L*Random();
  x2p = L*Random();
  y2p = L*Random();
  z2p = L*Random();*/

  // Eq. (5)
  for(wp=0; wp<Nwalkers; wp++) {
    // 1 -> 1', 1st particle of Walker wp is displaced -> 2'
    dx = W[wp].x[0] - x1p;
    dy = W[wp].y[0] - y1p;
    dz = W[wp].z[0] - z1p;
    R2_11 = FindNearestImage(&dx, &dy, &dz);
    // 2 -> 2'
    dx = W[wp].x[1] - x2p;
    dy = W[wp].y[1] - y2p;
    dz = W[wp].z[1] - z2p;
    R2_22 = FindNearestImage(&dx, &dy, &dz);

    // 1 -> 2', 1st particle of Walker wp is displaced -> 2'
    dx = W[wp].x[0] - x2p;
    dy = W[wp].y[0] - y2p;
    dz = W[wp].z[0] - z2p;
    R2_12 = FindNearestImage(&dx, &dy, &dz);
    // 2 -> 1'
    dx = W[wp].x[1] - x1p;
    dy = W[wp].y[1] - y1p;
    dz = W[wp].z[1] - z1p;
    R2_21 = FindNearestImage(&dx, &dy, &dz);

    weight += exp(-0.5*(R2_11+R2_22)/dt) - exp(-0.5*(R2_12+R2_21)/dt);
  }
  weight /= (DOUBLE) Nwalkers;
  //weight *= pow(2.*3.141592653589*dt,-3.);

  weight *= exp(0.5*R2diag/dt); // convert the initial Gaussian distribution to uniform

  /*W[w].x[0] = x1p;
  W[w].y[0] = y1p;
  W[w].z[0] = z1p;
  W[w].x[1] = x2p;
  W[w].y[1] = y2p;
  W[w].z[1] = z2p;
  ReduceWalkerToTheBox(&W[w]);*/

  W[w].R[0][0] = x1p;
  W[w].R[0][1] = y1p;
  W[w].R[0][2] = z1p;
  W[w].R[1][0] = x2p;
  W[w].R[1][1] = y2p;
  W[w].R[1][2] = z2p;

  if(weight<0) {
    //Warning("  negative weight! Killing walker\n");
    weight = 0.;
  }

  //Message("%e\n", weight);
  W[w].weight = weight;
  //W[w].E = Eo;
}


/**************************** DMC Pseudooptential Move *********************************/
// bosonic move
void DMCPseudopotentialMoveBosonsPBC(int w) {
  int i,j;
  DOUBLE dr[3],r2,rij_abs, rij0_abs;
  DOUBLE weight = 1.; // correction to the branching weight
  DOUBLE erf_arg, mult_arg, exp_arg;
  DOUBLE a1D;
  DOUBLE dx,dy,dz;

  a1D = a;

  // generate Gaussian distribution shift
  // R = R(i-1) + chi, f(chi) = const Exp(-R^2/2dt)
  for(i=0; i<N; i++) {
    RandomNormal3(&dx, &dy, &dz, 0., Sqrt(dt)); //(x,y,z, mu, sigma)
    CaseX(W[w].dR[i][0] = dx);
    CaseY(W[w].dR[i][1] = dy);
    CaseZ(W[w].dR[i][2] = dz);

    CaseX(W[w].R[i][0] = W[w].x[i]+dx);
    CaseY(W[w].R[i][1] = W[w].y[i]+dy);
    CaseZ(W[w].R[i][2] = W[w].z[i]+dz);
  }

  // correct for the interaction
  for(i=0; i<N; i++) { // two-body terms: particle-particle force
    for(j=i+1; j<N; j++) {
      CaseX(dr[0] = W[w].x[i] - W[w].x[j]);
      CaseY(dr[1] = W[w].y[i] - W[w].y[j]);
      CaseZ(dr[2] = W[w].z[i] - W[w].z[j]);

      rij0_abs = Sqrt(FindNearestImage(&dr[0], &dr[1], &dr[2]));

      CaseX(dr[0] = W[w].R[i][0] - W[w].R[j][0]);
      CaseY(dr[1] = W[w].R[i][1] - W[w].R[j][1]);
      CaseZ(dr[2] = W[w].R[i][2] - W[w].R[j][2]);
      rij_abs = Sqrt(FindNearestImage(&dr[0], &dr[1], &dr[2]));

      exp_arg = dt/(a1D*a1D) + 0.25*(W[w].dR[i][2]-W[w].dR[j][2])*(W[w].dR[i][2]-W[w].dR[j][2])/dt - (rij_abs+rij0_abs)/a1D;
      erf_arg = (rij_abs+rij0_abs)/Sqrt(4.*dt) - Sqrt(dt)/a1D;
      mult_arg = Sqrt(PI*dt)/a1D*Exp(exp_arg);

      weight *= 1. + mult_arg*erfc(erf_arg);

      if(weight<0) Warning("  negative weight!\n");
    }
  }

  for(i=0; i<N; i++) {
    ReduceToTheBox(W[w].R[i]);

    CaseX(W[w].x[i] = W[w].R[i][0]);
    CaseY(W[w].y[i] = W[w].R[i][1]);
    CaseZ(W[w].z[i] = W[w].R[i][2]);
  }

  W[w].weight = weight;
  //W[w].E = Eo;
}

/**********************************Gaussian Jump *****************************/
// R = R + chi, f(chi) = Exp(-R^2/2dt)
void GaussianJump(int w, DOUBLE sigma) {
  int i;
  DOUBLE dx, dy, dz;
#ifdef CENTER_OF_MASS_IS_NOT_MOVED
  DOUBLE CM[3];
#endif
#ifdef CENTER_OF_MASS_Z_IS_NOT_MOVED
  DOUBLE CM[3];
#endif
#ifdef CENTER_OF_MASS_IS_NOT_MOVED_TWO_PARTS
  DOUBLE CM[3];
#endif

#ifdef SECURE
  dx = dy = dz = 0.;
#endif

  for(i=0; i<N; i++) {
    RandomNormal3(&dx, &dy, &dz, 0., sigma); //(x,y,z, mu, sigma)

#ifdef SECURE
#ifdef TRIAL_2D
    if(fabs(dz)>1e-6) Error("Gaussian jump, non zero dz");
#endif

#ifdef TRIAL_1D  // in 1D movement in the radial direction is frozen
    if(fabs(dx)>1e-6) Error("Gaussian jump, non zero dx");
    if(fabs(dy)>1e-6) Error("Gaussian jump, non zero dy");
#endif
#endif

    CaseX(W[w].R[i][0] = W[w].x[i] + dx);
    CaseY(W[w].R[i][1] = W[w].y[i] + dy);
    CaseZ(W[w].R[i][2] = W[w].z[i] + dz);

    if(measure_SD) {
      CaseX(W[w].rreal[i][0] += dx);
      CaseY(W[w].rreal[i][1] += dy);
      CaseZ(W[w].rreal[i][2] += dz);
    }
#ifdef PARTICLES_NEVER_LEAVE_THE_BOX
    ReduceToTheBox(W[w].R[i]);
#endif
  }

#ifdef FIRST_PARTICLE_IS_NOT_MOVED // do not move the first particle
  CaseX(W[w].R[0][0] = W[w].x[0]);
  CaseY(W[w].R[0][1] = W[w].y[0]);
  CaseZ(W[w].R[0][2] = W[w].z[0]);
#endif

#ifdef CENTER_OF_MASS_IS_NOT_MOVED // put the center of mass to zero
  CM[0] = CM[1] = CM[2] = 0.;
  for(i=0; i<N; i++) {
    CaseX(CM[0] += W[w].R[i][0]);
    CaseY(CM[1] += W[w].R[i][1]);
    CaseZ(CM[2] += W[w].R[i][2]);
  }
  CM[0] /= N;
  CM[1] /= N;
  CM[2] /= N;
#ifdef BC_ABSENT
  for(i=0; i<N; i++) {
    CaseX(W[w].R[i][0] -= CM[0]);
    CaseY(W[w].R[i][1] -= CM[1]);
    CaseZ(W[w].R[i][2] -= CM[2]);
  }
#endif
#endif

#ifdef CENTER_OF_MASS_Z_IS_NOT_MOVED // put the center of mass to zero
  CM[0] = CM[1] = CM[2] = 0.;
  for(i=0; i<N; i++) {
    CaseZ(CM[2] += W[w].R[i][2]);
  }
  CM[2] /= N;
#ifdef BC_ABSENT
  for(i=0; i<N; i++) {
    CaseZ(W[w].R[i][2] -= CM[2]);
  }
#endif
#endif

#ifdef CENTER_OF_MASS_IS_NOT_MOVED_TWO_PARTS // put the center of mass to zero for the first half and the second half of the system
  CM[0] = CM[1] = CM[2] = 0.;
  for(i=0; i<N/2; i++) {
    CaseX(CM[0] += W[w].R[i][0]);
    CaseY(CM[1] += W[w].R[i][1]);
    CaseZ(CM[2] += W[w].R[i][2]);
  }
  CM[0] /= (N/2);
  CM[1] /= (N/2);
  CM[2] /= (N/2);
#ifdef BC_ABSENT
  for(i=0; i<N/2; i++) {
    CaseX(W[w].R[i][0] -= CM[0]);
    CaseY(W[w].R[i][1] -= CM[1]);
    CaseZ(W[w].R[i][2] -= CM[2]);
  }
#endif

  CM[0] = CM[1] = CM[2] = 0.;
  for(i=N/2; i<N; i++) {
    CaseX(CM[0] += W[w].R[i][0]);
    CaseY(CM[1] += W[w].R[i][1]);
    CaseZ(CM[2] += W[w].R[i][2]);
  }
  CM[0] /= (N/2);
  CM[1] /= (N/2);
  CM[2] /= (N/2);
#ifdef BC_ABSENT
  for(i=N/2; i<N; i++) {
    CaseX(W[w].R[i][0] -= (CM[0]-CMseparation));
    CaseY(W[w].R[i][1] -= (CM[1]-CMseparation));
    CaseZ(W[w].R[i][2] -= (CM[2]-CMseparation));
  }
#endif
#endif

#ifdef HARD_SPHERE
  if(Overlapping(W[w].R)) {
    tried++;
    overlaped++;
    GaussianJump(w, dt);
  }
#endif
}

/*************************** Gaussian Jump Walker *****************************/
// R = R + chi, f(chi) = Exp(-R^2/2dt)
void GaussianJumpWalker(int w, DOUBLE sigma) {
  int i;
  DOUBLE dx, dy, dz;

  for(i=0; i<N; i++) {
    RandomNormal3(&dx, &dy, &dz, 0., sigma); //(x,y,z, mu, sigma)

    CaseX(W[w].x[i] += dx);
    CaseY(W[w].y[i] += dy);
    CaseZ(W[w].z[i] += dz);

    if(measure_SD) {
      CaseX(W[w].rreal[i][0] += dx);
      CaseY(W[w].rreal[i][1] += dy);
      CaseZ(W[w].rreal[i][2] += dz);
    }
#ifdef PARTICLES_NEVER_LEAVE_THE_BOX
    ReduceWalkerToTheBox(&W[w]);
#endif
  }
}

/********************************** Drift Jump *******************************/
// R' = R + F(R) dt/2, i.e. is called with dt/2
void DriftJump(DOUBLE dt, int w) {
  int i;

#ifdef ONE_BODY_TRIAL_TERMS
#ifdef CHAINS_PINNED_TO_CENTER_OF_MASS
  CalculateAllChainsCenterOfMass(W[w].R, W[w].w);
#endif
#endif

  DriftForce(W[w].F, W[w].R, w);

  for(i=0; i<N; i++) {
    CaseX(W[w].Rp[i][0] = W[w].R[i][0] + dt*W[w].F[i][0]);
    CaseY(W[w].Rp[i][1] = W[w].R[i][1] + dt*W[w].F[i][1]);
    CaseZ(W[w].Rp[i][2] = W[w].R[i][2] + dt*W[w].F[i][2]);
#ifdef PARTICLES_NEVER_LEAVE_THE_BOX
    ReduceToTheBox(W[w].Rp[i]);
#endif

    if(SmartMC == DMC_LINEAR_DRIFT && measure_SD) {
      CaseX(W[w].rreal[i][0] += dt*W[w].F[i][0]);
      CaseY(W[w].rreal[i][1] += dt*W[w].F[i][1]);
      CaseZ(W[w].rreal[i][2] += dt*W[w].F[i][2]);
    }
  }
}

/********************************** Drift Jump 1 *****************************/
// R" = R + (F(R)+F(R') dt, dt = t/4, i.e. is called with dt/4
void DriftJump1(DOUBLE dt, int w) {
  int i;

#ifdef ONE_BODY_TRIAL_TERMS
#ifdef CHAINS_PINNED_TO_CENTER_OF_MASS
  CalculateAllChainsCenterOfMass(W[w].Rp, W[w].w);
#endif
#endif

  DriftForce(W[w].Fp, W[w].Rp, w);

  for(i=0; i<N; i++) {
    CaseX(W[w].Rpp[i][0] = W[w].R[i][0] + dt*(W[w].F[i][0]+W[w].Fp[i][0]));
    CaseY(W[w].Rpp[i][1] = W[w].R[i][1] + dt*(W[w].F[i][1]+W[w].Fp[i][1]));
    CaseZ(W[w].Rpp[i][2] = W[w].R[i][2] + dt*(W[w].F[i][2]+W[w].Fp[i][2]));
#ifdef PARTICLES_NEVER_LEAVE_THE_BOX
    ReduceToTheBox(W[w].Rpp[i]);
#endif
  }
}

/********************************** Drift Jump 2 ****************************/
/* R(i+1) = R(i) + F(R") dt */
void DriftJump2(DOUBLE dt, int w) {
  int i;

#ifdef ONE_BODY_TRIAL_TERMS
#ifdef CHAINS_PINNED_TO_CENTER_OF_MASS
  CalculateAllChainsCenterOfMass(W[w].Rpp, W[w].w);
#endif
#endif

  DriftForce(W[w].Fpp, W[w].Rpp, w);

  for(i=0; i<N; i++) {
    CaseX(W[w].Rppp[i][0] = W[w].R[i][0] + dt*W[w].Fpp[i][0]);
    CaseY(W[w].Rppp[i][1] = W[w].R[i][1] + dt*W[w].Fpp[i][1]);
    CaseZ(W[w].Rppp[i][2] = W[w].R[i][2] + dt*W[w].Fpp[i][2]);
    if(measure_SD) {
      CaseX(W[w].rreal[i][0] += dt*W[w].Fpp[i][0]);
      CaseY(W[w].rreal[i][1] += dt*W[w].Fpp[i][1]);
      CaseZ(W[w].rreal[i][2] += dt*W[w].Fpp[i][2]);
    }
#ifdef PARTICLES_NEVER_LEAVE_THE_BOX
    //ReduceToTheBox(W[w].Rppp[i]); // will be reduced in the energy measurement
#endif
  }
}

/********************************** Drift Force ******************************/
/* F = (\nabla \Psi) / \Psi */
void DriftForce(DOUBLE **F, DOUBLE **R, int w) {
  int i, j;
  DOUBLE dr[3] = {0,0,0};
  DOUBLE r, r2;
  DOUBLE Force, dF;
#ifdef CRYSTAL
#ifdef CRYSTAL_SYMMETRIC
  DOUBLE M1, e;
#endif
#endif
#ifdef SCALABLE
  int cc, ci, cj;
#endif
#ifdef CALCULATE_DERIVATIVES_NUMERICALLY
  DOUBLE epsilon; // accuracy
  DOUBLE coord_store; // store coordinates
  DOUBLE Umax, Umin; // moving by epsilon
#endif
#ifdef THREE_BODY_TERMS
  DOUBLE x2, y2, HR;
  DOUBLE CM12[3], xi[3], xj[3], xk[3], yk[3];
  int k;
#endif
#ifdef MANY_BODY_TRIAL_TRAP_CENTER_OF_MASS_NODE_Z
  DOUBLE z_CM = 0.;
#endif
#ifdef MANY_BODY_TRIAL_TRAP_CENTER_OF_MASS_HERMITE2_Z
  DOUBLE z2 = 0.;
#endif

  // set Force array elements to zero
  ArrayEmpty2D(F, i, N, j, 3, DOUBLE);

#ifdef CALCULATE_DERIVATIVES_NUMERICALLY
  epsilon = numeric_derivative_dr;

  for(i=0; i<N; i++) {
#ifdef MOVE_IN_X
    coord_store = R[i][0];
    R[i][0] += epsilon;
    Umax = Uvector(R, w);
    R[i][0] = coord_store - epsilon;
    Umin = Uvector(R, w);
    R[i][0] = coord_store;
    F[i][0] = (Umax - Umin) / (2.*epsilon);
#endif
#ifdef MOVE_IN_Y
    coord_store = R[i][1];
    R[i][1] += epsilon;
    Umax = Uvector(R, w);
    R[i][1] = coord_store - epsilon;
    Umin = Uvector(R, w);
    R[i][1] = coord_store;
    F[i][1] = (Umax - Umin) / (2.*epsilon);
#endif
#ifdef MOVE_IN_Z
    coord_store = R[i][2];
    R[i][2] += epsilon;
    Umax = Uvector(R, w);
    R[i][2] = coord_store - epsilon;
    Umin = Uvector(R, w);
    R[i][2] = coord_store;
    F[i][2] = (Umax - Umin) / (2.*epsilon);
#endif
  }
  //return;
  SaveMatrix("Fnum.dat", F, N, 3);
  ArrayEmpty2D(F, i, N, j, 3, DOUBLE);
#endif // end CALCULATE_DERIVATIVES_NUMERICALLY

#ifdef CRYSTAL // add contribution from the w.f. of crystal
#ifdef CRYSTAL_SYMMETRIC
   for(j=0; j<Crystal.size; j++) {
    W[w].Mo[j] = 0.; // fill Mo array, Mo[j] = w[j] sum_i e(i,j)
    for(i=0; i<N; i++) {
      CaseX(dr[0] = R[i][0] - Crystal.x[j]);
      CaseY(dr[1] = R[i][1] - Crystal.y[j]);
      CaseZ(dr[2] = R[i][2] - Crystal.z[j]);
      r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
      //if(CheckInteractionConditionWF(dr[0], dr[1], dr[2], r2)) {
        W[w].Mo[j] += Exp(-Crystal_dot_Rx_j*dr[0]*dr[0]-Crystal_dot_Ry_j*dr[1]*dr[1]-Crystal_dot_Rz_j*dr[2]*dr[2]);
      //}
    }
    W[w].Mo[j] *= Crystal_dot_weight_j;
  }

  for(i=0; i<N; i++) {
    for(j=0; j<Crystal.size; j++) {
      CaseX(dr[0] = R[i][0] - Crystal.x[j]);
      CaseY(dr[1] = R[i][1] - Crystal.y[j]);
      CaseZ(dr[2] = R[i][2] - Crystal.z[j]);
      r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
      //if(CheckInteractionConditionWF(dr[0], dr[1], dr[2], r2)) {
        e = Exp(-Crystal_dot_Rx_j*dr[0]*dr[0]-Crystal_dot_Ry_j*dr[1]*dr[1]-Crystal_dot_Rz_j*dr[2]*dr[2]);

        M1 = -2.*Crystal_dot_weight_j*e/W[w].Mo[j];
        CaseX(F[i][0] += M1*Crystal_dot_Rx_j*dr[0]);
        CaseY(F[i][1] += M1*Crystal_dot_Ry_j*dr[1]);
        CaseZ(F[i][2] += M1*Crystal_dot_Rz_j*dr[2]);
      //}
    }
  }
#endif
#endif

#ifdef ONE_BODY_TRIAL_TERMS
  for(i=0; i<N; i++) { // external field force
    OneBodyFp(&F[i][0], &F[i][1], &F[i][2], R[i][0], R[i][1], R[i][2], i, W[w].w);
  }
#endif

#ifdef SCALABLE
  for(cc=0; cc<Ncells; cc++) { // loop over cells
    for(ci=0; ci<W[w].c_Nlocal[cc]; ci++) { // loop over particles in the local cell
      i = W[w].c_index_local[cc][ci];
      for(cj=0; cj<W[w].c_Nall[cc]; cj++) { // loop over particles in the nearest cells avoiding DOUBLE counting
        j = W[w].c_index_all[cc][cj];
        if(j>i) {
#else
  for(i=0; i<N; i++) { // two-body terms: particle-particle force
    for(j=i+1; j<N; j++) {
#endif

      CaseX(dr[0] = R[i][0] - R[j][0]);
      CaseY(dr[1] = R[i][1] - R[j][1]);
      CaseZ(dr[2] = R[i][2] - R[j][2]);
      r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);

      if(CheckInteractionConditionWF(dr[0], dr[1], dr[2], r2)) {
        r = Sqrt(r2);
        Force = InterpolateFp(&G, r, W[w].spin[i], W[w].spin[j]) / r;
        CaseX(dF = Force * dr[0]; F[i][0] += dF; F[j][0] -= dF;)
        CaseY(dF = Force * dr[1]; F[i][1] += dF; F[j][1] -= dF;)
        CaseZ(dF = Force * dr[2]; F[i][2] += dF; F[j][2] -= dF;)
#ifdef BC_SPHERE // add effect of the curvature
        Rnow = (R[i][0]*R[j][0] + R[i][1]*R[j][1] + R[i][2]*R[j][2]) / (Rsphere*Rsphere);
        F[i][0] += Force*(-dr[0] + R[i][0]*Rnow - R[j][0]);
        F[i][1] += Force*(-dr[1] + R[i][1]*Rnow - R[j][1]);
        F[i][2] += Force*(-dr[2] + R[i][2]*Rnow - R[j][2]);
        F[j][0] += Force*( dr[0] + R[j][0]*Rnow - R[i][0]);
        F[j][1] += Force*( dr[1] + R[j][1]*Rnow - R[i][1]);
        F[j][2] += Force*( dr[2] + R[j][2]*Rnow - R[i][2]);
#endif
      }

#ifdef THREE_BODY_TERMS
      //  contains y[k] = r[i] - r[j]
      yk[0] = dr[0];
      yk[1] = dr[1];
      yk[2] = dr[2];
      y2 = r2; // y[k]

      // calculate CM12 = (r[i]+r[j])/2
      CM12[0] = R[i][0] - 0.5*dr[0];
      CM12[1] = R[i][1] - 0.5*dr[1];
      CM12[2] = R[i][2] - 0.5*dr[2];

      for(k=j+1; k<N; k++) {
        // calculate x[k] = 2./sqrt(3)*(r[k] - (r[i]+r[j])/2)
        xk[0] = R[k][0] - CM12[0];
        xk[1] = R[k][1] - CM12[1];
        xk[2] = R[k][2] - CM12[2];
        x2 = 4./3.*FindNearestImage(&xk[0], &xk[1], &xk[2]);
        xk[0] *= 2./sqrt(3.);
        xk[1] *= 2./sqrt(3.);
        xk[2] *= 2./sqrt(3.);
        HR = sqrt(x2 + y2); // hyperradius

        xi[0] = -0.5*xk[0] + 0.5*sqrt(3.)*yk[0];
        xi[1] = -0.5*xk[1] + 0.5*sqrt(3.)*yk[1];
        xi[2] = -0.5*xk[2] + 0.5*sqrt(3.)*yk[2];

        xj[0] = -0.5*xk[0] - 0.5*sqrt(3.)*yk[0];
        xj[1] = -0.5*xk[1] - 0.5*sqrt(3.)*yk[1];
        xj[2] = -0.5*xk[2] - 0.5*sqrt(3.)*yk[2];

        Force = Fp3(HR)*2./(sqrt(3.)*HR);

        F[i][0] += Force*xi[0];
        F[i][1] += Force*xi[1];
        F[i][2] += Force*xi[2];

        F[j][0] += Force*xj[0];
        F[j][1] += Force*xj[1];
        F[j][2] += Force*xj[2];

        F[k][0] += Force*xk[0];
        F[k][1] += Force*xk[1];
        F[k][2] += Force*xk[2];
      }
#endif
    }
  }
#ifdef SCALABLE
  }}
#endif

#ifdef MANY_BODY_TRIAL_TRAP_CENTER_OF_MASS_NODE_Z
  for(i=0; i<N; i++) z_CM += R[i][2];
  for(i=0; i<N; i++) F[i][2] += 1./z_CM;
#endif
#ifdef MANY_BODY_TRIAL_TRAP_CENTER_OF_MASS_HERMITE2_Z
  for(i=0; i<N; i++) z2 += (R[i][2] / (DOUBLE) N)*(R[i][2] / (DOUBLE) N);
  for(i=0; i<N; i++) F[i][2] += 2.*R[i][2] / ((DOUBLE)(N*N) * (z2 - beta));
#endif
  //SaveMatrix("Fanalyt.dat", F, N, 3);
}

/********************************** Drift Force ******************************/
/* F = (\nabla \Psi) / \Psi */
void CalculateDriftForceAndKineticEnergyNumerically(struct Walker *walker, DOUBLE *Ekin) {
  int i, j;

  DOUBLE epsilon, epsilon2; // accuracy
  DOUBLE coord_store; // store coordinates
  DOUBLE Uo, Umax, Umin; // moving by epsilon
  DOUBLE up; // U' = (dU/dx)
  DOUBLE upp; // U'' = (d^2U/dx^2)

  epsilon = numeric_derivative_dr;
  epsilon2 = epsilon*epsilon;

  Uo = U(*walker);

  for(i=0; i<N; i++) {
#ifdef MOVE_IN_X
    coord_store = walker->x[i];
    walker->x[i] = coord_store + epsilon;
    Umax = U(*walker);
    walker->x[i] = coord_store - epsilon;
    Umin = U(*walker);
    walker->x[i] = coord_store;
    up = (Umax - Umin) / (2.*epsilon);
    walker->F[i][0] = up;
    upp = (Umax -2.*Uo + Umin) / epsilon2;
    *Ekin -= upp + up*up;
#endif
#ifdef MOVE_IN_Y
    coord_store = walker->y[i];
    walker->y[i] = coord_store + epsilon;
    Umax = U(*walker);
    walker->y[i] = coord_store - epsilon;
    Umin = U(*walker);
    walker->y[i] = coord_store;
    up = (Umax - Umin) / (2.*epsilon);
    upp = (Umax -2.*Uo + Umin) / epsilon2;
    walker->F[i][1] = up;
    *Ekin -= upp + up*up;
#endif
#ifdef MOVE_IN_Z
    coord_store = walker->z[i];
    walker->z[i] = coord_store + epsilon;
    Umax = U(*walker);
    walker->z[i] = coord_store - epsilon;
    Umin = U(*walker);
    walker->z[i] = coord_store;
    up = (Umax - Umin) / (2.*epsilon);
    upp = (Umax -2.*Uo + Umin) / epsilon2;
    walker->F[i][2] = up;
    *Ekin -= upp + up*up;
#endif
  }
}

/********************************** Branching ********************************/
#ifdef BRANCHING_ALWAYS_STABLE // potentially unstable branching
void Branching(void) {
  int i, w, wnew;
  DOUBLE dE, weight, log_weight;
  int w_branching, Nbranching; // not all walkers have to be branched
  int Npop_new;
  int search;
  DOUBLE xmax, xmin, x;
  DOUBLE wmax = 0.;

  if(SmartMC == DMC_PSEUDOPOTENTIAL) {
    BranchingPseudopotential();
    return;
  }

  if(!branchng_present) {
    Eo = EnergyO();
    return;
  }

  // calculate weights
  w_branching = 0;
  for(w=0; w<NwalkersMax; w++) {
   if(W[w].status == ALIVE || W[w].status == VIRTUAL) {
      if(SmartMC == DMC_QUARTIC || SmartMC == DMC_PSEUDOPOTENTIAL || SmartMC == DMC_MOVE_ONE) {
        weight = W[w].weight;
#ifdef VIRTUAL_WALKERS
        log_weight = log(weight);
#endif
      }
      else { // DMC_QUADRATIC, DMC_LINEAR_GAUSSIAN, DMC_LINEAR_MIDDLE, DMC_LINEAR_METROPOLIS
        dE = Eo - 0.5*(W[w].E + W[w].Eold);
        log_weight = dt*dE;
        weight = Exp(log_weight);
      }

#ifdef VIRTUAL_WALKERS
      if(W[w].status == VIRTUAL) { // virtual walkers keep logarithmic weight
        branching_weight[w_branching] = log_weight;
        W[w].weight += log_weight; //accumulate weight
      } else
#endif
      {
        branching_weight[w_branching] = weight;
        W[w].weight = weight;
      }

      if(weight>wmax) wmax = weight; // find maximal weight

      w_branching++;
    }
  }
  Nbranching = w_branching;

#ifdef SECURE
  if(Nbranching == 0) Warning(" all walkers are dead");
#endif

  Npop_new = 0;
  for(w=0; w<Nbranching; w++) 
#ifdef VIRTUAL_WALKERS
    if(W[w].status == ALIVE) 
#endif
  {
    branching_xi[w] = Random(); // throw random numbers
    branching_multiplicity[w] = (int)(branching_weight[w]+branching_xi[w]); // define number of sons
    if(branching_multiplicity[w]<0) {
      Warning("Branching: huge multiplicity (larger than INT_MAX), walker %i, weight %lf, status %i\n", w, branching_weight[w], W[w].status);
      branching_multiplicity[w] = NwalkersMax;
    }
    Npop_new += branching_multiplicity[w]; // check if the new population size is within the proper limits
  }

  if(Npop_new<0) { // this can happen if the timestep is too large
    Warning("Branching: huge total multiplicity (larger than INT_MAX)\n");
    Npop_new = Npop_max+1;
  }
  // adjust the weight if necessary

  //Npop_max = Npop_min = Npop; // adjust always, results in additional bias
  if(Npop_new>Npop_max || Npop_new<Npop_min) { // adjust only outside of corridor
    search = ON;

    // find xmin and xmax so that Nmin < Npop < Nmax
    if(Npop_new>Npop_max) {
      xmax = 1.; // with no shift population is larger then Npop, can be unstable
      //xmax = (DOUBLE)Npop/wmax; // approximation of the upper bound using max weight
      xmin = 0.;
    }
    else { // Npop_new<Npop_min
      xmin = 1.; // with no shift population is smaller then Npop
      xmax = 2.;
      do {
        xmax *= 2.;
        Npop_new = 0;
        for(w=0; w<Nbranching; w++) if(W[w].status == ALIVE) Npop_new += (int)(min(xmax*branching_weight[w]+branching_xi[w],INT_MAX));
        if(Npop_new>Npop) search = OFF;
      }
      while(search);
    }

    // double check
    x = xmin;
    Npop_new = 0;
    for(w=0; w<Nbranching; w++) if(W[w].status == ALIVE) Npop_new += (int)(min(x*branching_weight[w]+branching_xi[w],INT_MAX));
    if(Npop_new>Npop) 
      Error("xmin = %lf Npop = %i, target = %i\n", x, Npop_new, Npop);
    //Message("xmin = %lf -> Npop = %i, ", x, Npop_new);
    x = xmax;
    Npop_new = 0;
    for(w=0; w<Nbranching; w++) if(W[w].status == ALIVE) Npop_new += (int)(min(x*branching_weight[w]+branching_xi[w],INT_MAX));
    if(Npop_new<Npop) Error("xmax = %lf Npop = %i, target = %i\n", x, Npop_new, Npop);
    //Message("xmax = %lf -> Npop = %i, target = %i\n", x, Npop_new, Npop);

    search = ON; // do interactions until proper weights are found
    while(search) {
      x = (xmin+xmax)/2.;
      Npop_new = 0;
      for(w=0; w<Nbranching; w++) {
        if(W[w].status == ALIVE) Npop_new += (int)(min(x*branching_weight[w]+branching_xi[w],INT_MAX));
      }

      if(Npop_new == Npop)
        search = OFF;
      else if(Npop_new>Npop)
        xmax = x;
      else
        xmin = x;
    }

    for(w=0; w<Nbranching; w++) if(W[w].status == ALIVE) branching_multiplicity[w] = (int)(x*branching_weight[w]+branching_xi[w]);
  }

  // do the branching
  w_branching = 0;
  for(w=0; w<NwalkersMax; w++) {
#ifdef VIRTUAL_WALKERS
    if(W[w].status == VIRTUAL) { // do not branch virtual (OBDM pure) walker
      //W[w].weight += branching_weight[w_branching]; // virtual walker stores logarithmic weight // logarithmic weight
      //W[w].weight *= branching_weight[w_branching]; // virtual walker stores logarithmic weight // linear weight
      if((iteration_global-1) % Nmeasure == 0) { // in case of OBDM pure measurement
        W[w].OBDMweight[W[w].OBDM_position] = W[w].weight; // store total weight
        W[w].OBDM_position++;
        if(W[w].OBDM_position>grid_pure_block) Warning("OBDM_position (%i) of Walker %i is out of range (%i), iter %i!\n", W[w].OBDM_position, w, grid_pure_block, iteration_global);
      }
    }
    else
#endif
    if(W[w].status == ALIVE) { // no virtual walkers
      if(branching_multiplicity[w_branching] == 0) {
        W[w].status = KILLED;
        W[w].weight = branching_weight[w_branching];
      }
      else {
        W[w].weight = branching_weight[w_branching]; // (DOUBLE) multiplicity;
        for(i=1; i<branching_multiplicity[w_branching]; i++) { // make replicae
          wnew = dead_walkers_storage[Ndead_walkers-1];
//#ifdef SECURE
          if(Ndead_walkers < 1) Error("  cannot branch walkers: no space (try  NwalkersMax = 3 Npop)\n");
          if(wnew < 0) Error("  negative walker index %i, Ndead_walkers = %i\n", wnew, Ndead_walkers-1);
//#endif
          CopyWalker(&W[wnew], &W[w]);
          W[wnew].status = REINCARNATION;
          W[wnew].weight = 0.;
          Ndead_walkers--;
        }
      }
    }
    w_branching++;
  }

  Eo = EnergyO();

  WalkersSort(); // put all walkers into order: ALIVE, update dead_walkers_storage
}
#elif defined(BRANCHING_FIXED_NUMBER_OF_WALKERS)
/********************************** Branching ********************************/
void Branching(void) {
  int w, i, Nreweight, go, wp, generate_new_random_number;
  int wneg = -1; // initiazized with negative value
  int wpos = -1;
#define MAXNWALKERS 10000
  int walker_sign[MAXNWALKERS];
  int index_of_killed_walker[MAXNWALKERS];
  DOUBLE walker_weight[MAXNWALKERS];
  DOUBLE Wglobal=0., Wpos=0., Wneg=0., Wtot;
  DOUBLE arrow_position;
  DOUBLE Eref;
  static int first_time = ON;

  if(first_time) {
    first_time = OFF;
    E_from_growth_estimator = (DOUBLE) Nwalkers;
  }

  Eref = (0.5-0.837308394272014);

  if(Nwalkers>MAXNWALKERS) Error("Branching: too many walkers, increase MAXNWALKERS = %i", MAXNWALKERS);

  for(w=0; w<Nwalkers; w++) { // calculate weights
    W[w].weight = Exp(dt*(Eref - 0.5*(W[w].E + W[w].Eold))); // weight
    //W[w].weight = Exp(dt*(Eref - W[w].E)); // weight
    //W[w].weight = Exp(dt*(Eo - 0.5*(W[w].E + W[w].Eold))); // weight
    //W[w].weight = Exp(dt*(Eo - W[w].E)); // weight
    Wglobal += W[w].weight; // calculate global weight
  }
  //Wglobal /= (DOUBLE) Npop;
  Wglobal /= Nwalkers;

  //E_from_growth_estimator /= Wglobal;

  E_from_growth_estimator = - log(Wglobal) / dt;

  for(w=0; w<Nwalkers; w++) { // calculated scaled weights
    W[w].weight /= Wglobal; // scaled weights
    if(W[w].weight > 1.) {
      walker_sign[w] = 1;
      Wpos += W[w].weight - 1.;
      walker_weight[w] = Wpos;
    }
    else {
      walker_sign[w] = -1;
      Wneg += 1. - W[w].weight;
      walker_weight[w] = Wneg;
    }
  }

  if(fabs(Wpos-Wneg)>1e-3*fabs(Wpos)) Warning("Branching incorrect weights: pos %f, neg %f, dif %e\n", Wpos, Wneg, Wpos-Wneg);

  // decide how many walkers should be reweighted
  Nreweight = (int) (Wpos + Random());
  accepted += Nreweight;
  rejected += Nwalkers - accepted;

  accepted / (accepted + rejected);

  Wtot = 0.5*(Wpos+Wneg);
  for(i=0; i<Nreweight; i++) { // do reweighting
    arrow_position = Wtot * Random();
    // find index of positive walker to replicate (wpos)
    go = ON;
    for(w=0; w<Nwalkers && go; w++) {
      if(walker_sign[w] == 1 && walker_weight[w] > arrow_position) {
        go = OFF;
        wpos = w;
      }
    }

    do{
      go = ON;
      arrow_position = Wtot * Random();
      generate_new_random_number = OFF;
      for(w=0; w<NwalkersMax && go; w++) {
        if(walker_sign[w] == -1 && walker_weight[w] > arrow_position) {
          go = OFF;
          for(wp=0; wp<i; wp++)
            if(w == index_of_killed_walker[wp]) {
              Warning("from %i negative walkers, walker %i (index %i) coincides walker %i and will be rethrown\n", Nreweight, i+1, w, wp+1);
              go = OFF;
              generate_new_random_number = ON;
            }
          wneg = w;
          index_of_killed_walker[i] = w;
        }
      }
    }
    while(generate_new_random_number);

    // kill walker wneg, relpicate walker wpos
    CopyWalker(&W[wneg], &W[wpos]);
    //Message("reweighting walkers: %i (%i) -> (%i)\n", Nreweight, wpos, wneg);
  }


/*  multiplicity = (int) (weight + Random());

#ifdef VIRTUAL_WALKERS
  if(W[w].status == VIRTUAL) {
    W[w].weight += log_weight; // virtual walker stores logarithmic weight
    if((iteration_global-1) % Nmeasure == 0) { // in case of OBDM pure measurement
      W[w].OBDMweight[W[w].OBDM_position] = W[w].weight; // store total weight
                                                         //if(W[w].OBDM_position) W[w].OBDMweight[W[w].OBDM_position] /= W[w].OBDMweight[W[w].OBDM_position-1]; // convert to relative weight
      W[w].OBDM_position++;
      if(W[w].OBDM_position>grid_pure_block) Warning("OBDM_position (%i) of Walker %i is out of range (%i), iter %i!\n", W[w].OBDM_position, w, grid_pure_block, iteration_global);
    }
    return 1; // do not branch virtual walker
  }
#endif
  */

  for(w=0; w<Nwalkers; w++) W[w].weight = 1.;

  Eo = EnergyO();

  //WalkersSort();
}
#else // branching not always stable
/******************************* Branching Walker ****************************/
int BranchingWalker(int w) {
  int multiplicity;
  int i, wp;
  DOUBLE dE, weight, log_weight;

#ifdef SECURE
  if(W[w].status != ALIVE && W[w].status != VIRTUAL) Error("Trying to replicate walker %i which is not alive (status %i)\n", w, W[w].status);
#endif

  if(SmartMC == DMC_MOVE_ONE || SmartMC == DMC_QUARTIC) {
    weight = W[w].weight;
  }
  else {
    dE = Eo - 0.5*(W[w].E + W[w].Eold);
    log_weight = dt*dE;
    weight = Exp(log_weight);
  }

  multiplicity = (int) (weight + Random());

#ifdef VIRTUAL_WALKERS
  if(W[w].status == VIRTUAL) {
    W[w].weight += log_weight; // virtual walker stores logarithmic weight
    if((iteration_global-1) % Nmeasure == 0) { // in case of OBDM pure measurement
      W[w].OBDMweight[W[w].OBDM_position] = W[w].weight; // store total weight
      //if(W[w].OBDM_position) W[w].OBDMweight[W[w].OBDM_position] /= W[w].OBDMweight[W[w].OBDM_position-1]; // convert to relative weight
      W[w].OBDM_position++;
      if(W[w].OBDM_position>grid_pure_block) Warning("OBDM_position (%i) of Walker %i is out of range (%i), iter %i!\n", W[w].OBDM_position, w, grid_pure_block, iteration_global);
    }
    return 1; // do not branch virtual walker
  }
#endif

  if(Nwalkers> Npop_max)
    multiplicity = (int) ((DOUBLE) multiplicity*reduce + Random());
  else if(Nwalkers < Npop_min)
    multiplicity = (int) ((DOUBLE) multiplicity*amplify + Random());

  if(multiplicity == 0) {
    W[w].status = KILLED;
    W[w].weight = weight;
  }
  else {
    W[w].weight = weight;// / (DOUBLE) multiplicity;
    for(i=1; i<multiplicity; i++) { // make replicae
      if(Ndead_walkers < 1) {
        if(multiplicity > Npop_max) {
          Warning("walker %i was killed because it had too many sons (%i)\n", w, multiplicity);
          multiplicity = 0;
          W[w].status = DEAD;
          return 0;
        }
        else
          Error("too many walkers\n Eo=%" LE "  E[%i]=%" LE "  sons %i\n", Eo, w, W[w].E, multiplicity);
      }
      wp = dead_walkers_storage[Ndead_walkers-1];
#ifdef SECURE
      if(W[wp].status == ALIVE || W[wp].status == REINCARNATION) Error("Branching, walker %i, status %i", wp, W[wp].status);
#endif
      CopyWalker(&W[wp], &W[w]);
      W[wp].status = REINCARNATION;
      W[wp].weight = 0.;
      Ndead_walkers--;
    }
  }
  return multiplicity;
}

/********************************** Branching ********************************/
void Branching(void) {
  int w;
  int Nwalkers_pure_check;

  //amplify = (DOUBLE) Npop_max / (DOUBLE) Nwalkers;
  //reduce  = (DOUBLE) Npop_min / (DOUBLE) Nwalkers;
  //Eo = EnergyO();

  // - df/dtau = (Eo - E_from_normalization)f
  //  f(0) = Nwalkers
  //  f(dt) = Npop
  E_from_growth_estimator = Eo + log((DOUBLE)Npop/(DOUBLE)Nwalkers) / dt;
  E_from_growth_estimator /= (DOUBLE) N;

  if(branchng_present)
    for(w=0; w<NwalkersMax; w++)
      if(W[w].status == ALIVE || W[w].status == VIRTUAL)
        BranchingWalker(w);

  Eo = EnergyO();

  WalkersSort();
}
#endif

/********************************* Copy To Walker ****************************/
void CopyVectorToWalker(struct Walker* W, DOUBLE** R, int *spin) {
  int i;

  for(i=0; i<N; i++) {
    CaseX(W->x[i] = R[i][0]);
    CaseY(W->y[i] = R[i][1]);
    CaseZ(W->z[i] = R[i][2]);
#ifdef SPINFULL
    W->spin[i] = spin[i];
#endif
  }
#ifdef CHAINS_PINNED_TO_CENTER_OF_MASS
  CalculateAllChainsCenterOfMassWalker(W);
#endif
  //W->weight = 1.;
}

/********************************* Copy To Vector ****************************/
void CopyWalkerToVector(DOUBLE** R, struct Walker W) {
  int i;

  for(i=0; i<N; i++) {
    CaseX(R[i][0] = W.x[i]);
    CaseY(R[i][1] = W.y[i]);
    CaseZ(R[i][2] = W.z[i]);
  }
}

/********************************* Copy Vector *******************************/
void CopyVector(DOUBLE** out, DOUBLE** in) {
  int i;

  for(i=0; i<N; i++) {
    CaseX(out[i][0] = in[i][0]);
    CaseY(out[i][1] = in[i][1]);
    CaseZ(out[i][2] = in[i][2]);
  }
}

/********************************** Energy O *********************************/
DOUBLE EnergyO(void) {
  int w;
  DOUBLE E;

  E = 0;
  EFF = 0.;
  Epot = 0.;
  Ekin = 0.;
  Eint = 0.;
  Eext = 0.;
  Nwalkersw = 0.;
#ifdef INTERACTION_WITH_DAMPING
  Edamping = 0.;
#endif
  for(w=0; w<NwalkersMax; w++) {
    if(W[w].status == ALIVE || W[w].status == KILLED || W[w].status == REINCARNATION) {
      E += W[w].weight*W[w].E;
      EFF += W[w].weight*W[w].EFF;
      Epot += W[w].weight*W[w].Epot;
      Ekin += W[w].weight*W[w].Ekin;
      Eint += W[w].weight*W[w].Eint;
      Eext += W[w].weight*W[w].Eext;

#ifdef INTERACTION_WITH_DAMPING
      Edamping += W[w].weight*W[w].Edamping;
#endif
      Nwalkersw += W[w].weight;
    }
  }
  E /= Nwalkersw;
  EFF /= Nwalkersw;
  Epot /= Nwalkersw;
  Ekin /= Nwalkersw;
  Eint /= Nwalkersw;
  Eext /= Nwalkersw;
#ifdef INTERACTION_WITH_DAMPING
  Edamping /= Nwalkersw;
#endif

  return E;
}

/********************************** Branching ********************************/
void BranchingPseudopotential(void) {
  int i, w, wnew;
  DOUBLE dE, weight;
  int w_branching, Nbranching; // not all walkers have to be branched
  int Npop_new;
  int search;
  DOUBLE xmax, xmin, x;
  DOUBLE wmax = 0.;
  DOUBLE total_weight = 0.;

  /*// Update coordinates
  for(w=0; w<NwalkersMax; w++) {
    W[w].x[0] = W[w].R[0][0];
    W[w].y[0] = W[w].R[0][1];
    W[w].z[0] = W[w].R[0][2];
    W[w].x[1] = W[w].R[1][0];
    W[w].y[1] = W[w].R[1][1];
    W[w].z[1] = W[w].R[1][2];
    ReduceWalkerToTheBox(&W[w]);
  }*/

  // Update walkers' coordinates
  for(w=0; w<Nwalkers; w++) {
    W[w].x[0] = W[w].Rpp[0][0];
    W[w].y[0] = W[w].Rpp[0][1];
    W[w].z[0] = W[w].Rpp[0][2];
    W[w].x[1] = W[w].Rpp[1][0];
    W[w].y[1] = W[w].Rpp[1][1];
    W[w].z[1] = W[w].Rpp[1][2];
    if(N == 3) {
      W[w].x[2] = W[w].Rpp[2][0];
      W[w].y[2] = W[w].Rpp[2][1];
      W[w].z[2] = W[w].Rpp[2][2];
    }
  }

  // calculate total weight
  for(w=0; w<NwalkersMax; w++) {
    if(W[w].status == ALIVE || W[w].status == VIRTUAL) {
      total_weight += W[w].weight;
    }
  }
  if(total_weight<1e-3) {
    Message("  small total weight = %e\n", total_weight);
  }

  // correct Eo in such a way that the total normalization is Npop
  // Weight' = Exp(Eo*dt) Weight = Npop
  Eo = log((DOUBLE)Npop/total_weight)/dt;
  //weight_adjust_factor = Exp(Eo*dt);
  //weight_adjust_factor = (DOUBLE)Npop/total_weight;

  // calculate weights
  w_branching = 0;
  for(w=0; w<NwalkersMax; w++) {
   if(W[w].status == ALIVE || W[w].status == VIRTUAL) {
      weight = W[w].weight*Exp(Eo*dt);

      branching_weight[w_branching] = weight; // assign weights array
      W[w].weight = weight;
      if(weight>wmax) wmax = weight; // find maximal weight
      w_branching++;
    }
  }
  Nbranching = w_branching;

  Npop_new = 0;
  for(w=0; w<Nbranching; w++) {
    branching_xi[w] = Random(); // throw random numbers
    branching_multiplicity[w] = (int)(branching_weight[w]+branching_xi[w]); // define number of sons
    if(branching_multiplicity[w]<0) {
      Warning("Branching: huge multiplicity (larger than INT_MAX)\n");
      branching_multiplicity[w] = NwalkersMax;
    }
    Npop_new += branching_multiplicity[w]; // check if the new population size is within the proper limits
  }

  if(Npop_new<0) { // this can happen if the timestep is too large
    Warning("Branching: huge total multiplicity (larger than INT_MAX)\n");
    Npop_new = Npop_max+1;
  }

  // do the branching
  w_branching = 0;
  for(w=0; w<NwalkersMax; w++) {
    if(W[w].status == ALIVE){ // not virtual walker
      if(branching_multiplicity[w_branching] == 0) {
        W[w].status = KILLED;
        W[w].weight = branching_weight[w_branching];
      }
      else {
        W[w].weight = branching_weight[w_branching]; // (DOUBLE) multiplicity;
        for(i=1; i<branching_multiplicity[w_branching]; i++) { // make replicae
          wnew = dead_walkers_storage[Ndead_walkers-1];
          CopyWalker(&W[wnew], &W[w]);
          W[wnew].status = REINCARNATION;
          W[wnew].weight = 0.;
          Ndead_walkers--;
          if(Ndead_walkers == 0) Error("Ndead_walkers = 0\n");
        }
      }
    }
    w_branching++;
  }

  WalkersSort();
}

/********************************** Branching ********************************/
void WalkersSort(void) {
  int w, Nwalkers_pure_check;

  Nwalkers = 0;
  Nwalkersw = 0;
  Nwalkers_weight_killed = 0.;
  for(w=0; w<NwalkersMax; w++) {
    if(W[w].status == REINCARNATION)
      W[w].status = ALIVE;
    else if(W[w].status == KILLED) {
      W[w].status = DEAD;
      Nwalkers_weight_killed += W[w].weight;
    }
  }

  for(w=0; w<NwalkersMax; w++) { // order all alive walkers
    if(W[w].status == ALIVE) {
      if(w != Nwalkers) { // walker should be moved to Nwalkers position
        if(W[Nwalkers].status == DEAD) { // overwrite
          CopyWalker(&W[Nwalkers], &W[w]);
          W[w].status = DEAD;
        }
        else { // exchange
          CopyWalker(&W[dead_walkers_storage[Ndead_walkers-1]], &W[Nwalkers]);
          W[dead_walkers_storage[Ndead_walkers-1]].status = W[Nwalkers].status;
          CopyWalker(&W[Nwalkers], &W[w]);
          W[Nwalkers].status = ALIVE;
          CopyWalker(&W[w], &W[dead_walkers_storage[Ndead_walkers-1]]);
          W[w].status = W[dead_walkers_storage[Ndead_walkers-1]].status;
          W[dead_walkers_storage[Ndead_walkers-1]].status = DEAD;
        }
      }
      Nwalkers++;
      Nwalkersw += W[w].weight;
    }
  }

  if(Nwalkers == 0) Error("All walkers have died");

#ifdef VIRTUAL_WALKERS
  Nwalkers_pure_check = 0;
  for(w=0; w<NwalkersMax; w++) { // order all virtual walkers
    if(W[w].status == VIRTUAL) {
      if(w != Nwalkers+Nwalkers_pure_check) {
        CopyWalker(&W[Nwalkers+Nwalkers_pure_check], &W[w]);
        W[w].status = DEAD;
        //Warning("  reordering virtual walkers: from position %i to %i\n", w, Nwalkers+Nwalkers_pure_check);
      }
      Nwalkers_pure_check++;
    }
  }
  //if(Nwalkers_pure_check != Nwalkers_pure) {
  //  Warning("  Branching: Nwalkers_pure is not self consistent (real %i, expected %i)\n", Nwalkers_pure_check, Nwalkers_pure);
  //  Nwalkers_pure = Nwalkers_pure_check;
  //}
  Nwalkers_pure = Nwalkers_pure_check;
#endif

  Ndead_walkers = 0;
  for(w=NwalkersMax-1; w>=Nwalkers+Nwalkers_pure; w--) {
    dead_walkers_storage[Ndead_walkers++] = w;
#ifdef SECURE
    if(W[w].status != DEAD) Error("Branching (2)");
#endif
  }

  /*for(w=0; w<NwalkersMax; w++) { // dump walkers
    if(W[w].status == ALIVE) Message("A");
    else if(W[w].status == DEAD) Message("D");
    else if(W[w].status == REINCARNATION) Message("R");
    else if(W[w].status == KILLED) Message("K");
    else if(W[w].status == VIRTUAL) Message("V");
  }
  Message("\n");*/
}
