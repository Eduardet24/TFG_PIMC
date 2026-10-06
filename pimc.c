/** @file pimc.c
 *  @brief Implementation of the PIMC algorithm.
*/
/*pimc.c*/
#include <stdio.h>
#include <limits.h>
#include <memory.h>
#include "main.h"
#include "randnorm.h"
#include "rw.h"
#include "compatab.h"
#include "memory.h"
#include "vmc.h"
#include "mymath.h"
#include "pimc.h"
#include "pigs.h"
#include MATHINCLUDE

/************************** PIMC Move One By One ****************************/
/** @brief Perform a single PIMC move for one bead.
 *  @param w Index of the bead to be moved.
*/
void PIMCMoveOneByOne(int w) {
  int i,j;
  int w1,w2; // index of the walker with [w-1] and [w+1]
  DOUBLE du;
  DOUBLE xi;
  DOUBLE dx,dy,dz;
  DOUBLE x,y,z;
  DOUBLE xp,yp,zp;
  DOUBLE dr[3];
  DOUBLE xm,ym,zm; // position of the middle point
  DOUBLE r2;
  DOUBLE tau;
  //DOUBLE sigma;

  w1 = w-1; // previous bead number
  w2 = w+1; // next bead number
  if(w1<0) w1 = Nwalkers-1; // cyclic condition
  if(w2==Nwalkers) w2 = 0;

  tau = 1./(T*(DOUBLE)Nwalkers); // inverse temperature on a single bead (epsilon)
  //sigma = 0.5*tau;

  // Move particle i of the bead
  for(i=0; i<N; i++) { // loop over particles
    du = 0.;

    x = W[w].x[i];
    y = W[w].y[i];
    z = W[w].z[i];
    xm = 0.5*(W[w1].x[i]+W[w2].x[i]); // middle point
    ym = 0.5*(W[w1].y[i]+W[w2].y[i]);
    zm = 0.5*(W[w1].z[i]+W[w2].z[i]);

    // Wp is the trial position
    RandomNormal3(&dx, &dy, &dz, 0., Sqrt(dt_vmc)); // generate a random number distributed according to Normal distribution

/** @brief Generate a random number distributed according to a Normal distribution.
 *  @param dx Pointer to the x-component of the random vector.
 *  @param dy Pointer to the y-component of the random vector.
 *  @param dz Pointer to the z-component of the random vector.
 *  @param mean The mean of the Normal distribution.
 *  @param sigma The standard deviation of the Normal distribution.
 */
#ifdef TRIAL_3D
    xp = x + dx;
    yp = y + dy;
    zp = z + dz;
#endif

#ifdef TRIAL_2D
    xp = x + dx;
    yp = y + dy;
    zp = 0.;
#endif

#ifdef TRIAL_1D
    xp = 0.;
    yp = 0.;
    zp = z + dz;
#endif
/** @brief Reduce the position to the simulation box.
 *  @param xp Pointer to the x-coordinate of the position.
 *  @param yp Pointer to the y-coordinate of the position.
 *  @param zp Pointer to the z-coordinate of the position.
 */
    ReduceToTheBoxXYZ(&xp, &yp, &zp); // in case of a box, put the new position to the box

    // external potential
    du -= tau*Vext(xp,yp,zp);
    du += tau*Vext(x,y,z);

    // new kinetic energy, (R' - Rm)^2
    dr[0] = xp - xm;
    dr[1] = yp - ym;
    dr[2] = zp - zm;

    r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
    //du -= r2/(2.*sigma);
    du -= r2/(2.*tau);

    // old kinetic energy
    dr[0] = x - xm;
    dr[1] = y - ym;
    dr[2] = z - zm;
    r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
    //du += r2/(2.*sigma);
    du += r2/(2.*tau);

    for(j=0; j<N; j++) {
      if(j != i) { // calculate new interaction energy
        dr[0] = xp - W[w].x[j];
        dr[1] = yp - W[w].y[j];
        dr[2] = zp - W[w].z[j];
        r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
        du -= tau*InteractionEnergy(sqrt(r2));

        // calculate old interaction energy
        dr[0] = x - W[w].x[j];
        dr[1] = y - W[w].y[j];
        dr[2] = z - W[w].z[j];
        r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
        du += tau*InteractionEnergy(sqrt(r2));

        /*// calculate new interaction energy
        dr[0] = xp - W[w1].x[j];
        dr[1] = yp - W[w1].y[j];
        dr[2] = zp - W[w1].z[j];
        r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
        du -= tau*InteractionEnergy(sqrt(r2));
        dr[0] = xp - W[w2].x[j];
        dr[1] = yp - W[w2].y[j];
        dr[2] = zp - W[w2].z[j];
        r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
        du -= tau*InteractionEnergy(sqrt(r2));

        // calculate old interaction energy
        dr[0] = x - W[w1].x[j];
        dr[1] = y - W[w1].y[j];
        dr[2] = z - W[w1].z[j];
        r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
        du += tau*InteractionEnergy(sqrt(r2));
        dr[0] = x - W[w2].x[j];
        dr[1] = y - W[w2].y[j];
        dr[2] = z - W[w2].z[j];
        r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);*/
        du += tau*InteractionEnergy(sqrt(r2));
      }
    }

    // The Metropolis code
    // f = Exp(u), p = Exp(up-u)
    if(du > 0.) {
      W[w].x[i] = xp;
      W[w].y[i] = yp;
      W[w].z[i] = zp;
      W[w].U -= du;
      accepted++;
    }
    else {
      xi = Random(); // generate a number between (0,1)
      if(xi<Exp(du)) {
        W[w].x[i] = xp;
        W[w].y[i] = yp;
        W[w].z[i] = zp;
        W[w].U -= du;
        accepted++;
      }
      else {
        rejected++;
      }
    }
  }
}


/*************************** PIGS pseudopotential *************************************/
/** @brief Perform a single PIMC move for the pseudopotential.
 *  @param w Index of the walker to be moved.
 */
void PIMCPseudopotentialMoveOneByOne(int w) {
  int i,j;
  int w1,w2; // index of the walker with [w-1] and [w+1]
  DOUBLE dS, dSkin, dSdelta, dSwf;
  DOUBLE dx,dy,dz;
  DOUBLE tau;
  int chain_edge = OFF;
  int up;
  DOUBLE x, y, z;
  DOUBLE xm, ym, zm; // position of the middle point
  //DOUBLE mA; // i.e. 2 mu
  //DOUBLE aA;
#ifdef EXTERNAL_POTENTIAL
  DOUBLE dEpot;
#endif
#ifdef SECURE
  DOUBLE dS_check;
  DOUBLE dSkin_check_old, dSpot_check_old, dSwf_check_old;
  DOUBLE dSkin_check_new, dSpot_check_new, dSwf_check_new;
#endif

  x = y = z = 0.;
  xm = ym = zm = 0.;

  //tau = dt; // imaginary time step
  tau = beta;

  w1 = w-1;
  w2 = w+1;
  if(w == 0) { // treat w=0 and w=Nwalkers-1 separately, according to the cyclic condition
    w1 = Nwalkers-1;
  }
  if(w == Nwalkers-1) { // treat w=0 and w=Nwalkers-1 separately
    w2 = 0;
  }

  // Move particle i of the bead
  for(i=0; i<N; i++) { // move i-th particle to the point rm
    dSkin = dSdelta = dSwf = 0.;
#ifdef EXTERNAL_POTENTIAL
    dEpot = 0.;
#endif

    CaseZ(Wp[w].z[i] = W[w].z[i]); // store coordinates

#ifdef SECURE
   dS_check = -PathIntegralPseudopotentialSdetailed(&dSkin_check_old, &dSpot_check_old, &dSwf_check_old);
#endif

    // Wp is the trial position
    RandomNormal3(&dx, &dy, &dz, 0., Sqrt(dt));
    W[w].z[i] += dz;
#ifdef BC_1DPBC
    ReduceToTheBoxXYZ(&W[w].x[i], &W[w].y[i], &W[w].z[i]);
#endif

  // contribution from the Laplacian, - (R'(i)-R(i+/-1))^2 / 2 + (R(i)-R(i+/-1))^2 / 2
#ifndef INFINITE_MASS_COMPONENT_UP // mA finite
    dSkin += PropagatorLaplacian(mA, tau, W[w].z[i] - W[w1].z[i]);
    dSkin += PropagatorLaplacian(mA, tau, W[w].z[i] - W[w2].z[i]);
    dSkin -= PropagatorLaplacian(mA, tau, Wp[w].z[i] - W[w1].z[i]);
    dSkin -= PropagatorLaplacian(mA, tau, Wp[w].z[i] - W[w2].z[i]);
#endif
#ifdef EXTERNAL_POTENTIAL
  dEpot += VextUp(W[w].x[i], W[w].y[i], W[w].z[i], i); // new external potential energy
  dEpot -= VextUp(Wp[w].z[i],Wp[w].y[i],Wp[w].z[i],i); // old external potential energy
#endif

    // contribution from Pseudopotential
    for(j=0; j<N; j++) if(j != i) { // interaction energy: AA
      dSdelta += PropagatorPseudopotenial1D_ij(mA, aA, tau, W[w].z[i],  W[w].z[j], W[w1].z[i], W[w1].z[j]);
      dSdelta += PropagatorPseudopotenial1D_ij(mA, aA, tau, W[w].z[i],  W[w].z[j], W[w2].z[i], W[w2].z[j]);
      dSdelta -= PropagatorPseudopotenial1D_ij(mA, aA, tau, Wp[w].z[i], W[w].z[j], W[w1].z[i], W[w1].z[j]);
      dSdelta -= PropagatorPseudopotenial1D_ij(mA, aA, tau, Wp[w].z[i], W[w].z[j], W[w2].z[i], W[w2].z[j]);
    }

#ifdef SECURE
    if(isnan(dSdelta)) Warning("NaN weight in PIGS move!\n");
#endif

    dS = dSkin + dSdelta + dSwf; // dSdelta - contribution from delta function, dSwf from the trial w.f.
#ifdef EXTERNAL_POTENTIAL
    dS -= tau*dEpot; // potential energy contribution
#endif
#ifdef SECURE
   dS_check += PathIntegralPseudopotentialSdetailed(&dSkin_check_new, &dSpot_check_new, &dSwf_check_new);
   if(fabs((dS-dS_check)/dS)>1e-5) {
     Message(" dS check (%lf %lf) mismatch dSkin %lf %lf, dSpot %lf %lf\n", dS, dS_check, dSkin, dSkin_check_new-dSkin_check_old, dSdelta, dSpot_check_new-dSpot_check_old);
     SaveCoordinates("check.dat");
   }
#endif

    if(Random() < Exp(dS)) { // The Metropolis code
      W[w].U -= dS;
      accepted_one++;
    }
    else {
      CaseX(W[w].x[i] = Wp[w].x[i]);
      CaseY(W[w].y[i] = Wp[w].y[i]);
      CaseZ(W[w].z[i] = Wp[w].z[i]);
      rejected_one++;
    }
  }
}

