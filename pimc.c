/** @file pimc.c
 *  @brief Implementation of the PIMC algorithm.
*/
/*pimc.c*/
#include <stdio.h>
#include <stdlib.h>
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
#include "utils.h"
#include MATHINCLUDE

#ifdef PIMC_HARMONIC_NORMAL_MODES
/******************* Exact Harmonic PIMC Normal-Mode Draw ******************/
void PIMCMoveHarmonicNormalModes(void) {
  int bead, mode;
  int beads = Nwalkers;
  DOUBLE tau, omega2, eigenvalue;
  DOUBLE dx, dy, dz;
  DOUBLE q0, q_cos, q_sin;
  DOUBLE amplitude, angle, step, cosine, sine;
  DOUBLE step_cosine, step_sine;

#ifdef TRIAL_1D
  int particle;
  if(N < 1 || beads < 1 || T <= 0. || omega_z2 <= 0.) {
    Error("Harmonic normal-mode sampling requires N>=1, T>0, and omega_z>0\n");
  }

  tau = 1. / (T * (DOUBLE)beads);
  omega2 = omega_z2;
  for(bead=0; bead<beads; bead++) {
    for(particle=0; particle<N; particle++) {
      W[bead].x[particle] = 0.;
      W[bead].y[particle] = 0.;
      W[bead].z[particle] = 0.;
    }
  }

  for(particle=0; particle<N; particle++) {
    // Eigenvalues of the primitive cyclic action's quadratic form.
    eigenvalue = tau * omega2;
    RandomNormal3(&dx, &dy, &dz, 0., 1. / Sqrt(eigenvalue));
    q0 = dz;
    amplitude = q0 / Sqrt((DOUBLE)beads);
    for(bead=0; bead<beads; bead++) W[bead].z[particle] += amplitude;

    for(mode=1; mode<=(beads-1)/2; mode++) {
      step = 2. * PI * (DOUBLE)mode / (DOUBLE)beads;
      step_cosine = Cos(step);
      step_sine = Sin(step);
      eigenvalue = 4. * Sin(0.5 * step) * Sin(0.5 * step) / tau + tau * omega2;
      RandomNormal3(&dx, &dy, &dz, 0., 1. / Sqrt(eigenvalue));
      q_cos = dz;
      RandomNormal3(&dx, &dy, &dz, 0., 1. / Sqrt(eigenvalue));
      q_sin = dz;

      cosine = 1.;
      sine = 0.;
      amplitude = Sqrt(2. / (DOUBLE)beads);
      for(bead=0; bead<beads; bead++) {
        W[bead].z[particle] += amplitude * (q_cos * cosine + q_sin * sine);
        angle = cosine * step_cosine - sine * step_sine;
        sine = sine * step_cosine + cosine * step_sine;
        cosine = angle;
      }
    }

    if(beads % 2 == 0) {
      mode = beads / 2;
      eigenvalue = 4. / tau + tau * omega2;
      RandomNormal3(&dx, &dy, &dz, 0., 1. / Sqrt(eigenvalue));
      q_cos = dz / Sqrt((DOUBLE)beads);
      for(bead=0; bead<beads; bead++) {
        W[bead].z[particle] += (bead % 2 == 0) ? q_cos : -q_cos;
      }
    }
  }
#elif defined(TRIAL_2D)
  int particle;
  if(N < 1 || beads < 1 || T <= 0. || omega_x2 <= 0. || omega_y2 <= 0.) {
    Error("Harmonic normal-mode sampling requires N>=1, T>0, and omega_x>0, omega_y>0\n");
  }

  tau = 1. / (T * (DOUBLE)beads);

#ifdef WORM_ALGORITHM
  if(worm_next != NULL && worm_open == 0) {
    int *visited = (int*) calloc(N, sizeof(int));
    int *ring_p = (int*) malloc(N * beads * sizeof(int));
    int *ring_w = (int*) malloc(N * beads * sizeof(int));
    DOUBLE *ring_x = (DOUBLE*) malloc(N * beads * sizeof(DOUBLE));
    DOUBLE *ring_y = (DOUBLE*) malloc(N * beads * sizeof(DOUBLE));

    if(visited != NULL && ring_p != NULL && ring_w != NULL && ring_x != NULL && ring_y != NULL) {
      for(particle = 0; particle < N; particle++) {
        if(visited[particle]) continue;

        /* Rastrear el anillo continuo de cuentas a lo largo del ciclo físico */
        int ring_beads = 0;
        int curr_p = particle;
        int curr_w = 0;
        do {
          if(curr_w == 0) visited[curr_p] = 1;
          ring_p[ring_beads] = curr_p;
          ring_w[ring_beads] = curr_w;
          ring_x[ring_beads] = 0.;
          ring_y[ring_beads] = 0.;
          ring_beads++;
          if(ring_beads >= N * beads) break;
          curr_p = worm_next[curr_w][curr_p];
          curr_w = (curr_w + 1) % beads;
        } while(!(curr_w == 0 && curr_p == particle));

        /* --- Muestreo X para el anillo de ring_beads --- */
        eigenvalue = tau * omega_x2;
        RandomNormal3(&dx, &dy, &dz, 0., 1. / Sqrt(eigenvalue));
        q0 = dx;
        amplitude = q0 / Sqrt((DOUBLE)ring_beads);
        for(int n=0; n<ring_beads; n++) ring_x[n] += amplitude;

        for(mode=1; mode<=(ring_beads-1)/2; mode++) {
          step = 2. * PI * (DOUBLE)mode / (DOUBLE)ring_beads;
          step_cosine = Cos(step);
          step_sine = Sin(step);
          eigenvalue = 4. * Sin(0.5 * step) * Sin(0.5 * step) / tau + tau * omega_x2;
          RandomNormal3(&dx, &dy, &dz, 0., 1. / Sqrt(eigenvalue));
          q_cos = dx;
          RandomNormal3(&dx, &dy, &dz, 0., 1. / Sqrt(eigenvalue));
          q_sin = dx;

          cosine = 1.;
          sine = 0.;
          amplitude = Sqrt(2. / (DOUBLE)ring_beads);
          for(int n=0; n<ring_beads; n++) {
            ring_x[n] += amplitude * (q_cos * cosine + q_sin * sine);
            angle = cosine * step_cosine - sine * step_sine;
            sine = sine * step_cosine + cosine * step_sine;
            cosine = angle;
          }
        }

        if(ring_beads % 2 == 0) {
          eigenvalue = 4. / tau + tau * omega_x2;
          RandomNormal3(&dx, &dy, &dz, 0., 1. / Sqrt(eigenvalue));
          q_cos = dx / Sqrt((DOUBLE)ring_beads);
          for(int n=0; n<ring_beads; n++) {
            ring_x[n] += (n % 2 == 0) ? q_cos : -q_cos;
          }
        }

        /* --- Muestreo Y para el anillo de ring_beads --- */
        eigenvalue = tau * omega_y2;
        RandomNormal3(&dx, &dy, &dz, 0., 1. / Sqrt(eigenvalue));
        q0 = dy;
        amplitude = q0 / Sqrt((DOUBLE)ring_beads);
        for(int n=0; n<ring_beads; n++) ring_y[n] += amplitude;

        for(mode=1; mode<=(ring_beads-1)/2; mode++) {
          step = 2. * PI * (DOUBLE)mode / (DOUBLE)ring_beads;
          step_cosine = Cos(step);
          step_sine = Sin(step);
          eigenvalue = 4. * Sin(0.5 * step) * Sin(0.5 * step) / tau + tau * omega_y2;
          RandomNormal3(&dx, &dy, &dz, 0., 1. / Sqrt(eigenvalue));
          q_cos = dy;
          RandomNormal3(&dx, &dy, &dz, 0., 1. / Sqrt(eigenvalue));
          q_sin = dy;

          cosine = 1.;
          sine = 0.;
          amplitude = Sqrt(2. / (DOUBLE)ring_beads);
          for(int n=0; n<ring_beads; n++) {
            ring_y[n] += amplitude * (q_cos * cosine + q_sin * sine);
            angle = cosine * step_cosine - sine * step_sine;
            sine = sine * step_cosine + cosine * step_sine;
            cosine = angle;
          }
        }

        if(ring_beads % 2 == 0) {
          eigenvalue = 4. / tau + tau * omega_y2;
          RandomNormal3(&dx, &dy, &dz, 0., 1. / Sqrt(eigenvalue));
          q_cos = dy / Sqrt((DOUBLE)ring_beads);
          for(int n=0; n<ring_beads; n++) {
            ring_y[n] += (n % 2 == 0) ? q_cos : -q_cos;
          }
        }

        /* Asignar las coordenadas calculadas a las cuentas físicas */
        for(int n=0; n<ring_beads; n++) {
          W[ring_w[n]].x[ring_p[n]] = ring_x[n];
          W[ring_w[n]].y[ring_p[n]] = ring_y[n];
          W[ring_w[n]].z[ring_p[n]] = 0.;
        }
      }

      free(visited);
      free(ring_p);
      free(ring_w);
      free(ring_x);
      free(ring_y);
      return;
    }
    if(visited) free(visited);
    if(ring_p) free(ring_p);
    if(ring_w) free(ring_w);
    if(ring_x) free(ring_x);
    if(ring_y) free(ring_y);
  }
#endif

  for(bead=0; bead<beads; bead++) {
    for(particle=0; particle<N; particle++) {
      W[bead].x[particle] = 0.;
      W[bead].y[particle] = 0.;
      W[bead].z[particle] = 0.;
    }
  }

  for(particle=0; particle<N; particle++) {
    // --- Muestreo X ---
    eigenvalue = tau * omega_x2;
    RandomNormal3(&dx, &dy, &dz, 0., 1. / Sqrt(eigenvalue));
    q0 = dx;
    amplitude = q0 / Sqrt((DOUBLE)beads);
    for(bead=0; bead<beads; bead++) W[bead].x[particle] += amplitude;

    for(mode=1; mode<=(beads-1)/2; mode++) {
      step = 2. * PI * (DOUBLE)mode / (DOUBLE)beads;
      step_cosine = Cos(step);
      step_sine = Sin(step);
      eigenvalue = 4. * Sin(0.5 * step) * Sin(0.5 * step) / tau + tau * omega_x2;
      RandomNormal3(&dx, &dy, &dz, 0., 1. / Sqrt(eigenvalue));
      q_cos = dx;
      RandomNormal3(&dx, &dy, &dz, 0., 1. / Sqrt(eigenvalue));
      q_sin = dx;

      cosine = 1.;
      sine = 0.;
      amplitude = Sqrt(2. / (DOUBLE)beads);
      for(bead=0; bead<beads; bead++) {
        W[bead].x[particle] += amplitude * (q_cos * cosine + q_sin * sine);
        angle = cosine * step_cosine - sine * step_sine;
        sine = sine * step_cosine + cosine * step_sine;
        cosine = angle;
      }
    }

    if(beads % 2 == 0) {
      eigenvalue = 4. / tau + tau * omega_x2;
      RandomNormal3(&dx, &dy, &dz, 0., 1. / Sqrt(eigenvalue));
      q_cos = dx / Sqrt((DOUBLE)beads);
      for(bead=0; bead<beads; bead++) {
        W[bead].x[particle] += (bead % 2 == 0) ? q_cos : -q_cos;
      }
    }

    // --- Muestreo Y ---
    eigenvalue = tau * omega_y2;
    RandomNormal3(&dx, &dy, &dz, 0., 1. / Sqrt(eigenvalue));
    q0 = dy;
    amplitude = q0 / Sqrt((DOUBLE)beads);
    for(bead=0; bead<beads; bead++) W[bead].y[particle] += amplitude;

    for(mode=1; mode<=(beads-1)/2; mode++) {
      step = 2. * PI * (DOUBLE)mode / (DOUBLE)beads;
      step_cosine = Cos(step);
      step_sine = Sin(step);
      eigenvalue = 4. * Sin(0.5 * step) * Sin(0.5 * step) / tau + tau * omega_y2;
      RandomNormal3(&dx, &dy, &dz, 0., 1. / Sqrt(eigenvalue));
      q_cos = dy;
      RandomNormal3(&dx, &dy, &dz, 0., 1. / Sqrt(eigenvalue));
      q_sin = dy;

      cosine = 1.;
      sine = 0.;
      amplitude = Sqrt(2. / (DOUBLE)beads);
      for(bead=0; bead<beads; bead++) {
        W[bead].y[particle] += amplitude * (q_cos * cosine + q_sin * sine);
        angle = cosine * step_cosine - sine * step_sine;
        sine = sine * step_cosine + cosine * step_sine;
        cosine = angle;
      }
    }

    if(beads % 2 == 0) {
      eigenvalue = 4. / tau + tau * omega_y2;
      RandomNormal3(&dx, &dy, &dz, 0., 1. / Sqrt(eigenvalue));
      q_cos = dy / Sqrt((DOUBLE)beads);
      for(bead=0; bead<beads; bead++) {
        W[bead].y[particle] += (bead % 2 == 0) ? q_cos : -q_cos;
      }
    }
  }
#endif
}
#endif

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
  DOUBLE sigma;

  w1 = w-1; // previous bead number
  w2 = w+1; // next bead number
  if(w1<0) w1 = Nwalkers-1; // cyclic condition
  if(w2==Nwalkers) w2 = 0;

  tau = 1./(T*(DOUBLE)Nwalkers); // inverse temperature on a single bead (epsilon)
  sigma = 0.5*tau;

  // Move particle i of the bead
  for(i=0; i<N; i++) { // loop over particles
    du = 0.;

    x = W[w].x[i];
    y = W[w].y[i];
    z = W[w].z[i];
#ifdef WORM_ALGORITHM
    int prev_p = (worm_prev != NULL) ? worm_prev[w][i] : i;
    int next_p = (worm_next != NULL) ? worm_next[w][i] : i;
    if(prev_p < 0) prev_p = i;
    if(next_p < 0) next_p = i;
#else
    int prev_p = i;
    int next_p = i;
#endif
    xm = 0.5*(W[w1].x[prev_p]+W[w2].x[next_p]); // middle point
    ym = 0.5*(W[w1].y[prev_p]+W[w2].y[next_p]);
    zm = 0.5*(W[w1].z[prev_p]+W[w2].z[next_p]);

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
    du -= r2/(2.*sigma);

    // old kinetic energy
    dr[0] = x - xm;
    dr[1] = y - ym;
    dr[2] = z - zm;
    r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
    du += r2/(2.*sigma);

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


#ifdef WORM_ALGORITHM
/***********************************************************************
 *  Algoritmo Worm (Prokof'ev-Svistunov / Dornheim 2014)
 *
 *  Movimientos del Worm en el sector G
 *
 *  - Enlaces:
 *    worm_next[w][i]: índice de partícula en el bead (w+1)%Nwalkers.
 *    worm_prev[w][i]: índice de partícula en el bead (w-1+Nwalkers)%Nwalkers.
 *    Si el enlace está abierto (Ira o Masha), vale -1.
 *
 *  - Sector Z (worm_open == 0):
 *    Todos los enlaces están cerrados.
 *    Movimiento: Open (corta un enlace aleatorio, creando Ira y Masha con gap = 1).
 *
 *  - Sector G (worm_open == 1):
 *    Ira en (worm_bead, worm_particle), Masha en (worm_bead_masha, worm_particle_masha).
 *    Gap: delta_w = (worm_bead_masha - worm_bead + Nwalkers) % Nwalkers.
 *    Movimientos:
 *      Advance/Recede cambian el gap en +/-1; Close reconecta cuando gap == 1.
 *      Para N>1 también se proponen intercambios Swap antes de Masha.
 *      IraSpatialMove desplaza espacialmente la cabeza Ira.
 ***********************************************************************/

static int worm_snapshot_seen[6] = {0, 0, 0, 0, 0, 0};
static const char *worm_snapshot_names[6] = {
  "Open", "Close", "Advance", "Recede", "Swap", "IraSpatialMove"
};

static int *worm_swap_candidates = NULL;
static DOUBLE *worm_swap_weights = NULL;

void FreeWorm(void) {
  int w;
  if(worm_swap_candidates != NULL) {
    free(worm_swap_candidates);
    worm_swap_candidates = NULL;
  }
  if(worm_swap_weights != NULL) {
    free(worm_swap_weights);
    worm_swap_weights = NULL;
  }
  if(worm_gap_histogram != NULL) {
    free(worm_gap_histogram);
    worm_gap_histogram = NULL;
  }
  if(worm_cycle_histogram != NULL) {
    free(worm_cycle_histogram);
    worm_cycle_histogram = NULL;
  }
  if(worm_next != NULL) {
    for(w = 0; w < NwalkersMax; w++) {
      if(worm_next[w]) free(worm_next[w]);
      if(worm_prev[w]) free(worm_prev[w]);
    }
    free(worm_next);
    free(worm_prev);
    worm_next = NULL;
    worm_prev = NULL;
  }
}

void InitWorm(void) {
  int w, i;
  if(worm_next != NULL) return;

  atexit(FreeWorm);

  {
    FILE *snapshot_file = fopen("outworm_snapshots.dat", "w");
    int move;
    for(move = 0; move < 6; move++) worm_snapshot_seen[move] = 0;
    if(snapshot_file != NULL) {
      fprintf(snapshot_file, "# worm_snapshot_v1\n");
      fclose(snapshot_file);
    }
  }

  worm_next = (int**) malloc(NwalkersMax * sizeof(int*));
  worm_prev = (int**) malloc(NwalkersMax * sizeof(int*));
  for(w = 0; w < NwalkersMax; w++) {
    worm_next[w] = (int*) malloc(N * sizeof(int));
    worm_prev[w] = (int*) malloc(N * sizeof(int));
    for(i = 0; i < N; i++) {
      worm_next[w][i] = i;
      worm_prev[w][i] = i;
    }
  }

  worm_swap_candidates = (int*) malloc(N * sizeof(int));
  worm_swap_weights = (DOUBLE*) malloc(N * sizeof(DOUBLE));

  worm_gap_histogram = (long long*) calloc(NwalkersMax + 1, sizeof(long long));
  worm_max_gap_observed = 0;

  worm_cycle_histogram = (long long*) calloc(N + 1, sizeof(long long));
  worm_cycle_samples = 0;
}

static int WormSnapshotWrite(int move, const char *phase) {
  FILE *snapshot_file;
  int bead, particle;

  snapshot_file = fopen("outworm_snapshots.dat", "a");
  if(snapshot_file == NULL) return 0;

  fprintf(snapshot_file, "EVENT %s %s %d %d %d %d %d %d %d\n",
          worm_snapshot_names[move], phase, Nwalkers, N, worm_open,
          worm_bead, worm_particle, worm_bead_masha, worm_particle_masha);
  for(bead = 0; bead < Nwalkers; bead++) {
    for(particle = 0; particle < N; particle++) {
      fprintf(snapshot_file, "BEAD %d %d %.17g %.17g %d %d\n",
              bead, particle, (double)W[bead].x[particle], (double)W[bead].y[particle],
              worm_next[bead][particle], worm_prev[bead][particle]);
    }
  }
  fprintf(snapshot_file, "END\n");
  fclose(snapshot_file);
  return 1;
}

static int WormSnapshotBefore(int move) {
  if(move < 0 || move >= 6 || worm_snapshot_seen[move]) return 0;
  return WormSnapshotWrite(move, "before");
}

static void WormSnapshotAfter(int move) {
  if(WormSnapshotWrite(move, "after")) worm_snapshot_seen[move] = 1;
}

#define WORM_MAX_GAP 6
#define WORM_PROB_CLOSE_GAP1     0.60
#define WORM_PROB_RECEDE_GAP1    0.25
#define WORM_PROB_SPATIAL_GAP1   0.15
#define WORM_PROB_ADVANCE_GAP_GT1 0.45
#define WORM_PROB_RECEDE_GAP_GT1  0.45
#define WORM_PROB_SPATIAL_GAP_GT1 0.10

#define WORM_PROB_CLOSE_N_GT1    0.40
#define WORM_PROB_SWAP_N_GT1     0.40
#define WORM_PROB_SPATIAL_N_GT1  0.20
#define WORM_PROB_ADVANCE_N_GT1  0.20
#define WORM_PROB_RECEDE_N_GT1   0.20
#define WORM_PROB_RECEDE_N_GAP1  0.20
#define WORM_PROB_SWAP_N_GAP1    0.20
#define WORM_PROB_SPATIAL_N_GAP1 0.20

static void WormOpenMove(void) {
  int w_cut, i_cut, j_cut, w_next;
  DOUBLE dx, dy, r2, tau_local, acc, p_close;

  worm_try_open++;

  w_cut = (int)(Random() * Nwalkers);
  if(w_cut >= Nwalkers) w_cut = Nwalkers - 1;
  i_cut = (int)(Random() * N);
  if(i_cut >= N) i_cut = N - 1;

  w_next = (w_cut + 1) % Nwalkers;
  j_cut = worm_next[w_cut][i_cut];
  if(j_cut < 0) return; /* ya abierto */

  dx = W[w_next].x[j_cut] - W[w_cut].x[i_cut];
  dy = W[w_next].y[j_cut] - W[w_cut].y[i_cut];
  r2 = dx * dx + dy * dy;

  tau_local = 1.0 / (T * (DOUBLE)Nwalkers);
  p_close = (N > 1) ? WORM_PROB_CLOSE_N_GT1 : WORM_PROB_CLOSE_GAP1;

  /* Aceptación Metropolis: ganancia del resorte cinético al cortar */
  acc = worm_C * Exp(r2 / (2.0 * tau_local))
    * p_close * (DOUBLE)Nwalkers * (DOUBLE)N;

  if(acc >= 1.0 || Random() < acc) {
    int record_snapshot = WormSnapshotBefore(0);
    worm_open = 1;
    worm_bead = w_cut;
    worm_particle = i_cut;
    worm.w = w_cut;
    worm.i = i_cut;
    worm.x_ira = W[w_cut].x[i_cut];
    worm.y_ira = W[w_cut].y[i_cut];

    worm_bead_masha = w_next;
    worm_particle_masha = j_cut;
    worm.x_masha = W[w_next].x[j_cut];
    worm.y_masha = W[w_next].y[j_cut];

    worm_next[w_cut][i_cut] = -1;
    worm_prev[w_next][j_cut] = -1;

    worm_acc_open++;
    if(record_snapshot) WormSnapshotAfter(0);
  }
}

static void WormCloseMove(void) {
  int w_masha_expected;
  DOUBLE dx, dy, r2, tau_local, acc, p_close;

  worm_try_close++;

  w_masha_expected = (worm_bead + 1) % Nwalkers;
  if(worm_bead_masha != w_masha_expected) {
    return; /* solo se puede cerrar si los beads son contiguos (gap == 1) */
  }

  dx = W[worm_bead_masha].x[worm_particle_masha] - W[worm_bead].x[worm_particle];
  dy = W[worm_bead_masha].y[worm_particle_masha] - W[worm_bead].y[worm_particle];
  r2 = dx * dx + dy * dy;

  tau_local = 1.0 / (T * (DOUBLE)Nwalkers);
  p_close = (N > 1) ? WORM_PROB_CLOSE_N_GT1 : WORM_PROB_CLOSE_GAP1;

  /* Aceptación Metropolis: coste del resorte cinético al reconectar */
  acc = (1.0 / worm_C) * Exp(-r2 / (2.0 * tau_local))
    / (p_close * (DOUBLE)Nwalkers * (DOUBLE)N);

  if(acc >= 1.0 || Random() < acc) {
    int record_snapshot = WormSnapshotBefore(1);
    worm_next[worm_bead][worm_particle] = worm_particle_masha;
    worm_prev[worm_bead_masha][worm_particle_masha] = worm_particle;

    worm_open = 0;
    worm_bead = -1;
    worm_particle = -1;
    worm_bead_masha = -1;
    worm_particle_masha = -1;

    worm_acc_close++;
    if(record_snapshot) WormSnapshotAfter(1);
  }
}

static void WormAdvanceMove(void) {
  int w_curr, w_next, p_idx, p_next, gap, particle, candidate_count;
  DOUBLE dx, dy, dummy_z, x_new, y_new, dV, acc, tau_local, proposal_ratio;
  DOUBLE p_advance, p_recede_reverse;

  worm_try_advance++;

  w_curr = worm_bead;
  p_idx = worm_particle;
  w_next = (w_curr + 1) % Nwalkers;

  /* No avanzar si ya estamos adyacentes a Masha (eso corresponde a Close) */
  if(w_next == worm_bead_masha) return;

  gap = (worm_bead_masha - worm_bead + Nwalkers) % Nwalkers;

  /* Recede leaves one detached bead on the slice immediately after Ira. */
  p_next = -1;
  candidate_count = 0;
  for(particle = 0; particle < N; particle++) {
    if(worm_prev[w_next][particle] < 0 && worm_next[w_next][particle] < 0) {
      p_next = particle;
      candidate_count++;
    }
  }
  if(candidate_count != 1) return;

  tau_local = 1.0 / (T * (DOUBLE)Nwalkers);

  /* Muestreo gaussiano del propagador cinético libre desde Ira */
  RandomNormal3(&dx, &dy, &dummy_z, 0., Sqrt(tau_local));
  x_new = W[w_curr].x[p_idx] + dx;
  y_new = W[w_curr].y[p_idx] + dy;

  /* Cambio en la acción potencial en el nuevo bead */
  dV = Vext(x_new, y_new, 0.);
  p_advance = (N > 1) ? WORM_PROB_ADVANCE_N_GT1 : WORM_PROB_ADVANCE_GAP_GT1;
  if(gap == 2) {
    p_recede_reverse = (N > 1) ? WORM_PROB_RECEDE_N_GAP1 : WORM_PROB_RECEDE_GAP1;
  } else {
    p_recede_reverse = (N > 1) ? WORM_PROB_RECEDE_N_GT1 : WORM_PROB_RECEDE_GAP_GT1;
  }
  proposal_ratio = p_recede_reverse / p_advance;
  acc = Exp(-tau_local * dV) * proposal_ratio;

  if(acc >= 1.0 || Random() < acc) {
    int record_snapshot = WormSnapshotBefore(2);
    W[w_next].x[p_next] = x_new;
    W[w_next].y[p_next] = y_new;

    worm_next[w_curr][p_idx] = p_next;
    worm_prev[w_next][p_next] = p_idx;

    worm_bead = w_next;
    worm_particle = p_next;
    worm.w = w_next;
    worm.i = p_next;
    worm.x_ira = x_new;
    worm.y_ira = y_new;

    worm_next[w_next][p_next] = -1;

    worm_acc_advance++;
    if(record_snapshot) WormSnapshotAfter(2);
  }
}

static void WormRecedeMove(void) {
  int w_curr, w_prev, p_idx, p_prev, gap;
  DOUBLE dV, acc, tau_local, proposal_ratio;
  DOUBLE p_recede, p_advance_reverse;

  worm_try_recede++;

  w_curr = worm_bead;
  p_idx = worm_particle;
  w_prev = (w_curr - 1 + Nwalkers) % Nwalkers;

  gap = (worm_bead_masha - worm_bead + Nwalkers) % Nwalkers;
    p_prev = worm_prev[w_curr][p_idx];
  /* Límite máximo de apertura (M_bar en Dornheim) y protección contra colisión con Masha */
    if(gap >= WORM_MAX_GAP || w_prev == worm_bead_masha || p_prev < 0
      || worm_next[w_prev][p_prev] != p_idx) return;

  tau_local = 1.0 / (T * (DOUBLE)Nwalkers);

  /* Al borrar el bead de Ira, se elimina su potencial */
  dV = Vext(W[w_curr].x[p_idx], W[w_curr].y[p_idx], 0.);
  p_advance_reverse = (N > 1) ? WORM_PROB_ADVANCE_N_GT1 : WORM_PROB_ADVANCE_GAP_GT1;
  if(gap == 1) {
    p_recede = (N > 1) ? WORM_PROB_RECEDE_N_GAP1 : WORM_PROB_RECEDE_GAP1;
  } else {
    p_recede = (N > 1) ? WORM_PROB_RECEDE_N_GT1 : WORM_PROB_RECEDE_GAP_GT1;
  }
  proposal_ratio = p_advance_reverse / p_recede;
  acc = Exp(+tau_local * dV) * proposal_ratio;

  if(acc >= 1.0 || Random() < acc) {
    int record_snapshot = WormSnapshotBefore(3);
    worm_next[w_prev][p_prev] = -1;
    worm_prev[w_curr][p_idx] = -1;

    worm_bead = w_prev;
    worm_particle = p_prev;
    worm.w = w_prev;
    worm.i = p_prev;
    worm.x_ira = W[w_prev].x[p_prev];
    worm.y_ira = W[w_prev].y[p_prev];

    worm_acc_recede++;
    if(record_snapshot) WormSnapshotAfter(3);
  }
}

static void WormSwapMove(void) {
  int w_curr, w_next, i_ira, k, p_old, cand_count;
  int *candidates = worm_swap_candidates;
  DOUBLE *weights = worm_swap_weights;
  DOUBLE sum_weights1, sum_weights2, r, cum, acc;
  DOUBLE tau_local, d2;
  int chosen_idx, chosen_k;

  worm_try_swap++;

  if(N < 2) return; /* Se requieren al menos 2 partículas para intercambio */
  if(candidates == NULL || weights == NULL) return;

  w_curr = worm_bead;
  i_ira = worm_particle;
  w_next = (w_curr + 1) % Nwalkers;
  tau_local = 1.0 / (T * (DOUBLE)Nwalkers);

  /* 1. Candidatos en slice w_next con enlace entrante válido (excluye Masha si está en w_next) */
  cand_count = 0;
  for(k = 0; k < N; k++) {
    if(worm_prev[w_next][k] >= 0) {
      candidates[cand_count++] = k;
    }
  }

  if(cand_count == 0) return;

  /* 2. Calcular chi_1(k) = exp(-|r_ira - r_{k, next}|^2 / (2*tau)) */
  sum_weights1 = 0.0;
  for(int c = 0; c < cand_count; c++) {
    k = candidates[c];
    DOUBLE dx = worm.x_ira - W[w_next].x[k];
    DOUBLE dy = worm.y_ira - W[w_next].y[k];
    d2 = dx * dx + dy * dy;
    weights[c] = Exp(-d2 / (2.0 * tau_local));
    sum_weights1 += weights[c];
  }

  if(sum_weights1 <= 0.0) return;

  /* 3. Elegir candidato k con probabilidad P(k) = chi_1(k) / Sigma_1 */
  r = Random() * sum_weights1;
  cum = 0.0;
  chosen_idx = 0;
  for(int c = 0; c < cand_count; c++) {
    cum += weights[c];
    if(r <= cum) {
      chosen_idx = c;
      break;
    }
  }
  chosen_k = candidates[chosen_idx];
  p_old = worm_prev[w_next][chosen_k];

  /* 4. Calcular Sigma_2 para el movimiento inverso desde la nueva cabeza Ira (w_curr, p_old) */
  sum_weights2 = 0.0;
  for(int c = 0; c < cand_count; c++) {
    k = candidates[c];
    DOUBLE dx = W[w_curr].x[p_old] - W[w_next].x[k];
    DOUBLE dy = W[w_curr].y[p_old] - W[w_next].y[k];
    d2 = dx * dx + dy * dy;
    sum_weights2 += Exp(-d2 / (2.0 * tau_local));
  }

  if(sum_weights2 <= 0.0) return;

  /* 5. Aceptación Metropolis: A = min(1, Sigma_1 / Sigma_2) */
  acc = sum_weights1 / sum_weights2;

  if(acc >= 1.0 || Random() < acc) {
    int record_snapshot = WormSnapshotBefore(4);
    /* Reconectar enlaces */
    worm_next[w_curr][i_ira] = chosen_k;
    worm_prev[w_next][chosen_k] = i_ira;
    worm_next[w_curr][p_old] = -1;

    /* Actualizar cabeza Ira */
    worm_bead = w_curr;
    worm_particle = p_old;
    worm.w = w_curr;
    worm.i = p_old;
    worm.x_ira = W[w_curr].x[p_old];
    worm.y_ira = W[w_curr].y[p_old];

    worm_acc_swap++;
    if(record_snapshot) WormSnapshotAfter(4);
  }
}

static void WormIraSpatialMove(void) {
  int w_curr, w_prev, p_idx;
  DOUBLE dx, dy, dummy_z, x_new, y_new, dV_kin, dV_ext, dS, acc, tau_local;
  DOUBLE dx_old, dy_old, dx_new, dy_new;

  worm_try_iramove++;

  w_curr = worm_bead;
  p_idx = worm_particle;
  w_prev = (w_curr - 1 + Nwalkers) % Nwalkers;

  tau_local = 1.0 / (T * (DOUBLE)Nwalkers);

  /* Desplazamiento gaussiano local */
  RandomNormal3(&dx, &dy, &dummy_z, 0., Sqrt(dt_vmc > 0 ? dt_vmc : tau_local));
  x_new = W[w_curr].x[p_idx] + dx;
  y_new = W[w_curr].y[p_idx] + dy;

  /* Cambio cinético con el bead anterior */
  int p_in = (worm_prev != NULL && worm_prev[w_curr][p_idx] >= 0) ? worm_prev[w_curr][p_idx] : p_idx;
  dx_old = W[w_curr].x[p_idx] - W[w_prev].x[p_in];
  dy_old = W[w_curr].y[p_idx] - W[w_prev].y[p_in];
  dx_new = x_new - W[w_prev].x[p_in];
  dy_new = y_new - W[w_prev].y[p_in];

  dV_kin = ((dx_new*dx_new + dy_new*dy_new) - (dx_old*dx_old + dy_old*dy_old)) / (2.0 * tau_local);
  dV_ext = tau_local * (Vext(x_new, y_new, 0.) - Vext(W[w_curr].x[p_idx], W[w_curr].y[p_idx], 0.));

  dS = -(dV_kin + dV_ext);
  acc = Exp(dS);

  if(acc >= 1.0 || Random() < acc) {
    int record_snapshot = WormSnapshotBefore(5);
    W[w_curr].x[p_idx] = x_new;
    W[w_curr].y[p_idx] = y_new;
    worm.x_ira = x_new;
    worm.y_ira = y_new;
    worm_acc_iramove++;
    if(record_snapshot) WormSnapshotAfter(5);
  }
}

void PIMCWormMoves(void) {
  int sub, gap;
  DOUBLE r;

  if(worm_next == NULL) {
    InitWorm();
  }

  /* Realizar varios intentos por iteración para una buena exploración ergódica */
  for(sub = 0; sub < 5; sub++) {
    if(worm_open == 0) {
      worm_count_Z++;
      WormOpenMove();
    } else {
      worm_count_G++;

      gap = (worm_bead_masha - worm_bead + Nwalkers) % Nwalkers;
      if(gap < NwalkersMax && worm_gap_histogram != NULL) {
        worm_gap_histogram[gap]++;
      }
      if(gap > worm_max_gap_observed) worm_max_gap_observed = gap;

      r = Random();
      if(N > 1) {
        if(gap == 1) {
          if(r < WORM_PROB_CLOSE_N_GT1) {
            WormCloseMove();
          } else if(r < WORM_PROB_CLOSE_N_GT1 + WORM_PROB_RECEDE_N_GAP1) {
            WormRecedeMove();
          } else if(r < WORM_PROB_CLOSE_N_GT1 + WORM_PROB_RECEDE_N_GAP1
                    + WORM_PROB_SWAP_N_GAP1) {
            WormSwapMove();
          } else {
            WormIraSpatialMove();
          }
        } else {
          if(r < WORM_PROB_ADVANCE_N_GT1) {
            WormAdvanceMove();
          } else if(r < WORM_PROB_ADVANCE_N_GT1 + WORM_PROB_RECEDE_N_GT1) {
            WormRecedeMove();
          } else if(r < WORM_PROB_ADVANCE_N_GT1 + WORM_PROB_RECEDE_N_GT1
                    + WORM_PROB_SWAP_N_GT1) {
            WormSwapMove();
          } else {
            WormIraSpatialMove();
          }
        }
      } else {
        /* Para N == 1: Advance, Recede, Spatial, Close */
        if(gap == 1) {
          if(r < WORM_PROB_CLOSE_GAP1) {
            WormCloseMove();
          } else if(r < WORM_PROB_CLOSE_GAP1 + WORM_PROB_RECEDE_GAP1) {
            WormRecedeMove();
          } else {
            WormIraSpatialMove();
          }
        } else {
          if(r < WORM_PROB_ADVANCE_GAP_GT1) {
            WormAdvanceMove();
          } else if(r < WORM_PROB_ADVANCE_GAP_GT1 + WORM_PROB_RECEDE_GAP_GT1) {
            WormRecedeMove();
          } else {
            WormIraSpatialMove();
          }
        }
      }
    }
  }
}

void MeasureWormPermutations(void) {
  int p, curr_p, curr_w, L;
  int *visited;

  if(worm_open != 0 || worm_next == NULL) return;

  visited = (int*) calloc(N, sizeof(int));
  if(visited == NULL) return;

  worm_cycle_samples++;

  for(p = 0; p < N; p++) {
    if(!visited[p]) {
      L = 0;
      curr_p = p;
      curr_w = 0;
      do {
        if(curr_w == 0) {
          visited[curr_p] = 1;
          L++;
        }
        curr_p = worm_next[curr_w][curr_p];
        curr_w = (curr_w + 1) % Nwalkers;
        if(L > N) break;
      } while(!(curr_w == 0 && curr_p == p));

      if(worm_cycle_histogram != NULL && L <= N) {
        worm_cycle_histogram[L] += L; /* ponderado por número de partículas en ciclo L */
      }
    }
  }

  free(visited);
}

void SaveWormPermutations(const char *filename) {
  FILE *f;
  int L;
  DOUBLE total_particles;

  if(worm_cycle_histogram == NULL || worm_cycle_samples == 0) return;
  f = fopen(filename, "w");
  if(f == NULL) return;

  total_particles = (DOUBLE)worm_cycle_samples * (DOUBLE)N;

  fprintf(f, "# cycle_length particle_count fraction\n");
  for(L = 1; L <= N; L++) {
    DOUBLE frac = total_particles > 0. ? (DOUBLE)worm_cycle_histogram[L] / total_particles : 0.;
    fprintf(f, "%d %lld %.6f\n", L, worm_cycle_histogram[L], frac);
  }
  fclose(f);
}

void WormDiagnostic(void) {
  long long total = worm_count_Z + worm_count_G;
  if(total == 0) total = 1;
  double frac_G = (double)worm_count_G / (double)total;
  double rate_open = worm_try_open > 0 ? (double)worm_acc_open / (double)worm_try_open * 100.0 : 0.0;
  double rate_close = worm_try_close > 0 ? (double)worm_acc_close / (double)worm_try_close * 100.0 : 0.0;
  double rate_adv = worm_try_advance > 0 ? (double)worm_acc_advance / (double)worm_try_advance * 100.0 : 0.0;
  double rate_rec = worm_try_recede > 0 ? (double)worm_acc_recede / (double)worm_try_recede * 100.0 : 0.0;
  double rate_swap = worm_try_swap > 0 ? (double)worm_acc_swap / (double)worm_try_swap * 100.0 : 0.0;
  double rate_ira = worm_try_iramove > 0 ? (double)worm_acc_iramove / (double)worm_try_iramove * 100.0 : 0.0;

  double p1 = 0., p2 = 0.;
  if(worm_cycle_samples > 0 && worm_cycle_histogram != NULL) {
    double tot = (double)worm_cycle_samples * (double)N;
    if(tot > 0.) {
      p1 = (double)worm_cycle_histogram[1] / tot;
      if(N >= 2) p2 = (double)worm_cycle_histogram[2] / tot;
    }
  }

  Message("Worm F4: Z=%lld G=%lld (fG=%.2f) | Op: %.1f%% | Cl: %.1f%% | Adv: %.1f%% | Rec: %.1f%% | Sw: %.1f%% | Ira: %.1f%% | P(1)=%.2f P(2)=%.2f\n",
          worm_count_Z, worm_count_G, frac_G,
          rate_open, rate_close, rate_adv, rate_rec, rate_swap, rate_ira,
          p1, p2);
}

void SaveWormGapHistogram(const char *filename) {
  FILE *f;
  int g;
  if(worm_gap_histogram == NULL) return;
  f = fopen(filename, "w");
  if(f == NULL) return;
  fprintf(f, "# gap count\n");
  for(g = 1; g <= WORM_MAX_GAP && g <= NwalkersMax; g++) {
    fprintf(f, "%d %lld\n", g, worm_gap_histogram[g]);
  }
  fclose(f);
}
#endif /* WORM_ALGORITHM */
