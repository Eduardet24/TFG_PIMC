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
 *  Algoritmo Worm (Prokof'ev-Svistunov) para estadística de Bose.
 *
 *  FASE 2: Movimientos Open y Close (Transiciones Z <-> G)
 *
 *  - Enlaces:
 *    worm_next[w][i]: índice de partícula en el bead (w+1)%Nwalkers.
 *    worm_prev[w][i]: índice de partícula en el bead (w-1+Nwalkers)%Nwalkers.
 *    Si el enlace está abierto (Ira o Masha), vale -1.
 *
 *  - Sector Z:
 *    Todos los enlaces están cerrados.
 *    Se propone Open: se elige un bead w y partícula i al azar,
 *    se corta el enlace hacia (w+1)%Nwalkers.
 *    Aceptación: acc = worm_C * exp(+r^2 / (2*tau))
 *
 *  - Sector G:
 *    El gusano está abierto con cabeza Ira en (worm_bead, worm_particle)
 *    y cola Masha en (worm_bead_masha, worm_particle_masha).
 *    Se propone Close: se verifica que Ira y Masha sean contiguos en tiempo.
 *    Aceptación: acc = (1.0 / worm_C) * exp(-r^2 / (2*tau))
 ***********************************************************************/

void InitWorm(void) {
  int w, i;
  if(worm_next != NULL) return;

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
}

static void WormOpenMove(void) {
  int w_cut, i_cut, j_cut, w_next;
  DOUBLE dx, dy, r2, tau, acc;

  worm_try_open++;

  /* Elegir aleatoriamente el bead y la partícula del enlace a abrir */
  w_cut = (int)(Random() * Nwalkers);
  if(w_cut >= Nwalkers) w_cut = Nwalkers - 1;
  i_cut = (int)(Random() * N);
  if(i_cut >= N) i_cut = N - 1;

  w_next = (w_cut + 1) % Nwalkers;
  j_cut = worm_next[w_cut][i_cut];
  if(j_cut < 0) return; /* ya abierto */

  /* Distancia espacial en 2D entre los dos extremos del enlace */
  dx = W[w_next].x[j_cut] - W[w_cut].x[i_cut];
  dy = W[w_next].y[j_cut] - W[w_cut].y[i_cut];
  r2 = dx * dx + dy * dy;

  tau = 1.0 / (T * (DOUBLE)Nwalkers);

  /* Metropolis: quitar el resorte cinético da un factor exp(+r^2 / (2*tau)) */
  acc = worm_C * Exp(r2 / (2.0 * tau));

  if(acc >= 1.0 || Random() < acc) {
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

    /* Abrir enlace */
    worm_next[w_cut][i_cut] = -1;
    worm_prev[w_next][j_cut] = -1;

    worm_acc_open++;
  }
}

static void WormCloseMove(void) {
  int w_masha_expected;
  DOUBLE dx, dy, r2, tau, acc;

  worm_try_close++;

  /* Condición para cerrar: Ira y Masha deben ser contiguos en tiempo */
  w_masha_expected = (worm_bead + 1) % Nwalkers;
  if(worm_bead_masha != w_masha_expected) {
    return;
  }

  dx = W[worm_bead_masha].x[worm_particle_masha] - W[worm_bead].x[worm_particle];
  dy = W[worm_bead_masha].y[worm_particle_masha] - W[worm_bead].y[worm_particle];
  r2 = dx * dx + dy * dy;

  tau = 1.0 / (T * (DOUBLE)Nwalkers);

  /* Metropolis: añadir el resorte cinético da un factor exp(-r^2 / (2*tau)) */
  acc = (1.0 / worm_C) * Exp(-r2 / (2.0 * tau));

  if(acc >= 1.0 || Random() < acc) {
    /* Reconectar enlace y cerrar gusano */
    worm_next[worm_bead][worm_particle] = worm_particle_masha;
    worm_prev[worm_bead_masha][worm_particle_masha] = worm_particle;

    worm_open = 0;
    worm_bead = -1;
    worm_particle = -1;
    worm_bead_masha = -1;
    worm_particle_masha = -1;

    worm_acc_close++;
  }
}

void PIMCWormMoves(void) {
  if(worm_next == NULL) {
    InitWorm();
  }

  if(worm_open == 0) {
    worm_count_Z++;
    WormOpenMove();
  } else {
    worm_count_G++;
    WormCloseMove();
  }
}

void WormDiagnostic(void) {
  long long total = worm_count_Z + worm_count_G;
  if(total == 0) total = 1;
  double frac_G = (double)worm_count_G / (double)total;
  double rate_open = worm_try_open > 0 ? (double)worm_acc_open / (double)worm_try_open * 100.0 : 0.0;
  double rate_close = worm_try_close > 0 ? (double)worm_acc_close / (double)worm_try_close * 100.0 : 0.0;

  Message("Worm F2: Z=%lld G=%lld (frac_G=%.3f) | Open: %lld/%lld (%.1f%%) | Close: %lld/%lld (%.1f%%) [C=%.2f]\n",
          worm_count_Z, worm_count_G, frac_G,
          worm_acc_open, worm_try_open, rate_open,
          worm_acc_close, worm_try_close, rate_close,
          worm_C);
}
#endif /* WORM_ALGORITHM */
