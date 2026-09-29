/*pigs.c*/
// contains path-integral implementation of Path Integram Monte Carlo and Path Integral Ground State methods

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

struct Worm worm;

#ifdef ONE_COMPONENT_CODE
  DOUBLE mA = 1.;
#endif

/***************************** Path Integral Move All ************************/
void PathIntegralMoveAll(void) {
// Moves a single particle in all beads
  int iloop,i,w,up;
  DOUBLE xi;
  DOUBLE dx,dy,dz;
  DOUBLE du;
  DOUBLE sigma;

  DOUBLE Sold, Snew, dS;
  DOUBLE dSkin, Skinold, Skin;
  DOUBLE dSdelta, Sdeltaold, Sdelta;
  DOUBLE dSwf, Swfold, Swf;

  // Move particle i of the bead
#ifdef ONE_COMPONENT_CODE
  for(iloop=0; iloop<N; iloop++) {
#else // two component code
#ifdef MOVE_ONLY_UP_PARTICLES
  for(iloop=0; iloop<Nup; iloop++) { // only particles up are moved
#else // move all particles
  for(iloop=0; iloop<Nup+Ndn; iloop++) { // move i-th particle to the point rm
#endif
#endif

    dSkin = dSdelta = dSwf = 0.;

#ifdef ONE_COMPONENT_CODE
    i = iloop; // index of the moved particle
#else
    if(iloop<Nup) { // move particle up
      up = ON;
      i = iloop;
    }
    else { // move particle down
      up = OFF;
      i = iloop-Nup;
    }
#endif

    Sold = PathIntegralPseudopotentialSdetailed(&Skinold, &Sdeltaold, &Swfold); // test change

#ifdef SECURE
   if(isnan(Sold)) Error(" NaN encountered is Sold\n");
#endif

    sigma = sqrt(0.5*dt_all);
    RandomNormal3(&dx, &dy, &dz, 0., sigma);

    //if(CheckCoordinatesDimensionality()) Message("dimensionality check not passed");

    //  Move particle np of the walker w
#ifdef ONE_COMPONENT_CODE
    for(w=0; w<Npop; w++) {
      CaseX(Wp[w].x[i] = W[w].x[i];) // store old coordinates
      CaseY(Wp[w].y[i] = W[w].y[i];)
      CaseZ(Wp[w].z[i] = W[w].z[i];)
      CaseX(W[w].x[i] += dx;) // move
      CaseY(W[w].y[i] += dy;)
      CaseZ(W[w].z[i] += dz;)
      ReduceToTheBoxXYZ(&W[w].x[i], &W[w].y[i], &W[w].z[i]);
    }
#else
    if(up) {
      for(w=0; w<Npop; w++) {
        CaseX(Wp[w].x[i] = W[w].x[i];) // store old coordinates
        CaseY(Wp[w].y[i] = W[w].y[i];)
        CaseZ(Wp[w].z[i] = W[w].z[i];)
        CaseX(W[w].x[i] += m_up_inv_sqrt*dx;) // move
        CaseY(W[w].y[i] += m_up_inv_sqrt*dy;)
        CaseZ(W[w].z[i] += m_up_inv_sqrt*dz;)
        ReduceToTheBoxXYZ(&W[w].x[i], &W[w].y[i], &W[w].z[i]);
      }
    }
    else {
      for(w=0; w<Npop; w++) {
        CaseX(Wp[w].xdn[i] = W[w].xdn[i];)
        CaseY(Wp[w].ydn[i] = W[w].ydn[i];)
        CaseZ(Wp[w].zdn[i] = W[w].zdn[i];)
        CaseX(W[w].xdn[i] += m_dn_inv_sqrt*dx;)
        CaseY(W[w].ydn[i] += m_dn_inv_sqrt*dy;)
        CaseZ(W[w].zdn[i] += m_dn_inv_sqrt*dz;)
        ReduceToTheBoxXYZ(&W[w].xdn[i], &W[w].ydn[i], &W[w].zdn[i]);
      }
    }
#endif

#ifdef CENTER_OF_MASS_IS_NOT_MOVED // propose such movement that CM position is not moved
    Error("  PIGS Move All does not support CENTER_OF_MASS_IS_NOT_MOVED\n");
#endif

    Snew = PathIntegralPseudopotentialSdetailed(&Skin, &Sdelta, &Swf);
    dS = Snew - Sold ; // test change
    dSkin = Skin - Skinold;
    dSdelta = Sdelta - Sdeltaold;
    dSwf = Swf - Swfold;

    if(Random() < Exp(dS)) { // The Metropolis code
      accepted_all++;
    }
    else {
      rejected_all++;
#ifdef ONE_COMPONENT_CODE
      for(w=0; w<Npop; w++) {
        CaseX(W[w].x[i] = Wp[w].x[i];) // restore old coordinates
        CaseY(W[w].y[i] = Wp[w].y[i];)
        CaseZ(W[w].z[i] = Wp[w].z[i];)
      }
#else
      if(up) {
        for(w=0; w<Npop; w++) {
          CaseX(W[w].x[i] = Wp[w].x[i];) // restore old coordinates
          CaseY(W[w].y[i] = Wp[w].y[i];)
          CaseZ(W[w].z[i] = Wp[w].z[i];)
        }
      }
      else {
        for(w=0; w<Npop; w++) {
          CaseX(W[w].xdn[i] = Wp[w].xdn[i];)
          CaseY(W[w].ydn[i] = Wp[w].ydn[i];)
          CaseZ(W[w].zdn[i] = Wp[w].zdn[i];)
        }
      }
#endif
    }
  }
}

/************************** PIMC Moce One By One ****************************/
double PropagatorPseudopotenial1D_ij(double m, double a, double dt, double xi0, double xj0, double xi, double xj) {
// dt - propagation time
// xij0 = x_i(0) - x_j(0), no absolute value
// xij = x_i(dt) - x_j(dt), no absolute value
// m = 2\mu
  double safe;
  double x;
  double xij,xij0,yij;
  double absxij, absxij0;
  double yerfcx, yexp;
  int particles_exchange_sign;

  xij0 = xi0 - xj0;
  xij = xi - xj;

#ifdef BC_1DPBC
  if(xij0> Lhalf) xij0 -= L;
  if(xij > Lhalf) xij -= L;
  if(xij0 < -Lhalf) xij0 += L;
  if(xij < -Lhalf) xij += L;
#ifdef SECURE
  if(xij0<-Lhalf || xij0>Lhalf)
    Warning(" xij0= %lf out of range\n", xij0);
  if(xij<-Lhalf || xij>Lhalf)
    Warning(" xij= %lf out of range\n", xij);
#endif
#endif
  absxij0 = fabs(xij0);
  absxij = fabs(xij);

  if(xij*xij0<0) // there was an exchange between particles
    particles_exchange_sign = -1;
  else // there was no exchange of particles
    particles_exchange_sign = 1;

  //if(absxij>Apar || absxij0>Apar) particles_exchange_sign = 1;
#ifdef BC_1DPBC
  //if(absxij>Lhalf-5.*sqrt(dt) || absxij0>Lhalf-5.*sqrt(dt)) particles_exchange_sign = 1;
  if(absxij>0.75*Lhalf || absxij0>0.75*Lhalf) particles_exchange_sign = 1;
  //if(absxij>5.*sqrt(dt) || absxij0>5.*sqrt(dt)) particles_exchange_sign = 1.;
#endif

  if(a<1e-6 && a>-1.1e-3) { // unitary gas
    //if(particles_exchange_sign == -1) {
    //  safe = log(1. - exp(-m/dt*Apar*Apar));
   // }
   // else 
{
      x = exp(-0.5*m/dt*(1.+(DOUBLE)particles_exchange_sign)*fabs(xij*xij0));
      if(fabs(x) > 1e-10) {
        safe = log(1. - x);
      }
      else { // return log(1-x) expansion
        safe = x + x*x / 2. + x*x*x / 3. + x*x*x*x / 4. + x*x*x*x*x / 5.;
      }
    }
    return safe;
  }

  // erfcx
  yij = sqrt(0.25*m/dt) * (absxij + absxij0) - sqrt(dt/m)/a;
  yerfcx = Erfcx(yij);
  yexp = yerfcx*exp(-0.5*m/dt*(1.+particles_exchange_sign)*fabs(xij*xij0));
  x = yexp * sqrt(PI*dt/m)/a; // a<0

  if(fabs(x) > 1e-10) {
    safe = log(1. + x);
  }
  else { // return log(1+x) expansion
    safe = x - x*x / 2. + x*x*x / 3. - x*x*x*x / 4. + x*x*x*x*x / 5.;
  }

  return safe;
}

/************************** Propagation Laplacian ***************************/
double PropagatorLaplacian(double m, double dt, double dr) {
// contribution from the Laplacian + (R(i)-R(i+/-1))^2 / 2

#ifdef BC_1DPBC
  dr = fabs(dr);
  if(dr>Lhalf) dr = L - dr;
#ifdef SECURE
  if(dr<0 || dr>Lhalf) Warning(" dr= %lf out of range\n", dr);
#endif
#endif

  return -0.5 * m * dr*dr / dt;
}

/*************************** Path Integral Pseudopotential ****************************/
DOUBLE PathIntegralPseudopotentialS(void) {
  DOUBLE dummy1, dummy2, dummy3;

  return PathIntegralPseudopotentialSdetailed(&dummy1, &dummy2, &dummy3);
}

DOUBLE PathIntegralPseudopotentialSdetailed(DOUBLE *Skincheck, DOUBLE *Sdeltacheck, DOUBLE *Swfcheck) {
  int i, j, w, w1; // index of the walker with [w-1] and [w+1]
  DOUBLE S, Skin, Sdelta, Swf;
  //DOUBLE mA, mB, mAB; // i.e. 2 mu
  DOUBLE tau;
#ifdef EXTERNAL_POTENTIAL
  DOUBLE Epot = 0.; // potential energy contribution
#endif

  tau = beta;

#ifdef INFINITE_MASS_COMPONENT_DN // 2\mu  = m
  mAB = 2.;
#endif

  if(SmartMC == PIGS_PSEUDOPOTENTIAL_OBDM) { // initialize with ira position
    W[worm.w].z[worm.i] = worm.z_ira;
  }

  Skin = Sdelta = Swf = 0.;
  for(w=0; w<Nwalkers; w++) {
    w1 = w-1;
    if(w == 0) {
      if(MC == PIMC) { // PIMC, cyclic conditions
        w1 = Npop-1;
      }
      else { // PIGS, skip head
        w = 1;
        w1 = 0;
      }
    }

#ifdef ONE_COMPONENT_CODE
    for(i=0; i<N; i++) Skin += PropagatorLaplacian(mA, tau, W[w1].z[i] - W[w].z[i]);
    for(i=0; i<N; i++) for(j=i+1; j<N; j++) Sdelta += PropagatorPseudopotenial1D_ij(mA,  aA, tau, W[w].z[i],   W[w].z[j],   W[w1].z[i],   W[w1].z[j]);

#ifdef EXTERNAL_POTENTIAL
    for(i=0; i<N; i++) Epot += Vext(W[w].x[i], W[w].y[i], W[w].z[i]);
#endif
#else // two component code
#ifdef FINITE_MASS_COMPONENT_UP // mA finite
    for(i=0; i<Nup; i++) Skin += PropagatorLaplacian(mA, tau, W[w1].z[i] - W[w].z[i]);
    for(i=0; i<Nup; i++) for(j=i+1; j<Nup; j++) Sdelta += PropagatorPseudopotenial1D_ij(mA,  aA, tau, W[w].z[i],   W[w].z[j],   W[w1].z[i],   W[w1].z[j]);
#endif
#ifdef FINITE_MASS_COMPONENT_DN // mB finite
    for(i=0; i<Ndn; i++) Skin += PropagatorLaplacian(mB, tau, W[w1].zdn[i] - W[w].zdn[i]);
    for(i=0; i<Ndn; i++) for(j=i+1; j<Ndn; j++) Sdelta += PropagatorPseudopotenial1D_ij(mB,  aB, tau, W[w].zdn[i], W[w].zdn[j], W[w1].zdn[i], W[w1].zdn[j]);
#endif
    for(i=0; i<Nup; i++) for(j=0; j<Ndn; j++)   Sdelta += PropagatorPseudopotenial1D_ij(mAB, a,  tau, W[w].z[i],   W[w].zdn[j], W[w1].z[i],   W[w1].zdn[j]);
#ifdef EXTERNAL_POTENTIAL
    for(i=0; i<Nup; i++) Epot += VextUp(W[w].x[i], W[w].y[i], W[w].z[i]);
    for(i=0; i<Ndn; i++) Epot += VextDn(W[w].xdn[i], W[w].ydn[i], W[w].zdn[i]);
#endif
#endif
  }

  if(MC == PIGS) {
    Swf += U(W[0]);
    Swf += U(W[Nwalkers-1]);
  }

  S = Skin + Sdelta + Swf; // dSdelta - contribution from delta function, dSwf from the trial w.f.
  *Skincheck = Skin;
  *Sdeltacheck = Sdelta;
  *Swfcheck = Swf;

  if(SmartMC == PIGS_PSEUDOPOTENTIAL_OBDM) {/*
    i = worm.i;
    w = worm.w;
    w2 = w + 1; // next
    Skin -= PropagatorLaplacian(mA, tau, worm.z_ira   - W[w2].z[worm.i]); // this link does not belong to ira
    Skin += PropagatorLaplacian(mA, tau, worm.z_masha - W[w2].z[worm.i]); // this link belongs to masha

#ifdef ONE_COMPONENT_CODE
#ifdef EXTERNAL_POTENTIAL // potential energy is divided in two pieces
    Epot -= 0.5*Vext(W[w].x[i], W[w].y[i], worm.z_ira);
    Epot += 0.5*Vext(W[w].x[i], W[w].y[i], worm.z_ira);
#endif
    for(j=0; j<N; j++) if(j != worm.i) {
#else // two component code
#ifdef EXTERNAL_POTENTIAL // potential energy is divided in two pieces
    Epot -= 0.5*VextUp(W[w].x[i], W[w].y[i], worm.z_ira);
    Epot += 0.5*VextUp(W[w].x[i], W[w].y[i], worm.z_ira);
#endif
    for(j=0; j<Nup; j++) if(j != worm.i) {
#endif
      Sdelta -= PropagatorPseudopotenial1D_ij(mA, aA, tau, worm.z_ira,   W[w].z[j], W[w2].z[i], W[w2].z[j]); // remove ira + w2
      Sdelta += PropagatorPseudopotenial1D_ij(mA, aA, tau, worm.z_masha, W[w].z[j], W[w2].z[i], W[w2].z[j]); // add masha + w2
    }*/
  }

#ifdef EXTERNAL_POTENTIAL
  S -= tau*Epot; // potential energy contribution
#endif

  return S;
}

/*************************** Path Integral Pseudopotential ****************************/
void PathIntegralPseudopotentialStaging(void) {
  DOUBLE xs,ys,zs,ms; // r*, position of the weighted center and its mass
  DOUBLE dx,dy,dz;
  DOUBLE sigma;
  int l; // // number of links in the staging chain (minimal two)
  int i, j, w, w1, w2; // index of the walker with [w-1] and [w+1]
  DOUBLE dS;
  DOUBLE tau;
  int staging_begin, staging_end;
  DOUBLE staging_head_z;
#ifdef EXTERNAL_POTENTIAL
  DOUBLE Epot = 0.; // potential energy contribution
#endif
#ifdef SECURE
  DOUBLE Sold, Snew;
  DOUBLE dS_check, dSpot_check;
  DOUBLE dSkin_check_old, dSpot_check_old, dSwf_check_old;
  DOUBLE dSkin_check_new, dSpot_check_new, dSwf_check_new;
#endif

  tau = beta;
  xs = ys = zs = 0.;//???

#ifdef INFINITE_MASS_COMPONENT_DN // 2\mu  = m
  mAB = 2.;
#endif

  if(SmartMC == PIGS_PSEUDOPOTENTIAL_OBDM) { // initialize with ira position
    Error("Staging is not implemented for OBDM");
  }

  // throw indices of the first and last beads
  staging_begin = rand() % Nwalkers;
  staging_end = rand() % Nwalkers;
  // through which particle is moved
  i = rand() % N;

  if(staging_end>Npop-2) return;//!!!

  if(staging_begin > staging_end) { // swap w1 and w2
    w = staging_begin;
    staging_begin = staging_end;
    staging_end = w;
  }
  l = staging_end - staging_begin;
  if(l<2) return; // nothing to be done

#ifdef SECURE
  Sold = PathIntegralPseudopotentialSdetailed(&dSkin_check_old, &dSpot_check_old, &dSwf_check_old);
#endif

  //for(w=staging_begin; w<staging_begin; w++) { // store coordinates
  for(w=0; w<Npop; w++) { // store coordinates!!!
    CaseX(Wp[w].x[i] = W[w].x[i]);
    CaseY(Wp[w].y[i] = W[w].y[i]);
    CaseZ(Wp[w].z[i] = W[w].z[i]);
  }

  staging_head_z = W[staging_end].z[i];
  if(staging_head_z - W[staging_begin].z[i] > Lhalf) {
    staging_head_z -= L;
  }
  if(staging_head_z - W[staging_begin].z[i] < -Lhalf) {
    staging_head_z += L;
  }

  // generate free path corresponding to Laplacian
  // j0 = w1
  // j = w - w1
  // l = w2 - w1
  //for(j=0; j<=l-2; j++) { // generate new path according to the exact kinetic propagator
  //w = staging_begin + 1 + j;
  for(w=staging_begin+1; w<staging_end; w++) { // generate new path according to the exact kinetic propagator
    j = w - staging_begin - 1;
    ms = mA*(DOUBLE)(l-j)/(DOUBLE)(l-j-1);
    sigma = sqrt(tau/ms);

    CaseZ(zs = (staging_head_z + W[w-1].z[i] * (DOUBLE)(l - j - 1))/(DOUBLE)(l-j));
    RandomNormal3(&dx, &dy, &dz, 0., sigma); //(x,y,z, mu, sigma)
    CaseX(W[w].x[i] = xs + dx;)
    CaseY(W[w].y[i] = ys + dy;)
    CaseZ(W[w].z[i] = zs + dz;)
  }

  for(w=staging_begin+1; w<staging_end; w++) ReduceToTheBoxXYZ(&W[w].x[i], &W[w].y[i], &W[w].z[i]); // put the new part into the box

  for(w=staging_begin+1; w<staging_end; w++) { // calculate dS
    w1 = w-1;
#ifdef ONE_COMPONENT_CODE
    for(j=0; j<N; j++) if(j!= i) dS += PropagatorPseudopotenial1D_ij(mA, aA, tau, W[w].z[i],  W[w].z[j], W[w1].z[i],  W[w1].z[j]);
    for(j=0; j<N; j++) if(j!= i) dS -= PropagatorPseudopotenial1D_ij(mA, aA, tau, Wp[w].z[i], W[w].z[j], Wp[w1].z[i], W[w1].z[j]);
#ifdef EXTERNAL_POTENTIAL
    Epot += Vext(W[w].x[i], W[w].y[i], W[w].z[i]); //!!! to be completed
    Epot -= Vext(Wp[w].x[i], Wp[w].y[i], Wp[w].z[i]);
#endif
#else // two component code
#ifdef FINITE_MASS_COMPONENT_UP // mA finite
    for(j=0; j<N; j++) if(j!= i) Sdelta += PropagatorPseudopotenial1D_ij(mA,  aA, tau, W[w].z[i],   W[w].z[j],   W[w1].z[i],   W[w1].z[j]);
#endif
#ifdef FINITE_MASS_COMPONENT_DN // mB finite
    for(j=0; j<N; j++) if(j!= i) Sdelta += PropagatorPseudopotenial1D_ij(mB,  aB, tau, W[w].zdn[i], W[w].zdn[j], W[w1].zdn[i], W[w1].zdn[j]);
#endif
    for(j=0; j<N; j++) if(j!= i)   Sdelta += PropagatorPseudopotenial1D_ij(mAB, a,  tau, W[w].z[i],   W[w].zdn[j], W[w1].z[i],   W[w1].zdn[j]);
#ifdef EXTERNAL_POTENTIAL
    Epot += VextUp(W[w].x[i], W[w].y[i], W[w].z[i]);
    Epot += VextDn(W[w].xdn[i], W[w].ydn[i], W[w].zdn[i]);
#endif
#endif
  }

  dS = 0.;
  w = staging_end; // !!! should be added to for()
  w1 = w-1;
  for(j=0; j<N; j++) if(j!= i) dS += PropagatorPseudopotenial1D_ij(mA, aA, tau, W[w].z[i],  W[w].z[j], W[w1].z[i],  W[w1].z[j]);
  for(j=0; j<N; j++) if(j!= i) dS -= PropagatorPseudopotenial1D_ij(mA, aA, tau, Wp[w].z[i], W[w].z[j], Wp[w1].z[i], W[w1].z[j]);

  //Swf += U(W[0]);
  //Swf += U(W[Nwalkers-1]);
  //S = Skin + Sdelta + Swf; // dSdelta - contribution from delta function, dSwf from the trial w.f.

  /*if(SmartMC == PIGS_PSEUDOPOTENTIAL_OBDM) {
    i = worm.i;
    w = worm.w;
    w2 = w + 1; // next
    Skin -= PropagatorLaplacian(mA, tau, worm.z_ira   - W[w2].z[worm.i]); // this link does not belong to ira
    Skin += PropagatorLaplacian(mA, tau, worm.z_masha - W[w2].z[worm.i]); // this link belongs to masha

#ifdef ONE_COMPONENT_CODE
#ifdef EXTERNAL_POTENTIAL // potential energy is divided in two pieces
    Epot -= 0.5*Vext(W[w].x[i], W[w].y[i], worm.z_ira);
    Epot += 0.5*Vext(W[w].x[i], W[w].y[i], worm.z_ira);
#endif
    for(j=0; j<N; j++) if(j != worm.i) {
#else // two component code
#ifdef EXTERNAL_POTENTIAL // potential energy is divided in two pieces
    Epot -= 0.5*VextUp(W[w].x[i], W[w].y[i], worm.z_ira);
    Epot += 0.5*VextUp(W[w].x[i], W[w].y[i], worm.z_ira);
#endif
    for(j=0; j<Nup; j++) if(j != worm.i) {
#endif
      Sdelta -= PropagatorPseudopotenial1D_ij(mA, aA, tau, worm.z_ira,   W[w].z[j], W[w2].z[i], W[w2].z[j]); // remove ira + w2
      Sdelta += PropagatorPseudopotenial1D_ij(mA, aA, tau, worm.z_masha, W[w].z[j], W[w2].z[i], W[w2].z[j]); // add masha + w2
    }
  }*/

#ifdef EXTERNAL_POTENTIAL
  dS -= tau*Epot; // potential energy contribution
#endif

  //dS = dSkin + dSdelta + dSwf; // dSdelta - contribution from delta function, dSwf from the trial w.f.
#ifdef EXTERNAL_POTENTIAL
  //dS -= tau*dEpot; // potential energy contribution
#endif

#ifdef SECURE
  Snew = PathIntegralPseudopotentialSdetailed(&dSkin_check_new, &dSpot_check_new, &dSwf_check_new);
  dSpot_check = dSpot_check_new-dSpot_check_old;
  if(fabs((dSpot_check - dS)/dS)>1e-3 && fabs(dS)>1e-6) {
    Message("dS %lf %lf\n", dS, dSpot_check_new-dSpot_check_old);
  }
#endif

  if(Random() < Exp(dS)) { // The Metropolis code
    W[w].U -= dS;
    accepted_all++;
  }
  else {
    for(w=staging_begin; w<staging_end; w++) { // restore coordinates
      CaseX(W[w].x[i] = Wp[w].x[i]);
      CaseY(W[w].y[i] = Wp[w].y[i]);
      CaseZ(W[w].z[i] = Wp[w].z[i]);
    }
    rejected_all++;
  }
}
