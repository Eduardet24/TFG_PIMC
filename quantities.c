/*quantities.c*/

#include <memory.h>
#include <stdio.h>
#include "main.h"
#include "crystal.h"
#include "quantities.h"
#include "utils.h"
#include "trial.h"
#include "vmc.h"
#include "dmc.h"
#include "randnorm.h"
#include "rw.h"
#include "trial.h"
#include "spline.h"
#include "speckles.h"
#include "memory.h"
#include "compatab.h"
#include "ewald.h"
#include "mymath.h"
#include MATHINCLUDE

/********************************* Energy ************************************/
static DOUBLE EnergyPIMC(DOUBLE *Epot, DOUBLE *Ekin, DOUBLE *Eff,
                         DOUBLE *Edamping, DOUBLE *Eint, DOUBLE *Eext) {
  DOUBLE spring = 0.;
  DOUBLE virial = 0.;
  DOUBLE centroid_z = 0.;
  DOUBLE potential = 0.;
  DOUBLE external = 0.;
  DOUBLE interaction = 0.;
  DOUBLE tau = 1. / (T * (DOUBLE) Nwalkers);
  int bead;
  int next_bead;
  int particle;
  int other_particle;
  int terms = 0;
  DOUBLE dr[3];
  DOUBLE r2;

  for(bead=0; bead<Nwalkers; bead++) {
    if(W[bead].status != ALIVE) continue;
    next_bead = bead + 1;
    if(next_bead == Nwalkers) next_bead = 0;
    terms++;
    centroid_z += W[bead].z[0];

    for(particle=0; particle<N; particle++) {
#ifdef TRIAL_1D
      spring += (W[bead].z[particle] - W[next_bead].z[particle]) *
                (W[bead].z[particle] - W[next_bead].z[particle]);
#elif defined(TRIAL_2D)
      spring += (W[bead].x[particle] - W[next_bead].x[particle]) *
                (W[bead].x[particle] - W[next_bead].x[particle]);
      spring += (W[bead].y[particle] - W[next_bead].y[particle]) *
                (W[bead].y[particle] - W[next_bead].y[particle]);
#else
      spring += (W[bead].x[particle] - W[next_bead].x[particle]) *
                (W[bead].x[particle] - W[next_bead].x[particle]);
      spring += (W[bead].y[particle] - W[next_bead].y[particle]) *
                (W[bead].y[particle] - W[next_bead].y[particle]);
      spring += (W[bead].z[particle] - W[next_bead].z[particle]) *
                (W[bead].z[particle] - W[next_bead].z[particle]);
#endif

#ifdef EXTERNAL_POTENTIAL
      external += Vext(W[bead].x[particle], W[bead].y[particle], W[bead].z[particle]);
#endif
    }

    for(particle=0; particle<N; particle++) {
      for(other_particle=particle+1; other_particle<N; other_particle++) {
        dr[0] = W[bead].x[particle] - W[bead].x[other_particle];
        dr[1] = W[bead].y[particle] - W[bead].y[other_particle];
        dr[2] = W[bead].z[particle] - W[bead].z[other_particle];
        r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
        interaction += InteractionEnergy(Sqrt(r2));
      }
    }
  }

  if(terms == 0) {
    *Epot = *Ekin = *Eff = *Edamping = *Eint = *Eext = 0.;
    return 0.;
  }

  *Eext = external / ((DOUBLE) terms * N);
  *Eint = interaction / ((DOUBLE) terms * N);
  *Epot = *Eext + *Eint;
  centroid_z /= (DOUBLE) terms;
#ifdef TRIAL_1D
  if(N == 1) {
    for(bead=0; bead<Nwalkers; bead++) {
      if(W[bead].status == ALIVE) {
        virial += (W[bead].z[0] - centroid_z) * omega_z2 * W[bead].z[0];
      }
    }
    *Ekin = D * T / 2. + virial / (2. * (DOUBLE) terms * N);
  }
  else {
    *Ekin = D / (2. * tau) - spring / (2. * (DOUBLE) terms * tau * tau * N);
  }
#else
  *Ekin = D / (2. * tau) - spring / (2. * (DOUBLE) terms * tau * tau * N);
#endif
  *Eff = *Edamping = 0.;
  potential = *Epot + *Ekin;
  return potential;
}

// returns energy of the system averaged over all walkers
DOUBLE Energy(DOUBLE *Epot, DOUBLE *Ekin, DOUBLE *Eff, DOUBLE *Edamping, DOUBLE *Eint, DOUBLE *Eext) {
  DOUBLE dEpot, dEkin, dEff, dEdamping, dEint, dEext;
  int w;
  int terms = 0;

  if(MC == PIMC) {
    return EnergyPIMC(Epot, Ekin, Eff, Edamping, Eint, Eext);
  }

  *Epot = *Ekin = *Eff = *Edamping = *Eint = *Eext = 0.;
  terms = 0;
  for(w=0; w<Nwalkers; w++) {
    if(W[w].status == ALIVE) {
      WalkerEnergy(&W[w], &dEpot, &dEkin, &dEff, &dEdamping, &dEint, &dEext);
      *Epot += dEpot;
      *Ekin += dEkin;
      *Eff  += dEff;
      *Eint += dEint;
      *Eext += dEext;
      *Edamping += dEdamping;
      terms++;
    }
  }

  *Epot /= (DOUBLE) (N * terms);
  *Ekin /= (DOUBLE) (N * terms);
  *Eff  /= (DOUBLE) (N * terms);
  *Edamping /= (DOUBLE) (N * terms);
  *Eint /= (DOUBLE) (N * terms);
  *Eext /= (DOUBLE) (N * terms);

  return *Epot + *Ekin;
}

/************************ Walker Energy **************************************/
DOUBLE WalkerEnergy0(struct Walker *W) {
  DOUBLE E;

  E = WalkerEnergy(W, &W->Epot, &W->Ekin, &W->EFF, &W->Edamping, &W->Eint, &W->Eext);

  W->Epot /= (DOUBLE) N;
  W->Ekin /= (DOUBLE) N;
  W->EFF /= (DOUBLE) N;
  W->Edamping /= (DOUBLE) N;
  W->Eint /= (DOUBLE) N;
  W->Eext /= (DOUBLE) N;

  return E;
}

DOUBLE WalkerEnerg_sum_of_two_wavefunctions(struct Walker *walker, DOUBLE *Epot, DOUBLE *Ekin, DOUBLE *EFF, DOUBLE *Edamping, DOUBLE *Eint, DOUBLE *Eext) {
  DOUBLE U1, U2, r;
  DOUBLE E1,Epot1,Ekin1,EFF1,Edamping1,Eint1,Eext1;
  DOUBLE E2,Epot2,Ekin2,EFF2,Edamping2,Eint2,Eext2;
  DOUBLE weight1, weight2;
  DOUBLE E;

  walker->status_sum_of_two_walkers_U = 1;
  walker->status_sum_of_two_walkers_E = 1;
  b = sum_of_two_wavefunctions_parameter1; // parameters of wavefunction 1
  U1 = U(*walker);
  E1 = WalkerEnergy(W, &Epot1, &Ekin1, &EFF1, &Edamping1, &Eint1, &Eext1);
  b = sum_of_two_wavefunctions_parameter2; // parameters of wavefunction 2
  U2 = U(*walker);
  E2 = WalkerEnergy(W, &Epot2, &Ekin2, &EFF2, &Edamping2, &Eint2, &Eext2);
  walker->status_sum_of_two_walkers_U = 0;
  walker->status_sum_of_two_walkers_E = 0;

  r = exp(U1-U2); // psi1 / psi2
  weight1 = sum_of_two_wavefunctions_alpha * r / (sum_of_two_wavefunctions_alpha * r + 1. - sum_of_two_wavefunctions_alpha);
  weight2 = (1. - sum_of_two_wavefunctions_alpha) / (sum_of_two_wavefunctions_alpha * r + 1. - sum_of_two_wavefunctions_alpha);
  *Epot = weight1*Epot1 + weight2*Epot2;
  *Ekin = weight1*Ekin1 + weight2*Ekin2;
  *EFF = weight1*EFF1 + weight2*EFF2;
  *Eint = weight1*Eint1 + weight2*Eint2;
  *Eext = weight1*Eext1 + weight2*Eext2;
  E = weight1*E1 + weight2*E2;
  return E;
}

/* notation: F = (\nabla \Psi) / \Psi */
DOUBLE WalkerEnergy(struct Walker *walker, DOUBLE *Epot, DOUBLE *Ekin, DOUBLE *EFF, DOUBLE *Edamping, DOUBLE *Eint, DOUBLE *Eext) {
  int i,j,k;
  DOUBLE dr[3] = {0,0,0};
  DOUBLE r, r2;
  DOUBLE Force, dF;
#ifdef CRYSTAL
#ifdef CRYSTAL_SYMMETRIC
  DOUBLE M1, e;
#endif
#endif
#ifdef THREE_BODY_TERMS
  DOUBLE x2, y2, HR;
  DOUBLE CM12[3], xi[3], xj[3], xk[3], yk[3];
#endif
#ifdef SCALABLE
  int cc, ci, cj;
  if(MC == DIFFUSION) UpdateScalableTable(walker);
#endif
#ifdef MANY_BODY_TRIAL_TRAP_CENTER_OF_MASS_NODE_Z
  DOUBLE z_CM = 0.;
#endif
#ifdef MANY_BODY_TRIAL_TRAP_CENTER_OF_MASS_HERMITE2_Z
  DOUBLE z2 = 0.;
#endif

#ifdef DO_NOT_CALCULATE_ENERGY_AND_U
  return 0.;
#endif

#ifdef HARD_SPHERE_EXCLUDED_VOLUME_CHANGE_OF_VARIABLES
  ExcludedVolumeCoordinatesWrap(*walker); // convert variables x_i -> x_i - i * a
#endif

  if(MC == VARIATIONAL && SmartMC == VMC_SUM_OF_TWO_WAVEFUNCTIONS && walker->status_sum_of_two_walkers_E == 0) {
    return WalkerEnerg_sum_of_two_wavefunctions(walker, Epot, Ekin, EFF, Edamping, Eint, Eext);
  }

  *Epot = *Ekin = *EFF = *Edamping = *Eint = *Eext = 0.;

  ArrayEmpty2D(walker->F, i, N, k, 3, DOUBLE); // set Force array elements to zero

#ifdef CALCULATE_DERIVATIVES_NUMERICALLY
  CalculateDriftForceAndKineticEnergyNumerically(walker, Ekin);
#endif

#ifdef CRYSTAL
#ifdef CRYSTAL_SYMMETRIC // symmetric crystal
  for(i=0; i<Crystal.size; i++) { // do not change i to j as "j" is used in a macro
#ifdef TRIAL_1D
    *Ekin += 2.*Crystal_dot_Rz_i;
#endif
#ifdef TRIAL_2D
    *Ekin += 2.*(Crystal_dot_Rx_i + Crystal_dot_Ry_i);
#endif
#ifdef TRIAL_3D
    *Ekin += 2.*(Crystal_dot_Rx_i + Crystal_dot_Ry_i + Crystal_dot_Rz_i);
#endif
  }

  for(j=0; j<Crystal.size; j++) {
    walker->Mo[j] = 0.;     // fill Mo array, Mo[j] = w[j] sum_i e(i,j)
    for(i=0; i<N; i++) {
      dr[0] = walker->x[i] - Crystal.x[j];
      dr[1] = walker->y[i] - Crystal.y[j];
      dr[2] = walker->z[i] - Crystal.z[j];
      r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
      walker->Mo[j] += Exp(-Crystal_dot_Rx_j*dr[0]*dr[0]-Crystal_dot_Ry_j*dr[1]*dr[1]-Crystal_dot_Rz_j*dr[2]*dr[2]);
    }
    walker->Mo[j] *= Crystal_dot_weight_j;
  }

  for(i=0; i<N; i++) {
    for(j=0; j<Crystal.size; j++) {
      dr[0] = walker->x[i] - Crystal.x[j];
      dr[1] = walker->y[i] - Crystal.y[j];
      dr[2] = walker->z[i] - Crystal.z[j];
      r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
      e = Exp(-Crystal_dot_Rx_j*dr[0]*dr[0]-Crystal_dot_Ry_j*dr[1]*dr[1]-Crystal_dot_Rz_j*dr[2]*dr[2]);

      M1 = -2.*Crystal_dot_weight_j*e/walker->Mo[j];
      walker->F[i][0] += M1*Crystal_dot_Rx_j*dr[0];
      walker->F[i][1] += M1*Crystal_dot_Ry_j*dr[1];
      walker->F[i][2] += M1*Crystal_dot_Rz_j*dr[2];

      *Ekin -= 4.*( Crystal_dot_Rx_j*Crystal_dot_Rx_j*dr[0]*dr[0]
                   +Crystal_dot_Ry_j*Crystal_dot_Ry_j*dr[1]*dr[1]
                   +Crystal_dot_Rz_j*Crystal_dot_Rz_j*dr[2]*dr[2]
                  )*Crystal_dot_weight_j*e/walker->Mo[j]*(1. - Crystal_dot_weight_j*e/walker->Mo[j]);
    }
  }
#endif// end crystal symmetric
#endif // end CRYSTAL

#ifdef ONE_BODY_TRIAL_TRAP_CENTER_OF_MASS_NODE_Z
  for(i=0; i<N; i++) z_CM += walker->z[i];
  *Ekin += (DOUBLE) N * 2./ (z_CM*z_CM);
  for(i=0; i<N; i++) walker->F[i][2] += 1./z_CM;
#endif

  for(i=0; i<N; i++) { // one-body terms
#ifdef EXTERNAL_POTENTIAL // add the external potential energy [units of hw]
    *Eext += Vext(walker->x[i], walker->y[i], walker->z[i]);
#endif
#ifdef EXTERNAL_POTENTIAL_CONSTANT_SPINFULL
    if(walker->spin == 0) *Eext += D;
#endif

#ifdef ONE_BODY_TRIAL_TERMS
#ifdef CALCULATE_DERIVATIVES_ANALYTICALLY
    OneBodyFp(&walker->F[i][0], &walker->F[i][1], &walker->F[i][2], walker->x[i], walker->y[i], walker->z[i], i, walker->w);
    *Ekin += OneBodyE(walker->x[i], walker->y[i], walker->z[i], i, walker->w);
#endif
#endif
  }

#ifdef INTERACTION_SOFT_SPHERE_HYPERRADIUS
  *Eint += Dpar*(DOUBLE) CheckWalkerHyperradiusOverlappingCount(walker);
#endif

#ifdef SCALABLE_POTENTIAL
  for(cc=0; cc<Ncells; cc++) { // loop over cells
    for(ci=0; ci<walker->c_Nlocal[cc]; ci++) { // loop over particles in the local cell
      i = walker->c_index_local[cc][ci];
      for(cj=0; cj<walker->c_Nall[cc]; cj++) { // loop over particles in the nearest cells avoiding DOUBLE counting
        j = walker->c_index_all[cc][cj];
        if(j>i) {
#else
  for(i=0; i<N; i++) { // two-body terms: particle-particle force
    for(j=i+1; j<N; j++) {
#endif
      CaseX(dr[0] = walker->x[i] - walker->x[j]);
      CaseY(dr[1] = walker->y[i] - walker->y[j]);
      CaseZ(dr[2] = walker->z[i] - walker->z[j]);
      r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
      r = Sqrt(r2);

#ifdef THREE_BODY_TERMS
#ifdef CALCULATE_DERIVATIVES_ANALYTICALLY
      //  contains y[k] = r[i] - r[j]
      yk[0] = dr[0];
      yk[1] = dr[1];
      yk[2] = dr[2];
      y2 = r2; // y[k]

      // calculate CM12 = (r[i]+r[j])/2
      CM12[0] = walker->x[i] - 0.5*dr[0];
      CM12[1] = walker->y[i] - 0.5*dr[1];
      CM12[2] = walker->z[i] - 0.5*dr[2];

      for(k=j+1; k<N; k++) {
        // calculate x[k] = 2./sqrt(3)*(r[k] - (r[i]+r[j])/2)
        xk[0] = walker->x[k] - CM12[0];
        xk[1] = walker->y[k] - CM12[1];
        xk[2] = walker->z[k] - CM12[2];
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

        walker->F[i][0] += Force*xi[0];
        walker->F[i][1] += Force*xi[1];
        walker->F[i][2] += Force*xi[2];

        walker->F[j][0] += Force*xj[0];
        walker->F[j][1] += Force*xj[1];
        walker->F[j][2] += Force*xj[2];

        walker->F[k][0] += Force*xk[0];
        walker->F[k][1] += Force*xk[1];
        walker->F[k][2] += Force*xk[2];

        *Ekin += 2.*E3(HR);
      }
#endif
#endif

      // two-body potential energy
#ifdef INTERACTION_WITHOUT_IMAGES // might be restricted to r < L/2
      if(CheckInteractionConditionPotential(dr[0], dr[1], dr[2], r2)) *Eint += InteractionEnergy_ij(r, walker->spin[i], walker->spin[j]);
#else 
      *Eint += InteractionEnergyWithImages_ij(dr[0], dr[1], dr[2]);
#endif

#ifdef INTERACTION_WITH_DAMPING // Imaginary part of interaction potential is treated as a perturbation
      if(CheckInteractionConditionPotential(dr[0], dr[1], dr[2], r2)) *Edamping += InteractionDamping(r);
#endif

#ifdef CALCULATE_DERIVATIVES_ANALYTICALLY
      if(CheckInteractionConditionWF(dr[0], dr[1], dr[2], r2)) {
        *Ekin += 2.*InterpolateE(&G, r, walker->spin[i], walker->spin[j]);
        Force = InterpolateFp(&G, r, walker->spin[i], walker->spin[j])/r;
        CaseX(walker->F[i][0] += Force * dr[0]);
        CaseX(walker->F[j][0] -= Force * dr[0]);
        CaseY(walker->F[i][1] += Force * dr[1]);
        CaseY(walker->F[j][1] -= Force * dr[1]);
        CaseZ(walker->F[i][2] += Force * dr[2]);
        CaseZ(walker->F[j][2] -= Force * dr[2]);
#ifdef BC_SPHERE // add effect of the curvature
        Rnow = (walker->x[i]*walker->x[j]+walker->y[i]*walker->y[j]+walker->z[i]*walker->z[j]) /(Rsphere*Rsphere); // (ri . rj) / R^2
        walker->F[i][0] += Force * (-dr[0] + walker->x[i]*Rnow - walker->x[j]);
        walker->F[i][1] += Force * (-dr[1] + walker->y[i]*Rnow - walker->y[j]);
        walker->F[i][2] += Force * (-dr[2] + walker->z[i]*Rnow - walker->z[j]);
        walker->F[j][0] += Force * ( dr[0] + walker->x[j]*Rnow - walker->x[i]);
        walker->F[j][1] += Force * ( dr[1] + walker->y[j]*Rnow - walker->y[i]);
        walker->F[j][2] += Force * ( dr[2] + walker->z[j]*Rnow - walker->z[i]);
        *Ekin -= 2.*(1. + 0.5*(1. + Rnow))*InterpolateE(&G, r, walker->spin[i], walker->spin[j]); // InterpolateE = f"
        *Ekin -= 2.*Force*(Rnow - r*r/(4.*Rsphere*Rsphere));
        *Ekin += 2.*Force*Force*(1.-Rnow*Rnow)*Rsphere*Rsphere;
#endif
      }
#endif
    }
  }
#ifdef SCALABLE_POTENTIAL
  }}
#endif
  // calculate modulus squared
  for(i=0; i<N; i++) {
    CaseX(*EFF += walker->F[i][0]*walker->F[i][0]);
    CaseY(*EFF += walker->F[i][1]*walker->F[i][1]);
    CaseZ(*EFF += walker->F[i][2]*walker->F[i][2]);
  }

#ifdef SPINFULL_TUNNELING // caclulate contribution from tunneling
  U_psiT_sigma_inv(*walker);
  for(i=0; i<N; i++) *Eint -= t_tunneling*exp(walker->psiT_sigma_inv[i]);
#endif

#ifdef MANY_BODY_TRIAL_TRAP_CENTER_OF_MASS_NODE_Z
  for(i=0; i<N; i++) z_CM += walker->z[i];
  for(i=0; i<N; i++) *Ekin -= (2./z_CM)*walker->F[i][2]; // Delta Psi = Delta f_J + (2/f_MB) (nabla f_J)/(f_J)
#endif

#ifdef MANY_BODY_TRIAL_TRAP_CENTER_OF_MASS_HERMITE2_Z
  for(i=0; i<N; i++) z2 += (walker->z[i] / (DOUBLE) N)*(walker->z[i] / (DOUBLE) N);
  for(i=0; i<N; i++) *Ekin -= 2. * (2.*walker->z[i] / ((DOUBLE)(N*N) * (z2 - beta))) * (walker->F[i][2]); // 2 (nabla f_J) (nabla f_MB)
  *Ekin -= 2. / ((DOUBLE)N * (z2 - beta));
#endif

#ifdef HARD_SPHERE_EXCLUDED_VOLUME_CHANGE_OF_VARIABLES
  ExcludedVolumeCoordinatesUnwrap(*walker); // convert variables x_i -> x_i - i * a
#endif

// in a homogeneous 3D system energy is measured in units of [h^2/2ma^2]
#ifdef UNITS_SET_2M_TO_1
  *Eint *= 2.;
  *Eext *= 2.;
#else
// otherwise in units of [h^2/mr^2] or [hw], i.e. rescale the kinetic term by 1/2
  *Ekin *= 0.5;
  *EFF  *= 0.5;
#endif
  *Epot = *Eext + *Eint; // total potential energy
#ifdef CALCULATE_DERIVATIVES_ANALYTICALLY
  *Ekin -= *EFF;
#endif
  *EFF += *Eint + *Eext;

  return *Epot + *Ekin;
}

/************************ Walker Potential energy ****************************/
DOUBLE WalkerPotentialEnergy(struct Walker *walker) {
  int i,j;
  DOUBLE dr[3] = {0,0,0}, r, r2;
  DOUBLE V = 0.;

  // one-body terms
  for(i=0; i<N; i++) V += Vext(walker->x[i], walker->y[i], walker->z[i]);

  // empty Force arrays
  for(i=0; i<N; i++) { 
    walker->F[i][0] = walker->F[i][1] = walker->F[i][2] = 0.;
  }

  // two-body terms
  for(i=0; i<N; i++) { 
    for(j=i+1; j<N; j++) {
      CaseX(dr[0] = walker->x[i] - walker->x[j]);
      CaseY(dr[1] = walker->y[i] - walker->y[j]);
      CaseZ(dr[2] = walker->z[i] - walker->z[j]);

#ifdef SECURITY
#ifdef TRIAL_1D
      dr[0] = dr[1] = 0.;
#endif
#ifdef TRIAL_2D
      dr[2] = 0.;
#endif
#endif
      r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
      r = Sqrt(r2);

#ifdef INTERACTION_WITHOUT_IMAGES // might be restricted to r < L/2
      if(CheckInteractionConditionPotential(dr[0], dr[1], dr[2], r2)) V += InteractionEnergy_ij(r, walker->spin[i], walker->spin[j]);
#else 
      V += InteractionEnergyWithImages_ij(dr[0], dr[1], dr[2]);
#endif
    }
  }
  return V;
}

/************************ Walker Potential energy ****************************/
DOUBLE WalkerPotentialEnergyTilde(struct Walker *W) {
  int i,j,k;
  DOUBLE dr[3] = {0,0,0}, r, r2;
  DOUBLE V = 0., dF, Force;

  // one-body terms
  for(i=0; i<N; i++) V += VextTilde(W->x[i], W->y[i], W->z[i]);

  // two-body terms
  for(i=0; i<N; i++) { 
    for(j=i+1; j<N; j++) {
      dr[0] = W->x[i] - W->x[j];
      dr[1] = W->y[i] - W->y[j];
      dr[2] = W->z[i] - W->z[j];
      r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
      r = Sqrt(r2);

      Force = InteractionEnergyTilde(r) / r; // returns (V(r))'
      for(k=0; k<3; k++) {
        dF = Force * dr[k];
        W->F[i][k] += dF;
        W->F[j][k] -= dF;
      }
    }
  }

  // two-body terms
  for(i=0; i<N; i++) V += W->F[i][0]*W->F[i][0] + W->F[i][1]*W->F[i][1] + W->F[i][2]*W->F[i][2];

  return V;
}

/************************ External Energy ************************************/
DOUBLE Vext(DOUBLE x, DOUBLE y, DOUBLE z) {
  DOUBLE E = 0.;

#ifdef SPECKLES
  E += SpecklesEpot(x, y);
#endif

#ifdef VEXT_SPLINE_1D
  E += exp(InterpolateSplineU(&GV, z));
#endif

#ifdef TRAP_POTENTIAL
  CaseX(E += 0.5*omega_x2*x*x);
  CaseY(E += 0.5*omega_y2*y*y);
  CaseZ(E += 0.5*omega_z2*z*z);
#endif

#ifdef EXTERNAL_POTENTIAL_Uo_R2
  E -= Uo / (x*x+y*y+z*z);
#endif

#ifdef VEXT_COS2 // lattice
#ifdef BC_1DPBC_Z
  E += Srecoil*Cos(kL*z)*Cos(kL*z);
//  E += Srecoil*(Sin(kL*x)*Sin(kL*x) + Syrecoil*Syrecoil*Sin(kL*y)*Sin(kL*y)); // b = omega_y / omega_z aspect ratio, 1D CIR calculation
#endif
#ifdef BC_2DPBC_SQUARE
  E += Srecoil*(Cos(kL*x)*Cos(kL*x) + Cos(kL*y)*Cos(kL*y));
#endif
#ifdef BC_3DPBC_CUBE // set UNITS_SET_2M_TO_1 to get output in recoil energies
  E += Srecoil*(Cos(kL*x)*Cos(kL*x) + Cos(kL*y)*Cos(kL*y) + Cos(kL*z)*Cos(kL*z));
#endif
#ifdef BC_ABSENT
  E += Srecoil*(Sin(kL*x)*Sin(kL*x) + Sin(kL*y)*Sin(kL*y));
#endif
#endif

#ifdef VEXT_SOLITON_FIXED_PHASE
 E += 2.*soliton_fixed_phase_V2*(1.-soliton_fixed_phase_V2)*(1.-soliton_fixed_phase_V2)
 /(-1. + 2.*soliton_fixed_phase_V2 + cosh(sqrt(2. - 2.*soliton_fixed_phase_V2)*(z-Lhalf)/soliton_fixed_phase_Xi))
 /(-1. + 2.*soliton_fixed_phase_V2 + cosh(sqrt(2. - 2.*soliton_fixed_phase_V2)*(z-Lhalf)/soliton_fixed_phase_Xi));
#endif

#ifdef VEXT_IMPURITY_3DSW
  if((x-Lhalf)*(x-Lhalf)+(y-Lhalf)*(y-Lhalf)+(z-Lhalf)*(z-Lhalf)<RoSW*RoSW)
    E -= trial_Vo;
#endif

#ifdef VEXT_GAUSSIAN
  E += gaussian_alpha * exp( -gaussian_beta *z*z);
#endif

#ifdef VEXT_TRAP_DOUBLE_WELL
   E += (0.5*z*z + alpha_latt*Exp(-0.5*z*z) - (Log(alpha_latt)+1.));
#endif

  return E;
}

/************************ External Energy Tilde ******************************/
// [V, [T, V]]
DOUBLE VextTilde(DOUBLE x, DOUBLE y, DOUBLE z) {

#ifdef TRAP_POTENTIAL // i.e. the anisotropic oscillator
#ifdef TRIAL_1D
  // for 0.5*z*z the commutator gives z^2
  return z*z;
#else // 2D 3D
  Error("  VextTilde not implemented\n");
  return 0.;
#endif
#else  // no oscillator
  return 0.;
#endif
}

#define EULERGAMMA 1.781072417990198
/************************ Interaction Energy *********************************/
// spinfull interaction energy
DOUBLE InteractionEnergy_ij(DOUBLE r, int spin1, int spin2) {

#ifdef INTERACTION_SPINFULL_DIPOLE_BILAYER
  if(spin1 == spin2)
    //return 0.;
    return D/(r*r*r);
  else
    return D*(r*r-2.*bilayer_width2)*pow(r*r+bilayer_width2, -2.5);
#endif

#ifdef INTERACTION_SPINFULL_COULOMB_BILAYER
  if(spin1 == spin2)
    return 1./r;
  else
    return -1./sqrt(r*r+bilayer_width2);
#endif

#ifdef INTERACTION_SPINFULL_GAUSSIAN
  if(spin1 == spin2)
    return 0.;
  else
    return -D*exp(-0.5*r*r);
#endif

#ifdef INTERACTION_SPINFULL_GAUSSIANS_3COMPONENT
  if(spin1 == 0 && spin2 == 0) return cR[1]*exp(-r*r); //AA
  if(spin1 == 0 && spin2 == 1) return cR[2]*exp(-r*r); //AB
  if(spin1 == 1 && spin2 == 0) return cR[2]*exp(-r*r); //BA
  if(spin1 == 0 && spin2 == 2) return cR[3]*exp(-r*r); //AC
  if(spin1 == 2 && spin2 == 0) return cR[3]*exp(-r*r); //CA
  if(spin1 == 1 && spin2 == 2) return cR[4]*exp(-r*r); //BC
  if(spin1 == 2 && spin2 == 1) return cR[4]*exp(-r*r); //CB
#endif

#ifdef INTERACTION_SPINFULL_DOUBLE_GAUSSIAN
  if(spin1 == spin2)
    return 0.;
  else
    return -gaussian_alpha*exp(-0.5*r*r) - gaussian_beta*exp(-0.25*r*r);
#endif

#ifdef INTERACTION_SPINFULL_DIPOLE_MULTILAYER
  if(spin1 == spin2)
    return 1./(r*r*r);
  else
    return (r*r-2.*bilayer_width2*((DOUBLE)(spin1-spin2)*(spin1-spin2)))*pow(r*r+bilayer_width2*((DOUBLE)(spin1-spin2)*(spin1-spin2)), -2.5);
#endif

  // default choice
  return InteractionEnergy(r);
}

// spinless interaction energy
DOUBLE InteractionEnergy(DOUBLE r) {
#ifdef INTERACTION_ABSENT
  return 0.;
#endif

#ifdef INTERACTION_SPINFULL_DIPOLE_BILAYER // only 12 term, is used for spline construction
  return D*(r*r-2.*bilayer_width2)*pow(r*r+bilayer_width2, -2.5);
#endif

#ifdef INTERACTION_SPINFULL_COULOMB_BILAYER // only 12 term, is used for spline construction
  return -1./(r*r+bilayer_width2);
#endif

#ifdef INTERACTION_LENNARD_JONES
  DOUBLE rt;

  rt = sigma_lj/r;
  rt *= rt;
  rt *= (rt*rt);

  return U_lj*(rt*rt - rt);
#endif

#ifdef INTERACTION_LENNARD_JONES10 // powers 10, 6
  DOUBLE rt;

  rt = 1./r;// -1
  rt *= rt; // -2

  return D*(rt*rt*rt*(rt*rt - 1.));
#endif

#ifdef INTERACTION_SUTHERLAND
  static int first_time = ON;
  static DOUBLE sqE,sqE2;

  if(first_time) {
    sqE = PI/L;
    sqE2 = sqE*sqE;
    first_time = OFF;
  }

  return D*(D-1.)*sqE2/(Sin(sqE*r)*Sin(sqE*r));
#endif

#ifdef INTERACTION_CALOGERO
#ifdef TRIAL_1D
  return D*(D-1.)/(r*r);
#endif
#ifdef TRIAL_2D
  return D/(r*r);
#endif
#ifdef TRIAL_3D
  return D*(D+1.)/(r*r);
#endif
#endif

#ifdef INTERACTION_YUKAWA
  //return 2./(D*D)*K0(r);; // D is de Boer parameter, see Magro et al. PRB 48, 411
  //return 2/(D*D)*Exp(-2.*r)/r; // D is de Boer parameter
  if(r == 0.) return 0.;
  return 2*D*Exp(-2.*r)/r; // D = M/m is mass ratio
  //return 13.33333333333*Exp(-2.*r)/r;
  /*if(r>0.5)
    return 13.33333333333*Exp(-2.*r)/r*(1.-0.5/r);
  else
    return 0.;*/
#endif

#ifdef INTERACTION_K0
  return K0(r)/(2*D*D); // D is de Boer parameter
#endif

#ifdef INTERACTION_PETROV
  DOUBLE ko,k1,arg;
  if(r<1) {
    r = 1.;
    arg = 2./EULERGAMMA;
    ko = K0(arg);
    k1 = K1(arg);

    return 4./r*(1./r-1-arg*ko*D/EULERGAMMA*(ko-arg*k1));
  }
  else {
    arg = 2.*r/EULERGAMMA;
    ko = K0(arg);
    k1 = K1(arg);

    return 8*D*ko/(EULERGAMMA*EULERGAMMA)*(arg*k1-ko);
  }
#endif

#ifdef INTERACTION_PETROV3D
  if(r<1.) r = 1.;
  return D*2./r*Exp(-2.*r)*(1. - 0.5/r);
#endif

#ifdef INTERACTION_COULOMB
  return D/r;
#endif

#ifdef INTERACTION_COULOMB_1D
    return - (digamma(r/L) + digamma(1.-r/L)) / L; //+ log((DOUBLE)N) / L
#endif

#ifdef INTERACTION_LINEAR_1D_COULOMB
    return r/a;
#endif

#ifdef INTERACTION_WITH_DAMPING
  /*int i;
  DOUBLE dx;

  if(r<=G.min) return G.ReV[0];

  i = (int) ((r-G.min)/G.step);
  if(i>= G.size) return 0;

  dx = r - G.x[i];
  return G.ReV[i]+(G.ReV[i+1] - G.ReV[i])*G.I_step * dx;*/
#endif

#ifdef INTERACTION_ROTON
  if(r<1)
    return 0.;
  else if(r<2.5)
    return -1e-3;
  else
    return 1./(r*r*r);
#endif

#ifdef INTERACTION_SOFT_SPHERE
#ifdef TRIAL_3D
  if(r<Apar)
    return Vo;
  else
    return 0;
#else
  if(r<Ro)
    return Vo;
  else
    return 0;
#endif
#endif

#ifdef INTERACTION_SQUARE_WELL
  if(r<RoSW)
    return -Vo;
  else
    return 0;
#endif

#ifdef INTERACTION_GAUSSIAN // 2/sqrt(pi*alpha) exp(-x^2/alpha)
  //return 2./sqrt(PI*D)*exp(-r*r/D);
  return a*exp(-b*r*r);
#endif

#ifdef INTERACTION_KURBAKOV
  //return D*exp(-0.5*r*r);
  return (0.1 + a/r)/(1. + r*r*r*r*r*r/(a*a*a*a*a*a));
#endif

#ifdef INTERACTION_KURBAKOV_GENERAL
  DOUBLE y,a1,A;
  y=0;A=0;//                                                      no potential
  if(Ipar==1){y=D*exp(-Apar*log(r));if(Apar>=6)A=1E-5;}//         power-law
// Kurbakov 8-LJD
  if(Ipar==2){a1=sqrt(r*r+Apar*Apar);y=2*D/(r*a1*(r+a1));}//      separated dipoles
  if(Ipar==3)y=D*exp(-r/Apar)/r;//                                Yukawa
// Kurbakov 8-LJD
  if(Ipar==4){a1=Apar/r;a1*=a1;a1*=a1*a1;y=D*a1*(a1-1);A=1E-5;}// Lennard-Jones
  if(Ipar==5){//                                                  Lennard-Jones dipoles
// ! Kurbakov 8-LJD !
    a1=1/(r*r*r);y=1*a1;a1*=Apar*Apar*Apar;
    a1*=a1;y+=D*a1*(a1-1);A=1E-5;
  }if(Ipar==6)y=D*exp(-r*r/(2*Apar*Apar));//                      Gaussian-shaped
// Kurbakov B-POT7
  if(Ipar==7)y=D*(Apar+1/r)/(1+r*r*r*r*r*r);//                    Rydberg-Coulomb
// Kurbakov C-HE
  if(Ipar==10){
    DOUBLE Heepsilon = 10.948; // in K
    DOUBLE HeC6 = 1.36745214; //K/A^6
    DOUBLE Herm = 2.963; //A
    DOUBLE HeC8 = 0.42123807; //K/A^8
    DOUBLE HeD = 1.4826; //
    DOUBLE HeC10 = 0.17473318; //K/A^10
    DOUBLE Healpha = 10.43329537; //A^-1
    DOUBLE HeA = 184431.01; //K
    DOUBLE Hebeta = -2.27965105; //A^-2
    DOUBLE HeF;
    r/=Herm;if(r<HeD)HeF=exp(-(HeD/r-1)*(HeD/r-1));else HeF=1;
    if(HeF)HeF*=(HeC6+(HeC8+HeC10/(r*r))/(r*r))/(r*r*r*r*r*r);
    y=Heepsilon*(HeA*exp(-Healpha*r+Hebeta*r*r)-HeF)/12.115916747700853;r*=Herm;//A=1E-5;
// Kurbakov 9-CORE
  }if(Kurbakov_softcore>0){
#ifdef KURBAKOV_SOFT_CORE
    a1=log(100/r)-gammaE;a1=exp(2*a1+100*log(Kurbakov_softcore/r));
    if(a1<1E308)y+=a1;else y=a1;
    a1=Kurbakov_softcore;if(DIMENSION==2)a1*=a1;if(DIMENSION==3)a1*=a1*a1;
    if(a1<1)A=a1*1E-5;else A=1E-5;
#endif
// Kurbakov 2-CUT
  }if(A){//                                                        cutting off
    if(y>1E308)y=1E308;A*=y;if(A>2E2)y/=exp(7/8.*log(A));
    else if(A>7E-3)y/=sqrt(sqrt(sqrt(1+A*A*A*A*A*A*A)));
  }return(y);
#endif

#ifdef INTERACTION_POWER
  return D*exp(-Apar*log(r));
#endif

#ifdef INTERACTION_DOUBLE_GAUSSIAN
  return -gaussian_alpha*exp(-0.5*r*r) - gaussian_beta*exp(-0.25*r*r);
#endif

#ifdef INTERACTION_R2
  return 1./(r*r);
#endif

#ifdef INTERACTION_DIPOLE
  return D/(r*r*r);
#endif

#ifdef INTERACTION_DIPOLE_TRUNCATED
  if(r<a)
    return 1./(a*a*a);
  else
    return 1./(r*r*r);
#endif

#ifdef INTERACTION_DIPOLE_SCREENED
    return 1./(r*r*r+a*a*a);
#endif

#ifdef INTERACTION_R4
  return 1./(r*r*r*r);
#endif

#ifdef INTERACTION_QUADRUPOLE
  return 1./(r*r*r*r*r);
#endif

#ifdef INTERACTION_RYDBERG
  return 1./(r*r*r*r*r*r);
#endif

#ifdef INTERACTION_RYDBERG_TRUNCATED
  if(r<a)
    return 1./(a*a*a*a*a*a);
  else
    return 1./(r*r*r*r*r*r);
#endif

#ifdef INTERACTION_RYDBERG_SCREENED
  return a/(r*r*r*r*r*r+b*b*b*b*b*b);
#endif

#ifdef INTERACTION_RYDBERG_WITH_DELTA
  //return -0.5+sqrt(C6/(r*r*r*r*r*r)+0.25);
  //return -0.5*delta_tilde+sqrt(1./(r*r*r*r*r*r)+0.25*delta_tilde*delta_tilde);
  return delta_tilde*(sqrt(0.25+1./(r*r*r*r*r*r*delta_tilde))-0.5);
#endif

#ifdef INTERACTION_R12
  return 1./(r*r*r*r*r*r*r*r*r*r*r*r);
#endif

#ifdef INTERACTION_Aziz
  DOUBLE epsilon = 10.948; // in K
  DOUBLE C6 = 1.36745214; //K/A^6
  DOUBLE rm = 2.963; //A
  DOUBLE C8 = 0.42123807; //K/A^8
  DOUBLE D = 1.4826; //
  DOUBLE C10 = 0.17473318; //K/A^10
  DOUBLE alpha = 10.43329537; //A^-1
  DOUBLE A = 1.8443101e5; //K
  DOUBLE beta = -2.27965105; //A^-2
  DOUBLE x;
  DOUBLE F;

  r /= rm;
  if(r<D)
    F = exp(-(D/r-1.)*(D/r-1.));
  else
    F = 1.;

  F *= (C6 + (C8 + C10/(r*r))/(r*r))/(r*r*r*r*r*r);

  return epsilon*( A*exp(-alpha*r + beta*r*r) -F)/energy_unit;
#endif

#ifdef INTERACTION_LOGARITHMIC
  return log(r);
#endif
  
#ifdef INTERACTION_HYDROGEN
  DOUBLE x;

  //return epsilon*( A*exp(-alpha*r + beta*r*r) -F)/energy_unit;
  if(r<10.5824278175887) {
    return InterpolateSplineF(&GV, r) / energy_unit;
  }
  else {
    x = r/10.5824278175887;
    return -3.340686042602373E-002/(x*x*x*x*x*x) / energy_unit;
  }
#endif

#ifdef INTERACTION_Q1D_DIPOLES
  return D*(-2.*r + sqrt(2*PI)*(1. + r*r)*Erfcx(r/sqrt(2.)));
#endif

#ifdef INTERACTION_CH2
  DOUBLE e;

  e = 1./cosh(r);
  return -Epar*e*e;
#endif

#ifdef INTERACTION_HARD_SPHERE_SQUARE_WELL
  if(r<Ro) return -trial_Vo;
#endif

#ifdef INTERACTION_GAPPED
  return D*D*(0.5-r/L)*(0.5-r/L);
#endif

  return 0.;
}

#ifdef INTERACTION_WITH_DAMPING
DOUBLE InteractionDamping(DOUBLE r) {
  int i;
  DOUBLE dx;

  if(r<=G.min) return G.ImV[0];

  i = (int) ((r-G.min)/G.step);
  if(i>= G.size) return 0;

  dx = r - G.x[i];
  return G.ImV[i]+(G.ImV[i+1] - G.ImV[i])*G.I_step * dx;
}
#endif

/************************ Interaction Energy *********************************/
// spinfull interaction energy
DOUBLE InteractionEnergyWithImages_ij(DOUBLE x, DOUBLE y, DOUBLE z) {
#ifdef INTERACTION_SUM_OVER_IMAGES 
  return InteractionEnergySumOverImages(x, y, z);
#endif
#ifdef INTERACTION_EWALD
  return EwaldSumPair(x, y, z);
#endif
#ifdef INTERACTION_SUM_OVER_IMAGES_SMOOTH_CUTOFF
  return SmoothCutoffSum(x, y, z);
#endif
  return 0;
}

/************************ Interaction Energy Tilde ***************************/
// [V, [T, V]]
DOUBLE InteractionEnergyTilde(DOUBLE r) {

#ifdef INTERACTION_ABSENT
  return 0.;
#endif

#ifdef INTERACTION_SUTHERLAND
  static int first_time = ON;
  static DOUBLE coef;
  DOUBLE V;

  if(first_time) {
    sqE = PI/L;
    sqE2 = sqE*sqE;

    coef = 2.*D*(D-1.)*sqE2*sqE;

    first_time = OFF;
  }

  V = coef/(tan(sqE*r)*sin(sqE*r)*sin(sqE*r));
  //return V*V;

  return V;
#endif

#ifdef INTERACTION_DIPOLE
  //return 9.*D*D/(r*r*r*r*r*r*r*r);
  return 3.*D/(r*r*r*r);
#endif

  Error("  Interaction energy not implemented!");
  return 0.;

}

/****************************** pair distribution *****************************/
void MeasurePairDistribution(void) {
  int w,i,j,k,n;
  DOUBLE r, r2, dr[3] = {0,0,0};
  DOUBLE R2;
  int PD_position;
  static int wait = 1;
#ifdef SPINFULL
  int m;
#endif
  PD.times_measured += Nwalkers;

#ifdef PAIR_DISTRBUTION_EMPTINESS_FORMATION_PROBABILITY_1D
  int empty;

  for(n=0; n<PD.size; n++) { // do the measurement for size PD.
    for(w=0; w<Nwalkers; w++) {
      empty = 1;
      for(i=0; i<N; i++) {
        if(W[w].z[i] < PD.step*(int) n) { 
         empty = 0;
         break;
        }
      }
      PD.N[n] += empty;
      PD.times_measured++;
    }
  }
  return;
#endif

  if(MC == DIFFUSION && wait == OFF) { // do a pure measurement
    //PD_pure.times_measured += Nwalkers;
    for(w=0; w<Nwalkers; w++) {
      PD_position = W[w].PD_position;
      //PD_pure.times_measured += W[w].weight;
      PD_pure.times_measured += 1;
      for(n=0; n<gridPD; n++) {
        //PD_pure.N[n] += W[w].PD[PD_position][n] * W[w].weight;
        PD_pure.N[n] += W[w].PD[PD_position][n];
        W[w].PD[PD_position][n] = 0;
#ifdef SPINFULL
        for(m=0; m<Nspin; m++) {
          PD_pure.PDSpin[m][n] += W[w].PDSpin[PD_position][m][n];
          W[w].PDSpin[PD_position][m][n] = 0;
        }
#endif
      }
      if(measure_PairDistr_R2) {
        PD_pure.r2 += W[w].PD_R2pure[PD_position];
      }
    }
  }

  for(w=0; w<Nwalkers; w++) {
    if(measure_PairDistr_R2) R2 = 0.;
    for(i=0; i<N; i++) {
      for(j=i+1; j<N; j++) {
        CaseX(dr[0] = W[w].x[i] - W[w].x[j]);
        CaseY(dr[1] = W[w].y[i] - W[w].y[j]);
        CaseZ(dr[2] = W[w].z[i] - W[w].z[j]);

#ifdef LATTICE_ZIGZAG
        dr[1] = 0.;
#endif
        r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
        if(measure_PairDistr_R2) R2 += r2;
        r = Sqrt(r2);
#ifdef HARD_SPHERE
        if(r2 < a2) Warning("Pair Distribution: Particles %i and %i overlap\n", i+1, j+1);
#endif
       if(CheckInteractionConditionPotential(dr[0], dr[1], dr[2], r2)) {
          n = (int) (r / PD.step);
          if(n<PD.size) {
            PD.N[n]++; // mixed estimator
            if(MC == DIFFUSION) W[w].PD[W[w].PD_position][n]++; // pure estimator

#ifdef SPINFULL
            m = (int)abs(W[w].spin[i] - W[w].spin[j]); // difference in spin
            PD.PDSpin[m][n]++;
            if(MC == DIFFUSION) W[w].PDSpin[W[w].PD_position][m][n]++; // pure estimator
#endif 
          }
        }

        if(measure_PairDistrMATRIX) {
#ifndef TRIAL_1D
          n = (int) ((dr[0]+PD_MATRIX_x) / (2.*PD_MATRIX.step_x));
          k = (int) ((dr[1]+PD_MATRIX_y) / (2.*PD_MATRIX.step_y));
          if(n>=0 && k>=0 && n<gridPD_MATRIX_x && k<gridPD_MATRIX_y) {
            PD_MATRIX.N[n][k]++;
            PD_MATRIX.times_measured++;
          }
          n = (int) ((-dr[0]+PD_MATRIX_x) / (2.*PD_MATRIX.step_x));
          k = (int) ((-dr[1]+PD_MATRIX_y) / (2.*PD_MATRIX.step_y));
          if(n>=0 && k>=0 && n<gridPD_MATRIX_x && k<gridPD_MATRIX_y) {
            PD_MATRIX.N[n][k]++;
            PD_MATRIX.times_measured++;
          }
#else // trapped 1D in z (n(r_1) n(r_2))
          n = (int) ((W[w].z[i]+PD_MATRIX_x) / (2.*PD_MATRIX.step_x));
          k = (int) ((W[w].z[j]+PD_MATRIX_x) / (2.*PD_MATRIX.step_y));
          if(n>=0 && k>=0 && n<gridPD_MATRIX_x && k<gridPD_MATRIX_y) {
            PD_MATRIX.N[n][k]++;
            PD_MATRIX.times_measured++;
            PD_MATRIX.N[k][n]++;
            PD_MATRIX.times_measured++;
          }
#endif
        }

        // measure g3(0): all distances (i,j), (i, k), (j,k) smaller than dx
        if(r < PD.step) { 
          for(k=j+1; k<N; k++) { // pair (i,k)
            CaseX(dr[0] = W[w].x[i] - W[w].x[k]);
            CaseY(dr[1] = W[w].y[i] - W[w].y[k]);
            CaseZ(dr[2] = W[w].z[i] - W[w].z[k]);
            if(FindNearestImage(&dr[0], &dr[1], &dr[2]) < PD.step*PD.step) {
              PD.g3++;
            }
          }
        }
      }
    }
    if(measure_PairDistr_R2) {
      R2 /= (DOUBLE) (N*(N-1)/2);
      PD.r2 += R2;
      if(MC == DIFFUSION) W[w].PD_R2pure[W[w].PD_position] = R2;
    }

    if(MC == DIFFUSION) { // move current position
      W[w].PD_position++;
      if(W[w].PD_position == grid_pure_block) W[w].PD_position = 0;
    }
  }

  if(MC == DIFFUSION && wait) { // skip the first block in PURE case
    wait++;
    if(wait == grid_pure_block +1) wait = OFF;
  }
}

/****************************** pair distribution *****************************/
void MeasureHyperradiusDistribution(void) {
  int w,i,j,k,n;
  DOUBLE HR, r2, dr[3] = {0,0,0};
  DOUBLE CM12[3] = {0,0,0};
  DOUBLE rho2;
  int HD_position;
  static int wait = 1;

  //HD.times_measured += Nwalkers;

  /*if(MC == DIFFUSION && wait == OFF) { // do a pure measurement
    HD_pure.times_measured += Nwalkers;
    for(w=0; w<Nwalkers; w++) {
      HD_position = W[w].HD_position;
      //PD_pure.times_measured += W[w].weight;
      PD_pure.times_measured += 1;
      for(n=0; n<gridHD; n++) {
        HD_pure.N[n] += W[w].HD[HD_position][n];
        W[w].HD[HD_position][n] = 0;
      }
    }
  }

  for(w=0; w<Nwalkers; w++) {
    for(i=0; i<N; i++) {
      for(j=i+1; j<N; j++) {
        // distance (1-2)
        dr[0] = W[w].x[i] - W[w].x[j];
        dr[1] = W[w].y[i] - W[w].y[j];
        dr[2] = W[w].z[i] - W[w].z[j];
        rho2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);

        // center of mass (1-2) position
        CM12[0] = W[w].x[i] - 0.5*dr[0];
        CM12[1] = W[w].y[i] - 0.5*dr[1];
        CM12[2] = W[w].z[i] - 0.5*dr[2];

        for(k=j+1; k<N; k++) {
          // distance (CM12 - 3)
          dr[0] = W[w].x[k] - CM12[0];
          dr[1] = W[w].y[k] - CM12[1];
          dr[2] = W[w].z[k] - CM12[2];
          r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);

          HR = sqrt(4./3.*r2 + rho2);

          n = (int) (HR / HD.step);
          if(n<HD.size) {
            HD.N[n]++; // mixed estimator
            //if(MC == DIFFUSION) W[w].HD[W[w].HD_position][n]++; // pure estimator
          }
        }
      }
    }
    //if(MC == DIFFUSION) { // move current position
    //  W[w].HD_position++;
    //  if(W[w].HD_position == grid_pure_block) W[w].HD_position = 0;
    // }
  }

  ///f(MC == DIFFUSION && wait) { // skip the first block in PURE case
  //  wait++;
  //  if(wait == grid_pure_block +1) wait = OFF;
  //}*/
}

/**************************** radial distribution *****************************/
void MeasureRadialDistribution(void) {
  int w;
  static int i=0;
  static int first_time = ON;

  //if(++i != Nmeasure) return;

  if(first_time == OFF && MC == DIFFUSION) SaveMeanR2DMC();
  i=0;

  RD.r2recent = RD.r4 = 0.;
  if(measure_RDz) RDz.r2recent = RDz.r4 = 0.;

  for(w=0; w<Nwalkers; w++) {
    MeasureRadialDistributionWalker(&W[w]);
  }

  RD.r2recent = RD.r2recent / (DOUBLE) (Nwalkers * N);
  //RD.r4 = RD.r4 / (DOUBLE) (Nwalkers * N) - RD.r2recent*RD.r2recent;
  //RD.r2recent = Sqrt(RD.r2recent);
  //RD.r4 = Sqrt(RD.r4);
  RD.r4 = RD.r4 / (DOUBLE) (Nwalkers * N);

  if(measure_RDz) {
    RDz.r2recent = RDz.r2recent / (DOUBLE) (Nwalkers * N);
    //RDz.r4 = RDz.r4 / (DOUBLE) (Nwalkers * N) - RDz.r2recent*RDz.r2recent;
    RDz.r4 = RDz.r4 / (DOUBLE) (Nwalkers * N);
    //RDz.r2recent = Sqrt(RDz.r2recent);
    //RDz.r4 = Sqrt(RDz.r4);
  }
}

/**************************** radial distribution walker **********************/
void MeasureRadialDistributionWalker(struct Walker *W) {
  int i, n, nx, ny;
  DOUBLE r2, r, z2, R2, Z2, R4, Z4;
  DOUBLE xCM, yCM, zCM;
  int calculate_CM = OFF;

  R2 = Z2 = 0.;
  R4 = Z4 = 0.;
  xCM = yCM = zCM = 0.;
#pragma omp atomic
  RD.times_measured++;

  if(measure_RDz) {
#pragma omp atomic
    RDz.times_measured++;
  }

  if(MC == DIFFUSION && W->RD_wait == OFF) { // do a pure measurement
#pragma omp atomic
    RD_pure.times_measured++;
    for(n=0; n<gridRD; n++) {
#pragma omp atomic
      RD_pure.N[n] += W->RD[W->RD_position][n]; // * W[w].weight;
      W->RD[W->RD_position][n] = 0;
    }

    if(measure_R2) {
      W->r2old = W->R2pure[W->RD_position];
      W->z2old = W->Z2pure[W->RD_position];
    }

    if(measure_RDz) {
#pragma omp atomic
      RDz_pure.times_measured++;
       for(n=0; n<gridRD; n++) {
#pragma omp atomic
        RDz_pure.N[n] += W->RDz[W->RD_position][n]; // * W[w].weight;
        W->RDz[W->RD_position][n] = 0;
      }
    }
  }

#ifdef BC_SPHERE
  calculate_CM = ON;
#else
  if(R2_subtract_CM) calculate_CM = ON;
#endif
  if(calculate_CM) {
    for(i=0; i<N; i++) {
      xCM += W->x[i];
      yCM += W->y[i];
      zCM += W->z[i];
    }
    xCM /= (DOUBLE) N;
    yCM /= (DOUBLE) N;
    zCM /= (DOUBLE) N;
  }

  for(i=0; i<N; i++) {
#ifdef BC_ABSENT // trap
    if(R2_subtract_CM) {
      r2 = (W->x[i]-xCM)*(W->x[i]-xCM) + (W->y[i]-yCM)*(W->y[i]-yCM) + (W->z[i]-zCM)*(W->z[i]-zCM);
    }
    else {
      r2 = W->x[i]*W->x[i] + W->y[i]*W->y[i] + W->z[i]*W->z[i];
    }
#else // hom. system
    r2 = W->x[i]*W->x[i]; // save x direction
#endif

    r = Sqrt(r2);

#ifndef BC_SPHERE // on a sphere measure CM r2
#pragma omp atomic
    RD.r2 += r2;
#pragma omp atomic
    RD.r4 += r2*r2;
#pragma omp atomic
    RD.r2recent += r2;
    R2 += r2;
    R4 += r2*r2;
#endif

#ifdef BC_ABSENT // trap
    if(fabs(W->z[i])< RD.width) {
#else // hom. system
    if(r<RD.width) {
#endif
      n = (int) (r/RD.step);
      if(n<RD.size) {
#pragma omp atomic
        RD.N[n]++;
        if(MC == DIFFUSION) W->RD[W->RD_position][n]++; // pure estimator
      }
    }

    z2 = (W->z[i]-zCM)*(W->z[i]-zCM); // save z component
    //z2 = W->z[i]*W->z[i]; // save z component

#ifdef TRIAL_2D // save y component in 2D hom. system
#ifdef BC_2DPBC
    z2 = W->y[i]*W->y[i];
#endif
#ifdef BC_1DPBC_X
    z2 = W->y[i]*W->y[i];
#endif
#endif
    Z2 += z2;
    Z4 += z2*z2;
#ifdef BC_ABSENT
    if(measure_RDz && r < RDz.width) {
#else
    if(measure_RDz) {
#endif
      //RDz.r2 += z2;
      //RDz.r2recent += z2;
      //RDz.r4 += z2*z2;
      n = (int) (Sqrt(z2)/RDz.step);
      if(n<RDz.size) {
#pragma omp atomic
        RDz.N[n]++;
        if(MC == DIFFUSION) W->RDz[W->RD_position][n]++; // pure estimator
      }
    }

#pragma omp atomic
    RDz.r2 += Z2;
#pragma omp atomic
    RDz.r2recent += Z2;
#ifdef TRIAL_1D // define freq. of oscillations as Omega2 = 4.*z2 / (z4-N*z2*z2)
#pragma omp atomic
    RDz.r4 += Z2*Z2;
#else
#pragma omp atomic
    RDz.r4 += Z4;
#endif
  }

#ifdef BC_SPHERE // on a sphere calculate mean square CM size
  r2 = xCM*xCM + yCM*yCM + zCM*zCM;
#pragma omp atomic
  RD.r2 += r2;
#pragma omp atomic
  RD.r4 += r2*r2;
#pragma omp atomic
  RD.r2recent += r2;
#pragma omp atomic
  R2 += (DOUBLE) N * r2; // will be divided by N
#pragma omp atomic
  R4 += (DOUBLE) N * r2*r2;
#endif

  if(measure_RadDistrMATRIX) {
    for(i=0; i<N; i++) {
      nx = (int) (W->x[i]/RD_MATRIX.step_x);
      ny = (int) (W->y[i]/RD_MATRIX.step_y);

      if(nx<gridRDx && ny<gridRDy) {
#pragma omp atomic
        RD_MATRIX.N[nx][ny]++;
#pragma omp atomic
        RD_MATRIX.times_measured++;
      }
    }
  }

  if(measure_R2) {
    //W->r2old = W->r2; are measured in pure part
    //W->z2old = W->z2;
    W->r2 = R2 / (DOUBLE) N;
    W->r4 = R4 / (DOUBLE) N;
    W->z2 = Z2 / (DOUBLE) N;

    if(MC == DIFFUSION) { // pure estimator
      W->R2pure[W->RD_position] = R2 / (DOUBLE) N;
      W->Z2pure[W->RD_position] = Z2 / (DOUBLE) N;
    }

#ifdef TRIAL_1D // define freq. of oscillations as Omega2 = 4.*z2 / (z4-N*z2*z2)
    W->z4 = Z2*Z2 / (DOUBLE) N;
#else
    W->z4 = Z4 / (DOUBLE) N;
#endif
  }

  if(MC == DIFFUSION) { // move current position
    W->RD_position++;
    if(W->RD_position == grid_pure_block) W->RD_position = 0;
    if(W->RD_wait) { // skip the first block in PURE case
      W->RD_wait++;
      if(W->RD_wait == grid_pure_block +1) 
        W->RD_wait = OFF;
    }
  }
}

/************************ Walker Phase Anyons ********************************/
void WalkerPhaseAnyons(struct Walker *walker, double *RePhase, double *ImPhase) {
  int i,j;
  double ReS = 1., ImS = 0., ReSnew, ImSnew;
  int signz;
  double phi, cosphi, sinphi;
  double norm;

  for(i=0; i<N; i++) {
    for(j=i+1; j<N; j++) {
      signz = (walker->z[i]<walker->z[j]) ? 1 : -1;
      phi = -_PI_2*alpha_R*signz;
      cosphi = cos(phi);
      sinphi = sin(phi);
      ReSnew = ReS*cosphi - ImS*sinphi;
      ImSnew = ReS*sinphi + ImS*cosphi;
      ReS = ReSnew;
      ImS = ImSnew;
    }
  }
  norm = sqrt(ReS*ReS+ImS*ImS);
  *RePhase = ReS/norm;
  *ImPhase = ImS/norm;
}

/******************************** OBDM ***************************************/
// McMillan point is generated randomly in full 1D, 2D, 3D geometry
// in quasi 2D it is generated as a shift in (x,y) plane
// McMillan calculation
// the contribution of particle i is: exp(u_m - u_im[i] - u_ij[i])
// u_ij(i) = sum_j u(|r_i - r_j|)
// u_im(i) = u(|r_M - r_i|)
// u_m = sum_j u(|r_m - r_j|)
int MeasureOBDM(int w) {
  int i,j,m, n;
  int McMillan_spin = 0; // type of the particle to move
  DOUBLE *u_ij, *u_mi; // define locally for parallelization
  DOUBLE u_m;
  int sign_u_m;
  DOUBLE xm = 0.;
  DOUBLE ym = 0.;
  DOUBLE zm = 0.;
  DOUBLE u, r, r2, dr[3] = {0,0,0};
#ifdef ONE_BODY_TRIAL_TERMS
  DOUBLE U_one_body_old = 0;
  DOUBLE U_one_body_new = 0;
#endif
#ifdef SECURE
  int i_check;
  DOUBLE u_check;
#endif
#ifdef OBDM_FERMIONS
  int s;
#endif
#ifdef OBDM_ANYONS
  double ReS, ImS, ReSnew, ImSnew, signz1, signz2;
#endif
#ifdef OBDM_EXCITED_BOSONS
  int sExcited;
#endif
#ifdef BC_ABSENT
  int N1, N2;
#endif
#ifdef HARD_SPHERE
  int overlap;
#ifdef SECURE
  int  overlap_check;
#endif
#endif
#ifdef CRYSTAL // symmetric crystal contribution
#ifdef CRYSTAL_SYMMETRIC
#  define OBDM_CALCULATE_CRYSTAL_SYMMETRIC
#endif
#endif
#ifdef TRIAL_1D
#ifndef BC_ABSENT
  if(OBDM.max <= Lhalf/N) return MeasureOBDMShortRange();
#endif
#endif

#ifdef BC_Q2D // e.g. speckles
#ifdef SPECKLES
  return MeasureOBDMorbital();
#endif
#endif

#ifdef THREE_BODY_TERMS
  return MeasureOBDMfromU();
#endif

#ifdef MEASURE_PURE_OBDM
  if(MC == DIFFUSION && (iteration_global) % (Nmeasure*grid_pure_block) == 0 && iteration_global) {
    if(iteration_global>1 && measure_Nk_pure == OFF) SaveOBDMpure();
    if(iteration_global>1 && measure_Nk_pure == ON) SaveMomentumDistributionPure();
    MeasureOBDMpure();
  }
#endif

#ifdef CRYSTAL // symmetric crystal contribution
#ifdef CRYSTAL_SYMMETRIC
  for(i=0; i<N; i++) {
    Wp[w].x[i] = W[w].x[i];
    Wp[w].y[i] = W[w].y[i];
    Wp[w].z[i] = W[w].z[i];
  }
  U_one_body_old = OneBodyUWalker(Wp[w]);
#endif
#endif

  u_ij = calloc(N, sizeof(DOUBLE));
  u_mi = calloc(N, sizeof(DOUBLE));
  //u_ij = (DOUBLE*) Calloc("u_ij", N, sizeof(DOUBLE));
  //u_mi = (DOUBLE*) Calloc("u_mi", N, sizeof(DOUBLE));
  //ArrayEmpty1D(u_ij, i, N); // initialize arrays with zeros

  for(i=0; i<N; i++) {
    for(j=i+1; j<N; j++) {
      CaseX(dr[0] = W[w].x[i] - W[w].x[j]);
      CaseY(dr[1] = W[w].y[i] - W[w].y[j]);
      CaseZ(dr[2] = W[w].z[i] - W[w].z[j]);
      r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
#ifdef HARD_SPHERE
      if(r2<a2) {
        Warning("Measure OBDM, overlapping configuration at input, exiting\n");
        return -1;
    }
#endif
      if(CheckInteractionConditionWF(dr[0], dr[1], dr[2], r2)) {
        u = InterpolateU(&G, Sqrt(r2), W[w].spin[i], W[w].spin[j]);
        u_ij[i] += u;
        u_ij[j] += u;
      }
  }
  }

  // Generate McMillan point McMillan_points times
  for(m=0; m<McMillan_points; m++) {
#ifdef HARD_SPHERE
    overlap = OFF;
#endif

#ifdef BC_ABSENT
    CaseX(xm = (1.-2.*Random())*OBDM.max);
    CaseY(ym = (1.-2.*Random())*OBDM.max);
    CaseZ(zm = (1.-2.*Random())*OBDM.max);
#else // PBC
    CaseX(xm = Random()*L);
    CaseY(ym = Random()*L);
    CaseZ(zm = Random()*L);
#endif
    u_m = 0;
    sign_u_m = 1;

    //#ifdef ONE_BODY_TRIAL_TERMS // effect of the external field, in general can depend on index
    //      u_m += OneBodyU(xm, ym, zm);
    //#endif

    // fill u_m and u_im
    for(i=0; i<N; i++) {
      CaseX(dr[0] = W[w].x[i] - xm);
      CaseY(dr[1] = W[w].y[i] - ym);
      CaseZ(dr[2] = W[w].z[i] - zm);
      r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
#ifdef HARD_SPHERE // begin check
      if(r2<a2) {
        if(overlap == OFF) // so far no overlap with other particles (good conf. for moving particle i)
          overlap = i+1;
        else // overlapping with two particles (bad. conf)
          overlap = -1;
    }
      else {
#endif

        if(CheckInteractionConditionWF(dr[0], dr[1], dr[2], r2)) {
          u_mi[i] = InterpolateU(&G, Sqrt(r2), W[w].spin[i], McMillan_spin); // same type of particle
          u_m += u_mi[i];
        }
        else {
          u_mi[i] = 0.;
        }
#ifdef HARD_SPHERE // end check
        }
#endif
      }

    for(i=0; i<N; i++)
#ifdef SPINFULL
      if(W[w].spin[i] == McMillan_spin) // in SPINFULL case measure OBDM of type 1
#endif
      {
#ifdef HARD_SPHERE // begin overlap loop
        if(overlap ==  OFF || overlap == i+1) { // no overlap, or overlapping with i-th particle
#endif
          CaseX(dr[0] = W[w].x[i] - xm);
          CaseY(dr[1] = W[w].y[i] - ym);
          CaseZ(dr[2] = W[w].z[i] - zm);
          r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
          r = Sqrt(r2);

#ifdef BC_Q2D // z direction is integrated out assuming quasi 2D geometry
          n = (int) (Sqrt(dr[0]*dr[0]+dr[1]*dr[1])/OBDM.step);
#else // standart case
          n = (int) (r/OBDM.step);
#endif

          // if(CheckInteractionConditionWF(dr[0], dr[1], dr[2], r2)) {
          // assume that u(r) = 0, r>L/2 in the case of scalable code
#ifdef ONE_BODY_TRIAL_TERMS // one-body terms are present
#ifdef OBDM_CALCULATE_CRYSTAL_SYMMETRIC // cannot be presented as a one-body contribution
          Wp[w].x[i] = xm;
          Wp[w].y[i] = ym;
          Wp[w].z[i] = zm;
          U_one_body_new = OneBodyUWalker(Wp[w]);
          Wp[w].x[i] = W[w].x[i];
          Wp[w].y[i] = W[w].y[i];
          Wp[w].z[i] = W[w].z[i];

          // calculate Jastrow and Crystal contribution
          u = Exp(u_m - u_mi[i] - u_ij[i] + U_one_body_new - U_one_body_old);
#else // can be presented as a one-body contribution
          u = Exp(u_m - u_mi[i] - u_ij[i] + OneBodyU(xm,ym,zm,i, W[w].w) - OneBodyU(W[w].x[i], W[w].y[i], W[w].z[i],i, W[w].w));
#endif
#else // one-body terms are absent
          u = Exp(u_m - u_mi[i] - u_ij[i]);
#endif

#ifdef SECURE // direct evaluation of the OBDM
          u_check = 0;
#ifdef HARD_SPHERE
          overlap_check = OFF;
#endif
          for(i_check=0; i_check<N; i_check++) {
            if(i_check != i) {

              dr[0] = W[w].x[i_check] - xm;
              dr[1] = W[w].y[i_check] - ym;
              dr[2] = W[w].z[i_check] - zm;
              r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
#ifdef HARD_SPHERE
              if(r2<a2) overlap_check = ON;
#endif
              if(CheckInteractionConditionWF(dr[0], dr[1], dr[2], r2)) u_check += InterpolateU(&G, Sqrt(r2), W[w].spin[i_check], W[w].spin[i_check]);

              dr[0] = W[w].x[i_check] - W[w].x[i];
              dr[1] = W[w].y[i_check] - W[w].y[i];
              dr[2] = W[w].z[i_check] - W[w].z[i];
#ifdef HARD_SPHERE
              if(r2<a2) overlap_check = ON;
#endif
              r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
              if(CheckInteractionConditionWF(dr[0], dr[1], dr[2], r2)) u_check -= InterpolateU(&G, Sqrt(r2), W[w].spin[i_check], W[w].spin[i_check]);
            }
          }
          u_check += OneBodyU(xm, ym, zm, i, w) - OneBodyU(W[w].x[i], W[w].y[i], W[w].z[i], i, w);
          u_check = Exp(u_check);
#ifdef HARD_SPHERE
          if(overlap_check == ON)
            Warning("SECURE measure OBDM: overlapping configuration\n");
#endif
          /*if(fabs((u_check-u)/u)>1e-6)
          Warning("SECURE measure OBDM check failed!\n");
          if(fabs(u)>1e3)
          Warning("SECURE measure OBDM check passed, but strange value obtained: %" LE " !\n", u);*/
          u = u_check;//!!!
#endif

#ifdef OBDM_FERMIONS
#ifdef OBDM_ANYONS
          ReS = 1.;
          ImS = 0.;
          //WalkerPhaseAnyons(&W[w], &ReS, &ImS);
          for(j=0; j<N; j++) {
            if(j==i) continue;

            // psi(x') / psi(x)
            signz1 = (zm<W[w].z[j])?(1.):(-1.); // calculate sign sign(z[i]-z[j])
            signz2 = (W[w].z[i]<W[w].z[j])?(1.):(-1.); // calculate sign sign(z[i]-z[j])
            ReSnew = ReS*cos(-_PI_2*alpha_R*(signz1-signz2)) - ImS*sin(-_PI_2*alpha_R*(signz1-signz2));
            ImSnew = ReS*sin(-_PI_2*alpha_R*(signz1-signz2)) + ImS*cos(-_PI_2*alpha_R*(signz1-signz2));
            ReS = ReSnew;
            ImS = ImSnew;
          }
          //u *= ReS;
          u=u;
#else
          s = 1;
          for(j=0; j<N; j++) {
            if(j!=i) {
              if(W[w].z[i]<W[w].z[j]) s *= -1;
              if(zm<W[w].z[j]) s *= -1;
            }
          }
          u *= s;
#endif
#endif

#ifdef OBDM_EXCITED_BOSONS // change the sign for the excited states
          sExcited = 1;
          for(j=0; j<N; j++) {
            if(j!=i) {
              //for(i=0; i<N; i++) {
              //for(j=i+1; j<N; j++) {
              //dr[2] = W[w].z[i] - W[w].z[j];
              //r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
              //sExcited *= InterpolateExactSign(Sqrt(r2));
              sExcited *= InterpolateExactSign(fabs(W[w].z[i] - W[w].z[j]));
              //dr[2] = zm - W[w].z[j];
              //r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
              //sExcited *= InterpolateExactSign(Sqrt(r2));
              sExcited *= InterpolateExactSign(fabs(zm - W[w].z[j]));
            }
          }
          u *= sExcited;
#endif

          if(n<OBDM.size && W[w].status == ALIVE) {
            if(W[w].status == ALIVE || W[w].status == KILLED) {
              // one-body correlation function
#pragma omp atomic
              OBDM.N[n] += W[w].weight;
#ifndef OBDM_ANYONS
#ifdef OBDM_EXCITED_BOSONS
#pragma omp atomic
              OBDM.f[n] += fabs(u) * W[w].weight * sExcited;
#else
#pragma omp atomic
              OBDM.f[n] += fabs(u) * W[w].weight;
#endif
#endif

#ifdef SECURE
              if(OBDM.f[n]/OBDM.N[n]>100)
                Warning("MeasureOBDM: suspicious value %" LE "  = %" LE " /%" LE " \n", OBDM.f[n]/OBDM.N[n],OBDM.f[n],OBDM.N[n]);
#endif
#ifdef OBDM_FERMIONS
#ifdef OBDM_ANYONS
              if(W[w].z[i] < zm) {
#pragma omp atomic
                OBDM.f[n] +=      ReS * u * W[w].weight*2.;
#pragma omp atomic
                OBDMfermi.f[n] += ImS * u * W[w].weight*2.;
            }
#else
#pragma omp atomic
              OBDMfermi.f[n] += u * W[w].weight;
#endif
#endif
#ifdef TRIAL_1D
              for(j=0; j<OBDM.Nksize; j++) {
#pragma omp atomic
                OBDM.Nk[j] += Cos(OBDM.k[j]*r)*fabs(u)*W[w].weight;
#ifdef OBDM_FERMIONS
#ifdef OBDM_ANYONS
                //OBDMfermi.Nk[j] += (Cos(OBDM.k[j]*r)*ReS - Sin(OBDM.k[j]*r)*ImS)*u*W[w].weight;
#pragma omp atomic
                if(W[w].z[i] < zm) OBDMfermi.Nk[j] += 2.*(Cos(OBDM.k[j]*r)*ReS + Sin(OBDM.k[j]*r)*ImS)*u*W[w].weight;
#else
#pragma omp atomic
                OBDMfermi.Nk[j] += Cos(OBDM.k[j]*r)*u*W[w].weight;
#endif
#endif
              }
#pragma omp atomic
              OBDM.Nktimes_measured++;
#endif // 1D

#ifdef SECURE
              if(MC == VARIATIONAL && W[w].weight != 1.) Error("wrong weight in variational calculation");
#endif

#ifdef BC_ABSENT // reduced one-body correlation function
              N1 = (int) (fabs(W[w].z[i]) / RDz.step);
              N2 = (int) (fabs(zm) / RDz.step);
              if(N1<RDz.size && N2<RDz.size) {
#pragma omp atomic
                OBDMtrap.N[n] += Sqrt(RDz.f[N1]*RDz.f[N2]);
#pragma omp atomic
                OBDMtrap.f[n] += u;
          }
#endif
            }
#ifdef HARD_SPHERE // end overlap loop
        }
#endif
      }
    }
  }

  free(u_ij);
  free(u_mi);
  return 0;
}

/******************************** OBDM ***************************************/
// uses only U() call
int MeasureOBDMfromU(void) {
  int w;

#pragma omp parallel for
  for(w=0; w<Nwalkers; w++) MeasureOBDMfromUWalker(w);

  return 0;
}

int MeasureOBDMfromUWalker(int w) {
  int i,j,m,n,sign_displacement;
  DOUBLE u_old, u_new, u;
  DOUBLE r;
#ifdef HARD_SPHERE_HYPERRADIUS
  int overlap;
#endif
#ifdef OBDM_FERMIONS
  int s;
#endif

  u_old = U(W[w]);
  for(i=0; i<N; i++) { // loop over particle which is moved
    for(m=0; m<McMillan_points; m++) { // loop over McMillan points
#ifdef HARD_SPHERE_HYPERRADIUS // begin overlap loop
    overlap = ON;
    while(overlap == ON) { // repeat until no overlap
#endif
      CopyWalkerCoord(&Wp[w], &W[w]);

      n = rand() % OBDM.size; // displacement integer
      r = n*OBDM.step; // displacement distance
      r /= sqrt(DIMENSION); // displacement amplitude

      // displace particle i
      CaseX(sign_displacement = 2*(rand() % 2)-1);
      CaseX(Wp[w].x[i] += (DOUBLE)sign_displacement * r);
      CaseY(sign_displacement = 2*(rand() % 2)-1);
      CaseY(Wp[w].y[i] += (DOUBLE)sign_displacement * r);
      CaseZ(sign_displacement = 2*(rand() % 2)-1);
      CaseZ(Wp[w].z[i] += (DOUBLE)sign_displacement * r);
      ReduceToTheBoxXYZ(&Wp[w].x[i], &Wp[w].y[i], &Wp[w].z[i]);

#ifdef HARD_SPHERE_HYPERRADIUS
      overlap = CheckWalkerHyperradiusOverlapping(Wp[w]);
    } // end overlap loop
#endif

      u = exp(U(Wp[w]) - u_old);

#ifdef OBDM_FERMIONS
#ifdef TRIAL_1D
      s = 1;
      for(j=0; j<N; j++) {
        if(j!=i) {
          if(W[w].z[i]<W[w].z[j]) s *= -1;
          if(Wp[w].z[i]<W[w].z[j]) s *= -1;
        }
      }
      u *= s;
#endif
#endif
      if(n<OBDM.size && W[w].status == ALIVE) {
        if(W[w].status == ALIVE || W[w].status == KILLED) {
          // one-body correlation function
          OBDM.N[n] += W[w].weight;
          OBDM.f[n] += fabs(u) * W[w].weight;
#ifdef SECURE
          if(OBDM.f[n]/OBDM.N[n]>100)
            Warning("MeasureOBDM: suspicious value %" LE "  = %" LE " /%" LE " \n", OBDM.f[n]/OBDM.N[n],OBDM.f[n],OBDM.N[n]);
#endif
#ifdef OBDM_FERMIONS
          OBDMfermi.f[n] += u * W[w].weight;
#endif
#ifdef TRIAL_1D
          for(j=0; j<OBDM.Nksize; j++) {
            OBDM.Nk[j] += Cos(OBDM.k[j]*r)*fabs(u)*W[w].weight;
#ifdef OBDM_FERMIONS
            OBDMfermi.Nk[j] += Cos(OBDM.k[j]*r)*u*W[w].weight;
#endif
          }
          OBDM.Nktimes_measured++;
#endif
        }
      }
    }
  }
  return 0;
}

/********************* Measure OBDM Pure *************************************/
int MeasureOBDMpure(void) {
#ifdef MEASURE_PURE_OBDM
  int w,j,w_v,i_moved;
  DOUBLE xm, ym, zm;
  DOUBLE u, r2, dr[3] = {0,0,0};
#ifdef OBDM_FERMIONS
  int s;
#endif
  int terms = 0;
  DOUBLE shift = 0.5; //0.45; // 0.25132 // in units of L
#ifndef TRIAL_1D // in 2D and 3D
  DOUBLE theta;
#endif
#ifdef TRIAL_3D
  DOUBLE phi;
#endif
#ifdef ONE_BODY_TRIAL_TERMS
  DOUBLE u_one_body_old, u_one_body_new;
#endif

  Nwalkers_pure = Npop_virtual;

  //for(w=0; w<Nwalkers && terms<Nwalkers_pure; w++) {
  //if(W[w].status != ALIVE)  Warning("  not ALIVE walker encountered in MeasureOBDMpure\n");
  //terms ++;
  for(terms=0; terms<Nwalkers_pure; terms++) {
    w = (int)(Random()*(DOUBLE)Nwalkers); // choose random walker to move

    xm = ym = zm = 0;
    u = 0; // initialize weight
    i_moved = (int)(Random()*(DOUBLE)N); // number of particle which will be moved r[i] -> rm

#ifdef TRIAL_3D
    if(measure_Nk_pure) {
      xm = Random()*L;
      ym = Random()*L;
      zm = Random()*L;
    }
    else {
      theta = Random()*2.*PI;
      phi = Random()*PI;

      xm = W[w].x[i_moved] + shift*L*Cos(theta)*Sin(phi);
      ym = W[w].y[i_moved] + shift*L*Sin(theta)*Sin(phi);
      zm = W[w].z[i_moved] + shift*L*Cos(phi);

      if(xm>L) xm -= L;
      if(xm<0) xm += L;
      if(ym>L) ym -= L;
      if(ym<0) ym += L;
      if(zm>L) zm -= L;
      if(zm<0) zm += L;

      if(xm<0 || xm>L || ym<0 || ym>L || zm<0 || zm>L) Error("  McMillan point out of the box!\n");
    }
#endif

#ifdef TRIAL_2D // generate 2D displacement
    if(measure_Nk_pure) {
      xm = Random()*L;
      ym = Random()*L;
      zm = 0;
    }
    else {
      theta = Random()*2.*pi;
      xm = W[w].x[i_moved] + shift*L*Cos(theta);
      ym = W[w].y[i_moved] + shift*L*Sin(theta);
      zm = W[w].z[i_moved];
      if(xm>L) xm -= L;
      if(xm<0) xm += L;
      if(ym>L) ym -= L;
      if(ym<0) ym += L;
    }
#endif

#ifdef TRIAL_1D // generate 1D displacement
    if(measure_Nk_pure) {
       CaseX(xm = 0);
       CaseY(ym = 0;)
       CaseZ(zm = Random()*L);
    }
    else {
      CaseX(xm = W[w].x[i]);
      CaseY(ym = W[w].y[i]);
      zm = W[w].z[i_moved] + shift*L;
      if(zm>L) zm -= L; // put back to the box
    }
#endif

#ifdef BC_ABSENT
    Message("Pure OBDM estimator is not implemented!\n");
#endif

    for(j=0; j<N; j++) {
      if(j != i_moved) {
        dr[0] = W[w].x[j] - xm;
        dr[1] = W[w].y[j] - ym;
        dr[2] = W[w].z[j] - zm;
        r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
        if(CheckInteractionConditionWF(dr[0], dr[1], dr[2], r2)) u += InterpolateU(&G, Sqrt(r2), W[w].spin[j], W[w].spin[i_moved]);

        dr[0] = W[w].x[j] - W[w].x[i_moved];
        dr[1] = W[w].y[j] - W[w].y[i_moved];
        dr[2] = W[w].z[j] - W[w].z[i_moved];
        r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
        if(CheckInteractionConditionWF(dr[0], dr[1], dr[2], r2)) u -= InterpolateU(&G, Sqrt(r2), W[w].spin[j], W[w].spin[i_moved]);
      }
    }
    u = Exp(u); // calculate Jastrow contribution

#ifdef ONE_BODY_TRIAL_TERMS // one-body terms
    for(j=0; j<N; j++) {
      Wp[w].x[j] = W[w].x[j];
      Wp[w].y[j] = W[w].y[j];
      Wp[w].z[j] = W[w].z[j];
    }
    u_one_body_old = OneBodyUWalker(Wp[w]);
    Wp[w].x[i_moved] = xm;
    Wp[w].y[i_moved] = ym;
    Wp[w].z[i_moved] = zm;
    u_one_body_new = OneBodyUWalker(Wp[w]);
    u *= Exp(u_one_body_new - u_one_body_old); // calculate one-body contribution
#endif

#ifdef OBDM_FERMIONS
    s = 1;
    for(j=0; j<N; j++) {
      if(j != i_moved) {
        if(W[w].z[i_moved] < W[w].z[j]) s *= -1;
        if(zm < W[w].z[j]) s *= -1;
      }
    }
    u *= s;
#endif

    // create a virtual walker
    w_v = dead_walkers_storage[Ndead_walkers-1]; // take it from dead walkers storage
    if(Ndead_walkers < 1) {
      Warning("  MeasureOBDMpure: too many walkers\n");
      return 1;
    }
    CopyWalker(&W[w_v], &W[w]); // inherit data

    W[w_v].x[i_moved] = xm; // correct the new partice position
    W[w_v].y[i_moved] = ym;
    W[w_v].z[i_moved] = zm;

    W[w_v].R_dmc_end_of_step[i_moved][0] = xm; // quadratic algorithm restores the coordinates
    W[w_v].R_dmc_end_of_step[i_moved][1] = ym;
    W[w_v].R_dmc_end_of_step[i_moved][2] = zm;

    W[w_v].OBDMpureB = fabs(u);
    W[w_v].OBDMpureF = u;
    W[w_v].OBDM_position = 0;
    W[w_v].status = VIRTUAL;
    W[w_v].weight = 0.; // virtual walker stores logarithmic weight

    dr[0] = W[w].x[i_moved] - xm;
    dr[1] = W[w].y[i_moved] - ym;
    dr[2] = W[w].z[i_moved] - zm;
    r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
    W[w_v].OBDMpure_r = Sqrt(r2);

    Ndead_walkers--;
  }
#endif
  return 0;
}

/******************************** OBDM Q2D ***********************************/
// transverse orbital should be explicitly normalized as
// (2 alpha / pi)^{1/4} exp(- alpha r^2)
int MeasureOBDMorbital(void) {
  int r_3D = OFF; // distance r in (x,y,z), otherwize in (x,y)
  int zm_same_as_z = ON; // keep z'=z, otherwise integrate over z and z'
#ifdef BC_Q2D
  DOUBLE Lintz = 3.; // integration [-Lintz, Lintz]
  //DOUBLE dummy, sigma = 0.5;
#endif
  int w,i,m,n;
  int i_moved;
  DOUBLE u_m;
  DOUBLE xm,ym,zm;
  DOUBLE xold,yold,zold;
  DOUBLE u;
  DOUBLE dr[3] = {0,0,0},r2;
#ifdef SPECKLES
  DOUBLE phir, phirp;
  DOUBLE xperiod, yperiod;
  int ix,iy;
  DOUBLE alatt;
#endif
#ifdef HARD_SPHERE
  int overlap;
#endif
  DOUBLE u_one_body_r, u_one_body_m;
  static int times = 0;

#ifdef CRYSTAL
  Error("MeasureOBDMorbital is not implemented\n");
#endif

  if(times == 0) {
#ifdef BC_Q2D
        Warning("OBDM Q2D parameters:\n");
    if(zm_same_as_z)
      Warning("  McMillan points are generated with z'=z\n");
    else
      Warning("  2D OBDM is obtained by integration over over z and z' in a box (%lf < z <%lf)\n", -Lintz, Lintz);

    if(r_3D)
      Warning("  distance r in OBDM is 3D distance in (x,y,z) variables\n");
    else
      Warning("  distance r in OBDM is 2D distance in (x,y) variables\n");
#else
    Warning("  Calling OBDM Q2D in a different geometry\n");
#endif
  }

  times++;
  for(w=0; w<Nwalkers; w++) {
    for(m=0; m<McMillan_points; m++) {
#ifdef HARD_SPHERE
      overlap = OFF;
#endif
      i_moved = rand() % N;

      xold = W[w].x[i_moved];
      yold = W[w].y[i_moved];
      zold = W[w].z[i_moved];

      xm = Random()*L;
      ym = Random()*L;

#ifdef BC_Q2D
      if(zm_same_as_z)
        zm = zold;
      else
        zm = (Random()-0.5)*2.*Lintz; // rectangular box integration
#else
        zm = zold;
#endif

      dr[0] = xm - xold;
      dr[1] = ym - yold;
      dr[2] = zm - zold;
      r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
     if(r_3D) {
        n = (int) (Sqrt(r2)/OBDM.step); // r is 3D distance
     }
     else {
       n = (int) (Sqrt(dr[0]*dr[0]+dr[1]*dr[1])/OBDM.step); // r is 2D distance in (x,y) plane
     }

      u_m = 0;

      // effect of the external field
      u_one_body_m = OneBodyU(xm, ym, zm, 0, W[w].w);
      u_one_body_r = OneBodyU(xold, yold, zold, 0, W[w].w);
      u_m += u_one_body_m - u_one_body_r;

      for(i=0; i<N; i++) {
        if(i != i_moved) { // new value of the w.f.
          dr[0] = W[w].x[i] - xm;
          dr[1] = W[w].y[i] - ym;
          dr[2] = W[w].z[i] - zm;
          r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
#ifdef HARD_SPHERE
          if(r2<a2) {
            overlap = ON;
          }
          else {
#endif
          if(CheckInteractionConditionWF(dr[0], dr[1], dr[2], r2)) {
            u_m += InterpolateU(&G, Sqrt(r2), W[w].spin[i], W[w].spin[i]);
          }
#ifdef HARD_SPHERE
          }
#endif

          dr[0] = W[w].x[i] - xold;
          dr[1] = W[w].y[i] - yold;
          dr[2] = W[w].z[i] - zold;
          r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
#ifdef HARD_SPHERE
          if(r2<a2) {
            overlap = ON;
          }
          else {
#endif
          if(CheckInteractionConditionWF(dr[0], dr[1], dr[2], r2)) {
            u_m -= InterpolateU(&G, Sqrt(r2), W[w].spin[i], W[w].spin[i]);
          }
#ifdef HARD_SPHERE
          }
#endif
        }
      }
      u = Exp(u_m);

      if(u>1e2 && verbosity) {
        Warning(" OBDM measurement: strange (large) value\n");
        u = 1e2;
      }

#ifdef HARD_SPHERE
      if((W[w].status == ALIVE || W[w].status == KILLED) && !overlap) {
#else
      if(W[w].status == ALIVE || W[w].status == KILLED) {
#endif

        // measure OBDM
        if(n<OBDM.size) {
          OBDM.N[n] += W[w].weight;
#ifndef BC_Q2D // usual case
          OBDM.f[n] += u * W[w].weight;
#else   // Quasi 2D
          if(zm_same_as_z)
            OBDM.f[n] += u * W[w].weight;
          else // rectangular box integration, explicit normalization of z orbita, \int phi^2(r) dr = sqrt(pi/2alpha)
            OBDM.f[n] += u*Lintz*sqrt(2.*alpha_z/PI)*W[w].weight;
#endif
        }

#ifdef SPECKLES
        if(measure_OBDM_MATRIX) {
          //phir   = exp(OneBodyU(xold, yold, zold));
          //phirp  = exp(OneBodyU(xm,ym,zm));
          phir   = exp(u_one_body_r);
          phirp  = exp(u_one_body_m);
          phir  /= exp(-(alpha_z-0.5)*zold*zold); // condensate orbital is assumed to be (-1/2 z^2)
          phirp /= exp(-(alpha_z-0.5)*zm*zm);
#ifdef TRIAL_2D
          OBDM.CF += u*phir*phirp*L*L;
#else // Quasi 2D or 3D 
          OBDM.CF += u*phir*phirp*L*L*2.*Lintz*sqrt(2.*alpha_z/PI); // rectangular box integration
          //OBDM.CF += u*phir*phirp*L*L*sigma*sqrt(2.*PI)/sqrt(0.5*PI/alpha_z); // Gaussian integration
#endif
          OBDM.CFN++;

          // measure condensate wave function as an iterated eigenvector
          alatt = L/sqrt((DOUBLE)N); //alatt = 1./sqrt(dens); // lattice period
          //xperiod = xold - alatt*(DOUBLE)((int)(xold/alatt));!!!
          //yperiod = yold - alatt*(DOUBLE)((int)(yold/alatt));
          xperiod = xm - alatt*(DOUBLE)((int)(xm/alatt));
          yperiod = ym - alatt*(DOUBLE)((int)(ym/alatt));

          ix = (int) (xperiod/OBDM_MATRIX.step+0.5); // reduce to period
          iy = (int) (yperiod/OBDM_MATRIX.step+0.5);
          if(ix==OBDM_MATRIX.size) ix = 0;
          if(iy==OBDM_MATRIX.size) iy = 0;
          if(ix>OBDM_MATRIX.size) 
            Error("OBDM MATRIX out of the box (x = %lf, reduced to period %lf, matrix size %lf), i= , size = !", xold, xperiod, OBDM_MATRIX.max, ix, OBDM_MATRIX.size);
          if(iy>OBDM_MATRIX.size) 
            Error("OBDM MATRIX out of the box (y = %lf, reduced to period %lf, matrix size %lf), i= %i, size = %i", yold, yperiod, OBDM_MATRIX.max, ix, OBDM_MATRIX.size);
          //OBDM_MATRIX.f[ix][iy] += u*phirp*sqrt(2.*alpha_z/PI);!!!
          OBDM_MATRIX.f[ix][iy] += u*phir*sqrt(2.*alpha_z/PI);
          //OBDM_MATRIX.N[ix][iy] ++;
          OBDM_MATRIX.times_measured++;
        }
#endif
      }
    }
  }

  return 0;
}

/********************* Measure OBDM Short Range ******************************/
int MeasureOBDMShortRange(void) {
// optimized calculation of the short range of the OBDM
  int w,i,m, n;
  int i_moved;
  DOUBLE u_m;
  DOUBLE zm;
  DOUBLE zold;
  DOUBLE u;
#ifdef OBDM_FERMIONS
  int s;
#endif

  for(w=0; w<Nwalkers; w++) {
    for(m=0; m<McMillan_points; m++) {
      i_moved = rand() % N;
      zold = W[w].z[i_moved];

      do {
        n = rand() % OBDM.size;
        zm = zold + n * OBDM.step*((rand()%2==0)?1:-1);
      }
      while(zm<0 || zm>L);

      if(zm<0 || zm>L) Error("OBDM short range error");

      u_m = 0;
#ifdef BC_ABSENT // effect of the external field
      u_m += OneBodyU(0, 0, zm, i_moved, w) - OneBodyU(0, 0, zold, i_moved, w);
#endif 
      for(i=0; i<N; i++) {
        if(i != i_moved) { // new value of the w.f.
          u_m += InterpolateU(&G, fabs(W[w].z[i] - zm), W[w].spin[i], W[w].spin[i]);
          u_m -= InterpolateU(&G, fabs(W[w].z[i] - zold), W[w].spin[i], W[w].spin[i]);
        }
      }
      u = Exp(u_m);

      if(W[w].status == ALIVE || W[w].status == KILLED) {
        OBDM.N[n] += W[w].weight;
        OBDM.f[n] += u * W[w].weight;
#ifdef OBDM_FERMIONS
        s = 1;
        for(i=0; i<N; i++) {
          if(i != i_moved) {
            if(zold<W[w].z[i]) s *= -1;
            if(zm < W[w].z[i]) s *= -1;
          }
        }
        OBDMfermi.f[n] += u * s * W[w].weight;
#endif
      }
    }
  }
  return 0;
}

// function orders particles and return the parity of the permutation
int PermutationParity(DOUBLE *r) {
  int i,j;
  int parity = 1;
  DOUBLE x;

  for(i=0; i<N; i++) {
    for(j=i+1; j<N; j++) {
      if(r[j]<r[i]) {
        x = r[i];
        r[i] = r[j];
        r[j] = x;
        parity *= -1;
      }
    }
  }

  return parity;
}

/******************************** OBDM ***************************************/
#ifdef TRIAL_1D
void MeasureOBDMMatrix(void) {
  int w,i,j,m,n1,n2;
  DOUBLE u_m;
  int sign_u_m;
  DOUBLE xm, ym, zm;
  DOUBLE u, r2, dr[3] = {0,0,0};
#ifdef OBDM_FERMIONS
  int s=1;
#endif
#ifdef OBDM_EXCITED_BOSONS
  int s1=1, s2;
#endif
  DOUBLE zCM = 0.;

  if(OBDM_subtract_CM) { //  calculate center of mass position
    for(i=0; i<N; i++) {
      zCM += W->z[i];
    }
    zCM /= (DOUBLE) N;
  }

#ifdef CRYSTAL
  Warning("  MeasureOBDMMatrix is not implemented for a crystal. Skipping.\n");
  return;
#endif

  for(w=0; w<Nwalkers; w++) { // fill u_ij
    for(i=0; i<N; i++) { // initialize arrays with zeros
      u_ij[i] = 0;
    }
    for(i=0; i<N; i++) {
      for(j=i+1; j<N; j++) {
        CaseX(dr[0] = W[w].x[i] - W[w].x[j]);
        CaseY(dr[1] = W[w].y[i] - W[w].y[j]);
        CaseZ(dr[2] = W[w].z[i] - W[w].z[j]);
        r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
        if(CheckInteractionConditionWF(dr[0], dr[1], dr[2], r2)) {
          u = InterpolateU(&G, Sqrt(r2), W[w].spin[i], W[w].spin[j]);
          u_ij[i] += u;
          u_ij[j] += u;
#ifdef OBDM_EXCITED_BOSONS // accumulate sign for excited states
          s1 *= InterpolateExactSign(Sqrt(r2));
#endif
        }
      }
    }

    // Generate McMillan point McMillan_points times
    for(m=0; m<OBDM_MATRIX.size; m++) {
      xm = 0;
      ym = 0;
      if(OBDM_subtract_CM)
        zm = OBDM_MATRIX.min + (DOUBLE)m *OBDM_MATRIX.step + zCM;
      else
        zm = OBDM_MATRIX.min + (DOUBLE)m *OBDM_MATRIX.step;

      sign_u_m = 1;
      u_m = OneBodyU(xm, ym, zm, 0, w);

      // fill u_m and u_im
      for(i=0; i<N; i++) {
        dr[2] = W[w].z[i] - zm;
        r2 = dr[2]*dr[2];
        if(CheckInteractionConditionWF(dr[0], dr[1], dr[2], r2)) {
          u_mi[i] = InterpolateU(&G, Sqrt(r2), W[w].spin[i], W[w].spin[i]);
          u_m += u_mi[i];
        }
      }

      for(i=0; i<N; i++) { // Move i -> m
#ifdef FIRST_PARTICLE_IS_NOT_MOVED
        if(i==0) i++; // skip the first particle
#endif
        if(R2_subtract_CM) {
          n1 = (int) ((W[w].z[i]-OBDM_MATRIX.min-zCM)/OBDM_MATRIX.step);
        }
        else {
          n1 = (int) ((W[w].z[i]-OBDM_MATRIX.min)/OBDM_MATRIX.step);
        }
        n2 = m;

        if(n1>=0 && n1<OBDM_MATRIX.size) {
          u = Exp(u_m - u_mi[i] - u_ij[i] - OneBodyU(W[w].x[i], W[w].y[i], W[w].z[i], i, w));
#ifdef OBDM_FERMIONS
          s = 1;
          for(j=0; j<N; j++) {
            if(j!=i) {
              if(W[w].z[i]<W[w].z[j]) s *= -1;
              if(zm<W[w].z[j]) s *= -1;
            }
          }
          //u *= s;
#endif
#ifdef OBDM_EXCITED_BOSONS // change the sign for the excited states
          s2 = 1;
          for(i=0; i<N; i++) {
            for(j=i+1; j<N; j++) {
              dr[2] = zm - W[w].z[j];
              r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
              s2 *= InterpolateExactSign(Sqrt(r2));
            }
          }
          u=u*s2;
#endif
          if(n1<0) Warning(" negative index n1 in g1(z1,z2)\n");
          if(n2<0) Warning(" negative index n2 in g1(z1,z2)\n");
          OBDM_MATRIX.f[n1][n2] += u;
          OBDM_MATRIX.f[n2][n1] += u;
#ifdef OBDM_FERMIONS
          OBDMfermi_MATRIX.f[n1][n2] += u*s;
          OBDMfermi_MATRIX.f[n2][n1] += u*s;
#endif
          OBDM_MATRIX.times_measured++;

          // measure momentum distribution
          for(j=0; j<OBDM.Nksize; j++) {
            OBDM.Nk[j] += Cos(OBDM.k[j]*(W[w].z[i]-zm))*u*W[w].weight;
#ifdef OBDM_FERMIONS
            OBDMfermi.Nk[j] += Cos(OBDM.k[j]*(W[w].z[i]-zm))*u*W[w].weight*s;
#endif
          }
          OBDM.Nktimes_measured++;
        }
      }
    }
  }
}
#else // TRIAL_2D
DOUBLE* mu_ij = NULL;
void MeasureOBDMMatrix(void) {
  int i, j, w, k, n;
  DOUBLE h, xp, yp, rp, smu_k, mu, r2, phi;
  DOUBLE alpha;

#ifdef BC_Q2D // e.g. speckles
  MeasureOBDMorbital();
  return;
#endif

  if(alpha_x != alpha_y) {
    Warning("  Measure OBDM matrix: alpha_x must be equal to alpha_y\n");
    return;
  }
  else {
    alpha = alpha_x;
  }

  if (!mu_ij) mu_ij = (DOUBLE*)malloc(N*sizeof(DOUBLE));
  h = OBDM_MATRIX.step;

  for(w=0; w<Nwalkers; w++) {
    OBDM_MATRIX.times_measured++;

    memset(mu_ij, 0x00, N*sizeof(DOUBLE));

    for(i=0; i<N; i++) {
      xp = W[w].x[i];
      yp = W[w].y[i];

      mu_p_k[i] = xp*xp+yp*yp;
      mu = -alpha*mu_p_k[i]; // mu = - alpha r^2
      mu_ij[i] += mu;
      mu_p_k[i] = Sqrt(mu_p_k[i]); // mu_p_k[i] = r(i)

      // mu_ij = u1(i) + \sum_{j \ne i} u2(ij)
      for(j=i+1; j<N; j++) {
        r2 = (xp - W[w].x[j])*(xp - W[w].x[j]) + (yp - W[w].y[j])*(yp - W[w].y[j]);
        mu = InterpolateU(&G, Sqrt(r2), W[w].spin[j], W[w].spin[j]);
        mu_ij[i] += mu;
        mu_ij[j] += mu;
      }
    }

    for(i=0; i<OBDM_MATRIX.size; i++) {
      phi = Random() * _2PI;

      //xp = yp = rp = h * (i+1);
      xp = yp = rp = h * ((DOUBLE)i+0.5);
      xp *= Cos(phi);
      yp *= Sin(phi);

      //smu_k = -alpha*h*h*(DOUBLE)(i*i);
      smu_k = -alpha*rp*rp;
      for(n=0; n<N; n++) {
        r2 = (xp - W[w].x[n])*(xp - W[w].x[n]) + (yp - W[w].y[n])*(yp - W[w].y[n]);
        mu_k[n] = InterpolateU(&G, Sqrt(r2), W[w].spin[n], W[w].spin[n]); // mu_k[j] = u2(i'j)
        smu_k += mu_k[n]; // smu_k = u1(i') + \sum_j u2(i'j)
      }
      for(n=0; n<N; n++) { //  particle n -> McMillan (rp)
        k = (int)(mu_p_k[n]/h); // bin to the left
        if(k < OBDM_MATRIX.size) {
          // mu_ij = u1(i) + \sum_{j \ne i} u2(ij)
          mu = Exp(smu_k - mu_k[n] - mu_ij[n])/(_2PI*_2PI*mu_p_k[n]*h); // psi'/psi dr/r, probability of this configuratioin is prop. to r
          //mu = Exp(smu_k - mu_k[n] - mu_ij[n])/Sqrt(mu_p_k[n]*rp);
          //mu = Exp(smu_k - mu_k[n] - mu_ij[n]);
          // [u1(i') + \sum_j u2(i'j)] - u2(i'j) - [u1(i) + \sum_{j \ne i} u2(ij)]
          // Exp(u) / r(i)
          OBDM_MATRIX.f[i][k] += mu;
          OBDM_MATRIX.f[k][i] += mu;
          OBDM_MATRIX.N[i][k] ++;
          OBDM_MATRIX.N[k][i] ++;
        }
      }
    }
  }
}
#endif

/******************************** TBDM ***************************************/
int MeasureTBDM(void) {
  DOUBLE CMshift[3] = {0.,0.,0.}; // shift of the center of the mass
  int i,j,w,nR;
  int index_up, index_down; // this pair will be moved
  DOUBLE dr[3] = {0,0,0},r2,CM2,u;
  int m,n1,n2;
#ifdef OBDM_FERMIONS
  int s;
#endif
#ifdef SECURE
  DOUBLE u_check;
#endif

#ifdef THREE_BODY_TERMS
  return MeasureTBDMfromU();
#endif

#ifdef MEASURE_PURE_TBDM
  if(MC == DIFFUSION && (iteration_global) % (Nmeasure*grid_pure_block) == 0 && iteration_global) {
    if(iteration_global>1) SaveOBDMpure(); // the same function as in save OBDM pure
    MeasureTBDMpure();
  }
#endif

  for(w=0; w<Nwalkers; w++) { // loop over walkers
    for(i=0; i<N; i++) { // initialize arrays with zeros
      u_tbdm_i[i] = 0.;
    }
    for(i=0; i<N; i++) {
      for(j=i+1; j<N; j++) {
        dr[0] = W[w].x[i] - W[w].x[j];
        dr[1] = W[w].y[i] - W[w].y[j];
        dr[2] = W[w].z[i] - W[w].z[j];
        r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
        //if(CheckInteractionConditionWF(dr[0], dr[1], dr[2], r2)) {
          u_tbdm_ij[i][j] = InterpolateU(&G, Sqrt(r2), W[w].spin[i], W[w].spin[j]); // store contribution of pair (i,j)
          u_tbdm_i[i] += u_tbdm_ij[i][j]; // store total contribution of particle i
          u_tbdm_i[j] += u_tbdm_ij[i][j];
        //}
      }
    }

    // Generate McMillan shift McMillan_points times
    for(m=0; m<McMillanTBDM_points; m++) {
      do { // calculate center of mass separation distance
#ifdef TRIAL_3D
        CMshift[0] = (2.*Random()-1.)*TBDM.max;
        CMshift[1] = (2.*Random()-1.)*TBDM.max;
        CMshift[2] = (2.*Random()-1.)*TBDM.max;
#endif
#ifdef TRIAL_2D
        CMshift[0] = (2.*Random()-1.)*TBDM.max;
        CMshift[1] = (2.*Random()-1.)*TBDM.max;
#endif
#ifdef TRIAL_1D
        CMshift[2] = (2.*Random()-1.)*TBDM.max;
#endif
        CM2 = FindNearestImage(&CMshift[0], &CMshift[1], &CMshift[2]);
        nR = (int) (Sqrt(CM2)/TBDM.step+0.);
      }
      while(CM2>=TBDM.max*TBDM.max);

      for(i=0; i<N; i++) { // calculate shifted position of each particle
        Wp[w].x[i] = W[w].x[i] + CMshift[0];
        Wp[w].y[i] = W[w].y[i] + CMshift[1];
        Wp[w].z[i] = W[w].z[i] + CMshift[2];
      }
      ReduceWalkerToTheBox(&Wp[w]);

      for(i=0; i<N; i++) u_tbdm_mcmillan_i[i] = 0.; // empty arrays

      for(i=0; i<N; i++) { // interaction of the new McMillan position
        for(j=0; j<N; j++) {
          if(j != i) {
            dr[0] = Wp[w].x[i] - W[w].x[j];
            dr[1] = Wp[w].y[i] - W[w].y[j];
            dr[2] = Wp[w].z[i] - W[w].z[j];
            r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
            //if(CheckInteractionConditionWF(dr[0], dr[1], dr[2], r2)) {
              u_tbdm_mcmillan_ij[i][j] = InterpolateU(&G, Sqrt(r2), W[w].spin[i], W[w].spin[j]); // pair contribution from McMillan point "i" and particle "j"
              u_tbdm_mcmillan_i[i] += u_tbdm_mcmillan_ij[i][j]; // total contribution from McMillan point "i"
            //}
          }
        }
      }

      // Shift particles (up) and (down) by McMillan shift
      TBDM_MATRIX.times_measured++;
      for(index_up=0; index_up<N; index_up++) {
        for(index_down=index_up+1; index_down<N; index_down++) {

          u = Exp( u_tbdm_mcmillan_i[index_up]  -u_tbdm_i[index_up]
                  +u_tbdm_mcmillan_i[index_down]-u_tbdm_i[index_down]
                   -u_tbdm_mcmillan_ij[index_up][index_down]-u_tbdm_mcmillan_ij[index_down][index_up]+2.*u_tbdm_ij[index_up][index_down]);

#ifdef SECURE // direct evaluation of the OBDM
          u_check = 0;
          for(i=0; i<N; i++) {
            if(i != index_up && i != index_down) {
              dr[0] = W[w].x[i] - Wp[w].x[index_up];
              dr[1] = W[w].y[i] - Wp[w].y[index_up];
              dr[2] = W[w].z[i] - Wp[w].z[index_up];
              r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
              if(CheckInteractionConditionWF(dr[0], dr[1], dr[2], r2)) u_check += InterpolateU(&G, Sqrt(r2), W[w].spin[i], W[w].spin[index_up]);
              dr[0] = W[w].x[i] - W[w].x[index_up];
              dr[1] = W[w].y[i] - W[w].y[index_up];
              dr[2] = W[w].z[i] - W[w].z[index_up];
              r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
              if(CheckInteractionConditionWF(dr[0], dr[1], dr[2], r2)) u_check -= InterpolateU(&G, Sqrt(r2), W[w].spin[i], W[w].spin[index_up]);
              dr[0] = W[w].x[i] - Wp[w].x[index_down];
              dr[1] = W[w].y[i] - Wp[w].y[index_down];
              dr[2] = W[w].z[i] - Wp[w].z[index_down];
              r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
              if(CheckInteractionConditionWF(dr[0], dr[1], dr[2], r2)) u_check += InterpolateU(&G, Sqrt(r2), W[w].spin[i], W[w].spin[index_down]);
              dr[0] = W[w].x[i] - W[w].x[index_down];
              dr[1] = W[w].y[i] - W[w].y[index_down];
              dr[2] = W[w].z[i] - W[w].z[index_down];
              r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
              if(CheckInteractionConditionWF(dr[0], dr[1], dr[2], r2)) u_check -= InterpolateU(&G, Sqrt(r2), W[w].spin[i], W[w].spin[index_down]);
            }
          }
          u_check = Exp(u_check);
          if(fabs((u_check-u)/u)>1e-6)
            Warning("SECURE measure OBDM check failed!\n");
          if(fabs(u)>1e3)
            Warning("SECURE measure OBDM check passed, but strange value obtained: %" LE " !\n", u);
#endif

          if(W[w].status == ALIVE || W[w].status == KILLED) {
            TBDM.f[nR] += u;
            TBDM.N[nR] += 1.;

            // store TBDM as a function of two arguments BOSONS
            //n1 = (int) ((fabs(W[w].z[index_up]-W[w].z[index_down])+0.)/TBDM_MATRIX.step); // size of a pair
            if(fabs(W[w].z[index_up]-W[w].z[index_down])<Lhalfwf)
              n1 = (int) ((fabs(W[w].z[index_up]-W[w].z[index_down]))/TBDM_MATRIX.step+0.); // size of a pair
            else
              n1 = (int) (((L-fabs(W[w].z[index_up]-W[w].z[index_down])))/TBDM_MATRIX.step+0.); // size of a pair
            n2 = (int) (Sqrt(CM2)/TBDM_MATRIX.step+0.);
            if(n1>=0 && n1<TBDM_MATRIX.size && n2>=0 && n2<TBDM_MATRIX.size) {
              TBDM_MATRIX.f[n1][n2] += u;
              TBDM_MATRIX.N[n1][n2] ++;
            }
#ifdef OBDM_FERMIONS
            s = 1;
            for(i=0; i<N; i++) {
              if(i!=index_up && i!=index_down) {
                if(W[w].z[index_up]<W[w].z[i]) s *= -1;
                if(W[w].z[index_down]<W[w].z[i]) s *= -1;
                //if(W[w].z[index_up]+CMshift[2]<W[w].z[i]) s *= -1;
                //if(W[w].z[index_down]+CMshift[2]<W[w].z[i]) s *= -1;
                if(Wp[w].z[index_up]<W[w].z[i]) s *= -1;
                if(Wp[w].z[index_down]<W[w].z[i]) s *= -1;
              }
            }
            TBDMfermi.f[nR] += s*u;
            // store TBDM as a function of two arguments FERMIONS
            if(n1>=0 && n1<TBDM_MATRIX.size && n2>=0 && n2<TBDM_MATRIX.size) {
              TBDMfermi_MATRIX.f[n1][n2] += u*s;
              TBDMfermi_MATRIX.N[n1][n2] ++;
            }
#endif
          }
        }
      }
    }
  }
  return 0;
}

/******************************** TBDM ***************************************/
// uses only U() call
int MeasureTBDMfromU(void) {
  int w;

#pragma omp parallel for
  for(w=0; w<Nwalkers; w++) MeasureTBDMfromUWalker(w);

  return 0;
}

/********************* Measure TBDM from U ***********************************/
int MeasureTBDMfromUWalker(int w) {
  int i,j,m,n;
  int sign_displacement;
  DOUBLE u_old, u_new, u;
  DOUBLE dr[3],r2,r;
#ifdef HARD_SPHERE_HYPERRADIUS
  int overlap;
#endif

  u_old = U(W[w]);
  for(i=0; i<N; i++) { // loop over particle which is moved
    for(m=0; m<McMillan_points; m++) { // loop over McMillan points
#ifdef HARD_SPHERE_HYPERRADIUS // begin overlap loop
      overlap = ON;
      while(overlap == ON) { // repeat until no overlap
#endif
        CopyWalkerCoord(&Wp[w], &W[w]);

        n = rand() % OBDM.size; // displacement integer
        r = n*OBDM.step; // displacement distance
        r /= sqrt(DIMENSION); // displacement amplitude

        // displace particles i and j
        //j = i + 1;
        j = i + rand() % (N-1);
        if(j==N) j = 0;
        CaseX(sign_displacement = 2*(rand() % 2)-1);
        CaseX(Wp[w].x[i] += (DOUBLE)sign_displacement * r);
        CaseX(Wp[w].x[j] += (DOUBLE)sign_displacement * r);
        CaseY(sign_displacement = 2*(rand() % 2)-1);
        CaseY(Wp[w].y[i] += (DOUBLE)sign_displacement * r);
        CaseY(Wp[w].y[j] += (DOUBLE)sign_displacement * r);
        CaseZ(sign_displacement = 2*(rand() % 2)-1);
        CaseZ(Wp[w].z[i] += (DOUBLE)sign_displacement * r);
        CaseZ(Wp[w].z[j] += (DOUBLE)sign_displacement * r);
        ReduceToTheBoxXYZ(&Wp[w].x[i], &Wp[w].y[i], &Wp[w].z[i]);
        ReduceToTheBoxXYZ(&Wp[w].x[j], &Wp[w].y[j], &Wp[w].z[j]);

#ifdef HARD_SPHERE_HYPERRADIUS
        overlap = CheckWalkerHyperradiusOverlapping(Wp[w]);
      } // end overlap loop
#endif

      u = exp(U(Wp[w]) - u_old);

      if(n<OBDM.size && W[w].status == ALIVE) {
        if(W[w].status == ALIVE || W[w].status == KILLED) {
          // one-body correlation function
          TBDM.N[n] += W[w].weight;
          TBDM.f[n] += fabs(u) * W[w].weight;
        }
      }
    }
  }
  return 0;
}

/********************* Measure TBDM Pure *************************************/
int MeasureTBDMpure(void) {
#ifdef MEASURE_PURE_TBDM
  int w,w_v,i_moved,j_moved;
  DOUBLE xmi, ymi, zmi, xmj, ymj, zmj; // McMillan point for [i] and [j] particles
  DOUBLE r_i_old[3], r_j_old[3]; // store old coordinates
  DOUBLE shift = 0.5; // in units of L
  DOUBLE u;
  int terms = 0;
#ifndef TRIAL_1D // in 2D and 3D
  DOUBLE theta;
#endif
#ifdef TRIAL_3D
  DOUBLE phi;
#endif
  int Nup, Ndn;
  DOUBLE Uold, Unew;

#ifndef SPINFULL
  Warning("  define SPINFULL for MeasureTBDMpure");
#endif 
#ifdef SPINFULL_TUNNELING
  Error("  MeasureTBDMpure not implemented for SPINFULL_TUNNELING");
#endif 

  Nwalkers_pure = Npop_virtual;
  Nup = N/2; // works only in balanced case
  Ndn = N/2;

  for(terms=0; terms<Nwalkers_pure; terms++) {
    w = (int)(Random()*(DOUBLE)Nwalkers); // choose random walker to move

    Uold = U(W[w]); // old walue of the w.f.

    xmi = ymi = zmi = 0;
    xmj = ymj = zmj = 0;

    i_moved = (int)(Random()*(DOUBLE)(Nup)); // number of particle which will be moved r[i] -> rm
    j_moved = Nup + (int)(Random()*(DOUBLE)(Ndn)); // number of particle which will be moved r[j] -> rm

    r_i_old[0] = W[w].x[i_moved]; // store old coordinates
    r_i_old[1] = W[w].y[i_moved];
    r_i_old[2] = W[w].z[i_moved];
    r_j_old[0] = W[w].x[j_moved];
    r_j_old[1] = W[w].y[j_moved];
    r_j_old[2] = W[w].z[j_moved];

#ifdef TRIAL_3D
    theta = Random()*2.*PI;
    phi = Random()*PI;
    xmi = W[w].x[i_moved] + shift*L*Cos(theta)*Sin(phi);
    ymi = W[w].y[i_moved] + shift*L*Sin(theta)*Sin(phi);
    zmi = W[w].z[i_moved] + shift*L*Cos(phi);
    xmj = W[w].x[j_moved] + shift*L*Cos(theta)*Sin(phi);
    ymj = W[w].y[j_moved] + shift*L*Sin(theta)*Sin(phi);
    zmj = W[w].z[j_moved] + shift*L*Cos(phi);
#endif

#ifdef TRIAL_2D // generate 2D displacement
    theta = Random()*2.*pi;
    xmi = W[w].x[i_moved] + shift*L*Cos(theta);
    ymi = W[w].y[i_moved] + shift*L*Sin(theta);
    zmi = W[w].z[i_moved];
    xmj = W[w].x[j_moved] + shift*L*Cos(theta);
    ymj = W[w].y[j_moved] + shift*L*Sin(theta);
    zmj = W[w].z[j_moved];
#endif

#ifdef TRIAL_1D // generate 1D displacement
    zmi = W[w].z[i_moved] + shift*L;
    zmj = W[w].z[j_moved] + shift*L;
#endif

    if(xmi>L) xmi -= L;
    if(xmi<0) xmi += L;
    if(ymi>L) ymi -= L;
    if(ymi<0) ymi += L;
    if(zmi>L) zmi -= L;
    if(zmi<0) zmi += L;
    if(xmj>L) xmj -= L;
    if(xmj<0) xmj += L;
    if(ymj>L) ymj -= L;
    if(ymj<0) ymj += L;
    if(zmj>L) zmj -= L;
    if(zmj<0) zmj += L;

    W[w].x[i_moved] = xmi;
    W[w].y[i_moved] = ymi;
    W[w].z[i_moved] = zmi;
    W[w].x[j_moved] = xmj;
    W[w].y[j_moved] = ymj;
    W[w].z[j_moved] = zmj;

    Unew = U(W[w]); // new walue of the w.f.
    u = exp(Unew-Uold);

    // create a virtual walker
    w_v = dead_walkers_storage[Ndead_walkers-1]; // take it from dead walkers storage
    if(Ndead_walkers < 1) {
      Warning("  MeasureTBDMpure: too many walkers\n");
      return 1;
    }
    CopyWalker(&W[w_v], &W[w]); // inherit data
    W[w_v].R_dmc_end_of_step[i_moved][0] = xmi; // quadratic algorithm restores the coordinates
    W[w_v].R_dmc_end_of_step[i_moved][1] = ymi;
    W[w_v].R_dmc_end_of_step[i_moved][2] = zmi;
    W[w_v].R_dmc_end_of_step[j_moved][0] = xmj;
    W[w_v].R_dmc_end_of_step[j_moved][1] = ymj;
    W[w_v].R_dmc_end_of_step[j_moved][2] = zmj;
    W[w_v].OBDMpureB = fabs(u);
    W[w_v].OBDMpureF = u;
    W[w_v].OBDM_position = 0;
    W[w_v].status = VIRTUAL;
    W[w_v].weight = 0.; // virtual walker stores logarithmic weight
    W[w_v].OBDMpure_r = shift*L;

    Ndead_walkers--;

    W[w].x[i_moved] = r_i_old[0]; // restore coordinates in the walkers
    W[w].y[i_moved] = r_i_old[1];
    W[w].z[i_moved] = r_i_old[2];
    W[w].x[j_moved] = r_j_old[0];
    W[w].y[j_moved] = r_j_old[1];
    W[w].z[j_moved] = r_j_old[2];
  }
#endif
  return 0;
}

/************************ Empty Distribution Arrays **************************/
void EmptyDistributionArrays(void) {
#ifndef MEMORY_CONTIGUOUS
  int i;
#endif

  PD.times_measured = 0;
  PD.r2 = 0.;
  PD.g3 = 0;
  ArrayEmpty1D(PD.N, i, PD.size);

  RD.times_measured = 0;
  RD.r2 = 0.;
  RD.r4 = 0.;
  ArrayEmpty1D(RD.N, i, RD.size);
  RD.r2recent = 0.;

  if(measure_RDz) {
    RDz.times_measured = 0;
    RDz.r2 = 0.;
    RDz.r4 = 0.;
    RDz.r2recent = 0.;
    ArrayEmpty1D(RDz.N, i, RD.size);
  }

  /*if(MC == DIFFUSION && measure_RadDistr) {
    RD_pure.times_measured = 0;
    RD_pure.r2 = 0.;
    ArrayEmpty1D(RD_pure.N, i, RD.size);

    if(measure_RDz) {
      RDz_pure.times_measured = 0;
      RDz_pure.r2 = 0.;
      ArrayEmpty1D(RDz_pure.N, i, RD.size);
    }
  }*/

  if(measure_OBDM) {
    ArrayEmpty1D(OBDM.f, i, OBDM.size);
    ArrayEmpty1D(OBDM.N, i, OBDM.size);
#ifdef BC_ABSENT
    ArrayEmpty1D(OBDMtrap.f, i, OBDM.size);
    ArrayEmpty1D(OBDMtrap.N, i, OBDM.size);
#endif
#ifdef OBDM_FERMIONS
    ArrayEmpty1D(OBDMfermi.f, i, OBDM.size);
#endif
#ifdef SPECKLES
    OBDM.CFN = 0;
    OBDM.CF = 0.;
#endif
  }

  if(measure_TBDM) {
    ArrayEmpty1D(TBDM.f, i, TBDM.size);
    ArrayEmpty1D(TBDM.N, i, TBDM.size);
#ifdef OBDM_FERMIONS
    ArrayEmpty1D(TBDMfermi.f, i, TBDM.size);
#endif
  }

  /*if(measure_PairDistr_pure) {
    PD_pure.times_measured = 0;
    PD_pure.r2 = 0.;
    ArrayEmpty1D(PD_pure.N, i, PD_pure.size);
  }*/

  OrderParameter.cos = 0;
  OrderParameter.sin = 0;
  OrderParameter.times_measured = 0;
}

/**************************** Static Structure Factor *************************/
void MeasureStaticStructureFactor(int w) {
  int i,n;
  DOUBLE cosF, sinF;
//#ifdef SPINFULL
  DOUBLE sinF1, cosF1, sinF2, cosF2;
//#endif
  DOUBLE Sk_contr;
  int Sk_position;
#ifdef TRIAL_1D
  DOUBLE k;
#else
  int nk;
  DOUBLE k1,k2,k3;
#endif

#pragma omp atomic
  Sk.times_measured++;

  if(MC == DIFFUSION && W[w].times_measured>grid_pure_block) { // do a pure measurement
#pragma omp atomic
    Sk_pure.times_measured++;
    Sk_position = W[w].Sk_position;
#ifdef SECURE
    if(Sk_position<0 || Sk_position>=grid_pure_block) Warning("Measure Sk: bad grid_pure_block index!\n");
#endif
    for(n=0; n<gridSk; n++) {
#pragma omp atomic
      Sk_pure.f[n] += W[w].Sk[Sk_position][n]; // * W[w].weight;
      W[w].Sk[Sk_position][n] = 0.;
//#ifndef TRIAL_1D // i.e. 2D and 3D
#pragma omp atomic
      Sk_pure.N[n] += W[w].Sk_N[Sk_position][n]; // * W[w].weight;
      W[w].Sk_N[Sk_position][n] = 0;
//#endif
    }
  }

#ifdef TRIAL_1D
#ifdef SPINFULL
  for(n=0; n<Sk.size; n++) {
    k = Sk.k[n];
    sinF1 = cosF1 = sinF2 = cosF2 = 0.;
    for(i=0; i<N; i++) {
      if(W[w].spin[i] == 0) {
        cosF1 += Cos(k*W[w].z[i]);
        sinF1 += Sin(k*W[w].z[i]);
      }
      else {
        cosF2 += Cos(k*W[w].z[i]);
        sinF2 += Sin(k*W[w].z[i]);
      }
}
#else
  for(n=0; n<Sk.size; n++) {
    k = Sk.k[n];
    sinF = cosF = 0;
    for(i=0; i<N; i++) {
      cosF += Cos(k*W[w].z[i]);
      sinF += Sin(k*W[w].z[i]);
    }
#pragma omp atomic
    Sk.cos[n] += cosF;
#pragma omp atomic
    Sk.sin[n] += sinF;
#endif
#else // 2D and 3D cases
  for(nk=0; nk<Sk.size; nk++) {
    sinF = cosF = 0.;
    cosF1 = cosF2 = sinF1 = sinF2 = 0.; //????

    k1 = Sk.kx[nk];
    k2 = Sk.ky[nk];
    k3 = Sk.kz[nk];

    for(i=0; i<N; i++) {
      cosF += Cos(k1*W[w].x[i] + k2*W[w].y[i] + k3*W[w].z[i]);
      sinF += Sin(k1*W[w].x[i] + k2*W[w].y[i] + k3*W[w].z[i]);
    }
    n = Sk.index[nk];
#pragma omp atomic
    Sk.N[n] += N;
#endif
#ifdef SPINFULL
    Sk_contr = cosF1 * cosF2 + sinF1 * sinF2;
#else
    Sk_contr = cosF*cosF + sinF*sinF;
#endif
#pragma omp atomic
    Sk.f[n] += Sk_contr;
    if(MC == DIFFUSION) { // accumulation of the pure estimator
      W[w].Sk[W[w].Sk_position][n] += Sk_contr;
#ifndef TRIAL_1D // 2D and 3D
      W[w].Sk_N[W[w].Sk_position][n] += N;
#endif
    }
  }

  if(MC == DIFFUSION) { // move current position
    W[w].Sk_position++;
    if(W[w].Sk_position == grid_pure_block) W[w].Sk_position = 0;
  }
}

/************************** Measure Superfluid Density ************************/
void MeasureSuperfluidDensity(void) {
  int w,i;
  DOUBLE CM[3], dr[3] = {0,0,0};
  static int SD_position = 0;
  static int wait = 1;
  int inew, iold;

  for(w=0; w<Nwalkers; w++) {
    // empty array with center of the mass coordinate
    CaseX(CM[0] = 0.);
    CaseY(CM[1] = 0.);
    CaseZ(CM[2] = 0.);

    for(i=0; i<N; i++) { // find position of the CM
      CaseX(CM[0] += W[w].rreal[i][0]);
      CaseY(CM[1] += W[w].rreal[i][1]);
      CaseZ(CM[2] += W[w].rreal[i][2]);
    }
    CaseX(CM[0] /= (DOUBLE)N);
    CaseY(CM[1] /= (DOUBLE)N);
    CaseZ(CM[2] /= (DOUBLE)N);

    CaseX(W[w].CM[0][SD_position] = CM[0]);
    CaseY(W[w].CM[1][SD_position] = CM[1]);
    CaseZ(W[w].CM[2][SD_position] = CM[2]);

    if(!wait) {
      inew = SD_position;
      // accumulate diffusion
      for(i=0; i<SD.size-1; i++) {
        iold = inew - i - 1;
        if(iold<0) iold += SD.size;
        CaseX(dr[0] = W[w].CM[0][inew] - W[w].CM[0][iold]);
        CaseY(dr[1] = W[w].CM[1][inew] - W[w].CM[1][iold]);
        CaseZ(dr[2] = W[w].CM[2][inew] - W[w].CM[2][iold]);
        CaseXSD(SD.CM2[i] += dr[0]*dr[0]);
        CaseYSD(SD.CM2[i] += dr[1]*dr[1]);
        CaseZSD(SD.CM2[i] += dr[2]*dr[2]);
        //SD.N[i]++;
      }
      CaseXSD(SD.times_measured++);
      CaseYSD(SD.times_measured++);
      CaseZSD(SD.times_measured++);
    }
  }

  SD_position++;
  if(SD_position == SD.size) SD_position = 0;

  if(wait) {
    wait++;
    if(wait == SD.size)
      wait = OFF;
  }
}

/****************************** order parameter *****************************/
// Measures the order parameter: OP_tmp = (\psi*)\psi,
// where \psi = 1/N <\sum Exp(i 2\pi n r_i)>
void MeasureOrderParameter(void) {
#ifdef TRIAL_1D
  DOUBLE tmp,tmp1;
  int w, i;

  tmp1  = 2. * PI * n;
  for(w=0; w<Nwalkers; w++) {
    for(i=0; i<N; i++) {
      tmp = tmp1 * W[w].z[i];
      OrderParameter.cos += Cos(tmp);
      OrderParameter.sin += Sin(tmp);
    }
  }
  OrderParameter.times_measured += Nwalkers;
#endif
#ifdef TRIAL_2D // peak in a crystal
  DOUBLE kx,ky;
  int w, i;
  DOUBLE cosFx, sinFx;
  DOUBLE cosFy, sinFy;
  static int wait = 1;
  int position;

  kx  = 4. * PI / crystal_side_x;
  ky  = 4. * PI / crystal_side_y;
  for(w=0; w<Nwalkers; w++) {
    cosFx = sinFx = cosFy = sinFy = 0;
    for(i=0; i<N; i++) {
      cosFx += Cos(kx*W[w].x[i]);
      sinFx += Sin(kx*W[w].x[i]);
      cosFy += Cos(ky*W[w].y[i]);
      sinFy += Sin(ky*W[w].y[i]);
    }
    OrderParameter.cos += (cosFx*cosFx+sinFx*sinFx)/(DOUBLE) (N*N);
    OrderParameter.sin += (cosFy*cosFy+sinFy*sinFy)/(DOUBLE) (N*N);
    if(MC == DIFFUSION) {// pure estimator
      position = W[w].Sk_position;
      if(wait == OFF) { // save data and erase old data
        OrderParameter.cos_pure += W[w].Skx[position];
        OrderParameter.sin_pure += W[w].Sky[position];
        W[w].Skx[position] = 0;
        W[w].Sky[position] = 0;
        OrderParameter.times_measured_pure++;
      }
      W[w].Skx[position] += (cosFx*cosFx+sinFx*sinFx)/(DOUBLE) (N*N);
      W[w].Sky[position] += (cosFy*cosFy+sinFy*sinFy)/(DOUBLE) (N*N);
      W[w].Sk_position++; // move current position
      if(W[w].Sk_position == grid_pure_block ) W[w].Sk_position = 0;
    }
  }
  OrderParameter.times_measured += Nwalkers;

  if(MC == DIFFUSION && wait) { // skip the first block in PURE case
    wait++;
    if(wait == grid_pure_block+1) wait = 0;
  }
#endif
}

/****************************** Initialize Form Factor ************************/
void InitializeFormFactor(void) {
#ifdef TRIAL_1D
  DOUBLE dk_min;

#ifdef BC_ABSENT // trap
  dk_min = 0.;

  Message("  Form factor spacing: SKT_dk = %lf\n", SKT_dk);
  Message("  range covered %lf<k<%lf\n", SKT_dk, SKT_dk*(DOUBLE)(gridSKT_k));
#else // hom. system
  if(SKT_dk == 0) { // in units of pi*n
    SKT_dk = 2./(DOUBLE)gridSKT_k; // 0<k<2pi
  }

#ifdef VEXT_COS2 // lattice
  dk_min = 2.*PI/L; // dk
  dk_min /= (PI*n); // in units of Fermi momentum
#else // no lattice
  dk_min = 2./(DOUBLE) Ndens; // dk
  if(SKT_dk<dk_min) SKT_dk = 2./(DOUBLE)N; // adjust if k is too small
#endif

  if(fabs(SKT_dk/dk_min-(DOUBLE)((int)(SKT_dk/dk_min+1e-6)))>1e-6) {
    Warning("wrong SKT_dk = %lf, should be multiples of %lf\n", SKT_dk, dk_min);
  }
  Message("  Form factor spacing: SKT_dk = %lf [in units of pi/N]\n", SKT_dk);
  Message("  range covered %lf<k<%lf [kF]\n", SKT_dk, SKT_dk*gridSKT_k);
#endif

#else
  SKT_dk = 2.*PI/L;
#endif
}

/****************************** Laplace transform of Form factor **************/
// S(k,tau) = \int S(k,omega) Exp(-omega tau) d omega
void MeasureFormFactor(int w) {
  int i;
  DOUBLE rhoIm, rhoRe;
  int k;
  //static int index_of_current_position = -1;
  int index_displaced;
#ifndef TRIAL_1D //2D,3D
  DOUBLE k1,k2,k3;
#else //1D
  DOUBLE sign,kF,dk_min,k1;
  int j;

  kF = PI*n;
#endif

  if(MC != DIFFUSION) return;

  W[w].Sk_tau_position++;
  if(W[w].Sk_tau_position == gridSKT_t) W[w].Sk_tau_position = 0;

  for(k=0; k<gridSKT_k; k++) {
    rhoIm = rhoRe = 0.;
#ifdef BC_ABSENT
    DOUBLE r2accum = 0.;
#endif
    for(i=0; i<N; i++) {
#ifdef MEASURE_FORM_FACTOR_FOR_SUM_RULES // use A = x_i and A = x_i^2 operators S(tau) = <A(tau)A(0)>
#ifdef TRIAL_1D
      if(k==0) // x, can be positive and negative
        rhoRe += W[w].z[i] / (DOUBLE) N;
      else if(k%2 == 0) // k = 2, 4, i.e. other even powers k+1=1,3,5,... like x^3, sign should be taken into account
        rhoRe += pow(fabs(W[w].z[i]), (DOUBLE)(k+1)) * ((W[w].z[i]<0)?(-1.):(1.)) / (DOUBLE) N;
      else // k = 1,3, ... i.e. even powers k+1=2,4,... like x^2, no sign issues
        rhoRe += pow(fabs(W[w].z[i]), (DOUBLE)(k+1)) / (DOUBLE) N;
#endif
#ifdef TRIAL_2D
      DOUBLE xCM = 0;

      for(int j=0; j<N; j++) xCM += W[w].x[j];
      xCM /= (DOUBLE) N;

      if(k==0) // x
        rhoRe += W[w].x[i]-xCM;
      else if(k%2 == 0) // powers k+1=1,3,5,...
        rhoRe += pow(fabs(W[w].x[i]-xCM), (DOUBLE)(k+1)) * ((W[w].x[i]<0)?(-1.):(1.));
      else // powers k+1=2,4,...
        rhoRe += pow(fabs(W[w].x[i]-xCM), (DOUBLE)(k+1));
#endif
#else // use density operator S(k,tau) : <(exp(ikr_i(tau)) (exp(-ikr_j(0))>
#ifdef TRIAL_1D
#ifdef BC_ABSENT // trap
      if(SKT_kmin>0)
        k1 = (DOUBLE)(k)*SKT_dk + SKT_kmin;
      else // avoid k=0
        k1 = (DOUBLE)(k+1)*SKT_dk + SKT_kmin;
#else // hom. system
      k1 = (DOUBLE)(k+1)*SKT_dk*kF;
#endif
      rhoRe += Cos(k1*W[w].z[i]);
      rhoIm -= Sin(k1*W[w].z[i]);
#else // 2D and 3D cases
      k1 = Sk.kx[k+1]; // skip trivial [0,0,0] term
      k2 = Sk.ky[k+1];
      k3 = Sk.kz[k+1];
      rhoRe += Cos(k1*W[w].x[i] + k2*W[w].y[i] + k3*W[w].z[i]);
      rhoIm -= Sin(k1*W[w].x[i] + k2*W[w].y[i] + k3*W[w].z[i]);
      // will be normalized to k=0
      //rhoRe += Cos((DOUBLE)k*SKT_dk*W[w].x[i]);
      //rhoIm -= Sin((DOUBLE)k*SKT_dk*W[w].x[i]);
      //rhoRe += Cos((DOUBLE)k*SKT_dk*W[w].y[i]);
      //rhoIm -= Sin((DOUBLE)k*SKT_dk*W[w].y[i]);
#endif
#endif // MEASURE_FORM_FACTOR_FOR_SUM_RULES
    }

#ifndef MEASURE_FORM_FACTOR_FOR_SUM_RULES
#ifdef BC_ABSENT
    DOUBLE dr[3];
    int j;
    for(i=0; i<N; i++) {
      for(j=i+1; j<N; j++) {
        dr[0] = W[w].x[i] - W[w].x[j];
        dr[1] = W[w].y[i] - W[w].y[j];
        dr[2] = W[w].z[i] - W[w].z[j];
        r2accum += FindNearestImage(&dr[0], &dr[1], &dr[2]);
        rhoRe = pow(r2accum, (DOUBLE)k);
        rhoIm = 0.;
      }
    }
#endif
#endif

    W[w].rhoRe[W[w].Sk_tau_position][k] = rhoRe;
    W[w].rhoIm[W[w].Sk_tau_position][k] = rhoIm;
  }

    /* enable for calculation of Psi
    W[w].U = U((W[w]));

#ifdef TRIAL_1D // define sign
    sign = 1.;
    for(i=0; i<N; i++) {
      for(j=0; j<i; j++) {

        //dr = W[w].z[i]-W[w].z[j];
        //if(dr < Lhalf) dr += L;
        //if(dr > Lhalf) dr -= L;

        //if(dr<0) sign *= -1;
        sign *=  W[w].z[i]-W[w].z[j];
      }
    }
    if(sign>0)
      sign = 1;
    else
      sign = -1;
    W[w].psi[index_of_current_position][0] = exp(W[w].U)*sign;
#else
    W[w].psi[index_of_current_position][0] = exp(W[w].U);
#endif*/

  if(W[w].times_measured>=gridSKT_t) {
    int index_oldest = (W[w].Sk_tau_position + 1) % gridSKT_t;
    for(i=0; i<gridSKT_t; i++) {
      for(k=0; k<gridSKT_k; k++) {
        index_displaced = W[w].Sk_tau_position - i; // mixed
        if(index_displaced<0) index_displaced += gridSKT_t; // mixed
#pragma omp atomic // mixed
        SKT.f[i][k] += W[w].rhoRe[W[w].Sk_tau_position][k]*W[w].rhoRe[index_displaced][k] + W[w].rhoIm[W[w].Sk_tau_position][k]*W[w].rhoIm[index_displaced][k];

//        index_displaced = (index_oldest + i) % gridSKT_t; // pure
//#pragma omp atomic
//        SKT.f[i][k] += W[w].rhoRe[index_displaced][k]*W[w].rhoRe[index_oldest][k] + W[w].rhoIm[index_displaced][k]*W[w].rhoIm[index_oldest][k];

        //SKT.N[i][k]++;
      }
      //PhiTau.f[i][0] += W[w].psi[index_displaced][0]/W[w].psi[index_of_current_position][0];
    }
#pragma omp atomic
    SKT.times_measured++;
  }
}

/**************************** Lindemann ratio *********************************/
void MeasureLindemannRatio(void) {
#ifdef CRYSTAL
  int w,i,n;
  DOUBLE dx,dy,dz,r2;
  int recent_position, position;
  static int wait = 1;
  // Lindemann - Lozovik parameter
  int ns,_nns,z;
  double ds,ds0;
  static int initialized = OFF; // initialize neighbour matrix
  static int **mns; //matrix of neighbour sites [Ncryst][6]
  int nns,s;

  if(measure_Lindemann_Lozovik && !initialized) {

    mns = (int**) Calloc("mns ", N, sizeof(int*));
    for(i=0; i<Crystal.size; i++) mns[i] = (int*) Calloc("mns ", 6, sizeof(int));


    for(nns=0; nns<6; nns++) {
      for(s=0; s<Crystal.size; s++) {
        for(ds0=1E300,ns=0; ns<Crystal.size; ns++) {
          z=1;
          for(_nns=0;_nns<nns;_nns++) if(ns==mns[s][_nns]) z=0;
          if(s!=ns&&z) {
            ds=sqrt((Crystal.x[ns]-Crystal.x[s])*(Crystal.x[ns]-Crystal.x[s])+(Crystal.y[ns]-Crystal.y[s])*(Crystal.y[ns]-Crystal.y[s]));
            if(ds<ds0){ds0=ds;mns[s][nns]=ns;}
            ds=sqrt((Crystal.x[ns]-Crystal.x[s]+Lx)*(Crystal.x[ns]-Crystal.x[s]+Lx)+(Crystal.y[ns]-Crystal.y[s])*(Crystal.y[ns]-Crystal.y[s]));
            if(ds<ds0){ds0=ds;mns[s][nns]=ns;}
            ds=sqrt((Crystal.x[ns]-Crystal.x[s]-Lx)*(Crystal.x[ns]-Crystal.x[s]-Lx)+(Crystal.y[ns]-Crystal.y[s])*(Crystal.y[ns]-Crystal.y[s]));
            if(ds<ds0){ds0=ds;mns[s][nns]=ns;}
            ds=sqrt((Crystal.x[ns]-Crystal.x[s])*(Crystal.x[ns]-Crystal.x[s])+(Crystal.y[ns]-Crystal.y[s]+Ly)*(Crystal.y[ns]-Crystal.y[s]+Ly));
            if(ds<ds0){ds0=ds;mns[s][nns]=ns;}
            ds=sqrt((Crystal.x[ns]-Crystal.x[s]+Lx)*(Crystal.x[ns]-Crystal.x[s]+Lx)+(Crystal.y[ns]-Crystal.y[s]+Ly)*(Crystal.y[ns]-Crystal.y[s]+Ly));
            if(ds<ds0){ds0=ds;mns[s][nns]=ns;}
            ds=sqrt((Crystal.x[ns]-Crystal.x[s]-Lx)*(Crystal.x[ns]-Crystal.x[s]-Lx)+(Crystal.y[ns]-Crystal.y[s]+Ly)*(Crystal.y[ns]-Crystal.y[s]+Ly));
            if(ds<ds0){ds0=ds;mns[s][nns]=ns;}
            ds=sqrt((Crystal.x[ns]-Crystal.x[s])*(Crystal.x[ns]-Crystal.x[s])+(Crystal.y[ns]-Crystal.y[s]-Ly)*(Crystal.y[ns]-Crystal.y[s]-Ly));
            if(ds<ds0){ds0=ds;mns[s][nns]=ns;}
            ds=sqrt((Crystal.x[ns]-Crystal.x[s]+Lx)*(Crystal.x[ns]-Crystal.x[s]+Lx)+(Crystal.y[ns]-Crystal.y[s]-Ly)*(Crystal.y[ns]-Crystal.y[s]-Ly));
            if(ds<ds0){ds0=ds;mns[s][nns]=ns;}
            ds=sqrt((Crystal.x[ns]-Crystal.x[s]-Lx)*(Crystal.x[ns]-Crystal.x[s]-Lx)+(Crystal.y[ns]-Crystal.y[s]-Ly)*(Crystal.y[ns]-Crystal.y[s]-Ly));
            if(ds<ds0){ds0=ds;mns[s][nns]=ns;}
          }
        }
      }
    }
    initialized = ON;
  }

  if(wait == OFF) { // i.e. pure measurement can be done
    for(w=0; w<Nwalkers; w++) {
      recent_position = W[w].Lind_position;
      for(n=0; n<grid_pure_block; n++) {
        position = recent_position - n - 1;
        if(position < 0) position += grid_pure_block;
        LindemannRatio.Fpure[n] += W[w].LindF[position]; // * W[w].weight;
        LindemannRatio.Npure[n] += W[w].LindN[position]; // * W[w].weight;
      }
      W[w].LindF[recent_position] = 0.;
      W[w].LindN[recent_position] = 0;
    }
  }

  for(w=0; w<Nwalkers; w++) {
    for(i=0; i<N; i++) {
      if(measure_Lindemann_Lozovik == 0) {
        dx = W[w].x[i] - Crystal.x[i];
        dy = W[w].y[i] - Crystal.y[i];
        dz = W[w].z[i] - Crystal.z[i];
        r2 = FindNearestImage(&dx, &dy, &dz);
        LindemannRatio.F += r2;
        LindemannRatio.N++;
        if(MC == DIFFUSION) { // pure estimator
          W[w].LindF[W[w].Lind_position] += r2;
          W[w].LindN[W[w].Lind_position]++;
        }
      }
      else { // calculation of the Lindemann-Lozovik ratio
        for(nns=0; nns<6; nns++) {
          dx = W[w].x[i]-Crystal.x[i]-W[w].x[mns[i][nns]]+Crystal.x[mns[i][nns]];
          dy = W[w].y[i]-Crystal.y[i]-W[w].y[mns[i][nns]]+Crystal.y[mns[i][nns]];
          dz = 0.;
          r2 = FindNearestImage(&dx, &dy, &dz);
          LindemannRatio.F += r2;
          LindemannRatio.N++; // counts number of measurements times number of neighbours
          if(MC == DIFFUSION) { // pure estimator
            W[w].LindF[W[w].Lind_position] += r2;
            W[w].LindN[W[w].Lind_position]++;
          }
        }
      }
    }

    // move current position
    if(MC == DIFFUSION) {
      W[w].Lind_position++;
      if(W[w].Lind_position == grid_pure_block) W[w].Lind_position = 0;
    }
  }

  if(MC == DIFFUSION && wait) { // skip the first block in PURE case
    wait++;
    if(wait == grid_pure_block +1) wait = OFF;
  }
#else
  //Warning("  Cannot measure Lindemann ratio in absence of a crystal");
#endif
}

/****************** Measure Wave Function Projection ***************************/
void MeasureWaveFunctionProjection(void) {
  static int first_time = ON;
  static double k;
  double dr[3] = {0,0,0};
  double r2,r;
  double psi1, psi2;
  double ratio, norm;
  int i,j,w;

  if(first_time == ON) {
    Warning("  Measuring reweighting, instead of energy projection <psi_1|psi_2> is saved\n");
    k = PI / L;
  }

  ratio = norm = 0.;
  for(w=0; w<Nwalkers; w++) {
    //psi1 = psi2 = 1;
    psi1 = psi2 = 0;
    for(i=0; i<N; i++) {
      for(j=i+1; j<N; j++) {
        dr[0] = W[w].x[i] - W[w].x[j];
        dr[1] = W[w].y[i] - W[w].y[j];
        dr[2] = W[w].z[i] - W[w].z[j];
        r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
        r = Sqrt(r2);

        //psi1 *= cos(k*r); // Tonks-Girardeau
        //psi2 *= exp(-r/a); // McGuire
        //psi1 += log(cos(k*r)); // Tonks-Girardeau
        //psi2 -= r/a; // McGuire
        psi1 += log(r); // Tonks-Girardeau in trap
        if(j==1) 
          psi2 -= r/a;
        else
          psi2 += log(r);
      }
      psi1 -= 0.5*W[w].z[i]*W[w].z[i];
      psi2 -= 0.5*W[w].z[i]*W[w].z[i];
    }

#ifdef MEASURE_PROJECTION_PARTICLES_12
//#ifdef TRIAL_MCGUIRE
    ratio += exp(psi1-psi2);
    norm += exp(-2.*(psi2));
#else
    //Message(" psi1 = %lf, U = %lf", psi1,  U(W[w]));
    ratio += exp(psi2-psi1);
    norm += exp(-2.*psi1);
#endif
  }
  ratio /= (double) Nwalkers;
  norm /= (double) Nwalkers;

  SaveEnergyVMC(ratio, norm);
  //SaveEnergyCls(norm);

  first_time = OFF;
}

/****************************** hyper radius **********************************/
void MeasureHyperRadius(void) {
  int w,i,j,k,n;
  DOUBLE r, r2, dr[3];
  int HR_position;
  static int wait = 1;
  DOUBLE x1,y1,z1,x2,y2,z2,x3,y3,z3;
  DOUBLE m1,m2,m3,mur1,mur2,mur3,l1,l2,l3;
  DOUBLE m_up, m_dn,mu,M;
  DOUBLE HR21, HR22, HR23;
  DOUBLE mass_ratio = 1.;

  //HR.times_measured += Nwalkers;

  if(MC && wait == OFF) { // do a pure measurement
    HR_pure.times_measured += Nwalkers;
    for(w=0; w<Nwalkers; w++) {
      HR_position = W[w].HR_position;
      for(n=0; n<gridg3; n++) {
        HR_pure.N[n] += W[w].HR[HR_position][n]; // * W[w].weight;
        W[w].HR[HR_position][n] = 0;
      }
    }
  }

  m_dn = 0.5*(1.+mass_ratio);   // m_dn = M/(2mu)
  m_up = 0.5*(1.+1./mass_ratio);// m_up = m/(2mu)
  m1 = m_up;
  m2 = m_dn;
  m3 = m_dn;

  M = m1+m2+m3;
  mu = sqrt(m1*m2*m3/M);
  mur1= m2*m3/(m2+m3);
  mur2= m1*m3/(m1+m3);
  mur3= m1*m2/(m1+m2);
  l1 = Sqrt(mu/mur1);
  l2 = Sqrt(mu/mur2);
  l3 = Sqrt(mu/mur3);

  for(w=0; w<Nwalkers; w++) {
    for(i=0; i<N; i++) { // 1 - light 1
      x1 = W[w].x[i];
      y1 = W[w].y[i];
      z1 = W[w].z[i];
      for(j=i+1; j<N; j++) { // 2 - heavy 1
        x2 = W[w].x[j];
        y2 = W[w].y[j];
        z2 = W[w].z[j];
        for(k=j+1; k<N; k++) { // 3 - heavy 2
          x3 = W[w].x[k];
          y3 = W[w].y[k];
          z3 = W[w].z[k];

          dr[0] = x1-x2;
          dr[1] = y1-y2;
          dr[2] = z1-z2;
          FindNearestImage3D(r2, dr[0], dr[1], dr[2]);
          //r2 = dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2];
          HR23 = r2/(l3*l3);

          dr[0] = x2-x3;
          dr[1] = y2-y3;
          dr[2] = z2-z3;
          FindNearestImage3D(r2, dr[0], dr[1], dr[2]);
          //r2 = dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2];
          HR21 = r2/(l1*l1);

          dr[0] = x1-x3;
          dr[1] = y1-y3;
          dr[2] = z1-z3;
          FindNearestImage3D(r2, dr[0], dr[1], dr[2]);
          //r2 = dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2];
          HR22 = r2/(l2*l2);

          dr[0] = x1 - (m2*x2+m3*x3)/(m2+m3);
          dr[1] = y1 - (m2*y2+m3*y3)/(m2+m3);
          dr[2] = z1 - (m2*z2+m3*z3)/(m2+m3);
          FindNearestImage3D(r2, dr[0], dr[1], dr[2]);
          //r2 = dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2];
          HR21 += r2*l1*l1;

          dr[0] = x2 - (m1*x1+m3*x3)/(m1+m3);
          dr[1] = y2 - (m1*y1+m3*y3)/(m1+m3);
          dr[2] = z2 - (m1*z1+m3*z3)/(m1+m3);
          FindNearestImage3D(r2, dr[0], dr[1], dr[2]);
          //r2 = dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2];
          HR22 += r2*l2*l2;

          dr[0] = x3 - (m2*x2+m1*x1)/(m2+m1);
          dr[1] = y3 - (m2*y2+m1*y1)/(m2+m1);
          dr[2] = z3 - (m2*z2+m1*z1)/(m2+m1);
          FindNearestImage3D(r2, dr[0], dr[1], dr[2]);
          //r2 = dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2];
          HR23 += r2*l3*l3;

          r = sqrt(HR21);
          r = sqrt(HR22);
          r = sqrt(HR23);

          // convert to same definition as used by Petrov
          r *= sqrt(sqrt(3.)/2.);

          n = (int) (r / HR.step);
          //if(n>=HR.size) n = HR.size-1;
          if(n<HR.size) {
            //if(fabs(HR21/HR22-1.)>1e-3)
            //    Message("%lf\n%lf\n", sqrt(HR21), sqrt(HR22), sqrt(HR23));
            HR.N[n]++; // mixed estimator
            //HR.f[n] += 1./(r*r*r*r*r);
            HR.f[n] += 1./(r*r*r);
          }
          if(MC) W[w].HR[W[w].HR_position][n]++; // pure estimator
          HR.times_measured++;
        }
      }
    }
  }
  /*  if(MC) { // move current position
      W[w].HR_position++;
      if(W[w].HR_position == grid_pure_block) W[w].HR_position = 0;
    }
  }*/
  if(MC && wait) { // skip the first block in PURE case
    wait++;
    if(wait == grid_pure_block +1) wait = OFF;
  }
}

/******************* Measure Effective Potential ******************************/
void MeasureEffectivePotential(void) {
  int w,i,n;
  DOUBLE r, r2, dr[3] = {0,0,0};
  DOUBLE Rp[3]; // r' - position where the field is measured

  Veff.times_measured += Nwalkers;

  for(n=0; n<Veff.size; n++) {
    Rp[0] = Veff.step * ((DOUBLE)n+0.5);
    Rp[1] = 0;
    Rp[2] = 0;
    for(w=0; w<Nwalkers; w++) {
      for(i=0; i<N; i++) { // sum over spins
        dr[0] = W[w].x[i] - Rp[0];
        dr[1] = W[w].y[i] - Rp[1];
        dr[2] = W[w].z[i] - Rp[2];
        r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
        r = Sqrt(r2);
        Veff.N[n]++; // mixed estimator
        Veff.f[n] += InteractionEnergy_ij(r, W[w].spin[i], Nspin/2);
      }
    }
  }
}

/************************** Measure Pure Potential Energy *********************/
void MeasurePurePotentialEnergy(void) {
  int w,i;
  static int pure_position = 0;
  static int wait = 1;
  int inew, iold;

  for(w=0; w<Nwalkers; w++) {
    W[w].Epot_pure[pure_position] = W[w].Epot; // update the recent value

    if(!wait) {
      inew = pure_position;
      for(i=0; i<Epot_pure.size; i++) { // loop over time
        iold = inew - i; // position of a previous measurement
        if(iold<0) iold += Epot_pure.size; // check bounds

        Epot_pure.f[i] += W[w].Epot_pure[iold];
      }
      Epot_pure.times_measured++;
    }
  }

  pure_position++;
  if(pure_position == grid_pure_block) pure_position = 0;

  if(wait) {
    wait++;
    if(wait == Epot_pure.size) wait = OFF;
  }
}

/*********************************** Spin Polarization *************************/
// calculates <P^2>
DOUBLE SpinPolarization(void) {
  int i, w;
  DOUBLE polarization_walker, sigma2;

  sigma2 = 0.;
  for(w=0; w<Nwalkers; w++) {
    polarization_walker = 0.;
    for(i=0; i<N; i++) {
      polarization_walker += (DOUBLE)W[w].spin[i]; // spin [0,1], polarization_walker[0,N]
    }
    polarization_walker /= (DOUBLE) N; // polarization_walker [0,1]
    polarization_walker -= 0.5; // polarization_walker [-0.5,0.5]
    polarization_walker *= 2; // polarization_walker [-1,1]
    sigma2 += polarization_walker*polarization_walker;
  }
  return sigma2 / (DOUBLE) Nwalkers;
}

/******************************* Measure Hessian Matrix ************************/
void MeasureHessianMatrix(void) {
  int i,j, N_long_vector, index;
  DOUBLE V_plus_dxi;
  DOUBLE V_minus_dxi;
  DOUBLE V_plus_dxi_plus_dxj;
  DOUBLE V_plus_dxi_minus_dxj;
  DOUBLE V_minus_dxi_plus_dxj;
  DOUBLE V_minus_dxi_minus_dxj;
  DOUBLE xi,yi,zi,xj,yj,zj;
  DOUBLE dx;
  DOUBLE *r, **FirstDerivative, **Hessian;

  N_long_vector = 0;
  CaseX(N_long_vector += N);
  CaseY(N_long_vector += N);
  CaseZ(N_long_vector += N);

  Message("  Saving the Hessian matrix\n");
  r = (DOUBLE*) Calloc("vecor r", N_long_vector, sizeof(DOUBLE));
  ArrayCalloc2D(Hessian, "Hessian ", i, N_long_vector, N_long_vector, DOUBLE, "DOUBLE");
  ArrayCalloc2D(FirstDerivative, "FirstDerivative ", i, N_long_vector, 2, DOUBLE, "DOUBLE");

//DOUBLE WalkerPotentialEnergy(struct Walker *walker)
  index = 0; // copy Walker coordinates to the long vector r[]
  CaseX(for(i=0;i<N;i++) r[index++] = W[0].x[i]);
  CaseY(for(i=0;i<N;i++) r[index++] = W[0].y[i]);
  CaseZ(for(i=0;i<N;i++) r[index++] = W[0].z[i]);

  dx = 1e-3*L; // dx ~ 1e-3

  // calculate vector with first derivatives
  for(i=0; i<N_long_vector; i++) {
    CopyRawVectorToWalkerDisplace_i_j(r, &W[0], i, +dx, 0, 0);
    V_plus_dxi = WalkerPotentialEnergy(&W[0]);
    CopyRawVectorToWalkerDisplace_i_j(r, &W[0], i, -dx, 0, 0);
    V_minus_dxi = WalkerPotentialEnergy(&W[0]);

    FirstDerivative[i][0] = r[i]; // vector with coordinates
    FirstDerivative[i][1] = (V_plus_dxi-V_minus_dxi) / (2.*dx); // vector with the first derivative
  }
  SaveMatrix("FirstDerivative.dat", FirstDerivative, N_long_vector, 2);

  // calculate Hessian matrix
  for(i=0; i<N_long_vector; i++) {
    for(j=0; j<N_long_vector; j++) {
      CopyRawVectorToWalkerDisplace_i_j(r, &W[0], i, +dx, j, +dx);
      V_plus_dxi_plus_dxj = WalkerPotentialEnergy(&W[0]);
      CopyRawVectorToWalkerDisplace_i_j(r, &W[0], i, +dx, j, -dx);
      V_plus_dxi_minus_dxj = WalkerPotentialEnergy(&W[0]);
      CopyRawVectorToWalkerDisplace_i_j(r, &W[0], i, -dx, j, +dx);
      V_minus_dxi_plus_dxj = WalkerPotentialEnergy(&W[0]);
      CopyRawVectorToWalkerDisplace_i_j(r, &W[0], i, -dx, j, -dx);
      V_minus_dxi_minus_dxj = WalkerPotentialEnergy(&W[0]);
      Hessian[i][j] = (V_plus_dxi_plus_dxj-V_plus_dxi_minus_dxj-V_minus_dxi_plus_dxj+V_minus_dxi_minus_dxj) / (4.*dx*dx);
    }
  }
  SaveMatrix("Hessian.dat", Hessian, N_long_vector, N_long_vector);
}

/***************************** Kinetic Energy Numerically ********************/
// F = (\nabla \Psi) / \Psi
// f"/f = u" + (u')^2
DOUBLE MeasureKineticEnergyNumerically(struct Walker *walker) {
  int i, j;

  DOUBLE epsilon, epsilon2; // accuracy
  DOUBLE coord_store; // store coordinates
  DOUBLE Uo, Umax, Umin; // moving by epsilon
  DOUBLE Up; // U' = (dU/dx)
  DOUBLE Upp; // U" = (d2U/dx2)
  DOUBLE Ekin = 0;

  DOUBLE Rmax = 0.;

  for(i=0; i<N; i++) { // find largest position
    CaseX(if(fabs(W->x[i]) > Rmax) Rmax = fabs(W->x[i]));
    CaseY(if(fabs(W->y[i]) > Rmax) Rmax = fabs(W->y[i]));
    CaseZ(if(fabs(W->z[i]) > Rmax) Rmax = fabs(W->z[i]));
  }
  epsilon = Rmax * 1e-5; // goal relative error of 1e-10
  epsilon2 = epsilon*epsilon;

  Uo = U(*walker);

  for(i=0; i<N; i++) {
#ifdef MOVE_IN_X
    coord_store = walker->x[i];
    walker->x[i] += epsilon;
    Umax = U(*walker);
    walker->x[i] = coord_store - epsilon;
    Umin = U(*walker);
    walker->x[i] = coord_store;
    Up = (Umax - Umin) / (2.*epsilon);
    Upp = (Umax - 2.*Uo + Umin) / epsilon2;
    Ekin -= Upp + Up*Up; 
#endif
#ifdef MOVE_IN_Y
    coord_store = walker->y[i];
    walker->y[i] += epsilon;
    Umax = U(*walker);
    walker->y[i] = coord_store - epsilon;
    Umin = U(*walker);
    walker->y[i] = coord_store;
    Up = (Umax - Umin) / (2.*epsilon);
    Upp = (Umax - 2.*Uo + Umin) / epsilon2;
    Ekin -= Upp + Up*Up; 
#endif
#ifdef MOVE_IN_Z
    coord_store = walker->z[i];
    walker->z[i] += epsilon;
    Umax = U(*walker);
    walker->z[i] = coord_store - epsilon;
    Umin = U(*walker);
    walker->z[i] = coord_store;
    Up = (Umax - Umin) / (2.*epsilon);
    Upp = (Umax - 2.*Uo + Umin) / epsilon2;
    Ekin -= Upp + Up*Up; 
#endif
  }

  return 0.5*Ekin;
}