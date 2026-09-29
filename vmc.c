/*vmc.c*/

#include <string.h>
#include <float.h>
#include "vmc.h"
#include "main.h"
#include "randnorm.h"
#include "utils.h"
#include "trial.h"
#include "dmc.h"
#include "crystal.h"
#include "spline.h"
#include "speckles.h"
#include "compatab.h"
#include "memory.h"
#include MATHINCLUDE
#include "ewald.h"
#include "mymath.h"

#define INFTY 1e100

#ifdef HARD_SPHERE
/********************************** Overlapping *******************************/
int Overlapping(DOUBLE **R) {
  int i,j;
  DOUBLE dr[3] = {0,0,0};
  DOUBLE r2;

  for(i=0; i<N; i++) {
    for(j=i+1; j<N; j++) {
      dr[0] = R[i][0] - R[j][0];
      dr[1] = R[i][1] - R[j][1];
      dr[2] = R[i][2] - R[j][2];
      r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
      if(r2 <= a2) 
        return 1;
    }
  }

  return 0;
}

/********************************** Overlapping Walker ************************/
int OverlappingWalker(const struct Walker* W) {
  int i,j;
  DOUBLE dr[3] = {0,0,0};
  DOUBLE r2;

  for(i=0; i<N; i++) {
    for(j=i+1; j<N; j++) {
      dr[0] = W->x[i] - W->x[j];
      dr[1] = W->y[i] - W->y[j];
      dr[2] = W->z[i] - W->z[j];
      r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
      if(r2 <= a2) return 1;
    }
  }

  return 0;
}

/************************ Check Overlapping ********************************/
int CheckOverlapping(void) {

  int w,i,j;
  DOUBLE dr[3] = {0,0,0};
  DOUBLE r2;

  for(w=0; w<Nwalkers; w++) {
    for(i=0; i<N; i++) {
      for(j=i+1; j<N; j++) {
        dr[0] = W[w].x[i] - W[w].x[j];
        dr[1] = W[w].y[i] - W[w].y[j];
        dr[2] = W[w].z[i] - W[w].z[j];

        r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);

        if(r2 <= a2)
          return w+1;
      }
      if(W[w].weight<0.5 && W[w].weight>1.5) {
        Error("wrong weight W[%i] -> " "%" LF, w, W[w].weight);
      }
    }
  }

  return 0;
}

/************************ Check Walker Overlapping ***************************/
int CheckWalkerOverlapping(const struct Walker W) {
  int i,j;
  DOUBLE dr[3] = {0,0,0};
  DOUBLE r2;

  for(i=0; i<N; i++) {
    for(j=i+1; j<N; j++) {
      dr[0] = W.x[i] - W.x[j];
      dr[1] = W.y[i] - W.y[j];
      dr[2] = W.z[i] - W.z[j];

      r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);

      if(r2 < a2)
        return 1;
    }
  }

  return 0;
}
#endif // end HS

/***************************** Move One By One *******************************/
void VMCMoveOneByOne(int w) {
  int i,k;
  DOUBLE du;
  DOUBLE xi;
  DOUBLE Wp[3] = {0,0,0};
  DOUBLE dr[3] = {0,0,0};
  DOUBLE dx,dy,dz;
  DOUBLE x = 0.;
  DOUBLE y = 0.;
  DOUBLE z = 0.;
  DOUBLE r2;
#ifdef HARD_SPHERE
  int overlap;
#endif
#ifdef SPINFULL
  int spin_new;
#endif
  int move_accepted;
#ifdef TRIAL_3D
  DOUBLE scale;
  scale = 1./Sqrt(lambda);
#endif

#ifdef SCALABLE
  int cc;
  int j;

  UpdateScalableTable(&W[w]);
#endif

#ifdef CHAINS_PINNED_TO_CENTER_OF_MASS
  int j;
  int Nchain;
#endif 

#ifdef CENTER_OF_MASS_IS_NOT_MOVED
  Error("CENTER_OF_MASS_IS_NOT_MOVED is not implemented for VMCMoveOneByOne\n");
#endif

  // Move particle i of the walker w
  for(i=0; i<N; i++) {
    du = 0.;

#ifdef FIRST_PARTICLE_IS_NOT_MOVED
    CaseX(W[w].x[0] = 0.);
    CaseY(W[w].y[0] = 0.);
    CaseZ(W[w].z[0] = 0.);
    if(i==0) i++; // skip the first particle
#endif

#ifdef HARD_SPHERE
    overlap = OFF;
#endif

    CaseX(x = W[w].x[i]);
    CaseY(y = W[w].y[i]);
    CaseZ(z = W[w].z[i]);

#ifdef SPINFULL
#ifdef SPINFULL_TUNNELING
    spin_new = rand()%2; // allow spin flip
#else
    spin_new = W[w].spin[i]; // do not change the spin
#endif
#endif

    // Wp is the trial position
    RandomNormal3(&dx, &dy, &dz, 0., Sqrt(dt_vmc));

#ifdef SCALABLE
    if(dx>Lhalfwf) {
      Warning("  Move One By One: decrease time step, particle jumps out of the cell!\n");
      Warning("  Decreasing the time step %" LF " -> %" LF "\n", dt_vmc, dt_vmc*0.1);
      dt_vmc *= 0.1;
      Warning("  Increasing acceptance rate %" LF " -> %" LF "\n", acceptance_rate, 0.5*(acceptance_rate+1));
    }
#endif

#ifdef SCALE_VMC_AMPLITUDE_WITH_ALPHA
    if(alpha_Rx>1e-15) dx /= sqrt(alpha_Rx);
    if(alpha_Ry>1e-15) dy /= sqrt(alpha_Ry);
    if(alpha_Rz>1e-15) dz /= sqrt(alpha_Rz);
#endif

    CaseX(Wp[0] = x + dx);
    CaseY(Wp[1] = y + dy);
    CaseZ(Wp[2] = z + dz);

    ReduceToTheBox(Wp);

//#ifdef ONE_BODY_TRIAL_TERMS // One body terms
//    du += OneBodyU(Wp[0], Wp[1], Wp[2], i) - OneBodyU(x, y, z, i);
//#endif
#ifdef ONE_BODY_TRIAL_TERMS // One body terms
#ifdef CHAINS_PINNED_TO_CENTER_OF_MASS
    // calculate old  one body terms of one chain
    Nchain = i / Nspin;
    for (j = Nspin * Nchain; j < Nspin * (Nchain + 1); j++) {
      du -= OneBodyU(W[w].x[j], W[w].y[j], W[w].z[j], j, W[w].w);
    }
    //calculate new chain center of mass
    CalculateOneChainCenterOfMass(Wp[0], Wp[1], Wp[2], i, W[w].w);

    // calculate new  one body terms of one chain
    for (j = Nspin * Nchain; j < Nspin * (Nchain + 1); j++) {
      if (j != i) {
        du += OneBodyU(W[w].x[j], W[w].y[j], W[w].z[j], j, W[w].w);
      }
      else {
        du += OneBodyU(Wp[0], Wp[1], Wp[2], j, W[w].w);
      }
    }
#else
    du += OneBodyU(Wp[0], Wp[1], Wp[2], i, W[w].w) - OneBodyU(x, y, z, i, W[w].w);
#endif
#endif

#ifdef MANY_BODY_TRIAL_TRAP_CENTER_OF_MASS_NODE_Z
    DOUBLE z_CM = 0.; // old
    
    for(k=0; k<N; k++) z_CM += W[w].z[k];
    du -= log(fabs(z_CM)); // (z1 + ... +zN)
    z_CM += Wp[2] - W[w].z[i]; // new
    du += log(fabs(z_CM)); // (z1 + ... +zN)
#endif

#ifdef MANY_BODY_TRIAL_TRAP_CENTER_OF_MASS_HERMITE2_Z
    DOUBLE z2 = 0.; // old
    
    for(k=0; k<N; k++) z2 += (W[w].z[k] / (DOUBLE) N)*(W[w].z[k] / (DOUBLE) N);
    du -= log(fabs(z2 - beta));
    z2 += (Wp[2] / (DOUBLE) N)*(Wp[2] / (DOUBLE) N) - (W[w].z[i] / (DOUBLE) N)*(W[w].z[i] / (DOUBLE) N); // new
    du += log(fabs(z2 - beta));
#endif

#ifdef CRYSTAL // crystal
#ifdef CRYSTAL_SYMMETRIC
    Error("  Cannot use VMCMoveOneByOne for symmetrized crystal w.f.\n");
#endif
#endif

#ifdef SCALABLE
    cc = FindCell(x, y, z); // current cell
    for(j=0; j<W[w].c_Nall[cc]; j++) { // loop over particles in the nearest cells
      k = W[w].c_index_all[cc][j]; // index of a neighbour particle
      if(k != i) {
        // calculate new value of the wavefunction
        dr[0] = Wp[0] - W[w].x[k];
        dr[1] = Wp[1] - W[w].y[k];
        dr[2] = Wp[2] - W[w].z[k];
        r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
        if(CheckInteractionConditionWF(dr[0], dr[1], dr[2], r2)) du += InterpolateU(&G, Sqrt(r2), W[w].spin[i], W[w].spin[k]);

        // calculate old value of the wavefunction
        dr[0] = x - W[w].x[k];
        dr[1] = y - W[w].y[k];
        dr[2] = z - W[w].z[k];
        r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
        if(CheckInteractionConditionWF(dr[0], dr[1], dr[2], r2)) du -= InterpolateU(&G, Sqrt(r2), W[w].spin[i], W[w].spin[k]);
      }
    }
#else
    for(k=0; k<N; k++) {
      if(k != i) {
        // calculate new value of the wavefunction
        CaseX(dr[0] = Wp[0] - W[w].x[k]);
        CaseY(dr[1] = Wp[1] - W[w].y[k]);
        CaseZ(dr[2] = Wp[2] - W[w].z[k]);
        r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);

#ifdef HARD_SPHERE
        if(r2<a2) {
          overlap = ON;
        }
        else {
#endif

#ifndef MEASURE_PROJECTION_PARTICLES_12
        if(CheckInteractionConditionWF(dr[0], dr[1], dr[2], r2)) 
#ifdef SPINFULL
          du += InterpolateU(&G, Sqrt(r2), spin_new, W[w].spin[k]);
#else
          du += InterpolateU(&G, Sqrt(r2), W[w].spin[i], W[w].spin[k]);
#endif
#else
        if(i==0 && k==1) 
          du -= Sqrt(r2)/a;
        else
          du += InterpolateU(&G, Sqrt(r2), W[w].spin[i], W[w].spin[k]);
#endif

#ifdef HARD_SPHERE
        }
#endif

        // calculate old value of the wavefunction
        CaseX(dr[0] = x - W[w].x[k]);
        CaseY(dr[1] = y - W[w].y[k]);
        CaseZ(dr[2] = z - W[w].z[k]);
        r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
#ifndef MEASURE_PROJECTION_PARTICLES_12
        if(CheckInteractionConditionWF(dr[0], dr[1], dr[2], r2)) du -= InterpolateU(&G, Sqrt(r2), W[w].spin[i], W[w].spin[k]);
#else
        if(i==0 && k==1) 
          du += Sqrt(r2)/a;
        else
          du -= InterpolateU(&G, Sqrt(r2), i, k);
#endif
      }
    }
#endif

    // The Metropolis code
    // f = Exp(u), p = Exp(up-u)
    move_accepted = OFF;
#ifdef HARD_SPHERE
    if(overlap) {
      rejected++;
    }
    else {
#endif
    if(du > 0) {
      move_accepted = ON;
    }
    else {
      xi = Random();
      if(xi<Exp(2.*du)) {
        move_accepted = ON;
      }
      else {
        rejected++;
      }
    }
#ifdef HARD_SPHERE // Metropolis end
    }
#endif

    if(move_accepted) {
      CaseX(W[w].x[i] = Wp[0]);
      CaseY(W[w].y[i] = Wp[1]);
      CaseZ(W[w].z[i] = Wp[2]);
      W[w].U -= du;
      accepted++;
#ifdef SPINFULL_TUNNELING
      W[w].spin[i] = spin_new;
#endif

#ifdef SECURE
       if(CheckWalkerOverlapping(W[w])) Error("MoveOneByOne error");
#endif
    }
    else {
#ifdef CHAINS_PINNED_TO_CENTER_OF_MASS
      CalculateOneChainCenterOfMass(W[w].x[i], W[w].y[i], W[w].z[i], i, w);
#endif
  }
}
}

/***************************** Move All **************************************/
void VMCMoveAll(int w) {
  int i;
  DOUBLE xi;
  DOUBLE dx,dy,dz;
  DOUBLE x,y,z;
  DOUBLE sigma;
  int overlap;
  DOUBLE scale;
  int accept;

#ifdef SPINFULL_TUNNELING
  Error("spin flip not implemented for VMCMoveAll, use VMCMoveOne instead\n");
#endif

  scale = 1./Sqrt(lambda);
  overlap = OFF;
  sigma = Sqrt(dt_vmc);

  // Move all particles of the walker w
  for(i=0; i<N; i++) {
    x = W[w].x[i];
    y = W[w].y[i];
    z = W[w].z[i];

    RandomNormal3(&dx, &dy, &dz, 0., sigma);

#ifdef HARD_SPHERE_HYPERRADIUS
    tried++;
#endif

#ifdef SECURE
#ifndef BC_ABSENT
#ifdef MOVE_IN_X
    if(dx*dx>L2) {
      Warning("Gaussian jump: large time step (x)! Exiting to avoid possible problems! \n");
      return;
    }
#endif
#ifdef MOVE_IN_Y
    if(dy*dy>L2) {
      Warning("Gaussian jump: large time step (y)! Exiting to avoid possible problems! \n");
      return;
    }
#endif
#ifdef MOVE_IN_Z
    if(dz*dz>L2) {
      Warning("Gaussian jump: large time step (z)! Exiting to avoid possible problems! \n");
      return;
    }
#endif
#endif
#endif

    // wp is the trial position
    CaseX(Wp[w].x[i] = x + dx);
    CaseY(Wp[w].y[i] = y + dy);
    CaseZ(Wp[w].z[i] = z + dz);
#ifdef SPINFULL
    Wp[w].spin[i] = W[w].spin[i];
#endif
  }

#ifdef CENTER_OF_MASS_IS_NOT_MOVED // propose such movement that CM position is not moved
  AdjustCenterOfMassWalker(&Wp[w]);
#endif
#ifdef CENTER_OF_MASS_Z_IS_NOT_MOVED // propose such movement that CM position is not moved
  AdjustCenterOfMassWalker(&Wp[w]);
#endif
#ifdef CENTER_OF_MASS_IS_NOT_MOVED_TWO_PARTS
  AdjustCenterOfMassWalker(&Wp[w]);
#endif
#ifdef FIRST_PARTICLE_IS_NOT_MOVED // actually particles are moved such that the first particle is in the origin
  AdjustCenterOfMassWalker(&Wp[w]);
#endif

  ReduceWalkerToTheBox(&Wp[w]);

#ifdef SCALABLE
  UpdateScalableTable(&Wp[w]);
#endif
  Wp[w].U = U(Wp[w]);
  Wp[w].weight = 1;

  //-infinity is returned if particles overlap in new configuration
  if(Wp[w].U == -INFTY) {
    overlaped++;
    rejected++;
    return;
  }

#ifdef SECURE
  if(CheckWalkerOverlapping(Wp[w]))
    Error("MoveAll error");
#endif

  if(overlap == ON) {
    overlaped++;
    rejected++;
  }
  else {// The Metropolis code
    if(Wp[w].U > W[w].U) {// f^2 = Exp(2u), p = Exp(2(up-u))
      accept = ON;
    }
    else {
      xi = Random();
      if(xi<Exp(2.*(Wp[w].U-W[w].U))) {
        accept = ON;
      }
      else {
        rejected++;
        return;
      }
    }
  }

#ifdef HARD_SPHERE_HYPERRADIUS
  if(CheckWalkerHyperradiusOverlapping(Wp[w])) {
    overlaped++;
    accept = OFF;
  }
#endif

  if(accept) {
    accepted++;
    CopyWalker(&W[w], &Wp[w]);
    return;
  }
}

/***************************** Smart MC move *********************************/
//   R' = R + dt F(R) + dt chi
//   p_chi(dr,dt) ~ Exp(-0.5 dr^2/dt)
//   p = min [1, Exp(w)]
//   w = 2(U(R')-U(R)) - 1/2 [F(R')+F(R)] [(F(R')-F(R))dt+2(R'-R)]
//
//   equivalently (with R' not reduced to the box)
//   w = 2[U(R')-U(R)] + 1/(2dt) [(R'-R-F(R)dt)^2-(R-R'-F(R')dt)^2]
void VMCMoveDrift(int w) {
  int i;
  DOUBLE xi;
  DOUBLE dx, dy, dz;
  DOUBLE dexp = 0.;
  DOUBLE dt_move;

  dt_move = dt_vmc;

#ifdef SCALABLE
  UpdateScalableTable(&W[w]);
#endif
  //W[w].U = U(W[w]);
  CopyWalkerToVector(W[w].R, W[w]);
  DriftForce(W[w].F, W[w].R, w);

  for(i=0; i<N; i++) {
    RandomNormal3(&dx, &dy, &dz, 0., Sqrt(dt_move)); //(x,y,z, mu, sigma)
    // wp is the trial position

    /*if(measure_SD) {
      CaseX(rreal[i][0] += dx);
      CaseY(rreal[i][1] += dy);
      CaseZ(rreal[i][2] += dz);
    }*/

    CaseX(W[w].dR[i][0] = dt_move * W[w].F[i][0] + dx);
    CaseY(W[w].dR[i][1] = dt_move * W[w].F[i][1] + dy);
    CaseZ(W[w].dR[i][2] = dt_move * W[w].F[i][2] + dz);

    CaseX(W[w].Rp[i][0] = W[w].x[i] + W[w].dR[i][0]);
    CaseY(W[w].Rp[i][1] = W[w].y[i] + W[w].dR[i][1]);
    CaseZ(W[w].Rp[i][2] = W[w].z[i] + W[w].dR[i][2]);

    //ReduceToTheBox(W[w].Rp[i]);
  }

  CopyVectorToWalker(&Wp[w], W[w].Rp, W[w].spin);
  ReduceWalkerToTheBox(&Wp[w]);
  Wp[w].U = U(Wp[w]);
  DriftForce(W[w].Fp, W[w].Rp, w);

  // w = 2(U(R')-U(R)) - 1/2 [F(R')+F(R)] [(F(R')-F(R))dt+2(R'-R)]
  for(i=0; i<N; i++) {
    CaseX(dexp += (W[w].Fp[i][0]+W[w].F[i][0])*((W[w].Fp[i][0]-W[w].F[i][0])*dt_move+2.*W[w].dR[i][0]));
    CaseY(dexp += (W[w].Fp[i][1]+W[w].F[i][1])*((W[w].Fp[i][1]-W[w].F[i][1])*dt_move+2.*W[w].dR[i][1]));
    CaseZ(dexp += (W[w].Fp[i][2]+W[w].F[i][2])*((W[w].Fp[i][2]-W[w].F[i][2])*dt_move+2.*W[w].dR[i][2]));
  }
  dexp = -0.5*dexp + 2.*(Wp[w].U - W[w].U); // with detailed balance

  /*// w = 2[U(R')-U(R)] + 1/(2dt) [(R'-R-F(R)dt)^2-(R-R'-F(R')dt)^2]
  // NB: can be used only if R' was not reduced to the box!
  for(i=0; i<N; i++) {
    CaseZ(dexp += ( (W[w].Rp[i][2]-W[w].z[i]-W[w].F[i][2]*dt)*(W[w].Rp[i][2]-W[w].z[i]-W[w].F[i][2]*dt)
                   -(W[w].z[i]-W[w].Rp[i][2]-W[w].Fp[i][2]*dt)*(W[w].z[i]-W[w].Rp[i][2]-W[w].Fp[i][2]*dt)
                  )/(2.*dt));
  }
  dexp += 2.*(Wp[w].U - W[w].U); // with detailed balance*/

  // The Metropolis code
  if(dexp > 0.) {
    accepted++;
    CopyWalkerCoord(&W[w], &Wp[w]);
    W[w].U = Wp[w].U;
    return;
  }
  else {
    xi = Random();
    if(xi<Exp(dexp)) {
      accepted++;
      CopyWalkerCoord(&W[w], &Wp[w]);
      W[w].U = Wp[w].U;
      return;
    }
    else {
      rejected++;
      //if(MC == DIFFUSION) VMCMoveDrift(w);
      return;
    }
  }
}

/***************************** Smart MC move *********************************/
//   R' = R + dt F(R) + dt chi
//   p_chi(dr) ~ Exp(-0.5/dr^2)
//   p = min [1, Exp(w)]
//   w = 2(U(r')-U(r)) - 0.5/dt [F(r')+F(r)] [(F(r')-F(r))dt+2(r'-r)]
void VMCMoveDriftOneByOne(int w) {
  int i,k,np;
  DOUBLE xi;
  DOUBLE dx, dy, dz;
  DOUBLE x,y,z;
  DOUBLE dexp = 0.;
  DOUBLE dt;
  DOUBLE r,r2,Force,du;
  DOUBLE F[3],Fp[3];
  DOUBLE Wp[3] = {0,0,0};
  DOUBLE dr[3] = {0,0,0};

#ifdef SCALABLE
  int cc;
  int j;

  UpdateScalableTable(&W[w]);
#endif

  dt = dt_vmc;

  // Move particle np of the walker w
  for(np=0; np<N; np++) {
    F[0] = F[1] = F[2] = 0.;
    Fp[0] = Fp[1] = Fp[2] = 0.;
    du = 0.;
    x = W[w].x[np];
    y = W[w].y[np];
    z = W[w].z[np];

#ifdef SCALABLE // Calculate Force and weight of the old position
    cc = FindCell(x, y, z); // current cell
    for(j=0; j<W[w].c_Nall[cc]; j++) { // loop over particles in the nearest cells
      i = W[w].c_index_all[cc][j]; // index of a neighbour particle
#else
    for(i=0; i<N; i++) {
#endif
      if(i != np) {
        dr[0] = x - W[w].x[i];
        dr[1] = y - W[w].y[i];
        dr[2] = z - W[w].z[i];
        r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
        if(CheckInteractionConditionWF(dr[0], dr[1], dr[2], r2)) {
          r = Sqrt(r2);
          du -= InterpolateU(&G, r, W[w].spin[np], W[w].spin[i]);
          Force = InterpolateFp(&G, r, W[w].spin[np], W[w].spin[i]) / r;
          for(k=0; k<3; k++) {
            F[k] += Force * dr[k];
          }
        }
      }
    }

    RandomNormal3(&dx, &dy, &dz, 0., Sqrt(dt)); // Gaussian jump

    dx += dt * F[0]; // drift move
    dy += dt * F[1];
    dz += dt * F[2];

#ifdef SCALABLE // adjust timestep
    if(dx>Lhalfwf) {
      Warning("  Move One By One: decrease timestep, particle jumps out of the cell!\n");
      Warning("  Decreasing the timestep %" LF " -> %" LF "\n", dt_vmc, dt_vmc*0.1);
      dt_vmc *= 0.1;
      dt *= 0.1;
      Warning("  Increasing accepance rate %" LF " -> %" LF "\n", acceptance_rate, 0.5*(acceptance_rate+1));
    }
#endif

    Wp[0] = x + dx; // Wp is the trial position
    Wp[1] = y + dy;
    Wp[2] = z + dz;

    ReduceToTheBox(Wp);

#ifdef SCALABLE // calculate new value of the wavefunction and new drift
    cc = FindCell(x, y, z); // current cell
    for(j=0; j<W[w].c_Nall[cc]; j++) { // loop over particles in the nearest cells
      i = W[w].c_index_all[cc][j]; // index of a neighbour particle
#else
    for(i=0; i<N; i++) {
#endif
      if(i != np) {
        dr[0] = Wp[0] - W[w].x[i];
        dr[1] = Wp[1] - W[w].y[i];
        dr[2] = Wp[2] - W[w].z[i];
        r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
        r = Sqrt(r2);
        if(CheckInteractionConditionWF(dr[0], dr[1], dr[2], r2)) {
          du += InterpolateU(&G, r, W[w].spin[np], W[w].spin[i]);
          Force = InterpolateFp(&G, r, W[w].spin[np], W[w].spin[i]) / r;
          for(k=0; k<3; k++) {
            Fp[k] += Force * dr[k];
          }
        }
      }
    }

    dexp =  (Fp[0]+F[0])*((Fp[0]-F[0])*dt+2.*dx);
    dexp += (Fp[1]+F[1])*((Fp[1]-F[1])*dt+2.*dy);
    dexp += (Fp[2]+F[2])*((Fp[2]-F[2])*dt+2.*dz);
    dexp = -0.5*dexp + 2.*du;

    // The Metropolis code
    // f = Exp(u), p = Exp(up-u)
    if(dexp > 0) {
      W[w].x[np] = Wp[0];
      W[w].y[np] = Wp[1];
      W[w].z[np] = Wp[2];
      W[w].U -= du;
      accepted++;
    }
    else {
      xi = Random();
      if(xi<Exp(dexp)) {
        W[w].x[np] = Wp[0];
        W[w].y[np] = Wp[1];
        W[w].z[np] = Wp[2];
        W[w].U -= du;
        accepted++;
      }
      else {
        rejected++;
      }
    }
  }
}

/***************************** VMCMoveSumOfTwoWavefunctions **************************************/
void VMCMoveSumOfTwoWavefunctions(int w) {
  int i;
  DOUBLE xi;
  DOUBLE dx,dy,dz;
  DOUBLE x,y,z;
  DOUBLE sigma;
  int overlap;
  DOUBLE scale;
  int accept;
  DOUBLE U1, U2;

  sigma = Sqrt(dt_vmc);

  //if(iteration_global == 1) { // initialize first time
    b = sum_of_two_wavefunctions_parameter1; // wavefunction 1
    U1 = U(W[w]);
    b = sum_of_two_wavefunctions_parameter2; // wavefunction 2
    U2 = U(W[w]);
    W[w].U = U2 + log(sum_of_two_wavefunctions_alpha*exp(U1-U2) + 1. - sum_of_two_wavefunctions_alpha); // sum of two wave functions
  //}

   if(fabs(W[w].U - U2)>1e-3)
    Message("oh");

  // Move all particles of the walker w
  for(i=0; i<N; i++) {
    x = W[w].x[i];
    y = W[w].y[i];
    z = W[w].z[i];

    RandomNormal3(&dx, &dy, &dz, 0., sigma);

    // wp is the trial position
    CaseX(Wp[w].x[i] = x + dx);
    CaseY(Wp[w].y[i] = y + dy);
    CaseZ(Wp[w].z[i] = z + dz);
#ifdef SPINFULL
    Wp[w].spin[i] = W[w].spin[i];
#endif
  }

#ifdef CENTER_OF_MASS_IS_NOT_MOVED // propose such movement that CM position is not moved
  AdjustCenterOfMassWalker(&Wp[w]);
#endif
#ifdef CENTER_OF_MASS_Z_IS_NOT_MOVED // propose such movement that CM position is not moved
  AdjustCenterOfMassWalker(&Wp[w]);
#endif
#ifdef CENTER_OF_MASS_IS_NOT_MOVED_TWO_PARTS
  AdjustCenterOfMassWalker(&Wp[w]);
#endif
#ifdef FIRST_PARTICLE_IS_NOT_MOVED // propose such movement that the first particle is not moved
  AdjustCenterOfMassWalker(&Wp[w]);
#endif

  ReduceWalkerToTheBox(&Wp[w]);

  b = sum_of_two_wavefunctions_parameter1; // wavefunction 1
  U1 = U(Wp[w]);
  b = sum_of_two_wavefunctions_parameter2; // wavefunction 2
  U2 = U(Wp[w]);
  // sum of two wave functions
  //Psi = sum_of_two_wavefunctions_alpha * exp(U1) + (1.-sum_of_two_wavefunctions_alpha) * exp(U2);
  //Utotal = log(Psi);
  //Psi = exp(U1)*(sum_of_two_wavefunctions_alpha + (1.-sum_of_two_wavefunctions_alpha) * exp(U2-U1));
  //Utotal = U1 + log((sum_of_two_wavefunctions_alpha + (1.-sum_of_two_wavefunctions_alpha) * exp(U2-U1)));
  Wp[w].U = U2 + log(sum_of_two_wavefunctions_alpha*exp(U1-U2) + 1. - sum_of_two_wavefunctions_alpha);
  if(isnan(U1) || isinf(U1)) {
    Warning("NaN encountered in Metropolis move (U1), exiting\n");
    return;
  }
  if(isnan(U2) || isinf(U2)) {
    Warning("NaN encountered in Metropolis move (U2), exiting\n");
    return;
  }
  if(isnan(Wp[w].U) || isinf(Wp[w].U)) {
    Warning("NaN encountered in Metropolis move (U), exiting\n");
    return;
  }
  Wp[w].weight = 1;

  //-infinity is returned if particles overlap in new configuration
  if(Wp[w].U == -INFTY) {
    overlaped++;
    rejected++;
    return;
  }

#ifdef SECURE
  if(CheckWalkerOverlapping(Wp[w]))
    Error("MoveAll error");
#endif

  // The Metropolis code
  if(Wp[w].U > W[w].U) {// f^2 = Exp(2u), p = Exp(2(up-u))
    accept = ON;
  }
  else {
    xi = Random();
    if(xi<Exp(2.*(Wp[w].U-W[w].U))) {
      accept = ON;
    }
    else {
      rejected++;
      return;
    }
  }

  if(accept) {
    accepted++;
    CopyWalker(&W[w], &Wp[w]);
    return;
  }
}

/********************************* U *****************************************/
DOUBLE U_sum_of_two_wavefunctions(struct Walker *W) {
  DOUBLE U1, U2, r;
  DOUBLE Utotal1, Utotal2, Utotal3;

  W->status_sum_of_two_walkers_U = 1;
  b = sum_of_two_wavefunctions_parameter1; // parameters of wavefunction 1
  U1 = U(*W);
  b = sum_of_two_wavefunctions_parameter2; // parameters of wavefunction 2
  U2 = U(*W);
  W->status_sum_of_two_walkers_U = 0;

  r = exp(U1-U2);
  /*if(r<1e-6) {
    r = sum_of_two_wavefunctions_alpha/(1.-sum_of_two_wavefunctions_alpha)*r;
    // ln(1+x) = x-x*x/2+...
    return U2 + log(1. - sum_of_two_wavefunctions_alpha) + r - 0.5*r*r;
  }*/

  if(r<DBL_EPSILON) // 2e-16
    Warning("r out of range");

  //Utotal1 = log(sum_of_two_wavefunctions_alpha*exp(U1) + (1. - sum_of_two_wavefunctions_alpha)*exp(U2));
  //Utotal2 = U2 + log(sum_of_two_wavefunctions_alpha*r + 1. - sum_of_two_wavefunctions_alpha);
  Utotal3 = U1 + log(sum_of_two_wavefunctions_alpha + (1. - sum_of_two_wavefunctions_alpha)/r);

  return Utotal3; // sum of two wave functions

  //return 
  //return U2 + log(sum_of_two_wavefunctions_alpha*r + 1. - sum_of_two_wavefunctions_alpha); // sum of two wave functions
}

/*********************** Excluded Volume Coordinates Wrap ********************/
// change z[i] -> z[i] - i*Apar
void ExcludedVolumeCoordinatesWrap(struct Walker Walker) {
  int i,j;
  int sorted = ON;
  DOUBLE zmin;

#ifndef TRIAL_TONKS_SIMPLE /* (1D) Tonks w.f. |sin(sqrt{Bpar}*(x-a))|, Rpar=0, a=0 means TG */
  Error("use TRIAL_TONKS_SIMPLE wave function with a=0 and use Apar for defining the scattering length");
#endif

  // check if z positions are already sorted
  for(i=1; i<N; i++) {
    if(Walker.z[i] < Walker.z[i-1]) {
      sorted = OFF;
      break;
    }
  }

  if(sorted == OFF) {
    qsort(Walker.z, N, sizeof(DOUBLE), compare_double);
  }

  sorted = ON;
  // check if z positions are already sorted
  for(i=1; i<N; i++) {
    if(Walker.z[i] < Walker.z[i-1]) {
      sorted = OFF;
      break;
    }
  }
  if(sorted == OFF) 
    Message("strange: sorted %i, ", sorted);

   // apply the Girardeau's excluded volume transformation
  for(i=0; i<N; i++) {
    Walker.z[i] -= Apar*(DOUBLE) i;
  }

  // adjust w.f
  // (1D) Tonks w.f. |sin(sqrt{Bpar}*(x-a))|, Rpar=0, a=0 means TG
  // return Log(fabs(Sin(sqE*(x-a))));
  L = (DOUBLE) N / n;
  L = L - (DOUBLE)(N-1)*Apar;
  Lz = L;
  L_inv_z = 1./Lz;
  Lhalf = L/2.;
  L2 = L*L;
  sqE = PI / L;
  kL = sqE;
  a = 0.;
  a2 = 0.;
}

/*********************** Excluded Volume Coordinates Wrap ********************/
// change z[i] -> z[i] + i*Apar
void ExcludedVolumeCoordinatesUnwrap(struct Walker Walker) {
  int i;

  // undo the Girardeau's excluded volume transformation
  for(i=0; i<N; i++) {
    Walker.z[i] += Apar*(DOUBLE) i;
  }

  // adjust w.f
  // (1D) Tonks w.f. |sin(sqrt{Bpar}*(x-a))|, Rpar=0, a=0 means TG
  // return Log(fabs(Sin(sqE*(x-a))));
  L = (DOUBLE) N / n;
  Lz = L;
  L_inv_z = 1./Lz;
  Lhalf = L/2.;
  L2 = L*L;
  sqE = PI / L;
  kL = sqE;
  a = Apar;
  a2 = Apar*Apar;
}

DOUBLE U(struct Walker Walker) {
  int i,j;
  DOUBLE u = 0.;
  DOUBLE dr[3] = {0,0,0};
  DOUBLE x,y,z;
  DOUBLE r2;

#ifdef SCALABLE
  int cc,ci,cj;
#endif 
#ifdef THREE_BODY_TERMS
  DOUBLE CM12[3]; // CM position
  DOUBLE rho2;
  DOUBLE HR; // hyperradius
  int k;
#endif

#ifdef DO_NOT_CALCULATE_ENERGY_AND_U
  return 0.;
#endif

#ifdef HARD_SPHERE_EXCLUDED_VOLUME_CHANGE_OF_VARIABLES
  if(CheckWalkerOverlapping(Walker)) return -INFTY;
  ExcludedVolumeCoordinatesWrap(Walker); // convert variables x_i -> x_i - i * a
#endif

  if(MC == VARIATIONAL && SmartMC == VMC_SUM_OF_TWO_WAVEFUNCTIONS && Walker.status_sum_of_two_walkers_U == 0) {
    return U_sum_of_two_wavefunctions(&Walker);
  }

// Effect of the external field
#ifdef ONE_BODY_TRIAL_TERMS
  u = OneBodyUWalker(Walker);
#endif

// 2-body Jastrow contribution
#ifdef SCALABLE
  for(cc=0; cc<Ncells; cc++) { // loop over cells
    for(ci=0; ci<Walker.c_Nlocal[cc]; ci++) { // loop over particles in the local cell
      i = Walker.c_index_local[cc][ci];
      x = Walker.x[i];
      y = Walker.y[i];
      z = Walker.z[i];
      for(cj=0; cj<Walker.c_Nall[cc]; cj++) { // loop over particles in the nearest cells avoiding DOUBLE counting
        j = Walker.c_index_all[cc][cj];
        if(j>i) {
          dr[0] = x - Walker.x[j];
          dr[1] = y - Walker.y[j];
          dr[2] = z - Walker.z[j];
#ifdef TRIAL_1D
          dr[0] = dr[1] = 0.;
#endif
#ifdef TRIAL_2D
          dr[2] = 0.;
#endif
          r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
#ifdef HARD_SPHERE
          if(r2 <= a2) return -INFTY;
#endif
          if(r2<G.max2 && CheckInteractionConditionWF(dr[0], dr[1], dr[2], r2)) {
            u += InterpolateU(&G, Sqrt(r2), Walker.spin[i], Walker.spin[j]);
          }
        }
      }
    }
  }
#else
  for(i=0; i<N; i++) {
    CaseX(x = Walker.x[i]);
    CaseY(y = Walker.y[i]);
    CaseZ(z = Walker.z[i]);
    for(j=i+1; j<N; j++) {
      CaseX(dr[0] = x - Walker.x[j]);
      CaseY(dr[1] = y - Walker.y[j]);
      CaseZ(dr[2] = z - Walker.z[j]);
      r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
#ifdef HARD_SPHERE
     if(r2 <= a2) {
#ifdef HARD_SPHERE_EXCLUDED_VOLUME_CHANGE_OF_VARIABLES // undo changes
         ExcludedVolumeCoordinatesUnwrap(Walker); // convert variables x_i -> x_i - i * a
#endif
         return -INFTY;
      }
#endif
      if(CheckInteractionConditionWF(dr[0], dr[1], dr[2], r2)) {
        u += InterpolateU(&G, Sqrt(r2), Walker.spin[i], Walker.spin[j]);
      }
    }
  }
#endif

#ifdef THREE_BODY_TERMS
  for(i=0; i<N; i++) {
    for(j=i+1; j<N; j++) {
      // distance (1-2)
      CaseX(dr[0] = Walker.x[i] - Walker.x[j]);
      CaseY(dr[1] = Walker.y[i] - Walker.y[j]);
      CaseZ(dr[2] = Walker.z[i] - Walker.z[j]);
      rho2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);

      // center of mass (1-2) position
      CaseX(CM12[0] = Walker.x[i] - 0.5*dr[0]);
      CaseY(CM12[1] = Walker.y[i] - 0.5*dr[1]);
      CaseZ(CM12[2] = Walker.z[i] - 0.5*dr[2]);

      for(k=j+1; k<N; k++) {
        // distance (CM12 - 3)
        CaseX(dr[0] = Walker.x[k] - CM12[0]);
        CaseY(dr[1] = Walker.y[k] - CM12[1]);
        CaseZ(dr[2] = Walker.z[k] - CM12[2]);
        r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);

        HR = sqrt(4./3.*r2 + rho2);

        u += u3(HR);
      }
    }
  }
#endif

#ifdef HARD_SPHERE_EXCLUDED_VOLUME_CHANGE_OF_VARIABLES
  ExcludedVolumeCoordinatesUnwrap(Walker); // convert variables x_i -> x_i - i * a
#endif

  return u;
}

/********************************* U vector **********************************/
// needed for numerical calculation of derivatives
DOUBLE Uvector(DOUBLE **R, int w) {
  int i,j;
  DOUBLE u = 0.;
  DOUBLE dr[3] = {0,0,0};
  DOUBLE x,y,z;
  DOUBLE r2;
#ifdef THREE_BODY_TERMS
  DOUBLE CM12[3]; // CM position
  DOUBLE rho2;
  DOUBLE HR; // hyperradius
  int k;
#endif

  for(i=0; i<N; i++) {
    x = R[i][0];
    y = R[i][1];
    z = R[i][2];
#ifdef ONE_BODY_TRIAL_TERMS
    u += OneBodyU(x, y, z, i, W[w].w); // Effect of the external field
#endif
    for(j=i+1; j<N; j++) {
      CaseX(dr[0] = x - R[j][0]);
      CaseY(dr[1] = y - R[j][1]);
      CaseZ(dr[2] = z - R[j][2]);
      r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
#ifdef HARD_SPHERE
      if(r2 <= a2) {
         return -INFTY;
      }
#endif
      if(CheckInteractionConditionWF(dr[0], dr[1], dr[2], r2)) {
        u += InterpolateU(&G, Sqrt(r2), W[w].spin[i], W[w].spin[j]);
      }
    }
  }

#ifdef THREE_BODY_TERMS
  for(i=0; i<N; i++) {
    for(j=i+1; j<N; j++) {
      // distance (1-2)
      dr[0] = R[i][0] - R[j][0];
      dr[1] = R[i][1] - R[j][1];
      dr[2] = R[i][2] - R[j][2];
      rho2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);

      // center of mass (1-2) position
      CM12[0] = R[i][0] - 0.5*dr[0];
      CM12[1] = R[i][1] - 0.5*dr[1];
      CM12[2] = R[i][2] - 0.5*dr[2];

      for(k=j+1; k<N; k++) {
        // distance (CM12 - 3)
        dr[0] = R[k][0] - CM12[0];
        dr[1] = R[k][1] - CM12[1];
        dr[2] = R[k][2] - CM12[2];
        r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);

        HR = sqrt(4./3.*r2 + rho2);

        u += u3(HR);
      }
    }
  }
#endif

  return u;
}

/********************************* U *****************************************/
// for a given walker, calculate the effect of spin inversion psi_T(R,-sigma_i) / psi_T(R,sigma_i)
// and store in psiT_sigma_inv[]
// two component system is assumed
void U_psiT_sigma_inv(struct Walker Walker) {
#ifdef SPINFULL_TUNNELING
  int i,j;
  DOUBLE u = 0.;
  DOUBLE dr[3] = {0,0,0};
  DOUBLE x,y,z;
  DOUBLE r2;

// Effect of the external field
//#ifdef ONE_BODY_TRIAL_TERMS
//  u = OneBodyUWalker(Walker);
//#endif

  for(i=0; i<N; i++) { // invert spin of particle i
    Walker.psiT_sigma_inv[i] = 0.;
    CaseX(x = Walker.x[i]);
    CaseY(y = Walker.y[i]);
    CaseZ(z = Walker.z[i]);
    for(j=0; j<N; j++) {
      if(j != i) {
        CaseX(dr[0] = x - Walker.x[j]);
        CaseY(dr[1] = y - Walker.y[j]);
        CaseZ(dr[2] = z - Walker.z[j]);
        r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
        if(CheckInteractionConditionWF(dr[0], dr[1], dr[2], r2)) {
          Walker.psiT_sigma_inv[i] += InterpolateU(&G, Sqrt(r2), !Walker.spin[i], Walker.spin[j]) - InterpolateU(&G, Sqrt(r2), Walker.spin[i], Walker.spin[j]);
        }
      }
    }
  }
#endif
}

/************************ Check Particles in the box *************************/
int CheckParticlesInTheBox(const struct Walker W) {
  int i;

  for(i=0; i<N; i++) {
    if(W.x[i] <0 || W.x[i]>Lx || W.y[i] <0 || W.y[i]>Ly || W.z[i] <0 || W.z[i]>Lz) 
      Error("Particle out side Particle %i (%" LF " %" LF " %" LF ")", i, W.x[i], W.y[i], W.z[i]);
  }

  return 0;
}

/******************** Check Interaction Condition ****************************/
int CheckInteractionConditionWF(const DOUBLE x, const DOUBLE y, const DOUBLE z, const DOUBLE r2) {

#ifdef BC_3DPBC // Homogeneous 3D system
  if(r2 < Lhalfwf2)
    return 1;
  else
    return 0;
#endif

#ifdef BC_2DPBC
  if(r2 < Lhalfwf2)
    return 1;
  else
    return 0;
#ifdef SECURE
  if(z != 0) Error("2D system: non zero z component is encountered!");
#endif
#endif

#ifdef  BC_2DPBC_NON_ORTHOGONAL
#error "CheckInteractionCondition not implemented"
#endif

#ifdef BC_1DPBC_Z // 1D tube
  if(fabs(z)<Lhalfwf)
    return 1;
  else
    return 0;
#endif

#ifdef BC_1DPBC_X // 1D tube
  if(fabs(x)<Lhalfwf)
    return 1;
  else
    return 0;
#endif

#ifdef BC_SPHERE
  return 1;
#endif

#ifdef BC_ABSENT
  return 1;
#endif
}

int CheckInteractionConditionPotential(const DOUBLE x, const DOUBLE y, const DOUBLE z, const DOUBLE r2) {

#ifdef TRIAL_2D // 2D hom.
#ifdef SECURE
  if(z != 0) Error("2D system: non zero z component is encountered!");
#endif
#endif

#ifdef BC_3DPBC // Homogeneous 3D system
  if(r2 < Lcutoff_pot2)
    return 1;
  else
    return 0;
#endif

#ifdef BC_2DPBC // 2D tube
  if(r2 < Lcutoff_pot2)
    return 1;
  else
    return 0;
#endif

#ifdef BC_1DPBC_Z // 1D tube
  if(fabs(z)<Lcutoff_pot)
    return 1;
  else
    return 0;
#endif

#ifdef BC_1DPBC_X // 1D tube
  if(fabs(x)<Lcutoff_pot)
    return 1;
  else
    return 0;
#endif

#ifdef BC_SPHERE
  return 1;
#endif

#ifdef BC_ABSENT
  return 1;
#endif
}

/***************************** Update Scalable table **************************/
void UpdateScalableTable(struct Walker *W) {
  int i,j,cc;
  int s[27], Nnear;

  for(i=0; i<Ncells; i++) W->c_Nall[i] = W->c_Nlocal[i] = 0;

  for(i=0; i<N; i++) {
    cc = FindCell(W->x[i], W->y[i], W->z[i]);
//#ifdef SECURE
   if(cc>=Ncells)
     Error("problems in Update Scalable Table (2)");

   if(W->c_Nlocal[cc]<0)
     Error("problems in Update Scalable Table (3)");
//#endif

    W->c_index_local[cc][W->c_Nlocal[cc]] = i;
    W->c_Nlocal[cc]++;

    Nnear = FindNearestCells(W->x[i],W->y[i],W->z[i], s);
    for(j=0; j<Nnear; j++) {
      W->c_index_all[s[j]][W->c_Nall[s[j]]] = i;
      W->c_Nall[s[j]]++;
    }
#ifdef SECURE
    if(W->c_Nall[cc]>N) 
      Error("impossible situation in UpdateScalableTable");
#endif
  }
#ifdef SECURE
  for(cc=0; cc<Ncells; cc++) { // loop over cells
    for(j=0; j<W->c_Nlocal[cc]; j++) { // loop over particles in the local cell
      i = W->c_index_local[cc][j];
      if(cc != FindCell(W->x[i], W->y[i], W->z[i]))
        Error("Update Scalable table - no self consistency!");
    }
  }
#endif
}

/***************************** Find Cell **************************************/
// point to the cell in which the particle is contained
int FindCell(DOUBLE x, DOUBLE y, DOUBLE z) {
#ifdef TRIAL_1D
  return (int) (z/Lcell);
#endif
#ifdef TRIAL_2D
  #ifdef SECURE
    int i1,i2;
    i1 = (int) (x/Lcell);
    i2 = (int) (y/Lcell);
    if(i1<0) {
      Warning("  paricle outside the box (1)!\n");
      i1 += Nscal;
    }
    if(i2<0) {
      Warning("  paricle outside the box (2)!\n");
      i2 += Nscal;
    }
    if(i1>=Nscal) {
      Warning("  paricle outside the box (3)!\n");
      i1 -= Nscal;
    }
    if(i2>=Nscal) {
      Warning("  paricle outside the box (4)!\n");
      i2 -= Nscal;
    }
    return i1 + i2*Nscal;
  #else
    return (int) (x/Lcell) + (int) (y/Lcell)*Nscal;
  #endif
#endif
#ifdef TRIAL_3D
  return (int) (x/Lcell) + (int) (y/Lcell)*Nscal + (int) (z/Lcell)*Nscal*Nscal;
#endif
}

/***************************** Find Cell **************************************/
int FindNearestCells(DOUBLE x, DOUBLE y, DOUBLE z, int *s) {
// points to the current and nearest cells
// returns number of nearest cells
#ifdef TRIAL_1D
  int sz, i, index, sznear;
  sz = (int) (z/Lcell);
  index = 0;
  for(i=-1; i<=1; i++) {
    sznear = sz + i;
    if(sznear == -1) sznear = Nscal-1;
    if(sznear == Nscal) sznear = 0;
    s[index] = sznear;
    index++;
  }
  if(Nscal>=3)
    return 3;
  else
    return Nscal;
#endif

#ifdef TRIAL_2D
  int sx, sy, i,j,index,sxnear,synear;
  sx = (int) (x/Lcell);
  sy = (int) (y/Lcell);
  index = 0;
  for(i=-1; i<=1; i++) {
    sxnear = sx + i;
    if(sxnear == -1) sxnear = Nscal-1;
    if(sxnear == Nscal) sxnear = 0;
    for(j=-1; j<=1; j++) {
      synear = sy + j;
      if(synear == -1) synear = Nscal-1;
      if(synear == Nscal) synear = 0;
      s[index] = sxnear + Nscal*synear;
#ifdef SECURE
      if(s[index]>Ncells) 
        Error("scalable code - impossible cell index (1)");
#endif
      index++;
    }
  }
  if(Nscal>=3)
    return 9;
  else
    return 1;
#endif

#ifdef TRIAL_3D
  int sx, sy, sz, i,j,index,sxnear,synear,sznear;
  sx = (int) (x/Lcell);
  sy = (int) (y/Lcell);
  sz = (int) (z/Lcell);
  index = 0;
  for(i=-1; i<=1; i++) {
    sxnear = sx + i;
    if(sxnear == -1) sxnear = Nscal-1;
    if(sxnear == Nscal) sxnear = 0;
    for(i=-1; i<=1; i++) {
      synear = sy + i;
      if(synear == -1) synear = Nscal-1;
      if(synear == Nscal) synear = 0;
      for(j=-1; j<=1; j++) {
        sznear = sz + j;
        if(sznear == -1) sznear = Nscal-1;
        if(sznear == Nscal) sznear = 0;
        s[index] = sxnear + Nscal*synear + Nscal*Nscal*sznear;
#ifdef SECURE
        if(s[index]>Ncells) 
          Error("scalable code - impossible cell index (2)");
#endif
        index++;
      }
    }
  }
  if(Nscal>=3)
    return 27;
  else
    return 1;
#endif
}

/***************************** ClassicalMoveOneByOne *******************************/
void ClassicalMoveOneByOne(int w) {
  int i,j;
  DOUBLE xi;
  DOUBLE Wp[3] = {0,0,0};
  DOUBLE dr[3] = {0,0,0};
  DOUBLE r2;
  DOUBLE x = 0.;
  DOUBLE y = 0.;
  DOUBLE z = 0.;
  DOUBLE dx,dy,dz;
  DOUBLE dE;
  int overlap;

#ifdef CENTER_OF_MASS_IS_NOT_MOVED // propose such movement that CM position is not moved
  AdjustCenterOfMassWalker(&W[w]);
#endif
#ifdef FIRST_PARTICLE_IS_NOT_MOVED
  AdjustCenterOfMassWalker(&W[w]);
#endif

  // Move particle np of the walker w
  for(i=0; i<N; i++) {
    dE = 0.;

    overlap = OFF;

    CaseX(x = W[w].x[i]);
    CaseY(y = W[w].y[i]);
    CaseZ(z = W[w].z[i]);

    // Wp is the trial position
    RandomNormal3(&dx, &dy, &dz, 0., Sqrt(dt_vmc));

    CaseX(Wp[0] = x + dx);
    CaseY(Wp[1] = y + dy);
    CaseZ(Wp[2] = z + dz);

    ReduceToTheBox(Wp);

#ifdef HARD_SPHERE // Check overlapping
#ifdef BC_1DPBC_Z
    for(j=0; j<N; j++) {
      if(j != i) {
        if(fabs(W[w].z[j] - Wp[2])<a) overlap = ON;
        if(fabs(W[w].z[j] - Wp[2]-L)<a) overlap = ON;
        if(fabs(W[w].z[j] - Wp[2]+L)<a) overlap = ON;
        //if(fabs(W[w].z[j] - (DOUBLE)j*a - (Wp[2]-(DOUBLE)i*a))<a) overlap = ON;
      }
    }
#endif
#endif

#ifdef EXTERNAL_POTENTIAL // Effect of the trapping
    //du -= alpha_x*dx*(2.*x+dx) + alpha_y*dy*(2.*y+dy) + alpha_z*dz*(2.*z+dz);
    //dE += alpha_x*dx*(2.*x+dx) + alpha_y*dy*(2.*y+dy) + alpha_z*dz*(2.*z+dz);
    dE += Vext(Wp[0], Wp[1], Wp[2]) - Vext(x, y, z);
#endif

    for(j=0; j<N; j++) {
      if(j != i) {
        // calculate new value of energy
        CaseX(dr[0] = Wp[0] - W[w].x[j]);
        CaseY(dr[1] = Wp[1] - W[w].y[j]);
        CaseZ(dr[2] = Wp[2] - W[w].z[j]);
        r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
#ifdef INTERACTION_WITHOUT_IMAGES // might be restricted to r < L/2
        if(CheckInteractionConditionPotential(dr[0], dr[1], dr[2], r2)) dE += InteractionEnergy_ij(Sqrt(r2), W[w].spin[i], W[w].spin[j]);
#else 
        dE += InteractionEnergyWithImages_ij(dr[0], dr[1], dr[2]);
#endif

        // calculate old value of the wavefunction
        CaseX(dr[0] = x - W[w].x[j]);
        CaseY(dr[1] = y - W[w].y[j]);
        CaseZ(dr[2] = z - W[w].z[j]);
        r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
#ifdef INTERACTION_WITHOUT_IMAGES // might be restricted to r < L/2
        if(CheckInteractionConditionPotential(dr[0], dr[1], dr[2], r2)) dE -= InteractionEnergy_ij(Sqrt(r2), W[w].spin[i], W[w].spin[j]);
#else 
        dE -= InteractionEnergyWithImages_ij(dr[0], dr[1], dr[2]);
#endif

      }
    }

    // The Metropolis code
    if(dE <= 0. && overlap == OFF) {
      CaseX(W[w].x[i] = Wp[0]);
      CaseY(W[w].y[i] = Wp[1]);
      CaseZ(W[w].z[i] = Wp[2]);
      W[w].E += dE;
      accepted++;
    }
    else {
      xi = Random();
      if(xi<Exp(-dE/T) && overlap == OFF) {
        CaseX(W[w].x[i] = Wp[0]);
        CaseY(W[w].y[i] = Wp[1]);
        CaseZ(W[w].z[i] = Wp[2]);
        W[w].E += dE;
        accepted++;
      }
      else {
        rejected++;
      }
    }
  }
}

/***************************** Check All Walkers Are Alive *************************/
int CheckWalkersAlive(void) {
  int status = 0;
  int w;

  for(w=0; w<Nwalkers; w++) {
    if(!W[w].status == ALIVE) {
      Warning("Check Walkers Alive: walker W[%i] has status %i\n", w, W[w].status);
      status = 1;
    }
  }
  return status;
}

/***************************** Check Walkers in the box ****************************/
int CheckWalkersInTheBox(void) {
  int status = OFF;
  int w,i;

  for(w=0; w<Nwalkers; w++) {
    for(i=0; i<N; i++) {
      CaseX(if(W[w].x[i]>L || W[w].x[i]<0) status = ON);
      CaseY(if(W[w].y[i]>L || W[w].y[i]<0) status = ON);
      CaseZ(if(W[w].z[i]>L || W[w].z[i]<0) status = ON);
      if(status) {
        Warning(" Check Walkers in the box failed for w=%w, i=%i, ", w, i, W[w].x[i], W[w].y[i], W[w].z[i]);
        return 1;
      }
    }
  }

  return status;
}

/********************************* U One Body ********************************/
// returns w.f. of given walker in form f = Exp(U)
DOUBLE OneBodyUWalker(struct Walker Walker) {

#ifndef ONE_BODY_TRIAL_TERMS // homogeneous liquid phase
  return 0;
#else
  int i;
  DOUBLE u = 0.;

#ifdef CRYSTAL
#ifdef CRYSTAL_SYMMETRIC
  DOUBLE dr[3] = {0,0,0};
  DOUBLE r2;
  DOUBLE ucr;
  int j;
#endif
#endif

#ifdef MANY_BODY_TRIAL_TRAP_CENTER_OF_MASS_NODE_Z
  DOUBLE z_CM = 0.;

  for(i=0; i<N; i++) z_CM += Walker.z[i];
  u += log(fabs(z_CM)); // (z1 + ... +zN)
#endif

#ifdef MANY_BODY_TRIAL_TRAP_CENTER_OF_MASS_HERMITE2_Z
  DOUBLE z2 = 0.;

  for(i=0; i<N; i++) z2 += (Walker.z[i] / (DOUBLE) N)*(Walker.z[i] / (DOUBLE) N);
  u += log(fabs(z2 - beta));
#endif

  for(i=0; i<N; i++) 
    u += OneBodyU(Walker.x[i], Walker.y[i], Walker.z[i], i, Walker.w);

#ifdef CRYSTAL
#ifdef CRYSTAL_SYMMETRIC
  for(i=0; i<Crystal.size; i++) { // sum over lattice sites
    ucr = 0.;
    for(j=0; j<N; j++) { // sum over particles
      dr[0] = Walker.x[j] - Crystal.x[i];
      dr[1] = Walker.y[j] - Crystal.y[i];
      dr[2] = Walker.z[j] - Crystal.z[i];
      r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
      //if(CheckInteractionConditionWF(dr[0], dr[1], dr[2], r2))  
      //restriction in the summation leads to strange results in Crystal.R = 0 case
        //ucr += Exp(-Crystal.R * r2);
        ucr += Crystal_dot_weight_j*Exp(-Crystal_dot_Rx_i*dr[0]*dr[0]-Crystal_dot_Ry_i*dr[1]*dr[1]-Crystal_dot_Rz_i*dr[2]*dr[2]);
    }
    u += Log(ucr);
  }
#endif
#endif

  return u;
#endif
}

/********************************* One Body U ********************************/
DOUBLE OneBodyU(DOUBLE x, DOUBLE y, DOUBLE z, int i, int w) {
#ifndef ONE_BODY_TRIAL_TERMS // i.e. homogeneous liquid system
  return 0;
#else
  DOUBLE u = 0.;

#ifdef CRYSTAL
#ifdef CRYSTAL_NONSYMMETRIC
  DOUBLE dr[3],r2;
#endif
#endif

#ifdef ONE_BODY_IMPURITY
  int j;
  DOUBLE dr[3],r2;
#endif

#ifdef ONE_BODY_TRIAL_TRAP // i.e. in a trap
#ifdef TRIAL_3D
  u -= alpha_x * x*x + alpha_y * y*y + alpha_z * z*z;
#endif
#ifdef TRIAL_2D
  u -= alpha_x * x*x + alpha_y * y*y;
#endif
#ifdef TRIAL_1D
  u -= alpha_z * z*z;
#endif
#endif

#ifdef ONE_BODY_TRIAL_TRAP_POWER_LAW
  u += -alpha_latt*(x*x + y*y + z*z) + 0.5*beta_latt*log(x*x + y*y + z*z);
  //u += -alpha_latt*(x*x + y*y + z*z) + beta_latt*log(sqrt(x*x + y*y + z*z));
#endif

#ifdef INTERPOLATE_SPLINE_ONE_BODY_Z_WF
  u += InterpolateSplineU(&G1, z);
#endif

#ifdef SPECKLES
  u += SpecklesU(x, y);
#endif

#ifdef ONE_BODY_IMPURITY // impurity
  for(j=0; j<Crystal.size; j++) {
    dr[0] = x - Crystal.x[j];
    dr[1] = y - Crystal.y[j];
    dr[2] = z - Crystal.z[j];
    r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
    u += ImpurityInterpolateU(Sqrt(r2));
  }
#endif

#ifdef ONE_BODY_TRIAL_SIN
#ifdef BC_3DPBC_CUBE
  u += alpha_latt*log(sin(kL*x)*sin(kL*x) + sin(kL*y)*sin(kL*y) + sin(kL*z)*sin(kL*z));
#endif
#ifdef BC_2DPBC_SQUARE
  u += alpha_latt*log(sin(kL*x)*sin(kL*x) + sin(kL*y)*sin(kL*y));
#endif
#ifdef BC_1DPBC_Z
  u += log(beta_latt+pow(fabs(sin(kL*z)),alpha_latt));
#endif
#endif

#ifdef ONE_BODY_TRIAL_ZERO_BOUNDARY_CONDITION
  u += beta_z * log(fabs(sin(PI*z/L)));
#endif

#ifdef ONE_BODY_TRIAL_COS_SERIES 
 u += log(1. + alpha_latt*cos(2.*kL*z) + beta_latt*cos(4.*kL*z) + gamma_latt*cos(6.*kL*z));
#endif 

#ifdef ONE_BODY_TRIAL_ZBC
  u += log(fabs(sin(kL*z)));
#endif

#ifdef ONE_BODY_TRIAL_ABC
  u += alpha_latt*log(fabs(sin(kL*z)));
  //u += log(fabs(sin(kL*z)));
#endif

#ifdef ONE_BODY_TRIAL_1D_WAVE_COS
  u += log(fabs(cos(one_body_momentum_k*z)));
#endif

#ifdef ONE_BODY_TRIAL_1D_WAVE_SIN
  if(one_body_momentum_k>0) u += log(fabs(sin(one_body_momentum_k*z)));
#endif

#ifdef ONE_BODY_TRIAL_1D_WAVE_COS_FINITE
  u += log(fabs(1.+Srecoil*cos(one_body_momentum_k*z)));
#endif

#ifdef ONE_BODY_SOLITON
  u += 0.5*log(solitonV2 + (1.-solitonV2)*tanh(solitonk*(z-Lhalf))*tanh(solitonk*(z-Lhalf)));
#endif

#ifdef CRYSTAL
#ifndef CRYSTAL_SYMMETRIC
#ifdef CHAINS_PINNED_TO_CENTER_OF_MASS
  dr[0] = x - W[w].Crystal_x[i];
  dr[1] = y - W[w].Crystal_y[i];
  dr[2] = z - W[w].Crystal_z[i];
  r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
  u -= Crystal_dot_Rx_i*dr[0]*dr[0] + Crystal_dot_Ry_i*dr[1]*dr[1] + Crystal_dot_Rz_i*dr[2]*dr[2];
#else
  dr[0] = x - Crystal.x[i];
  dr[1] = y - Crystal.y[i];
  dr[2] = z - Crystal.z[i];
  r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
  u -= Crystal_dot_Rx_i * dr[0]*dr[0] + Crystal_dot_Ry_i * dr[1]*dr[1] + Crystal_dot_Rz_i * dr[2]*dr[2];
#endif
#endif
#endif

#ifdef ONE_BODY_TRIAL_DOUBLE_GAUSSIAN
  u += log(exp(-alpha_z * (z-R0par)*(z-R0par)) + exp(-alpha_z * (z+R0par)*(z+R0par)));
#endif

#ifdef ONE_BODY_TRIAL_LATTICE
  DOUBLE f=0, dr[3]={0.,0.,0.},r2;
  int j;

  for(j=0; j<Crystal.size; j++) {
    CaseX(dr[0] = x - Crystal.x[j]);
    CaseY(dr[1] = y - Crystal.y[j]);
    CaseZ(dr[2] = z - Crystal.z[j]);
    r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
    f += exp(-alpha_latt*r2);
  }
  u += log(beta_latt + gamma_latt*f);
#endif

#ifdef ONE_BODY_TRIAL_TRAP_DOWN_ANTISYMMETRIC_Z
  if(W[w].spin[i]) u += log(fabs(z));
#endif

#ifdef ONE_BODY_TRIAL_TRAP_DOWN_HERMITE2_Z
  //if(!up) u += log(fabs(z*z-0.5));
  if(W[w].spin[i]) u += log(fabs(z*z-beta));
#endif

  return u;
#endif
}

/********************************* One Body Fp *******************************/
// Fx = [dF / dx] / F
// Fy = [dF / dy] / F
// Fz = [dF / dz] / F
void OneBodyFp(DOUBLE *Fx, DOUBLE *Fy, DOUBLE *Fz, DOUBLE x, DOUBLE y, DOUBLE z, int i, int w) {
#ifdef CRYSTAL
#ifdef CRYSTAL_NONSYMMETRIC // nonsymmetric crystal
  DOUBLE dr[3]={0.,0.,0},r2;
#endif
#endif
#ifdef ONE_BODY_IMPURITY // add contribution from the w.f. of impurities
  int j;
  DOUBLE dri[3]={0.,0.,0}, ri2, ri, Force;
#endif

#ifdef ONE_BODY_TRIAL_TRAP
  CaseX(*Fx -= two_alpha_x * x);
  CaseY(*Fy -= two_alpha_y * y);
  CaseZ(*Fz -= two_alpha_z * z);
#endif

#ifdef ONE_BODY_TRIAL_TRAP_POWER_LAW
  *Fx += (-2.*alpha_latt + beta_latt/(x*x + y*y + z*z))*x;
  *Fy += (-2.*alpha_latt + beta_latt/(x*x + y*y + z*z))*y;
  *Fz += (-2.*alpha_latt + beta_latt/(x*x + y*y + z*z))*z;
#endif

#ifdef SPECKLES
  SpecklesFp(Fx, Fy, x, y);
#endif

#ifdef INTERPOLATE_SPLINE_ONE_BODY_Z_WF
  *Fz += InterpolateSplineFp(&G1, z);
#endif

#ifdef ONE_BODY_IMPURITY // add contribution from the w.f. of impurities
  for(j=0; j<Crystal.size; j++) {
    dri[0] = x - Crystal.x[j];
    dri[1] = y - Crystal.y[j];
    dri[2] = z - Crystal.z[j];
    ri2 = FindNearestImage(&dri[0], &dri[1], &dri[2]);
    if(CheckInteractionConditionWF(dri[0], dri[1], dri[2], ri2)) {
      ri = Sqrt(ri2);
      //*Ekin += ImpurityInterpolateE(r);
      Force = ImpurityInterpolateFp(ri)/ri;
      *Fx += Force * dri[0];
      *Fy += Force * dri[1];
      *Fz += Force * dri[2];
    }
  }
#endif

#ifdef CRYSTAL
#ifndef CRYSTAL_SYMMETRIC
#ifdef CHAINS_PINNED_TO_CENTER_OF_MASS
  CaseX(dr[0] = x - W[w].Crystal_x[i]);
  CaseY(dr[1] = y - W[w].Crystal_y[i]);
  CaseZ(dr[2] = z - W[w].Crystal_z[i]);

  r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);

  //if(CheckInteractionConditionWF(dr[0], dr[1], dr[2], r2)) {
  CaseX(*Fx -= 2.*Crystal_dot_Rx_i * dr[0]);
  CaseY(*Fy -= 2.*Crystal_dot_Ry_i * dr[1]);
  CaseZ(*Fz -= 2.*Crystal_dot_Rz_i * dr[2]);
#else
  CaseX(dr[0] = x - Crystal.x[i]);
  CaseY(dr[1] = y - Crystal.y[i]);
  CaseZ(dr[2] = z - Crystal.z[i]);

  r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);

  //if(CheckInteractionConditionWF(dr[0], dr[1], dr[2], r2)) {
    CaseX(*Fx -= 2.*Crystal_dot_Rx_i * dr[0]);
    CaseY(*Fy -= 2.*Crystal_dot_Ry_i * dr[1]);
    CaseZ(*Fz -= 2.*Crystal_dot_Rz_i * dr[2]);
  //}
#endif
#endif
#endif

#ifdef ONE_BODY_TRIAL_SIN
#ifdef BC_3DPBC_CUBE
  *Fx += alpha_latt*kL*sin(2.*kL*x)/(sin(kL*x)*sin(kL*x) + sin(kL*y)*sin(kL*y) + sin(kL*z)*sin(kL*z));
  *Fy += alpha_latt*kL*sin(2.*kL*y)/(sin(kL*x)*sin(kL*x) + sin(kL*y)*sin(kL*y) + sin(kL*z)*sin(kL*z));
  *Fz += alpha_latt*kL*sin(2.*kL*z)/(sin(kL*x)*sin(kL*x) + sin(kL*y)*sin(kL*y) + sin(kL*z)*sin(kL*z));
#endif
#ifdef BC_2DPBC_SQUARE
  *Fx += alpha_latt*kL*sin(2.*kL*x)/(sin(kL*x)*sin(kL*x) + sin(kL*y)*sin(kL*y));
  *Fy += alpha_latt*kL*sin(2.*kL*y)/(sin(kL*x)*sin(kL*x) + sin(kL*y)*sin(kL*y));
#endif
#ifdef BC_1DPBC_Z
  //*Fz += alpha_latt*kL*cos(kL*z)*pow(fabs(sin(kL*z)),alpha_latt-1.)/(beta_latt+pow(fabs(sin(kL*z)),alpha_latt));
  //*Fz += alpha_latt*kL*pow(fabs(sin(kL*z)),alpha_latt)/tan(kL*z)/(beta_latt+pow(fabs(sin(kL*z)),alpha_latt));
  // cotan(x) = tan(M_PI_2 - x)
  *Fz += alpha_latt*kL*pow(fabs(sin(kL*z)),alpha_latt)*tan(_PI_2 - kL*z)/(beta_latt+pow(fabs(sin(kL*z)),alpha_latt));
#endif
#endif

#ifdef ONE_BODY_TRIAL_ZERO_BOUNDARY_CONDITION
  *Fz += beta_z*(PI/L)/tan((PI/L)*z);
#endif

#ifdef ONE_BODY_TRIAL_COS_SERIES 
  z = z - (int)(z); // attention!, z is reduces to [0, 1]
  *Fz -= kL*(2.*alpha_latt*sin(2.*kL*z) + 4.*beta_latt*sin(4.*kL*z) + 6.*gamma_latt*sin(6.*kL*z)) / (1. + alpha_latt*cos(2.*kL*z) + beta_latt*cos(4.*kL*z) + gamma_latt*cos(6.*kL*z)) ;
#endif 

#ifdef ONE_BODY_TRIAL_ZBC
  *Fz += kL/tan(kL*z);
#endif

#ifdef ONE_BODY_TRIAL_ABC
  *Fz += alpha_latt*kL/tan(kL*z);
#endif

#ifdef ONE_BODY_TRIAL_1D_WAVE_COS
  *Fz -= one_body_momentum_k*tan(one_body_momentum_k*z);
#endif

#ifdef ONE_BODY_TRIAL_1D_WAVE_SIN
  if(one_body_momentum_k>0) *Fz += one_body_momentum_k/tan(one_body_momentum_k*z);
#endif

#ifdef ONE_BODY_TRIAL_1D_WAVE_COS_FINITE
  *Fz -= one_body_momentum_k*Srecoil*sin(one_body_momentum_k*z) / (1.+Srecoil*cos(one_body_momentum_k*z));
#endif

#ifdef ONE_BODY_SOLITON
  *Fz += 2.*solitonk*(1.-solitonV2)*tanh(solitonk*(z-Lhalf)) /(-1. + 2.*solitonV2 + cosh(2.*solitonk*(z-Lhalf)));
#endif

#ifdef ONE_BODY_TRIAL_DOUBLE_GAUSSIAN
  *Fz += -2*alpha_z*(z - R0par*tanh(2*alpha_z*R0par*z));
#endif

#ifdef ONE_BODY_TRIAL_LATTICE
  DOUBLE f=0., fx=0., fy=0., fz=0., dr[3]={0.,0.,0.},r2;
  int j;

  for(j=0; j<Crystal.size; j++) {
    CaseX(dr[0] = x - Crystal.x[j]);
    CaseY(dr[1] = y - Crystal.y[j]);
    CaseZ(dr[2] = z - Crystal.z[j]);
    r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
    f += exp(-alpha_latt*r2);
    CaseX(fx += dr[0] * exp(-alpha_latt*r2));
    CaseY(fy += dr[1] * exp(-alpha_latt*r2));
    CaseZ(fz += dr[2] * exp(-alpha_latt*r2));
  }
  f = beta_latt + gamma_latt * f;
  CaseX(*Fx -= fx *2.*alpha_latt*gamma_latt / f);
  CaseY(*Fy -= fy *2.*alpha_latt*gamma_latt / f);
  CaseZ(*Fz -= fz *2.*alpha_latt*gamma_latt / f);
#endif

#ifdef ONE_BODY_TRIAL_TRAP_DOWN_ANTISYMMETRIC_Z
  if(W[w].spin[i]) *Fz += 1./z;
#endif

#ifdef ONE_BODY_TRIAL_TRAP_DOWN_HERMITE2_Z
  if(W[w].spin[i]) *Fz += 2.*z/(z*z-beta);
#endif

}

/********************************* One Body Eloc *****************************/
// Eloc = [-f_x"/f+(f_x'/f)^2] + [-f_y"/f+(f_y'/f)^2] + [-f_z"/f+(f_z'/f)^2]
DOUBLE OneBodyE(DOUBLE x, DOUBLE y, DOUBLE z, int i, int w) {
  DOUBLE E = 0.;

#ifdef ONE_BODY_IMPURITY // add contribution from the w.f. of impurities
  DOUBLE dri[3]={0.,0.,0}, ri2, ri;
  int j;
#endif

#ifdef ONE_BODY_TRIAL_TRAP
  //*Ekin +=  N *(alpha_x + alpha_y + alpha_z); // * 0.5
  CaseX(E += two_alpha_x);
  CaseY(E += two_alpha_y);
  CaseZ(E += two_alpha_z);
#endif

#ifdef ONE_BODY_TRIAL_TRAP_POWER_LAW
  E += 6.*alpha_latt - beta_latt/(x*x + y*y + z*z);
#endif

#ifdef SPECKLES
  E += SpecklesFpp(x,y);
#endif

#ifdef INTERPOLATE_SPLINE_ONE_BODY_Z_WF
  E += InterpolateSplineE(&G1, z);
#endif

#ifdef ONE_BODY_IMPURITY // add contribution from the w.f. of impurities
  for(j=0; j<Crystal.size; j++) {
    dri[0] = x - Crystal.x[j];
    dri[1] = y - Crystal.y[j];
    dri[2] = z - Crystal.z[j];
    ri2 = FindNearestImage(&dri[0], &dri[1], &dri[2]);
    if(CheckInteractionConditionWF(dri[0], dri[1], dri[2], ri2)) {
      ri = Sqrt(ri2);
      E += ImpurityInterpolateE(ri);
    }
  }
#endif

#ifdef CRYSTAL
#ifndef CRYSTAL_SYMMETRIC
#ifdef TRIAL_1D
  E += 2.*Crystal_dot_Rz_i;
#endif
#ifdef TRIAL_2D
  E += 2.*(Crystal_dot_Rx_i + Crystal_dot_Ry_i);
#ifdef CHAINS_PINNED_TO_CENTER_OF_MASS
  E -= 2.*(Crystal_dot_Rx_i + Crystal_dot_Ry_i)/(DOUBLE)Nspin;
#endif
#endif
#ifdef TRIAL_3D
  E += 2.*(Crystal_dot_Rx_i + Crystal_dot_Ry_i + Crystal_dot_Rz_i);
#endif
#endif
#endif

#ifdef ONE_BODY_TRIAL_SIN
#ifdef BC_3DPBC_CUBE
  E += (4.*kL*kL*alpha_latt*(3. - 3.*cos(2.*kL*x) + cos(2.*kL*(x - y)) - 3.*cos(2.*kL*y) + cos(2.*kL*(x + y)) 
        + cos(2.*kL*(x - z)) + cos(2.*kL*(y - z)) - 3.*cos(2.*kL*z) + cos(2.*kL*(x + z)) + cos(2.*kL*(y + z))))
       /(cos(2.*kL*x) + cos(2.*kL*y) + cos(2.*kL*z) - 3.)
       /(cos(2.*kL*x) + cos(2.*kL*y) + cos(2.*kL*z) - 3.);
#endif
#ifdef BC_2DPBC_SQUARE
  E += 4.*kL*kL*alpha_latt*(-2.*cos(2.*kL*x) + cos(2.*kL*(x - y)) + cos(2.*kL*(x + y)) + 4.*sin(kL*y)*sin(kL*y))
       /(cos(2.*kL*x) + cos(2.*kL*y)-2.)
       /(cos(2.*kL*x) + cos(2.*kL*y)-2.);
#endif
#ifdef BC_1DPBC_Z // in order to treat negative values of sin correctly change sin^alpha(x) -> [sin^2(x)]^(alpha/2)
  E += -0.5*(alpha_latt*kL*kL*pow(sin(kL*z)*sin(kL*z),0.5*alpha_latt-1.)*(beta_latt*(-2. + alpha_latt + alpha_latt*cos(2.*kL*z)) 
           - 2.*pow(sin(kL*z)*sin(kL*z),0.5*alpha_latt)))
           / (beta_latt + pow(sin(kL*z)*sin(kL*z),0.5*alpha_latt))
           / (beta_latt + pow(sin(kL*z)*sin(kL*z),0.5*alpha_latt));
#endif
#endif

#ifdef ONE_BODY_TRIAL_ZERO_BOUNDARY_CONDITION
  E += beta_z*(PI/L)*(PI/L)/(sin((PI/L)*z)*sin((PI/L)*z));
#endif

#ifdef ONE_BODY_TRIAL_COS_SERIES
  z = z - (int)(z); // attention!, z is reduces to [0, 1]
  E += 
  kL*kL*(4.*alpha_latt*cos(2.*kL*z) + 16.*beta_latt*cos(4.*kL*z) + 36.*gamma_latt*cos(6.*kL*z)) / (1. + alpha_latt*cos(2.*kL*z) + beta_latt*cos(4.*kL*z) + gamma_latt*cos(6.*kL*z)) 
  + (kL*(2.*alpha_latt*sin(2.*kL*z) + 4.*beta_latt*sin(4.*kL*z) + 6.*gamma_latt*sin(6.*kL*z)) / (1. + alpha_latt*cos(2.*kL*z) + beta_latt*cos(4.*kL*z) + gamma_latt*cos(6.*kL*z)))
   *(kL*(2.*alpha_latt*sin(2.*kL*z) + 4.*beta_latt*sin(4.*kL*z) + 6.*gamma_latt*sin(6.*kL*z)) / (1. + alpha_latt*cos(2.*kL*z) + beta_latt*cos(4.*kL*z) + gamma_latt*cos(6.*kL*z)));
#endif 

#ifdef ONE_BODY_TRIAL_ZBC
  E += kL*kL/(sin(kL*z)*sin(kL*z));
#endif

#ifdef ONE_BODY_TRIAL_ABC
  E += kL*kL*alpha_latt/(sin(kL*z)*sin(kL*z));
#endif

#ifdef ONE_BODY_TRIAL_1D_WAVE_COS
  // -f_z"/f+(f_z'/f)^2
  //E += one_body_momentum_k*one_body_momentum_k*(1.+tan(one_body_momentum_k*z)*tan(one_body_momentum_k*z));
  E += one_body_momentum_k*one_body_momentum_k/(cos(one_body_momentum_k*z)*cos(one_body_momentum_k*z));
#endif

#ifdef ONE_BODY_TRIAL_1D_WAVE_SIN
  // -f_z"/f+(f_z'/f)^2
  if(one_body_momentum_k>0) E += one_body_momentum_k*one_body_momentum_k/(sin(one_body_momentum_k*z)*sin(one_body_momentum_k*z));
#endif

#ifdef ONE_BODY_TRIAL_1D_WAVE_COS_FINITE
  E += one_body_momentum_k*one_body_momentum_k*Srecoil*(Srecoil+cos(one_body_momentum_k*z)) / ((1.+Srecoil*cos(one_body_momentum_k*z))*(1.+Srecoil*cos(one_body_momentum_k*z)));
#endif

#ifdef ONE_BODY_SOLITON
  E += solitonk*solitonk*(1.-solitonV2)
    *(1. - 4.*solitonV2 - 2.*cosh(2.*solitonk*(z-Lhalf)) + cosh(4.*solitonk*(z-Lhalf)))
    /cosh(solitonk*(z-Lhalf))
    /cosh(solitonk*(z-Lhalf))
    /(-1. + 2.*solitonV2 + cosh(2.*solitonk*(z-Lhalf)))
    /(-1. + 2.*solitonV2 + cosh(2.*solitonk*(z-Lhalf)));
#endif

#ifdef ONE_BODY_TRIAL_DOUBLE_GAUSSIAN
  E += 2.*alpha_z*(1. - 2*alpha_z*R0par*R0par/(cosh(2.*alpha_z*R0par*z)*cosh(2.*alpha_z*R0par*z)));
#endif

#ifdef ONE_BODY_TRIAL_LATTICE
  DOUBLE f=0., fx=0., fy=0., fz=0., fxx=0., fyy=0., fzz=0., dr[3]={0.,0.,0.},r2;
  int j;

  for(j=0; j<Crystal.size; j++) {
    CaseX(dr[0] = x - Crystal.x[j]);
    CaseY(dr[1] = y - Crystal.y[j]);
    CaseZ(dr[2] = z - Crystal.z[j]);
    r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
    f += exp(-alpha_latt*r2);
    CaseX(fx += dr[0] * exp(-alpha_latt*r2));
    CaseY(fy += dr[1] * exp(-alpha_latt*r2));
    CaseZ(fz += dr[2] * exp(-alpha_latt*r2));
    CaseX(fxx += exp(-alpha_latt*r2) * (2.*alpha_latt*dr[0]*dr[0]-1.));
    CaseY(fyy += exp(-alpha_latt*r2) * (2.*alpha_latt*dr[1]*dr[1]-1.));
    CaseZ(fzz += exp(-alpha_latt*r2) * (2.*alpha_latt*dr[2]*dr[2]-1.));
  }
  f = beta_latt + gamma_latt * f;
  CaseX(fx *= -2.*alpha_latt*gamma_latt / f);
  CaseY(fy *= -2.*alpha_latt*gamma_latt / f);
  CaseZ(fz *= -2.*alpha_latt*gamma_latt / f);
  CaseX(fxx *= 2.*alpha_latt *gamma_latt / f);
  CaseX(fyy *= 2.*alpha_latt *gamma_latt / f);
  CaseX(fzz *= 2.*alpha_latt *gamma_latt / f);
  // Eloc = [-f_x"/f+(f_x'/f)^2] + [-f_y"/f+(f_y'/f)^2] + [-f_z"/f+(f_z'/f)^2]
  E += -fxx + fx*fx -fyy + fy*fy - fzz +fz*fz;
#endif

#ifdef ONE_BODY_TRIAL_TRAP_DOWN_ANTISYMMETRIC_Z
  if(W[w].spin[i]) E += 1./(z*z);
#endif

#ifdef ONE_BODY_TRIAL_TRAP_DOWN_HERMITE2_Z
  if(W[w].spin[i]) E += 2.*(beta+z*z)/((beta-z*z)*(beta-z*z));
#endif

  return E;
}

/************************ Check Walker Hyperradius Overlapping ***************************/
// for a given three particles (1,2,3) the hyperradius coordinates are constructed as CM of (1,2) pair and the distance between this point and 3
int CheckWalkerHyperradiusOverlapping(const struct Walker W) {
  int i,j,k;
  DOUBLE CM12[3] = {0,0,0};
  DOUBLE dr[3] = {0,0,0};
  DOUBLE rho2,r2;

  for(i=0; i<N; i++) {
    for(j=i+1; j<N; j++) {
      // distance (1-2)
      dr[0] = W.x[i] - W.x[j];
      dr[1] = W.y[i] - W.y[j];
      dr[2] = W.z[i] - W.z[j];
      rho2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);

      // center of mass (1-2) position
      CM12[0] = W.x[i] - 0.5*dr[0];
      CM12[1] = W.y[i] - 0.5*dr[1];
      CM12[2] = W.z[i] - 0.5*dr[2];

      for(k=j+1; k<N; k++) {
        // distance (CM12 - 3)
        dr[0] = W.x[k] - CM12[0];
        dr[1] = W.y[k] - CM12[1];
        dr[2] = W.z[k] - CM12[2];
        r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);

        if(4./3.*r2 + rho2 < R3*R3)
          return 1;
      }
    }
  }

  return 0;
}

int CheckWalkerHyperradiusOverlappingCount(struct Walker *W) {
  int i,j,k;
  DOUBLE CM12[3] = {0,0,0};
  DOUBLE dr[3] = {0,0,0};
  DOUBLE rho2,r2;
  int overlapped = 0;

  for(i=0; i<N; i++) {
    for(j=i+1; j<N; j++) {
      // distance (1-2)
      CaseX(dr[0] = W->x[i] - W->x[j]);
      CaseY(dr[1] = W->y[i] - W->y[j]);
      CaseZ(dr[2] = W->z[i] - W->z[j]);
      rho2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);

      // center of mass (1-2) position
      CaseX(CM12[0] = W->x[i] - 0.5*dr[0]);
      CaseY(CM12[1] = W->y[i] - 0.5*dr[1]);
      CaseZ(CM12[2] = W->z[i] - 0.5*dr[2]);

      for(k=j+1; k<N; k++) {
        // distance (CM12 - 3)
        CaseX(dr[0] = W->x[k] - CM12[0]);
        CaseY(dr[1] = W->y[k] - CM12[1]);
        CaseZ(dr[2] = W->z[k] - CM12[2]);
        r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);

        if(4./3.*r2 + rho2 < R3*R3)
          overlapped++;
      }
    }
  }

  return overlapped;
}

/************************ Check Vector Hyperradius Overlapping ***************************/
// for a given three particles (1,2,3) the hyperradius coordinates are constructed as CM of (1,2) pair and the distance between this point and 3
int CheckVectorHyperradiusOverlapping(DOUBLE **R) {
  int i,j,k;
  DOUBLE CM12[3] = {0,0,0};
  DOUBLE dr[3] = {0,0,0};
  DOUBLE rho2,r2;

  for(i=0; i<N; i++) {
    for(j=i+1; j<N; j++) {
      // distance (1-2)
      dr[0] = R[i][0] - R[j][0];
      dr[1] = R[i][1] - R[j][1];
      dr[2] = R[i][2] - R[j][2];
      rho2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);

      // center of mass (1-2) position
      CM12[0] = R[i][0] - 0.5*dr[0];
      CM12[1] = R[i][1] - 0.5*dr[1];
      CM12[2] = R[i][2] - 0.5*dr[2];

      for(k=j+1; k<N; k++) {
        // distance (CM12 - 3)
        dr[0] = R[k][0] - CM12[0];
        dr[1] = R[k][1] - CM12[1];
        dr[2] = R[k][2] - CM12[2];
        r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);

        if(4./3.*r2 + rho2 < R3*R3)
          return 1;
      }
    }
  }

  return 0;
}

/************************ Adjust Center Of Mass Walker ***********************************/
void AdjustCenterOfMassWalker(struct Walker *Walker) {
  int i;
  DOUBLE CM[3],M;

#ifdef CENTER_OF_MASS_IS_NOT_MOVED // put the center of mass to zero
  CM[0] = CM[1] = CM[2] = 0.;
  for(i=0; i<N; i++) {
    CaseX(CM[0] += Walker->x[i]);
    CaseY(CM[1] += Walker->y[i]);
    CaseZ(CM[2] += Walker->z[i]);
  }
  M = (double) N;
#ifdef BC_ABSENT
  CM[0] = CM[0]/M;
  CM[1] = CM[1]/M;
  CM[2] = CM[2]/M;
#else
  CM[0] = CM[0]/M - Lhalf;
  CM[1] = CM[1]/M - Lhalf;
  CM[2] = CM[2]/M - Lhalf;
#endif

  for(i=0; i<N; i++) {
    CaseX(Walker->x[i] -= CM[0]);
    CaseY(Walker->y[i] -= CM[1]);
    CaseZ(Walker->z[i] -= CM[2]);
  }
#endif

#ifdef CENTER_OF_MASS_Z_IS_NOT_MOVED // put the center of mass to zero
  CM[0] = CM[1] = CM[2] = 0.;
  for(i=0; i<N; i++) {
    CaseZ(CM[2] += Walker->z[i]);
  }
  M = (double) N;
#ifdef BC_ABSENT
  CM[2] = CM[2]/M;
#else
  CM[2] = CM[2]/M - Lhalf;
#endif

  for(i=0; i<N; i++) {
    CaseZ(Walker->z[i] -= CM[2]);
  }
#endif

#ifdef CENTER_OF_MASS_IS_NOT_MOVED_TWO_PARTS // put the center of mass to zero for the first half and the second half of the system
  CM[0] = CM[1] = CM[2] = 0.;
  for(i=0; i<N/2; i++) {
    CaseX(CM[0] += Walker->x[i]);
    CaseY(CM[1] += Walker->y[i]);
    CaseZ(CM[2] += Walker->z[i]);
  }
  CM[0] /= (N/2);
  CM[1] /= (N/2);
  CM[2] /= (N/2);
#ifdef BC_ABSENT
  for(i=0; i<N/2; i++) {
    CaseX(Walker->x[i] -= CM[0]);
    CaseY(Walker->y[i] -= CM[1]);
    CaseZ(Walker->z[i] -= CM[2]);
  }
#endif

  CM[0] = CM[1] = CM[2] = 0.;
  for(i=N/2; i<N; i++) {
    CaseX(CM[0] += Walker->x[i]);
    CaseY(CM[1] += Walker->y[i]);
    CaseZ(CM[2] += Walker->z[i]);
  }
  CM[0] /= (N/2);
  CM[1] /= (N/2);
  CM[2] /= (N/2);
#ifdef BC_ABSENT
  for(i=N/2; i<N; i++) {
    CaseX(Walker->x[i] -= (CM[0]-CMseparation));
    CaseY(Walker->y[i] -= (CM[1]-CMseparation));
    CaseZ(Walker->z[i] -= (CM[2]-CMseparation));
  }
#endif
#endif

  //ReduceWalkerToTheBox(&W[w]);

#ifdef FIRST_PARTICLE_IS_NOT_MOVED
  for(i=N-1; i>=0; i--) {
    CaseX(Walker->x[i] -= Walker->x[0]);
    CaseY(Walker->y[i] -= Walker->y[0]);
    CaseZ(Walker->z[i] -= Walker->z[0]);
  }
#endif
}