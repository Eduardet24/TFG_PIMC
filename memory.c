/*memory.c*/

#include <malloc.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "memory.h"
#include "trial.h"
#include "utils.h"
#include "main.h"
#include "compatab.h"

#define ON 1
#define OFF 0

/*************************** Allocate Grid **********************************/
void AllocateWFGrid(struct Grid *G, unsigned long int size) {

  if(verbosity) Message("\nMemory allocation : allocating grid of size %li...", size);

  G->x   = (DOUBLE*) Calloc("G->x   ", size, sizeof(DOUBLE));
  G->fp  = (DOUBLE*) Calloc("G->fp  ", size, sizeof(DOUBLE));
  G->lnf = (DOUBLE*) Calloc("G->lnf",  size, sizeof(DOUBLE));
  G->E   = (DOUBLE*) Calloc("G->E",    size, sizeof(DOUBLE));
#ifdef INTERACTION_WITH_DAMPING
  G->ReV = (DOUBLE*) Calloc("G->ReV   ", size, sizeof(DOUBLE));
  G->ImV = (DOUBLE*) Calloc("G->ImV  ", size, sizeof(DOUBLE));
#endif
  G->f   = (DOUBLE*) Calloc("G->f  ",  size+1, sizeof(DOUBLE));
  G->fpp = (DOUBLE*) Calloc("G->fpp",  size+1, sizeof(DOUBLE));

  G->size = size;

  if(verbosity) Message(" done\n");
}

/*************************** Allocate Walkers *******************************/
void AllocateSingleWalker(struct Walker* W) {
  int i,j;

  // initialize first scalar elements of the structure. Arrays are initialized later
  W->weight = 1.;
  W->times_measured = 0;
  W->PD_position = 0;
  W->Sk_position = 0;
  W->Sk_tau_position = 0;
#ifdef MEASURE_CORRELATION_FUNCTIONS_IN_THE_MIDDLE_OF_DRIFT
  W->status_end_of_step_initialized = OFF;
#endif
  W->status_sum_of_two_walkers_U = 0;
  W->status_sum_of_two_walkers_E = 0;

  // data which is copied by CopyWalker
  // 1)
#ifndef NPARTICLES // i.e. allocate dynamically
  W->x = (DOUBLE*) Calloc("W->x", N, sizeof(DOUBLE));
  W->y = (DOUBLE*) Calloc("W->y", N, sizeof(DOUBLE));
  W->z = (DOUBLE*) Calloc("W->z", N, sizeof(DOUBLE));
  W->spin = (int*) Calloc("W->spin", N, sizeof(int));
#ifdef SPINFULL
    W->psiT_sigma_inv = (DOUBLE*) Calloc("W->psiT_sigma_inv", N, sizeof(DOUBLE));
#endif
#ifdef CHAINS_PINNED_TO_CENTER_OF_MASS
    W->Crystal_x = (DOUBLE*) Calloc("W->Crystal_x", N, sizeof(DOUBLE));
    W->Crystal_y = (DOUBLE*) Calloc("W->Crystal_y", N, sizeof(DOUBLE));
    W->Crystal_z = (DOUBLE*) Calloc("W->Crystal_z", N, sizeof(DOUBLE));
#endif
#endif

#ifdef CRYSTAL_SYMMETRIC
    W->Mo = (DOUBLE*) Calloc("Mo", N, sizeof(DOUBLE));
#endif

  if(MC == DIFFUSION) { // pure estimators
    // 2)
#ifdef MEASURE_CORRELATION_FUNCTIONS_IN_THE_MIDDLE_OF_DRIFT
    if(SmartMC == DMC_QUADRATIC) { // distributions are measured in the middle of the step, so last coordinates should be stored
      ArrayCalloc2D(W->R_dmc_end_of_step, "dR_drift ", i, N, 3, DOUBLE, "DOUBLE");
    }
#endif

    // 3)
    if(measure_energy) {
     W->Epot_pure = (DOUBLE*) Calloc("W->Epot_pure", grid_pure_block, sizeof(DOUBLE));
    }

    // 4)
    if(measure_PairDistr) {
      ArrayCalloc2D(W->PD, "W->PD ", i, grid_pure_block, gridPD, int, "int");
#ifdef SPINFULL
      ArrayCalloc3D(W->PDSpin, "W->PDSpin ", i, grid_pure_block, j, Nspin, gridPD, int, "int");
#endif
      if(measure_PairDistr_R2) W->PD_R2pure = (DOUBLE*) Calloc("W->PD_R2pure", grid_pure_block, sizeof(DOUBLE));
    }

    // 5)
    if(measure_RadDistr) {
      W->RD_wait = 1;
      ArrayCalloc2D(W->RD, "W->RD ", i, grid_pure_block, gridRD, int, "int");
      ArrayCalloc2D(W->RDz, "W->RDz ", i, grid_pure_block, gridRD, int, "int");
      W->R2pure = (DOUBLE*) Calloc("W->R2pure", grid_pure_block, sizeof(DOUBLE));
      W->Z2pure = (DOUBLE*) Calloc("W->Z2pure", grid_pure_block, sizeof(DOUBLE));
    }

    // 6)
    if(measure_SD) {
      ArrayCalloc2D(W->rreal, "W->rreal ", i, N, 3, DOUBLE, "DOUBLE");
      ArrayCalloc2D(W->CM, "W->CM ", i, 3, SD.size, DOUBLE, "DOUBLE");
    }

    // 7)
    if(measure_Sk) {
      ArrayCalloc2D(W->Sk, "W->Sk ", i, grid_pure_block, gridSk, DOUBLE, "DOUBLE");
      ArrayCalloc2D(W->Sk_N, "W->Sk_N ", i, grid_pure_block, gridSk, unsigned long int, "unsigned long int");
    }

    // 8)
    if(measure_FormFactor) {
      ArrayCalloc2D(W->rhoRe, "W->rhoRe ", i, gridSKT_t, gridSKT_k, DOUBLE, "DOUBLE");
      ArrayCalloc2D(W->rhoIm, "W->rhoIm ", i, gridSKT_t, gridSKT_k, DOUBLE, "DOUBLE");
      ArrayCalloc2D(W->psi, "W->rhoIm ", i, gridSKT_t, gridSKT_k, DOUBLE, "DOUBLE");
    }

    // 9)
    if(measure_pure_coordinates) {
      W->x_pure = (DOUBLE*) Calloc("W->x_pure", N, sizeof(DOUBLE));
      W->y_pure = (DOUBLE*) Calloc("W->y_pure", N, sizeof(DOUBLE));
      W->z_pure = (DOUBLE*) Calloc("W->z_pure", N, sizeof(DOUBLE));
    }

    // 10)
    if(measure_Lind) {
      W->LindF = (DOUBLE*) Calloc("W->LindF", grid_pure_block, sizeof(DOUBLE));
      W->LindN = (int*) Calloc("W->LindN", grid_pure_block, sizeof(int));
      W->Lind_position = 0;
    }

    // 11)
#ifdef VIRTUAL_WALKERS
    if(measure_OBDM || measure_TBDM) {
      W->OBDMweight = (DOUBLE*) Calloc("W->OBDMweight", grid_pure_block, sizeof(DOUBLE));
    }
#endif

    // 12)
    if(measure_g3) {
      ArrayCalloc2D(W->HR, "W->HR ", i, grid_pure_block, gridg3, DOUBLE, "DOUBLE");
    }

    // 13)
    if(measure_OP) {
      W->Skx = (DOUBLE*) Calloc("W->Skx", grid_pure_block , sizeof(DOUBLE));
      W->Sky = (DOUBLE*) Calloc("W->Sky", grid_pure_block , sizeof(DOUBLE));
    }
  } // DIFFUSION

  // data which is NOT copied by CopyWalker
  ArrayCalloc2D(W->F, "F", i, N, 6, DOUBLE, "DOUBLE");

  if(MC == DIFFUSION || (MC == VARIATIONAL && SmartMC == VMC_MOVE_DRIFT_ALL)) { // in case of DMC allocate arrays for jumps
    ArrayCalloc2D(W->R, "R", i, N, 3, DOUBLE, "DOUBLE");
    ArrayCalloc2D(W->Rp, "Rp", i, N, 3, DOUBLE, "DOUBLE");
    ArrayCalloc2D(W->Rpp, "Rpp", i, N, 3, DOUBLE, "DOUBLE");
    ArrayCalloc2D(W->Rppp, "Rppp", i, N, 3, DOUBLE, "DOUBLE");
    ArrayCalloc2D(W->Fp, "Fp", i, N, 3, DOUBLE, "DOUBLE");
    ArrayCalloc2D(W->Fpp, "Fpp", i, N, 3, DOUBLE, "DOUBLE");
    if(MC == DIFFUSION && (SmartMC == DMC_LINEAR_METROPOLIS || SmartMC == DMC_QUADRATIC_METROPOLIS)) {
      ArrayCalloc2D(W->dR_drift, "dR_drift ", i, N, 3, DOUBLE, "DOUBLE");
      ArrayCalloc2D(W->dR_drift_old, "dR_drift_old ", i, N, 3, DOUBLE, "DOUBLE");
    }
    if(MC == DIFFUSION && (SmartMC == DMC_LINEAR_METROPOLIS || SmartMC == DMC_QUADRATIC_METROPOLIS || SmartMC == DMC_PSEUDOPOTENTIAL)) {
      ArrayCalloc2D(W->dR, "dR", i, N, 3, DOUBLE, "DOUBLE");
    }
    //if(measure_spectral_weight) {
    //}
  } // DIFFUSION

#ifdef SCALABLE
  W->c_Nlocal = (int*) Calloc("W->c_Nlocal", Ncells, sizeof(int));
  ArrayCalloc2D(W->c_index_local, "W->c_index_local ", i, Ncells, N, int, "int");
  W->c_Nall = (int*) Calloc("W->c_Nall", Ncells, sizeof(int));
  ArrayCalloc2D(W->c_index_all, "W->c_index_all ", i, Ncells, N, int, "int");
  W->c_Npositive = (int*) Calloc("W->c_Npositive", Ncells, sizeof(int));
  ArrayCalloc2D(W->c_index_positive, "W->c_index_positive ", i, Ncells, N, int, "int");
#endif
}

/*************************** Copy Walker ************************************/
void CopyWalker(struct Walker *out, const struct Walker *in) {
  int i, j, k;
#ifdef SPINFULL
  int js, ks;
#endif

  out->E = in->E;
  out->EFF = in->EFF;
  out->Ekin = in->Ekin;
  out->Epot = in->Epot;
  out->Eint = in->Eint;
  out->Eext = in->Eext;
  out->Eold = in->Eold;
  out->U = in->U;
  out->status = in->status;
  out->weight = in->weight;
  out->r2 = in->r2;
  out->z2 = in->z2;
  out->r2pd = in->r2pd;
  out->r2old  = in->r2old;
  out->z2old  = in->z2old;
#ifdef INTERACTION_WITH_DAMPING
  out->Edamping = in->Edamping;
#endif
  out->times_measured = out->times_measured

  // 1)
  CaseX(ArrayCopy1D(out->x, in->x, i, N));
  CaseY(ArrayCopy1D(out->y, in->y, i, N));
  CaseZ(ArrayCopy1D(out->z, in->z, i, N));
#ifdef SPINFULL
  ArrayCopy1D(out->spin, in->spin, i, N);
#endif
#ifdef SPINFULL_TUNNELING
  ArrayCopy1D(out->psiT_sigma_inv, in->psiT_sigma_inv, i, N);
#endif
#ifdef CHAINS_PINNED_TO_CENTER_OF_MASS
  CaseX(ArrayCopy1D(out->Crystal_x, in->Crystal_x, i, N));
  CaseY(ArrayCopy1D(out->Crystal_y, in->Crystal_y, i, N));
  CaseZ(ArrayCopy1D(out->Crystal_z, in->Crystal_z, i, N));
#endif

#ifdef MEASURE_CORRELATION_FUNCTIONS_IN_THE_MIDDLE_OF_DRIFT
  out->status_end_of_step_initialized = in->status_end_of_step_initialized;

  if(MC == DIFFUSION && SmartMC == DMC_QUADRATIC) {
    ArrayCopy2D(out->R_dmc_end_of_step, in->R_dmc_end_of_step, i, N, j, 3);
  }
#endif

  if(MC == DIFFUSION && SmartMC == DMC_QUADRATIC_METROPOLIS) {
    ArrayCopy2D(out->dR_drift_old, in->dR_drift_old, i, N, j, 3);
  }

  if(MC == DIFFUSION && SmartMC == DMC_LINEAR_METROPOLIS) {
    ArrayCopy2D(out->dR_drift_old, in->dR_drift_old, i, N, j, 3);
  }

  // 2)
  //if(MC == DIFFUSION && (SmartMC == DMC_LINEAR_METROPOLIS || SmartMC == DMC_LINEAR_METROPOLIS_MIDDLE)) {
  //  ArrayCopy2D(out->dR_drift_old, in.dR_drift_old, i, Nup, j, 6);
   // out->status_end_of_step_initialized = in.status_end_of_step_initialized;
  //}

  if(MC == DIFFUSION) {
    // 2)
    if(SmartMC == DMC_QUADRATIC_METROPOLIS) { // distributions are measured in the middle of the step, so last coordinates should be stored
#ifdef MEASURE_CORRELATION_FUNCTIONS_IN_THE_MIDDLE_OF_DRIFT
      out->status_end_of_step_initialized = in->status_end_of_step_initialized;
      ArrayCopy2D(out->R_dmc_end_of_step, in->R_dmc_end_of_step, i, N, j, 3);
#endif
    }

    // 3)
    if(measure_energy) {
     ArrayCopy1D(out->Epot_pure, in->Epot_pure, i, grid_pure_block);
    }

    // 4)
    if(measure_PairDistr) {
      out->PD_position = in->PD_position;
      out->g3_store = in->g3_store;
      ArrayCopy2D(out->PD, in->PD, i, grid_pure_block, j, gridPD);
#ifdef SPINFULL
      ArrayCopy3D(out->PDSpin, in->PDSpin, i, grid_pure_block, js, Nspin, ks, gridPD);
#endif
      if(measure_PairDistr_R2) ArrayCopy1D(out->PD_R2pure, in->PD_R2pure, i, grid_pure_block);
    }

    // 5)
    if(measure_RadDistr) {
      out->RD_position = in->RD_position;
      out->RD_wait = in->RD_wait;
      ArrayCopy2D(out->RD, in->RD, i, grid_pure_block, j, gridRD);
      ArrayCopy2D(out->RDz, in->RDz, i, grid_pure_block, j, gridRD);
      if(measure_R2) {
        ArrayCopy1D(out->R2pure, in->R2pure, i, grid_pure_block);
        ArrayCopy1D(out->Z2pure, in->Z2pure, i, grid_pure_block);
      }
    }

    // 6)
    if(measure_SD) {
      ArrayCopy2D(out->CM, in->CM, i, 3, j, SD.size);
      ArrayCopy2D(out->rreal, in->rreal, i, N, j, 3);
    }

    // 7)
    if(measure_Sk) {
      out->Sk_count = in->Sk_count;
      out->Sk_position = in->Sk_position;
      ArrayCopy2D(out->Sk, in->Sk, i, grid_pure_block, j, gridSk);
      ArrayCopy2D(out->Sk_N, in->Sk_N, i, grid_pure_block, j, gridSk);
    }

    // 8)
    if(measure_FormFactor) {
      out->Sk_tau_position = in->Sk_tau_position;
      ArrayCopy2D(out->rhoRe, in->rhoRe, i, gridSKT_t, j, gridSKT_k);
      ArrayCopy2D(out->rhoIm, in->rhoIm, i, gridSKT_t, j, gridSKT_k);
      ArrayCopy2D(out->psi, in->psi, i, gridSKT_t, j, gridSKT_k);
    }

    // 9)
    if(measure_pure_coordinates) {
      CaseX(ArrayCopy1D(out->x_pure, in->x_pure, i, N));
      CaseY(ArrayCopy1D(out->y_pure, in->y_pure, i, N));
      CaseZ(ArrayCopy1D(out->z_pure, in->z_pure, i, N));
    }

    // 10)
    if(measure_Lind) {
      out->Lind_position = in->Lind_position;
      ArrayCopy1D(out->LindF, in->LindF, i, grid_pure_block);
      ArrayCopy1D(out->LindN, in->LindN, i, grid_pure_block);
    }

    // 11)
#ifdef VIRTUAL_WALKERS
    if(measure_OBDM || measure_TBDM) {
      out->OBDMpureB = in->OBDMpureB;
      out->OBDMpureF = in->OBDMpureF;
      out->OBDMpure_r = in->OBDMpure_r;
      out->OBDM_position = in->OBDM_position;
      ArrayCopy1D(out->OBDMweight, in->OBDMweight, i, in->OBDM_position);
    }
#endif

    // 12)
    //if(measure_g3) {
    //}

    // 13)
    if(measure_OP) {
      out->Skxy_count = in->Skxy_count;
      out->Sk_position = in->Sk_position;
      ArrayCopy1D(out->Skx, in->Skx, i, grid_pure_block);
      ArrayCopy1D(out->Sky, in->Sky, i, grid_pure_block);
    }
  }

#ifdef SCALABLE
  ArrayCopy1D(out->c_Nall, in->c_Nall, i, Ncells);
  ArrayCopy1D(out->c_Nlocal, in->c_Nlocal, i, Ncells);
  ArrayCopy1D(out->c_Npositive, in->c_Npositive, i, Ncells);
  for(i=0; i<Ncells; i++) {
    for(j=0; j<in->c_Nall[i]; j++) out->c_index_all[i][j] = in->c_index_all[i][j];
    for(j=0; j<in->c_Nlocal[i]; j++) out->c_index_local[i][j] = in->c_index_local[i][j];
    for(j=0; j<in->c_Npositive[i]; j++) out->c_index_positive[i][j] = in->c_index_positive[i][j];
  }
#endif
}

/*************************** Allocate Walkers *******************************/
void AllocateWalkers(void) {
  int w, i;
#ifdef SPINFULL
  int j;
#endif

  if(verbosity) Message("Memory allocation : allocating walkers ...");

  walkers_storage = (int*) Calloc("Walkers storage", NwalkersMax, sizeof(int));
  dead_walkers_storage = (int*) Calloc("Dead walkers storage", NwalkersMax, sizeof(int));
  branching_weight = (DOUBLE*) Calloc("branching weight", NwalkersMax, sizeof(DOUBLE));
  branching_xi = (DOUBLE*) Calloc("branching xi", NwalkersMax, sizeof(DOUBLE));
  branching_multiplicity = (int*) Calloc("branching multiplicity", NwalkersMax, sizeof(int));

#ifdef VIRTUAL_WALKERS
  if((measure_OBDM || measure_TBDM) && MC == DIFFUSION) {
    Warning("  Increasing maximal number of walkers by number of virtual walkers, %i -> %i\n", NwalkersMax, NwalkersMax+Npop_virtual);
    NwalkersMax = NwalkersMax + Npop_virtual;
  }
#else
  Npop_virtual = 0;
#endif

  W = (struct Walker*) Calloc("W", NwalkersMax, sizeof(struct Walker));
  for(w=0; w<NwalkersMax; w++) {
    AllocateSingleWalker(&W[w]);
    W[w].w = w;
  }

  Wp = (struct Walker*) Calloc("Wp", NwalkersMax, sizeof(struct Walker));
  for(w=0; w<NwalkersMax; w++) {
    AllocateSingleWalker(&Wp[w]);
    Wp[w].w = w;
  }

  for(w=0; w<NwalkersMax; w++) { // initialize spin
    for(i=0; i<N; i++) {
#ifdef SPINLESS
      W[w].spin[i] = 0;
      Wp[w].spin[i] = 0;
#else
      W[w].spin[i] = (int) (i-Nspin*(i / Nspin));
      Wp[w].spin[i] = (int) (i-Nspin*(i / Nspin));
#endif
    }
  }

  if(verbosity) Message(" done\n");
}

/*************************** Allocate Values ********************************/
void AllocateGrids(void) {
  int i;

  if(verbosity) Message("Memory allocation : allocating values ...");

  OBDM.f = (DOUBLE*) Calloc("OBDM->f", gridOBDM, sizeof(DOUBLE));
  OBDM.N = (DOUBLE*) Calloc("OBDM->N", gridOBDM, sizeof(DOUBLE));
  OBDM.size = gridOBDM;
  if(OBDM.max <= 0) OBDM.max = Lcutoff_pot;
  OBDM.min = 0;
  OBDM.step = OBDM.max/gridOBDM;
  OBDM.kmin *= PI*n;
  OBDM.kmax *= PI*n;
  if(OBDM.Nksize) {
    OBDM.k  = (DOUBLE*) Calloc("OBDM->k", OBDM.Nksize, sizeof(DOUBLE));
    OBDM.Nk = (DOUBLE*) Calloc("OBDM->Nk", OBDM.Nksize, sizeof(DOUBLE));

#ifndef BC_ABSENT // PBC
    //if(OBDM.kmin == 0) OBDM.kmin = 2.*PI/L // zero momentum is a valid value for the momentum distribution
    Message("  OBDM.kmin = %" LE " \n", OBDM.kmin);
    if(OBDM.kmax == 0) OBDM.kmax = OBDM.kmin + (N-1)*2.*PI/L;
    Message("  OBDM.kmax = %" LE " \n", OBDM.kmax);
    for(i=0; i<OBDM.Nksize; i++) { // initialize the values of the momenta (1D)
      OBDM.k[i] = OBDM.kmin + (OBDM.kmax-OBDM.kmin)/(N-1) * i;
    }
#else // trap
    //OBDM.kmin = NkMaxTrap/(DOUBLE) (OBDM.Nksize-1);
    OBDM.kmin = NkMaxTrap / (DOUBLE)(1+2*(OBDM.Nksize-1)); // k_0 = dk/2
    //OBDM.kmin = NkMaxTrap / (DOUBLE)(OBDM.Nksize); // k_0 = dk
    Message("  OBDM.kmin = %" LE " \n", OBDM.kmin);
    OBDM.kmax = NkMaxTrap;
    Message("  OBDM.kmax = %" LE " \n", OBDM.kmax);
    for(i=0; i<OBDM.Nksize; i++) { // initialize the values of the momenta (1D)
      //OBDM.k[i] = OBDM.kmin + (OBDM.kmax-OBDM.kmin)/(OBDM.Nksize-1) * i;
      OBDM.k[i] = OBDM.kmin*(DOUBLE)(1+2*i);
      //OBDM.k[i] = OBDM.kmin*(DOUBLE)(i+1);
    }
#endif
  }

  if(measure_TBDM) {
    TBDM.f = (DOUBLE*) Calloc("TBDM->f", gridOBDM, sizeof(DOUBLE));
    TBDM.N = (DOUBLE*) Calloc("TBDM->N", gridOBDM, sizeof(DOUBLE));
    TBDM.size = gridOBDM;
    if(TBDM.max == 0) TBDM.max = OBDM.max;
    TBDM.min = 0;
    TBDM.step = TBDM.max/gridOBDM;

    ArrayCalloc2D(TBDM_MATRIX.f, "TBDM_MATRIX->f ", i, gridOBDM_MATRIX, gridOBDM_MATRIX, DOUBLE, "DOUBLE");
    ArrayCalloc2D(TBDM_MATRIX.N, "TBDM_MATRIX->N ", i, gridOBDM_MATRIX, gridOBDM_MATRIX, int, "int");
    TBDM_MATRIX.times_measured = 0;
    TBDM_MATRIX.size = gridOBDM_MATRIX;
    if(TBDM_MATRIX.max == 0) TBDM_MATRIX.max = 0.5*Lcutoff_pot;
    TBDM_MATRIX.min = TBDM_MATRIX.Norma = 0;
#ifdef BC_ABSENT
    TBDM_MATRIX.step = 2*TBDM_MATRIX.max/(gridOBDM_MATRIX-1.);
#else
    TBDM_MATRIX.step = TBDM_MATRIX.max/(gridOBDM_MATRIX);
#endif

#ifdef OBDM_FERMIONS
    TBDMfermi.f = (DOUBLE*) Calloc("TBDMfermi->f", gridOBDM, sizeof(DOUBLE));

    ArrayCalloc2D(TBDMfermi_MATRIX.f, "TBDMfermi_MATRIX->f ", i, gridOBDM_MATRIX, gridOBDM_MATRIX, DOUBLE, "DOUBLE");
    ArrayCalloc2D(TBDMfermi_MATRIX.N, "TBDMfermi_MATRIX->N ", i, gridOBDM_MATRIX, gridOBDM_MATRIX, int, "int");
    TBDMfermi_MATRIX.size = gridOBDM_MATRIX;
    if(TBDMfermi_MATRIX.max == 0) TBDMfermi_MATRIX.max = Lcutoff_pot;
    TBDMfermi_MATRIX.min = TBDMfermi_MATRIX.Norma = 0;
#ifdef BC_ABSENT
    TBDMfermi_MATRIX.step = 2*TBDMfermi_MATRIX.max/(gridOBDM_MATRIX-1);
#else
    TBDMfermi_MATRIX.step = TBDMfermi_MATRIX.max/(gridOBDM_MATRIX);
#endif
#endif
  }

  PD.N = (int*) Calloc("PD", gridPD, sizeof(int));
  PD.size = gridPD;
  PD.width = 0.1;
  PD.min = 0;
  if(PD.max == 0) PD.max = Lcutoff_pot;
  PD.r2 = 0.;
  PD.step = PD.max/gridPD;
  PD.times_measured = 0;

 #ifdef SPINFULL
   ArrayCalloc2D(PD.PDSpin, "PD->PDSpin ", i, Nspin, gridPD, int, "int");  
 #endif 

  HR.N = (int*) Calloc("HR", gridg3, sizeof(int));
  HR.f = (DOUBLE*) Calloc("HR", gridg3, sizeof(DOUBLE));
  HR.size = gridg3;
  HR.width = 0.1;
  HR.min = 0;
  if(HR.max == 0) HR.max = Lcutoff_pot;
  HR.r2 = 0.;
  HR.step = HR.max/gridg3;
  HR.times_measured = 0;

  if(measure_OBDM_MATRIX == ON) {
    ArrayCalloc2D(OBDM_MATRIX.f, "OBDM_MATRIX->f ", i, gridOBDM_MATRIX, gridOBDM_MATRIX, DOUBLE, "DOUBLE");
    ArrayCalloc2D(OBDM_MATRIX.N, "OBDM_MATRIX->N ", i, gridOBDM_MATRIX, gridOBDM_MATRIX, int, "int");
#ifdef OBDM_FERMIONS
    ArrayCalloc2D(OBDMfermi_MATRIX.f, "OBDMfermi_MATRIX->f ", i, gridOBDM_MATRIX, gridOBDM_MATRIX, DOUBLE, "DOUBLE");
    ArrayCalloc2D(OBDMfermi_MATRIX.N, "OBDMfermi_MATRIX->N ", i, gridOBDM_MATRIX, gridOBDM_MATRIX, int, "int");
#endif

    ArrayCalloc2D(Nk_MATRIX.f, "Nk_MATRIX->f ", i, gridSk, gridSk, DOUBLE, "DOUBLE");
    ArrayCalloc2D(Nk_MATRIX.N, "Nk_MATRIX->N ", i, gridSk, gridSk, int, "int");
#ifdef OBDM_FERMIONS
    ArrayCalloc2D(Nkfermi_MATRIX.f, "Nkfermi_MATRIX->f ", i, gridSk, gridSk, DOUBLE, "DOUBLE");
    ArrayCalloc2D(Nkfermi_MATRIX.N, "Nkfermi_MATRIX->N ", i, gridSk, gridSk, int, "int");
#endif
  }

  if(MC == DIFFUSION && measure_FormFactor) {
    ArrayCalloc2D(SKT.f, "SKT->f ", i, gridSKT_t, gridSKT_k, DOUBLE, "DOUBLE");
    ArrayCalloc2D(SKTre.f, "SKTre->f ", i, gridSKT_t, gridSKT_k, DOUBLE, "DOUBLE");
    ArrayCalloc2D(SKTim.f, "SKTim->f ", i, gridSKT_t, gridSKT_k, DOUBLE, "DOUBLE");
    ArrayCalloc2D(SKT.N, "SKT->N ", i, gridSKT_t, gridSKT_k, int, "int");
    ArrayCalloc2D(PhiTau.f, "PhiTau->f ", i, gridSKT_t, gridSKT_k, DOUBLE, "DOUBLE");
  }

  OBDM_MATRIX.size = gridOBDM_MATRIX;
#ifdef SPECKLES
  if(OBDM_MATRIX.max == 0) OBDM_MATRIX.max = L/sqrt(N); // reduce to one period
#else
  if(OBDM_MATRIX.max == 0) OBDM_MATRIX.max = L;
#endif
  OBDM_MATRIX.Norma = 0;
#ifdef BC_ABSENT
  OBDM_MATRIX.min = - OBDM_MATRIX.max;
  OBDM_MATRIX.step = (OBDM_MATRIX.max-OBDM_MATRIX.min)/(OBDM_MATRIX.size-1.); // in [-L;L]
#else
  OBDM_MATRIX.min = 0.;
  OBDM_MATRIX.step = (OBDM_MATRIX.max-OBDM_MATRIX.min)/(OBDM_MATRIX.size); // in [0;L)
#endif

#ifdef BC_ABSENT
  OBDMtrap.f = (DOUBLE*) Calloc("OBDMtrap->f", gridOBDM, sizeof(DOUBLE));
  OBDMtrap.N = (DOUBLE*) Calloc("OBDMtrap->N", gridOBDM, sizeof(DOUBLE));
  OBDMtrap.size = gridOBDM;
  if(OBDMtrap.max == 0) OBDMtrap.max = Lcutoff_pot;
  OBDMtrap.min = 0;
  OBDMtrap.step = OBDMtrap.max/gridOBDM;
#endif

  if(measure_PairDistrMATRIX) {
    ArrayCalloc2D(PD_MATRIX.f, "PD_MATRIX->f ", i, gridPD_MATRIX_x, gridPD_MATRIX_y, DOUBLE, "DOUBLE");
    ArrayCalloc2D(PD_MATRIX.N, "PD_MATRIX->N ", i, gridPD_MATRIX_x, gridPD_MATRIX_y, int, "int");
    if(PD_MATRIX_x == 0) PD_MATRIX_x = L_half_x;
    if(PD_MATRIX_y == 0) PD_MATRIX_y = L_half_y;
    PD_MATRIX.step_x = PD_MATRIX_x/(gridPD_MATRIX_x-1.);
    PD_MATRIX.step_y = PD_MATRIX_y/(gridPD_MATRIX_y-1.);
  }

  if(measure_RadDistrMATRIX) {
    ArrayCalloc2D(RD_MATRIX.f, "RD_MATRIX->f ", i, gridRDx, gridRDy, DOUBLE, "DOUBLE");
    ArrayCalloc2D(RD_MATRIX.N, "RD_MATRIX->N ", i, gridRDx, gridRDy, int, "int");
    if(RD_MATRIX.maxx<1e-8) {
      Warning(" Setting maximal x size of RD_MATRIX to Lx\n");
      RD_MATRIX.maxx = Lx;
    }
    if(RD_MATRIX.maxy<1e-8) {
      Warning(" Setting maximal x size of RD_MATRIX to Ly\n");
      RD_MATRIX.maxy = Ly;
    }
    RD_MATRIX.step_x = RD_MATRIX.maxx/(DOUBLE) gridRDx;
    RD_MATRIX.step_y = RD_MATRIX.maxy/(DOUBLE) gridRDy;
  }

  if(MC == DIFFUSION && measure_PairDistr) {
    PD_pure.N = (int*) Calloc("PD_pure", gridPD, sizeof(int));
    PD_pure.size = gridPD;
    PD_pure.width = 0.1;
    PD_pure.max = PD.max;
    if(PD_pure.max == 0) PD_pure.max = Lcutoff_pot;
    PD_pure.step = PD_pure.max/gridPD;
#ifdef SPINFULL
    ArrayCalloc2D(PD_pure.PDSpin, "PD_pure->PDSpin ", i, Nspin, gridPD, int, "int");
#endif 
  }

  if(MC == DIFFUSION) { 
    Epot_pure.x = (DOUBLE*) Calloc("Epot_pure", grid_pure_block, sizeof(DOUBLE));
    Epot_pure.f = (DOUBLE*) Calloc("Epot_pure", grid_pure_block, sizeof(DOUBLE));
    Epot_pure.size = grid_pure_block;
    Epot_pure.min = 0;
    Epot_pure.step = dt * (DOUBLE) Nmeasure;
    Epot_pure.max = Epot_pure.step*(Epot_pure.size-1.);
    for(i=0; i<Epot_pure.size; i++) Epot_pure.x[i] = Epot_pure.step * (DOUBLE) i;
  }

  if(measure_RadDistr) {
    RD.N = (int*) Calloc("RD", gridRD, sizeof(int));
    RD.size = gridRD;
    if(RD.max == 0) RD.max = L;
    RD.step = RD.max/gridRD;
  }

  if(MC == DIFFUSION && measure_RadDistr) {
    RD_pure.N = (int*) Calloc("RD_pure", gridRD, sizeof(int));
    RD_pure.size = gridRD;
    RD_pure.max = RD.max;
    RD_pure.step = RD.step;
  }

  RDz.N = (int*) Calloc("RDz.N", gridRD, sizeof(int));
  RDz.f = (DOUBLE*) Calloc("RDz.f", gridRD, sizeof(DOUBLE));

  RDz.size = gridRD;
  RDz.min = 0;
  if(RDz.max == 0) RDz.max = Lz;
#ifdef TRIAL_2D // save y component in 2D hom. system
#ifdef BC_2DPBC
  if(RDz.max == 0) RDz.max = Ly;
#endif
#ifdef BC_1DPBC_X
  if(RDz.max == 0) RDz.max = Ly;
#endif
#endif
  RDz.step = RDz.max/gridRD;

  if(MC == DIFFUSION) {
    RDz_pure.N = (int*) Calloc("RDz_pure", gridRD, sizeof(int));
    RDz_pure.size = gridRD;
    RDz_pure.step = RDz.step;
    RDz_pure.max = RDz.max;
  }

  if(measure_g3 && MC == DIFFUSION) {
    HR_pure.N = (int*) Calloc("HR_pure", gridg3, sizeof(int));
    HR_pure.size = gridRD;
    HR_pure.step = HR.max/gridg3;
  }

  if(measure_effective_potential) {
    Veff.size = gridRD;
    Veff.step = RD.max/gridRD;
    Veff.N = (int*) Calloc("Veff", Veff.size, sizeof(int));
    Veff.f = (DOUBLE*) Calloc("Veff", Veff.size, sizeof(DOUBLE));
  }

  if(measure_Sk) {
    Sk.f = (DOUBLE*) Calloc("Sk->f", gridSk, sizeof(DOUBLE));
    Sk.cos = (DOUBLE*) Calloc("Sk->cos", gridSk, sizeof(DOUBLE));
    Sk.sin = (DOUBLE*) Calloc("Sk->sin", gridSk, sizeof(DOUBLE));
    Sk.N = (unsigned long int*) Calloc("Sk->N", gridSk, sizeof(unsigned long int));
    Sk.size = gridSk;
    Sk.k =  (DOUBLE*) Calloc("SK->k", gridSk, sizeof(DOUBLE));
    Sk.kx = (DOUBLE*) Calloc("SK->kx", gridSk, sizeof(DOUBLE));
    Sk.ky = (DOUBLE*) Calloc("SK->ky", gridSk, sizeof(DOUBLE));
    Sk.kz = (DOUBLE*) Calloc("SK->kz", gridSk, sizeof(DOUBLE));
    Sk.index = (int*) Calloc("SK->index", gridSk, sizeof(int));
    Sk.degeneracy = (int*) Calloc("SK->degeneracy", gridSk, sizeof(int));

    if(Sk.L == 0) {
      Sk.L = L;
      Message("  Sk max = %" LE " \n", Sk.L);
    }
    else {
      Sk.L = 2.*PI*Sk.size/Sk.L;
    }

    for(i=0; i<Sk.size; i++) { // initialize the values of the momenta (1D)
      Sk.k[i] = 2.*PI/Sk.L * (i+1.);
    }
  }

  if(measure_Sk && MC == DIFFUSION) {
    Sk_pure.f = (DOUBLE*) Calloc("Sk_pure.f", gridSk, sizeof(DOUBLE));
    Sk_pure.N = (unsigned long int*) Calloc("Sk_pure.N", gridSk, sizeof(unsigned long int));
    Sk_pure.times_measured = 0;
    Sk_pure.size = gridSk;
    Sk_pure.L = Sk.L;
  }

  if(measure_Lind && MC == DIFFUSION) {
     LindemannRatio.Fpure = (DOUBLE*) Calloc("LindemannRatio.Fpure", grid_pure_block, sizeof(DOUBLE));
     LindemannRatio.Npure = (int*) Calloc("LindemannRatio.Npure", grid_pure_block, sizeof(int));
  }

  if(measure_SD) {
    SD.CM2 = (DOUBLE*) Calloc("SD->CM2", SD.size, sizeof(DOUBLE));
    SD.N   = (int*) Calloc("SD->N", SD.size, sizeof(int));
    SD.times_measured = 0;
    //SD.size = 100;
    //SD.spacing = 1;
  }

  u_mi = (DOUBLE*) Calloc("u_mi", N, sizeof(DOUBLE));
  u_ij = (DOUBLE*) Calloc("u_ij", N, sizeof(DOUBLE));
#ifdef OBDM_FERMIONS
  sign_u_mi = (int*) Calloc("sign_u_mi", N, sizeof(int));
  sign_u_ij = (int*) Calloc("sign_u_ij", N, sizeof(int));
  order = (DOUBLE*) Calloc("order", N, sizeof(DOUBLE));
  OBDMfermi.f = (DOUBLE*) Calloc("OBDMfermi->f", gridOBDM, sizeof(DOUBLE));
  OBDMfermi.Nk = (DOUBLE*) Calloc("OBDMfermi->f", OBDM.Nksize, sizeof(DOUBLE));
#endif
  mu_k = (DOUBLE*) Calloc("mu_k", N, sizeof(DOUBLE));
  mu_p_k =  (DOUBLE*) Calloc("mu_p_k", N, sizeof(DOUBLE));
  mu_pp_k =  (DOUBLE*) Calloc("mu_pp_k", N, sizeof(DOUBLE));

  if(measure_TBDM) {
    u_tbdm_i = (DOUBLE*) Calloc("u_tbdm_i", N, sizeof(DOUBLE));
    u_tbdm_mcmillan_i = (DOUBLE*) Calloc("u_tbdm_mcmillan_i", N, sizeof(DOUBLE));
    ArrayCalloc2D(u_tbdm_ij, "u_tbdm_ij ", i, N, N, DOUBLE, "DOUBLE");
    ArrayCalloc2D(u_tbdm_mcmillan_ij, "u_tbdm_mcmillan_ij ", i, N, N, DOUBLE, "DOUBLE");
  }
}
/************************** Calloc ******************************************/
void* Calloc(const char* name, unsigned length, size_t size) {
   void *p = calloc(length, size);
   if(p == NULL) Error("Not enough memory for %s", name);

   return p;
}

/************************* CallocContiguous2D ********************************/
void* CallocContiguous2D(const char* name, unsigned dim1, unsigned dim2, char *type) {
// Allocate 2D array of size [dim1][dim2]
// type is a string: "DOUBLE", "double", "float", "int"
  DOUBLE **pDOUBLE;
  double **pdouble;
  float **pfloat;
  int **pint;
  unsigned i;

  if(strcmp(type, "DOUBLE") == 0) {
    pDOUBLE = (DOUBLE**) Calloc(name, dim1, sizeof(DOUBLE*));
    pDOUBLE[0] = (DOUBLE*) Calloc(name, dim1*dim2, sizeof(DOUBLE)); // allocate contiguous memory
    for(i=1; i<dim1; i++) pDOUBLE[i] = pDOUBLE[0] + i*dim2; // set pointers
    return pDOUBLE;
  }
  else if(strcmp(type, "double") == 0) {
    pdouble = (double**) Calloc(name, dim1, sizeof(double*));
    pdouble[0] = (double*) Calloc(name, dim1*dim2, sizeof(double)); // allocate contiguous memory
    for(i=1; i<dim1; i++) pdouble[i] = pdouble[0] + i*dim2; // set pointers
    return pdouble;
  }
  else if(strcmp(type, "float") == 0) {
    pfloat = (float**) Calloc(name, dim1, sizeof(float*));
    pfloat[0] = (float*) Calloc(name, dim1*dim2, sizeof(float)); // allocate contiguous memory
    for(i=1; i<dim1; i++) pfloat[i] = pfloat[0] + i*dim2; // set pointers
    return pfloat;
  }
  else { //if(strcmp(type, "int") == 0) {
    pint = (int**) Calloc(name, dim1, sizeof(int*));
    pint[0] = (int*) Calloc(name, dim1*dim2, sizeof(int)); // allocate contiguous memory
    for(i=1; i<dim1; i++) pint[i] = pint[0] + i*dim2; // set pointers
    return pint;
  }
}

/************************* CallocContiguous3D ********************************/
void* CallocContiguous3D(const char* name, unsigned dim1, unsigned dim2, unsigned dim3, char *type) {
// Allocate 3D array of size [dim1][dim2][dim3]
// type is a string: "DOUBLE", "double", "float", "int"
  DOUBLE ***pDOUBLE;
  double ***pdouble;
  float ***pfloat;
  int ***pint;
  unsigned i,j;

  if(strcmp(type, "DOUBLE") == 0) {
    pDOUBLE = (DOUBLE***) Calloc(name, dim1, sizeof(DOUBLE**));
    pDOUBLE[0] = (DOUBLE**) Calloc(name, dim1*dim2, sizeof(DOUBLE*));
    pDOUBLE[0][0] = (DOUBLE*) Calloc(name, dim1*dim2*dim3, sizeof(DOUBLE));
    for(i=0; i<dim1; i++) {
      pDOUBLE[i] = pDOUBLE[0] + i*dim2;
      pDOUBLE[i][0] = pDOUBLE[0][0] + i*(dim2*dim3);
      for(j=1; j<dim2; j++) {
        pDOUBLE[i][j] = pDOUBLE[i][0] + j*dim3;
      }
    }
    return pDOUBLE;
  }
  else if(strcmp(type, "double") == 0) {
    pdouble = (double***) Calloc(name, dim1, sizeof(double**));
    pdouble[0] = (double**) Calloc(name, dim1*dim2, sizeof(double*));
    pdouble[0][0] = (double*) Calloc(name, dim1*dim2*dim3, sizeof(double));
    for(i=0; i<dim1; i++) {
      pdouble[i] = pdouble[0] + i*dim2;
      pdouble[i][0] = pdouble[0][0] + i*(dim2*dim3);
      for(j=1; j<dim2; j++) {
        pdouble[i][j] = pdouble[i][0] + j*dim3;
      }
    }
    return pdouble;
  }
  else if(strcmp(type, "float") == 0) {
    pfloat = (float***) Calloc(name, dim1, sizeof(float**));
    pfloat[0] = (float**) Calloc(name, dim1*dim2, sizeof(float*));
    pfloat[0][0] = (float*) Calloc(name, dim1*dim2*dim3, sizeof(float));
    for(i=0; i<dim1; i++) {
      pfloat[i] = pfloat[0] + i*dim2;
      pfloat[i][0] = pfloat[0][0] + i*(dim2*dim3);
      for(j=1; j<dim2; j++) {
        pfloat[i][j] = pfloat[i][0] + j*dim3;
      }
    }
    return pfloat;
  }
  else { //if(strcmp(type, "int") == 0) {
    pint = (int***) Calloc(name, dim1, sizeof(int**));
    pint[0] = (int**) Calloc(name, dim1*dim2, sizeof(int*));
    pint[0][0] = (int*) Calloc(name, dim1*dim2*dim3, sizeof(int));
    for(i=0; i<dim1; i++) {
      pint[i] = pint[0] + i*dim2;
      pint[i][0] = pint[0][0] + i*(dim2*dim3);
      for(j=1; j<dim2; j++) {
        pint[i][j] = pint[i][0] + j*dim3;
      }
    }
    return pint;
  }
}


/********************************* Copy Walker *******************************/
void CopyWalkerCoord(struct Walker *out, const struct Walker *in) {

#ifndef MEMORY_CONTIGUOUS
  int i;
#endif
  CaseX(ArrayCopy1D(out->x, in->x, i, N));
  CaseY(ArrayCopy1D(out->y, in->y, i, N));
  CaseZ(ArrayCopy1D(out->z, in->z, i, N));
#ifdef SPINFULL
  ArrayCopy1D(out->spin, in->spin, i, N);
#endif
#ifdef SPINFULL_TUNNELING
  ArrayCopy1D(out->psiT_sigma_inv, in->psiT_sigma_inv, i, N);
#endif
}

/************************** Copy Raw Vector to Walker ************************/
void CopyRawVectorToWalkerDisplace_i_j(DOUBLE *r, struct Walker *W, int index1, DOUBLE dx1, int index2, DOUBLE dx2) {
// copies vector r[] of size (3N) to the Walker
// shifts poistions i and j by dx1 and dx2
  int i, index;

  index = 0;
  //CaseX(for(i=0;i<N;i++) W.x[i] = r[index++]);
  CaseX(for(i=0;i<N;i++) {W->x[i] = r[index]; if(index == index1) W->x[i]+=dx1; if(index == index2) W->x[i]+=dx2; index++;});
  CaseY(for(i=0;i<N;i++) {W->y[i] = r[index]; if(index == index1) W->y[i]+=dx1; if(index == index2) W->y[i]+=dx2; index++;});
  CaseZ(for(i=0;i<N;i++) {W->z[i] = r[index]; if(index == index1) W->z[i]+=dx1; if(index == index2) W->z[i]+=dx2; index++;});
}

/********************************* Walker Compare ****************************/
int WalkerCompare(const struct Walker x, const struct Walker y) {
  int i;
  int ret = 0;

  for(i=0; i<N; i++) {
    CaseX(if(x.x[i] != y.x[i]) ret = 1);
    CaseY(if(x.y[i] != y.y[i]) ret = 1);
    CaseZ(if(x.z[i] != y.z[i]) ret = 1);
  }
  if((x.E != y.E) && !isnan(x.E) && !isnan(y.E)) ret = 1; // comparison is not succesfull when E is infinite
  if((x.U != y.U) && !isnan(x.U) && !isnan(y.U)) ret = 1;

  if(ret) {
    for(i=0; i<N; i++) {
      Warning("%i   %" LF " %" LF " %" LF "  %" LF " %" LF " %" LF "\n", i+1, x.x[i], x.y[i], x.z[i], y.x[i], y.y[i], y.z[i]);
    }
    Warning("%" LF " %" LF "\n", x.E, y.E);
    Warning("%" LF " %" LF "\n", x.U, y.U);
  }

  return ret;
}

/************************ Vector to Walker Compare ***************************/
int VectorToWalkerCompare(const DOUBLE **x, const struct Walker y) {
  int i;
  int ret = 0;

  for(i=0; i<N; i++) {
    CaseX(if(x[i][0] != y.x[i]) ret = 1);
    CaseY(if(x[i][1] != y.y[i]) ret = 1);
    CaseZ(if(x[i][2] != y.z[i]) ret = 1);
  }

  if(ret) {
    for(i=0; i<N; i++) {
      Warning("%i   %" LF " %" LF " %" LF "  %" LF " %" LF " %" LF "\n", i+1, x[i][0], x[i][1], x[i][2], y.x[i], y.y[i], y.z[i]);
    }
  }

  return ret;
}