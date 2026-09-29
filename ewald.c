#include "main.h"
#include "vmc.h"
#include "utils.h"
#include "compatab.h"
#include MATHINCLUDE
#include "memory.h"
#include "ewald.h"
#include <stdlib.h>
#include "mymath.h"

DOUBLE Ewald_alpha = 0.5;
DOUBLE Ewald_alpha2,Ewald_alpha3,Ewald_alpha4,Ewald_alpha5,Ewald_alpha6;
int Ewald_Nx = 3; // number of terms in real space
int Ewald_Nk = 6; // number of terms in momentum space
int Ewald_Nk_max = 9; // used to calclate "exact" value and allocate memory
int Ewald_Nx_check_convergence = 9; // 1< Nx <Ewald_Nk_max for doing test
int Ewald_Nk_check_convergence = 9; // 1< Nx <Ewald_Nk_max

#ifdef TRIAL_2D
#ifdef INTERACTION_COULOMB
#  define Ewald_unit_Lk(L) (L)
#  define Ewald_rho(n,n2) (erfc(Ewald_alpha*n)/n)
#  define Ewald_kappa(n,n2) (erfc(PI*n/Ewald_alpha)/n)
#  define Ewald_C1 (-2.*sqrt(PI)/Ewald_alpha)
#  define Ewald_C2 (-2.*Ewald_alpha/sqrt(PI))
#endif

#ifdef INTERACTION_DIPOLE
#  define Ewald_unit_Lk(L) (L*L*L)
#  define Ewald_rho(n,n2) (erfc(Ewald_alpha*n)/(n*n2) + 2.*Ewald_alpha/(sqrt(PI)*n2)*exp(-Ewald_alpha2*n2))
#  define Ewald_kappa(n,n2) (4.*(sqrt(PI)*Ewald_alpha*exp(-PI2*n*n/Ewald_alpha2)-PI2*n*erfc(PI*n/Ewald_alpha)))
#  define Ewald_C1 (4.*sqrt(PI)*Ewald_alpha)
#  define Ewald_C2 (- 4.*Ewald_alpha3/(3.*sqrt(PI)))
#endif

#ifdef INTERACTION_R4
#  define Ewald_unit_Lk(L) (L*L*L*L)
#  define Ewald_rho(n,n2) (exp(-Ewald_alpha2*n2)*(Ewald_alpha2*n2+1)/(n2*n2))
#  define Ewald_kappa(n,n2) (PI*Ewald_alpha4*(exp(-PI2*n2/Ewald_alpha2)-PI2*n2/Ewald_alpha2*E1(PI2*n2/Ewald_alpha2)))
#  define Ewald_C1 (PI*Ewald_alpha2)
#  define Ewald_C2 (-Ewald_alpha4/2.)
#endif

#ifdef INTERACTION_R5
#  define Ewald_unit_Lk(L) (L*L*L*L*L)
#  define Ewald_rho(n,n2) ((erfc(Ewald_alpha*n)+4.*exp(-Ewald_alpha2*n2)/(3.*sqrt(PI))*(3.*Ewald_alpha*n/2.+Ewald_alpha3*n*n2))/(n*n2*n2))
#  define Ewald_kappa(n,n2) (8.*sqrt(PI)*Ewald_alpha3/9.*(exp(-PI2*n2/Ewald_alpha2)*(1.-2.*PI2*n2/Ewald_alpha2)+2.*PI3*sqrt(PI)*n2*n/Ewald_alpha3*erfc(PI*n/Ewald_alpha)))
#  define Ewald_C1 (8.*sqrt(PI)*Ewald_alpha3/9.)
#  define Ewald_C2 (-8.*Ewald_alpha5/(15.*sqrt(PI)))
#endif

#ifdef INTERACTION_RYDBERG
#  define Ewald_unit_Lk(L) (L*L*L*L*L*L)
#  define Ewald_rho(n,n2) ((Ewald_alpha4/(2.*n2)+Ewald_alpha2/(n2*n2)+1./(n2*n2*n2))*exp(-Ewald_alpha2*n2))
#  define Ewald_kappa(n,n2) (PI*Ewald_alpha4/4.*(exp(-PI2*n2/Ewald_alpha2)*(1.-PI2*n2/Ewald_alpha2)+PI4*n2*n2/Ewald_alpha4*E1(PI2*n2/Ewald_alpha2)))
#  define Ewald_C1 (0.25*PI*Ewald_alpha4)
#  define Ewald_C2 (- Ewald_alpha6/6.)
#endif

#ifdef INTERACTION_R12
#  define Ewald_unit_Lk(L) (L*L*L*L*L*L*L*L*L*L*L*L)
#  define Ewald_rho(n,n2) (exp(-Ewald_alpha2*n2)*(1.+Ewald_alpha2*n2+Ewald_alpha4*n2*n2/2.+Ewald_alpha2*n2*Ewald_alpha2*n2*Ewald_alpha2*n2/6.+Ewald_alpha2*n2*Ewald_alpha2*n2*Ewald_alpha2*n2*Ewald_alpha2*n2/24.+Ewald_alpha2*n2*Ewald_alpha2*n2*Ewald_alpha2*n2*Ewald_alpha2*n2*Ewald_alpha2*n2/120.)/(n2*n2*n2*n2*n2*n2))
#  define Ewald_C1 (PI*Ewald_alpha5*Ewald_alpha5/600.)
#  define Ewald_C2 (-Ewald_alpha5*Ewald_alpha5*Ewald_alpha2/720.)
#  define Ewald_kappa(n, n2) exp(-(PI2*n2/Ewald_alpha2))*(1./5.-(PI2*n2/Ewald_alpha2)/20.+(PI2*n2/Ewald_alpha2)*(PI2*n2/Ewald_alpha2)/60.-(PI2*n2/Ewald_alpha2)*(PI2*n2/Ewald_alpha2)*(PI2*n2/Ewald_alpha2)/120.+(PI2*n2/Ewald_alpha2)*(PI2*n2/Ewald_alpha2)*(PI2*n2/Ewald_alpha2)*(PI2*n2/Ewald_alpha2)/120.-(PI2*n2/Ewald_alpha2)*(PI2*n2/Ewald_alpha2)*(PI2*n2/Ewald_alpha2)*(PI2*n2/Ewald_alpha2)*(PI2*n2/Ewald_alpha2)/120.*E1((PI2*n2/Ewald_alpha2)))
#endif
#endif // end 2D

#ifdef TRIAL_3D
#ifdef INTERACTION_COULOMB
#  define Ewald_unit_Lk(L) (L)
#  define Ewald_rho(n,n2) (erfc(Ewald_alpha*n)/n)
#  define Ewald_kappa(n,n2) (exp(-PI2*n2/Ewald_alpha2)/(PI*n2))
#  define Ewald_C1 (-PI/Ewald_alpha2)
#  define Ewald_C2 (-2.*Ewald_alpha/sqrt(PI))
#endif

#ifdef INTERACTION_R2
#  define Ewald_unit_Lk(L) (L*L)
#  define Ewald_rho(n,n2) (exp(-Ewald_alpha2*n2)/n2)
#  define Ewald_kappa(n,n2) (PI*erfc(PI*n/Ewald_alpha)/n)
#  define Ewald_C1 (-2.*PI*sqrt(PI)/Ewald_alpha)
#  define Ewald_C2 (-Ewald_alpha2)
#endif

#ifdef INTERACTION_CALOGERO
#  define Ewald_unit_Lk(L) (L*L)
#  define Ewald_rho(n,n2) (D*exp(-Ewald_alpha2*n2)/n2)
#  define Ewald_kappa(n,n2) (D*PI*erfc(PI*n/Ewald_alpha)/n)
#  define Ewald_C1 (-D*2.*PI*sqrt(PI)/Ewald_alpha)
#  define Ewald_C2 (-D*Ewald_alpha2)
#endif

#ifdef INTERACTION_R4
#  define Ewald_unit_Lk(L) (L*L*L*L)
#  define Ewald_rho(n,n2) (exp(-Ewald_alpha2*n2)*(Ewald_alpha2*n2+1)/(n2*n2))
#  define Ewald_kappa(n,n2) (2.*PI*(sqrt(PI)*Ewald_alpha*exp(-PI2*n2/Ewald_alpha2)-PI2*n*erfc(PI*n/Ewald_alpha)))
#  define Ewald_C1 (2.*PI*(sqrt(PI)*Ewald_alpha))
#  define Ewald_C2 (-Ewald_alpha4/2.)
#endif

#ifdef INTERACTION_R5
#  define Ewald_unit_Lk(L) (L*L*L*L*L)
#  define Ewald_rho(n,n2) ((erfc(Ewald_alpha*n)+4.*exp(-Ewald_alpha2*n2)/(3.*sqrt(PI))*(3.*Ewald_alpha*n/2.+Ewald_alpha3*n*n2))/(n*n2*n2))
#  define Ewald_kappa(n,n2) (4.*PI*Ewald_alpha2/3.*(exp(-PI2*n2/Ewald_alpha2) - PI2*n2/Ewald_alpha2*E1(PI2*n2/Ewald_alpha2)))
#  define Ewald_C1 (4.*PI*Ewald_alpha2/3.)
#  define Ewald_C2 (-8.*Ewald_alpha5/(15.*sqrt(PI)))
#endif

#ifdef INTERACTION_RYDBERG
#  define Ewald_unit_Lk(L) (L*L*L*L*L*L)
#  define Ewald_rho(n,n2) ((Ewald_alpha4/(2.*n2)+Ewald_alpha2/(n2*n2)+1./(n2*n2*n2))*exp(-Ewald_alpha2*n2))
#  define Ewald_kappa(n,n2) (PI*sqrt(PI)*Ewald_alpha3/3.*(exp(-PI2*n2/Ewald_alpha2)*(1.-2.*PI2*n2/Ewald_alpha2)+2.*PI3*sqrt(PI)*n*n2/Ewald_alpha3*erfc(PI*n/Ewald_alpha)))
#  define Ewald_C1 (PI*sqrt(PI)*Ewald_alpha3/3.)
#  define Ewald_C2 (- Ewald_alpha6/6.)
#endif
#endif // end 3D

#ifndef Ewald_C1 // dummy
#  define Ewald_unit_Lk(L) 1.
#  define Ewald_rho(n,n2) 1.
#  define Ewald_C1 1.
#  define Ewald_C2 1.
#  define Ewald_kappa(n, n2) 1.
#endif

DOUBLE Ewald_Xi;
DOUBLE Ewald_Lk, Ewald_L;

#ifdef LATTICE_TRIANGULAR
DOUBLE **Ewald_cos_coef; // [Nk][Nk] coefficients of the sum in momentum space
DOUBLE Ewald_Lx_Ly, Ewald_Ly_Lx;
DOUBLE Ewald_Lx_Ly_sq, Ewald_Ly_Lx_sq;
#else
DOUBLE *Ewald_cos_coef; // [Nk*Nk] coefficients of the sum in momentum space
#endif

DOUBLE *Ewald_Exact; // [Nwalkers] 'exact' energy
DOUBLE Ewald_Exact_mean; // averaged over walkers

/**************************** Ewald Check ************************************/
int EwaldCheck(void) {
  DOUBLE dr[3], r2, r;
  int i, j, k, w, Nx, Nx_store, Nk, Nk_store;
  int max_direct = 5; // max number of cells in direct summation
  DOUBLE Echeck; // store energy
  clock_t time_check; // and time
  DOUBLE Nx_num_store, alpha_store;
  int nx = 0;
  int ny = 0;
  int nz = 0;
  int kmax;
  DOUBLE Eself;
  DOUBLE EcheckEwaldW0; // Ewald energy for the first walker

  if(fabs(erfc(1.)-0.15729920705028513)>1e-6) Error("  erfc() is not implemented");
  if(fabs(E1(1.)-0.21938393439552026)>1e-6) Error("  E1() is not implemented");

  // test consistency of EwaldSumWalker() and EwaldSumPair()
  Message("\n  Checking Ewald (energy per particle): compare \n");

  // initialize Ewald summation
  //EwaldInit();

  // store parameters
  Nx_store = Ewald_Nx;
  Nk_store = Ewald_Nk;
  alpha_store = Ewald_alpha;

  // calculate "Exact energy"
  Ewald_Nx = Ewald_Nk_max;
  Ewald_Nk = Ewald_Nk_max;
  Ewald_Exact_mean = EwaldSumAllWalkers();
  Message("\n  Exact Ewald energy is %.15lf\n\n", Ewald_Exact_mean/(double)N);

  // self energy convergence
  Message("  Convergence of the Self energy per particle\n");
  Message("  Xi %.15lf\n", 0.5*Ewald_Xi*Ewald_Lk);
  InteractionSelfEnergy();
  Message("  direct sum:\n");
  Message("  cell E\n");

  if(DIMENSION == 2)
    kmax = 10000;
  else
    kmax = 100;

  for(k=10; k<kmax; k*=10) {
    Eself = 0.;
#ifdef MOVE_IN_X
    for(nx=-k; nx<=k; nx++)
#endif
    {
#ifdef MOVE_IN_Y
      for(ny=-k; ny<=k; ny++)
#endif
      {
#ifdef MOVE_IN_Z
        for(nz=-k; nz<=k; nz++)
#endif
        {
          if(nx!=0 || ny!=0 || nz!=0) {
            Eself += InteractionEnergy(Sqrt((nx*Lx)*(nx*Lx)+(ny*Ly)*(ny*Ly)+(nz*Lz)*(nz*Lz)));
          }
        }
      }
    }
   Eself /= (double) (N-1); // self energy 
   Message("  %i %lf\n", k, 0.5*Eself*(double)(N-1));
  }

  time_check = -clock();
  Echeck = EwaldSumWalker(&W[0]);
  time_check += clock();
  EcheckEwaldW0 = Echeck/(double)N;
  Message("  E = %lf, time %lf (standard sum)\n", EcheckEwaldW0,(DOUBLE)time_check/(DOUBLE)CLOCKS_PER_SEC);

  Echeck = 0;
  time_check = -clock();
  for(i=0; i<N; i++) { // two-body terms: particle-particle force
    for(j=i+1; j<N; j++) {
      dr[0] = W[0].x[i] - W[0].x[j];
      dr[1] = W[0].y[i] - W[0].y[j];
      dr[2] = W[0].z[i] - W[0].z[j];
      Echeck += EwaldSumPair(dr[0],dr[1],dr[2]);
    }
  }
  time_check += clock();
  Message("  E = %lf, time %lf (summation by pairs)\n", Echeck/(double)N,(DOUBLE)time_check/(DOUBLE)CLOCKS_PER_SEC);

  // test convergence with Nx
  Message("  Checking convergence with number of cells\n");
  Message("    Nx Nk Nx_cells Ny_cells E dE/E time\n");

  for(Nx=1; Nx<=Ewald_Nx_check_convergence; Nx++) {
    Ewald_Nx = Nx;
    for(Nk=1; Nk<=1; Nk++) {
      Ewald_Nk = Nk;
      time_check = -clock();
      Echeck = EwaldSumAllWalkers();
      time_check += clock();
      Message("    %3i %3i %3i %3i %.15e %.15e %e\n", Nx, Nk, (2*Nx-1)*(2*Nx-1), (2*Nk-1)*(2*Nk-1), Echeck/(double)N, Echeck/Ewald_Exact_mean-1., (DOUBLE)time_check/(DOUBLE)CLOCKS_PER_SEC);
    }
  }

  for(Nx=1; Nx<=1; Nx++) {
    Ewald_Nx = Nx;
    for(Nk=1; Nk<=Ewald_Nk_check_convergence; Nk++) {
      Ewald_Nk = Nk;
      time_check = -clock();
      Echeck = EwaldSumAllWalkers();
      time_check += clock();
      Message("    %3i %3i %3i %3i %.15e %.15e %e\n", Nx, Nk, (2*Nx-1)*(2*Nx-1), (2*Nk-1)*(2*Nk-1), Echeck/(double)N, Echeck/Ewald_Exact_mean-1., (DOUBLE)time_check/(DOUBLE)CLOCKS_PER_SEC);
    }
  }

  Message("  Checking Ewald summation against summation by pairs\n");
  Echeck = 0;
  time_check = -clock();
  for(i=0; i<N; i++) { // two-body terms: particle-particle force
    for(j=i+1; j<N; j++) {
      dr[0] = W[0].x[i] - W[0].x[j];
      dr[1] = W[0].y[i] - W[0].y[j];
      dr[2] = W[0].z[i] - W[0].z[j];
      r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
      r = Sqrt(r2);
      Echeck += EwaldSumPair(dr[0],dr[1],dr[2]);
    }
  }
  Echeck/= (double)N;
  time_check += clock();

  Message("     E = %lf, dE/E = %e, time %lf (Ewald summation)\n", Echeck, Echeck/EcheckEwaldW0-1., (DOUBLE)time_check/(DOUBLE)CLOCKS_PER_SEC);

  // first time self energy is calculated
  InteractionEnergySumOverImages(1,1,1);

  Message("\n   Direct summation (no Ewald):\n");
  Message("     cutoff  E  dE/E  time\n");
  Nx_num_store = sum_over_images_cut_off_num;
  for(k=1; k<max_direct; k++) {
    sum_over_images_cut_off_num = (double) k; // g2 =1 summation
    Ewald_Nx = k; // direct summation
    Echeck = 0.;
    time_check = -clock();
    for(w=0; w<1; w++) {
      for(i=0; i<N; i++) { // two-body terms: particle-particle force
        for(j=i+1; j<N; j++) {
          dr[0] = W[0].x[i] - W[0].x[j];
          dr[1] = W[0].y[i] - W[0].y[j];
          dr[2] = W[0].z[i] - W[0].z[j];
          r2 = FindNearestImage(&dr[0], &dr[1], &dr[2]);
          r = Sqrt(r2);
          Echeck += InteractionEnergySumOverImages(dr[0],dr[1],dr[2]);
        }
      }
    }
    time_check += clock();
    //Echeck /= (DOUBLE) Nwalkers;
    Echeck /= (DOUBLE)N;
    //Message("     %i %lf %e %e\n", k, Echeck/(double)N, Echeck/Ewald_Exact_mean-1., (DOUBLE)time_check/(DOUBLE)CLOCKS_PER_SEC);
    Message("     %i %lf %e %e\n", k, Echeck, Echeck/EcheckEwaldW0-1., (DOUBLE)time_check/(DOUBLE)CLOCKS_PER_SEC);
  }

  //EwaldOptimizeInitExact();
  //Ewald_alpha = EwaldOptimize(1, 1);
  //Ewald_alpha = EwaldOptimize(1, 2);

  Message("  Ewald parameters: alpha=%lf Nx=%i Nk=%i\n", Ewald_alpha, Ewald_Nx, Ewald_Nk);

  if(fabs(Echeck/EcheckEwaldW0-1.) > 1e-3)
    Error("  Ewald summation test failed (direct: %e, Ewald: %e). Exiting.\n", Echeck, EcheckEwaldW0);
  else
    Message("  Ewald summation test passed.\n");

  Ewald_Nx = Nx_store;
  Ewald_Nk = Nk_store;
  Ewald_alpha = alpha_store;
  sum_over_images_cut_off_num = Nx_num_store;
  EwaldInitCoef();

  return 0;
}

/**************************** Ewald Init *************************************/
void EwaldInit(void) {
// initialize arrays
#ifdef LATTICE_TRIANGULAR
  int i;

  Ewald_Exact = (DOUBLE*) Calloc("Ewald_Exact", Nwalkers, sizeof(DOUBLE));
  Ewald_Ly_Lx = Ly/Lx;
  Ewald_Lx_Ly = Lx/Ly;
  Ewald_Ly_Lx_sq = sqrt(Ly/Lx);
  Ewald_Lx_Ly_sq = sqrt(Lx/Ly);
  Ewald_L = 1./sqrt(Lx*Ly);
  Ewald_Lk = Ewald_unit_Lk(Ewald_L);
  ArrayCalloc2D(Ewald_cos_coef, "Ewald cos coeff ", i, Ewald_Nk_max+1, Ewald_Nk_max+1, DOUBLE, "DOUBLE");
#else
  Ewald_Exact = (DOUBLE*) Calloc("Ewald_Exact", Nwalkers, sizeof(DOUBLE));
  Ewald_L = 1./L;
  Ewald_Lk = Ewald_unit_Lk(Ewald_L);
  Ewald_cos_coef = (DOUBLE*) Calloc("Ewald cos coeff", (Ewald_Nk_max+1)*(Ewald_Nk_max+1), sizeof(DOUBLE));
#endif

  Ewald_alpha2 = Ewald_alpha*Ewald_alpha;
  Ewald_alpha3 = Ewald_alpha2*Ewald_alpha;
  Ewald_alpha4 = Ewald_alpha2*Ewald_alpha2;
  Ewald_alpha5 = Ewald_alpha4*Ewald_alpha;
  Ewald_alpha6 = Ewald_alpha3*Ewald_alpha3;

  EwaldInitCoef();
  Ewald_Xi = EwaldXi();

  EwaldPsi(Lx/2.*Ewald_L, Ly/2.*Ewald_L, 0);
}

/**************************** Ewald Init Coef ********************************/
void EwaldInitCoef(void) {
  int nx = 0;
  int ny = 0;
  int nz = 0;
  DOUBLE n;
#ifdef LATTICE_TRIANGULAR
  DOUBLE n2;
#else
  int n2;
#endif

  // sum over momentum space
#ifdef MOVE_IN_X
  for(nx=-Ewald_Nk_max; nx<=Ewald_Nk_max; nx++)
#endif
  {
#ifdef MOVE_IN_Y
    for(ny=-Ewald_Nk_max; ny<=Ewald_Nk_max; ny++)
#endif
    {
#ifdef MOVE_IN_Z
      for(nz=-Ewald_Nk_max; nz<=Ewald_Nk_max; nz++)
#endif
      {
#ifdef LATTICE_TRIANGULAR
        n2 = nx*nx*Ewald_Ly_Lx + ny*ny*Ewald_Lx_Ly;
#else
        n2 = nx*nx + ny*ny + nz*nz;
#endif

        n = sqrt((DOUBLE)n2);
#ifdef LATTICE_TRIANGULAR
        if(n2 == 0)
          Ewald_cos_coef[0][0] = 0;
        else
          Ewald_cos_coef[abs(nx)][abs(ny)] = Ewald_kappa(n,(DOUBLE)n2);
#else
        if(n2 == 0)
          Ewald_cos_coef[0] = 0;
        else
          Ewald_cos_coef[n2] = Ewald_kappa(n,(DOUBLE)n2);
#endif
      }
    }
  }
}

/**************************** Ewald Xi ***************************************/
DOUBLE EwaldXi(void) {
  DOUBLE xi = 0.;
  DOUBLE n;
  int nx = 0;
  int ny = 0;
  int nz = 0;
  int nmax = 100;
#ifdef LATTICE_TRIANGULAR
  DOUBLE n2;

  // sum over real space
  for(nx=-nmax; nx<=nmax; nx++) {
    for(ny=-nmax; ny<=nmax; ny++) {
      n2 = nx*nx*Ewald_Lx_Ly + ny*ny*Ewald_Ly_Lx;
      n = sqrt((DOUBLE)n2);
      if(n2>0) xi += Ewald_rho(n,(DOUBLE)n2);
    }
  }

  // sum over real space
  for(nx=-nmax; nx<=nmax; nx++) {
    for(ny=-nmax; ny<=nmax; ny++) {
      n2 = nx*nx*Ewald_Ly_Lx + ny*ny*Ewald_Lx_Ly;
      n = sqrt((DOUBLE)n2);
      if(n2>0) xi += Ewald_kappa(n,(DOUBLE)n2);
    }
  }
#else
  int n2;

  // sum over real and momentum space
#ifdef MOVE_IN_X
  for(nx=-nmax; nx<=nmax; nx++)
#endif
  {
#ifdef MOVE_IN_Y
    for(ny=-nmax; ny<=nmax; ny++)
#endif
    {
#ifdef MOVE_IN_Z
      for(nz=-nmax; nz<=nmax; nz++)
#endif
      {
        n2 = nx*nx + ny*ny + nz*nz;
        n = sqrt((DOUBLE)n2);
        if(n2>0) {
          xi += Ewald_kappa(n,(DOUBLE)n2);
          xi += Ewald_rho(n,(DOUBLE)n2);
        }
      }
    }
  }
#endif
  xi += Ewald_C1 + Ewald_C2;

  return xi;
}

/**************************** Ewald Psi **************************************/
DOUBLE EwaldPsi(DOUBLE x, DOUBLE y, DOUBLE z) {
// x, y, z are reduced units, i.e x/L, y/L
  int nx = 0;
  int ny = 0;
  int nz = 0;
#ifndef LATTICE_TRIANGULAR
  int n2;
#endif
  DOUBLE psi = 0.;
  DOUBLE r_n_x=0., r_n_y=0., r_n_z=0., r_n, r_n2;

  // sum over real space
#ifdef MOVE_IN_X
  for(nx=-Ewald_Nx; nx<=Ewald_Nx; nx++)
#endif
  {
#ifdef MOVE_IN_Y
    for(ny=-Ewald_Nx; ny<=Ewald_Nx; ny++)
#endif
    {
#ifdef MOVE_IN_Z
      for(nz=-Ewald_Nx; nz<=Ewald_Nx; nz++)
#endif
      {
#ifdef LATTICE_TRIANGULAR
        r_n_x = Ewald_Lx_Ly_sq*(DOUBLE) nx + x;
        r_n_y = Ewald_Ly_Lx_sq*(DOUBLE) ny + y;
#else
        r_n_x = (DOUBLE) nx + x;
        r_n_y = (DOUBLE) ny + y;
        r_n_z = (DOUBLE) nz + z;
#endif
        r_n2 = r_n_x*r_n_x + r_n_y*r_n_y + r_n_z*r_n_z;
        r_n = sqrt(r_n2);
        psi += Ewald_rho(r_n,r_n2);
      }
    }
  }

  // sum over momentum space
#ifdef MOVE_IN_X
  for(nx=-Ewald_Nk; nx<=Ewald_Nk; nx++)
#endif
  {
#ifdef MOVE_IN_Y
    for(ny=-Ewald_Nk; ny<=Ewald_Nk; ny++)
#endif
    {
#ifdef MOVE_IN_Z
      for(nz=-Ewald_Nk; nz<=Ewald_Nk; nz++)
#endif
      {
#ifdef LATTICE_TRIANGULAR
        psi += Ewald_cos_coef[abs(nx)][abs(ny)]*cos(_2PI*((DOUBLE) nx*x*Ewald_Ly_Lx_sq + (DOUBLE) ny*y*Ewald_Lx_Ly_sq));
#else
        n2 = nx*nx + ny*ny + nz*nz;
        psi += Ewald_cos_coef[n2]*cos(_2PI*((DOUBLE) nx * x + (DOUBLE) ny * y + (DOUBLE) nz * z));
#endif
      }
    }
  }
  psi += Ewald_C1;

  return psi;
}

/**************************** Ewald Sum Walker *******************************/
DOUBLE EwaldSumWalker(struct Walker *Walker) {
  int i,j;
  DOUBLE E = 0;
  DOUBLE Epsi = 0;
  DOUBLE Exi = 0;
  DOUBLE x,y,z,dr[3];

  Ewald_alpha2 = Ewald_alpha*Ewald_alpha;
  Ewald_alpha3 = Ewald_alpha2*Ewald_alpha;
  Ewald_alpha4 = Ewald_alpha2*Ewald_alpha2;
  Ewald_alpha5 = Ewald_alpha4*Ewald_alpha;
  Ewald_alpha6 = Ewald_alpha3*Ewald_alpha3;

  for(i=0; i<N; i++) {
    x = Walker->x[i];
    y = Walker->y[i];
    z = Walker->z[i];
    for(j=i+1; j<N; j++) {
      dr[0] = x - Walker->x[j];
      dr[1] = y - Walker->y[j];
      dr[2] = z - Walker->z[j];
      //FindNearestImage(&dr[0], &dr[1], &dr[2]);
      Epsi += EwaldPsi(dr[0]*Ewald_L, dr[1]*Ewald_L, dr[2]*Ewald_L);
    }
  }

  Exi = 0.5*Ewald_Xi;

  // total energy
  Exi *= (double)N*Ewald_Lk;
  Epsi *= Ewald_Lk;

  E = Epsi + Exi;

  return E;
}

/**************************** Ewald Sum All Walkers **************************/
DOUBLE EwaldSumAllWalkers(void) {
  DOUBLE E=0;
  int w;

  for(w=0; w<Nwalkers; w++) E += EwaldSumWalker(&W[w]);

  return E / (double) Nwalkers;
}

/************************ Interaction Self Energy ****************************/
// for dipoles and quadrupoles only, otherwise use general forumla
#ifdef INTERACTION_DIPOLE
# define OPTIMIZED_SELF_ENERGY_SUMMATION
#endif
#ifdef INTERACTION_QUADRUPOLE
# define OPTIMIZED_SELF_ENERGY_SUMMATION
#endif

#ifndef OPTIMIZED_SELF_ENERGY_SUMMATION
DOUBLE InteractionSelfEnergy(void) {
  int nx = 0;
  int ny = 0;
  int nz = 0;
  DOUBLE Eself1, Eself2, Eself;
  int nmax;

  // single size
  Eself = 0.;
  nmax = (int) sum_over_images_cut_off_num;
#ifdef MOVE_IN_X
  for(nx=-nmax; nx<=nmax; nx++)
#endif
#ifdef MOVE_IN_Y
    for(ny=-nmax; ny<=nmax; ny++)
#endif
#ifdef MOVE_IN_Z
      for(nz=-nmax; nz<=nmax; nz++)
#endif
      if(nx!=0 || ny!=0 || nz!=0)
        Eself += InteractionEnergy(Sqrt((nx*Lx)*(nx*Lx)+(ny*Ly)*(ny*Ly)+(nz*Lz)*(nz*Lz)));
  Eself *= 0.5*(double) N; // total self energy
  //Eself /= (double) N; // energy per particle
  Eself /= (double) (N*(N-1)/2); // to be called N(N-1)/2 times
  Eself2 = Eself;

  // double size
  nmax *= 2;
  Eself = 0.;
#ifdef MOVE_IN_X
  for(nx=-nmax; nx<=nmax; nx++)
#endif
#ifdef MOVE_IN_Y
    for(ny=-nmax; ny<=nmax; ny++)
#endif
#ifdef MOVE_IN_Z
      for(nz=-nmax; nz<=nmax; nz++)
#endif
      if(nx!=0 || ny!=0 || nz!=0)
        Eself += InteractionEnergy(Sqrt((nx*Lx)*(nx*Lx)+(ny*Ly)*(ny*Ly)+(nz*Lz)*(nz*Lz)));

  Eself *= 0.5*(double) N; // total self energy
  //Eself /= (double) N; // energy per particle
  Eself /= (double) (N*(N-1)/2); // to be called N(N-1)/2 times
  Eself1 = Eself;

  // extrapolation to infinite size
  Eself = 2.*Eself1-Eself2;
  Message("  Interaction self Energy Sum Over Images: %lf (per particle)\n", Eself*0.5*(double)(N-1));
  
  //Eself = Ewald_Xi*Ewald_Lk/((double)(N-1));

  return Eself;
}
#else // (dipoles) assumes g2(r)=1 outside of the direct summation radius // Kurbakov
DOUBLE InteractionSelfEnergy(void) { // Kurbakov
  int nx,ny; // Kurbakov
  DOUBLE rr2,rr2_,Eself; // Kurbakov
  DOUBLE sum_over_images_cut_off_num_d; // Kurbakov
  Eself=0; // Kurbakov
  sum_over_images_cut_off_num_d = sum_over_images_cut_off_num * Sqrt(Lx*Ly); // Kurbakov
  for(nx=-(int)(sum_over_images_cut_off_num_d/Lx)-1;nx<=(int)(sum_over_images_cut_off_num_d/Lx)+1;nx++) // Kurbakov
    for(ny=-(int)(sum_over_images_cut_off_num_d/Ly*1.5)-1;ny<=(int)(sum_over_images_cut_off_num_d/Ly*1.5)+1;ny++){ // Kurbakov
      rr2_=nx*Lx*nx*Lx+ny*Ly*ny*Ly*0.5;rr2=nx*Lx*nx*Lx+ny*Ly*ny*Ly; // Kurbakov
      if(0<rr2_&&rr2_<sum_over_images_cut_off_num_d*sum_over_images_cut_off_num_d){ // Kurbakov
#ifdef INTERACTION_QUADRUPOLE
      Eself+=1./rr2/rr2/Sqrt(rr2); // for quadrupoles // Kurbakov
#endif
#ifdef INTERACTION_DIPOLE
      Eself+=1./rr2/Sqrt(rr2); // for dipoles // Kurbakov
#endif
      } // Kurbakov
  } // Kurbakov

#ifdef INTERACTION_QUADRUPOLE
  Eself+=2.*PI/3./sum_over_images_cut_off_num_d/sum_over_images_cut_off_num_d/sum_over_images_cut_off_num_d/Lx/Ly*0.66312316693288855; // for quadrupoles // Kurbakov
#endif
#ifdef INTERACTION_DIPOLE
  Eself+=2.*PI/sum_over_images_cut_off_num_d/Lx/Ly*0.85984660010223786; // Kurbakov
#endif

  Eself *= 0.5*(double) N; // total self energy // Kurbakov
  Eself /= (double) (N*(N-1)/2); // to be called N(N-1)/2 times // Kurbakov
  return Eself; // Kurbakov
} // Kurbakov
#endif // Kurbakov

/***************** Interaction Self Energy Sum Over Images *******************/
// direct sum of the potential energy over images
#ifndef OPTIMIZED_SELF_ENERGY_SUMMATION
DOUBLE InteractionEnergySumOverImages(DOUBLE x, DOUBLE y, DOUBLE z) {
  DOUBLE E = 0.;
  static int initialized = OFF;
  static DOUBLE Eself = 0.;
  int nx, ny, nz;

  nx = ny = nz = 0;

  if(initialized == OFF) { // first time calculate self energy
    Eself = InteractionSelfEnergy();
    initialized = ON;
  }

#ifdef MOVE_IN_X
  for(nx=-(int)sum_over_images_cut_off_num; nx<=(int)sum_over_images_cut_off_num; nx++)
#endif
#ifdef MOVE_IN_Y
    for(ny=-(int)sum_over_images_cut_off_num; ny<=(int)sum_over_images_cut_off_num; ny++)
#endif
#ifdef MOVE_IN_Z
      for(nz=-(int)sum_over_images_cut_off_num; nz<=(int)sum_over_images_cut_off_num; nz++)
#endif
      E += InteractionEnergy(Sqrt((x+nx*Lx)*(x+nx*Lx)+(y+ny*Ly)*(y+ny*Ly)+(z+nz*Lz)*(z+nz*Lz)));

  E += Eself;

  return E;
}
#else // (dipoles) assumes g2(r)=1 outside of the direct summation radius
DOUBLE InteractionEnergySumOverImages(DOUBLE x, DOUBLE y, DOUBLE z) {
  int nx,ny;
  DOUBLE rr2,rr2_,rrx,rry,rr2x,rr2y,s; // Kurbakov
  static DOUBLE sum_over_images_cut_off_num_d;
  static int initialized = OFF;
  static DOUBLE Eself,Etail;
  static int Nnx,Nny;
  static DOUBLE Lnx,Lny;

  if(initialized == OFF) { // first time calculate self energy

    initialized = ON; // Kurbakov
    if(sum_over_images_cut_off_num < 1.2) { // Kurbakov
      Warning(" Cannot do summation with sum_over_images_cut_off_num<1.2, setting it to 1.2\n"); // Kurbakov
      sum_over_images_cut_off_num = 1.2; // Kurbakov
    } // Kurbakov
    Eself = InteractionSelfEnergy(); // Kurbakov
    sum_over_images_cut_off_num_d = sum_over_images_cut_off_num * Sqrt(Lx*Ly);
    Nnx=(int)(sum_over_images_cut_off_num_d/Lx)+1;Lnx=-Nnx*Lx; // Kurbakov
    Nny=(int)(sum_over_images_cut_off_num_d/Ly*1.5)+1;Lny=-Nny*Ly; // Kurbakov
#ifdef INTERACTION_QUADRUPOLE
  Etail=2.*PI*0.66312316693288855/(3.*sum_over_images_cut_off_num_d*sum_over_images_cut_off_num_d*sum_over_images_cut_off_num_d*Lx*Ly)+Eself; // for quadrupoles // Kurbakov
#endif
#ifdef INTERACTION_DIPOLE
  Etail=2.*PI*0.85984660010223786/(sum_over_images_cut_off_num_d*Lx*Ly)+Eself; // Kurbakov
#endif
  }
  s=0;
  for(rrx=Lnx+x,nx=-Nnx; nx<=Nnx; nx++,rrx+=Lx) { // Kurbakov
    rr2x=rrx*rrx; // Kurbakov
    for(rry=Lny+y,ny=-Nny; ny<=Nny; ny++,rry+=Ly) { // Kurbakov
    rr2y=rry*rry; // Kurbakov
      rr2_=rr2x+rr2y*0.5; // Kurbakov
      if(rr2_<sum_over_images_cut_off_num_d*sum_over_images_cut_off_num_d) { // Kurbakov
        rr2=rr2x+rr2y; // Kurbakov
#ifdef INTERACTION_QUADRUPOLE
        s += 1./(rr2*rr2*Sqrt(rr2)); // for quadrupoles // Kurbakov
#endif
#ifdef INTERACTION_DIPOLE
        s += 1./(rr2*Sqrt(rr2)); // Kurbakov
#endif
      } // Kurbakov
    } // Kurbakov
  } // Kurbakov
  return(s+Etail);
}
#endif

/**************************** Ewald Sum Pair *********************************/
// contribution to the energy from a single pair of particles plus self energy
//DOUBLE EwaldSumPair(DOUBLE x, DOUBLE y, DOUBLE z) {
//  return (EwaldPsi(x*Ewald_L, y*Ewald_L, z*Ewald_L) + Ewald_Xi/(double) (N-1))*Ewald_Lk;;
//}

/**************************** Ewald Optimize *********************************/
// returns best accuracy that can be achieved with given Nx and Nk
DOUBLE EwaldOptimize(int Nx, int Nk) {
  DOUBLE Ewald_alpha_store;
  DOUBLE E, sigma;
  DOUBLE sigma_min = 1, alpha_opt;
  clock_t time_check;

  Ewald_alpha_store = Ewald_alpha;

  Ewald_Nx = Nx;
  Ewald_Nk = Nk;

  Message("# Nx = %i Nk = %i (1-alpha, 2-Nx, 3-Nk, 4-E, 5-sigma, 6-time)\n", Nx, Nk);
  //Message("alpha sigma\n");
  for(Ewald_alpha = Ewald_alpha_store*0.5; Ewald_alpha<Ewald_alpha_store*2; Ewald_alpha += Ewald_alpha_store*0.1) {
    Ewald_alpha2 = Ewald_alpha*Ewald_alpha;
    Ewald_alpha3 = Ewald_alpha2*Ewald_alpha;
    Ewald_alpha4 = Ewald_alpha2*Ewald_alpha2;
    Ewald_alpha5 = Ewald_alpha4*Ewald_alpha;
    Ewald_alpha6 = Ewald_alpha3*Ewald_alpha3;

    EwaldInitCoef();

    time_check = -clock();
    EwaldOptimizeSigma(&E, &sigma);
    time_check += clock();

    Message("%e %i %i %e %e %e\n", Ewald_alpha, Nx, Nk, E/(double)N, sigma/Ewald_Exact_mean, (DOUBLE)time_check/(DOUBLE)CLOCKS_PER_SEC);

    if(sigma/Ewald_Exact_mean<sigma_min) {
      sigma_min = sigma/Ewald_Exact_mean;
      alpha_opt = Ewald_alpha;
    }
  }
  Message("\n\n");

  Ewald_alpha = Ewald_alpha_store;
  Ewald_alpha2 = Ewald_alpha*Ewald_alpha;
  Ewald_alpha3 = Ewald_alpha2*Ewald_alpha;
  Ewald_alpha4 = Ewald_alpha2*Ewald_alpha2;
  Ewald_alpha5 = Ewald_alpha4*Ewald_alpha;
  Ewald_alpha6 = Ewald_alpha3*Ewald_alpha3;
  //Message("best sigma = %lf with alpha = %lf\n", sigma_min, alpha_opt);

  return alpha_opt;
}

/**************************** Ewald Optimize *********************************/
void EwaldOptimizeInitExact(void) {
  int w;
  int Nk_store, Nx_store;

  Nx_store = Ewald_Nx;
  Nk_store = Ewald_Nk;
  Ewald_Nx = Ewald_Nk_max;
  Ewald_Nk = Ewald_Nk_max;

  for(w=0; w<Nwalkers; w++) Ewald_Exact[w] = EwaldSumWalker(&W[w]);

  Ewald_Nx = Nx_store;
  Ewald_Nk = Nk_store;
}

/**************************** Ewald Optimize *********************************/
void EwaldOptimizeSigma(DOUBLE *Emean, DOUBLE *sigma) {
  int w;
  DOUBLE E;

  *sigma = 0.;
  *Emean = 0.;

  for(w=0; w<Nwalkers; w++) {
    E = EwaldSumWalker(&W[w]);
    *sigma += (E-Ewald_Exact[w])*(E-Ewald_Exact[w]);
    *Emean += E;
  }
  *sigma = sqrt(*sigma/(double) Nwalkers);
  *Emean /= (double) Nwalkers;
}

/******************************* Smooth cutoff *******************************/
#ifdef INTERACTION_SUM_OVER_IMAGES_SMOOTH_CUTOFF
double Etail_smooth_cutoff;
double Rc_smooth_cutoff2, _Rc_smooth_cutoff;
double Rrx_smooth_cutoff, Rry_smooth_cutoff, Rrz_smooth_cutoff;
int Nx_smooth_cutoff, Ny_smooth_cutoff, Nz_smooth_cutoff;

double SmoothCutoffSum(double x, double y, double z) {
// must be -Lx/2 <= x < Lx/2, -Ly/2 <= y < Ly/2 and -Lz/2 <= z < Lz/2
// for 2D case one must be z=0
  CaseX(int nx; double rrx;)
  CaseY(int ny; double rry;)
  CaseZ(int nz; double rrz;)
  double rr2,rr,th,s=0.;

  // calculation of energy (self-enegry corresponds to x=y=z=0)
  CaseZ(for(nz=-Nz_smooth_cutoff, rrz=z-Rrz_smooth_cutoff; nz<=Nz_smooth_cutoff; nz++, rrz+=Lz)) 
  {
    CaseY(for(ny=-Ny_smooth_cutoff, rry=y-Rry_smooth_cutoff; ny<=Ny_smooth_cutoff; ny++, rry+=Ly)) 
    {
      CaseX(for(nx=-Nx_smooth_cutoff, rrx=x-Rrx_smooth_cutoff; nx<=Nx_smooth_cutoff; nx++, rrx += Lx)) 
      {
        rr2 = CaseX(rrx*rrx) CaseY(+ rry*rry) CaseZ(+ rrz*rrz);
        if(rr2 > 0. && rr2 < Rc_smooth_cutoff2) {
          // calculation of smoothing
          rr = sqrt(rr2);
          th = 1. - rr*_Rc_smooth_cutoff;
          th = 1. - th*th;
          th *= th;
          th *= th;
          th *= th*th;
          th = 1. - th;
          th *= th;
          // calculation of interaction
          s += th*InteractionEnergy(rr);// general case
        }
      }
    }
  }

#ifdef INTERACTION_CALOGERO
  return D*(s+Etail_smooth_cutoff);
#else
  return (s+Etail_smooth_cutoff);
#endif
}
#endif

int SmoothCutoffInit(void) {
#ifdef INTERACTION_SUM_OVER_IMAGES_SMOOTH_CUTOFF
  double u,rr,rr2,deriv,th,s,s1,s2,th1,th2,RC,n_D;
  int k,km,A,a,jx,jy,jz,Jx,Jy,Jz,nn2;
  static int first_time = ON;

  Message("  Smooth cut off check:\n");
  if(first_time == ON) { // calculation of tail
#ifdef TRIAL_1D
    RC=100.;
    n_D= 1./n;
#endif
#ifdef TRIAL_2D
    RC=40.;
    n_D = 1./sqrt(n);
#endif
#ifdef TRIAL_3D
    RC=25.;
    n_D = 1./pow(n, 1./3.);
#endif
    Message("  Automatically converting  Rc_smooth_cutoff %lf -> %lf\n", Rc_smooth_cutoff, Rc_smooth_cutoff*n_D);
    Rc_smooth_cutoff *= n_D;
    RC *= n_D;
    first_time = OFF;
  }

  Rc_smooth_cutoff2=Rc_smooth_cutoff*Rc_smooth_cutoff;
  _Rc_smooth_cutoff=1/Rc_smooth_cutoff;

  CaseX(Nx_smooth_cutoff=(int)(Rc_smooth_cutoff*L_inv_x+0.5);)
  CaseY(Ny_smooth_cutoff=(int)(Rc_smooth_cutoff*L_inv_y+0.5);)
  CaseZ(Nz_smooth_cutoff=(int)(Rc_smooth_cutoff*L_inv_z+0.5);)

  CaseX(Rrx_smooth_cutoff=Nx_smooth_cutoff*Lx;)
  CaseY(Rry_smooth_cutoff=Ny_smooth_cutoff*Ly;)
  CaseZ(Rrz_smooth_cutoff=Nz_smooth_cutoff*Lz;)
/*
  // calculation of tail
#ifdef TRIAL_2D
#ifdef INTERACTION_COULOMB
  Etail_smooth_cutoff=-4.283844763680682*Rc_smooth_cutoff/Lx/Ly;// Coulomb
#endif
#ifdef INTERACTION_CALOGERO
  Etail_smooth_cutoff=(2.498663137395704-log(Rc_smooth_cutoff)*2.*PI)/Lx/Ly;// Calogero
#endif
#ifdef INTERACTION_DIPOLE
  Etail_smooth_cutoff=9.50400881303096/Rc_smooth_cutoff/Lx/Ly;// dipoles
#endif
#ifdef INTERACTION_R5
  Etail_smooth_cutoff=8.126834372198/Rc_smooth_cutoff/Rc_smooth_cutoff/Rc_smooth_cutoff/Lx/Ly;// quardupoles
#endif
#endif // TRIAL_2D

#ifdef TRIAL_3D
#ifdef INTERACTION_COULOMB
  Etail_smooth_cutoff=-2.999820701306844*Rc_smooth_cutoff*Rc_smooth_cutoff/Lx/Ly/Lz; // Coulomb
#endif
#ifdef INTERACTION_CALOGERO
  Etail_smooth_cutoff=-8.567689527361368*Rc_smooth_cutoff*L_inv_x*L_inv_y*L_inv_z; // Calogero
#endif
#ifdef INTERACTION_RYDBERG
  Etail_smooth_cutoff=16.25366874439/Rc_smooth_cutoff/Rc_smooth_cutoff/Rc_smooth_cutoff/Lx/Ly/Lz; // Rydberg
#endif
#endif // TRIAL_3D
*/

// calculation of tail in a general case of interaction
  km=256;A=6;
  s1=0.;
  for(k=0;k<km;k++){
    u=(k+.5)/km*PI/2;
    rr=RC;for(a=0;a<A;a++)rr*=tan(u);
    deriv=rr*2*A/sin(2*u);
    if(rr<RC)th = 1.-rr/RC;else th=0;
    th = 1. - th*th;
    th *= th;
    th *= th;
    th *= th;
    th *= th*th;
    th = 1. - th;
    th *= th;
#ifdef TRIAL_1D
    if(th<1.)s1+= InteractionEnergy(rr)*(1.-th)*2.*deriv;
#endif
#ifdef TRIAL_2D
    if(th<1.)s1+= InteractionEnergy(rr)*(1.-th)*2.*pi*rr*deriv;
#endif
#ifdef TRIAL_3D
    if(th<1.)s1+= InteractionEnergy(rr)*(1.-th)*4.*pi*rr*rr*deriv;
#endif
  }
  s1*=PI/2/km;km*=2;
  s2=0.;
  for(k=0;k<km;k++){
    u=(k+.5)/km*PI/2;
    rr=RC;for(a=0;a<A;a++)rr*=tan(u);
    deriv=rr*2*A/sin(2*u);
    if(rr<RC)th = 1.-rr/RC;else th=0;
    th = 1. - th*th;
    th *= th;
    th *= th;
    th *= th;
    th *= th*th;
    th = 1. - th;
    th *= th;
#ifdef TRIAL_1D
    if(th<1.)s2+= InteractionEnergy(rr)*(1-th)*2.*deriv;
#endif
#ifdef TRIAL_2D
    if(th<1.)s2+= InteractionEnergy(rr)*(1-th)*2.*pi*rr*deriv;
#endif
#ifdef TRIAL_3D
    if(th<1.)s2+= InteractionEnergy(rr)*(1-th)*4.*pi*rr*rr*deriv;
#endif
  }
  s2*=PI/2./km;km/=2;
  if(fabs(s1/s2-1.)<1E-4){
    Jx=Jy=Jz=(int)(RC/n_D+.5);s=s2;
#ifdef TRIAL_1D
    Jy=Jz=0;
#endif
#ifdef TRIAL_2D
    Jz=0;
#endif
    for(jx=-Jx;jx<=Jx;jx++)for(jy=-Jy;jy<=Jy;jy++)for(jz=-Jz;jz<=Jz;jz++)if(jx||jy||jz){
      rr=sqrt(jx*jx+jy*jy+jz*jz+0.)*n_D;
      if(rr<Rc_smooth_cutoff)th = 1.-rr/Rc_smooth_cutoff;else th=0;
      th = 1. - th*th;
      th *= th;
      th *= th;
      th *= th*th;
      th = 1. - th;
      th *= th;
      th1 = th;
      if(rr<RC)th = 1.-rr/RC;else th=0;
      th = 1. - th*th;
      th *= th;
      th *= th;
      th *= th;
      th *= th*th;
      th = 1. - th;
      th *= th;
      th2 = th;
      s += InteractionEnergy(rr)*(th2-th1)/n;
    }
  }else{
    Jx=Jy=Jz=(int)(Rc_smooth_cutoff/n_D+1);s=0;
#ifdef TRIAL_1D
    Jy=Jz=0;
#endif
#ifdef TRIAL_2D
    Jz=0;
#endif
    for(jx=-Jx;jx<=Jx;jx++)for(jy=-Jy;jy<=Jy;jy++)for(jz=-Jz;jz<=Jz;jz++)if(jx||jy||jz){
      rr=sqrt(jx*jx+jy*jy+jz*jz+0.)*n_D;
      if(rr<Rc_smooth_cutoff)th = 1.-rr/Rc_smooth_cutoff;else th=0;
      th = 1. - th*th;
      th *= th;
      th *= th;
      th *= th*th;
      th = 1. - th;
      th *= th;
      s -= InteractionEnergy(rr)*th/n;
    }
  }
  s*=N/(N-1.);
  Etail_smooth_cutoff=s;
#ifdef TRIAL_1D
  Etail_smooth_cutoff/=Lx;
#endif
#ifdef TRIAL_2D
  Etail_smooth_cutoff/=Lx*Ly;
#endif
#ifdef TRIAL_3D
  Etail_smooth_cutoff/=Lx*Ly*Lz;
#endif

// calculation of self-energy (x=y=z=0)
  CaseX(int nx; double rrx;)
  CaseY(int ny; double rry;)
  CaseZ(int nz; double rrz;)
  s=0.;
  CaseZ(for(nz=-Nz_smooth_cutoff, rrz=0-Rrz_smooth_cutoff; nz<=Nz_smooth_cutoff; nz++, rrz+=Lz)) 
  {
    CaseY(for(ny=-Ny_smooth_cutoff, rry=0-Rry_smooth_cutoff; ny<=Ny_smooth_cutoff; ny++, rry+=Ly)) 
    {
      CaseX(for(nx=-Nx_smooth_cutoff, rrx=0-Rrx_smooth_cutoff; nx<=Nx_smooth_cutoff; nx++, rrx += Lx)) 
      {
        rr2 = CaseX(rrx*rrx) CaseY(+ rry*rry) CaseZ(+ rrz*rrz);
        nn2 = CaseX(nx*nx) CaseY(+ny*ny) CaseZ(+nz*nz);
        if(nn2 && rr2 < Rc_smooth_cutoff2) {
          rr = sqrt(rr2);
          th = 1. - rr*_Rc_smooth_cutoff;
          th = 1. - th*th;
          th *= th;
          th *= th;
          th *= th*th;
          th = 1. - th;
          th *= th;
          s += th*InteractionEnergy(rr);
        }
      }
    }
  }
  Etail_smooth_cutoff += s/(N-1.);

  Message("  tail energy: (E/N)_tail=%g\n", Etail_smooth_cutoff);
#endif
  return 0;
}