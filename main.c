/*main.c written by G.E. Astrakharchik 2001-2017*/

#include <stdio.h>
#include <time.h>
#include "main.h"
#include "randnorm.h"
#include "memory.h"
#include "utils.h"
#include "rw.h"
#include "vmc.h"
#include "quantities.h"
#include "dmc.h"
#include "gencoord.h"
#include "display.h"
#include "crystal.h"
#include "spline.h"
#include "speckles.h"
#include "optimiz.h"
#include "ewald.h"
#include "pimc.h"
#include "pigs.h"
#include "tvmc.h"
#include "compatab.h"
#include MATHINCLUDE

#ifdef MPI
#  include <mpi.h>
#  include "parallel.h"
#endif
#ifdef _OPENMP
#  include <omp.h>
#  include "parallel.h"
#endif

int grid_trial = 1000;
int gridOBDM = 601;
int gridOBDM_MATRIX = 60;
int gridRD = 601;
int gridPD = 601;
int gridg3 = 601;
int gridRDx = 1;
int gridRDy = 1;
int grid_pure_block  = 100;
int gridSk = 10;
int gridPD_MATRIX_x = 10;
int gridPD_MATRIX_y = 10;
int gridSKT_t = 100;
int gridSKT_k = 30;
DOUBLE SKT_dk = 0;
DOUBLE SKT_kmin = 0; 

unsigned long int Niter=100, Nmeasure=10, iteration_global=0;
DOUBLE Rpar=1, Apar=1, Bpar=20, Cpar=0.1, Dpar = 1, R0par =1, Epar = 0, Ipar = 0, Kpar = 1, Mpar = 1, alpha_R = -PI, beta_z = 1.;
DOUBLE Kpar11=1, Kpar12=1, aA=1;
DOUBLE alpha_x=0., alpha_x2, two_alpha_x, alpha_Rx= -PI;
DOUBLE alpha_y=0., alpha_y2, two_alpha_y, alpha_Ry= -PI;
DOUBLE alpha_z=0., alpha_z2, two_alpha_z, alpha_Rz= -PI;
DOUBLE omega_x = 1, omega_y = 1, omega_z = 1; // parameters of the trap
DOUBLE omega_x2, omega_y2, omega_z2;
DOUBLE alpha_latt = 1, alpha_Jastrow = 1, beta_latt = 1, gamma_latt = 1;
DOUBLE D = 1, n = 0.01, dt = 1e-3, dt_vmc = 1e-3, dt_one = 1., dt_all = 1e-3, beta = 1.;
DOUBLE acceptance_rate = 0.;
int N=16, Nimp=0, NCrystal, Ndens=16, Nspin = 1;
int Nwalkers=100, NwalkersMax=100;
int Nwalkers_pure = 0; // number of walkers used for pure OBDM measurement
DOUBLE Nwalkersw=100., Nwalkers_weight_killed = 0.;
int blck=10, blck_heating=0;
int verbosity=0;
DOUBLE E, EFF, Eo, Epot, Ekin, Elat, Eck, Edamping, Eint, Eext;
DOUBLE Srecoil=1, Syrecoil=1;
int Nlattice, Nlattice3;
DOUBLE L, Lhalf, L2, Lhalf2, Lmax, Lm1, Lwf, Lhalfwf, Lhalfwf2, Lcutoff_pot, Lcutoff_pot2;
DOUBLE Lx, Ly, Lz, L_half_x, L_half_y, L_half_z, L_inv_x, L_inv_y, L_inv_z, kL;
DOUBLE a=1, a2, b=1;
DOUBLE lambda=1, lambda2, lambda4;
DOUBLE Rzigzag=0.;
DOUBLE T=0.; 
DOUBLE Tannealing = 1;
DOUBLE PD_MATRIX_x = 1.;
DOUBLE PD_MATRIX_y = 1.;
DOUBLE NkMaxTrap = 10;
DOUBLE inwf_power = 1.; // load density.in file and change the power to inwf_power
DOUBLE m_mu = 0.5; // reduced mass
DOUBLE Rsphere = 1., Rnow; // sphere radius size
int generate_crystal_coordinates = ON; 

char file_particles[25] = "in3Dprev.in";
char file_energy[25] = "oute.dat";
char file_OBDM[25] = "outdr.dat";
char file_OBDM_MATRIX[25] = "outobdm.dat";
char file_PD[25] = "outpd.dat";
char file_PD_pure[25] = "outpdp.dat";
char file_PDz[25] = "outpdz.dat";
char file_RD[25] = "outrd.dat";
char file_RD_pure[25] = "outrdp.dat";
char file_RDz[25] = "outrdz.dat";
char file_RDz_pure[25] = "outrdzp.dat";
char file_R2[25] = "outr2.dat";
char file_z2[25] = "outz2.dat";
char file_R2_pure[25] = "outr2pur.dat";
char file_z2_pure[25] = "outz2pur.dat";
char file_wf[25]= "inwf.in";
char file_wf1b[25]= "inwf1b.in";
char file_Sk[25] = "outsk.dat";
char file_Sk_pure[25] = "outskp.dat";
char file_SD[25] = "outsd.dat";

int measure_energy = ON;
int measure_OBDM = OFF;
int measure_TBDM = OFF;
int measure_OBDM_MATRIX = OFF;
int measure_SD = OFF;
int measure_RadDistr = ON;
int measure_RadDistrMATRIX = OFF;
int measure_PairDistrMATRIX = OFF;
int measure_PairDistr = ON;
int measure_PairDistr_R2 = OFF;
int measure_Sk = OFF;
int measure_R2 = OFF;
int measure_OP = OFF;
int measure_Lind = OFF;
int measure_Nk_pure = OFF;
int measure_energy_barrier = OFF;
int measure_SkMATRIX = OFF;
int measure_FormFactor = OFF;
int video = OFF;
int branchng_present = ON;
int square_box = ON;
int file_particles_append = OFF;
int file_append = OFF;
int optimization = OFF;
int measure_Lindemann_Lozovik = OFF;
int measure_wavefunction_projection = OFF;
int measure_RDz = OFF;
int measure_g3 = OFF;
int R2_subtract_CM = OFF; // used for the calculation of the breathing mode
int OBDM_subtract_CM = OFF; 
int measure_effective_potential = OFF; // effective potential with respect to center of mass position
int measure_pure_coordinates = OFF;
int measure_tVMC = OFF;
int measure_Hessian_matrix = OFF;
int check_kinetic_energy = ON;

struct Walker *W, *Wp;
struct Grid G; // Jastrow terms of w.f.
struct Grid G1; // One=body terms of w.f.
struct Grid GV; // Vext(z)

struct sOBDM OBDM, OBDMtrap, OBDMfermi, TBDM, TBDMfermi;
struct sSk Sk, Sk_pure;
struct sMATRIX OBDM_MATRIX, PD_MATRIX, RD_MATRIX, TBDM_MATRIX, OBDMfermi_MATRIX, TBDMfermi_MATRIX, SKT, SKTre, SKTim, PhiTau, Nk_MATRIX, Nkfermi_MATRIX;
struct Distribution PD, PDz, PD_pure, RD, RDz, RD_pure, RDz_pure, HR, HR_pure, Veff;
struct sSD SD;
struct sOrderParameter OrderParameter;
struct sLindemannRatio LindemannRatio;
struct sPureEstimator Epot_pure;

DOUBLE *u_mi, *u_ij;
DOUBLE *u_tbdm_i, **u_tbdm_ij, *u_tbdm_mcmillan_i, **u_tbdm_mcmillan_ij;
DOUBLE *order;
int *sign_u_mi,*sign_u_ij;

DOUBLE *mu_k, *mu_p_k, *mu_pp_k;
DOUBLE lattice_length = 0;
DOUBLE Etail;
DOUBLE Vo = 1; // Square Well height
DOUBLE Ro = 1; // Square Well width
DOUBLE delta_tilde = 1; // parameter of Rydberg potential
DOUBLE solitonV2 = 0.1; // soliton w.f. parameters
DOUBLE solitonXi = 1;
DOUBLE solitonk;
DOUBLE soliton_fixed_phase_V2 = 0.1; // soliton effective potential parameters
DOUBLE soliton_fixed_phase_Xi = 1;
DOUBLE R3 = 1; // hyperradius 

int accepted = 0;
int rejected = 0;
int accepted_one = 0;
int rejected_one = 0;
int accepted_all = 0;
int rejected_all = 0;
int overlaped = 0;
int tried = 0;
int Nmeasured = 0;
int McMillan_points = 1;
int McMillanTBDM_points = 100;
int Npop = 50, Npop_max, Npop_min, Npop_virtual = 20;
int generate_new_coordinates = 0;
int generate_new_coordinates_from_lattice = 0;
int R2times_measured = 0;
DOUBLE Evar = 0,  Evar2 = 0.;
DOUBLE EFFvar=0., EFFvar2 = 0.;
DOUBLE Epotvar=0., Epotvar2 = 0.;
DOUBLE Ekinvar=0., Ekinvar2 = 0.;

DOUBLE xR2 = 0., zR2 = 0.;

DOUBLE reduce = 0.9, amplify = 1.1;

clock_t time_start = 0;
clock_t time_end = 0;
clock_t time_simulating = 0;
clock_t time_measuring = 0;
clock_t time_communicating = 0;
time_t time_era_begin, time_era_end;

int MC = ON;
int SmartMC = 0;
int *walkers_storage;        /* contains indices to the alive walkers */
int *dead_walkers_storage;   /* contains indices to empty walkers */
int Ndead_walkers = 0;
DOUBLE *branching_weight; // [NwalkersMax]
DOUBLE *branching_xi;
int *branching_multiplicity;
int Nscal=1; // replicate system times in each direction
int Ncells=1; //number of cells 
DOUBLE Lcell, Lwf_store;
int Nwf_store;
DOUBLE sum_over_images_cut_off_num = 2.;
DOUBLE Rc_smooth_cutoff = 100.;
int Niter_store = 100; // store configuration for reweighting
DOUBLE bilayer_width = 0, bilayer_width2 = 0;
DOUBLE Rpar11, Rpar12 = 0.77;
DOUBLE Apar12, Bpar12, Cpar12;
DOUBLE t_tunneling = 0.;
DOUBLE Uo = 1; // Vext(r) = - Uo / r^2
DOUBLE gaussian_alpha = 1;
DOUBLE gaussian_beta = -1;
DOUBLE RoSW = 1.; // square-well size
DOUBLE Epotcutoff = -1; // cut-off energy in solving numerically the two-body problem
DOUBLE cR[10], cI[10], tvmcNpar = 3, tvmcNobs = 5; // tVMC variational parameters, real and imaginary parts
DOUBLE CMseparation = 1;
DOUBLE hard_core_diameter = 0.;
DOUBLE dt_Kurbakov, Kurbakov_softcore = 0.;
int Kurbakov_BC = 0;
DOUBLE one_body_momentum_i = 0;
DOUBLE one_body_momentum_k;
DOUBLE E_from_growth_estimator;
DOUBLE sum_of_two_wavefunctions_parameter1 = 1.;
DOUBLE sum_of_two_wavefunctions_parameter2 = 1.;
DOUBLE sum_of_two_wavefunctions_alpha = 0.5;
DOUBLE numeric_derivative_dr = 0.001;

#ifdef TRIAL_TWO_BODY_SCATTERING_NUMERICAL_SPINFULL
  struct Grid **GridSpin; // [spin1][spin2] array
#endif

int main(int argn, char **argv) {
  DOUBLE sigma, error;
  DOUBLE Eblck = 0;
  int w,i;
  int block, iter;
  unsigned long int measured_blck = 0;
  clock_t time_begin;

  time_start = clock();
  time(&time_era_begin);

  Message(SRC_CODE_DATE);
#ifdef MPI
  ParallelInitialize(argn, argv);
  ParallelRandInit();
#else
#ifdef _OPENMP
  ParallelInitialize(argn, argv);
  ParallelRandInit();
#else
  RandInit();
#endif
#endif

#ifdef SECURE
  Warning("Secure/Debug option enabled\n");
#endif
#ifdef _DEBUG
  Warning("  Code was compiled with DEBUG option and might be not efficient.\n\n"); 
#ifdef __INTEL_COMPILER
  Message("  Intel compiler is used.\n");
#endif
#ifdef _MSC_VER
  Message("  Microsoft visual studio compiler is used.\n");
#endif
#endif
#ifdef DATATYPE_FLOAT
  Message("  Goal precision: float (8 digits after comma)\n");
#endif
#ifdef DATATYPE_DOUBLE
  Message("  Goal precision: double (16 digits after comma)\n");
#endif
#ifdef DATATYPE_LONG_DOUBLE
  Message("  Goal precision: long double (32 digits after comma)\n");
#endif

#ifdef BRANCHING_ALWAYS_STABLE
  Message("  Branching version: avoids system collapse\n");
#else
  Message("  Branching version: permits system collapse\n");
#endif

  CheckMantissa();

#ifdef IMAGES_FLOOR
  Message("floor()");
#endif
#ifdef IMAGES_IF_IF
  Message("if() if()");
#endif
#ifdef IMAGES_ABS_IF
  Message("abs() if()");
#endif

#ifdef CALCULATE_DERIVATIVES_NUMERICALLY
  Message("  Derivatives will be calculated numerically (kinetic energy, drift force)\n");
#endif
#ifdef CALCULATE_DERIVATIVES_ANALYTICALLY
  Message("  Derivatives will be calculated analytically (kinetic energy, drift force)\n");
#endif

  LoadParameters();

  Message("\n\n  particles move in ");
  CaseX(Message("X ");)
  CaseY(Message("Y ");)
  CaseZ(Message("Z ");)
  Message("directions\n");
  if(measure_SD && MC == DIFFUSION) {
    Message("  superfluid response is measured in ");
    CaseXSD(Message("X ");)
    CaseYSD(Message("Y ");)
    CaseZSD(Message("Z ");)
    Message("directions\n");
  }

  if(energy_unit != 1.) Warning(" unit of energy is different from 1 and equals %" LG "\n", energy_unit);

  if(boundary == NO_BOUNDARY_CONDITIONS && L<0.) {
    L = 100;
  }
  else if(boundary == ONE_BOUNDARY_CONDITION) {
    L = (DOUBLE) Ndens / n;
  }
  else if(boundary == TWO_BOUNDARY_CONDITIONS) {
#ifdef BC_2DPBC_HEXAGON
    Ly = 2.*sqrt((DOUBLE) (Ndens) / n *2./(3.*sqrt(3.)));
    Lx = sqrt(3.)/2.*Ly;
    Message("  Box size: Lx= %lf, Ly= %lf\n", Lx, Ly);
    L = Lx;
#else // cube
    L = sqrt((DOUBLE) (Ndens) / n);
#endif
  }
  else if(boundary == THREE_BOUNDARY_CONDITIONS) {
#ifdef BC_3DPBC_TRUNCATED_OCTAHEDRON
    L = pow((DOUBLE) (2*Ndens) / n, (1./3.));
#else // cube
    //L = pow((DOUBLE) (Ndens) / n, (1./3.));
    L = exp(log(Ndens / n) / 3.);
#endif
  }
  Eck = 1e10;

#ifdef CENTER_OF_MASS_IS_NOT_MOVED
  Message("  Center of mass is not moved\n");
#endif
#ifdef CENTER_OF_MASS_Z_IS_NOT_MOVED
  Message("  Center of mass Z is not moved\n");
#endif
#ifdef CENTER_OF_MASS_IS_NOT_MOVED_TWO_PARTS
  Message("  Center of mass is not moved for the first and second halves of the system\n");
#endif
#ifdef FIRST_PARTICLE_IS_NOT_MOVED
  Message("  Position of the first particle is not moved\n");
#endif

#ifdef MEASURE_FORM_FACTOR_FOR_SUM_RULES
  if(measure_FormFactor) Message("  Form Factor: S(tau) = (1/N) <A(tau)A(0)> with A = (x_i)^k  operators\n");
#else
  if(measure_FormFactor) Message("  Form Factor: S(k,tau) = (1/N) <(exp(ikr_i(tau)) (exp(-ikr_j(0))>, i.e. density operator\n");
#endif

  Message("  Interaction potential (not included in Eloc part): ");
#ifdef INTERACTION_ABSENT
  Message("ABSENT\n");
#endif
#ifdef INTERACTION_DIPOLE
  Message("DIPOLE\n");
#endif
#ifdef INTERACTION_QUADRUPOLE
  Message("QUADRUPOLE\n");
#endif
#ifdef INTERACTION_COULOMB
  Message("COULOMB\n");
#endif
#ifdef INTERACTION_COULOMB_1D
  Message("COULOMB 1D\n");
#endif
#ifdef INTERACTION_LINEAR_1D_COULOMB
  Message("LINEAR COULOMB 1D (|x|)\n");
#endif
#ifdef INTERACTION_LENNARD_JONES
  Message("LENNARD JONES\n");
#endif
#ifdef INTERACTION_YUKAWA
  Message("YUKAWA + CORRECTION\n");
#endif
#ifdef INTERACTION_YUKAK0
  Message("Ko(r) function\n");
#endif
#ifdef INTERACTION_PETROV
  Message("'Petrov' function\n");
#endif
#ifdef INTERACTION_PETROV_CONST
  Message("'Petrov' + const function\n");
#endif
#ifdef INTERACTION_PETROV_D
  Message(" D*2./r*Exp(-2.*r)*(1. - 0.5/r) - 'Petrov 3D' function\n");
#endif
#ifdef INTERACTION_SQUARE_WELL
  Message("Square well\n");
  ConstructSqWellDepth();
  FindSWBindingEnergy();
#endif
#ifdef INTERACTION_SPINFULL_DIPOLE_BILAYER
  Message(" Spinfull dipole bilayer: 1/r^3 for the same spin, (r*r-2.*bilayer_width2)*pow(r*r+bilayer_width2, -2.5) for different\n"); 
#endif
#ifdef INTERACTION_SPINFULL_COULOMB_BILAYER
  Message(" Spinfull Coulomb bilayer: 1/r for the same spin, -1/sqrt(r^2-bilayer_width2) for different\n"); 
#endif
#ifdef SPINFULL_TUNNELING
  Message("  Tunneling is allowed\n");
#endif
#ifdef INTERACTION_HYDROGEN
  Message("  Spline interpolation will be used for Vint(r), loading from inVint.in file\n");
  LoadInteractionPotential(&GV, "inVint.in");
  SplineConstruct(&GV, 0., 0e30);
  SaveWaveFunctionSpline(&GV, "Vintspline.dat");
#endif
#ifdef INTERACTION_Q1D_DIPOLES
  Message("Q1D Dipoles, D*(-2 x + sqrt(2*pi) (1 + x^2) Erfcx(x/sqrt(2))), with strength D = %lf\n", D);
#endif

#ifdef UNITS_SET_2M_TO_1
  Message("  Unit of energy: h^2/2ma^2\n");
#else
  Message("  Unit of energy: h^2/ma^2\n");
#endif

  Warning("Maximal allowed energy Eck = %" LE " \n", Eck);

  Message("\n  L=%" LE " \n", L);
  Lhalf = 0.5*L;
  Lhalf2 = Lhalf*Lhalf;
  L2 = L*L;
  Lm1 = 1. / L;
#ifdef BOX_EQUAL_SIDES
  Lx = Ly = Lz = L;
#endif

#ifdef CRYSTAL
#ifdef CRYSTAL_NONSYMMETRIC
  Message("  Solid phase simulation: non symmetrized trial wave function.\n");
#else
  Message("  Solid phase simulation: symmetrized trial wave function.\n");
#endif

  // Lx Lx n = Ndens
  ConstructGridCrystal(&Crystal, NCrystal, n * (DOUBLE) NCrystal / (DOUBLE) Ndens);

  if(generate_crystal_coordinates == OFF) LoadCrystalCoordinates();
#endif

#ifdef LATTICE_TRIANGULAR
#ifdef CRYSTAL
#ifdef TRIAL_2D
  Message("  Simulation is done in an asymmetric rectangular box\n");

#endif
#endif
#else
#ifdef CRYSTAL
  Message("  Simulation is done in a symmetric (square/cubic) box\n");
#endif
#endif

#ifdef ONE_BODY_TRIAL_LATTICE
  ConstructGridCrystalOneBody();
#endif

  L_half_x = 0.5*Lx;
  L_half_y = 0.5*Ly;
  L_half_z = 0.5*Lz;
  L_inv_x = 1./Lx;
  L_inv_y = 1./Ly;
  L_inv_z = 1./Lz;

#ifdef BC_3DPBC_TRUNCATED_OCTAHEDRON
  Lcutoff_pot = L*sqrt(3.)/4.;
#else
#ifdef BC_2DPBC_HEXAGON
  Lcutoff_pot = Lx/2.;
#else // standard case
  Lcutoff_pot = L/2.;
#endif
#endif
  Message("  Potential energy, corr functions, etc. are calculated up to the cut-off distance Rc= %lf\n", Lcutoff_pot);
  Lcutoff_pot2 = Lcutoff_pot*Lcutoff_pot;

  AllocateGrids(); // allocate memory for distribution arrays

#ifndef TRIAL_1D // 2D and 3D
  if(measure_Sk) LoadParametersSk();
#endif

#ifdef BC_ABSENT // in the trap load distribution for calculation of OBDM
  if(measure_OBDM) LoadRadialZDistribution();
#endif

  if(measure_FormFactor && MC == DIFFUSION) InitializeFormFactor();

#ifdef SCALABLE
  Lcell = L/(DOUBLE)Nscal;
  Nwf_store = N;
  Lwf = Lcell; // w.f. is constructed using Lcell
  N = N/Ncells;
#else
  Lwf = L;
#endif

#ifdef BC_3DPBC_TRUNCATED_OCTAHEDRON
  Lwf = L*sqrt(3.)/2.;
#endif

#ifdef BC_2DPBC_HEXAGON
  Lwf = Lx/2.;
#endif

  Lhalfwf = 0.5*Lwf;
  Lhalfwf2 = Lhalfwf*Lhalfwf;

#ifdef INTERPOLATE_SPLINE_ONE_BODY_Z_WF
  Message("  Spline interpolation will be used for one-body terms of w.f.\n");
  LoadTrialWaveFunction(&G1, file_wf1b, ON);
  SplineConstruct(&G1, 0., 0e30);

  SaveWaveFunction(&G1, "wf1b.dat", 0., 2.);
  SaveWaveFunctionSpline(&G1, "wf1bspline.dat");

  CheckTrialWaveFunction(&G1);
#endif

#ifdef VEXT_SPLINE_1D
  Message("  Spline interpolation will be used for Vext(z), loading from inVext.in file\n");
  LoadTrialWaveFunction(&GV, "inVext.in", OFF);
  SplineConstruct(&GV, 0., 0e30);

  SaveWaveFunctionSpline(&GV, "Vextspline.dat");
#endif

#ifdef ONE_BODY_TRIAL_SIN
#ifdef BC_3DPBC_CUBE
  Message("  one-body term: [sin(kL*x)^2 + sin(kL*y)^2 + sin(kL*z)^2]^alpha_R with kL = PI*Nlattice/L.\n");
#endif
#ifdef BC_2DPBC_SQUARE
  Message("  one-body term: [sin(kL*x)^2 + sin(kL*y)^2]^alpha_R with kL = PI*Nlattice/L.\n");
#endif
#ifdef BC_1DPBC_Z
  Message("  one-body term: sin(kL*z)^(2*alpha_R) with kL = PI*Nlattice/L.\n");
#endif
  if(kL<0) { // kL is used in w.f.
    kL = PI*Nlattice/L;
    Message("  Defining kL from Nlattice, kL = %lf\n", kL);
    if(fabs(kL)<1e-6) Error(" zero kL\n");
  }
#endif

#ifdef ONE_BODY_TRIAL_ZERO_BOUNDARY_CONDITION
  Message("  Zero boundary condition: one-body term: sin(kL*z) with kL = PI/L\n");
#endif

#ifdef ONE_BODY_TRIAL_COS_SERIES 
  Message("  one-body term (1D) 1 + alpha_latt*cos(2*k*z) + beta_latt*cos(4*k*z)+ gamma_latt*cos(6*k*z)\n");
  if(kL<0) { // kL is used in w.f.
    kL = PI*Nlattice/L;
    Message("  Defining kL from Nlattice, kL = %lf\n", kL);
  }
#endif 

#ifdef VEXT_COS2 // kL is used in Vext
  if(kL<0) {
    kL = PI*Nlattice/L;
    Message("  Defining kL from Nlattice, kL = %lf\n", kL);
  }
#endif

#ifdef ONE_BODY_TRIAL_ZBC
  Message("  one-body term ZBC: sin(kL*z) with kL = PI/L\n");
  kL = PI/L;
#endif

#ifdef ONE_BODY_TRIAL_ABC
  Message("  one-body term ABC: sin(kL*z)^alpha_latt with kL = PI/L\n");
  kL = PI/L;
#endif

#ifdef ONE_BODY_TRIAL_1D_WAVE_COS
  Message("  One-body term component: |cos(one_body_momentum_i pi z_j/L)|, where one_body_momentum_i = 0,1,2,... and j=1,...,N\n");
  one_body_momentum_k = PI/L*one_body_momentum_i;
  Message("  Momentum %"LE" corresponding to excitation of level %i\n", one_body_momentum_k, one_body_momentum_i);
#endif

#ifdef ONE_BODY_TRIAL_1D_WAVE_SIN
  Message("  One-body term component: |sin(one_body_momentum_i pi z_j/L)|, where one_body_momentum_i = 0,1,2,... and j=1,...,N\n");
  one_body_momentum_k = PI/L*one_body_momentum_i;
  Message("  Momentum %"LE" corresponding to excitation of level %i\n", one_body_momentum_k, one_body_momentum_i);
#endif

#ifdef ONE_BODY_TRIAL_1D_WAVE_COS_FINITE
  Message("  One-body term component: |1 + S cos(one_body_momentum_i pi z/L)|, where one_body_momentum_i = 0,1,2,...\n");
  one_body_momentum_k = 2.*PI/L*(DOUBLE) one_body_momentum_i;
  Message("  Momentum %"LE" corresponding to excitation of level %i\n", one_body_momentum_k, one_body_momentum_i);
#endif

#ifdef ONE_BODY_TRIAL_TRAP_DOWN_ANTISYMMETRIC_Z
  Message("  Antisymmetric state in z_dn\n");
#endif

#ifdef ONE_BODY_TRIAL_TRAP_DOWN_HERMITE2_Z
  Message("  Hermitite n=2 polynomial in z_dn, position of the node: z^2 = beta = %lf\n", beta);
#endif

  SaveOneBodyWaveFunction(); // afer one-body terms were defined

#ifdef MANY_BODY_TRIAL_TRAP_CENTER_OF_MASS_NODE_Z
  Message("  Hermite n=1 polynomial in center of mass z1 + ... + zN coordinate, position of the node z_{CM} = 0\n");
#endif

#ifdef MANY_BODY_TRIAL_TRAP_CENTER_OF_MASS_HERMITE2_Z
  Message("  Symmetrized Hermite n=2 polynomial, position of the node z_{CM}^2 = %lf\n", beta);
#endif

#ifdef VEXT_COS2
  //kL = PI*Nlattice/L;
  //if(Nlattice == 0) Error("Nlattice must be defined\n");
  Message("  External lattice with height %" LF " Erec and kL = %e is present\n", Srecoil, kL);

#ifndef BC_2DPBC_SQUARE // i.e. NOT in Q2D case
#ifndef UNITS_SET_2M_TO_1
  Warning("  Changing lattice height to recoil energies %lf [Er] -> %lf [Eo], Er = %lf\n", Srecoil, Srecoil*(kL*kL/2.),kL*kL/2.);
  Srecoil *= 0.5*kL*kL;
  Message("  set UNITS_SET_2M_TO_1 to get output in recoil energies\n");
#endif
#endif

#ifdef VEXT_SOLITON_FIXED_PHASE
  Message("  External effective potential: fixed phase soliton\n");
#endif

  Message("  Lattice size is %i, height %lf\n", Nlattice, Srecoil);

#ifdef BC_ABSENT // in a trapped case one can compare to harm. oscl. frequency
  Message(" Effective HO frequency omega = %lf\n", kL*sqrt(2.*Srecoil)); // S kL^2 x^2 = 1/2 omega^2 x^2
#endif

  SaveLattice();
#endif

#ifdef SPECKLES
  Message("\n  Speckle potential will be used.\n");
  SpecklesLoad();
#endif

//#ifdef IMPURITY // allocate memory and define
//  ConstructGridCrystal(&Crystal, NCrystal, n * (DOUBLE) NCrystal / (DOUBLE) Ndens);
//#endif

#ifdef VEXT_IMPURITY_3DSW
  ConstructVextSqWellDepth();
#endif

#ifdef INTERACTION_SOFT_SPHERE
  trial_Vo = ConstructSoftSphereHeight(Ro, a);
#endif

#ifdef INTERACTION_HARD_SPHERE_SQUARE_WELL
  trial_Vo = ConstructHardSphereSquareWellDepth(a, Ro, Apar);
#endif

#ifdef ONE_BODY_IMPURITY
  Message("  Impurities are present in the system.\n");
  ConstructGridImpurity(&Crystal, Rpar, Nimp);
  SaveImpurityWaveFunction();
#endif

#ifdef HARD_SPHERE_HYPERRADIUS
  Message("  Hard sphere hyperradius of size R3 = %lf is used (DMCor MoveAll)\n", R3);
#endif

  Message("\nConstructing particle-particle wavefunction ...\n  ");
#ifdef TRIAL_HS
  ConstructGridHS(&G);
#endif
#ifdef TRIAL_PSEUDOPOTENTIAL_PHONON
  ConstructGridPseudopotentialPhonon(&G);
#endif
#ifdef TRIAL_HS_SIMPLE
  ConstructGridHSSimple(&G);
#endif
#ifdef TRIAL_HS_SW_SIMPLE
  ConstructGridHSSWSimple(&G);
#endif
#ifdef TRIAL_PSEUDOPOTENTIAL_SIMPLE
  ConstructGridPseudopotentialSimple(&G);
#endif
#ifdef TRIAL_SS
  ConstructGridSS(&G, a, Apar, Rpar, Bpar, grid_trial, 0, Lwf);
#endif
#ifdef TRIAL_SS_LARGE
  ConstructGridSSlarge(&G);
#endif
#ifdef TRIAL_SS_EXP
  ConstructGridSSExp(&G);
#endif
#ifdef TRIAL_LIM
  ConstructGridLim(&G, 0, Lwf, grid_trial);
#endif
#ifdef TRIAL_POWER
  ConstructGridPower(&G, 0, Lwf, grid_trial);
#endif
#ifdef TRIAL_HS1D
  ConstructGridHS1D(&G, Rpar, Bpar, a, grid_trial, 0, Lwf);
#endif
#ifdef TRIAL_TONKS
  if(Rpar<1e-8) {
    Bpar = PI*PI/((L-2.*a)*(L-2.*a)); // k^2 for which f'(L/2) = 0
    Warning(" Bpar adjusted to the Tonks-Girardeau gas energy %" LE " \n\n", Bpar);
  }
  ConstructGridFree(&G, 1e100, Bpar, a,  grid_trial, 0, Lwf);
#endif
#ifdef TRIAL_TONKS_CUTOFF
  sqE = 0.5*PI/Rpar;
  sqE2 = sqE*sqE;
#endif
#ifdef TRIAL_LIEB
  ConstructGridLieb(&G, grid_trial);
#endif
#ifdef TRIAL_PHONON
  ConstructGridPhonons(&G);
#endif
#ifdef TRIAL_COS
  sqE = Apar*PI/L;
#endif
#ifdef TRIAL_SIN_EXCITED_BOSONS
  sqE = Apar*PI/L;
#endif
#ifdef TRIAL_COS_OFFSET
  sqE = Cpar*PI/L;
#endif
#ifdef TRIAL_PHONON_LUTTINGER
  ConstructGridPhononLuttinger(&G);
#endif
#ifdef TRIAL_PHONON_LUTTINGER_LATTICE
  ConstructGridPhononLuttinger(&G);
#endif
#ifdef TRIAL_PHONON_LUTTINGER_PIECEWISE
  ConstructGridPhononLuttinger(&G);
#endif
#ifdef TRIAL_PHONON_LUTTINGER_AUTO
  ConstructGridPhononLuttingerAuto(&G);
#endif
#ifdef TRIAL_PHONON_LUTTINGER_AUTO11_SPINFULL
  ConstructGridPhononLuttingerAuto11(&G);
#endif
#ifdef TRIAL_PHONON_LUTTINGER_AUTO12_SPINFULL
  ConstructGridPhononLuttingerAuto12(&G);
#endif
#ifdef TRIAL_LIEB_LUTTINGER_PHONON
  ConstructGridLiebLuttingerPhonon(&G, Rpar);
#endif
#ifdef TRIAL_LIEB_EXPONENTIAL_DECAY
  ConstructGridLiebExponentialDecay(&G);
#endif
#ifdef TRIAL_SS_PHONON_LUTINGER
  ConstructGridSoftSpherePhonon(&G);
#endif
#ifdef TRIAL_SS_PHONON_LUTINGER_AUTO
  ConstructGridSoftSpherePhononAuto(&G);
#endif
#ifdef TRIAL_SCHIFF_VERLET // with or without PHONONS
  ConstructGridSchiffVerlet(&G, Apar, Bpar, n);
#endif
#ifdef TRIAL_DIPOLES_PHONON
  ConstructGridDipolePhonon(&G);
#endif
#ifdef TRIAL_POWER_Ko
  ConstructGridPowerKo(&G);
#endif
#ifdef TRIAL_DIPOLE_Ko
  ConstructGridDipoleKo(&G);
#endif
#ifdef TRIAL_DIPOLE_Ko11
  ConstructGridDipoleKo11(&G);
#endif
#ifdef TRIAL_QUADRUPOLE_Ko
  ConstructGridQuadrupoleKo(&G);
#endif
#ifdef TRIAL_DIPOLE_Ko_TRAP
  Message(" TRIAL_DIPOLE_Ko_TRAP\n");
#endif
#ifdef TRIAL_DIPOLE_PHONON_FINITE_D
  ConstructGridDipolePhononFiniteD(&G);
#endif
#ifdef TRIAL_YUKAWA_2D
  ConstructGridYukawa2D(&G, Apar, Bpar, n);
#endif
#ifdef TRIAL_YUKAWA_3D_SYMM
  ConstructGridYukawa3D(&G, Apar, Bpar, Cpar);
#endif
#ifdef TRIAL_CSM_3D_PLASMON
  ConstructGridCSM3DPlasmon(&G, Rpar);
#endif
#ifdef TRIAL_YUKAWA_3D
  Message("  Yukawa 3D trial wavefunction (non symmetric)\n");
#endif
#ifdef TRIAL_SR3D
  Rpar = Lhalf;
  ConstructGridShortRange3D(&G, 1., Rpar);
#endif
#ifdef TRIAL_DECAY_PHONON
  ConstructGridDecayPhonon(&G, Apar);
#endif
#ifdef TRIAL_SUTHERLAND
  ConstructGridSutherland(&G);
#endif
#ifdef TRIAL_CALOGERO
  ConstructGridCalogero(&G);
#endif
#ifdef TRIAL_CALOGERO_2D
  ConstructGridCalogero2D(&G, Apar);
#endif
#ifdef TRIAL_COULOMB_1D_PHONON
  ConstructGridCoulomb1DPhonon(&G);
#endif
#ifdef TRIAL_DIPOLE_1D_PHONON
  ConstructGridDipole1DPhonon(&G);
#endif
#ifdef TRIAL_COULOMB_3D_PHONON
  ConstructGridCoulombPhonon(&G);
#endif
#ifdef TRIAL_DIPOLE_1D_TRAP
  Message("  Jastrow term: dipole 1D trap\n");
#endif
#ifdef TRIAL_COULOMB_1D_TRAP
  ConstructGridCoulomb1DTrap(&G);
#endif
#ifdef TRIAL_COULOMB_2D_TRAP
  ConstructGridCoulomb2DTrap(&G);
#endif
#ifdef TRIAL_EXP_PHONON
  ConstructGridExpPhonon(&G);
#endif
#ifdef TRIAL_HS2D_FINITE_ENERGY
  ConstructGridHS2DFiniteEnergy(&G);
#endif
#ifdef TRIAL_EXTERNAL
  LoadTrialWaveFunction(&G, file_wf, ON);
#endif
#ifdef TRIAL_TWO_BODY_SCATTERING_NUMERICAL
#ifdef JASTROW_DIMENSIONALITY_FROM_DPAR
  Message("  Dpar = %lf will be used as an effective dimensionality for constructing Jastrow term\n", Dpar);
#endif
  ConstructGridSpline2Body(&G, &InteractionEnergy);
#ifdef INTERPOLATE_LINEAR_JASTROW_WF
  LinearInterpolationConstruct(&G);
#endif
#endif
#ifdef TRIAL_TWO_BODY_SCATTERING_NUMERICAL_SPINFULL
  GridSpin = (struct Grid**) Calloc("GridSpin", Nspin, sizeof(struct Grid*)); // allocate [Nspin][Nspin] array of Jastrow terms
  for(i=0; i<Nspin; i++) GridSpin[i] = (struct Grid*) Calloc("GridSpin", Nspin, sizeof(struct Grid));

  //ConstructGridSpline2Body_ij(&GridSpin[0][2], 0, 2);
  for(i=0; i<Nspin; i++) for(w=0; w<Nspin; w++) ConstructGridSpline2Body_ij(&GridSpin[i][w], i, w);
  for(i=0; i<Nspin; i++) for(w=0; w<Nspin; w++) SplineConstruct(&GridSpin[i][w], 0., 0e30);
#endif

#ifdef TRIAL_KURBAKOV_NUMERIC
  ConstructGridKurbakovNumeric(&G);
#endif
#ifdef TRIAL_TWO_BODY_SCATTERING_PHONONS_NUMERICAL
  ConstructGridSpline2BodyPhonons(&G);
#endif
#ifdef TRIAL_HARD_RODS_3D
  ConstructGridHardRods3D(&G);
#endif
#ifdef TRIAL_SQ_WELL_ZERO_CONST
  ConstructGridSquareWellZeroConst();
#endif
#ifdef TRIAL_ZERO_RANGE_UNITARY
  ConstructGridZeroRangeUnitaryPhonons(&G);
#endif
#ifdef TRIAL_ZERO_RANGE_3D
  ConstructGridZeroRange3D(&G);
#endif
#ifdef TRIAL_SQUARE_WELL_FS
  ConstructGridSquareWellFreeState(&G);
#endif
#ifdef TRIAL_SQUARE_WELL_FS_PHONONS
  ConstructGridSquareWellFreeStatePhonons(&G);
#endif
#ifdef TRIAL_SQUARE_WELL_BS_PHONONS
  ConstructGridSquareWellBoundStatePhonons(&G);
#endif
#ifdef TRIAL_SQUARE_WELL_FS_TRAP
  ConstructGridSquareWellTrap();
#endif
#ifdef TRIAL_SQUARE_WELL_FS_TRAP_MCGUIRE1D
  ConstructGridSquareWellTrapMcGuire1D();
#endif
#ifdef TRIAL_SQUARE_WELL_BS_TRAP
  ConstructGridSquareWellBoundStateTrap(&G);
#endif
#ifdef ORBITAL_SQ_WELL_ZERO_ENERGY
  ConstructGridBCSSquareWellZeroEnergy();
#endif
#ifdef TRIAL_ZERO_RANGE_K0_TRAP
  Atrial = 2.*exp(-0.5772156649015329)/a;
#endif
#ifdef TRIAL_ZERO_RANGE_K0_EXP_TAIL_TRAP
  Atrial = 2.*exp(-0.5772156649015329)/a;
#endif
#ifdef THREE_BODY_ZERO_ENERGY_SOLUTION_1D_PHONONS
  ConstructGrid3BodyZeroEnergyPhonons();
#endif
#ifdef TRIAL_MCGUIRE_PHONON
  ConstructGridMcGuirePhonon();
#endif
#ifdef TRIAL_SPHERE_ZERO_RANGE_K0_SCREENING
  Message("  Jastrow: Apar,b: Ko(k*r) exp(-r*r/(b*(Apar+r))) where k = 2e^(-gamma)/a");
#endif
#ifdef TRIAL_SPHERE_ZERO_RANGE_LOG_SCREENING
  Message("  Jastrow: ln(a/r) * exp(-r/b)\n");
  if(a<R) Error(" scattering length must be larger than the radius, a>R\n");
#endif
#ifdef TRIAL_SOFT_DISKS_2D_ZERO_ENERGY_PHONONS
  ConstructGridSoftDisks2DZeroEnergyPhonons(&G);
#endif
#ifdef TRIAL_SS_1D
  ConstructSS_Simple(&G);
#endif

  Message("done\n");

#ifdef SCALABLE
  N = Nwf_store;
#endif

#ifdef SCALE_VMC_AMPLITUDE_WITH_ALPHA
  Message(" VMC displacement amplitude will be scaled by 1/sqrt(alpha_Rx, alpha_Ry, alpha_Rz) = (%lf, %lf, %lf)\n", 1./sqrt(alpha_Rx), 1./sqrt(alpha_Ry), 1./sqrt(alpha_Rz));
#endif

#ifdef LOAD_WF_FROM_FILE
  LoadTrialWaveFunction(&G, file_wf, ON);
#endif

#ifdef SAMPLE_WF
  Message("  Distributing w.f. on grids for a better performance...\n");
  SampleTrialWaveFunction(&G, 0., Lhalfwf, grid_trial, &InterpolateExactU, &InterpolateExactFp, &InterpolateExactE);
#endif

#ifdef INTERPOLATE_SPLINE_JASTROW_WF
  Message("  Spline interpolation will be used for Jastrow terms\n");
#ifdef JASTROW_LEFT_BOUNDARY_PSEUDOPOTENTIAL
  SplineConstruct(&G, -G.f[0]/a, 0e30); //  f'(0) = -f(0)/a
#else 
  SplineConstruct(&G, 0., 0e30); // f'(0) = 0
#endif
  SaveWaveFunctionSpline(&G, "wfspline.dat");
#endif

#ifdef SAMPLE_WF
  CheckTrialWaveFunctionInterpolation(&G, &InterpolateExactU, &InterpolateExactFp, &InterpolateExactE);
#endif

#ifdef SYMMETRIZE_TRIAL_WF
  Message("  An explicit symmetrization of w.f. will be used according to the rule u(r) -> u(r)+u(L-r)\n");
#endif

#ifdef POWER_TRIAL_WF
  Message("  An explicit change of the Jastrow terms will be done according to the rule u(r) -> alpha_Jastrow u(r) with alpha=%e\n", alpha_Jastrow);
#endif

  //SaveWaveFunction(&G, "wf.dat", 1e-6, Lhalfwf);
  SaveWaveFunction(&G, "wf.dat", 0., Lhalfwf);

  CheckTrialWaveFunction(&G);

  // be sure that the interaction energy vanishes at r=L/2
#ifdef INTERACTION_ABSENT
  Message("  Interaction potential (not included in Eloc part): ABSENT\n");
#endif
#ifdef INTERACTION_DIPOLE
  Message("  Interaction potential (not included in Eloc part): DIPOLE\n");
#endif
#ifdef INTERACTION_DIPOLE_TRUNCATED
  Message("  Interaction potential (not included in Eloc part): 1/r^3, r>a and const otherwise\n");
#endif
#ifdef INTERACTION_LENNARD_JONES
  Message("  Interaction potential (not included in Eloc part): LENNARD JONES\n");
#endif
#ifdef INTERACTION_GAUSSIAN
    Message("  Interaction potential (not included in Eloc part): Gaussian\n");
#endif
#ifdef INTERACTION_RYDBERG
  Message("  Interaction potential (not included in Eloc part): 1/r^6\n");
#endif
#ifdef INTERACTION_RYDBERG_TRUNCATED
  Message("  Interaction potential (not included in Eloc part): 1/r^6, r>a and const otherwise\n");
#endif

  Message("  Summation of interaction potential over all images: ");
#ifdef INTERACTION_SUM_OVER_IMAGES
  Message("direct summation\n");
#endif
#ifdef INTERACTION_EWALD
  Message("Ewald summation\n");
#endif
#ifdef INTERACTION_SUM_OVER_IMAGES_SMOOTH_CUTOFF
  Message("smooth cutoff summation\n");
#endif
#ifdef INTERACTION_WITHOUT_IMAGES
  Message("disabled\n");
#endif

  Etail = 0.5*n*IntegrateSimpson(Lhalf, Lhalf*10, &InterpolateEG, 100000L);
  Message("  Tail energy %" LE " \n", energy_unit*Etail);
  Message("  Mean-field energy (g2=1) %" LE " \n", energy_unit*0.5*n*IntegrateSimpson(1e-8*Lhalf, Lhalf*10, &InterpolateEG, 100000L));

#ifdef TRIAL_1D
  Message("  Tonks energy %" LE " \n", energy_unit*PI*PI*n*n/6.);
#endif

#ifdef THREE_BODY_TERMS
  Message("  Three-body terms in the wave function are included in the form of:\n");
#ifdef THREE_BODY_HYPERRADIUS _QUADRATIC
  Message("    f_3(R): 0, R<R3; 1 - (R - Apar)^2 / (R3 - Apar)^2 , R3<R < Rpar; 1, R>Apar\n");
#endif
#ifdef THREE_BODY_ZERO_ENERGY_SOLUTION
  Message("    f_3(R): 0, R<R3; 1 - (R3/R)^4\n");
#endif
#ifdef THREE_BODY_ZERO_ENERGY_SOLUTION_SYMMETRIZED 
  Message("    f_3(R): f_3(R):  1 - (R3/R)^4 - (R3/(L-R))^4\n");
#endif
#ifdef THREE_BODY_ZERO_ENERGY_SOLUTION_1D_SYMMETRIZED
  Message("    f_3(R): 0, R<R3; log(R/R3)+log((L-R)/R3)\n");
#endif
#ifdef THREE_BODY_ZERO_ENERGY_SOLUTION_1D_PHONONS
  Message("    (1D,2D,3D) f_3(R): 0, R<R3; log(R/R3), r<Rpar; B exp(-C/r + D/r^2)\n");
#endif
  SaveThreeBodyWaveFunction();
#endif

#ifdef JASTROW_RIGHT_BOUNDARY_BOUND_STATE
  Message("    Jastrow term: large distance is taken from bound state asymptotic\n");
#endif 

#ifdef TRAP_POTENTIAL
  CaseX(Message("  Harmonic potential X: %lf x^2\n", 0.5*omega_x2));
  CaseY(Message("  Harmonic potential Y: %lf y^2\n", 0.5*omega_y2));
  CaseZ(Message("  Harmonic potential Z: %lf z^2\n", 0.5*omega_z2));
#endif

  if(optimization == 2) {
    Optimization();
    return 0;
  }

#ifdef EXTERNAL_POTENTIAL
  SaveVext();
#endif

#ifdef INTERACTION_EWALD // Initialize Ewald summation before calculating energy for particle coordinates
  EwaldInit();
#endif

#ifdef INTERACTION_SUM_OVER_IMAGES_SMOOTH_CUTOFF // Initialize and check smooth summation
  SmoothCutoffInit();
#endif
  if(generate_new_coordinates) {
    AllocateWalkers();
    Nwalkers = Npop;
    GenerateCoordinates();
  }
  else {
    LoadParticleCoordinates();
  }

  if(measure_Hessian_matrix) MeasureHessianMatrix();

  if(MC == DIFFUSION && measure_pure_coordinates) CopyPureCoordinates();

#ifdef MEMORY_CONTIGUOUS
  Message("  \nMemory allocation: contiguous arrays (i.e. optimized)\n");
#else
  Message("  \nMemory allocation: scattered arrays (i.e. non optimized)\n");
#endif
  CheckMemcpy();

#ifdef INTERACTION_EWALD // Initialize and check Ewald summation
  EwaldCheck();
#endif

#ifdef CRYSTAL
  if(measure_energy_barrier) SaveDefectBarrier();
#endif

  if(MC == DIFFUSION) {
    Npop_max = Npop + Npop / 10;
    Npop_min = Npop - Npop / 10;

    reduce  = 0.5*(DOUBLE) (Npop_min + Npop_max) / (DOUBLE) (Npop_max);
    amplify = 0.5*(DOUBLE) (Npop_min + Npop_max) / (DOUBLE) (Npop_min);
    Ndead_walkers = 0;
    for(i=0; i<Nwalkers; i++) {
      W[i].status = ALIVE;
      walkers_storage[i] = i;
    }
    for(i=Nwalkers; i<NwalkersMax; i++) {
      W[i].status = DEAD;
      dead_walkers_storage[i-Nwalkers] = i;
      Ndead_walkers++;
    }
    if(measure_SD) {
      for(w=0; w<Nwalkers; w++) {
        for(i=0; i<N; i++) {
          W[w].rreal[i][0] = W[w].x[i];
          W[w].rreal[i][1] = W[w].y[i];
          W[w].rreal[i][2] = W[w].z[i];
        }
      }
    }
  }

  if((MC == DIFFUSION || MC == VARIATIONAL) && check_kinetic_energy) CheckEnergyAndDerivatives();

#ifdef MPI
  ParallelInitialize2(); // allocate arrays
#endif

  Eo = 0;
  for(w=0; w<Nwalkers; w++) {
#ifdef SCALABLE
    UpdateScalableTable(&W[w]);
#endif
#ifdef CHAINS_PINNED_TO_CENTER_OF_MASS
    CalculateAllChainsCenterOfMassWalker(&W[w]);
#endif 
    W[w].Eold = W[w].E = WalkerEnergy0(&W[w]);
    W[w].U = U(W[w]);
    Eo += W[w].E;
  }
  Eo /= (DOUBLE) Nwalkers;

  E = Energy(&Epot, &Ekin, &EFF, &Edamping, &Eint, &Eext);
  Message("\nEnergy %" LG "  = %i * %" LG " , EFF = %" LG "\n", E*N, N, E, EFF);

  Message("Starting heating...\n");
  time_communicating = clock() - time_start;
  dt_Kurbakov=dt;

  for(block=-blck_heating; block<blck; block++) {
    time_begin = clock();
    accepted = accepted_one = accepted_all = 0; // control rate initialization
    rejected = rejected_one = rejected_all = 0;
    overlaped = 0;
    tried = 0;
    Eblck = 0.;
    measured_blck = 0;

    if(block == 0) Message("Starting measurements...\n");
    EmptyDistributionArrays();

    if(verbosity>1) Message("Block No %i\n", block+1);

    for(iter=1; iter<=Niter; iter++) {
      if(MC == DIFFUSION && generate_new_coordinates && block<0)dt=dt_Kurbakov*exp(-80.*exp(20.*log((.5-(DOUBLE)block*(DOUBLE)Niter-(DOUBLE)iter)/(DOUBLE)Niter/(DOUBLE)blck_heating)));
#ifdef MPI
      if(iter == 1 && block==-blck_heating) {
        DistributeAllWalkers();
        time_done_in_parallel -= MPI_Wtime();  // initialize counter of the parallel job
      }
#endif
      iteration_global++;

#pragma omp parallel for
      for(w=0; w<Nwalkers+Nwalkers_pure; w++) {
        if(MC == VARIATIONAL) { // Variational Monte Carlo
          if(SmartMC == VMC_MOVE_ALL) {
            VMCMoveAll(w);
          }
          else if(SmartMC == VMC_MOVE_ONE) {
            VMCMoveOneByOne(w);
          }
          else if(SmartMC == VMC_MOVE_DRIFT_ALL) {
            VMCMoveDrift(w);
          }
          else if(SmartMC == VMC_MOVE_DRIFT_ONE) {
            VMCMoveDriftOneByOne(w);
          }
          else { // VMC_SUM_OF_TWO_WAVEFUNCTIONS 
            //VMCMoveSumOfTwoWavefunctions(w);
            VMCMoveAll(w);
          }
        }
        else if(MC == DIFFUSION) { // Diffusion Monte Carlo
          if(SmartMC == DMC_QUADRATIC) {
            DMCQuadraticMove(w);
          }
          else if(SmartMC == DMC_LINEAR_DRIFT) {
            DMCLinearMoveDrift(w);
          }
          else if(SmartMC == DMC_LINEAR_GAUSSIAN) {
            DMCLinearMoveGaussian(w);
          }
          else if(SmartMC == DMC_LINEAR_METROPOLIS) {
            DMCLinearMetropolisMove(w);
          }
          else if(SmartMC == DMC_QUADRATIC_METROPOLIS) {
            DMCQuadraticMetropolisMove(w);
          }
          else if(SmartMC == DMC_QUARTIC) {
            DMCQuarticMove(w);
          }
          else if(SmartMC == DMC_PSEUDOPOTENTIAL) {
            if(w==0) DMCPseudopotentialMove(w);
          }
          else if(SmartMC == DMC_MOVE_ONE) {
            DMCLinearMoveOne(w);
          }
        }
        else if(MC == CLASSICAL) { // classical Monte Carlo
          ClassicalMoveOneByOne(w);
        }
        else if(MC == PIMC) {
          PIMCMoveOneByOne(w);
          //PIMCPseudopotentialMoveOneByOne(w);
          //PathIntegralMoveAll();
          //PathIntegralPseudopotentialStaging();
        }

        if(MC == CLASSICAL && Tannealing != 1.) {
          T /= Tannealing;
          //dt_vmc /= Tannealing;
        }
      } // OMP for pragma
//#pragma omp single

      if(MC == DIFFUSION) { // Branching
#ifdef MPI
        //BranchingMPI();

        time_done_in_parallel += MPI_Wtime();  // update counter
        CollectAllWalkers();
        if(myid  == MASTER) Branching();
        time_done_in_parallel -= MPI_Wtime();  // update counter

        E = Eo / (DOUBLE) N;
        if(measure_energy && block >= 0 && iter % Nmeasure == 0) {
          SaveEnergyDMC(E);
        }
        if(iter%Niter==0) {
          SaveCoordinates(file_particles);
          Message("saving coordinates\n");
        }
        Evar += E;
        Evar2 += E*E;

        DistributeAllWalkers();
#else
        Branching();
#endif
      }

      if(MC == DIFFUSION && block>= 0 && measure_SD && iter % SD.spacing == 0) MeasureSuperfluidDensity();
      //if(Niter_store && block>=0 && iter % Niter_store == 0) SaveCoordinatesOptimization();

      if(block >= 0 && iter % Nmeasure == 0) {
        time_simulating += clock() - time_begin;
        Measure(block, iter);
        time_begin = clock();
        if(video) Picture();
      }
    }

    time_simulating += clock() - time_begin;

    Message("%i ", block+1); // print acceptance rates
    time(&time_era_end);
    MessageTimeElapsed((int)(time_era_end-time_era_begin));

    if(acceptance_rate != 0. && (MC != DIFFUSION)) {
      Message("dt= %.2g ", dt_vmc);
      if(accepted) // adjust the timestep
        dt_vmc *= (DOUBLE) accepted/((DOUBLE) (accepted+rejected)*acceptance_rate);
      else
        dt_vmc *= 1e-2;
    }
    if(accepted+rejected) Message("(%.2" LF "%%) ", (DOUBLE) accepted/(DOUBLE)(accepted+rejected)*100.);
    if(accepted_one+rejected_one) Message("(one %.2" LF "%%) ", (DOUBLE) accepted_one/(DOUBLE)(accepted_one+rejected_one)*100.);
    if(accepted_all+rejected_all) Message("(all %.2" LF "%%) ", (DOUBLE) accepted_all/(DOUBLE)(accepted_all+rejected_all)*100.);
    if(overlaped) Message("[%.2" LF "%%] ", (DOUBLE) overlaped/(DOUBLE)(tried)*100.);
    if(MC == CLASSICAL && Tannealing != 1.) Message("T=%.1" LE " ", T);

    Message("<E>=%3g E=%3g ", energy_unit*Evar/(DOUBLE) Nmeasured, energy_unit*E);
    Message("EFF=%3g Ekin=%3g Epot=%" LG "  ", energy_unit*EFF, energy_unit*Ekin, energy_unit*Epot);
    if(MC == DIFFUSION) Message("Nw=%i ", Nwalkers);

#ifdef BC_ABSENT
    Message(" <%.2" LF ",%.2" LF ">", RD.r2recent, RDz.r2recent);
#endif

    Message("\n");
  }

#ifdef MPI
  time_done_in_parallel += MPI_Wtime();  // update counter
  CollectAllWalkers();
  SaveCoordinates(file_particles);
#endif

  if(video) CloseGraph();

  if(Nmeasured == 0) 
    Warning("No measurements done\n");
  else {
    Evar /= (DOUBLE) Nmeasured;
    Evar2 /= (DOUBLE) Nmeasured;
    sigma = Sqrt(Evar2-Evar*Evar);
    error = sigma/(Sqrt((DOUBLE)(Nmeasured-1)));
    Message("Evar = %.15" LE " +/- %" LE "\n", energy_unit*Evar, energy_unit*error);

    EFFvar /= (DOUBLE) Nmeasured;
    EFFvar2 /= (DOUBLE) Nmeasured;
    sigma = Sqrt(EFFvar2-EFFvar*EFFvar);
    error = sigma/(Sqrt((DOUBLE)(Nmeasured-1)));
    Message("EFF  = %" LE " +/- %" LE "\n", energy_unit*EFFvar, energy_unit*error);

    Epotvar /= (DOUBLE) Nmeasured;
    Epotvar2 /= (DOUBLE) Nmeasured;
    sigma = Sqrt(Epotvar2-Epotvar*Epotvar);
    error = sigma/(Sqrt((DOUBLE)(Nmeasured-1)));
    Message("Epot = %" LE " +/- %" LE "\n", energy_unit*Epotvar, energy_unit*error);

    Ekinvar /= (DOUBLE) Nmeasured;
    Ekinvar2 /= (DOUBLE) Nmeasured;
    sigma = Sqrt(Ekinvar2-Ekinvar*Ekinvar);
    error = sigma/(Sqrt((DOUBLE)(Nmeasured-1)));
    Message("Ekin = %" LE " +/- %" LE "\n", energy_unit*Ekinvar, energy_unit*error);

    if(boundary == NO_BOUNDARY_CONDITIONS) Message("Eint = %" LF "\n", 2.*(Epotvar-Ekinvar));
    Message("  Total energy %.15" LE "\n", energy_unit*(Evar + Etail));

    if(R2times_measured) Message("<r2>^0.5 = %" LF "\n", Sqrt(xR2/(DOUBLE)R2times_measured));
    if(R2times_measured) Message("<z2>^0.5 = %" LF "\n", Sqrt(zR2/(DOUBLE)R2times_measured));
 }

  time_end = clock();
  time(&time_era_end);
  Message("\ntime: Ttotal   %.2" LF " sec = %.2" LF " hours = %.2" LF" days\n", (DOUBLE)(time_era_end-time_era_begin), (DOUBLE)(time_era_end-time_era_begin)/3600., (DOUBLE)(time_era_end-time_era_begin)/3600/24.);
  Message("time: Tprecise %.2" LF " sec = %.2" LF " hours = %.2" LF" days\n", (DOUBLE)(time_end-time_start)/(DOUBLE)CLOCKS_PER_SEC, (DOUBLE)(time_end-time_start)/(DOUBLE)CLOCKS_PER_SEC/3600., (DOUBLE)(time_end-time_start)/(DOUBLE)CLOCKS_PER_SEC/3600./24.);
  Message("time used for doing simulation    %.2" LF " seconds\n", (DOUBLE) time_simulating/(DOUBLE)CLOCKS_PER_SEC);
  Message("time used for doing measurements  %.2" LF " seconds\n", (DOUBLE) time_measuring/(DOUBLE)CLOCKS_PER_SEC);
  Message("time used for doing communication %.2" LF " seconds\n", (DOUBLE) time_communicating/(DOUBLE)CLOCKS_PER_SEC);

#ifdef MPI

  ParallelStop();
#endif

  Message("done\n");
  return 0;
}

/*********************************** Measure *********************************/
void Measure(int block, int iter) {
  int w;
  clock_t time_begin;

  time_begin = clock();
  Nmeasured++;

  if(measure_OBDM_MATRIX) MeasureOBDMMatrix();
  if(measure_PairDistr) MeasurePairDistribution();
  if(measure_g3) MeasureHyperRadius();
  if(measure_OP) MeasureOrderParameter();
  if(measure_Lind) MeasureLindemannRatio();
  if(measure_wavefunction_projection) MeasureWaveFunctionProjection();
  if(measure_tVMC) MeasuretVMC(iter);

  if(measure_energy && MC == VARIATIONAL) { //  || SmartMC == 2
    Eo = E = Energy(&Epot, &Ekin, &EFF, &Edamping, &Eint, &Eext);
    E = Epot + Ekin;
    // reweightig procedure
    //if(reweighting) MeasureEnergyReweighting();
    SaveEnergyVMC(E, EFF);
  }

  //if(measure_OBDM) MeasureOBDM();
  if(measure_TBDM) MeasureTBDM();
  if(measure_effective_potential) MeasureEffectivePotential();

#pragma omp parallel for
  for(w=0; w<Nwalkers; w++) {
    if(measure_RadDistr) MeasureRadialDistributionWalker(&W[w]);
    if(measure_Sk) MeasureStaticStructureFactor(w);
    if(measure_FormFactor) MeasureFormFactor(w);
    if(measure_OBDM) MeasureOBDM(w);
  }

  if(MC == VARIATIONAL) {
    if(verbosity > 1) {
      if(boundary == THREE_BOUNDARY_CONDITIONS) {
        Message("%i %li E=%3e EFF=%3e Ekin=%" LE "  Nw=%i\n", block+1, iter, E, EFF, Ekin, Nwalkers);
      }
      else if(boundary == ONE_BOUNDARY_CONDITION) {
        Message("%i %li E=%3lf EFF=%3lf Ep=%3lf Ekin=%" LF " <%.2" LF "> Nw=%i\n", 
          block+1, iter, E, EFF, Epot, Ekin, RD.r2 / (DOUBLE) (RD.times_measured * N), Nwalkers);
      }
      else {
#ifdef LATTICE
        Message("%i %li E=%3lf EFF=%3lf Ep=%3lf Ekin=%" LF " Elat=%" LF " Nw=%i\n", block+1, iter, E, EFF, Epot, Ekin, Elat, Nwalkers);
#else
        Message("%i %li E=%3lf EFF=%3lf Ep=%3lf Ekin=%" LF " <%.2" LF ",%.2" LF "> Nw=%i\n", 
          block+1, iter, E, EFF, Epot, Ekin, RD.r2 / (DOUBLE) (RD.times_measured * N),
          RDz.r2 / (DOUBLE) (RDz.times_measured * N), Nwalkers);
#endif
      }
    }
    Evar += E;
    Evar2 += E*E;
    EFFvar += EFF;
    EFFvar2 += EFF*EFF;
    Epotvar += Epot;
    Epotvar2 += Epot*Epot;
    Ekinvar += Ekin;
    Ekinvar2 += Ekin*Ekin;
    if(measure_R2 && iter % Nmeasure == 0) SaveMeanR2(); // uses same counter as Radial Distribution
#ifdef SPINFULL
    SaveSpinPolarization();
#endif
  }
  else if(MC == DIFFUSION) { // Diffusion Monte Carlo method
#ifndef MPI
    E = Eo / (DOUBLE) N;

    if(measure_energy) {
      MeasurePurePotentialEnergy();
      SaveEnergyDMC(E);
    }
    Evar += E;
    Evar2 += E*E;
    if(verbosity>1) Message("%i %li E = %" LE "  walkers: %i\n", block+1, iter, E, Nwalkers);
#endif

#ifdef SPINFULL
    SaveSpinPolarization();
#endif

    if(measure_R2 && iter % Nmeasure == 0) SaveMeanR2DMC();

     if(measure_pure_coordinates && iter % (grid_pure_block*Nmeasure) == 0) SavePureCoordinates();
  }
  else if(MC == CLASSICAL) { // Classical Monte Carlo method

    if(measure_RadDistr) MeasureRadialDistributionWalker(&W[0]);
    if(measure_Sk) MeasureStaticStructureFactor(0);
#pragma omp parallel for
    for(w=1; w<Nwalkers; w++) {
      if(measure_RadDistr) MeasureRadialDistributionWalker(&W[w]);
      if(measure_Sk) MeasureStaticStructureFactor(w);
    }

    if(measure_R2 && measure_RadDistr) SaveMeanR2();

    Eo = E = Energy(&Epot, &Ekin, &EFF, &Edamping, &Eint, &Eext);

    for(w=0; w<Nwalkers; w++) {
      if(measure_energy) {
        WalkerEnergy(&W[w], &Epot, &Ekin, &EFF, &Edamping, &Eint, &Eext);
        Epot /= (DOUBLE) N;
        SaveEnergyCls(Epot);
      }
    }
  }
  else if(MC == PIMC) {
    if(measure_RadDistr) MeasureRadialDistribution();
    if(measure_R2 && measure_RadDistr) SaveMeanR2();

    Eo = E = Energy(&Epot, &Ekin, &EFF, &Edamping, &Eint, &Eext);

    if(measure_energy && iter % Nmeasure == 0) {
      SaveEnergyCls(E);
    }
  }

  for(w=0; w<Nwalkers; w++) W[w].times_measured++; // should be executed last

  if(iter % Niter == 0) SaveBlockData(block);

  time_measuring += clock() - time_begin;
}
