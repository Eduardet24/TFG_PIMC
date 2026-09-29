/*rw.c*/

#include <stdio.h>
#include <string.h>
#include "rw.h"
#include "main.h"
#include "memory.h"
#include "utils.h"
#include "vmc.h"
#include "quantities.h"
#include "trial.h"
#include "spline.h"
#include "speckles.h"
#include "crystal.h"
#include "compatab.h"
#include "dmc.h"
#include MATHINCLUDE
#ifdef MPI
#  include "parallel.h"
#endif

#define CHECK(text, type, address) {if(strcmp(ch,text)==0) {    fscanf(in, type, address); Message("\n  %s" type, text, *address);    }  }
#define CHECKs(text, type, address) {if(strcmp(ch,text)==0) {    fscanf(in, type, address);      Message("\n  %s" type, text, address);    }  }

/******************************** Load Config ******************************/
int LoadParameters(void) {
  FILE *in;
  char ch[512]; //cambio de 100 a 512 por core dumped al compilar con --> gcc -g -fsanitize=address -fopenmp *.c -o MC -lm

  Message("Loading parameters ... ");
  in = fopen(INPATH "MC.cfg", "r");

  if(in == NULL) {
    perror("\nError: can't load MC.cfg file,\nexiting ...");
    Error("\nError: can't load " INPATH "MC.cfg file,\nexiting ...");
  }
  RD.width = RDz.width = 1000;
  RD.max = RDz.max = 5;
  OBDM.Nksize = 100;
  OBDM.max = -1.;
  SD.spacing = 1;
  NCrystal = -1;
  Ndens = -1;
  kL = -1;
  SKT_dk = -1.;
  L = -1.;

  RD_MATRIX.maxx = RD_MATRIX.maxy = 0.;
  Dpar = (DOUBLE) DIMENSION;

  while(fscanf(in, "%s\n", ch) != EOF) {
    CHECK("MC=","%i", &MC);
    CHECK("SmartMC=","%i", &SmartMC);
    CHECK("N=","%i", &N);
    CHECK("Ndens=","%i", &Ndens);
    CHECK("NCrystal=", "%i", &NCrystal);
    CHECK("Nimp=","%i", &Nimp);
    CHECK("Nspin=","%i", &Nspin);
    CHECK("n=", "%" LF "", &n);
    CHECK("a=", "%" LF "", &a);
    CHECK("b=", "%" LF "", &b);
    CHECK("Niter=", "%li", &Niter);
    CHECK("blck_heating=", "%i", &blck_heating);
    CHECK("blck=", "%i", &blck);
    CHECK("dt=", "%" LF "", &dt);
    CHECK("dt_one=", "%" LF "", &dt_one);
    CHECK("dt_all=", "%" LF "", &dt_all);
    CHECK("dt_vmc=", "%" LF "", &dt_vmc);
    CHECK("Rpar=", "%" LE, &Rpar);
    CHECK("Apar=", "%" LF "", &Apar);
    CHECK("Bpar=", "%" LF "", &Bpar);
    CHECK("Cpar=", "%" LF "", &Cpar);
    CHECK("Dpar=", "%" LF "", &Dpar);
    CHECK("R0par=", "%" LF "", &R0par);
    CHECK("Epar=", "%" LF "", &Epar);
    CHECK("Kpar=", "%" LE, &Kpar);
    CHECK("Mpar=", "%" LE, &Mpar);
    CHECK("Rpar12=", "%" LE, &Rpar12);
    CHECK("Apar12=", "%" LF "", &Apar12);
    CHECK("Bpar12=", "%" LF "", &Bpar12);
    CHECK("Cpar12=", "%" LF "", &Cpar12);
    CHECKs("file_particles=", "%s", file_particles);
    CHECKs("file_wf=", "%s", file_wf);
    CHECKs("file_wf1b=", "%s", file_wf1b);
    CHECKs("file_energy=", "%s", file_energy);
    CHECKs("file_OBDM=", "%s", file_OBDM);
    CHECKs("file_OBDM_MATRIX=", "%s", file_OBDM_MATRIX);
    CHECKs("file_PD=", "%s", file_PD);
    CHECKs("file_PD_pure=", "%s", file_PD_pure);
    CHECKs("file_RD=", "%s", file_RD);
    CHECKs("file_PDz=", "%s", file_PDz);
    CHECKs("file_RDz=", "%s", file_RDz);
    CHECKs("file_R2=", "%s", file_R2);
    CHECKs("file_z2=", "%s", file_z2);
    CHECKs("file_R2_pure=", "%s", file_R2_pure);
    CHECKs("file_z2_pure=", "%s", file_z2_pure);
    CHECKs("file_Sk=", "%s", file_Sk);
    CHECKs("file_Sk_pure=", "%s", file_Sk_pure);
    CHECKs("file_SD=", "%s", file_SD);
    CHECK("verbosity=", "%i", &verbosity);
    CHECK("Nmeasure=", "%li", &Nmeasure);
    CHECK("grid_trial=", "%i", &grid_trial);
    CHECK("gridOBDM=", "%i", &gridOBDM);
    CHECK("gridNk=", "%i", &OBDM.Nksize);
    CHECK("gridOBDM_MATRIX=", "%i", &gridOBDM_MATRIX);
    CHECK("gridRD=", "%i", &gridRD);
    CHECK("gridRDx=", "%i", &gridRDx);
    CHECK("gridRDy=", "%i", &gridRDy);
    CHECK("gridPD=", "%i", &gridPD);
    CHECK("gridg3=", "%i", &gridg3);
    CHECK("measure_g3=", "%i", &measure_g3);
    CHECK("McMillan_points=", "%i", &McMillan_points);
    CHECK("McMillanTBDM_points=", "%i", &McMillanTBDM_points);
    CHECK("NwalkersMax=", "%i", &NwalkersMax);
    CHECK("Npop=", "%i", &Npop);
    CHECK("Npop_virtual=", "%i", &Npop_virtual);
    CHECK("generate_new_coordinates=", "%i", &generate_new_coordinates);
    CHECK("generate_new_coordinates_from_lattice=", "%i", &generate_new_coordinates_from_lattice);
    CHECK("measure_energy=", "%i", &measure_energy);
    CHECK("measure_OBDM=", "%i", &measure_OBDM);
    CHECK("measure_TBDM=", "%i", &measure_TBDM);
    CHECK("measure_OBDM_MATRIX=", "%i", &measure_OBDM_MATRIX);
    CHECK("measure_RadDistr=", "%i", &measure_RadDistr);
    CHECK("measure_RadDistrMATRIX=", "%i", &measure_RadDistrMATRIX);
    CHECK("measure_PairDistr=", "%i", &measure_PairDistr);
    CHECK("measure_PairDistr_R2=", "%i", &measure_PairDistr_R2);
    CHECK("acceptance_rate=", "%" LF "", &acceptance_rate);
    CHECK("video=", "%i", &video);
    CHECK("lambda=", "%" LF "", &lambda);
    CHECK("alpha_x=", "%" LF "", &alpha_x);
    CHECK("alpha_y=", "%" LF "", &alpha_y);
    CHECK("alpha_z=", "%" LF "", &alpha_z);
    CHECK("alpha_R=", "%" LF "", &alpha_R);
    CHECK("alpha=", "%" LF "", &alpha_R);
    CHECK("alpha_Rx=", "%" LF "", &alpha_Rx);
    CHECK("alpha_Ry=", "%" LF "", &alpha_Ry);
    CHECK("alpha_Rz=", "%" LF "", &alpha_Rz);
    CHECK("measure_R2=", "%i", &measure_R2);
    CHECK("RDwidth=", "%" LF "", &RD.width);
    CHECK("RDzwidth=", "%" LF "", &RDz.width);
    CHECK("RDmax=", "%" LF "", &RD.max);
    CHECK("RDzmax=", "%" LF "", &RDz.max);
    CHECK("PDmax=", "%" LF "", &PD.max);
    CHECK("Skmax=", "%" LF "", &Sk.L);
    CHECK("Nkmin=", "%" LF "", &OBDM.kmin);
    CHECK("Nkmax=", "%" LF "", &OBDM.kmax);
    CHECK("OBDMmax=", "%" LF "", &OBDM.max);
    CHECK("OBDM_MATRIXmax=", "%" LF "", &OBDM_MATRIX.max);
    CHECK("gridPD_pure_block=", "%i", &grid_pure_block);
    CHECK("grid_pure_block=", "%i", &grid_pure_block);
    CHECK("measure_Sk=", "%i", &measure_Sk);
    CHECK("gridSk=", "%i", &gridSk);
    CHECK("gridSD=", "%i", &SD.size);
    CHECK("SDspacing=", "%i", &SD.spacing);
    CHECK("measure_SD=", "%i", &measure_SD);
    CHECK("measure_OP=", "%i", &measure_OP);
    CHECK("D=", "%" LF "", &D);
    CHECK("Q=", "%" LF "", &D);
    CHECK("measure_Lind=", "%i", &measure_Lind);
    CHECK("Nscal=", "%i", &Nscal);
    CHECK("Rzigzag=", "%" LF "", &Rzigzag);
    CHECK("measure_Nk_pure=", "%i", &measure_Nk_pure);
    CHECK("T=", "%" LF "", &T);
    CHECK("Tannealing=", "%" LF "", &Tannealing);
    CHECK("measure_PairDistrMATRIX=", "%i", &measure_PairDistrMATRIX);
    CHECK("PD_MATRIX_x=", "%" LF "", &PD_MATRIX_x);
    CHECK("PD_MATRIX_y=", "%" LF "", &PD_MATRIX_y);
    CHECK("gridPD_MATRIX_x=", "%i", &gridPD_MATRIX_x);
    CHECK("gridPD_MATRIX_y=", "%i", &gridPD_MATRIX_y);
    CHECK("gridSKT_k=", "%i", &gridSKT_k);
    CHECK("gridSKT_t=", "%i", &gridSKT_t);
    CHECK("SKT_dk=", "%" LF "", &SKT_dk);
    CHECK("SKT_kmin=", "%" LF "", &SKT_kmin);
    CHECK("sum_over_images_cut_off_num=", "%" LF "", &sum_over_images_cut_off_num);
    CHECK("file_append=", "%i", &file_append);
    CHECK("file_particles_append=", "%i", &file_particles_append);
    CHECK("measure_energy_barrier=", "%i", &measure_energy_barrier);
    CHECK("RDxmax=", "%" LF, &RD.max);
    CHECK("RDymax=", "%" LF, &RDz.max);
    CHECK("RadDistrMATRIXmaxx=", "%" LF, &RD_MATRIX.maxx);
    CHECK("RadDistrMATRIXmaxy=", "%" LF, &RD_MATRIX.maxy);
    CHECK("optimization=","%i", &optimization);
    CHECK("Niter_store=","%i", &Niter_store);
    CHECK("measure_SkMATRIX=","%i", &measure_SkMATRIX);
    CHECK("measure_FormFactor=","%i", &measure_FormFactor);
    CHECK("measure_Lindemann_Lozovik=","%i", &measure_Lindemann_Lozovik);
    CHECK("NkMaxTrap=", "%" LF, &NkMaxTrap);
    CHECK("generate_crystal_coordinates=","%i", &generate_crystal_coordinates);
    CHECK("measure_wavefunction_projection=","%i", &measure_wavefunction_projection);
    CHECK("inwf_power=", "%" LF, &inwf_power);
    CHECK("Vo=", "%" LF, &Vo);
    CHECK("Ro=", "%" LF, &Ro);
    CHECK("alpha_latt=", "%" LF, &alpha_latt);
    CHECK("beta=", "%" LF, &beta);
    CHECK("beta_latt=", "%" LF, &beta_latt);
    CHECK("gamma_latt=", "%" LF, &gamma_latt);
    CHECK("alpha_Jastrow=", "%" LF, &alpha_Jastrow);
    CHECK("delta_tilde=", "%" LF, &delta_tilde);
    CHECK("omega_x=", "%" LF, &omega_x);
    CHECK("omega_y=", "%" LF, &omega_y);
    CHECK("omega_z=", "%" LF, &omega_z);
    CHECK("S=", "%" LF "", &Srecoil);
    CHECK("Sy=", "%" LF "", &Syrecoil);
    CHECK("kL=", "%" LF, &kL);
    CHECK("solitonV=", "%" LF, &solitonV2);
    CHECK("solitonXi=", "%" LF, &solitonXi);
    CHECK("soliton_fixed_phase_V=",  "%" LF, &soliton_fixed_phase_V2);
    CHECK("soliton_fixed_phase_Xi=", "%" LF, &soliton_fixed_phase_Xi);
    CHECK("R2_subtract_CM=","%i", &R2_subtract_CM);
    CHECK("R3=", "%" LF, &R3);
    CHECK("bilayer_width=", "%" LF, &bilayer_width);
    CHECK("t_tunneling=", "%" LF, &t_tunneling);
    CHECK("Uo=", "%" LF, &Uo);
    CHECK("measure_effective_potential=", "%i", &measure_effective_potential);    
    CHECK("gaussian_alpha=", "%" LF, &gaussian_alpha);
    CHECK("gaussian_beta=", "%" LF, &gaussian_beta);
    CHECK("Rc_smooth_cutoff=", "%" LF, &Rc_smooth_cutoff);    
    CHECK("measure_pure_coordinates=", "%i", &measure_pure_coordinates);
    CHECK("measure_tVMC=", "%i", &measure_tVMC);    
    CHECK("beta_z=", "%" LF "", &beta_z);
    CHECK("RoSW=", "%" LF "", &RoSW);
    CHECK("Epotcutoff=", "%" LF "", &Epotcutoff);
    CHECK("c0R=", "%" LF "", &cR[0]);
    CHECK("c1R=", "%" LF "", &cR[1]);
    CHECK("c2R=", "%" LF "", &cR[2]);
    CHECK("c3R=", "%" LF "", &cR[3]);
    CHECK("c4R=", "%" LF "", &cR[4]);
    CHECK("c5R=", "%" LF "", &cR[5]);
    CHECK("c6R=", "%" LF "", &cR[6]);
    CHECK("c7R=", "%" LF "", &cR[7]);
    CHECK("c8R=", "%" LF "", &cR[8]);
    CHECK("c9R=", "%" LF "", &cR[9]);
    CHECK("c10R=", "%" LF "", &cR[10]);
    CHECK("c0I=", "%" LF "", &cI[0]);
    CHECK("c1I=", "%" LF "", &cI[1]);
    CHECK("c2I=", "%" LF "", &cI[2]);
    CHECK("c3I=", "%" LF "", &cI[3]);
    CHECK("c4I=", "%" LF "", &cI[4]);
    CHECK("c5I=", "%" LF "", &cI[5]);
    CHECK("c6I=", "%" LF "", &cI[6]);
    CHECK("c7I=", "%" LF "", &cI[7]);
    CHECK("c8I=", "%" LF "", &cI[8]);
    CHECK("c9I=", "%" LF "", &cI[9]);
    CHECK("c10I=", "%" LF "", &cI[10]);
    CHECK("tvmcNpar=", "%" LF "", &tvmcNpar);
    CHECK("tvmcNobs=", "%" LF "", &tvmcNobs);
    CHECK("Kpar11=", "%"LF, &Kpar11);
    CHECK("Kpar12=", "%"LF, &Kpar12);
    CHECK("aA=", "%"LF, &aA);
    CHECK("CMseparation=", "%"LF, &CMseparation);
    CHECK("measure_Hessian_matrix=", "%i", &measure_Hessian_matrix);
    CHECK("hard_core_diameter=", "%"LF, &hard_core_diameter);
    CHECK("Kurbakov_BC=", "%i", &Kurbakov_BC);
    CHECK("Ipar=", "%"LF, &Ipar);
    CHECK("check_kinetic_energy=", "%i", &check_kinetic_energy);
    CHECK("one_body_momentum_i=", "%"LF, &one_body_momentum_i);
    CHECK("L=", "%" LF "", &L);
    CHECK("measure_RDz=", "%i", &measure_RDz);    
    CHECK("Rsphere=", "%" LF "", &Rsphere);
    CHECK("sum_of_two_wavefunctions_parameter1=", "%" LF "", &sum_of_two_wavefunctions_parameter1);
    CHECK("sum_of_two_wavefunctions_parameter2=", "%" LF "", &sum_of_two_wavefunctions_parameter2);
    CHECK("sum_of_two_wavefunctions_alpha=", "%" LF "", &sum_of_two_wavefunctions_alpha);
    CHECK("numeric_derivative_dr=", "%" LF "", &numeric_derivative_dr);
  }

  a2 = a*a;
  lambda2 = lambda*lambda;
  lambda4 = lambda2*lambda2;
  alpha_x2 = alpha_x * alpha_x;
  alpha_y2 = alpha_y * alpha_y;
  alpha_z2 = alpha_z * alpha_z;
  two_alpha_x = 2. * alpha_x;
  two_alpha_y = 2. * alpha_y;
  two_alpha_z = 2. * alpha_z;
  omega_x2 = omega_x*omega_x;
  omega_y2 = omega_y*omega_y;
  omega_z2 = omega_z*omega_z;
  bilayer_width2 = bilayer_width*bilayer_width;

#ifdef TRAP_POTENTIAL
#ifdef TRIAL_3D
  Message("\n\n  Frequencies of the harmonic trap:\n");
  Message("    omega_x = %lf\n", omega_x);
  Message("    omega_y = %lf\n", omega_y);
  Message("    omega_z = %lf\n\n", omega_z);
#endif
#endif

  Message("\n\n  Type of the algorithm: ");
  if(MC == DIFFUSION) {
    if(SmartMC == DMC_QUADRATIC) Message("DMC: QUADRATIC\n");
    if(SmartMC == DMC_LINEAR_DRIFT) Message("DMC: LINEAR energy is calculated after the DRIFT\n");
    if(SmartMC == DMC_LINEAR_GAUSSIAN) Message("DMC: LINEAR energy is calculated after the GAUSSIAN jump\n");
    if(SmartMC == DMC_LINEAR_METROPOLIS) Message("DMC: LINEAR METROPOLIS\n");
    if(SmartMC == DMC_QUADRATIC_METROPOLIS) Message("DMC: QUADRATIC METROPOLIS\n");
    if(SmartMC == DMC_MOVE_ONE) Message("DMC: moving particles one by one\n");
    if(SmartMC == DMC_QUARTIC) Message("DMC: QUARTIC\n");
    if(SmartMC == DMC_PSEUDOPOTENTIAL) Message("DMC: PSEUDIOPOTENTIAL\n");
  }
  else if(MC == VARIATIONAL) { // VMC
    if(SmartMC == VMC_MOVE_ALL) Message("VMC: moving ALL particles\n");
    if(SmartMC == VMC_MOVE_ONE) Message("VMC: moving particles one by one\n");
    if(SmartMC == VMC_MOVE_DRIFT_ALL) Message("VMC: all particles are displaced by drift and a Gaussian jump\n");
    if(SmartMC == VMC_MOVE_DRIFT_ONE) Message("VMC: particles are displaced one by one by drift and a Gaussian jump\n");
    if(SmartMC == VMC_SUM_OF_TWO_WAVEFUNCTIONS) Message("VMC: sum of two wavefunctions\n");
  }
  else if(MC == CLASSICAL) {
    Message("classical\n");
  }
  else if(MC == PIMC) {
    Message("PIMC\n");
  }

  //Print type of the trial pair wave function
  Message("  Trial pair wavefunctions: ");
#ifdef TRIAL_EXTERNAL
  Message("loaded from file\n");
#endif
#ifdef TRIAL_HS
  Message("(3D HS+HS solution) + (exponential decay) -> hard spheres\n");
#endif
#ifdef TRIAL_HS_SIMPLE
  Message("(3D) HS + HS solution, no variational parameters \n");
#endif
#ifdef TRIAL_HS_SW_SIMPLE
  Message("(3D) HS SW solution, no variational parameters \n");
#endif
#ifdef TRIAL_SS
#ifdef TRIAL_3D
  Message("(3D SS+SS solution) + (exponential decay) -> hard spheres\n");
#else
  Message("(2D SS+SS solution), enable symmetrization!\n");
#endif
#endif
#ifdef TRIAL_SS_LARGE
  Message("(3D SS) + (exponential decay), with SS size Apar > Rpar\n");
#endif
#ifdef TRIAL_HS1D
  Message("(1D HS+HS solution) + (exponential decay) -> hard rodes\n");
#endif
#ifdef TRIAL_TONKS
  Message("Tonks wave function (1D) f = |Sin(sqrt{par}*(x-a))|\n");
#endif
#ifdef TRIAL_LIEB
  Message("Lieb wave function (1D) f = Cos(k(|z|-Rpar)) with f'(0)/f(0) = - 1/a_{1D}, for Rpar=0 -> Rpar = L/2\n");
#endif
#ifdef TRIAL_LIEB_LUTTINGER_PHONON
  Message("Lieb Luttinger phonon (1D) wave function, f = sin(Atrial + Btrial r + Ctrial r^2)^(1/Rpar), where f(0) = Apar and Rpar is equivalent to the Luttinger parameter K\n");
#endif
#ifdef TRIAL_LIEB_EXPONENTIAL_DECAY
  Message("Lieb exponential decay (1D) Atrial (r - a); 1 - Btrial Exp[-Apar r] - Btrial Exp[-Apar (L - r)], parameters Rpar and Apar\n");
#endif
#ifdef TRIAL_LIM
  Message("f = 1 - a/r\n");
#endif
#ifdef TRIAL_PP_FS_TRAP_MCGUIRE1D
  Message("f(r) = exp(1 - (Apar**2*b - a*r**2)/(a*Apar*b + a*b*r))/r, i.e. a: 3D scattering length, b: 1D scattering length, Apar - var. par\n");
#endif
#ifdef TRIAL_POWER
  Message("f = 1 - a/r^Apar\n");
#endif
#ifdef TRIAL_PHONON
  Message("phonon wavefunction (1D) A Cos(k(x-B)) (x<L*Rpar), sin^a(pi*x/L), i.e. Lieb + phonons matching, parameters: Rpar\n");
#endif
#ifdef TRIAL_PHONON_LUTTINGER
  Message("phonon Luttinger wavefunction (1D) A Cos(k(x-B)) (x<L*Rpar), sin^a(pi*x/L), i.e. Lieb + phonons matching, parameters: Kpar, Rpar\n");
#endif
#ifdef TRIAL_PHONON_LUTTINGER_LATTICE
  Message("phonon Luttinger lattice wavefunction (1D) Kpar, Rpar + Apar12, Rpar12, Bpar12. Note that zero boundary condition might be changed\n");
#endif
#ifdef TRIAL_PHONON_LUTTINGER_PIECEWISE
  Message("phonon Luttinger wavefunction piecewise\n");
#endif
#ifdef TRIAL_TONKS_TRAP
  Message("Tonks-Girardeau wavefunction in trap (1D) f(z) = z - a\n");
#endif
#ifdef TRIAL_TONKS_TRAP11
  Message("11: Tonks-Girardeau wavefunction in trap f(r) = |r - b|\n");
#endif

#ifdef HARD_SPHERE
  Warning("Hard-sphere check is enabled\n");
#ifdef HARD_SPHERE_KILL_OVERLAPPED
  Message("When hard-spheres overlap in DMC: KILL the walker\n");
#else
  Message("When hard-spheres overlap in DMC: RETHROW the walker\n");
#endif
#else
  Warning("Hard-sphere check is NOT enabled\n");
#endif

  // Print number of dimensions
  Message("  Number of dimensions: ");
#ifdef TRIAL_1D
  Message("1D\n");
#endif
#ifdef TRIAL_2D
  Message("2D\n");
#endif
#ifdef TRIAL_3D
  Message("3D\n");
#endif

  Message("  Boundary conditions: ");
#ifdef  BC_ABSENT
  Message("absent (external confinement)\n");
#endif
#ifdef BC_1DPBC_Z
  Message("1D Periodic boundary conditions (in Z direction)\n");
#endif
#ifdef BC_1DPBC_X
  Message("1D Periodic boundary conditions (in X direction)\n");
#endif
#ifdef BC_2DPBC
  Message("2D Periodic boundary conditions\n");
#endif
#ifdef BC_2DPBC_HEXAGON
  Message("box geometry: HEXAGON\n");
#endif
#ifdef BC_3DPBC
  Message("3D Periodic boundary conditions\n");
#endif
#ifdef BC_3DPBC_TRUNCATED_OCTAHEDRON
  Message("box geometry: TRUNCATED OCTAHEDRON\n");
#endif

  Message("done\n");

  acceptance_rate /= 100.;

#ifdef TRIAL_SUTHERLAND
  Message("Calogero-sutherland wavefunction\n");
#endif
#ifdef TRIAL_CALOGERO
  Message("Calogero wavefunction\n");
#endif

#ifndef BC_TRAP
  if(alpha_y < 0) { // beta = 1 / 4 lambda <z2>
#ifdef TRIAL_TONKS
    alpha_y = 1 / 4*N;
    Warning("Automaticaly setting TONKS value for parameter    alpha_y = %" LF "\n", alpha_y);
#else
    alpha_y = 5./4. * lambda*pow(3*(N-1)*lambda*a, -2./3.);
    if(alpha_y>0.5) alpha_y = 0.5;
    Warning("Automaticaly setting MEAN-FIELD value for parameter\n  alpha_y = %" LF "\n", alpha_y);
#endif
  }
#endif

  if(Npop == 1) {
    branchng_present = OFF;
    Warning(" Branchnig is absent.\n");
  }

  if(Npop>NwalkersMax) {
    if(MC == DIFFUSION) {
      Warning("Npop > NwalkersMax (%i>%i), adjusting NmalkersMax to %i\n", Npop, NwalkersMax, Npop*2);
      NwalkersMax = 2*Npop;
    }
    else {
      Warning("Npop > NwalkersMax (%i>%i), adjusting NmalkersMax to %i\n", Npop, NwalkersMax, Npop);
      NwalkersMax = Npop;
    }
  }

  if(boundary == NO_BOUNDARY_CONDITIONS || boundary == ONE_BOUNDARY_CONDITION) {
#ifdef TRIAL_HS
    Warning("Automaticaly converting Rpar=%" LF "[a] => Rpar=%" LF "[aosc]\n", Rpar, Rpar*a);
    Rpar *= a;
#endif
#ifdef TRIAL_SS
    if(Apar > a) Warning("WARNING: Apar is larger than a\n");
    Warning("Automaticaly converting Apar=%" LF "[a] => Apar=%" LF "[aosc]\n", Apar, Apar*a);
    Apar *= a;
    Warning("Automaticaly converting Rpar=%" LF "[a] => Rpar=%" LF "[aosc]\n", Rpar, Rpar*a);
    Rpar *= a;
#endif
  }

#ifdef BC_ABSENT
  if(boundary != NO_BOUNDARY_CONDITIONS) Error("Cannot simulate free system: define BC_ABSENT and recompile the program\n");
//  Warning("Automatically converting dt=%" LE " [a] => dt=%" LE " [aosc]\n", dt, dt*a*a);
//  dt *= a*a;
#ifdef TRIAL_2D
  if(lambda != 0) {
    Warning("Automatically setting lambda = 0\n");
    lambda = 0;
  }
#endif
#endif

#ifdef BC_1DPBC_ANY
  if(boundary != ONE_BOUNDARY_CONDITION) Error("Cannot simulate system with 1D PBC: define BC_1DPBC_X or BC_1DPBC_Y and recompile the program\n");
  if(lambda != 1.) {
    Warning("Forcing lambda = 1\n");
    lambda = 1.;
  }
#endif

#ifdef BC_2DPBC
  if(boundary != TWO_BOUNDARY_CONDITIONS) Error("Cannot simulate system with 2D PBC: define BC_2DPBC and recompile the program\n");
  if(lambda != 1.) {
    Warning("Forcing lambda = 1\n");
    lambda = 1.;
  }
#endif

#ifdef BC_3DPBC
  if(boundary != THREE_BOUNDARY_CONDITIONS) Error("Cannot simulate system with 3D PBC: define BC_3DPBC and recompile the program\n");
  if(a != 1.) {
    Warning("Note that a differs from 1\n");
    //a = 1.;
  }
#endif

#ifdef TRIAL_1D
  if(SKT_dk<0) {
    SKT_dk = 2./(DOUBLE)(N);
    Message("  Adjusting SKT_dk to %lf\n", SKT_dk);
  }
#endif

#ifndef __WATCOMC__
#ifndef LINUX_GRAPHICS
#ifndef __MVC__
  video = OFF;
#endif
#endif
#endif

#ifdef SCALABLE
  Ncells = (int)(pow((DOUBLE)Nscal,DIMENSION)+0.1);
  Message("  Scalability is enabled: number of cells %i\n", Ncells);
  Warning("  Changing number of particles %i -> %i\n", N, N*Ncells);
  N *= Ncells;
#ifdef SCALABLE_POTENTIAL
  Message("  Potential energy is treated in a scalable way\n");
#else
  Message("  Potential energy is treated in an exact way (scalability is not full)\n");
#endif
#endif

#ifdef LATTICE
#ifdef TRIAL_2D
  Nlattice = (pow(NCrystal, 1./2.)+1e-12);
  Nlattice3 = Nlattice*Nlattice;
  if(Nlattice3 != NCrystal) Error("Cannot construct square lattice with %i lattice sites!", NCrystal);
  Message("  Size of the square lattice is %ix%i\n", Nlattice, Nlattice);
#endif
#ifdef TRIAL_3D
  Nlattice = (int)(pow(NCrystal, 1./3.)+1e-12);
  Nlattice3 = Nlattice*Nlattice*Nlattice;
  if(Nlattice3 != NCrystal) Error("Cannot construct cubic lattice with %i lattice sites!", NCrystal);
  Message("  Size of the cubic lattice is %ix%ix%i\n", Nlattice, Nlattice, Nlattice);
#endif
#endif

#ifdef VEXT_COS2
#ifdef BC_3DPBC_CUBE
  Nlattice = (int)(pow(NCrystal, 1./3.)+1e-12);
  Nlattice3 = Nlattice*Nlattice*Nlattice;
  if(Nlattice3 != NCrystal) Error("Cannot construct cubic lattice with %i lattice sites!", NCrystal);
  Message("  Size of the cubic lattice is %ix%ix%i\n", Nlattice, Nlattice, Nlattice);
#endif
#ifdef BC_2DPBC_SQUARE
  Nlattice = (pow(NCrystal, 1./2.)+1e-12);
  Nlattice3 = Nlattice*Nlattice;
  if(Nlattice3 != NCrystal) Error("Cannot construct square lattice with %i lattice sites!", NCrystal);
  Message("  Size of the square lattice is %ix%i\n", Nlattice, Nlattice);
#endif
#ifdef BC_1DPBC_Z
  if(NCrystal<0)
    Nlattice = N;
  else
    Nlattice = NCrystal;
  Message("  Size of the square lattice is %i\n", Nlattice);
#endif
#endif

#ifdef VEXT_TRAP_DOUBLE_WELL
  Message("  Double well potential:\n");
  Message("    alpha_latt = %"LF"\n", alpha_latt);
#endif  

#ifdef CRYSTAL
#ifdef CRYSTAL_NONSYMMETRIC
  Message("  Simulation of the nonsymmetric crystal phase\n");
#endif
#ifdef CRYSTAL_SYMMETRIC
  Message("  Simulation of the symmetric crystal phase\n");
#endif
  
  if(NCrystal != -1) Message("  Crystal.size is %i\n", NCrystal);
  if(NCrystal == -1) NCrystal = N;
#endif

  if(Ndens == -1) {
    Ndens = N;
    Message("  Size of the box will be calculated using Ndens = %i\n", Ndens);
  }

  if(measure_Nk_pure == ON && measure_OBDM && MC == DIFFUSION) Message("  Pure measurement: momentum distribution has been chosen\n");
  if(measure_Nk_pure == OFF && measure_OBDM && MC == DIFFUSION) Message("  Pure measurement: OBDM has been chosen\n");

#ifdef INTERACTION_LENNARD_JONES10
  Warning("  changing Apar to satisfy boundary condition -> %lf\n", sqrt(D)/4.);
  Apar = sqrt(D)/4.;
#endif

#ifdef INTERACTION_LOGARITHMIC
  Message("  Interaction potential logarithmic\n");
#endif

#ifndef TRIAL_1D
  if(measure_FormFactor && !measure_Sk) {
    Warning("  Enabling measure_Sk in order to proceed with measurements of Form Factor\n");
    measure_Sk = ON;
  }
#endif

#ifdef MEASURE_CORRELATION_FUNCTIONS_IN_THE_MIDDLE_OF_DRIFT
  Message("  Correlation functions are measured: in the middle point of the drift (Rpp)\n");
#else
  Message("  Correlation functions are measured: at the end of the drift (Rppp)\n");
#endif

#ifdef MEASURE_PURE_OBDM
  Message("  Measurement of pure OBDM or nk enabled\n");
#endif
#ifdef MEASURE_PURE_TBDM
  Message("  Measurement of pure TBDM enabled\n");
#endif

#ifdef ONE_BODY_SOLITON
  Message("  soliton velocity %lf [in units of speed of sound]\n", solitonV2);
  solitonV2 = solitonV2*solitonV2;
  solitonk = (sqrt(1.-solitonV2)/(sqrt(2.)*solitonXi));
#endif

#ifdef VEXT_SOLITON_FIXED_PHASE
  Message("  Vext: soliton velocity %lf [in units of speed of sound]\n", soliton_fixed_phase_V2);
  soliton_fixed_phase_V2 = soliton_fixed_phase_V2*soliton_fixed_phase_V2;
#endif

#ifdef VEXT_IMPURITY_3DSW
  Message("  Vext: Square Well impurity");
#endif

#ifdef SPINLESS
  Nspin = 1;
#endif

  if(Nspin>1)
    Message("  Number of spin components: %i\n", Nspin);
  else
    Message("  One component system (no spin)\n");

#ifdef TRIAL_1D
 if(measure_RadDistr) measure_RDz = ON;
#endif

  if(MC != DIFFUSION) measure_SD = OFF;

  fclose(in);

  return 0;
}

/************************* Load Particle Coordinates *************************/
int LoadParticleCoordinates(void) {
  FILE *in;
  int w, i;
  DOUBLE Ew;

  if(verbosity) Warning("Loading walker coordinates ...\n  ");

  in = fopen(INPATH "in3Dprev.in", "r");
  if(in == NULL) {
    perror("\nError:");
    Error("\ncan't read %s", file_particles);
  }

  // Define number of walkers
  fscanf(in, "%i\n", &Nwalkers);
  if(Nwalkers>NwalkersMax) {
    Warning("  Increasing number of maximal walkers in order to be able to load saved coordinates! %i -> %i\n", NwalkersMax, Nwalkers);
    NwalkersMax = Nwalkers;
  }

  AllocateWalkers();

  E = 0.;
  for(w=0; w<Nwalkers; w++) {
    fscanf(in, "%" LF "\n", &Ew);

#ifdef RESCALE_COORDINATES_TO_OSC_UNITS // rescale to oscillator units
    Ew /= (2.*a*a);
#endif
    //Ew /= (DOUBLE) N;
    W[w].E = Ew;
    E += Ew;
    for(i=0; i<N; i++) {
      fscanf(in, "%" LF " %" LF " %" LF  "\n", &W[w].x[i], &W[w].y[i], &W[w].z[i]);

#ifdef RESCALE_COORDINATES_TO_OSC_UNITS // Coordinates are stored in scattering length distances
      W[w].x[i] *= a;
      W[w].y[i] *= a;
      W[w].z[i] *= a;
#endif

#ifdef TRIAL_1D
      if(W[w].x[i] != 0 || W[w].y[i] != 0) {
        Warning("1D simulation: forcing coordinates %i %i " LF LF " be equal to zero\n", w+1, i+1, W[w].x[i], W[w].y[i]);
        W[w].x[i] = W[w].y[i] = 0.;
      }
#endif

#ifdef TRIAL_2D
      if(W[w].z[i] != 0) {
        Warning("2D simulation: forcing coordinates %i %i %" LF " be equal to zero\n", w+1, i+1, W[w].z[i]);
        W[w].z[i] = W[w].z[i] = 0.;
      }
#endif
      // inside of the box check
      if((boundary == ONE_BOUNDARY_CONDITION && (W[w].z[i]>Lz || W[w].z[i]<0))
      || (boundary == THREE_BOUNDARY_CONDITIONS && (W[w].x[i]>Lx || W[w].x[i]<0 
          || W[w].y[i]>Ly || W[w].y[i]<0 || W[w].z[i]>Lz || W[w].z[i]<0))
          ||(boundary == TWO_BOUNDARY_CONDITIONS && (W[w].x[i]>Lx || W[w].x[i]<0 
          || W[w].y[i]>Ly || W[w].y[i]<0))) {
        Warning("Walker %i, particle %i\n", w+1, i+1);
        Warning("coords %" LF " %" LF " %" LF ", L = %" LF "\n", W[w].x[i], W[w].y[i], W[w].z[i], L);
        Error("bad particle coordinate\n");
      }
    }
    if(CheckWalkerOverlapping(W[w])) 
      Error("Bad initial configuration, particles overlap, check %s file", file_particles);

  }

  E /= (DOUBLE) Nwalkers;
  Nwalkersw = (DOUBLE) Nwalkers;
  fclose(in);

#ifdef SPINFULL // load spins
  Message("  Loading spins ... ");
  in = fopen(INPATH "in3Dspin.in", "r");
  if(in == NULL) {
    perror("\nError:");
    Error("\ncan't read %s", "in3Dspin.in");
  }
  else{
    for(w=0; w<Nwalkers; w++) {
      for(i=0; i<N; i++) {
        fscanf(in, "%i\n", &W[w].spin[i]);
        if(W[w].spin[i] > Nspin-1) Error("  spin=%i of particle i=%i is out of Nspin = %i range!", W[w].spin[i]+1, i, Nspin);
      }
    }
    fclose(in);
  }
  Message("done\n");
#endif

#pragma omp parallel for
  for(w=0; w<Nwalkers; w++) { // once spins are defined, calculate U
    W[w].U = U(W[w]);
  }

  /*if(measure_SD) {
    in = fopen("inSDprev.in", "r");
    if(in == NULL) {
      Warning("  Cannot restore SD information for the particles\n");
    }
    else {
      for(w=0; w<Nwalkers; w++) {
        for(i=0; i<N; i++) fscanf(in, "%" LF " %" LF " %" LF  "\n", &W[w].rreal[i][0], &W[w].rreal[i][1], &W[w].rreal[i][2]);
        for(i=0; i<SD.size; i++) fscanf(in, "%" LF " %" LF " %" LF  "\n", &W[w].CM[0][i], &W[w].CM[1][i], &W[w].CM[2][i]);
      }
      fclose(in);
    }
  }*/

#ifdef LATTICE
  // loading and converting!
  //Warning("changing number of particles %i -> %i\n", N, N/2);
  //N = N/2;
  //SaveCoordinates();
#endif

#ifdef HARD_SPHERE
  for(w=0; w<Nwalkers; w++) if(CheckWalkerOverlapping(W[w])) Error("  particles overlap, walker %i\n", w);
#endif

#ifdef HARD_SPHERE_HYPERRADIUS
  for(w=0; w<Nwalkers; w++) if(CheckWalkerHyperradiusOverlapping(W[w])) Error("  particles hyperradius overlap, walker %i\n", w);
#endif

  if(verbosity) Warning("done\n");
  return 0;
}

/****************************** Save Coordinates ************************/
int SaveCoordinates(char *file_particles) {
  FILE *out;
  int w, i;
  static int times_saved = 0;
  char filename[20];
#ifdef SPINFULL
  FILE *file_spin, *file_coord_spin;
#endif

#ifdef SECURE
  if(CheckOverlapping()) Error("Save Coordinates : particles overlap");
#endif

  if(file_particles_append<2)
    out = fopen(file_particles, file_particles_append?"a":"w");
  else {
    sprintf(filename, "in3Dprev.%i", times_saved++);
    out = fopen(filename, "w");
  }

  if(out == NULL) {
    perror("\nError:");
    Warning("can't write to %s file\n", file_particles);
    return 1;
  }
#ifdef SPINFULL
  file_spin = fopen("in3Dspin.in", file_particles_append?"a":"w");
  file_coord_spin = fopen("coordspin.in", file_particles_append?"a":"w");
#endif

  if(file_particles_append != 1) fprintf(out, "%i\n", Nwalkers);

  for(w=0; w<Nwalkers; w++) {
#ifdef RESCALE_COORDINATES_TO_OSC_UNITS // rescale to oscillator units
    fprintf(out, "%.15" LE "\n", 2.*a*a*WalkerEnergy0(&W[w]));
#else
    fprintf(out, "%.15" LE "\n", WalkerEnergy0(&W[w]));
#endif

    for(i=0; i<N; i++) {
#ifdef RESCALE_COORDINATES_TO_OSC_UNITS // rescale to scattering length
      fprintf(out, "%.15" LF " %.15" LF " %.15" LF "\n", W[w].x[i]/a, W[w].y[i]/a, W[w].z[i]/a);
#else  
      fprintf(out, "%.15" LE " %.15" LE " %.15" LE "\n", W[w].x[i], W[w].y[i], W[w].z[i]);
#ifdef SPINFULL
      fprintf(file_spin, "%i\n", W[w].spin[i]);
      fprintf(file_coord_spin, "%.15" LE " %.15" LE " %.15" LE " %i\n", W[w].x[i], W[w].y[i], W[w].z[i], W[w].spin[i]);
#endif
#endif
    }
#ifdef SPINFULL
    fprintf(file_spin, "\n");
    fprintf(file_coord_spin, "\n\n");
#endif
  }
  fclose(out);
#ifdef SPINFULL
  fclose(file_spin);
  fclose(file_coord_spin);
#endif

  /*if(measure_SD) {
    out = fopen("inSDprev.in", "w");
    if(out == NULL) {
      perror("\nError:");
      Warning("can't save SD data\n");
      return 1;
    }
    for(w=0; w<Nwalkers; w++) {
      for(i=0; i<N; i++) fprintf(out, "%.15" LE " %.15" LE " %.15" LE "\n", W[w].rreal[i][0], W[w].rreal[i][1], W[w].rreal[i][2]);
      for(i=0; i<SD.size; i++) fprintf(out, "%.15" LE " %.15" LE " %.15" LE "\n", W[w].CM[0][i], W[w].CM[1][i], W[w].CM[2][i]);
    }
    fclose(out);
  }*/

  return 0;
}

/******************************** Save Energy VMC ************************/
int SaveEnergyVMC(DOUBLE E, DOUBLE EFF) {
  FILE *out;
  static int first_time = ON;

#ifdef MPI
  ParallelSaveEnergyVMC();
  return 0;
#endif

  out = fopen(file_energy, (file_append || !first_time)?"a":"w");
  if(out == NULL) {
    perror("\nError:");
    Warning("can't write to energy file %s\n", file_energy);
    return 1;
  }

  if(first_time && !file_append) {
    fprintf(out, "#E EFF Epot Ekin Eint"); // NB for delta-pseudopotential interaction Ekin = EFF-Epot, Edelta = E-EFF
#ifdef INTERACTION_WITH_DAMPING
    fprintf(out, " Edamping");
#endif
#ifdef SPINFULL
    fprintf(out, " spin");
#endif
#ifdef EXTERNAL_POTENTIAL
    fprintf(out, " Eext");
#endif
    fprintf(out, "\n");
    first_time = OFF;
  }

  fprintf(out, "%.15" LE " %.15" LE " %.15" LE " %.15" LE " %.15" LE, E*energy_unit, EFF*energy_unit, Epot*energy_unit, Ekin*energy_unit, Eint*energy_unit);
#ifdef INTERACTION_WITH_DAMPING
  fprintf(out, " " %.15" LE, Edamping*energy_unit);
#endif
#ifdef SPINFULL
  fprintf(out, " %lf", SpinPolarization());
#endif
#ifdef EXTERNAL_POTENTIAL
  fprintf(out, " %.15" LE, Eext*energy_unit);
#endif

  fprintf(out, "\n");

  fclose(out);

  return 0;
}

/******************************** Save Energy Cls ************************/
int SaveEnergyCls(DOUBLE E) {
  FILE *out;
  static long int i = 0;  // variable i tells how many times the energy was saved 

  out = fopen(file_energy, (file_append || i++)?"a":"w");
  if (out == NULL) {
    perror("\nError:");
    Warning("can't write to energy file %s\n", file_energy);
    return 1;
  }

  fprintf(out, "%.15" LE "\n", E*energy_unit);

  fclose(out);
  return 0;
}

/******************************** Save Energy DMC ************************/
int SaveEnergyDMC(DOUBLE E) {
  FILE *out;
  static int first_time = ON;

#ifdef MPI
//  ParallelSaveEnergyDMC(E);
//  return 0;
  if(myid) return 0;
#endif

  out = fopen(file_energy, (file_append || !first_time)?"a":"w");
  if(out == NULL) {
    perror("\nError:");
    Warning("can't write to energy file %s\n", file_energy);
    return 1;
  }

  if(first_time && !file_append) {
    fprintf(out, "#E ");
    if(branchng_present == OFF) fprintf(out, " EFF");
    fprintf(out, " Nwalkers Epot Ekin Eint");
 
#ifdef EXTERNAL_POTENTIAL
    fprintf(out, " Eext");
#endif
#ifdef INTERACTION_WITH_DAMPING
    fprintf(out, " Edamping");
#endif
#ifdef MEASURE_ENERGY_FROM_GROWTH_ESTIMATOR
    fprintf(out, " Efrombranching");
#endif
#ifdef SPINFULL
    fprintf(out, " spin");
#endif
    fprintf(out, "\n");
    first_time = OFF;
  }

  fprintf(out, "%.15e", energy_unit*E);
  if(branchng_present == OFF) fprintf(out, " %.15e", energy_unit*EFF);
  fprintf(out, " %i %e %e %e", Nwalkers, energy_unit*Epot, energy_unit*Ekin, energy_unit*Eint);
#ifdef EXTERNAL_POTENTIAL
  fprintf(out, " %e", energy_unit*Eext);
#endif
#ifdef INTERACTION_WITH_DAMPING
  fprintf(out, "%.15" LE " ", Edamping*energy_unit);
#endif
#ifdef MEASURE_ENERGY_FROM_GROWTH_ESTIMATOR
  fprintf(out, " %e", E_from_growth_estimator);
#endif

#ifdef SPINFULL
   fprintf(out, " %lf", SpinPolarization());
#endif
   fprintf(out, "\n");

  fclose(out);

  return 0;
}

/*********************** Save Pair Distribution ******************************/
// 4.*pi/3.*((r+dr)^3-r^3)
// 0,
// 2 pi r dr
// pi*((r+dr)^2-r^2)
int SavePairDistribution(void) {
  int i;
  DOUBLE r,dr;
  FILE *out, *out2 = NULL;
  int save_pure = ON;
  DOUBLE normalization, normalization_pure;
  static int initialized = OFF;
  DOUBLE V, combinations;

#ifdef SPINFULL
  int m;
  char PD_dat_spinfull[100];
  char PD_dat_spinfull_pure[100];
  FILE *out3, *out4 = NULL;
#endif 

  dr = PD.step;

  if(PD.times_measured == 0) {
    Warning("Attempt to save pair distribution function without any measurement\n");
    return 1;
  }

  if(measure_PairDistrMATRIX) SavePairDistribtuionMatrix();
  if(measure_PairDistr_R2) SavePairDistribtuionR2();

  out = fopen(file_PD, (initialized || file_append)?"a":"w");
  if(out == NULL) {
    perror("\nError:");
    Warning("can't save distribution function to file %s\n", file_PD);
    return 1;
  }

  if(MC != DIFFUSION) { // VMC or Classical MC
    save_pure = OFF;
  }

  if(save_pure) {
    if(PD_pure.times_measured == 0) {
      Warning("Attempt to save pure pair distribution function without any measurement\n");
      save_pure = OFF;
    }
    out2 = fopen(file_PD_pure, (initialized || file_append)?"a":"w");
    if(out2 == NULL) {
      perror("\nError:");
      Warning("can't save distribution function to file %s\n", file_PD_pure);
      return 1;
    }
  }

  V = Ndens / n;
  combinations = (DOUBLE) (((N-1)*N) / 2);

#ifdef BC_ABSENT
  normalization = 1. / ((DOUBLE)(PD.times_measured*N)*PD.step);
  if(MC == DIFFUSION) normalization_pure = 1. / ((DOUBLE)(PD_pure.times_measured*N)*PD.step);
#else // no trap
#ifdef TRIAL_1D
  //normalization = 1./((DOUBLE)(PD.times_measured)*combinations/V*PD.step);
  normalization = 1. / ((DOUBLE)(PD.times_measured*N)*PD.step);
  //if(MC == DIFFUSION) normalization_pure = 1./((DOUBLE)(PD_pure.times_measured)*combinations/V*PD.step);
  if(MC == DIFFUSION) normalization_pure = 1./((DOUBLE)(PD_pure.times_measured*N)*PD.step);
#endif

#ifdef TRIAL_2D // cell contents corresponds to average over (R^2-(R-dr)^2)
  normalization = 1./((DOUBLE)(PD.times_measured)*combinations/V*2.*PI*PD.step);
  if(MC == DIFFUSION) normalization_pure = 1./((DOUBLE)(PD_pure.times_measured)*combinations/V*2.*PI*PD.step);
#endif

#ifdef TRIAL_3D
  normalization = 1./((DOUBLE)(PD.times_measured)*combinations/V*4.*PI*PD.step);
  if(MC == DIFFUSION) normalization_pure = 1./((DOUBLE)(PD_pure.times_measured)*combinations/V*4.*PI*PD.step);
#endif
#endif

#ifdef LATTICE_ZIGZAG // 1D g2(r)
  normalization = 1. / ((DOUBLE)(PD.times_measured*N)*n*PD.step);
  if(MC == DIFFUSION) normalization_pure = 1. / ((DOUBLE)(PD_pure.times_measured*N)*n*PD.step);
#endif

   for(i=0; i<PD.size; i++) {
     r = PD.min + (DOUBLE)(i + 0.5) * PD.step;
#ifdef BC_ABSENT
     fprintf(out, "%" LF " %.15" LE "\n", r, normalization*(DOUBLE)PD.N[i]);
     if (save_pure) fprintf(out2, "%" LF " %.15" LE "\n", r, normalization_pure*(DOUBLE)PD_pure.N[i]);
#else // no trap
#ifdef TRIAL_1D
     fprintf(out, "%" LF " %.15" LE "\n", r*n, normalization*(DOUBLE)PD.N[i]);
     if (save_pure) fprintf(out2, "%" LF " %.15" LE "\n", r*n, normalization_pure*(DOUBLE)PD_pure.N[i]);
#endif

#ifdef TRIAL_2D
#ifndef LATTICE_ZIGZAG // 2D
     fprintf(out, "%" LF " %.15" LE "\n", r, normalization*(DOUBLE)PD.N[i] / r);
     if (save_pure) fprintf(out2, "%" LF " %.15" LE "\n", r, normalization_pure*(DOUBLE)PD_pure.N[i] / r);
#else // quasi 2D
      fprintf(out, "%" LF " %.15" LE "\n", r*n, normalization*(DOUBLE)PD.N[i]);
      if (save_pure) fprintf(out2, "%" LF " %.15" LE "\n", r*n, normalization_pure*(DOUBLE)PD_pure.N[i]);
#endif
#endif

#ifdef TRIAL_3D
      //fprintf(out, "%" LF " %.15" LE "\n", r, normalization*(DOUBLE)PD.N[i]/(r*r));
      //if(save_pure) fprintf(out2, "%" LF " %.15" LE "\n", r, normalization_pure*(DOUBLE)PD_pure.N[i]/(r*r));
      fprintf(out, "%" LF " %.15" LE "\n", r, normalization*(DOUBLE)PD.N[i] / (r*r));
      if(save_pure) fprintf(out2, "%" LF " %.15" LE "\n", r, normalization_pure*(DOUBLE)PD_pure.N[i] / (r*r));
      //fprintf(out, "%" LF " %.15" LE "\n", r+dr/2., normalization*4*PI*dr/(4.*PI/3.*((r+dr)*(r+dr)*(r+dr)-r*r*r))*(DOUBLE)PD.N[i]);
      //if(save_pure) fprintf(out2, "%" LF " %.15" LE "\n", r+dr/2., normalization_pure*4*PI*dr/(4.*PI/3.*((r+dr)*(r+dr)*(r+dr)-r*r*r))*(DOUBLE)PD_pure.N[i]);
#endif
#endif

#ifdef HARD_SPHERE
     if (r < a2 && PD.N[i]>0)  Warning("  unreal pair distribution function\n");
#endif

      PD.N[i] = 0;
      if(save_pure) PD_pure.N[i] = 0;
   }
    fclose(out);
    if(save_pure) fclose(out2);

#ifdef SPINFULL
    for(m=0; m<Nspin; m++) {
      if(m == 0) {
        combinations = (DOUBLE)(N/Nspin*(N/Nspin-1)/2.);
      }
      else { 
        combinations = (DOUBLE)(N / Nspin*(N / Nspin)); 
      }
#ifdef BC_ABSENT
      normalization = 1. / ((DOUBLE)(PD.times_measured*N)*PD.step);
      if(MC == DIFFUSION) normalization_pure = 1. / ((DOUBLE)(PD_pure.times_measured*N)*PD.step);
#else // no trap
#ifdef TRIAL_1D
      normalization = 1. / ((DOUBLE)(PD.times_measured)*combinations / V*PD.step);
      if(MC == DIFFUSION) normalization_pure = 1. / ((DOUBLE)(PD_pure.times_measured)*combinations / V*PD.step);
#endif
#ifdef TRIAL_2D // cell contents corresponds to average over (R^2-(R-dr)^2)
      normalization = 1. / ((DOUBLE)(PD.times_measured)*combinations / V*2.*PI*PD.step);
      if(MC == DIFFUSION) normalization_pure = 1. / ((DOUBLE)(PD_pure.times_measured)*combinations / V*2.*PI*PD.step);
#endif
#ifdef TRIAL_3D
      normalization = 1. / ((DOUBLE)(PD.times_measured)*combinations / V*4.*PI*PD.step);
      if(MC == DIFFUSION) normalization_pure = 1. / ((DOUBLE)(PD_pure.times_measured)*combinations / V*4.*PI*PD.step);
#endif
#endif
      sprintf(PD_dat_spinfull, "outpdspin%i.dat", m);
      out3 = fopen(PD_dat_spinfull, (initialized || file_append) ? "a" : "w");
      if(out3 == NULL) {
        perror("\nError:");
        Warning("can't save distribution SPINFULL function to file %s\n", PD_dat_spinfull);
        return 1;
      }

      if(save_pure) {
        sprintf(PD_dat_spinfull_pure, "outpdpspin%i.dat", m);
        out4 = fopen(PD_dat_spinfull_pure, (initialized || file_append) ? "a" : "w");
        if(out4 == NULL) {
          perror("\nError:");
          Warning("can't save distribution SPINFULL function to file %s\n", PD_dat_spinfull_pure);
          return 1;
        }
      }

      for(i=0; i<PD.size; i++) {
        r = PD.min + (DOUBLE)(i + 0.5) * PD.step;
#ifdef TRIAL_1D
        fprintf(out3, "%" LF " %.15" LE "\n", r, normalization*(DOUBLE)PD.PDSpin[m][i]);
        if(save_pure) fprintf(out4, "%" LF " %.15" LE "\n", r, normalization_pure*(DOUBLE)PD_pure.PDSpin[m][i]);
#endif
#ifdef TRIAL_2D
        fprintf(out3, "%" LF " %.15" LE "\n", r, normalization*(DOUBLE)PD.PDSpin[m][i] / r);
        if(save_pure) fprintf(out4, "%" LF " %.15" LE "\n", r, normalization_pure*(DOUBLE)PD_pure.PDSpin[m][i] / r);
#endif
#ifdef TRIAL_3D
        fprintf(out3, "%" LF " %.15" LE "\n", r, normalization*(DOUBLE)PD.PDSpin[m][i] / (r*r));
        if(save_pure) fprintf(out4, "%" LF " %.15" LE "\n", r, normalization_pure*(DOUBLE)PD_pure.PDSpin[m][i] / (r*r));
#endif
        PD.PDSpin[m][i] = 0;
        if(save_pure) PD_pure.PDSpin[m][i] = 0;
      }
      fclose(out3);
      if(save_pure) fclose(out4);
    }
#endif

#ifdef TRIAL_1D // measure g3(0)
    // save g3(0)
    out = fopen("outg3.dat", (initialized || file_append)?"a":"w");
    if(out == NULL) {
      perror("\nError:");
      Warning("can't save distribution function to file %s\n", "outg3.dat");
      return 1;
    }

#ifdef BC_1DPBC_ANY
    fprintf(out, "%" LG " %.15" LE "\n", PD.step, 6*(DOUBLE)PD.g3 /((DOUBLE) PD.times_measured*N*n*n*4*PD.step*PD.step));
#endif

#ifdef BC_2DPBC
    fprintf(out, "%.15" LE "\n", 6*(DOUBLE)PD.g3 /((DOUBLE) PD.times_measured*N*n*n*PI*PI*PD.step*PD.step*PD.step*PD.step));
#endif

#ifdef BC_3DPBC
    fprintf(out, "%.15" LE "\n", 6*(DOUBLE)PD.g3 /((DOUBLE) PD.times_measured*N*n*n*16./9.*PI*PI*PD.step*PD.step*PD.step * PD.step*PD.step*PD.step));
#endif

    PD.g3 = 0;
    fclose(out);
#endif

  PD.times_measured = 0;
  if(MC == DIFFUSION) PD_pure.times_measured = 0;
  initialized = ON;

  return 0;
}

/*********************** Save Piar Distribution R2 ***************************/
int SavePairDistribtuionR2(void) {
  int i;
  int save_pure = ON;
  FILE *out = NULL;
  FILE *out2 = NULL;
  static int initialized = OFF;

  if(PD.times_measured == 0) {
    Warning("Attempt to save pair distribution function without any measurement\n");
    return 1;
  }

  out = fopen("outpdr2.dat", (initialized || file_append)?"a":"w");
  if(out == NULL) {
    perror("\nError:");
    Warning("can't save distribution function to file %s\n", "outpdr2.dat");
    return 1;
  }

  if(MC != DIFFUSION) { // VMC or Classical MC
    save_pure = OFF;
  }

  if(save_pure) {
    if(PD_pure.times_measured == 0) {
      Warning("Attempt to save pure pair distribution function without any measurement\n");
      save_pure = OFF;
    }
    out2 = fopen("outpdr2p.dat", (initialized || file_append)?"a":"w");
    if(out2 == NULL) {
      perror("\nError:");
      Warning("can't save distribution function to file %s\n", "outpdr2p");
      return 1;
    }
  }

  fprintf(out, "%.15" LE "\n", PD.r2 / (DOUBLE) PD.times_measured);
  PD.r2 = 0.;
  fclose(out);

  if(save_pure) {
    fprintf(out2, "%.15" LE "\n", PD_pure.r2 / (DOUBLE) PD_pure.times_measured);
    PD_pure.r2 = 0.;
    fclose(out2);
  }

  initialized = ON;
  return 0;
}

/*********************** Save Radial Distribution ****************************/
int SaveRadialDistribution(void) {
  int i;
  DOUBLE r;
  FILE *out, *out2 = NULL;
  int save_pure = ON;
  char file_name[40], file_name_pure[40];
  DOUBLE normalization, normalization2;
#ifdef BC_ABSENT
#ifdef TRIAL_3D
  DOUBLE norm,f;
#endif
#endif
  static int initialized = OFF;

#ifdef BC_1DPBC_Z // not measured
  return 1;
#endif

  if(RD.times_measured == 0) {
    Warning("Attempt to save radial distribution function without any measurement\n");
    return 1;
  }

#ifdef BC_ABSENT
  strcpy(file_name, file_RD);
#else
  strcpy(file_name, "outrdx.dat");
#endif

  out = fopen(file_name, (initialized || file_append)?"a":"w");
  if(out == NULL) {
    perror("\nError:");
    Warning("can't save radial distribution function to file %s\n", file_name);
    return 1;
  }

  if(MC != DIFFUSION) save_pure = OFF;

  if(save_pure) {
    if(RD_pure.times_measured == 0) {
      Warning("Attempt to save pure radial distribution function without any measurement\n");
      save_pure = OFF;
    }
#ifdef BC_ABSENT
    strcpy(file_name_pure, file_RD_pure);
#else
    strcpy(file_name_pure, "outrdxp.dat");
#endif

    out2 = fopen(file_name_pure, (initialized || file_append)?"a":"w");
    if(out2 == NULL) {
      perror("\nError:");
      Warning("can't save pure radial distribution function to file %s\n", file_name_pure);
      return 1;
    }
  }

#ifdef BC_ABSENT // Save distribution in a trap
#ifdef TRIAL_1D
  normalization = 1./(DOUBLE)(2.*RD.step*RD.times_measured*N);
  if(save_pure) normalization2 = 1./(DOUBLE)(2.*RD.step*RD_pure.times_measured*N);
#endif
#ifdef TRIAL_2D
  normalization = 1./(DOUBLE)(2*PI*RD.step*RD.times_measured*N);
  if(save_pure) normalization2 = 1./(DOUBLE)(2*PI*RD.step*RD_pure.times_measured*N);
#endif
#ifdef TRIAL_3D
  normalization = 1./(DOUBLE)(4.*PI*RD.step*RD.times_measured*N);
  if(save_pure) normalization2 = 1./(DOUBLE)(4.*PI*RD.step*RD_pure.times_measured*N);
#endif

  for(i=0; i<RD.size; i++) {
    r = RD.min + (DOUBLE) (i+0.5) * RD.step;
#ifdef TRIAL_1D
    fprintf(out, "%" LF " %.15" LE "\n", r, normalization*(DOUBLE)RD.N[i]);
    if(save_pure) fprintf(out2, "%" LF " %.15" LE "\n", r, normalization2*(DOUBLE)RD_pure.N[i]);
#endif
#ifdef TRIAL_2D
    fprintf(out, "%" LF " %.15" LE "\n", r, normalization*(DOUBLE)RD.N[i]/r);
    if(save_pure) fprintf(out2, "%" LF " %.15" LE "\n", r, normalization2*(DOUBLE)RD_pure.N[i]/r);
#endif
#ifdef TRIAL_3D
    fprintf(out, "%" LF " %.15" LE "\n", r, normalization*(DOUBLE)RD.N[i]/(r*r));
    if(save_pure) fprintf(out2, "%" LF " %.15" LE "\n", r, normalization2*(DOUBLE)RD_pure.N[i]/(r*r));
#endif
}

#else // Save distribution n(x) in a homogeneous system
  normalization = 1./((DOUBLE) (RD.times_measured)*N*RD.step/RD.max);
  if(save_pure) normalization2 =  1./((DOUBLE) (RD_pure.times_measured)*N*RD.step/RD.max);

  for(i=0; i<RD.size; i++) {
    r = RD.min + (DOUBLE) (i+0.5) * RD.step;
    fprintf(out, "%" LF " %.15" LE "\n", r, normalization*(DOUBLE)RD.N[i]);
    if(save_pure) fprintf(out2, "%" LF " %.15" LE "\n", r, normalization2*(DOUBLE)RD_pure.N[i]);
  }
#endif

  if(save_pure) {
    for(i=0; i<RD.size; i++) RD_pure.N[i] = 0;
    RD_pure.times_measured = 0;
  }

  fclose(out);
  if(save_pure) fclose(out2);

  initialized = ON;

  return 0;
}

// for anisotropic trap or tube save the z direction data
int SaveRadialZDistribution(void) {
  int i;
  DOUBLE r;
  FILE *out, *out2 = NULL;
  char file_name[40], file_name_pure[40];
  static int initialized = OFF;
  DOUBLE normalization;
  DOUBLE normalization2;
  int save_pure = ON;

#ifdef BC_ABSENT
    strcpy(file_name, file_RDz);
#else
    strcpy(file_name, "outrdz.dat");
#endif

#ifdef TRIAL_3D  // 3D trap
#ifdef BC_2DPBC  // quasi 2D
    strcpy(file_name, "outrdz.dat");
#endif
#endif

  out = fopen(file_name, (initialized || file_append)?"a":"w");
  if(out == NULL) {
    perror("\nError:");
    Warning("can't save distribution function to file %s\n", file_name);
    return 1;
  }

  if(MC != DIFFUSION) save_pure = OFF;

  if(save_pure) {
#ifdef BC_ABSENT
      strcpy(file_name_pure, file_RDz_pure);
#else
      strcpy(file_name_pure, "outrdzp.dat");
#endif

    out2 = fopen(file_name_pure, (initialized || file_append)?"a":"w");
    if(out2 == NULL) {
      perror("\nError:");
      Warning("can't save distribution function to file %s\n", file_name_pure);
      return 1;
    }
  }
  // Save distribution n(z)
#ifdef TRIAL_1D // 1D
#ifdef BC_ABSENT // trap

  normalization = 1./(2*RDz.step*(DOUBLE) (RDz.times_measured * N));
  for(i=0; i<RDz.size; i++) {
    r = (DOUBLE) (i+0.5) * RDz.step;
    fprintf(out, "%" LF " %.15" LE "\n", r, normalization*(DOUBLE)RDz.N[i]);
  }

  if(save_pure) {
    if(RDz_pure.times_measured == 0) {
      Warning("  Attempt to save pure RDz distribution without measurements\n");
    }
    else {
      //normalization = 1./(2*RDz.step*(DOUBLE) (RDz.times_measured * N));
      normalization2 =  1./ (2.*RDz.step * (DOUBLE) (RDz_pure.times_measured*N));
      for(i=0; i<RDz.size; i++) {
        r = (DOUBLE) (i+0.5) * RDz.step;
        //fprintf(out2, "%" LF " %.15" LE "\n", r, normalization*(DOUBLE)RDz.N[i]);
        fprintf(out2, "%" LF " %.15" LE "\n", r, normalization2*(DOUBLE)RDz_pure.N[i]);
      }
    }
  }
#else // 1D PBC
  normalization = 1./(RDz.step * (DOUBLE) (RDz.times_measured * N));
  if(save_pure) normalization2 =  1./ (RDz.step * (DOUBLE) (RDz_pure.times_measured*N));
  for(i=0; i<RDz.size; i++) {
    r = (DOUBLE) (i+0.5) * RDz.step;
#ifdef VEXT_COS2 // in lattice use lattice units
    fprintf(out, "%" LF " %.15" LE "\n", r, normalization*(DOUBLE)RDz.N[i]*N/n);
    if(save_pure) fprintf(out2, "%" LF " %.15" LE "\n", r, normalization2*(DOUBLE)RDz_pure.N[i]*N/n);
#else // otherwise scale length with density
    fprintf(out, "%" LF " %.15" LE "\n", r*n, normalization*(DOUBLE)RDz.N[i]*N/n);
    if(save_pure) fprintf(out2, "%" LF " %.15" LE "\n", r*n, normalization2*(DOUBLE)RDz_pure.N[i]*N/n);
#endif
  }
#endif
#endif

#ifdef TRIAL_2D // 2D
#ifdef BC_ABSENT // trap
  normalization = 1/Sqrt(2.*PI * RDz.step * RDz.width * RDz.width * (DOUBLE) (RDz.times_measured * N));
  for(i=0; i<RDz.size; i++) {
    r = (DOUBLE) (i+0.5) * RDz.step;
    fprintf(out, "%" LF " %.15" LE "\n", r, normalization*Sqrt((DOUBLE)RDz.N[i]));
  }
#else // 2D PBC
  normalization =  1./ (2.*RDz.step * (DOUBLE) (RDz.times_measured*N));

  if(save_pure) normalization2 =  1./ (2.*RDz.step * (DOUBLE) (RDz_pure.times_measured*N));
  for(i=0; i<RDz.size; i++) {
    r = (DOUBLE) (i+0.5) * RDz.step;
    fprintf(out, "%" LF " %.15" LE "\n", r, normalization*(DOUBLE)RDz.N[i]);
    if(save_pure) fprintf(out2, "%" LF " %.15" LE "\n", r, normalization2*(DOUBLE)RDz_pure.N[i]);
  }
#endif
#endif

#ifdef TRIAL_3D // 3D trap
#ifdef BC_ABSENT 
  // save square root of the density
  //normalization = 1/Sqrt(2.*PI * RDz.step * RDz.width * RDz.width * (DOUBLE) (RDz.times_measured * N));
  normalization = 1./(RDz.step*(DOUBLE) (RDz.times_measured * N));
  if(save_pure) normalization2 = 1./(RDz.step*(DOUBLE) (RDz_pure.times_measured * N));
#else // homogeneous
  normalization = 1./((DOUBLE) (RDz.times_measured)*N*RDz.step/RDz.max);
  if(save_pure) normalization2 =  1./((DOUBLE) (RDz_pure.times_measured)*N*RDz.step/RDz.max);
#endif

  for(i=0; i<RDz.size; i++) {
    r = ((DOUBLE) (i)+0.5) * RDz.step;
    fprintf(out, "%" LF " %.15" LE "\n", r, normalization*(DOUBLE)RDz.N[i]);
    if(save_pure) fprintf(out2, "%" LF " %.15" LE "\n", r, normalization2*(DOUBLE)RDz_pure.N[i]);
  }
#endif

  if(save_pure) {
    for(i=0; i<RDz.size; i++) RDz_pure.N[i] = 0;
    RDz_pure.times_measured = 0;
  }

  fclose(out);
  if(save_pure) fclose(out2);

  initialized = ON;

  return 0;
}

int LoadRadialZDistribution(void) {
  int i;
  DOUBLE r, dummy;
  FILE *in;

  /* for anisotropic trap or tube save the z direction data */
  in = fopen("rdz.dat", "r");
  if (in == NULL) {
    perror("\nError:");
    Warning("can't load radial z distribution function from file rdz.dat\n");
    return 1;
  }

  for(i=0; i<RDz.size; i++) {
    fscanf(in, "%" LF " %" LF " %" LF  "\n", &r, &RDz.f[i], &dummy);
  }
  if(fabs(RDz.max-r) >RDz.step) {
    Warning("  Loadind radial z distribution: input file has different z scale\n");
    Warning("  OBDM will be measured incorrectly !\n");
  }
  fclose(in);

  return 0;
}

/****************************** Save OBDM ************************************/
int SaveOBDM(void) {
  int i;
  DOUBLE r;
  FILE *out;
#ifdef BC_ABSENT // i.e. in the trap
  FILE *out2;
#endif
#ifdef OBDM_FERMIONS
  FILE *out3;
#endif
  static int initialized = OFF;

  out = fopen(file_OBDM, (initialized || file_append)?"a":"w");
  if(out == NULL) {
    perror("\nError:");
    Warning("can't save one body density matrix to file %s\n", file_OBDM);
    return 1;
  }

#ifdef BC_ABSENT // i.e. in the trap
  out2 = fopen("outdrn.dat", (initialized || file_append)?"a":"w");
  if(out2 == NULL) {
    perror("\nError:");
    Warning("can't save one body density matrix to file outdrn.dat\n");
    return 1;
  }
#endif

#ifdef OBDM_FERMIONS
  out3 = fopen("outdrf.dat", (initialized || file_append)?"a":"w");
  if(out3 == NULL) {
    perror("\nError:");
    Warning("can't save one body density matrix to file outdrf.dat\n");
    return 1;
  }
#endif
  if(OBDM.f[0]/OBDM.N[0]<0.1) {
     OBDM.f[0] = 0.;
     Warning("OBDM goes to zero as r -> 0 !\n");
  }

  for(i=0; i<OBDM.size; i++) {
    r = (DOUBLE) (i+0.5) * OBDM.step;

//#ifdef BC_1DPBC_ANY
//    r *= n; // rescale by density
//#endif

    if(OBDM.f[i]/OBDM.N[i]>100) Warning("Strange (large) value found in OBDM!\n");

    if(OBDM.N[i] > 0)
      fprintf(out, "%" LF " %.15" LE "\n", r, OBDM.f[i]/OBDM.N[i]);
    else
      fprintf(out, "%" LF " 0.0\n", r);

#ifdef BC_ABSENT // i.e. in the trap
    if(OBDMtrap.N[i] > 0) 
      fprintf(out2, "%" LF " %.15" LE "\n", r, OBDMtrap.f[i]/OBDMtrap.N[i]);
    else 
      fprintf(out2, "%" LF " 0.0\n", r);
#endif

#ifdef OBDM_FERMIONS
    if(OBDM.N[i] > 0) 
      fprintf(out3, "%" LF " %.15" LE "\n", r, OBDMfermi.f[i]/OBDM.N[i]);
    else
      fprintf(out3, "%" LF " 0.0\n", r);
#endif
  }

  fclose(out);
#ifdef BC_ABSENT // i.e. in the trap
  fclose(out2);
#endif
#ifdef OBDM_FERMIONS
  fclose(out3);
#endif

#ifdef SPECKLES
  out = fopen("outcf.dat", (initialized || file_append)?"a":"w");
  if(OBDM.CFN>0)
    fprintf(out, " %.15" LE "\n", OBDM.CF/(DOUBLE)OBDM.CFN);
  else
    fprintf(out,"0\n");
  fclose(out);
#endif

  initialized = ON;

  return 0;
}

/****************************** Save Momentum distribution *****************/
int SaveMomentumDistribution(void) {
#ifdef TRIAL_1D
  FILE *out;
  static int initialized = OFF;
  DOUBLE normalization;

  int k;
  // Save momentum distribution
  // n(k) = 2 \int_0^\infty Cos(k x) g_1(x) dx
  out = fopen("outnk.dat", (initialized || file_append)?"a":"w");
  if(out == NULL) {
    perror("\nError:");
    Warning("can't save distribution function to file outnk.dat\n");
    return 1;
  }

  normalization = 2.*OBDM.max/(DOUBLE)OBDM.Nktimes_measured; // 2 \int_0^\inf n(k) dk/(2pi) = 1
  for(k=0; k<OBDM.Nksize; k++) {
#ifdef BC_ABSENT
    fprintf(out, "%" LF " %.15" LE "\n", OBDM.k[k], normalization*OBDM.Nk[k]);
#else // save in units of the density
    fprintf(out, "%" LF " %.15" LE "\n", OBDM.k[k]/n, OBDM.Nk[k]*(DOUBLE)N/(DOUBLE)OBDM.Nktimes_measured);
#endif
    OBDM.Nk[k] = 0.;
  }
  fclose(out);

#ifdef OBDM_FERMIONS
  out = fopen("outnkf.dat", (initialized || file_append)?"a":"w");
  if(out == NULL) {
    perror("\nError:");
    Warning("can't save one body density matrix to file outnkf.dat\n");
    return 1;
  }

  for(k=0; k<OBDM.Nksize; k++) {
    //nk = 0.;
    //for(i=0; i<OBDM.size; i++) {
    //  r = (DOUBLE) (i+0.5) * OBDM.step;
    //  if(OBDM.N[i] > 0) nk += Cos(OBDM.k[k]*r)*OBDMfermi.f[i]/OBDM.N[i];
    //}
    // save in units of the density
    //fprintf(out3, "%" LF " %.15" LE "\n", OBDM.k[k]/n, nk*2*OBDM.step*n);
#ifdef BC_ABSENT
    fprintf(out, "%" LF " %.15" LE "\n", OBDM.k[k], normalization*OBDMfermi.Nk[k]);
#else // save in units of the density
    fprintf(out, "%" LF " %.15" LE "\n", OBDM.k[k]/n, OBDMfermi.Nk[k]*(DOUBLE) N/(DOUBLE)OBDM.Nktimes_measured);
#endif
    OBDMfermi.Nk[k] = 0.;
  }
  fclose(out);
  OBDM.Nktimes_measured = 0;
#endif

  initialized = ON;
#endif

  return 0;
}

/****************************** Save OBDM pure *******************************/
int SaveOBDMpure(void) {
#ifdef VIRTUAL_WALKERS // can be called also instead of Save TBDM pure
  int i,w;
  DOUBLE t, weight;
  static int initialized_pure = OFF;
  FILE *out, *out2;
  DOUBLE drw, drp;
  int Nvirtual;
#ifdef OBDM_FERMIONS
  FILE *out3;
  DOUBLE drpf;
#endif

  out = fopen("outdrw.dat", (initialized_pure || file_append)?"a":"w");
  if(out == NULL) {
    perror("\nError:");
    Warning("can't save one body density matrix to file %s\n", "outdrw.dat");
    return 1;
  }
  out2 = fopen("outdrp.dat", (initialized_pure || file_append)?"a":"w");
  if(out2 == NULL) {
    perror("\nError:");
    Warning("can't save one body density matrix to file %s\n", "outdrp.dat");
    return 1;
  }
#ifdef OBDM_FERMIONS
  out3 = fopen("outdrfp.dat", (initialized_pure || file_append)?"a":"w");
  if(out3 == NULL) {
    perror("\nError:");
    Warning("can't save one body density matrix to file %s\n", "outdrfp.dat");
    return 1;
  }
#endif

  for(i=0; i<W[Nwalkers].OBDM_position; i++) { // save data at imaginary time "t"
    t = (i+1)*dt*(DOUBLE)Nmeasure;
    drp = drw = 0.;
#ifdef OBDM_FERMIONS
    drpf = 0.;
#endif
    Nvirtual = 0;
    for(w=Nwalkers; w<Nwalkers+Nwalkers_pure; w++) { // average over all virtual walkers
      if(W[w].status != VIRTUAL) {
        Warning("  SaveOBDMpure: skipping some not pure walkers\n");
      }
      else {
        Nvirtual++;
        weight = Exp(W[w].OBDMweight[i]); // weight is logarithmic
        //weight = W[w].OBDMweight[i]; // weight is linear
        drw += weight;
        drp += W[w].OBDMpureB*weight;
#ifdef OBDM_FERMIONS
        drpf += W[w].OBDMpureF*weight;
#endif
      }
    }
    fprintf(out,"%" LE "  %.15" LE "\n", t, drw/(DOUBLE)Nvirtual);
    fprintf(out2,"%" LE "  %.15" LE "\n", t, drp/(DOUBLE)Nvirtual);
#ifdef OBDM_FERMIONS
    fprintf(out3,"%" LE " x %.15" LE "\n", t, drpf/(DOUBLE)Nvirtual);
#endif
  }
  for(w=Nwalkers; w<Nwalkers+Nwalkers_pure; w++) W[w].status = DEAD;

  WalkersSort();

  fclose(out);
  fclose(out2);
#ifdef OBDM_FERMIONS
  fclose(out3);
#endif
  initialized_pure = ON;
#endif
  return 0;
}

/****************************** Save Momentum Distribution Pure **************/
// Save PURE momentum distribution
// n(k) = 2 \int_0^\infty Cos(k x) g_1(x) dx
int SaveMomentumDistributionPure(void) {
#ifdef MEASURE_PURE_OBDM
  int i,w,k;
  DOUBLE t, weight;
  static int initialized = OFF;
  FILE *out, *out2, *out3;
  DOUBLE nkB, drw;
  int Nvirtual;
#ifdef OBDM_FERMIONS
  DOUBLE nkF;
  FILE *out4;
#endif

#ifndef TRIAL_1D
  Warning(" Save momentum distribtion is not implemented!\n");
  return 1;
#endif

  out = fopen("outdrw.dat", (initialized || file_append)?"a":"w");
  if(out == NULL) {
    perror("\nError:");
    Warning("can't save pure momentum distribution to file %s\n", "outdrw.dat");
    return 1;
  }

  out3 = fopen("outnkp.dat", (initialized || file_append)?"a":"w");
  if(out3 == NULL) {
    perror("\nError:");
    Warning("can't save distribution function to file outnkp.dat\n");
    return 1;
  }

#ifdef OBDM_FERMIONS
  out4 = fopen("outnkfp.dat", (initialized || file_append)?"a":"w");
  if(out4 == NULL) {
    perror("\nError:");
    Warning("can't save pure momentum distribution to file %s\n", "outdrfp.dat");
    return 1;
  }
#endif

  for(i=0; i<W[Nwalkers].OBDM_position; i++) {
    t = (i+1)*dt*(DOUBLE)Nmeasure;
    drw = 0.;
    Nvirtual = 0;
    for(w=Nwalkers; w<Nwalkers+Nwalkers_pure; w++) {
      if(W[w].status != VIRTUAL) {
        Warning("  SaveMomentumDistributionPure: skipping some not pure walkers\n");
      }
      else {
        Nvirtual++;
        weight = Exp(W[w].OBDMweight[i]);
        drw += weight;
      }
    }
    fprintf(out,"%" LE "  %.15" LE "\n", t, drw/(DOUBLE)Nvirtual);
  }

  for(i=0; i<W[Nwalkers].OBDM_position; i++) { // time "t"
    t = (i+1)*dt*(DOUBLE)Nmeasure;
    for(k=0; k<OBDM.Nksize; k++) { // loop over k
      Nvirtual = 0;
      nkB = 0.;
#ifdef OBDM_FERMIONS
      nkF = 0.;
#endif
      for(w=Nwalkers; w<Nwalkers+Nwalkers_pure; w++) { // average over virtual walkers
        Nvirtual++;
        weight = Exp(W[w].OBDMweight[i]);
        //drp += W[w].OBDMpureB*weight;
        nkB += Cos(OBDM.k[k]*W[w].OBDMpure_r)*weight*W[w].OBDMpureB;
#ifdef OBDM_FERMIONS
        nkF += Cos(OBDM.k[k]*W[w].OBDMpure_r)*weight*W[w].OBDMpureF;
#endif
      }
      fprintf(out3,"%.2" LE " ", nkB/(DOUBLE)Nvirtual);
#ifdef OBDM_FERMIONS
      fprintf(out4,"%.2" LE " ", nkF/(DOUBLE)Nvirtual);
#endif
    }
    fprintf(out3,"\n");
#ifdef OBDM_FERMIONS
    fprintf(out4,"\n");
#endif
  }

  out2 = fopen("outnkgridt.dat", "w"); // save time vector
  if(out2 == NULL) {
    perror("\nError:");
    Warning("can't save pure momentum distribution to file %s\n", "outnkgridt.dat");
    return 1;
  }

  for(i=0; i<W[Nwalkers].OBDM_position; i++) fprintf(out2,"%" LE "\n", (i+1)*dt*(DOUBLE)Nmeasure);
  fclose(out2);

  out2 = fopen("outnkgridk.dat", "w"); // save k vector
  if(out2 == NULL) {
    perror("\nError:");
    Warning("can't save pure momentum distribution to file %s\n", "outnkgridk.dat");
    return 1;
  }

  for(k=0; k<OBDM.Nksize; k++) fprintf(out2,"%" LE "\n", OBDM.k[k]);

  fclose(out2);

  for(w=Nwalkers; w<Nwalkers+Nwalkers_pure; w++) W[w].status = KILLED;

  fclose(out);
  fclose(out3);
#ifdef OBDM_FERMIONS
  fclose(out4);
#endif
  initialized = ON;
#endif
  return 0;
}

/****************************** Save TBDM ************************************/
int SaveTBDM(void) {
  int i;
  DOUBLE r;
  FILE *out;
  static int initialized = OFF;

  out = fopen("outtr.dat", (initialized || file_append)?"a":"w");
  if (out == NULL) {
    perror("\nError:");
    Warning("can't save two body density matrix to file %s\n", "outtr.dat");
    return 1;
  }

  if(TBDM.f[0]/TBDM.N[0]<0.1) {
     TBDM.f[0] = 0.;
     Warning("TBDM goes to zero as r -> 0 !\n");
  }

  for(i=0; i<TBDM.size; i++) {
    r = (DOUBLE) (i+0.5) * TBDM.step;
    if(TBDM.N[i] > 0) 
      fprintf(out, "%" LF " %.15" LE "\n", r, TBDM.f[i]/TBDM.N[i]);
    else
      fprintf(out, "%" LF " 0.0\n", r);
  }
  fclose(out);

#ifdef OBDM_FERMIONS
  out = fopen("outtrf.dat", (initialized || file_append)?"a":"w");
  if (out == NULL) {
    perror("\nError:");
    Warning("can't save one body density matrix to file %s\n", "outtrf.dat");
    return 1;
  }

  if(TBDMfermi.f[0]/TBDM.N[0]<0.1) {
     TBDMfermi.f[0] = 0.;
     Warning("TBDM (fermions) goes to zero as r -> 0 !\n");
  }

  for(i=0; i<TBDM.size; i++) {
    r = (DOUBLE) (i+0.5) * TBDM.step;
    if(TBDM.N[i] > 0) 
      fprintf(out, "%" LF " %.15" LE "\n", r, TBDMfermi.f[i]/TBDM.N[i]);
    else
      fprintf(out, "%" LF " 0.0\n", r);
  }
  fclose(out);
#endif

  if(measure_OBDM_MATRIX == OFF) SaveTBDMMatrix();

  initialized = ON;

  return 0;
}

/****************************** Save OBDM Matrix *****************************/
int SaveOBDMMatrix(void) {
  int i, j;
  FILE *out, *out1;
  static int initialized = OFF;
  DOUBLE normalization;

  out = fopen(file_OBDM_MATRIX, (initialized || file_append)?"a":"w");
  if(out == NULL) {
    perror("\nError:");
    Warning("can't save one body density matrix to file %s\n", file_OBDM_MATRIX);
    return 1;
  }

  normalization = 1./(DOUBLE) OBDM_MATRIX.times_measured;
  for(i=0; i<OBDM_MATRIX.size; i++) {
    for(j=0; j<OBDM_MATRIX.size; j++) {
      if(OBDM_MATRIX.size>100) // save space
        fprintf(out, "%.3" LE " ", normalization*OBDM_MATRIX.f[i][j]);
      else
        fprintf(out, "%.15" LE " ", normalization*OBDM_MATRIX.f[i][j]);
    }
    fprintf(out, "\n");
  }
  fclose(out);

  for(i=0; i<OBDM_MATRIX.size; i++) {
    for(j=0; j<OBDM_MATRIX.size; j++) {
      OBDM_MATRIX.f[i][j] = 0;
      OBDM_MATRIX.N[i][j] = 0;
    }
  }

#ifdef OBDM_FERMIONS
  out = fopen("outobdmf.dat", (initialized || file_append)?"a":"w");
  for(i=0; i<OBDM_MATRIX.size; i++) {
    for(j=0; j<OBDM_MATRIX.size; j++) {
      if(OBDM_MATRIX.size>100) // save space
        fprintf(out, "%.3" LE " ", normalization*OBDMfermi_MATRIX.f[i][j]);
      else
        fprintf(out, "%.15" LE " ", normalization*OBDMfermi_MATRIX.f[i][j]);
      OBDMfermi_MATRIX.f[i][j] = 0;
    }
    fprintf(out, "\n");
  }
  fclose(out);
#endif

  OBDM_MATRIX.times_measured = 0;

  out1 = fopen("outobdmh.dat", "w");
  fprintf(out1, "%" LG "\n", OBDM_MATRIX.step);
  fprintf(out1, "%i\n", OBDM_MATRIX.size);
  fprintf(out1, "%i\n", N);
  fprintf(out1, "%" LG "\n", D);
  fclose(out1);

  initialized = ON;

  return 0;
}

/****************************** Save Momentum Distribution Matrix *****************************/
int SaveMomentumDistributionMatrix(void) {
  int i, j;
  FILE *out, *out1;
  static int initialized = OFF;
  DOUBLE normalization;
  
  out = fopen("outnkxy.dat", (initialized || file_append)?"a":"w");
  if(out == NULL) {
    perror("\nError:");
    Warning("can't save one body density matrix to file %s\n", "outnkxyb.dat");
    return 1;
  }

  normalization = 1./(DOUBLE) OBDM_MATRIX.times_measured;
  for(i=0; i<gridSk; i++) {
    for(j=0; j<gridSk; j++) {
      fprintf(out, "%.15" LE " ", normalization*Nk_MATRIX.f[i][j]);
      Nk_MATRIX.f[i][j] = 0;
    }
    fprintf(out, "\n");
  }
  fclose(out);
  
#ifdef OBDM_FERMIONS
  out = fopen("outnkxyf.dat", (initialized || file_append)?"a":"w");
  for(i=0; i<gridSk; i++) {
    for(j=0; j<gridSk; j++) {
      fprintf(out, "%.15" LE " ", normalization*Nkfermi_MATRIX.f[i][j]);
      Nkfermi_MATRIX.f[i][j] = 0;
    }
    fprintf(out, "\n");
  }
  fclose(out);
#endif

  out1 = fopen("outnkxyh.dat", "w");
  fprintf(out1, "%" LG "\n", Nk_MATRIX.step);
  fprintf(out1, "%i", gridSk);
  fclose(out1);

  initialized = ON;

  return 0;
}

/****************************** Save OBDM Matrix *****************************/
int SaveTBDMMatrix(void) {
  int i, j;
  FILE *out, *out1;
  static int initialized = OFF;

  out = fopen("outtbdm.dat", (initialized || file_append)?"a":"w");
  if (out == NULL) {
    perror("\nError:");
    Warning("can't save one body density matrix to file %s\n", "outtbdm.dat");
    return 1;
  }

  for(i=0; i<TBDM_MATRIX.size; i++) {
    for(j=0; j<TBDM_MATRIX.size; j++) {
      if(TBDM_MATRIX.N[i][j] > 0)
        fprintf(out, "%.15" LE " ", TBDM_MATRIX.f[i][j]/(DOUBLE) TBDM_MATRIX.times_measured/(2*TBDM_MATRIX.step*2*TBDM_MATRIX.step*(1-1./(DOUBLE)N)));
      else
        fprintf(out, "0 ");
      TBDM_MATRIX.f[i][j] = 0;
      TBDM_MATRIX.N[i][j] = 0;
    }
    fprintf(out, "\n");
  }
  fclose(out);

  out1 = fopen("outtbdmh.dat", "w");
  fprintf(out1, "%" LG "\n", TBDM_MATRIX.step);
  fprintf(out1, "%i", TBDM_MATRIX.size);
  fclose(out1);

#ifdef OBDM_FERMIONS
  out = fopen("outtbdmf.dat", (initialized || file_append)?"a":"w");
  if (out == NULL) {
    perror("\nError:");
    Warning("can't save one body density matrix to file %s\n", "outtbdm.dat");
    return 1;
  }

  for(i=0; i<TBDM_MATRIX.size; i++) {
    for(j=0; j<TBDM_MATRIX.size; j++) {
      if(TBDMfermi_MATRIX.N[i][j] > 0)
        fprintf(out, "%.15" LE " ", TBDMfermi_MATRIX.f[i][j]/(DOUBLE) TBDM_MATRIX.times_measured/(2*TBDM_MATRIX.step*2*TBDM_MATRIX.step*(1-1./(DOUBLE)N)));
      else
        fprintf(out, "0 ");
      TBDMfermi_MATRIX.f[i][j] = 0;
      TBDMfermi_MATRIX.N[i][j] = 0;
    }
    fprintf(out, "\n");
  }
  fclose(out);
#endif

  TBDM_MATRIX.times_measured = 0;
  initialized = ON;

  return 0;
}

/****************************** Save Pair Distribution Matrix *****************************/
int SavePairDistribtuionMatrix(void) {
  int i, j;
  DOUBLE x, y;
  FILE *out, *out1;
  static int initialized = OFF;
  DOUBLE normalization;

  out = fopen("outpdm.dat", (initialized || file_append)?"a":"w");
  out1 = fopen("outpdh.dat", "w");
  if (out == NULL) {
    perror("\nError:");
    Warning("can't save pair distribution matrix to file %s\n", "outpdm.dat");
    return 1;
  }

  normalization = 1./ ((DOUBLE)(PD_MATRIX.times_measured)*PD_MATRIX.step_x*PD_MATRIX.step_y);
  for(i=0; i<gridPD_MATRIX_x; i++) {
    x = (DOUBLE) (i+0.5) * PD_MATRIX.step_x;
    for(j=0; j<gridPD_MATRIX_y; j++) {
      y = (DOUBLE) (j+0.5) * PD_MATRIX.step_y;

      if(PD_MATRIX.N[i][j] > 0)
        fprintf(out, "%.2" LE " ", normalization*PD_MATRIX.N[i][j]);
      else
        fprintf(out, "0 ");
      PD_MATRIX.N[i][j] = 0;
    }
    fprintf(out, "\n");
  }
  PD_MATRIX.times_measured = 0;
  fclose(out);
  fprintf(out1, "%i\n", gridPD_MATRIX_x);
  fprintf(out1, "%i\n", gridPD_MATRIX_y);
  fprintf(out1, "%" LG "\n", PD_MATRIX.step_x);
  fprintf(out1, "%" LG "\n", PD_MATRIX.step_y);
  fclose(out1);

  initialized = ON;

  return 0;
}

/****************************** Save R2 **************************************/
// R2 - <R^2>
// R4 - <R^4>
int SaveR2(DOUBLE f1, DOUBLE f2) {
  FILE *out;
  static int initialized = OFF;

  out = fopen(file_R2, (initialized || file_append)?"a":"w");
  if(out == NULL) {
    perror("\nError:");
    Warning("can't save R2 to file %s\n", file_R2);
    return 1;
  }

  fprintf(out, "%" LG " %" LG "\n", f1, f2);

  fclose(out);
  initialized = ON;

  return 0;
}

/****************************** Save Z2 **************************************/
int SaveZ2(DOUBLE f1, DOUBLE f2) {
  FILE *out;
  static int initialized = OFF;

  out = fopen(file_z2, (initialized || file_append)?"a":"w");
  if(out == NULL) {
    perror("\nError:");
    Warning("can't save Z2 to file %s\n", file_z2);
    return 1;
  }

  fprintf(out, "%" LG " %" LG "\n", f1, f2);

  fclose(out);
  initialized = ON;

  return 0;
}

/****************************** Save Pure R2 *********************************/
int SavePureR2(DOUBLE f) {
  FILE *out;
  static int initialized = OFF;

  out = fopen(file_R2_pure, (initialized || file_append)?"a":"w");
  if(out == NULL) {
    perror("\nError:");
    Warning("can't save R2 to file %s\n", file_R2_pure);
    return 1;
  }

  if(initialized) fprintf(out, "%" LG "\n", f);

  fclose(out);
  initialized = ON;

  return 0;
}

/****************************** Save Pure Z2 *********************************/
int SavePureZ2(DOUBLE f) {
  FILE *out;
  static int initialized = OFF;

  out = fopen(file_z2_pure, (initialized || file_append)?"a":"w");
  if(out == NULL) {
    perror("\nError:");
    Warning("can't save R2 to file %s\n", file_z2_pure);
    return 1;
  }

  if(initialized) fprintf(out, "%" LG "\n", f);

  fclose(out);
  initialized = ON;

  return 0;
}

/**************************** Save Mean R2 ************************************/
void SaveMeanR2(void) {
  xR2 += RD.r2recent;
  SaveR2(RD.r2recent/(DOUBLE)RD.times_measured, RD.r4/(DOUBLE)RD.times_measured);

  zR2 += RDz.r2recent;
  SaveZ2(RDz.r2recent/(DOUBLE)RDz.times_measured, RDz.r4/(DOUBLE)RDz.times_measured);

  RD.r2recent = RD.r4 = 0.;
  RDz.r2recent = RDz.r4 = 0.;
  //RD.times_measured = RDz.times_measured = 0;

  R2times_measured++;
}

/**************************** Save Mean R2 DMC ********************************/
void SaveMeanR2DMC(void) {
  DOUBLE R2=0., Z2=0.;
  DOUBLE R4=0., Z4=0.;
  DOUBLE R2pure=0., Z2pure=0.;
  DOUBLE weight, norm = 0.;
  int w;

  for(w=0; w<Nwalkers; w++) {
    weight = W[w].weight;
    norm += weight;
    R2 += W[w].r2 * weight;
    Z2 += W[w].z2 * weight;
    R4 += W[w].r4 * weight;
    Z4 += W[w].z4 * weight;
    //if(boundary == NO_BOUNDARY_CONDITIONS) {
    R2pure += W[w].r2old * weight;
    Z2pure += W[w].z2old * weight;
    //}
  }

  norm = 1./ (DOUBLE) Nwalkers;

  R2 *= norm;
  xR2 += R2;
  //R2 = Sqrt(R2);
  R4 *= norm;
  //R4 = Sqrt(R4);
  SaveR2(R2, R4);
  RD.r2recent = R2;

  R2pure = norm*R2pure;
  if(W[0].RD_wait == OFF) SavePureR2(R2pure);

  //if(boundary == NO_BOUNDARY_CONDITIONS) {
  Z2 *= norm;
  zR2 += Z2;
  //Z2 = Sqrt(Z2);
  Z4 *= norm;
  //Z4 = Sqrt(Z4);
  SaveZ2(Z2, Z4);
  RDz.r2recent = Z2;

  Z2pure = norm*Z2pure;
  if(W[0].RD_wait == OFF) SavePureZ2(Z2pure);
  //}
}

/****************************** Save RD Matrix *****************************/
int SaveRadialDistributionMatrix(void) {
  int i, j;
  //DOUBLE x, y;
  FILE *out, *out1;
  static int initialized = OFF;
  DOUBLE normalization;

  out = fopen("outrdm.dat", (initialized || file_append)?"a":"w");
  out1 = fopen("outrdmh.dat", "w");
  if(out == NULL) {
    perror("\nError:");
    Warning("can't save pair distribution matrix to file %s\n", "outrdm.dat");
    return 1;
  }

  //normalization = (DOUBLE) (N)/((DOUBLE)(RD_MATRIX.times_measured)*RD_MATRIX.step_x*RD_MATRIX.step_y);
  normalization = 1./((DOUBLE)(RD_MATRIX.times_measured)*RD_MATRIX.step_x*RD_MATRIX.step_y);
  for(i=0; i<gridRDx; i++) {
    //x = ((DOUBLE) i + 0.5)* RD_MATRIX.step_x;
    for(j=0; j<gridRDy; j++) {
      //y = ((DOUBLE) j + 0.5) * RD_MATRIX.step_y;

      if(RD_MATRIX.N[i][j] > 0)
        fprintf(out, "%.15" LE " ", normalization*RD_MATRIX.N[i][j]);
      else
        fprintf(out, "0 ");
      RD_MATRIX.N[i][j] = 0;
    }
    fprintf(out, "\n");
  }
  fprintf(out, "\n");
  fclose(out);
  fprintf(out1, "%i %i %" LG "  %" LG " ", gridRDx, gridRDy, RD_MATRIX.step_x, RD_MATRIX.step_y);
  fclose(out1);

  RD_MATRIX.times_measured = 0;

  initialized = ON;

  return 0;
}

/****************************** SaveWaveFunction *****************************/
void SaveWaveFunction(struct Grid *G, char *file_name, DOUBLE min, DOUBLE max) {
  int i;
  FILE *out;
  char wf_dat_spinfull[100];
  DOUBLE dx, x;
  DOUBLE f, lnf, fp, Eloc, E;
  DOUBLE testNaN;
  int N = 1000;
  int spin1, spin2;

  if(grid_trial>3) N = grid_trial;

  if(max<min) {
    Warning("Save Wavefunction: 'min' argument is larger than the 'max' one\n");
    max = 10.*min;
  }

#ifdef HARD_SPHERE
  if(min<a) Warning("Save wavefunction: 'min' is smaller than the size of the HS\n");
  min = a;
#endif

  //min = 1.;
  //max = min + G->step * 10;

  dx = (max-min) / (N-1);
  dx *= 1.2; // go above L/2

#ifdef INTERPOLATE_TYPE_SPLINE
  min = G->min;
#endif

#ifdef HARD_SPHERE
    min = a;
#endif

  for(spin1=0; spin1<Nspin; spin1++) {
    for(spin2=0; spin2<Nspin; spin2++) {

      Message("Saving wavefunction (%i - %i) components\n", spin1+1, spin2+1);

#ifdef SPINFULL
      sprintf(wf_dat_spinfull, "wf%i%i.dat", spin1+1, spin2+1);
#else
      strcpy(wf_dat_spinfull, file_name);
#endif
      out = fopen(wf_dat_spinfull, "w");

      fprintf(out, "#N= %i\n", N);
      fprintf(out, "#1x 2f 3lnf 4fp 5Eloc 6E\n");
      fprintf(out, "#1)x 2)f 3)lnf 4)fp 5)-f''/f - 2/r f'/f 6)-f''/f - (d-1)/r f'/f + (f'/f)^2 7)-f''/f - (d-1)/r f'/f + V 8) V\n"); // Eloc = -f''/f - (d-1)/r f'/f + (f'/f)^2
      fprintf(out, "# 1) x                 2) f                   3) lnf                  4) fp                   5) -f''/f - (d-1)/r f'/f    6) -f''/f - (d-1)/r f'/f + (f'/f)^2 7)-f''/f - (d-1)/r f'/f + V 8) V\n");

      for(i=0; i<N; i++) {
        //fprintf(out, "%.15" LE " %.15" LE " %.15" LE " %.15" LE " %.15" LE " %.15" LE "\n", x, Exp(InterpolateU(&G, x)), InterpolateU(&G, x),InterpolateFp(&G, x), InterpolateE(&G, x)-InterpolateFp(&G, x)*InterpolateFp(&G, x)+InteractionEnergy(x), InterpolateE(&G, x));
        x = min+dx*(DOUBLE) i;
        testNaN = exp(InterpolateU(G, x, spin1, spin2));
        f = exp(InterpolateU(G, x, spin1, spin2));
#ifdef OBDM_EXCITED_BOSONS
        f *= InterpolateExactSign(x);
#endif
        if(testNaN == testNaN && testNaN != 2*testNaN) { // check for NaN by comparison to itself and INF
          fprintf(out, "%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n", x, //1
            f, //2
            InterpolateU(G, x, spin1, spin2), //3
            InterpolateFp(G, x, spin1, spin2), //4
            InterpolateE(G, x, spin1, spin2) - InterpolateFp(G, x, spin1, spin2)*InterpolateFp(G, x, spin1, spin2), //5
            InterpolateE(G, x, spin1, spin2), //6
            InterpolateE(G, x, spin1, spin2) - InterpolateFp(G, x, spin1, spin2)*InterpolateFp(G, x, spin1, spin2) + InteractionEnergy_ij(x, spin1, spin2), //7
            energy_unit*InteractionEnergy_ij(x, spin1, spin2)); //8
        }
      }
     fclose(out);
    }
  }

/*#ifdef INTERPOLATE_SPLINE_ONE_BODY_Z_WF
  out = fopen(wf_dat, "w"); // harmonic oscillator
  fprintf(out, "#N= %i\n", N);
  fprintf(out, "#x f lnf fp Eloc E\n");
  for(i=0; i<N; i++) {
    x = min+dx * (DOUBLE) i;

    lnf = -alpha_z*x*x;
    fp = -2.*alpha_z*x;
    Eloc = alpha_z;

    lnf = InterpolateSplineU(&G1, x);
    fp = InterpolateSplineFp(&G1, x);
    Eloc = InterpolateSplineE(&G1, x);

    f = Exp(lnf);
    E = Eloc;
    
    fprintf(out, "%.15" LE " %.15" LE " %.15" LE " %.15" LE " %.15" LE " %.15" LE "\n", x, f, lnf, fp, Eloc, E);
  }
  fclose(out);
#endif*/
}

/****************************** SaveWaveFunction *****************************/
void SaveOneBodyWaveFunction(void) {
#ifdef ONE_BODY_TRIAL_TERMS
   int i;
   FILE* out;
   DOUBLE dx, x, f, testNaN, max, min;
   int Ngrid = 1000;

   out = fopen("wf1.dat", "w");

   Message("  Saving one-body wavefunction.\n");
   max = L;
   min = 0.;

   dx = (max - min) / (Ngrid - 1);

   for(i=0; i<Ngrid; i++) {
      x = min + dx * (DOUBLE)i;
#ifdef MOVE_IN_X
      testNaN = exp(OneBodyU(x, 0, 0, i, 0));
#endif
#ifdef MOVE_IN_Y
      testNaN = exp(OneBodyU(0, x, 0, i, 0));
#endif
#ifdef MOVE_IN_Z
      testNaN = exp(OneBodyU(0, 0, x, i, 0));
#endif
      //if(testNaN == testNaN && testNaN != 2 * testNaN) 
      fprintf(out, "%.15e %.15e\n", x, testNaN);
   }
   fclose(out);
#endif
}

/****************************** SaveWaveFunction *****************************/
void SaveThreeBodyWaveFunction(void) {
#ifdef THREE_BODY_TERMS
  int i;
  FILE *out;
  DOUBLE dx, x;
  DOUBLE f, lnf, fp, Eloc, E;
  DOUBLE testNaN;
  DOUBLE max,min;
  int N = 1000;

  out = fopen("wf3.dat", "w");

  Message("  Saving three-body wavefunction.\n");
  max = Lhalf;
  min = 0.;

  dx = (max-min) / (N-1);

  fprintf(out, "#N= %i\n", N);
  fprintf(out, "#1)x 2)f 3)lnf \n");
  fprintf(out, "#1x 2f 3lnf 4fp 5Eloc 6E\n");
  fprintf(out, "#1)x 2)f 3)lnf 4)fp 5)-f''/f - 5/r f'/f 6)-f''/f - 5/r f'/f + (f'/f)^2\n"); // Eloc = -f''/f - 2/r f'/f + (f'/f)^2

  for(i=0; i<N; i++) {
    x = min+dx*(DOUBLE) i;
    testNaN = exp(u3(x));
    if(testNaN == testNaN && testNaN != 2*testNaN) { // check for NaN by comparison to itself and INF
      fprintf(out, "%.15e %.15e %.15e %.15e %.15e %.15e\n", x,//1
        exp(u3(x)), //2
        u3(x),//3
        Fp3(x), //4
        E3(x) - Fp3(x)*Fp3(x), // 5
        E3(x)); // 6
    }
  }
  fclose(out);
#endif
}

/****************************** SaveWaveFunction *****************************/
void SaveImpurityWaveFunction(void) {
#ifdef ONE_BODY_IMPURITY
  int i;
  FILE *out;
  DOUBLE dx, x;
  DOUBLE f, lnf, fp, Eloc, E;
  DOUBLE testNaN;
  DOUBLE max,min;
  int N = 1000;

  out = fopen("wf1.dat", "w");

  Message("Saving one-body impurity wavefunction.\n");
  max = Lhalf;
  min = 0.;

  dx = (max-min) / (N-1);

  fprintf(out, "#N= %i\n", N);
  fprintf(out, "#1)x 2)f 3)lnf \n");
  fprintf(out, "#1x 2f 3lnf 4fp 5Eloc 6E\n");
  fprintf(out, "#1)x 2)f 3)lnf 4)fp 5)-f''/f - 5/r f'/f 6)-f''/f - 5/r f'/f + (f'/f)^2\n"); // Eloc = -f''/f - 2/r f'/f + (f'/f)^2

  for(i=0; i<N; i++) {
    x = min+dx*(DOUBLE) i;
    testNaN = exp(ImpurityInterpolateU(x));
    if(testNaN == testNaN && testNaN != 2*testNaN) { // check for NaN by comparison to itself and INF
      fprintf(out, "%.15e %.15e %.15e %.15e %.15e %.15e\n", x,//1
        exp(ImpurityInterpolateU(x)), //2
        ImpurityInterpolateU(x),//3
        ImpurityInterpolateFp(x), //4
        ImpurityInterpolateE(x) - ImpurityInterpolateFp(x)*ImpurityInterpolateFp(x), // 5
        ImpurityInterpolateE(x)); // 6
    }
  }
  fclose(out);
#endif
}

/****************************** SaveSk ***************************************/
int SaveStaticStructureFactor(void) {
  int i;
  FILE *out;
  DOUBLE normalization;
  static int initialized = OFF;
  DOUBLE rhok;
  static int add_column_with_rho_k_subtracted = OFF;

  // Save mixed/variational estimator
  out = fopen(file_Sk, (initialized || file_append)?"a":"w");
  if(out == NULL) {
    perror("\nError:");
    Warning("can't save S(k) function to file %s\n", file_Sk);
    return 1;
  }

  // check if the function was measured
  if(Sk.times_measured == 0) {
    Warning("Attempt to save static structre factor without any measurement\n");
    return 1;
  }

#ifdef BC_ABSENT
  normalization = 1./ ((DOUBLE)(Sk.times_measured));
#endif

#ifdef BC_1DPBC_ANY
  normalization = 1./ ((DOUBLE)(Sk.times_measured));
#endif

#ifdef BC_2DPBC
  normalization = 1./ ((DOUBLE)(Sk.times_measured));
  //normalization = N*Sk.N[0]/Sk.f[0];
#endif

#ifdef BC_3DPBC
  //normalization = N*Sk.N[0]/Sk.f[0]/(DOUBLE)Sk.times_measured;
  normalization = 1./(DOUBLE)Sk.times_measured;
#endif

#ifdef BC_SPHERE
  normalization = 1./(DOUBLE)Sk.times_measured;
#endif

  // first time, check if <rho_k> contribution is different from zero for any value of k
  if(initialized == OFF) {
    for(i=0; i<Sk.size; i++) {
      rhok = (normalization*Sk.cos[i]*normalization*Sk.cos[i])/(DOUBLE)N+(normalization*Sk.sin[i]*normalization*Sk.sin[i])/(DOUBLE)N;
      if(rhok>1e-3) {
        Warning("  <rho_k> = %lf contribution to S(k = %lf) is different from zero\n", rhok, Sk.k[i]);
        add_column_with_rho_k_subtracted = ON;
      }
    }
  }

  for(i=0; i<Sk.size; i++) {
#ifdef BC_ABSENT // BC absent
    if(add_column_with_rho_k_subtracted) { // Save <rho_{k} rho_{-k}> and also subtract <rho_k>^2
      rhok = (normalization*Sk.cos[i]*normalization*Sk.cos[i])/(DOUBLE)N+(normalization*Sk.sin[i]*normalization*Sk.sin[i])/(DOUBLE)N;
      fprintf(out, "%" LF " %.15" LE " %.15" LE "\n", Sk.k[i]/n, normalization*Sk.f[i]/(DOUBLE)N-(normalization*Sk.cos[i]*normalization*Sk.cos[i])/(DOUBLE)N-(normalization*Sk.sin[i]*normalization*Sk.sin[i])/(DOUBLE)N, normalization*Sk.f[i]/(DOUBLE)N);
    }
    else {
      fprintf(out, "%" LF " %.15" LE "\n", Sk.k[i]/n, normalization*Sk.f[i]/(DOUBLE)N-(normalization*Sk.cos[i]*normalization*Sk.cos[i])/(DOUBLE)N-(normalization*Sk.sin[i]*normalization*Sk.sin[i])/(DOUBLE)N);
    }
#endif

#ifdef BC_1DPBC_ANY // rescale by density 1D PBC
    if(measure_SkMATRIX) { // measure (kx, ky) matrix
      fprintf(out, " %" LE, normalization*Sk.f[i]/(DOUBLE)N);
      if((i+1)%(int)sqrt(Sk.size) == 0) fprintf(out, "\n");
    }
    else {
      if(add_column_with_rho_k_subtracted) { // Save <rho_{k} rho_{-k}> and also subtract <rho_k>^2
        rhok = (normalization*Sk.cos[i]*normalization*Sk.cos[i])/(DOUBLE)N+(normalization*Sk.sin[i]*normalization*Sk.sin[i])/(DOUBLE)N;
        fprintf(out, "%" LF " %.15" LE " %.15" LE "\n", Sk.k[i]/n, normalization*Sk.f[i]/(DOUBLE)N-(normalization*Sk.cos[i]*normalization*Sk.cos[i])/(DOUBLE)N-(normalization*Sk.sin[i]*normalization*Sk.sin[i])/(DOUBLE)N, normalization*Sk.f[i]/(DOUBLE)N);
      }
      else {
        fprintf(out, "%" LF " %.15" LE "\n", Sk.k[i]/n, normalization*Sk.f[i]/(DOUBLE)N-(normalization*Sk.cos[i]*normalization*Sk.cos[i])/(DOUBLE)N-(normalization*Sk.sin[i]*normalization*Sk.sin[i])/(DOUBLE)N);
      }
    }
#endif

#ifdef BC_2DPBC // rescale by density 2D PBC
    if(measure_SkMATRIX) { // measure (kx, ky) matrix
      //if(i == 0) Sk.f[0] = 0.;
      fprintf(out, "%.15" LF " ", Sk.N[i]?(normalization*Sk.f[i]/(Sk.N[i])):0);
      if(i<Sk.size-1 && Sk.ky[i+1]==0) fprintf(out,"\n");
    }
    else {
      //if(Sk.degeneracy[i] && Sk.k[i]>0) fprintf(out, "%" LF " %.15" LE "\n", Sk.k[i], Sk.N[i]?(normalization*Sk.f[i]/(Sk.N[i])):0);
      if(Sk.degeneracy[i] && Sk.k[i]>0) fprintf(out, "%" LF " %.15" LE "\n", Sk.k[i], normalization*Sk.f[i]/(N*Sk.degeneracy[i]));
    }
    Sk.N[i] = 0;
#endif

#ifdef BC_3DPBC // rescale by density 3D PBC
    if(Sk.degeneracy[i] && Sk.k[i]>0) fprintf(out, "%" LF " %.15" LE "\n", Sk.k[i], normalization*Sk.f[i]/(N*Sk.degeneracy[i]));
    Sk.N[i] = 0;
#endif
    Sk.f[i] = 0.;
    Sk.cos[i] = 0.;
    Sk.sin[i] = 0.;
  }
  Sk.times_measured = 0;
  fclose(out);

  // Save pure estimator
  if(measure_Sk && MC == DIFFUSION) {
    out = fopen(file_Sk_pure, (initialized || file_append)?"a":"w"); // open file
    if(out == NULL) {
      perror("\nError:");
      Warning("can't save S(k) function to file %s\n", file_Sk_pure);
      initialized = ON;
      return 1;
    }
    // check if the function was measured
    if(Sk_pure.times_measured == 0) {
      Warning("Attempt to save pure static structure factor without any measurement\n");
      initialized = ON;
      return 1;
    }
    // calculate the normalization factor
#ifdef BC_ABSENT
    normalization = 1./((DOUBLE)(Sk_pure.times_measured));
#endif

#ifdef BC_1DPBC_ANY
    normalization = 1./((DOUBLE)(Sk_pure.times_measured));
#endif

#ifdef BC_2DPBC
    //normalization = N*Sk_pure.N[0]/Sk_pure.f[0];
    normalization = 1./((DOUBLE)(Sk_pure.times_measured));
#endif

#ifdef BC_3DPBC
    //normalization = N*Sk_pure.N[0]/Sk_pure.f[0]/(DOUBLE)Sk_pure.times_measured;
    normalization = 1./(DOUBLE)Sk_pure.times_measured;
#endif

    // save data
    for(i=0; i<Sk.size; i++) {
#ifdef BC_ABSENT // PBC absent
      //fprintf(out, "%" LF " %.15" LE "\n", Sk.k[i], normalization*Sk_pure.f[i]/N-(normalization*Sk_pure.cos[i]*normalization*Sk_pure.cos[i])/N);
      fprintf(out, "%" LF " %.15" LE "\n", Sk.k[i], normalization*Sk_pure.f[i]/N);
#endif

#ifdef BC_1DPBC_ANY // rescale by density 1D PBC
      //fprintf(out, "%" LF " %.15" LE "\n", Sk.k[i]/n, normalization*Sk_pure.f[i]/N);
      fprintf(out, "%" LF " %.15" LE "\n", Sk.k[i]/n, normalization*Sk_pure.f[i]/N - N*(sin(Sk.k[i] * Lhalf) / (Sk.k[i] * Lhalf)) * (sin(Sk.k[i] * Lhalf) / (Sk.k[i] * Lhalf)));
#endif

#ifdef BC_2DPBC // rescale by density 2D PBC
      if(Sk.degeneracy[i] && Sk.k[i]>0) 
        //fprintf(out, "%" LF " %.15" LE "\n", Sk.k[i], Sk_pure.N[i]?(normalization*Sk_pure.f[i]/(Sk_pure.N[i])):0);
        fprintf(out, "%" LF " %.15" LE "\n", Sk.k[i], normalization*Sk_pure.f[i]/(N*Sk.degeneracy[i]));
      Sk_pure.N[i] = 0;
#endif

#ifdef BC_3DPBC // rescale by density 3D PBC
      if(Sk.degeneracy[i] && Sk.k[i]>0) 
        fprintf(out, "%" LF " %.15" LE "\n", Sk.k[i], normalization*Sk_pure.f[i]/(N*Sk.degeneracy[i]));

      Sk_pure.N[i] = 0;
#endif
      Sk_pure.f[i] = 0.;
    }
    fclose(out);
    Sk_pure.times_measured = 0;
  }

  initialized = ON;
  return 0;
}

/***************************** Save Superfluid Density *********************/
int SaveSuperfluidDensity(void) {
  FILE *out;
  int i;
  DOUBLE Dreal;
  static int initialized = OFF;

  // check if the function was measured
  if(SD.times_measured == 0) {
    Warning("Attempt to save superfluid fraction without any measurement\n");
    return 1;
  }

  out = fopen(file_SD, (initialized || file_append)?"a":"w");
  initialized = ON;
  if(ferror(out)) {
    Warning("can't save superfluid density\n");
    return 1;
  }

  for(i=0; i<SD.size-1; i++) {
    Dreal = (DOUBLE) N * SD.CM2[i] / ((DOUBLE)(i+1) * dt * SD.spacing * SD.times_measured);
    fprintf(out, "%" LG " %.15" LE "\n", (DOUBLE)(i+1)*dt*SD.spacing, Dreal);
    //if(SD.times_measured != SD.N[i]) Warning("Save SD   %i %i\n", SD.times_measured, SD.N[i]);
  }

  ArrayEmpty1D(SD.CM2, i, SD.size);
  SD.times_measured = 0;

  fclose(out);
  return 0;
}

/***************************** Save Superfluid Density *********************/
int SaveSuperfluidDensityArray(void) {
  FILE *out;
  int i;

  out = fopen("outsdarr.dat","w");
  if(ferror(out)) {
    Warning("can't save superfluid density array\n");
    return 1;
  }

  fprintf(out, "%i\n", SD.times_measured);
  for(i=0; i<SD.size; i++) fprintf(out, "%.15" LE " %i\n", SD.CM2[i], SD.N[i]);

  fclose(out);
  return 0;
}

/***************************** Load Superfluid Density *********************/
int LoadSuperfluidDensityArray(void) {
  FILE *in;
  int i;

  in = fopen("outsdarr.dat","r");
  if(in == NULL) {
    Warning("  Can't load superfluid density array. Skipping.\n");
    for(i=0; i<SD.size; i++) {
      SD.CM2[i] = 0.;
      SD.N[i] = 0;
    }
    return 1;
  }
  else {
    fscanf(in, "%i\n", &SD.times_measured);
    for(i=0; i<SD.size; i++) fscanf(in, "%" LF " %i\n", &SD.CM2[i], &SD.N[i]);
    fclose(in);
    return 0;
  }
}

/***************************** Save Order Parameter ************************/
void SaveOrderParameter(void) {
  FILE *out;
  static int initialized = OFF;
  DOUBLE norm;

  out = fopen("outop.dat", (initialized || file_append)?"a":"w");
  if(ferror(out)) {
    Warning("can't save order parameter\n");
  }

#ifdef TRIAL_1D
  norm = N*OrderParameter.times_measured;
  fprintf(out, "%.15" LE "\n", Sqrt(OrderParameter.cos*OrderParameter.cos+OrderParameter.sin*OrderParameter.sin)/norm);
  fclose(out);
#else
  norm = 1./(DOUBLE)OrderParameter.times_measured;
  fprintf(out, "%.15" LE " %.15" LE "\n", OrderParameter.cos*norm, OrderParameter.sin*norm);
  fclose(out);
  if(MC == DIFFUSION) {
    out = fopen("outopp.dat", (initialized || file_append)?"a":"w");
    initialized = ON;
    if(ferror(out)) {
      Warning("can't save order parameter\n");
    }
    norm = 1./(DOUBLE)OrderParameter.times_measured_pure;
    if(OrderParameter.times_measured_pure) fprintf(out, "%.15" LE " %.15" LE "\n", OrderParameter.cos_pure*norm, OrderParameter.sin_pure*norm);
    fclose(out);
  }
#endif
  OrderParameter.cos = 0.;
  OrderParameter.sin = 0.;
  OrderParameter.times_measured = 0;
  if(MC == DIFFUSION) {
    OrderParameter.cos_pure = 0.;
    OrderParameter.sin_pure = 0.;
    OrderParameter.times_measured_pure = 0;
  }
  initialized = ON;
}

/************************* Load Parameters Sk ********************************/
void LoadParametersSk(void) {
  FILE *in,*out  = NULL;
  int i,j;
  int number;
  //int nx,ny,nz,n2;
  DOUBLE nx,ny,nz,n2;
  int index;
  DOUBLE k,kold;

  if(verbosity) Warning("Loading Sk parameters...\n  ");

  in = fopen(INPATH "in3Dkas.in", "r");
  if(in == NULL) {
    perror("\nError:");
    Error("\ncan't read %s", "in3Dkas.in");
  }

  index = 0; // index pointing to the first degenerate momenta
  kold = -1;
  Sk.size = gridSk;
  if(measure_SkMATRIX) { // do not load k from file, but generate a matrix
    Warning("  Changing Sk into matrix with spacing (2 pi/Lx, 0.015 2 pi /aho)\n");
    index = 0;
    number = (int) sqrt((int)(Sk.size+0.1));
    for(i=0; i<number; i++) {
      for(j=0; j<number; j++) {
        Sk.kx[index] = 2.*PI/Lx*(double)(i);
        Sk.ky[index] = 0.015*2.*PI*(double)(j);
        Sk.kz[index] = 0.;
        Sk.index[index] = index;
        Sk.k[index] = Sqrt(Sk.kx[index]*Sk.kx[index]+Sk.ky[index]*Sk.ky[index]);
        Message("  %i\t%" LG " \t%" LG "\n", index+1, Sk.kx[index], Sk.ky[index]);
        index++;
      }
    }
  }
  else {
#ifdef BC_1DPBC_X // HO in y direction
    Warning("  momentum in x direction is in inverse interparticle units, in y direction in osc. units\n");
    Message("  N\tkx\tky\n");
#endif

#ifndef TRIAL_1D
    Fopen(out, "outskindex.dat", "w", "sk index");
#endif

    for(i=0; i<Sk.size; i++) {
      fscanf(in, "%i", &number);
      fscanf(in, "%" LF "", &nx);
      fscanf(in, "%" LF "", &ny);
      fscanf(in, "%" LF "", &nz);
      fscanf(in, "%" LF "\n", &n2);
#ifdef BC_1DPBC_X // HO in y direction
      Sk.kx[i] = 2.*PI/(Lx/(DOUBLE) N)*nx;
      Sk.ky[i] = ny;
      Sk.index[i] = i;
      Sk.k[i] = Sqrt(Sk.kx[i]*Sk.kx[i]+Sk.ky[i]*Sk.ky[i]);
      Message("  %i\t%" LG " \t%" LG "\n", i+1, Sk.kx[i], Sk.ky[i]);
#else
      if(fabs(nx)>1e-8) Sk.kx[i] = PI/L_half_x*nx;
      if(fabs(ny)>1e-8) Sk.ky[i] = PI/L_half_y*ny;
      if(fabs(nz)>1e-8) Sk.kz[i] = PI/L_half_z*nz;
      k = Sqrt(Sk.kx[i]*Sk.kx[i]+Sk.ky[i]*Sk.ky[i]+Sk.kz[i]*Sk.kz[i]);
      Sk.k[i] = k;

      Message("%i  %lf %lf %lf\n", i+1, nx, ny, nz);
      if(k>kold)   // comment this line to skip degenerated values of momenta
        index = i;

      Sk.index[i] = index;
      Sk.degeneracy[index]++;
      kold = k;
#endif

#ifndef TRIAL_1D
      if(i) fprintf(out, "%.15"LE "\n", Sk.k[i]);
#endif
    }
    fclose(in);
  }
#ifndef TRIAL_1D
  fclose(out);
#endif

  if(verbosity) Warning("done\n");
}

/***************************** Save Lindeman Ratio *************************/
void SaveLindemannRatio(void) {
  FILE *out;
  static int initialized = OFF;
  static int initializedp = OFF;
  int n;

  if(lattice_length == 0) {
    Warning("  Cannot save Lindemann ratio as lattice length is not defined\n");
    return;
  }

  if(LindemannRatio.N == 0) {
    Warning("  Trying to save pure Lindemann ratio without any measurement\n");
    return;
  }

  out = fopen("outlind.dat", (initialized || file_append)?"a":"w");
  initialized = ON;
  if(ferror(out)) {
    Warning("can't save Lindemann ratio\n");
  }

  if(measure_Lindemann_Lozovik == OFF)
    fprintf(out, "%.15" LE "\n", Sqrt(LindemannRatio.F/LindemannRatio.N/(lattice_length*lattice_length)));
  else
    fprintf(out, "%.15" LE "\n", Sqrt(LindemannRatio.F/LindemannRatio.N/(lattice_length*lattice_length)));

  LindemannRatio.F = 0.;
  LindemannRatio.N = 0;

  fclose(out);

  if(MC == DIFFUSION) { // save pure estimator
    if(LindemannRatio.Npure[0] == 0) {
      Warning("  Trying to save pure Lindemann ratio without any measurement\n");
      return;
    }
    out = fopen("outlindp.dat", (initialized || file_append)?"a":"w");
    initializedp = ON;
    if(ferror(out)) {
      Warning("can't save pure Lindemann ratio\n");
    }

    for(n=0; n<grid_pure_block; n++) {
      fprintf(out, "%.15" LE " ", Sqrt(LindemannRatio.Fpure[n]/LindemannRatio.Npure[n]/(lattice_length*lattice_length)));

      LindemannRatio.Fpure[n] = 0.;
      LindemannRatio.Npure[n] = 0;
    }
  fprintf(out, "\n");

  fclose(out);
  }
}

/************************* Load Crystal Coordinates *************************/
int LoadCrystalCoordinates(void) {
#ifdef ONE_BODY_IMPURITY
  Crystal.size = 1;
  Crystal.x[0] = 0;
  Crystal.y[0] = 0;
  Crystal.z[0] = 0.5*L;
  Warning("  Impurity is positioned at L/2");
#endif

#ifdef CRYSTAL
  FILE *in;
//  DOUBLE r;
  int i;

  Message("  Loading crystal coordinates ...\n");

  in = fopen(INPATH "in3Dcryst.in", "r");
  if(in == NULL) {
    perror("\nError:");
    Error("\ncan't read %s", "in3Dcryst.in");
  }

#ifndef BC_ABSENT
  fscanf(in, "    Lx=%" LF " Ly=%" LF " Lz=%" LF "\n", &Lx, &Ly, &Lz);
  Warning("Box size is set from file to Lx=%lf, Ly=%lf, Lz=%lf, density n=%lf\n", Lx, Ly, Lz, (DOUBLE)N/(Lx*Ly*Lz));
  if(fabs((DOUBLE)N/(Lx*Ly*Lz)/n-1)>1e-3) Error("  wrong density!\n");
#endif

  Crystal.size = NCrystal;
  for(i=0; i<NCrystal; i++) {
#ifdef CRYSTAL_WIDTH_ARRAY
    fscanf(in, "%" LF " %" LF " %" LF " %" LF " %" LF " %" LF " %" LF "\n", &Crystal.x[i], &Crystal.y[i], &Crystal.z[i], &Crystal.Rx[i], &Crystal.Ry[i], &Crystal.Rz[i], &Crystal.weight[i]);
#else
    fscanf(in, "%" LF " %" LF " %" LF "\n", &Crystal.x[i], &Crystal.y[i], &Crystal.z[i]);
#endif
    //r = Sqrt(Crystal.x[i]*Crystal.x[i]+Crystal.y[i]*Crystal.y[i]);
    //Crystal.x[i] *= Rpar;
    //Crystal.y[i] *= Rpar;
  }
  fclose(in);

#ifdef CRYSTAL_WIDTH_ARRAY
  Crystal.Rx_av = 0.;
  Crystal.Ry_av = 0.;
  Crystal.Rz_av = 0.;
  for(i=0; i<NCrystal; i++) {
    Message("  Crystal [%i]:", i);
    Message(" x = %" LE , Crystal.x[i]);
    Message(" y = %" LE , Crystal.y[i]);
    Message(" z = %" LE , Crystal.z[i]);
    Message(" Rx = %" LE , Crystal.Rx[i]);
    Message(" Ry = %" LE , Crystal.Ry[i]);
    Message(" Rz = %" LE , Crystal.Rz[i]);
    Message(" weight %" LE " \n", Crystal.weight[i]);

    Crystal.Rx_av += Crystal.Rx[i];
    Crystal.Ry_av += Crystal.Ry[i];
    Crystal.Rz_av += Crystal.Rz[i];
  }
  Crystal.Rx_av /= (double) NCrystal;
  Crystal.Ry_av /= (double) NCrystal;
  Crystal.Rz_av /= (double) NCrystal;
#endif

  Message("  done\n");
#endif

  return 0;
}

/******************************** Save SKT ****************************/
int SaveFormFactor(void) {
  FILE *out;
  static int first_time = ON;  // variable i tells how many times the energy was saved
  int i,k;

  out = fopen("outsktau.dat", (first_time == ON && file_append == OFF)?"w":"a");
  if(out == NULL) {
    perror("\nError:");
    Warning("can't write to file outsktau.dat\n");
    return 1;
  }

  for(i=0; i<gridSKT_t; i++) {
    fprintf(out, "%.15e ", i*dt*(DOUBLE)Nmeasure);
#ifdef TRIAL_1D
    for(k=0; k<gridSKT_k; k++) fprintf(out, "%.15" LE " ", SKT.f[i][k] / (DOUBLE) (SKT.times_measured*N));
#else
    fprintf(out, "%lf ", (DOUBLE)(i+1)*SKT_dk*PI*n);
    //fprintf(out, "0 "); // zero momentum
    for(k=0; k<gridSKT_k-1; k++) { // skip 0 point
      fprintf(out, "%.15" LE " ", SKT.f[i][k]/(DOUBLE)(SKT.times_measured*N));
      //fprintf(out, "%.15lf %.15lf\n", Sk.k[k], SKT.f[i][k] / (DOUBLE) (SKT.times_measured*N));
      //fprintf(out, "%.15lf %.15lf\n", Sk.k[k], N*SKT.f[i][k]/SKT.f[0][0]);
    }
#endif
    fprintf(out, "\n");

    for(k=0; k<gridSKT_k; k++) SKT.f[i][k] = 0.;
    for(k=0; k<gridSKT_k; k++) PhiTau.f[i][k] = 0.;
  }
  SKT.times_measured = 0;
  fclose(out);

  if(first_time == ON) {
    out = fopen("outsktauh.dat","w");
#ifdef TRIAL_1D
    fprintf(out, "%i %i %.15" LE " %.15" LE "\n", gridSKT_t, gridSKT_k, dt*(DOUBLE)Nmeasure, SKT_dk);
#else
    fprintf(out, "%i %i %.15" LE " %.15" LE "\n", gridSKT_t, gridSKT_k-1, dt*(DOUBLE)Nmeasure, SKT_dk);
#endif
    fclose(out);
    first_time = OFF;
  }

  first_time = OFF;

  return 0;
}

/************************* Save Defect Barrier *************************/
int SaveDefectBarrier(void) {
  Error("SaveDefectBarrier not implemented");

/*#ifdef CRYSTAL
  FILE *out;
  int i, j;
  DOUBLE x, y, E, Emin, u, u0;
  DOUBLE dx, dy;
  int M = 100; // matrix size M x M

  // for a fixed density
  //dx = dy = 3*M/sqrt(n);
  // to plot complete sample
  dx = Lx/(DOUBLE) M;
  dy = Ly/(DOUBLE) M;

  Message("Saving defect barrier profile to file barrier.dat ...\n  ");

  out = fopen(INPATH "barrier.dat", "w");
  if(out == NULL) {
    perror("\nError:");
    Error("\ncan't save");
  }

  for(i=0; i<NCrystal; i++) {
    wp.x[i] = Crystal.x[i]; 
    wp.y[i] = Crystal.y[i]; 
    wp.z[i] = Crystal.z[i]; 
  }

  // find Emin
  Emin = 1e10;
  for(i=0; i<M; i++) {
    x = (DOUBLE)(i)*dx;
    for(j=0; j<M; j++) {
      y = (DOUBLE)(j)*dy;

      wp.x[N-1] = x;
      wp.y[N-1] = y;
      wp.z[N-1] = 0;

      E = WalkerEnergy0(&Wp[w]);
      if(E<Emin) Emin = E;
    }
  }

  u0 = U(wp);
  for(i=0; i<M; i++) {
    x = (DOUBLE)(i)*dx;
    for(j=0; j<M; j++) {
      y = (DOUBLE)(j)*dy;

      wp.x[N-1] = x;
      wp.y[N-1] = y;
      wp.z[N-1] = 0;

      E = WalkerEnergy0(&Wp[w]);
      u = U(wp);

      fprintf(out, "%" LF " %" LF " %" LF " %" LF " %" LF  "\n", x, y, E-Emin, E, u-u0);
    }
    fprintf(out, "\n");
  }
  fclose(out);

  Message("done\n");
#endif*/

  return 0;
}

/*********************************** Save Lattice ****************************/
void SaveLattice(void) {
#ifdef BC_2DPBC_SQUARE

  FILE *out;
  int i, j, Nmax;
  DOUBLE x1,x2,x3,Fx,Fy,Fz;

  Message("\n  Saving lattice and w.f. ... ");
  out = fopen("lattice.dat","w");
  fprintf(out, "x y Vpot lnf Fx Fy Eloc\n");
  if(out == NULL) Error("\nError: Cannot open the file.");
  Lmax = Lx;
  Nmax = 100;
  for(i=0; i<Nmax; i++) {
    x1 = ((DOUBLE) i+1.)/(DOUBLE) (Nmax) * Lmax;
    for(j=0; j<Nmax; j++) {
      x2 = ((DOUBLE) j+1.)/(DOUBLE) (Nmax) * Lmax;

      Fx = Fy = Fz = 0;
      x3 = 0.;
      OneBodyFp(&Fx, &Fy, &Fz, x1, x2, x3,0,0);

      fprintf(out,"%.16" LE " %.16" LE " %.16" LE " %.4" LG " %.4" LG " %.4" LG " %.4" LG "\n", x1, x2, Vext(x1,x2,x3), OneBodyU(x1,x2,x3,0,0), Fx, Fy, OneBodyE(x1,x2,x3,0,0));
    }
    fprintf(out,"\n");
  }
  fclose(out);
#endif

#ifdef BC_1DPBC_Z
  FILE *out;
  int i, Nmax;
  DOUBLE x,y,z,Fx,Fy,Fz;

  x=y=z=0;
  Lmax = 2.; // show only the first two periods

  Message("\n  Saving lattice and w.f. ... ");
  out = fopen("lattice.dat","w");
  fprintf(out, "#z Vpot lnf Fz Eloc\n");
  if(out == NULL) Error("\nError: Cannot open the file.");
  Nmax = 100;
  for(i=0; i<Nmax; i++) {
    z = ((DOUBLE) i+1.)/(DOUBLE) (Nmax) * Lmax;
    Fx = 0;
    Fy = 0;
    Fz = 0;
    OneBodyFp(&Fx, &Fy, &Fz, x, y, z, 0, 0);

    fprintf(out,"%.8" LE " %.8" LE " %.8" LG " %.8" LG " %.8" LG "\n", z, Vext(x,y,z), OneBodyU(x,y,z,0,0), Fz, OneBodyE(x,y,z,0,0));
  }
  fclose(out);
#endif
}

/*********************************** Save Vext *******************************/
void SaveVext(void) {
  FILE *out;
  int i, Nmax;
  DOUBLE x,y,z;

  x=y=z=0;

  Message("\n  Saving external potential ... ");
  out = fopen("Vext.dat","w");
  fprintf(out, "#z Vext\n");
  if(out == NULL) Error("\nError: Cannot open the file.");
  Nmax = 100;
  for(i=0; i<Nmax; i++) {
    z = ((DOUBLE) i+1.)/(DOUBLE) (Nmax) * L;
    fprintf(out,"%.8" LE " %.8" LG "\n", z, Vext(x,y,z));
  }
  fclose(out);
}

/*********************************** Save Spin Polarization ******************/
void SaveSpinPolarization(void) {
  FILE *out;
  static int first_time = ON; 

  out = fopen("outspin.dat", (first_time == ON && file_append == OFF)?"w":"a");
  fprintf(out, "%lf\n", SpinPolarization());
  fclose(out);

  first_time = OFF;
}

/******************************** Save Hyper Radius *******************/
int SaveHyperRadius(void) {
  int i;
  DOUBLE r;
  FILE *out;
  DOUBLE normalization;
  static int initialized = OFF;

  if(HR.times_measured == 0) {
    Warning("Attempt to save hyper radius without any measurement\n");
    return 1;
  }

  Fopen(out, "outhr.dat", (initialized || file_append)?"a":"w", "hyper radius");

  //normalization = 1./(n*n/(2.*PI)*(DOUBLE)(HR.times_measured)*HR.step);
  normalization = 1./((DOUBLE)(HR.times_measured)*HR.step);

  for(i=0; i<HR.size; i++) {
    r = HR.min + (DOUBLE) (i+0.5) * HR.step;

    fprintf(out, "%f %.15e\n", r, normalization*HR.f[i]);
    HR.N[i] = 0;
    HR.f[i] = 0.;
  }
  HR.times_measured = 0;
  fclose(out);

  HR.times_measured = 0;
  initialized = ON;

  return 0;
}


/*********************************** Save Data *******************************/
void SaveBlockData(int block) { // save distributions
#ifdef MPI
  time_done_in_parallel += MPI_Wtime();

  CollectAllWalkersData();

  if(myid == MASTER) {
#else // do not save, as it is parallel
  SaveCoordinates(file_particles);
#endif

  if(measure_OBDM) SaveOBDM();
  if(measure_OBDM_MATRIX) SaveOBDMMatrix();
  if(measure_TBDM) SaveTBDM();
  if((measure_OBDM_MATRIX || measure_OBDM) && DIMENSION == 1) SaveMomentumDistribution();
  if(measure_PairDistr) SavePairDistribution();
  if(measure_g3) SaveHyperRadius();
  if(measure_RadDistr) SaveRadialDistribution();
  if(measure_RadDistr && measure_RDz) SaveRadialZDistribution();
  if(measure_RadDistrMATRIX) SaveRadialDistributionMatrix();
  if(measure_Sk) SaveStaticStructureFactor();
  if(MC == DIFFUSION && measure_SD) SaveSuperfluidDensity();
  //if(MC == DIFFUSION && measure_SD) SaveSuperfluidDensityArray();
  if(measure_OP) SaveOrderParameter();
  if(measure_Lind) SaveLindemannRatio();
  if(MC == DIFFUSION && measure_FormFactor) SaveFormFactor();
  if(measure_effective_potential) SaveEffectivePotential();
  if(MC == DIFFUSION) SaveEpotPure();

#ifdef MPI
  }
  EmptyDistributionArrays();
  time_done_in_parallel -= MPI_Wtime();
#endif
}

/*********************** Save Effective Potential ****************************/
int SaveEffectivePotential(void) {
  int i;
  DOUBLE r;
  FILE *out, *out2;
  int save_pure = ON;
  char file_name[40] = "outveff.dat";
  DOUBLE normalization;
#ifdef BC_ABSENT
#ifdef TRIAL_3D
  DOUBLE norm,f;
#endif
#endif
  static int initialized = OFF;

  if(Veff.times_measured == 0) {
    Warning("Attempt to save effective potential without any measurement\n");
    return 1;
  }

  out = fopen(file_name, (initialized || file_append)?"a":"w");
  if(out == NULL) {
    perror("\nError:");
    Warning("can't save effective potential to file %s\n", file_name);
    return 1;
  }

  for(i=0; i<Veff.size; i++) {
    r = Veff.min + (DOUBLE) (i+0.5) * Veff.step;
    fprintf(out, "%" LF " %.15" LE "\n", r, (Veff.N[i]>0)?(Veff.f[i]/(DOUBLE)Veff.N[i]):(0.));
    Veff.N[i] = 0;
    Veff.f[i] = 0.;
  }
  Veff.times_measured = 0;
  fclose(out);

  initialized = ON;

  return 0;
}

/****************************** Save Pure Coordinates ************************/
void CopyPureCoordinates(void) {
  int w, i;

  for(w=0; w<Nwalkers; w++) {
    for(i=0; i<N; i++) {
      CaseX(W[w].x_pure[i] = W[w].x[i]);
      CaseY(W[w].y_pure[i] = W[w].y[i]);
      CaseZ(W[w].z_pure[i] = W[w].z[i]);
    }
  }
}

int SavePureCoordinates(void) {
  FILE *out;
  int w, i;
  static int first_time = ON;

  Fopen(out, "in3Dpure.dat", first_time?"w":"a", "coordinates");

  if(MC == DIFFUSION) { // save 1st walker
    w = 0;
    fprintf(out, "%g\n", W[w].E);
    for(i=0; i<N; i++) fprintf(out, "%.15f %.15f %.15f\n", W[w].x_pure[i], W[w].y_pure[i], W[w].z_pure[i]);
  }
  else  { // VARIATIONAL, save all walkers
    for(w=0; w<Nwalkers; w++) { // Save all walkers
      fprintf(out, "%g\n", W[w].E);
      for(i=0; i<N; i++) fprintf(out, "%.15f %.15f %.15f\n", W[w].x_pure[i], W[w].y_pure[i], W[w].z_pure[i]);
    }
  }
  fclose(out);

  CopyPureCoordinates();

  first_time = OFF;

  return 0;
}

/****************************** Save Epot Pure *******************************/
int SaveEpotPure(void) {
  int i;
  FILE *out;
  static int first_time = ON;
  DOUBLE normalization;

  Fopen(out, "outepot.dat", first_time?"w":"a", "pure Epot");
  //normalization = energy_unit / ((DOUBLE)Epot_pure.times_measured * N);
  normalization = energy_unit / ((DOUBLE)Epot_pure.times_measured);

  for(i=0; i<Epot_pure.size; i++) fprintf(out, "%.15e %.15e\n", Epot_pure.x[i], Epot_pure.f[i] * normalization);

  fclose(out);

  for(i=0; i<Epot_pure.size; i++) Epot_pure.f[i] = 0.;
  Epot_pure.times_measured = 0;

  first_time = OFF;

  return 0;
}
