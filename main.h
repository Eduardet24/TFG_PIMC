#ifndef __MAIN_H_
#define __MAIN_H_

#define SRC_CODE_DATE "Monte Carlo for bosons. Source code of 23/IX/2026\n\n"

/* define SECURE to enable all self-consistent checks */
//#define SECURE

//#define MPI // define MPI for parallel computations
//#define NPARTICLES 10 // in parallel or static case N must be specified explicitly
#define OPENMP_MAX_NPROC 272 // 2048 maximal number of nodes which can be used in Open MP; compile with /openmp option

//#define DATATYPE_FLOAT
#define DATATYPE_DOUBLE
//#define DATATYPE_LONG_DOUBLE

#define IMAGES_FLOOR
//#define IMAGES_IF_IF
//#define IMAGES_ABS_IF

//#define CALCULATE_DERIVATIVES_NUMERICALLY // only log of w.f. is needed

//#define MEMORY_CONTIGUOUS // use static allocation for walkers, otherwise dynamic array
//#define BRANCHING_ALWAYS_STABLE // force branching that cannot collapse
//#define BRANCHING_FIXED_NUMBER_OF_WALKERS
#define MEASURE_CORRELATION_FUNCTIONS_IN_THE_MIDDLE_OF_DRIFT
//#define MEASURE_CORRELATION_FUNCTIONS_AT_THE_END_OF_DRIFT
//#define MEASURE_PURE_OBDM // add virtual walkers for pure OBDM or pure nk measurements
//#define MEASURE_PURE_TBDM // add virtual walkers for pure OBDM or pure nk measurements
//#define MEASURE_ENERGY_FROM_GROWTH_ESTIMATOR // additional energy estimator for DMC algorithm
#define MEASURE_FORM_FACTOR_FOR_SUM_RULES
//#define DO_NOT_CALCULATE_ENERGY_AND_U // use for VMC for the calculation of the correlation functions only

/* define ONLY ONE type of the interaction potential (which is not included in Eloc!)*/
//#define HARD_SPHERE // Check |r|>a condition, can be avoided if w.f. is very small for |r|<a 
//#define HARD_SPHERE_KILL_OVERLAPPED // DMC chose to kill or rethrow an overlapping walker
//#define HARD_SPHERE_HYPERRADIUS // kill all walkers with r3<R3, is implemented for VMC Move All, DMC quadratic
//#define KURBAKOV_SOFT_CORE // U(r)+=(a/r)^100*(100*exp(-gammaE)/r)^2 // Kurbakov 9-CORE
//#define HARD_SPHERE_EXCLUDED_VOLUME_CHANGE_OF_VARIABLES // use Girardeau's trick for simulation of hard rods (VMC MoveAll) works with TRIAL_TONKS_SIMPLE and ONE_BODY_TRIAL_ZBC

//#define SPINFULL // define if multiple components are used
//#define SPINFULL_TUNNELING // allow spins to tunnel
//#define SPINFULL_GENERATE_COORDENATES_IN_CHAIN 
//#define OBDM_FERMIONS
//#define OBDM_ANYONS

#define INTERACTION_ABSENT
//#define INTERACTION_Aziz // uses energy_unit inside InteractionEnergy call, to have [K] as local energy unit
//#define INTERACTION_HYDROGEN // uses energy_unit inside InteractionEnergy call, to have [K] as local energy unit
//#define INTERACTION_YUKAWA
//#define INTERACTION_LENNARD_JONES
//#define INTERACTION_LENNARD_JONES10
//#define INTERACTION_SUTHERLAND
//#define INTERACTION_CALOGERO // 1D: D(D-1) / r^2, 2D: D / r^2, 3D: D(D+1) / r^2
//#define INTERACTION_K0
//#define INTERACTION_PETROV
//#define INTERACTION_PETROV3D
//#define INTERACTION_COULOMB_1D // - (digamma(r/L) + digamma(1.-r/L)) / L
//#define INTERACTION_LINEAR_1D_COULOMB // r/a;
//#define INTERACTION_WITH_DAMPING // Imaginary part of interaction potential is treated as a perturbation
//#define INTERACTION_ROTON
//#define INTERACTION_SOFT_SPHERE // Vo, Ro
//#define INTERACTION_SQUARE_WELL // of size RoSW
//#define INTERACTION_HARD_SPHERE_SQUARE_WELL // hard sphere size (a), Ro, square well size (depth Vo), Apar - scattering length
//#define INTERACTION_GAUSSIAN // a*exp(-b*r*r)
//#define INTERACTION_KURBAKOV //
//#define INTERACTION_KURBAKOV_GENERAL // uses D, Apar, Ipar
//#define INTERACTION_DOUBLE_GAUSSIAN // -gaussian_alpha*exp(-0.5*r*r) - gaussian_beta*exp(-0.25*r*r)
//#define INTERACTION_LOGARITHMIC
//#define INTERACTION_POWER // D/r^Apar, Apar>2
//#define INTERACTION_COULOMB // 1/r
//#define INTERACTION_R2 // 1/r^2
//#define INTERACTION_DIPOLE // D/r^3
//#define INTERACTION_QUADRUPOLE  // 1/r^5
//#define INTERACTION_SPINFULL_DIPOLE_BILAYER
//#define INTERACTION_SPINFULL_DIPOLE_MULTILAYER
//#define INTERACTION_SPINFULL_GAUSSIAN // exp(-x^2 / 2)
//#define INTERACTION_SPINFULL_GAUSSIANS_3COMPONENT // strength a11,a12,a22 exp(-x^2)
//#define INTERACTION_SPINFULL_DOUBLE_GAUSSIAN
//#define INTERACTION_DIPOLE_TRUNCATED
//#define INTERACTION_DIPOLE_SCREENED // 1/(r^3+a^3)
//#define INTERACTION_R4 // 1/r^4
//#define INTERACTION_RYDBERG // 1/r^6
//#define INTERACTION_RYDBERG_TRUNCATED
//#define INTERACTION_RYDBERG_SCREENED // 1/(r^6+a^6)
//#define INTERACTION_RYDBERG_WITH_DELTA //-1/2+sqrt(C6/r^6+1/2);
//#define INTERACTION_R12 // 1/r^12
//#define INTERACTION_Q1D_DIPOLES // D*(-2 x + sqrt(2*pi) (1 + x^2) Erfcx(x/sqrt(2)))
//#define INTERACTION_SOFT_SPHERE_HYPERRADIUS // all walkers with r3<R3 get energy +1
//#define INTERACTION_CH2 // -Epar/ch^2(r), modified Poschl Teller potential, Epar=2 corresponds to unitarity
//#define INTERACTION_GAPPED

//#define INTERACTION_SUM_OVER_IMAGES // enable to sum up to convergence all PBC images, sums up to sum_over_images_cut_off_num
//#define INTERACTION_EWALD // enable to use Ewald summation
//#define INTERACTION_SUM_OVER_IMAGES_SMOOTH_CUTOFF

//#define SCALABLE // define for duplicating the system Nc times in all directions
//#define SCALABLE_POTENTIAL	 // define for truncation of the potential
//#define SCALE_VMC_AMPLITUDE_WITH_ALPHA // VMC displacement amplitude will be scaled by 1/sqrt(alpha_Rx, alpha_Ry, alpha_Rz)

// put w.f. on a grid for a better performance (define to enable)
//s#define LOAD_WF_FROM_FILE
//#define SAMPLE_WF // resample w.f. to grid for a better performance

//#define INTERPOLATE_SPLINE_JASTROW_WF // spline interpolation, in case of HS set 0 to 1e-15
//#define INTERPOLATE_LOG // use ln f(r) rather than f(r)
//#define INTERPOLATE_LOG_ONE_BODY // use ln f1(r) for interpolation rather than f1(r) (always enable for speckles)
//#define INTERPOLATE_LINEAR_JASTROW_WF // linear interpolation

//#define SYMMETRIZE_TRIAL_WF // enable to use u(r) -> u(r)+u(L-r)
//#define POWER_TRIAL_WF // enable to use u(r) -> alpha_Jastrow u(r) (not compatible with SYMMETRIZE)

/* define ONLY ONE type of the trial wave function*/
#define  TRIAL_ABSENT11 // Ideal Bose gas
#define  TRIAL_ABSENT12 // Ideal Bose gas
#define  TRIAL_ABSENT22 // Ideal Bose gas
//#define  TRIAL22_SAME_AS_11 // if 22 is not implemented, take it as a copy of 11
//#define  TRIAL12_SAME_AS_11 // if 22 is not implemented, take it as a copy of 11
//#define  TRIAL11_SAME_AS_12 // if 11 is not implemented, take it as a copy of 12
//#define  TRIAL22_SAME_AS_12 // if 22 is not implemented, take it as a copy of 12
//#define  TRIAL_EXTERNAL /* Load trial wave function from 'wf.in', 'wfpar.in' files*/
//#define  TRIAL_TWO_BODY_SCATTERING_NUMERICAL // find two body scattering solution and interpolate it with a spline
//#define  TRIAL_TWO_BODY_SCATTERING_NUMERICAL_SPINFULL // same as before (max distance Rpar12) for 11,12,22, etc, define SPINFULL
//#define  TRIAL_KURBAKOV_NUMERIC  // uses 0<Rpar<1 matching distance, Epar scattering energy, Dpar effective dimensionality, Mpar effective mass, a<0 scattering length or a>0 hard core diameter, Kurbakov_BC

//#define  JASTROW_RIGHT_BOUNDARY_BOUND_STATE // use bound state (spinless case)
//#define  JASTROW_RIGHT_BOUNDARY_SPINFULL_SAME_SPIN_BOUND_STATE // use bound state asymptotics for same spin
//#define  JASTROW_RIGHT_BOUNDARY_SPINFULL_DIFFERENT_SPIN_BOUND_STATE // use bound state asymptotics for different spin
//#define  JASTROW_RIGHT_BOUNDARY_CONSTANT
//#define  JASTROW_RIGHT_BOUNDARY_PHONONS // (2D,3D) Rpar, Apar 
//#define  JASTROW_LEFT_BOUNDARY_ZERO // f(0) = 0 
//#define  JASTROW_LEFT_BOUNDARY_U_MINUS_INF_ZERO// u(t->-inf)/dt = 0, f(r) = sinh u (asinh (ln (r)))
//#define  JASTROW_LEFT_BOUNDARY_ZERO_DERIVATIVE // f'(0) = 0
//#define  JASTROW_LEFT_BOUNDARY_PSEUDOPOTENTIAL // (1D) f'(0)/f(0) = -1/a

//#define  JASTROW_DIMENSIONALITY_FROM_DPAR // use Dpar parameter instead of D=1,2,3
//#define  JASTROW_EFFECTIVE_REPULSIVE_INTERACTION_FROM_R0PAR //solve two-body problem spin1=spin2 with R0PAR/r^3 instead of 1/r^3
//#define  TRIAL_TWO_BODY_SCATTERING_PHONONS_NUMERICAL // find two body scattering solution with (2D) exp(-Apar11*(1./r+1./(L-r))], r>Rpar11
//#define  TRIAL_DIPOLE_PHONON_FINITE_D /* (2D) Ko(2sqrt(D/r)), r<Rpar, exp{-Apar*(1./r+1./(Lwf-r))}, r>Rpar */
//#define  TRIAL_SCHIFF_VERLET_TRAP /*use SYMMETRIZE_TRIAL_WF for PBC (1D-3D):(1D) exp(-(Apar/(n|z|))^Bpar) (2D, hom) exp{-[Apar/(r/L)]^Bpar} (3D) exp(-Apar/r^Bpar)*/
//#define  TRIAL_SCHIFF_VERLET_TRAP_PHONON2D /*use SYMMETRIZE_TRIAL_WF for PBC (2D) exp{-[Apar/(r/L)]^Bpar - Cpar*(1/r+1/(L-r))}*/
//#define  TRIAL_SCHIFF_VERLET_PHONON /* (1D-3D):(1D) exp(-(Apar/(n|z|))^Bpar) (2D, hom) exp{-(Apar/r)^Bpar} (e^-Cr+e^-C(r-L)) (3D) exp(-Apar/r^Bpar), r<Rpar Lhalf, phonons otherwise*/
//#define  TRIAL_DIPOLES_PHONON /* (2D):  Atrial*exp(-2./sqrt(r)), r<Rpar*Lhalf, Btrial*exp(-Ctrial*(1./r+1./(L-r))) r>Rpar*Lhalf */
//#define  TRIAL_POWER_Ko /*(2D) Atrial*K0^-alpha_1, r<Rpar*L/2, r<Rpar*Lhalf, Btrial*exp(-Ctrial*(1./r+1./(L-r))) r>Rpar*Lhalf, cutoff at small distances */
//#define  TRIAL_DIPOLE_Ko /*(2D) Atrial Ko(2/sqrt(r)), r<Rpar*Lhalf; Btrial*exp(-Ctrial*(1./r+1./(L-r))) r>Rpar*Lhalf*/
//#define  TRIAL_QUADRUPOLE_Ko /*(2D) ... ; Btrial*exp(-Ctrial*(1./r+1./(L-r))) r>Rpar*Lhalf*/
//#define  TRIAL_DIPOLE_Ko11 /*(2D) Atrial Ko(2/sqrt(r)), r<Rpar*Lhalf; Btrial*exp(-Ctrial*(1./r+1./(L-r))) r>Rpar*Lhalf*/
//#define  TRIAL_DIPOLE_Ko11_TRAP /*(2D)11 Ko(2*sqrt(D/r))
//#define  TRIAL_DIPOLE_Ko_TRAP /*(2D) Ko(2*sqrt(D/r))
//#define  TRIAL_HS /*(3D) 2-body HS solution + symmetrized phonons, Rpar[L/2], Bpar - scatt. momentum k [2pi/L]: sin[k(x-a)]/x; A*Exp(-B/x^2-B/(L-x)^2) */
//#define  TRIAL_PSEUDOPOTENTIAL_PHONON    /*(3D) 2-body HS solution + symmetrized phonons, Rpar[L/2], Bpar - scatt. momentum k [2pi/L]: sin[k(x-a)]/x; A*Exp(-B/x^2-B/(L-x)^2) */
//#define  TRIAL_HS_SIMPLE  /*(3D,2D) 3D HardSphere + A1SW/x * sin(sqE*(x-a)), no variational parameters, r>a */
//#define  TRIAL_HS_SW_SIMPLE  /*(3D) 3D HardSphere + Square Well A1SW/x * sin(sqE*(x-a)), no variational parameters, r>a */
//#define  TRIAL_PSEUDOPOTENTIAL_SIMPLE /*(3D) 3D HardSphere + A1SW/x * sin(sqE*(x-a)), no variational parameters, r>0, node at r=a, works for positive a */
//#define  TRIAL_HARD_RODS_3D /*(1D,2D,3D) sin[pi (r-1)/(L-2)]^Apar */
//#define  TRIAL_EXP_PHONON /* (2D):  exp{-[Apar/(r/L)]^Bpar}, r<Rpar*Lhalf, exp(-Ctrial*(1./r+1./(L-r))) r>Rpar*Lhalf */

//#define  TRIAL_TONKS /* (1D) Tonks w.f. |sin(sqrt{Bpar}*(x-a))|, Rpar=0, a=0 means TG */
//#define  TRIAL_TONKS_SIMPLE /* (1D) Tonks w.f. |sin(sqE*x)| with sqE = pi/L or pi/(L-N*Apar)*/
//#define  TRIAL_TONKS_TRAP11 //(1D,2D,3D) f(r) = |r-b|, enable HS11 to have f(r)=0 for r<b
//#define  TRIAL_TONKS_TRAP //(1D,2D,3D) f(r) = |r-a|, enable HS to have f(r)=0 for r<a
//#define  TRIAL_TONKS_TRAP_tVMC // (1D) f(x) = 1 - x/a  + (c3R + i c3I) x^2 + (c4R + i c4I) x^3 + (c5R + i c5I) x^4
//#define  TRIAL_CH2_3D_TRAP // (3D) th(r)/r, unitarity free state 
//#define  TRIAL_CH2_1D_TRAP_FS // (1D) th(r), excited state with one node
//#define  TRIAL_CH2_1D_TRAP_BS // (1D) 1/ch(r), bound state
//#define  TRIAL_TONKS_CUTOFF // (1D) |sin(pi x/Rpar)|, a<Rpar; 1,x>Rpar
//#define  TRIAL_SUTHERLAND /* (1D) sin^a r(pi*x/L), a=1/2 (1+(1+D)^1/2) */
//#define  TRIAL_SS_1D // (1D) exp(k*x) + Apar*exp(-k*x)) x<Ro, sin(pi*x/L)^(1/Kpar) x>Ro, var. parameters Apar, Kpar
//#define  TRIAL_MCGUIRE // (1D,3D) 1D: f(r) = exp(-r/a), 3D: f(r) = exp(-r/a)/r, bound state of a pseudopotential
//#define  TRIAL_MCGUIRE_PHONON // (1D) exp(x/b); (x<Rpar12), sin^a(pi*x/L) otherwise
//#define  TRIAL_CALOGERO /*(1D,2D,3D) |r|^lambda, lambda = D from cfg
//#define  TRIAL_CALOGERO_2D /*(2D PBC) Atrial r^2 + Btrial/r^2, r<Rpar*Lhalf; exp(Ctrial-Apar*(1./r+1./(L-r))) r>Rpar*Lhalf*/
//#define  TRIAL_COULOMB_1D_PHONON /*(1D) r^1/2I_1(2r^1/2), r<Rpar, sin^beta(pi r/L), r>Rpar*/
//#define  TRIAL_DIPOLE_1D_PHONON /*(1D) D=1, r^1/2 I_1(2/r^1/2), r<Rpar, sin^beta(pi r/L), r>Rpar*/
//#define  TRIAL_DIPOLE_1D_TRAP /* (1D) any D, (r/D)^1/2 I_1(2 (D/r)^1/2) */
//#define  TRIAL_COULOMB_1D_TRAP /*(1D) r^1/2I_1(2r^1/2)*/
//#define  TRIAL_COULOMB_2D_TRAP /*(2D) I_0[2(Dr)^1/2]*/
//#define  TRIAL_COULOMB_3D_PHONON /*(3D) I_1(2r^1/2)/r^1/2, r<Rpar, exp(-c/r^2-c/(L-r)^2), r>Rpar*/
//#define  TRIAL_HS1D  /* (1D) HS + HS solution + exponent in 1D (hard rode) */
//#define  TRIAL_LIEB  /* (1D) cos(k(|z|-Rpar)), with f'(0)/f(0) = - 1/a_{1D}, for Rpar=0 -> Rpar = L/2 */
//#define  TRIAL_LIEB_RENORMALIZED
//#define  TRIAL_PHONON_NEGATIVE /* (1D) cos(k(|z|-Rpar)) (x<L*Rpar), sin^a(pi*x/L) a<0, |a| is the node*/
//#define  TRIAL_PHONON /* (1D) Rpar: A cos(k(x-B)) (x<L*Rpar), sin^a(pi*x/L), i.e. Lieb + phonons mat, if 0<Rpar<1 treated as in units of L/2, otherwise is taken as it is*/
//#define  TRIAL_PHONON_LUTTINGER /* (1D) Kpar, Rpar: A cos(k(x-B)) (x<L*Rpar), |sin(pi*x/L)|^{1/Kpar}, i.e. Lieb + phonons, Rpar is taken as it is */
//#define  TRIAL_PHONON_LUTTINGER_LATTICE /* Same as above times (1 + Apar12 (tanh[Bpar12 (x - Rpar12)] + tanh[Bpar12 (-x - Rpar12)])), changes boundary condition, i.e. not for delta interaction */
//#define  TRIAL_PHONON_LUTTINGER_PIECEWISE /* */
//#define  TRIAL_LIEB_LUTTINGER_PHONON /* sin(Atrial + Btrial r + Ctrial r^2)^(1/Rpar), where f(0) = Apar and Rpar is equivalent to the Luttinger parameter K*/
//#define  TRIAL_LIEB_EXPONENTIAL_DECAY /* (1D) Atrial (r - a); 1 - Btrial Exp[-Apar r] - Btrial Exp[-Apar (L - r)], parameters Rpar and Apar */
//#define  TRIAL_PHONON_LUTTINGER_AUTO // (1D) Kpar (automatic Rpar): A cos(k(x-B)) (x<L*Rpar), sin^a(pi*x/L), i.e. Lieb + phonons
//#define  TRIAL_PHONON_LUTTINGER_AUTO11_SPINFULL // (1D) aA,  Kpar11 (automatic Rpar11): A cos(k(x-B)) (x<L*Rpar11), sin^a(pi*x/L), i.e. Lieb + phonons
//#define  TRIAL_PHONON_LUTTINGER_AUTO12_SPINFULL // (1D) aAB, Kpar12 (automatic Rpar12): A cos(k(x-B)) (x<L*Rpar11), sin^a(pi*x/L), i.e. Lieb + phonons
//#define  TRIAL_SS_PHONON_LUTINGER // (1D) Rpar,Kpar: A ch(kappa x) (x<Ro), B sin(kx+delta) (Ro<x<Rpar), sin^a(pi*x/L), a=1/Kpar
//#define  TRIAL_SS_PHONON_LUTINGER_AUTO // (1D) A ch(k x) (x<Rpar), sin^a(pi*x/L)
//#define  TRIAL_GAPPED /* (1D) gap D: Reatto Chester for gapped system, only f(x), not f'(x); f''(x); f(x) = exp{-sum_i L*L*D/(2.*pi*pi*N*i*i)*cos(2.*pi*x/L)} */
//#define  TRIAL_COS /* (1D) cos(Apar*pi*x/L)*/
//#define  TRIAL_SIN_EXCITED_BOSONS /* (1D) sin(Apar*pi*x/L + Bpar)^(1/Kpar)*/
//#define  TRIAL_COS_OFFSET /* (1D) Apar + Bpar*cos(Cpar*pi*x/L)*/
//#define  TRIAL_DECAY /*(2D) exp(-Apar/(1+Bpar*r^2)), for Ko(r) interaction (Magro et al.) */
//#define  TRIAL_DECAY_PHONON /*(2D) exp(-Apar/(1+Bpar*r^2)),  r<Rpar*Lhalf, Btrial*exp(-Ctrial*(1./r+1./(L-r))) r>Rpar*Lhalf */
//#define  TRIAL_YUKAWA_2D  /*(2D) */
//#define  TRIAL_HS2D  /* (2D) log(exp(-0.5*(1+Apar)*r)*(exp(Apar*r)-exp(Apar)))+shift;  0<Apar<1, Apar is close to 1 */
//#define  TRIAL_HS2D_XING /*(2D) th(Apar(r^Bpar-1)) */
//#define  TRIAL_HS2D_ZERO_ENERGY /* (2D) ln(r)*/
//#define  TRIAL_HS2D_FINITE_ENERGY /* (2D) parameter: a<Rpar<L/2 - matching distance to 1, use Rpar=0 for setting Rpar=L/2*/
//#define  TRIAL_2D_ZERO_RANGE_BOUND_STATE // (2D) Ko(k*r) where k = 2e^(-gamma)/a
//#define  TRIAL_2D_ZERO_RANGE_K0_SCREENING // (2D) Apar,b: Ko(k*r) exp(-r*r/(b*(Apar+r))) where k = 2e^(-gamma)/a
//#define  TRIAL_SPHERE_ZERO_RANGE_BOUND_STATE // (sphere) Ko(k*r) where k = 2e^(-gamma)/a, Eloc(r) = +f"(r)/f(r)
//#define  TRIAL_SPHERE_ZERO_RANGE_FREE_STATE // (sphere) ln(a/r), Eloc(r) = +f"(r)/f(r)
//#define  TRIAL_SPHERE_ZERO_RANGE_LOG_MATCH // (sphere) Rpar: ln(r/a), r<Rpar; ln(r/a) * exp((r/Rpar-1) / ln(r/a)), Eloc(r) = +f"(r)/f(r)
//#define  TRIAL_SPHERE_ZERO_RANGE_K0_SCREENING // (sphere) Apar,b: Ko(k*r) exp(-r*r/(b*(Apar+r))) where k = 2e^(-gamma)/a, Eloc(r) = +f"(r)/f(r)
//#define  TRIAL_SPHERE_ZERO_RANGE_LOG_SCREENING // (sphere) b: ln(a/r) * exp(-r/b), Eloc(r) = +f"(r)/f(r), a>R

//#define  TRIAL_MATRIX   /* Load trial wavefuncion (matrix) from 'wf2D.in file'*/

//#define  TRIAL_LIM   /* (1D,2D,3D) 1 - a/r  i.e. limit of 3D HS + HS solution */
//#define  TRIAL_POWER /* (1D,3D,3D) 1 - a/r^Apar */
//#define  TRIAL_PP_FS_TRAP_MCGUIRE1D // (3D) exp(1 - (Apar**2*b - a*r**2)/(a*Apar*b + a*b*r))/r, i.e. a<0: 3D scattering length, b: 1D scattering length, Apar - var. par
//#define  TRIAL_SS    /* (3D,2D): (3D PBC) SS + SS solution + exponent in 3D (hard sphere), (3D trap): sinh(Kr)/r, B-A/r D) symmetrize */
//#define  TRIAL_SS_LARGE /*(3D) 2b body SS + phonons, SS size Apar > Rpar*/
//#define  TRIAL_SS_EXP /*(3D) exp(-Apar r^2) + phonons,  Rpar*/

//#define  TRIAL_SQUARE_WELL_BS_PHONONS /* (3D) SS+phonons at Rpar, C1 sin(Kappa r)/r, C2 exp(-k r)/r, C3 exp(-c/r^2-c/(L-r)^2)*/
//#define  TRIAL_SQUARE_WELL_FS_PHONONS /* (3D) SS+SS(E>0) A sin(Kappa r)/r, B sin(k r + delta): Apar, Rpar */
//#define  TRIAL_SQUARE_WELL_BS_TRAP // (3D) sin(Kappa r)/r, C2 exp(-k r)/r, SW width RoSW
//#define  TRIAL_SQUARE_WELL_FS_TRAP // (3D) free state, a<0 (3D) SS sin(K r)/r [r<RoSW], 1+a/r [r<Rpar]
//#define  TRIAL_SQUARE_WELL_FS_TRAP_MCGUIRE1D
//#define  TRIAL_SQUARE_WELL_FS /* (3D) PBC SS+SS(E>0) A sin(Kappa r)/r, B sin(k r + delta) */
//#define  TRIAL_ZERO_RANGE_UNITARY /* (3D) BC absent: 1/r PBC: 1/r, r<Rpar, otherwise phonons exp(-c/r^2-c/(L-r)^2)*/ 
//#define  TRIAL_ZERO_RANGE_TRAP /* (3D) (r-a)/r, define HARD_SPHERE for f(r) = 0, r<a
//#define  TRIAL_ZERO_RANGE_3D /* (3D) sin(k*r+delta)
//#define  TRIAL_ZERO_RANGE_K0_TRAP // 2D Ko(Atrial*r), Atrial = 2e^{-gamma}/a
//#define  TRIAL_ZERO_RANGE_LOG_TRAP // 2D |ln(r/a)|
//#define  TRIAL_ZERO_RANGE_LOG_EXP_TAIL_TRAP // 2D ln(a/r) exp(-r^2/(b (Apar + r)))
//#define  TRIAL_ZERO_RANGE_K0_EXP_TAIL_TRAP // 2D Ko(Atrial*r) exp(-r^2/(b (Apar + r))), Atrial = 2e^{-gamma}/a
//#define  TRIAL_ZERO_RANGE_LOG_NO_NODE_TRAP // 2D Rpar: ln(1+a/r) * exp(-r/Rpar)*(1 + r/Rpar)
//#define  TRIAL_ZERO_RANGE_LOG_MATCH_TRAP // 2D Rpar: ln(r/a), r<Rpar; ln(r/a) * exp((r/Rpar-1) / ln(r/a))
//#define  TRIAL_SR3D /*VMC only! (3D) A/r |sin(kr + B)|, no variational parameters */ 
//#define  TRIAL_YUKAWA_3D_SYMM /* (3D) exp( - Apar exp(-Bpar r)/r (1 - exp(-r/Cpar)))*(1-exp(-pow((L/2-r),Btrial))) */
//#define  TRIAL_YUKAWA_3D /* exp( - Apar exp(-Bpar r)/r (1 - exp(-r/Cpar)))*/
//#define  TRIAL_CSM_3D_PLASMON // Rpar, lambda: r^lambda(1 + const r^2); exp(-const/r^3 - const / r^2)
//#define  TRIAL_SOFT_DISKS_2D_ZERO_ENERGY_PHONONS /* a>0, (2D) SD(E=0) parameter 0<Rpar<1: Io(kappa r), r<Ro; A ln(r/aA), r<Rpar; B exp(-C/r + D/r^2) */

//#define  THREE_BODY_TERMS
//#define  THREE_BODY_HYPERRADIUS_QUADRATIC // f_3(R): 0, R<R3; 1 - (R - Apar)^2 / (R3 - Apar)^2 , R3<R < Rpar; 1, R>Apar
//#define  THREE_BODY_ZERO_ENERGY_SOLUTION_3D // f_3(R): 0, R<R3; 1 - (R3/R)^4
//#define  THREE_BODY_ZERO_ENERGY_SOLUTION_2D // f_3(R): 0, R<R3; 1 - (R3/R)^2
//#define  THREE_BODY_ZERO_ENERGY_SOLUTION_1D // f_3(R): 0, R<R3; log(R/R3)
//#define  THREE_BODY_ZERO_ENERGY_SOLUTION_3D_SYMMETRIZED // f_3(R):  1 - (R3/R)^4 - (R3/(L-R))^4
//#define  THREE_BODY_ZERO_ENERGY_SOLUTION_1D_SYMMETRIZED // f_3(R): ln[R*(L-R)/(R3*(L-R3))], define HARD_SPHERE_HYPER RADIUS for f_3(R) = 0, R<R3
//#define  THREE_BODY_ZERO_ENERGY_SOLUTION_1D_PHONONS // (1D) f_3(R): 0, R<R3; ln(R/R3), r<Rpar; B exp(-C/r + D/r^2)
//#define  THREE_BODY_BOUND_STATE_SOLUTION_1D

/* define the phase: CRYSTAL (simulation of solid), LIQUID/GAS otherwise */

//#define CRYSTAL // solid phase add localizing terms to the trial w.f. and define one of the following options:
//#define CRYSTAL_NONSYMMETRIC // nonsymmetric w.f. for the crystal (symmetric otherwise)
//#define CRYSTAL_SYMMETRIC     // symmetric w.f. prod_{j,Nc} sum_{i=1,N} exp(-alpha_Rx[j](x_i-x_j^c)-alpha_Ry[j](y_i-y_j^c)-alpha_Rz[j](z_i-z_j^c))
//#define CRYSTAL_WIDTH_ARRAY // define to treat Crystal.Rx[i] Ry[i] Rz[i] w[i] as parameters loaded from crystal coordinates file (x,y,z,Rx,Ry,Rz,w)

//#define LATTICE_TRIANGULAR // (2D) square/cubic otherwise
//#define LATTICE_3D_BCC // 3D body centered cubic lattice
//#define LATTICE_3D_FCC // 3D face centered cubic lattice
//#define LATTICE_3D_CUBIC // default 3D lattice
//#define LATTICE_CHECKERBOARD // square (2D) checker board Ncrystal/2 lattice sites
//#define LATTICE_DOMINO // square (2D) domino Ncrystal/2 lattice sites
//#define LATTICE_ZIGZAG // zig-zag (2D), linear in X, zigzag in Y
//#define MEASURE_PROJECTION_PARTICLES_12 // measure w.f. overlap in 1D

//#define CENTER_OF_MASS_IS_NOT_MOVED // only VMC Move All and DMC in trap (no PBC)
//#define CENTER_OF_MASS_Z_IS_NOT_MOVED // only VMC Move All and DMC in trap (no PBC)
//#define CENTER_OF_MASS_IS_NOT_MOVED_TWO_PARTS  // CMseparation is added to {x,y,z} coordinates of the second part of the system
//use R2_subtract_CM=1 instead of #define CENTER_OF_MASS_SUBTRACTED_IN_RADIAL_DISTRIBUTION
//#define CHAINS_PINNED_TO_CENTER_OF_MASS // DMC or VMC (VMC only for move one by one)
//#define FIRST_PARTICLE_IS_NOT_MOVED // pin position of 1st particlr
//#define PAIR_DISTRBUTION_EMPTINESS_FORMATION_PROBABILITY_1D // calculate probability of a gap of size x between in particles in 1D

// select dimensionality
#define TRIAL_1D
//#define TRIAL_2D
//#define TRIAL_3D

/* define ONLY ONE of following boundary conditions */
#define  BC_ABSENT /* free boundary condition */
//#define  BC_1DPBC_Z   /* boundary condition in Z direction */
//#define  BC_1DPBC_X   /* boundary condition in X direction */
//#define  BC_2DPBC_SQUARE     /* boundary condition in (X,Y) plane */
//#define  BC_2DPBC_HEXAGON
//#define  BC_2DPBC_NON_ORTHOGONAL /**/
//#define  BC_3DPBC_CUBE   /* PBC in all directions */
//#define  BC_3DPBC_TRUNCATED_OCTAHEDRON /* overlap between cube and |x-L/2|+|y-L/2|+|z-L/2| < 3L/4 */
//#define  BC_SPHERE   /* sphere of size Rsphere */

/* define if trap potential is present */
#define TRAP_POTENTIAL

/* if trap one-body terms are present in the trial w.f., enable ONE_BODY_TRIAL_TERMS*/
//#define ONE_BODY_TRIAL_TERMS
//#define ONE_BODY_TRIAL_TRAP // alpha_x * x*x + alpha_y * y*y + alpha_z * z*z;
//#define ONE_BODY_TRIAL_SIN // (1D,2D,3D) lattice: one-body term: (3D,2D)[sin(kL*x)^2 + sin(kL*y)^2 + sin(kL*z)^2]^alpha_R with kL = PI*Nlattice/L (1D) beta_latt + |sin(kL*z)|^alpha_latt
//#define ONE_BODY_TRIAL_ZERO_BOUNDARY_CONDITION // one-body term: |sin(pi/L*z)|^beta_z, default value beta_z=1, see also ONE_BODY_TRIAL_ABC, also can be used together with ONE_BODY_TRIAL_SIN
//#define ONE_BODY_TRIAL_ZBC // sin(pi*z/L)
//#define ONE_BODY_TRIAL_ABC // antiperiodic boundary condition: |sin(pi x /L)|^alpha_latt
//#define ONE_BODY_TRIAL_1D_WAVE_COS // down particles: |cos(one_body_momentum_i pi z/L)|, where one_body_momentum_i = 0,1,2,...
//#define ONE_BODY_TRIAL_1D_WAVE_SIN // down particles: |sin(one_body_momentum_i pi z/L)|, where one_body_momentum_i = 0,1,2,...
//#define ONE_BODY_TRIAL_COS_SERIES // one-body term (1D) 1 + alpha_latt*cos(2*k*z) + beta_latt*cos(4*k*z)+ gamma_latt*cos(6*k*z)
//#define ONE_BODY_TRIAL_1D_WAVE_COS_FINITE // 1 + S cos(one_body_momentum_i 2 pi z/L)
//#define ONE_BODY_SOLITON
//#define ONE_BODY_TRIAL_TRAP_POWER_LAW  // exp(-alpha_latt *r*r) * r^beta_latt
//#define ONE_BODY_TRIAL_LINEAR_CHAIN_UP
//#define ONE_BODY_IMPURITY // define and chose one of the following
//#define ONE_BODY_IMPURITY_1D
//#define ONE_BODY_IMPURITY_3DSW // Square Well
//#define ONE_BODY_TRIAL_DOUBLE_GAUSSIAN // exp(-alpha_z * (z-R0par)^2) + exp(-alpha_z * (z+R0par)^2)
//#define VEXT_IMPURITY_3DSW
//#define ONE_BODY_TRIAL_LATTICE // beta_latt + gamma_latt * sum exp(-alpha_latt*(r-r_i)^2), i=1,...,Ncrystal
//#define ONE_BODY_TRIAL_TRAP_DOWN_ANTISYMMETRIC_Z // (HO,n=1) multiplies one-body wave function by z_dn, so that there is a node at z_dn = 0
//#define ONE_BODY_TRIAL_TRAP_DOWN_HERMITE2_Z // (HO,n=2) multiplies one-body wave function by z_dn^2 - beta, default value beta=1/2, so that there is a node at |z_dn| = 1/sqrt(2)
//#define MANY_BODY_TRIAL_TRAP_CENTER_OF_MASS_NODE_Z // (z1+...+zN)
//#define MANY_BODY_TRIAL_TRAP_CENTER_OF_MASS_HERMITE2_Z // sum (z_i/N)^2 - beta

/* alternatively use splines */
//#define INTERPOLATE_SPLINE_ONE_BODY_Z_WF // use spline interpolation for 1-body terms in trial w.f.

/* define external confinement in a homogeneous system, define EXTERNAL_POTENTIAL */
//#define VEXT_COS2 // cos^2 external lattice potential
//#define VEXT_SOLITON_FIXED_PHASE
//#define VEXT_SPLINE_1D // use spline interpolation for Vext(z)
//#define EXTERNAL_POTENTIAL_Uo_R2
//#define EXTERNAL_POTENTIAL_CONSTANT_SPINFULL // D : potential energy of spin 0 
//#define VEXT_GAUSSIAN // gaussian_alpha * exp( -gaussian_beta *x*x)
//#define VEXT_TRAP_DOUBLE_WELL // (0.5*z*z + alpha_latt*Exp(-0.5*z*z) - (Log(alpha_latt)+1))

//#define SPECKLES // define for applying speckle potential

// Choose one type of interpolation method
//#define SPLINES_2D // 2D interpolation: use 2D splines 
//#define POINT9 // 2D interpolation: use interpolation between nearest 9 points
#define POINT16 // 2D interpolation: use interpolation between nearest 16 points
//#define POINT16CONTINUOUS // 2D interpolation: interp between nearest 16 points, 1st deriv continous
//#define POINT16CONTINUOUS_2deriv // 2D interpolation: interp between nearest 16 points, 1st, 2nd deriv continous enable togethe with POINT16CONTINUOUS
//#define POINT4 // 2D interpolation: 4 points and 4 2nd derivatives

#define energy_unit 1.
//#define energy_unit 16.079147954821714 // He3 units, i.e. [hbar^2 / (mHe3 A^2)] / (kB K) = 16.
//#define energy_unit 12.115916747700853 // He4 units
//#define energy_unit 12.485545424687677 // He critical mass for two-body bound state, 12.4855 u
//#define energy_unit 48.11317400458453 // H units, i.e. [hbar^2 / (mH A^2)] / (kB K)
//#define energy_unit 24.077826199121347 // Deuterium units
//#define energy_unit 16.07904559586791 // Tritium units
//#define energy_unit 15.64361051812288  // hydrogen, meff = 3.1u
//#define energy_unit 15.154747689431542 // hydrogen, meff = 3.2u
//#define energy_unit 14.93524884938404 // hydrogen, meff = 3.247u threshold
//#define energy_unit 14.921597724978747 // hydrogen, meff = 3.25u
//#define energy_unit 14.695512910963922 // hydrogen, meff = 3.3u
//#define energy_unit 14.263291942994396 // hydrogen, meff = 3.4u
//#define energy_unit 13.855769316051697 // hydrogen, meff = 3.5u
//#define energy_unit 12.123798151545234 // hydrogen, meff = 4u
//#define energy_unit 9.699038521236188 // hydrogen, meff = 5u

/* disable checks of the continuity of the trial wf.*/
//#define DONT_CHECK_CONTINUITY_OF_WF
//#define DONT_CHECK_CONTINUITY_OF_WF_AT_L_2
#define DONT_CHECK_CONTINUITY_OF_ELOC

// by default SD is measured in all accessible directions, unless the following is selected
//#define SD_MEASURE_ONLY_X_DIRECTION
//#define SD_MEASURE_ONLY_Y_DIRECTION
//#define SD_MEASURE_ONLY_Z_DIRECTION

//#define UNITS_SET_2M_TO_1 // unit of energy is h2/2ma2 instead of h2/ma2
//#define RESCALE_TO_RADIAL_OSC_UNITS // units of h w_perp (1D trapped system)
//#define RESCALE_COORDINATES_TO_OSC_UNITS

#ifdef MPI
//#  define EXEC_NAME "/home-kirk/astra/VMC/DMC"
//#  define INPATH "/disc1/users/astra/MC/"
//#  define INPATH  "/home/astra/MC/"
//#  define OUTPATH "/home/astra/MC/"
#  define INPATH  "./"
#  define OUTPATH "./"
#else
#  define INPATH  "./"
#  define OUTPATH "./"
#endif

#ifdef TRIAL_TONKS
#  define DONT_CHECK_CONTINUITY_OF_WF
#endif

#ifdef TRIAL_LIEB
#  define DONT_CHECK_CONTINUITY_OF_ELOC
#endif

#ifdef TRIAL_LIEB_LUTTINGER_PHONON
#  define DONT_CHECK_CONTINUITY_OF_ELOC
#endif

#ifdef TRIAL_SCHIFF_VERLET_TRAP
#  define DONT_CHECK_CONTINUITY_OF_WF
#  define DONT_CHECK_CONTINUITY_OF_ELOC
#endif
#ifdef TRIAL_SCHIFF_VERLET_TRAP_PHONON2D
#  define DONT_CHECK_CONTINUITY_OF_WF
#  define DONT_CHECK_CONTINUITY_OF_ELOC
#endif

#ifdef BC_SPHERE
#  define DONT_CHECK_CONTINUITY_OF_WF_AT_L_2
#endif

#ifdef BC_ABSENT
#  define boundary NO_BOUNDARY_CONDITIONS
//#  define EXTERNAL_POTENTIAL
//#  define BC_TRAP
#  define DONT_CHECK_CONTINUITY_OF_WF_AT_L_2
#endif

#ifdef TRIAL_1D
#  define DONT_CHECK_CONTINUITY_OF_WF_AT_L_2
#endif

#ifdef BC_TRAP
#  define EXTERNAL_POTENTIAL
#  define BC_TRAP
#  define ONE_BODY_TRIAL_TERMS
#  define SAVE_RADIAL_DISTR 1
//#ifndef INTERPOLATE_SPLINE_ONE_BODY_Z_WF // 1-body terms are described by splines
//#    define ONE_BODY_TRIAL_TRAP
//#endif
#else
#  define SAVE_RADIAL_DISTR 0
#endif

#ifdef TRAP_POTENTIAL 
#  define EXTERNAL_POTENTIAL
#endif

#ifdef SPECKLES
#  define EXTERNAL_POTENTIAL
#endif

#ifdef VEXT_SOLITON_FIXED_PHASE
#  define EXTERNAL_POTENTIAL
#endif

#ifdef VEXT_IMPURITY_3DSW
#  define EXTERNAL_POTENTIAL
#endif

#ifdef BC_1DPBC_Z
#  define boundary ONE_BOUNDARY_CONDITION
#  define BC_1DPBC_ANY
#endif

#ifdef BC_1DPBC_X
#  define boundary ONE_BOUNDARY_CONDITION
#  define BC_1DPBC_ANY
#endif

#ifdef BC_2DPBC_SQUARE
#  define BC_2DPBC
#  define boundary TWO_BOUNDARY_CONDITIONS
#endif

#ifdef BC_2DPBC_HEXAGON
#  define BC_2DPBC
#  define boundary TWO_BOUNDARY_CONDITIONS
#endif

#ifdef BC_3DPBC_CUBE
#  define boundary THREE_BOUNDARY_CONDITIONS
#  define BC_3DPBC
#endif

#ifdef BC_3DPBC_TRUNCATED_OCTAHEDRON
#  define boundary THREE_BOUNDARY_CONDITIONS
#  define BC_3DPBC
#endif

#ifdef BC_SPHERE
#  define boundary NO_BOUNDARY_CONDITIONS 
#endif

#ifdef LATTICE 
#  define EXTERNAL_POTENTIAL
#  //define CRYSTAL
#  undef LATTICE_TRIANGULAR
#endif

#ifdef VEXT_COS2
#  define EXTERNAL_POTENTIAL
#endif

#ifdef VEXT_SPLINE_1D
#  define EXTERNAL_POTENTIAL
#endif

#ifdef EXTERNAL_POTENTIAL_Uo_R2
#  define EXTERNAL_POTENTIAL
#endif

#ifdef VEXT_TRAP_DOUBLE_WELL
#  define EXTERNAL_POTENTIAL
#endif

#ifdef TRIAL_SUTHERLAND
#  define DONT_CHECK_CONTINUITY_OF_WF
#  define DONT_CHECK_CONTINUITY_OF_ELOC
#endif
#ifdef TRIAL_CALOGERO
#  define DONT_CHECK_CONTINUITY_OF_WF
#  define DONT_CHECK_CONTINUITY_OF_ELOC
#endif

#ifdef TRIAL_PHONON_NEGATIVE
#  define TRIAL_PHONON
#endif

#ifdef CRYSTAL
#  define ONE_BODY_TRIAL_TERMS
#endif

#ifdef ONE_BODY_IMPURITY
#  define ONE_BODY_TRIAL_TERMS
#endif

#ifdef ONE_BODY_TRIAL_SIN
#  define ONE_BODY_TRIAL_TERMS
#endif

#ifdef ONE_BODY_TRIAL_COS_SERIES
#  define ONE_BODY_TRIAL_TERMS
#endif

#ifdef ONE_BODY_SOLITON
#  define ONE_BODY_TRIAL_TERMS
#endif

#ifdef ONE_BODY_TRIAL_LATTICE
#  define ONE_BODY_TRIAL_TERMS
#endif

#ifdef BC_ABSENT
//#  define ONE_BODY_TRIAL_TERMS
//#ifndef INTERPOLATE_SPLINE_ONE_BODY_Z_WF // 1-body terms are described by splines
//#    define ONE_BODY_TRIAL_TRAP
//#endif
#endif

#ifdef CRYSTAL
#  ifndef CRYSTAL_NONSYMMETRIC
#    ifndef CRYSTAL_SYMMETRIC
//#        error "specify symmetry of the crystal"
#       define CRYSTAL_NONSYMMETRIC
#    endif
#  endif
#endif

#ifdef ONE_BODY_TRIAL_TRAP
#  define ONE_BODY_TRIAL_TERMS
#endif

#ifdef ONE_BODY_TRIAL_ZBC
#  define ONE_BODY_TRIAL_TERMS
#endif

#ifdef ONE_BODY_TRIAL_ABC
#  define ONE_BODY_TRIAL_TERMS
#endif

#ifdef ONE_BODY_TRIAL_1D_WAVE
#  define ONE_BODY_TRIAL_TERMS
#endif

#ifdef ONE_BODY_TRIAL_ZERO_BOUNDARY_CONDITION
#  define ONE_BODY_TRIAL_TERMS
#endif

#ifdef SPECKLES
#  define ONE_BODY_TRIAL_TERMS
#endif

#ifdef INTERPOLATE_SPLINE_ONE_BODY_Z_WF
#  define ONE_BODY_TRIAL_TERMS
#endif

#ifdef ONE_BODY_TRIAL_TRAP_POWER_LAW
#  define ONE_BODY_TRIAL_TERMS
#endif

// calculate OBDM and nk for fermions
#ifdef TRIAL_SUTHERLAND
#  define OBDM_FERMIONS
#endif
#ifdef TRIAL_CALOGERO
#  define OBDM_FERMIONS
#endif
#ifdef TRIAL_COULOMB_1D_PHONON
#  define OBDM_FERMIONS
#endif
#ifdef TRIAL_DIPOLE_1D_PHONON
#  define OBDM_FERMIONS
#endif
#ifdef INTERACTION_LINEAR_1D_COULOMB
#  define OBDM_FERMIONS
#endif
#ifdef TRIAL_SIN_EXCITED_BOSONS
#  define OBDM_FERMIONS
#  define OBDM_EXCITED_BOSONS // take into account a possible change in the sign of the w.f.
#endif
#ifdef TRIAL_SIN_EXCITED_BOSONS_KPAR
#  define OBDM_FERMIONS
#  define OBDM_EXCITED_BOSONS // take into account a possible change in the sign of the w.f.
#endif
#ifdef TRIAL_GAPPED
#  define OBDM_FERMIONS
#endif
#ifdef TRIAL_TONKS_TRAP // p-wave fermions
#  define OBDM_FERMIONS
#endif
#ifdef TRIAL_MCGUIRE // p-wave fermions
#  define OBDM_FERMIONS
#endif
#ifdef TRIAL_LIEB
#  define OBDM_FERMIONS
#endif
#ifdef OBDM_ANYONS
#  define OBDM_FERMIONS
#endif
#ifdef TRIAL_TONKS
#  define OBDM_FERMIONS
#endif

#define VARIATIONAL 0
#define DIFFUSION 1
#define CLASSICAL 2
#define PIMC 3
#define PIGS 4

#define NO_BOUNDARY_CONDITIONS 0
#define THREE_BOUNDARY_CONDITIONS 1
#define ONE_BOUNDARY_CONDITION 2
#define TWO_BOUNDARY_CONDITIONS 3

// SmartMC constants
#define DMC_QUADRATIC 0
#define DMC_LINEAR_DRIFT 1
#define DMC_LINEAR_GAUSSIAN 2
#define DMC_LINEAR_METROPOLIS 3
#define DMC_QUADRATIC_METROPOLIS 4
#define DMC_QUARTIC 5
#define DMC_PSEUDOPOTENTIAL 6
#define DMC_MOVE_ONE 7

#define VMC_MOVE_ALL 0
#define VMC_MOVE_ONE 1
#define VMC_MOVE_DRIFT_ALL 2
#define VMC_MOVE_DRIFT_ONE 3
#define VMC_SUM_OF_TWO_WAVEFUNCTIONS 4

#define PIGS_PRIMITIVE 0
#define PIGS_PSEUDOPOTENTIAL 1
#define PIGS_PSEUDOPOTENTIAL_OBDM 2

#ifdef TRIAL_1D
#  define DIMENSION 1
#endif
#ifdef TRIAL_2D
#  define DIMENSION 2
#endif
#ifdef TRIAL_3D
#  define DIMENSION 3
#endif

#ifdef TRIAL_EXTERNAL
#  define SUBSTITUTE_INTERPOLATED_WF
#endif

//#ifdef INTERPOLATE_LINEAR_JASTROW_WF
//#  define SUBSTITUTE_INTERPOLATED_WF
//#endif

#ifdef INTERPOLATE_SPLINE_JASTROW_WF
#define INTERPOLATE_TYPE_SPLINE
#endif

#ifdef INTERPOLATE_SPLINE_ONE_BODY_Z_WF 
#define INTERPOLATE_TYPE_SPLINE
#endif
#include <time.h>

#define RAN2_SEED 1023
int AllocateMemory(void);

#ifdef CRYSTAL
#  include "crystal.h"
#endif

#ifdef DATATYPE_FLOAT
typedef float DOUBLE;
#  define LF "f"
#  define LE "e"
#  define LG "g"
#  define Cos cosf
#  define Sin sinf
#  define Tg tanf
#  define Arctg atanf
#  define Log logf
#  define Exp expf
#  define Sqrt sqrtf
#endif

#ifdef DATATYPE_DOUBLE
//typedef double DOUBLE;
#  define DOUBLE double
#  define LF "lf"
#  define LE "le"
#  define LG "lg"
#  define Cos cos
#  define Sin sin
#  define Tg tan
#  define Arctg atan
#  define Log log
#  define Exp exp
#  define Sqrt sqrt
#endif

#ifdef DATATYPE_LONG_DOUBLE
#  define DOUBLE long double
#  define LF "Lf"
#  define LE "Le"
#  define LG "Lg"
#  define Cos cosl
#  define Sin sinl
#  define Tg tanl
#  define Arctg atanl
#  define Log logl
#  define Exp expl
#  define Sqrt sqrtl
#endif

#ifdef CRYSTAL_WIDTH_ARRAY
#  define Crystal_dot_Rx_i Crystal.Rx[i]
#  define Crystal_dot_Ry_i Crystal.Ry[i]
#  define Crystal_dot_Rz_i Crystal.Rz[i]
#  define Crystal_dot_weight_i Crystal.weight[i]
#  define Crystal_dot_Rx_j Crystal.Rx[j]
#  define Crystal_dot_Ry_j Crystal.Ry[j]
#  define Crystal_dot_Rz_j Crystal.Rz[j]
#  define Crystal_dot_weight_j Crystal.weight[j]
#else
#  define Crystal_dot_Rx_i Crystal.Rx
#  define Crystal_dot_Ry_i Crystal.Ry
#  define Crystal_dot_Rz_i Crystal.Rz
#  define Crystal_dot_weight_i 1.
#  define Crystal_dot_Rx_j Crystal.Rx
#  define Crystal_dot_Ry_j Crystal.Ry
#  define Crystal_dot_Rz_j Crystal.Rz
#  define Crystal_dot_weight_j 1.
#endif

#ifndef INTERACTION_SUM_OVER_IMAGES
#  ifndef INTERACTION_EWALD
#    ifndef INTERACTION_SUM_OVER_IMAGES_SMOOTH_CUTOFF
#      define INTERACTION_WITHOUT_IMAGES
#    endif
#  endif
#endif

#ifdef MPI
#  define MEMORY_CONTIGUOUS
#endif

#ifdef TRIAL_3D
#  define MOVE_IN_X
#  define MOVE_IN_Y
#  define MOVE_IN_Z
#endif

#ifdef TRIAL_2D
#  define MOVE_IN_X
#  define MOVE_IN_Y
#endif

#ifdef TRIAL_1D
#ifdef BC_1DPBC_X // in 1D trap move in z
#  define MOVE_IN_X
#else
#  define MOVE_IN_Z
#endif
#endif

#ifdef MOVE_IN_X
# define CaseX(arg) arg
#else
# define CaseX(arg)
#endif

#ifdef MOVE_IN_Y
# define CaseY(arg) arg
#else
# define CaseY(arg)
#endif

#ifdef MOVE_IN_Z
# define CaseZ(arg) arg
#else
# define CaseZ(arg)
#endif

#ifdef SD_MEASURE_ONLY_X_DIRECTION
#  define CaseXSD(arg) arg
#  define CaseYSD(arg)
#  define CaseZSD(arg)
#else
#ifdef SD_MEASURE_ONLY_Y_DIRECTION
#  define CaseXSD(arg)
#  define CaseYSD(arg) arg
#  define CaseZSD(arg)
#else
#ifdef MOVE_IN_X
#  define CaseXSD(arg) arg
#else
#  define CaseXSD(arg)
#endif
#ifdef MOVE_IN_Y
#  define CaseYSD(arg) arg
#else
#  define CaseYSD(arg)
#endif
#ifdef MOVE_IN_Z
#  define CaseZSD(arg) arg
#else
#  define CaseZSD(arg)
#endif
#endif
#endif

#ifdef SD_MEASURE_ONLY_Z_DIRECTION
#  define CaseXSD(arg)
#  define CaseYSD(arg)
#  define CaseZSD(arg) arg
#endif

#ifdef SPECKLES // spline interpolation fail if any coordinate is outside of the array limits
# define PARTICLES_NEVER_LEAVE_THE_BOX
#endif

#ifdef TRIAL_SCHIFF_VERLET_TRAP
# define TRIAL_SCHIFF_VERLET
#endif
#ifdef TRIAL_SCHIFF_VERLET_PHONON
# define TRIAL_SCHIFF_VERLET
#endif
#ifdef TRIAL_SCHIFF_VERLET_TRAP_PHONON2D
# define TRIAL_SCHIFF_VERLET
#endif

#ifdef TRIAL_3D // quasi 2D geometry
#ifdef BC_2DPBC
#define BC_Q2D
#endif
#endif

#ifndef BC_2DPBC_HEXAGON 
#  define BOX_EQUAL_SIDES // set Lx=Ly=Lz
#endif

#ifndef LATTICE_3D_BCC
#  ifndef LATTICE_3D_FCC
#    define LATTICE_3D_CUBIC
#  endif
#endif

#ifdef BC_1DPBC_Z
#  define BC_1DPBC
#endif
#ifdef BC_1DPBC_X
#  define BC_1DPBC
#endif
#ifdef BC_2DPBC_SQUARE
#  define BC_2DPBC
#endif
#ifdef BC_2DPBC_HEXAGON
#  define BC_2DPBC
#endif
#ifdef BC_2DPBC_NON_ORTHOGONAL
#  define BC_2DPBC
#endif
#ifdef BC_3DPBC_CUBE
#  define BC_3DPBC
#endif
#ifdef BC_3DPBC_TRUNCATED_OCTAHEDRON
#  define BC_3DPBC
#endif

#ifndef TRIAL_1D
#  ifndef TRIAL_2D
#    ifndef TRIAL_3D
#      error "define TRIAL_XD"
#    endif
#  endif
#endif

#ifndef CALCULATE_DERIVATIVES_NUMERICALLY
#  define CALCULATE_DERIVATIVES_ANALYTICALLY
#endif

#ifndef SPINFULL
#  define SPINLESS
#endif

//#ifdef SPINFULL
//extern struct Distribution **PDSpin; // PDSpin[spin1 - spin2][grid_PD] array
//#endif

#ifdef TRIAL22_SAME_AS_11
#  define InterpolateExactU22 InterpolateExactU11
#  define InterpolateExactFp22 InterpolateExactFp11
#  define InterpolateExactE22 InterpolateExactE11
#endif

#ifdef TRIAL12_SAME_AS_11
#  define InterpolateExactU InterpolateExactU11
#  define InterpolateExactFp InterpolateExactFp11
#  define InterpolateExactE InterpolateExactE11
#endif

#ifdef TRIAL11_SAME_AS_12
#  define InterpolateExactU11 InterpolateExactU
#  define InterpolateExactFp11 InterpolateExactFp
#  define InterpolateExactE11 InterpolateExactE
#endif

#ifdef TRIAL22_SAME_AS_12
#  define InterpolateExactU22 InterpolateExactU
#  define InterpolateExactFp22 InterpolateExactFp
#  define InterpolateExactE22 InterpolateExactE
#endif


#ifdef MEASURE_PURE_OBDM
#  ifdef MEASURE_PURE_TBDM
#    error "define either MEASURE_PURE_OBDM or MEASURE_PURE_TBDM, but not both of them"
#  endif
#endif

#ifdef MEASURE_PURE_OBDM
#  define VIRTUAL_WALKERS
#endif

#ifdef MEASURE_PURE_TBDM
#  define VIRTUAL_WALKERS
#endif

#ifdef TRIAL_TWO_BODY_SCATTERING_NUMERICAL_SPINFULL
  extern struct Grid **GridSpin; // [spin1][spin2] array
#undef INTERPOLATE_SPLINE_JASTROW_WF
#endif

#ifdef INTERACTION_SPINFULL_GAUSSIAN
#  define SAME_SPIN_CONSTANT_JASTROW
#endif

#ifdef INTERACTION_SPINFULL_DOUBLE_GAUSSIAN
#  define SAME_SPIN_CONSTANT_JASTROW
#endif

#ifdef JASTROW_RIGHT_BOUNDARY_SPINFULL_SAME_SPIN_BOUND_STATE
#  define JASTROW_RIGHT_BOUNDARY_BOUND_STATE
#endif

#ifdef JASTROW_RIGHT_BOUNDARY_SPINFULL_DIFFERENT_SPIN_BOUND_STATE
#  define JASTROW_RIGHT_BOUNDARY_BOUND_STATE
#endif

#define JASTROW_GOES_TO_ONE_AT_LARGEST_DISTANCE

#ifdef JASTROW_RIGHT_BOUNDARY_PHONONS
#  undef JASTROW_GOES_TO_ONE_AT_LARGEST_DISTANCE
#endif

#ifdef INTERACTION_LINEAR_1D_COULOMB
#  define DONT_CHECK_CONTINUITY_OF_WF
#endif

#ifndef JASTROW_LEFT_BOUNDARY_ZERO
#  ifndef JASTROW_LEFT_BOUNDARY_PSEUDOPOTENTIAL
#    define  JASTROW_LEFT_BOUNDARY_ZERO_DERIVATIVE
#  endif
#endif

struct Walker {
  /* coordinates */
#ifdef NPARTICLES
  DOUBLE x[NPARTICLES];
  DOUBLE y[NPARTICLES];
  DOUBLE z[NPARTICLES];
  int spin[NPARTICLES];
#ifdef SPINFULL
  DOUBLE psiT_sigma_inv[NPARTICLES]; // psi_T(R,-sigma_i) / psi_T(R,sigma_i)
#endif
#else
  DOUBLE *x;
  DOUBLE *y;
  DOUBLE *z;
  int *spin;
#ifdef SPINFULL
  DOUBLE *psiT_sigma_inv; // psi_T(R,-sigma_i) / psi_T(R,sigma_i)
#endif
#endif
  DOUBLE E; 
  DOUBLE EFF; 
  DOUBLE Ekin; 
  DOUBLE Epot; 
  DOUBLE Edamping; 
  DOUBLE Eint;
  DOUBLE Eext;
#ifdef CHAINS_PINNED_TO_CENTER_OF_MASS
  DOUBLE *Crystal_x;
  DOUBLE *Crystal_y;
  DOUBLE *Crystal_z;
#endif  
#ifdef VIRTUAL_WALKERS // Pure OBDM or TBDM
  DOUBLE OBDMpureB; // bosons
  DOUBLE OBDMpureF; // fermions
  DOUBLE OBDMpure_r;  // displacement in r for nk calculation
  DOUBLE *OBDMweight; 
  int OBDM_position;
#endif
#ifdef CRYSTAL_SYMMETRIC
  DOUBLE *Mo; // needed for the symm. crystal w.f. 
#endif
  DOUBLE Eold; 
  DOUBLE U; // f=exp(U)
  int status; //0-alive, 1-dead, 2-reincarnation, 3-killed, 4-virtual(OBDM)
#ifdef MEASURE_CORRELATION_FUNCTIONS_IN_THE_MIDDLE_OF_DRIFT
  int status_end_of_step_initialized;
#endif
  int status_sum_of_two_walkers_U;
  int status_sum_of_two_walkers_E;
  DOUBLE weight;
  int w; // Walker index

  DOUBLE r2, z2, r4, z4, r2pd; // storing information of the current block
  DOUBLE r2old, z2old; // storing information of the previous block

  int times_measured; // needed pure estimators, for the arrays to be filly filled

  int **PD;   // Pure pair distribution [grid_pure_block][gridPD]
  DOUBLE **rhoRe;  // form factor measurement [gridSKT_t][gridSKT_k]
  DOUBLE **rhoIm;  // form factor measurement [gridSKT_t][gridSKT_k]
  DOUBLE **psi;  // one-body density matrix <psi(tau)psi(0)> [gridSKT_t][gridSKT_k]
  int PD_position;

  int **RD;   // Pure radial distribution [grid_pure_block][gridRD]
  int **RDz;  // Pure radial distribution [grid_pure_block][gridRD]
  DOUBLE **HR;   // 3body correlation function in hyperradius [grid_pure_block][gridg3]
  DOUBLE *R2pure; // pure r2
  DOUBLE *Z2pure; // pure r2
  DOUBLE *PD_R2pure; // pure r2 calculated over pair distribution
  int RD_wait;
  int RD_position;
  int HR_position;
  int g3_store; // g3(0)

  DOUBLE **Sk; // Pure static structure factor
  unsigned long int **Sk_N; // Pure static structure factor
  DOUBLE *Skx, *Sky; // check for a crystal structure
  int Sk_position;
  int Sk_tau_position;
  int Sk_count;
  int Skxy_count;

  DOUBLE *LindF;   // Pure Lindemann ratio [grid_pure_block]
  int *LindN;   // Pure Lindemann ratio [grid_pure_block]
  int Lind_position;

  // calculation of the superfluid density
  DOUBLE **CM; // [(x,y,z)][SD.size] center of mass position
  DOUBLE **rreal; // [N][3]

  DOUBLE *Epot_pure; // Pure potential energy [grid_pure_block]
  DOUBLE *x_pure; // Pure coordinates [N]
  DOUBLE *y_pure; // Pure coordinates [N]
  DOUBLE *z_pure; // Pure coordinates [N]
  
  // scalable code
  int *c_Nlocal; // [Nc] number of particles in the current cell only
  int **c_index_local; // [Nc][N_local] indices of particles belonging to the current cell
  int *c_Nall; // [Nc] number of particles in the current cell and in the nearest cells
  int **c_index_all; // [Nc][N_all] indices of particles belonging to the current and nearest cells
  int *c_Npositive; // [Nc] same as cN_all, but avoids DOUBLE counting
  int **c_index_positive; // [Nc][N_positive] same as c_index_all, but avoids DOUBLE counting

  // DMC or smart VMC drift variables 
  DOUBLE **F, **Fp, **Fpp; /*FF[N][3]*/
  DOUBLE **R, **Rp, **Rpp, **Rppp, **R_dmc_end_of_step;
  //DOUBLE **Fold, **Rold, **dR_drift_old, **dR_drift, **dR_gauss; // for Metropolis quadratic move
  DOUBLE **dR_drift, **dR_drift_old, **dR, **dR_gauss; // for Metropolis quadratic move

#ifdef SPINFULL
  int ***PDSpin;
#endif 
};

#define ALIVE 0
#define DEAD 1
#define REINCARNATION 2
#define KILLED 3
#define VIRTUAL 4

#define ONE_COMPONENT_CODE

/* define Walker structure BEFORE including "quantities.h" */
#include "quantities.h"

extern int MC, SmartMC;
extern struct Walker *W, *Wp;
extern struct Grid G;  // Jastrow
extern struct Grid G1; // One-body
extern struct Grid GV; // Vext

extern struct sOBDM OBDM, OBDMtrap, OBDMfermi, TBDM, TBDMfermi;
extern struct sMATRIX OBDM_MATRIX, TBDM_MATRIX, OBDMfermi_MATRIX, TBDMfermi_MATRIX, PD_MATRIX, RD_MATRIX, SKT, SKTre, SKTim, PhiTau, Nk_MATRIX, Nkfermi_MATRIX;
extern struct Distribution PD, PDz, PD_pure, RD, RDz, RD_pure, RDz_pure, HD, HR, HR_pure, Veff;
extern struct sSk Sk, Sk_pure;
extern struct sSD SD;
extern struct sOrderParameter OrderParameter;
extern struct sLindemannRatio LindemannRatio;
extern struct sPureEstimator Epot_pure;

extern DOUBLE L, Lhalf, L2, Lhalf2, Lmax, Lwf, Lhalfwf, Lhalfwf2, Lcutoff_pot, Lcutoff_pot2;
extern DOUBLE Lx, Ly, Lz, L_half_x, L_half_y, L_half_z, L_inv_x, L_inv_y, L_inv_z, kL;
extern int grid_trial;
extern int gridOBDM;
extern int gridOBDM_MATRIX;
extern int gridRD;
extern int gridPD;
extern int gridg3;
extern int gridRDx;
extern int gridRDy;
extern int grid_pure_block;
extern int gridRD_pure_block;
extern int gridSk;
extern int gridPD_MATRIX_x;
extern int gridPD_MATRIX_y;
extern int gridSKT_t, gridSKT_k;
extern DOUBLE SKT_dk, SKT_kmin;
extern DOUBLE PD_MATRIX_x;
extern DOUBLE PD_MATRIX_y;
extern unsigned long int Niter, Nmeasure, iteration_global;
extern DOUBLE Rpar, Apar, Bpar, Cpar, Dpar, R0par, Epar, Ipar, Kpar, Mpar, alpha_R, alpha_latt, alpha_Jastrow, beta_latt, gamma_latt, beta_z;
extern DOUBLE Rpar11, Rpar12, Apar12, Bpar12, Cpar12;
extern DOUBLE alpha_x, alpha_x2, two_alpha_x, alpha_Rx;
extern DOUBLE alpha_y, alpha_y2, two_alpha_y, alpha_Ry;
extern DOUBLE alpha_z, alpha_z2, two_alpha_z, alpha_Rz;
extern DOUBLE omega_x, omega_y, omega_z;
extern DOUBLE omega_x2, omega_y2, omega_z2;
extern DOUBLE n, a, a2, b, lambda, lambda2, lambda4;

extern DOUBLE dt, dt_vmc, dt_one, dt_all, acceptance_rate, beta;
extern int N, Nimp, Ndens, NCrystal, Nspin;
extern int Nwalkers, NwalkersMax, Nwalkers_pure;
extern DOUBLE Nwalkersw, Nwalkers_weight_killed;
extern int Npop, Npop_max, Npop_min, Npop_virtual;
extern int blck, blck_heating;
extern int McMillan_points;
extern int McMillanTBDM_points;

extern char file_particles[], file_energy[], file_wf[], file_wf1b[];
extern char file_OBDM[], file_OBDM_MATRIX[],  file_PDz[], file_RD[], file_RDz[];
extern char file_R2[], file_z2[], file_R2_pure[], file_z2_pure[];
extern char file_PD[], file_PD_pure[], file_RD_pure[], file_RDz_pure[];
extern char file_Sk[], file_Sk_pure[], file_SD[];

extern int measure_energy;
extern int measure_OBDM;
extern int measure_TBDM;
extern int measure_OBDM_MATRIX;
extern int measure_SD;
extern int measure_RadDistr;
extern int measure_RadDistr_pure;
extern int measure_RadDistrMATRIX;
extern int measure_PairDistr;
extern int measure_PairDistrMATRIX;
extern int measure_PairDistr_R2;
extern int measure_R2;
extern int measure_Sk;
extern int measure_OP;
extern int measure_RDz;
extern int measure_g3;
extern int measure_effective_potential;

extern int video;
extern int branchng_present;

extern int verbosity;
extern int accepted, rejected;
extern int accepted_one, rejected_one;
extern int accepted_all, rejected_all;
extern int overlaped;
extern int generate_new_coordinates, generate_new_coordinates_from_lattice;

extern DOUBLE *u_mi, *u_ij;
extern DOUBLE *u_tbdm_i, **u_tbdm_ij, *u_tbdm_mcmillan_i, **u_tbdm_mcmillan_ij;
extern DOUBLE *order;
extern int *sign_u_mi, *sign_u_ij;
extern DOUBLE *mu_k, *mu_p_k, *mu_pp_k;

extern DOUBLE E, EFF, Eo, Ekin, Epot, Elat, Evar, Eck, Edamping, Eint, Eext;
extern DOUBLE xR2, zR2;
extern DOUBLE Srecoil, Syrecoil; 
extern int Nlattice, Nlattice3;

extern clock_t time_start;
extern clock_t time_simulating;
extern clock_t time_measuring;
extern clock_t time_communicating;

extern int *dead_walkers_storage;
extern int *walkers_storage;
extern int *killed_walkers_storage;
extern int Ndead_walkers;
extern DOUBLE reduce, amplify;
extern DOUBLE *branching_weight; // [NwalkersMax]
extern DOUBLE *branching_xi;
extern int *branching_multiplicity;

extern int overlaped;
extern int tried, Nmeasured;
extern DOUBLE D;
extern int square_box;
extern int measure_Lind;
extern DOUBLE lattice_length;

extern int Nscal;
extern int Ncells;
extern DOUBLE Lcell;

extern int measure_Nk_pure;

extern DOUBLE Rzigzag;
extern DOUBLE T;
extern DOUBLE Tannealing;

extern int R2times_measured;
extern DOUBLE sum_over_images_cut_off_num;
extern DOUBLE Rc_smooth_cutoff;
extern DOUBLE NkMaxTrap;

extern int file_particles_append;
extern int measure_energy_barrier;
extern int file_append;
extern int optimization;
extern int Niter_store;
extern int measure_SkMATRIX;
extern int measure_FormFactor;
extern int measure_Lindemann_Lozovik;
extern int generate_crystal_coordinates;
extern int measure_wavefunction_projection;
extern int measure_pure_coordinates;
extern int measure_tVMC;
extern int R2_subtract_CM;
extern int OBDM_subtract_CM;
extern int measure_Hessian_matrix;

extern DOUBLE inwf_power;
extern DOUBLE Vo, Ro;
extern DOUBLE delta_tilde;

extern DOUBLE solitonV2, solitonXi, solitonk, soliton_fixed_phase_V2, soliton_fixed_phase_Xi;
extern DOUBLE R3; // hyperradius
extern DOUBLE bilayer_width, bilayer_width2;
extern DOUBLE t_tunneling;
extern DOUBLE Uo;
extern DOUBLE gaussian_alpha;
extern DOUBLE gaussian_beta;
extern DOUBLE RoSW;
extern DOUBLE Epotcutoff;
extern DOUBLE cR[10], cI[10], tvmcNpar, tvmcNobs;
extern DOUBLE Kpar11, Kpar12, aA;
extern DOUBLE CMseparation;
extern DOUBLE hard_core_diameter;
extern int Kurbakov_BC;
extern DOUBLE Kurbakov_softcore;
extern DOUBLE m_mu;
extern int check_kinetic_energy;
extern DOUBLE one_body_momentum_i;
extern DOUBLE one_body_momentum_k;
extern DOUBLE E_from_growth_estimator;
extern DOUBLE Rsphere;
extern DOUBLE Rnow;
extern DOUBLE sum_of_two_wavefunctions_parameter1, sum_of_two_wavefunctions_parameter2, sum_of_two_wavefunctions_alpha;
extern DOUBLE numeric_derivative_dr;

void Measure(int, int);

#endif