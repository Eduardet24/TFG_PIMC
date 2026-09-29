/*crystal.c*/

#include <stdio.h>
#include "main.h"
#include "crystal.h"
#include "utils.h"
#include "memory.h"
#include "vmc.h"
#include "compatab.h"
#include "randnorm.h"
#include "trial.h"
#include MATHINCLUDE

DOUBLE crystal_side_x= 1;
DOUBLE crystal_side_y= 1;

#ifdef CRYSTAL

struct Cryst Crystal;

/************************** Construct Grid Crystal ***************************/
void ConstructGridCrystal(struct Cryst *Crystal, int size, DOUBLE n) {
  int i;

#ifdef BC_3DPBC
#ifdef LATTICE_3D_BCC // BCC lattice
  DOUBLE basisx[2]={0.,0.5};
  DOUBLE basisy[2]={0.,0.5};
  DOUBLE basisz[2]={0.,0.5};
  int N_atoms_in_cell = 2;
#endif
#ifdef LATTICE_3D_FCC // FCC lattice
  DOUBLE basisx[4]={0.,0.5,0.0,0.5};
  DOUBLE basisy[4]={0.,0.5,0.5,0.0};
  DOUBLE basisz[4]={0.,0.0,0.5,0.5};
  int N_atoms_in_cell = 4;
#endif
#ifdef LATTICE_3D_CUBIC // FCC lattice
  DOUBLE basisx[4]={0};
  DOUBLE basisy[4]={0};
  DOUBLE basisz[4]={0};
  int N_atoms_in_cell = 1;
#endif
  DOUBLE lattice[3];

  DOUBLE side;
  int nc, j, k, ll;
  int index; 
#endif
#ifdef BC_2DPBC
  int index,j,k;
  DOUBLE basis[2][2]={{0,0},{0.5,0.5}};
  DOUBLE lattice[3];
  int nc;
  DOUBLE ratio;
  int Imax, Jmax;
  int invert;
#ifdef LATTICE_TRIANGULAR
  DOUBLE aspect_ratio, best_aspect_ratio;
  int i_best, j_best;
#endif
#ifdef SPINFULL
  int s;
#endif
#endif

#ifdef ONE_BODY_IMPURITY
  size = 1;
  NCrystal = 1;
#endif

  if(verbosity) Message("\nMemory allocation : allocating crystal grid of size %i...", size);

  Crystal->x   = (DOUBLE*) Calloc("Crystal->x   ", size, sizeof(DOUBLE));
  Crystal->y   = (DOUBLE*) Calloc("Crystal->y   ", size, sizeof(DOUBLE));
  Crystal->z   = (DOUBLE*) Calloc("Crystal->z   ", size, sizeof(DOUBLE));
  Crystal->size = size;

#ifdef ONE_BODY_IMPURITY
  Crystal->size = 1;
  Crystal->x[0] = 0;
  Crystal->y[0] = 0;
  Crystal->z[0] = 0.5*L;
  Warning("  Impurity is positioned at L/2");
  return;
#endif

#ifdef CRYSTAL_WIDTH_ARRAY
  Message("  Parameters Crystal->Rx[i], Ry[i], Rz[i] are independent and will be loaded from file\n");
  Crystal->Rx   = (DOUBLE*) Calloc("Crystal->Rx   ", size, sizeof(DOUBLE));
  Crystal->Ry   = (DOUBLE*) Calloc("Crystal->Ry   ", size, sizeof(DOUBLE));
  Crystal->Rz   = (DOUBLE*) Calloc("Crystal->Rz   ", size, sizeof(DOUBLE));
  Crystal->weight = (DOUBLE*) Calloc("Crystal->weight   ", size, sizeof(DOUBLE));

  for(i=0; i<size; i++) Crystal->weight[i] = 1.;
#else
  Message("  Parameters Crystal->Rx, Ry, Rz are same for all particles\n");
  if(alpha_R == -PI && (alpha_Rx == -PI && alpha_Ry == -PI && alpha_Rz == -PI)) Error("Crystal parameters (alpha_R) are not defined!\n");
  if(alpha_R == 0)  Warning("Suspicious variational parameter in Crystal alpha_R = 0 !\n");
  if(alpha_Rx == -PI) {
    Warning("  Setting alpha_Rx = alpha_R = %lf\n", alpha_R);
    alpha_Rx = alpha_R;
  }
  if(alpha_Ry == -PI) {
    Warning("  Setting alpha_Ry = alpha_R = %lf\n", alpha_R);
    alpha_Ry = alpha_R;
  }
  if(alpha_Rz == -PI) {
    Warning("  Setting alpha_Rz = alpha_R = %lf\n", alpha_R);
    alpha_Rz = alpha_R;
  }
#ifdef TRIAL_2D
  Warning("  Setting alpha_Rz = 0\n");
  alpha_Rz = 0.;
#endif

  Crystal->Rx = alpha_Rx;
  Crystal->Ry = alpha_Ry;
  Crystal->Rz = alpha_Rz;
#endif

#ifdef BC_1DPBC_Z
  for(i=0; i<size; i++) {
    Crystal->x[i] = 0;
    Crystal->y[i] = 0;
    Crystal->z[i] = L/size*(0.5+i);
  }
  Crystal->size = size;
#endif

#ifdef BC_2DPBC // 2D system
//� f(r_i) = Exp(-alpha (r_i-r^c_i)^2)
  Message("  One-body term of the trial wavefunction: Crystal->Rx = %" LE " \n", Crystal->Rx);
  Message("                                           Crystal->Ry = %" LE " \n", Crystal->Ry);
  Message("                                           Crystal->Rz = %" LE " \n", Crystal->Rz);

#ifndef LATTICE_TRIANGULAR // square crystal
#ifdef BC_2DPBC
  Message("  Square lattice was chosen\n");
  Message("  2D crystal is being constructed\n");
#endif
#ifdef BC_1DPBC_X
  Message("  Linear (in X direction) crystal is being constructed\n");
#ifdef LATTICE_ZIGZAG
  Message("  Zig-zag symmetry will be used. Displacement Rzigzag= %" LF "\n", Rzigzag);
  if(N%2==1) Error("  Number of particles must be even in order to use Zig zag symmetry!\n");
#endif
#endif
  square_box = ON;
#ifdef BC_2DPBC
  index = 0;

  if((int)(Sqrt(size))*(int)(Sqrt(size)) != size) Error("number of crystal sites is incompatible with square lattice");

  for(i=0; i<Sqrt(size); i++) {
    for(j=0; j<Sqrt(size); j++) {
      if(index>size) Error("crystal index out of range\n");
      Crystal->x[index] = L/Sqrt(size)*(0.5+i);
      Crystal->y[index] = L/Sqrt(size)*(0.5+j);
      Crystal->z[index] = 0;
      index++;
    }
  }

#endif
#ifdef BC_1DPBC_X
#ifdef LATTICE_ZIGZAG
  for(i=0; i<size; i++) {
    Crystal->x[i] = L/size*(0.5+i);
    Crystal->y[i] = (i%2==0)?(Rzigzag):(-Rzigzag);
    Crystal->z[i] = 0;
  }
  index = size;
#else // linear crystal
  for(i=0; i<size; i++) {
    Crystal->x[i] = L/size*(0.5+i);
    Crystal->y[i] = 0;
    Crystal->z[i] = 0;
  }
  index = size;
#endif

#endif
#ifdef LATTICE_CHECKERBOARD
  Warning("  Checkerboard crystal w.f. will be used!");

  index = 0;
  for(i=0; i<Sqrt(size); i++) {
    for(j=0; j<Sqrt(size); j++) {
      if((i+j)%2) {
        Crystal->x[index] = L/Sqrt(size)*(0.5+i);
        Crystal->y[index] = L/Sqrt(size)*(0.5+j);
        Crystal->z[index] = 0;
        index++;
      }
    }
  }
#endif
#ifdef LATTICE_DOMINO
  Warning("  \"Domino\" crystal w.f. will be used!");

  index = 0;
  for(i=0; i<Sqrt(size); i++) {
    for(j=0; j<Sqrt(size)/2; j++) {
      if((i+j)%2) {
        Crystal->x[index] = L/Sqrt(size)*(0.5+i);
        Crystal->y[index] = L/Sqrt(size)*(0.5+2*j);
        Crystal->z[index] = 0;
        index++;
        Crystal->x[index] = L/Sqrt(size)*(0.5+i);
        Crystal->y[index] = L/Sqrt(size)*(0.5+2*j+1);
        Crystal->z[index] = 0;
        index++;
      }
    }
  }
#endif
  size = index;
  Crystal->size = size;
  Message("  Particles will be localized close to %i sites\n", Crystal->size);

#else
  // "pseudo-elementary cell" contains 2 atoms: 2I in x direction, J in y direction. lx=3^1/2 ly
  // depending on the number of atoms x and y can be effectively exchanged, so that always Lx<Ly
  Message("  Triangular lattice has been chosen\n");
  square_box = OFF;
  nc = 2; // atoms in a "pseudo-elementary" rectangular cell
#ifdef BC_2DPBC_NON_ORTHOGONAL
  Warning("  Non-orthogonal basis will be used for finding nearest images.\n");
#endif

#ifdef SPINFULL
  Message("  Changing lattice size from %i to %i \n" , size, size/Nspin );
  size /= Nspin;
#endif

  // fix number of rows and columns
  //I = (int)(Sqrt(0.5*size/Sqrt(3.))+0.5);
  //J = (int)(Sqrt(3.)*(DOUBLE)I+0.5);
  //if(2*Imax*Jmax != size) Error("  Cannot construct triangular lattice! Adjust number of lattice sites NCrystal = %i -> %i\n", size, 2*Imax*Jmax);

  Message("  Constructing best triangular lattice...");
//  for(N=10; N<=1000; N+=2) {
  best_aspect_ratio = 1e10;
  for(i=1; i<NCrystal/2; i++) {
    for(j=1; j<NCrystal; j++) {
      //if(2*i*j == N) {
      if(2*i*j == size) {
        // suppose we have j lattice sites in X direction and (2i) lattice sites in Y direction
        aspect_ratio = (DOUBLE) i * Sqrt(3) / ((DOUBLE) j);
        if(fabs(aspect_ratio-1.)<best_aspect_ratio) {
          best_aspect_ratio = fabs(aspect_ratio-1.);
          i_best = i;
          j_best = j;
        }
      }
    }
  }
  Imax = i_best;
  Jmax = j_best;
  ratio = Sqrt(3)*(DOUBLE)Imax/(DOUBLE)Jmax;

  if(best_aspect_ratio == 1e10)
    Message("  Cannot construct triangular lattice!");
  else
    Message(" done\n  Lattice size: %i = %i x %i, ratio = %" LG "\n", N, 2*Imax, Jmax, ratio<1?ratio:1./ratio);
    //Message(" %3i %3i %3i %" LG "\n", N, J, 2*Imax, ratio<1?ratio:1./ratio);

  if(ratio>1) { // should be Lx / Ly < 1
    Message("  Size of the rectangle in x direction is larger than in y direction\n");
    invert = ON;
  }
  else {
    Message("  Size in x direction is smaller than in y direction, rotating the lattice\n");
    invert = OFF;
  }

  crystal_side_x = Sqrt(size*Nspin*ratio/n) / (DOUBLE)Imax; // i.e. assuming Lx<Ly
  crystal_side_y = crystal_side_x / Sqrt(3.);
  Ly = Sqrt(size*Nspin/(n*ratio));
  Lx = Sqrt(size*Nspin/n*ratio);

  index = 0;
  for(i=0; i<Imax; i++) {
    for(j=0; j<Jmax; j++) {
      for(k=0; k<nc; k++) {
        lattice[0] = crystal_side_x*((DOUBLE)i+0.25) + basis[k][0]*crystal_side_x;
        lattice[1] = crystal_side_y*((DOUBLE)j+0.25) + basis[k][1]*crystal_side_y;
        lattice[2] = 0;

        if(lattice[0]>Lx) Error(" lattice site outside Lx!\n");
        if(lattice[1]>Ly) Error(" lattice site outside Ly!\n");
        if(index>= size*Nspin) Error(" Failed to construct the grid!\n");
#ifdef SPINFULL
        for (s = 0; s < Nspin; s++)
#endif 

        {
          if (invert == OFF) {
            Crystal->x[index] = lattice[0];
            Crystal->y[index] = lattice[1];
            Crystal->z[index] = lattice[2];
          }
          else {
            Crystal->x[index] = lattice[1];
            Crystal->y[index] = lattice[0];
            Crystal->z[index] = lattice[2];
          }
          index++;
        }
      }
    }
  }
  if(invert == OFF) {
    Lx = Imax*crystal_side_x;
    Ly = Jmax*crystal_side_y;
#ifdef BC_2DPBC_NON_ORTHOGONAL
    Crystal->e1x = 0;
    Crystal->e1y = crystal_side_y*J;
    Crystal->e2x = crystal_side_y/2.*2*Imax;
    Crystal->e2y = crystal_side_y*Sqrt(3.)/2.*2*Imax;
#endif
  }
  else {
    Ly = Imax*crystal_side_x;
    Lx = Jmax*crystal_side_y;
#ifdef BC_2DPBC_NON_ORTHOGONAL
    Crystal->e1x = crystal_side_y*J;
    Crystal->e1y = 0;
    Crystal->e2x = crystal_side_y*Sqrt(3.)/2.*2*Imax;
    Crystal->e2y = crystal_side_y/2.*2*Imax;
#endif
  }
  L = Lx;
  Lhalf = 0.5*L;
  Lhalf2 = Lhalf*Lhalf;
  L2 = L*L;
  Message("  Size of simulation box Lx = %" LF ", Ly = %" LF ". Ratio: Lx/Ly = %" LF "\n", Lx, Ly, Lx/Ly);
  Message("  Density of particles:     N        / (Lx Ly)  = %" LF "\n", (DOUBLE) N / (Lx*Ly));
  Message("  Density of lattice sites: NCrystal / (Lx Ly)  = %" LF "\n", (DOUBLE) NCrystal / (Lx*Ly));
  Message("  Elementary cell side: %" LF " %" LF "\n",crystal_side_x, crystal_side_y);
  if(invert == OFF) {
    Message("Lattice size %i x %i\n", 2*Imax, Jmax);
    Message("Peaks in Sk should appear at (%i, %i) and (%i, %i)\n",2*Imax, 0, Imax, Jmax);
  }
  else {
    Message("Lattice size %i x %i\n", Jmax, 2*Imax);
    Message("Peaks in Sk should appear at (%i, %i) and (%i, %i)\n",0,2*Imax, Jmax, Imax);
  }

  lattice_length = crystal_side_y;
  //if(Lx != I * crystal_side_x) Error("  Inconsistent size of Lx %" LF " %" LF "\n", Lx, I*crystal_side_x);
  //if(Ly != J * crystal_side_y) Error("  Inconsistent size of Ly %" LF " %" LF "\n", Ly, J*crystal_side_y);
  if(fabs(size*Nspin/(Lx*Ly)-n)/n>1e-8) Error("  Wrong density was obtained (%" LF " instead of %" LF ")\n", size*Nspin/(Lx*Ly), n);

#ifdef BC_2DPBC_NON_ORTHOGONAL
  Crystal->e1abs_inv2 = 1./(Crystal->e1x*Crystal->e1x+Crystal->e1y*Crystal->e1y);
  Crystal->e2abs_inv2 = 1./(Crystal->e2x*Crystal->e2x+Crystal->e2y*Crystal->e2y);
  Message("  Volume " "%" LF " %" LF " %" LF  "\n", Lx*Ly, Crystal->e1x*Crystal->e2x+Crystal->e1y*Crystal->e2y);
#endif

  if(index != size*Nspin) Error("  Failed to construct the grid! Expecting %i lattice sites, found %i\n", size, index);
#endif
#endif

#ifdef BC_3DPBC
  Message("  One-body term of the trial wave function: Crystal->Rx = %" LE " \n", Crystal->Rx);
  Message("  One-body term of the trial wave function: Crystal->Ry = %" LE " \n", Crystal->Ry);
  Message("  One-body term of the trial wave function: Crystal->Rz = %" LE " \n", Crystal->Rz);

#ifdef LATTICE_3D_BCC // BCC lattice
  side=pow(2./n, 1./3.);
  lattice_length = side;
  Message("  BCC lattice has been chosen\n  Elementary cell side: %10.5f\n",side);

  nc = (int) (pow(size/2, 1./3.)+1e-3);

  if(2*nc*nc*nc != size) Error("  Inconsistent number of lattice sites for BCC lattice: NCrystal=2 %i^3\n", nc);
#endif
#ifdef LATTICE_3D_FCC // FCC lattice
  side=pow(4./n, 1./3.);
  lattice_length = side;
  Message("  FCC lattice has been chosen\n  Elementary cell side: %10.5f\n",side);

  nc = (int) pow(size*(0.25+1e-5), 1./3.);

  if(4*nc*nc*nc != size) Error("  Inconsistent number of lattice sites for FCC lattice: NCrystal=4 %i^3\n", nc);
#endif
#ifdef LATTICE_3D_CUBIC
  side=pow(1./n, 1./3.);
  lattice_length = side;
  Message("  FCC lattice has been chosen\n  Elementary cell side: %10.5f\n",side);

  nc = (int) pow(size*(1.+1e-5), 1./3.);

  if(nc*nc*nc != size) Error("  Inconsistent number of lattice sites for CUBIC lattice: NCrystal= %i^3\n", nc);
#endif
  Lx = Ly = Lz = pow((DOUBLE)N/n, 1./3.);

  index = 0;
  for(i=0; i<nc; i++) {
    for(j=0; j<nc; j++) {
      for(k=0; k<nc; k++) {
         for(ll=0; ll<N_atoms_in_cell; ll++) {
           lattice[0] = ((DOUBLE)i + basisx[ll] + 0.5)*side;
           lattice[1] = ((DOUBLE)j + basisy[ll] + 0.5)*side;
           lattice[2] = ((DOUBLE)k + basisz[ll] + 0.5)*side;

           ReduceToTheBox(lattice);

           Crystal->x[index] = lattice[0];
           Crystal->y[index] = lattice[1];
           Crystal->z[index] = lattice[2];

           index++;
         }
      }
    }
  }

  if(index != size) Warning("  Problems in construction of the crystal\n");
#endif

#ifndef CRYSTAL_WIDTH_ARRAY
  Crystal->Rx_av = Crystal->Rx;
  Crystal->Ry_av = Crystal->Ry;
  Crystal->Rz_av = Crystal->Rz;
#endif

  if(verbosity) Message(" done\n");
}

/************************** Construct Grid Impurity **************************/
void ConstructGridImpurity(struct Cryst *Crystal, const DOUBLE R, int size) {
#ifndef CRYSTAL_WIDTH_ARRAY
  DOUBLE x;
  DOUBLE xmin, xmax;
  DOUBLE y, ymin, ymax;
  DOUBLE precision = 1e-8;
  int search = ON;
  int i;
  DOUBLE alpha;
#ifdef ONE_BODY_IMPURITY
  FILE *out;
  DOUBLE dx;
#endif

#ifdef ONE_BODY_IMPURITY_1D
  Message("1D impurities\n");
#endif
#ifdef ONE_BODY_IMPURITY_3DSW
  static int first_time = ON;
  DOUBLE r;
  DOUBLE mu;
  DOUBLE a;
  Message("3D impurities with Square Well interactions\n");

  a = b;

  if(a < 0) {
    Warning("  Automatically changing scattering length to positive value\n");
    a = -a;
    b = -b;
  }
#endif

  if(verbosity) Message("\nMemory allocation : allocating impurity grid of size %i...", size);

  Crystal->x = (DOUBLE*) Calloc("Crystal->x   ", size, sizeof(DOUBLE));
  Crystal->y = (DOUBLE*) Calloc("Crystal->y   ", size, sizeof(DOUBLE));
  Crystal->z = (DOUBLE*) Calloc("Crystal->z   ", size, sizeof(DOUBLE));
  Crystal->size = size;

#ifdef ONE_BODY_IMPURITY_1D
  if(R == 0) {
    Crystal->Rz = Lhalfwf;
    Warning(" Variational parameter Rpar => %" LG "\n", Lhalfwf);
  }
  else
    Crystal->Rz = Lwf*R;
  Message("  One-body term of the trial wavefunction: Crystal->R = %" LE " \n", Crystal->Rz);
#endif

  for(i=0; i<size; i++) {
    Crystal->x[i] = 0.;
    Crystal->y[i] = 0.;
    Crystal->z[i] = 0.;
#ifdef ONE_BODY_IMPURITY_1D
    Crystal->z[i] = Lwf/size*(0.5+i);
#endif
#ifdef ONE_BODY_IMPURITY_3DSW
    if(i==0)
      Crystal->x[i] = Crystal->y[i] = Crystal->z[i] = 0.5*L;
    else {
      Crystal->x[i] = L*Random();
      Crystal->x[i] = L*Random();
      Crystal->x[i] = L*Random();
    }
#endif
  }

  Message("  Constructing wavefuncion...\n");

#ifdef ONE_BODY_IMPURITY_1D
  /* x tg xR = 1/|b| */
  alpha = 1./b;
  xmin = (1e-10)/(2.*Crystal->Rz);
  xmax = (PI-1e-10)/(2.*Crystal->Rz);
  ymin = xmin * Tg(xmin*Crystal->Rz) - alpha;
  ymax = xmax * Tg(xmax*Crystal->Rz) - alpha;

  if(ymin*ymax > 0) Exit(1, "Can't construct the grid : No solution found");

  while(search) {
    x = (xmin+xmax)/2.;
    y = x * Tg(x*Crystal->Rz) - alpha;

    if(y*ymin < 0) {
      xmax = x;
      ymax = y;
    } else {
      xmin = x;
      ymin = y;
    }

    if(fabs(y/alpha)<precision) {
      x = (xmax*ymin-xmin*ymax) / (ymin-ymax);
      y = x * Tg(x*Crystal->Rz) - alpha;
      search = OFF;
    }
  }

  Crystal->k = x;
  Message("  k = %" LG "\n", Crystal->k);
  Crystal->k2 = x*x;
#endif

#ifdef ONE_BODY_IMPURITY_3DSW
  // reduced mass
  mu = 1.;

  // kappa = sqrt(2 mu Vo ) / hbar
  trial_Kappa = sqrt(2.*mu*trial_Vo);
  trial_Kappa2 = 2.*mu*trial_Vo;

  if(first_time) {
    Message("  BCS square well zero energy + const, var. par.: Rpar12 and Bpar12\n");
    Warning("  Scattering length must be negative in order to use this w.f.\n");
    Warning("  Sin(K r)/r [r<RoSW], 1+a/r [r<Rpar12], 1. + A*(Exp(-Bpar12*r)+Exp(Bpar12*(r-L))) [r<L/2]\n");
    Warning("  Changing Rpar12 %"LE" -> %"LE"\n", Rpar12, Rpar12*Lhalf);
    first_time = OFF;
    Rpar12 *= Lhalf;
  }

  if(RoSW > Rpar12) Error(" RoSW = %lf > Rpar12 = %lf = %lf L/2, must be otherwise, i.e. Rpar12/[L/2] > %lf\n", RoSW, Rpar12, Rpar12/Lhalf, RoSW/Lhalf);

  A3SW = 1./(1.+a/Rpar12*(1.-1./(Rpar12*Bpar12)*(Exp(Bpar12*(L-2*Rpar12))+1.)/(Exp(Bpar12*(L-2*Rpar12))-1.)));
  r = RoSW;
  A1SW = A3SW*(1.+a/r) / (Sin(trial_Kappa*r)/r);
  A4SW = A3SW*a/(Rpar12*Rpar12*Bpar12)/(Exp(-Bpar12*Rpar12)-Exp(Bpar12*(Rpar12-L)));

  //Message("  Check a = %"LE", goal value %"LE"\n", RoSW*(1.-tg(kappa*R)/(kappa*R)), a);

  r = RoSW;
  r = (trial_Kappa / tan(trial_Kappa*r) - 1. / r) / (-a / (r*r) / (1. + a / r)); // testing continuity of the first derivative
  if(fabs(r-1)>0.1) Warning("  ConstructGridImpurity: first derivative is discontinuous");
#endif

  if(verbosity) Message(" done\n");

#ifdef ONE_BODY_IMPURITY
  out = fopen("wfimp.dat", "w"); // Save wavefunction
  Message("Saving wavefunction... ");

  //xmax = Crystal->Rz*1.5;
  xmax = L/2.;
  xmin = xmax/N;
  dx = (xmax-xmin) / N;

  fprintf(out, "x f lnf fp Eloc E\n");

  for(i=1; i<1000; i++) {
    x = xmin + dx * (DOUBLE) (i);

    fprintf(out, "%.15" LE " %.15" LE " %.15" LE " %.15" LE " %.15" LE " %.15" LE "\n", x, Exp(ImpurityInterpolateU(x)), ImpurityInterpolateU(x),
        ImpurityInterpolateFp(x), ImpurityInterpolateE(x)-ImpurityInterpolateFp(x)*ImpurityInterpolateFp(x), ImpurityInterpolateE(x));
  }
  fclose(out);
  Message("done.\n");
#endif
#endif
}

#ifdef ONE_BODY_IMPURITY
/************************ Interpolate U ************************************/
// U = ln(f)
#ifdef ONE_BODY_IMPURITY_1D
DOUBLE ImpurityInterpolateU(DOUBLE x) {
  if(x<Crystal.Rz)
    return Log(Cos(Crystal.k*(x-Crystal.Rz)));
  else
    return 0;
}

/************************ InterpolateFp ************************************/
// fp = f'/f
DOUBLE ImpurityInterpolateFp(DOUBLE x) {
  if(x<Crystal.Rz)
    return -Crystal.k * Tg(Crystal.k*(x-Crystal.Rz));
  else
    return 0;
}

/************************ InterpolateE ************************************/
// E =  -(f"/f - (f'/f)^2);
DOUBLE ImpurityInterpolateE(DOUBLE x) {
  DOUBLE c;

  if(x<Crystal.Rz) {
    c = Tg(Crystal.k*(x-Crystal.Rz));
    return Crystal.k2*(1. + c*c);
  }
  else
    return 0;
}
#endif
#endif

#ifdef ONE_BODY_IMPURITY_3DSW
DOUBLE ImpurityInterpolateU(DOUBLE r) {
  if(r<RoSW)
    return Log(A1SW*Sin(trial_Kappa*r)/r);
  else if(r<Rpar12)
    //return Log(A3SW*(1.+a/r));
    return Log(A3SW*(1.+b/r));
  else if(r<Lhalf)
    return Log(1. + A4SW*(Exp(-Bpar12*r)+Exp(Bpar12*(r-L))));
  else
    return 0;
}

/************************ InterpolateFp ************************************/
// fp = f'/f
DOUBLE ImpurityInterpolateFp(DOUBLE r) {
  if(r<RoSW)
    return trial_Kappa/tan(trial_Kappa*r) - 1./r;
  else if(r<Rpar12)
    //return -a/(r*r)/(1.+a/r);
    return -b/(r*r)/(1.+b/r);
  else if(r<Lhalf)
    return A4SW*Bpar12*(-Exp(-Bpar12*r)+Exp(Bpar12*(r-L)))/(1. + A4SW*(Exp(-Bpar12*r)+Exp(Bpar12*(r-L))));
  else
    return 0.;
}

/************************ InterpolateE ************************************/
// E =  -(f"/f - (f'/f)^2);
DOUBLE ImpurityInterpolateE(DOUBLE r) {
  DOUBLE fp;

  fp = ImpurityInterpolateFp(r);
  if(r<RoSW)
    return trial_Kappa2 + fp*fp;
  else if(r<Rpar12)
    return fp*fp;
  else if(r<Lhalf)
    return -((Bpar12*A4SW*(Exp(Bpar12*L)*(-2 + Bpar12*r) + Exp(2*Bpar12*r)*(2 + Bpar12*r)))/((Exp(Bpar12*(L + r)) + A4SW*(Exp(Bpar12*L) + Exp(2*Bpar12*r)))*r)) + fp*fp;
  else
    return 0;
}
#endif


/************************** Chains pinned to center of mass, Calculate center of mass ***************************/
int ChainsPinnedCenterMass(void) {
#ifdef CHAINS_PINNED_TO_CENTER_OF_MASS
  Message("  \nChains pinned to center of mass\n");
  int w, i, j;
  for (w = 0; w < Nwalkers; w++) {
    for (i = 0; i < N / Nspin; i++) { // i = Number of chain
      CaseX(W[w].Crystal_x[i * Nspin] = 0.);
      CaseY(W[w].Crystal_y[i * Nspin] = 0.);
      CaseZ(W[w].Crystal_z[i * Nspin] = 0.);
      for (j = i * Nspin; j < Nspin * (i + 1); j++) { // j  Chain index
        CaseX(W[w].Crystal_x[i * Nspin] += W[w].x[j]);
        CaseY(W[w].Crystal_y[i * Nspin] += W[w].y[j]);
        CaseZ(W[w].Crystal_z[i * Nspin] += W[w].z[j]);
      }
      CaseX(W[w].Crystal_x[i * Nspin] /= (DOUBLE)Nspin);
      CaseY(W[w].Crystal_y[i * Nspin] /= (DOUBLE)Nspin);
      CaseZ(W[w].Crystal_z[i * Nspin] /= (DOUBLE)Nspin);
      for (j = i * Nspin + 1; j < Nspin * (i + 1); j++) {
        CaseX(W[w].Crystal_x[j] = W[w].Crystal_x[i * Nspin]);
        CaseY(W[w].Crystal_y[j] = W[w].Crystal_y[i * Nspin]);
        CaseZ(W[w].Crystal_z[j] = W[w].Crystal_z[i * Nspin]);
      }
    }
  }
#endif
  return 0;
}

/************************** Calculate All Chains Center Of Mass ***************************/
#ifdef CHAINS_PINNED_TO_CENTER_OF_MASS
void CalculateAllChainsCenterOfMass(DOUBLE **R, int w) {
  int i, j;
  DOUBLE dr[3] = { 0.,0.,0 }, r2;
    for (i = 0; i < N / Nspin; i++) { // i = Number of chain
      CaseX(W[w].Crystal_x[i*Nspin] = 0.);
      CaseY(W[w].Crystal_y[i*Nspin] = 0.);
      CaseZ(W[w].Crystal_z[i*Nspin] = 0.);
      for (j = i * Nspin; j < Nspin * (i + 1); j++) { // j  Chain index

        CaseX(dr[0] = R[i*Nspin][0] - R[j][0]);
        CaseY(dr[1] = R[i*Nspin][1] - R[j][1]);
        CaseZ(dr[2] = R[i*Nspin][2] - R[j][2]);

        r2 = FindNearestImageCenterOfMass(&dr[0], &dr[1], &dr[2]);

        CaseX(W[w].Crystal_x[i*Nspin] += R[j][0] + dr[0]);
        CaseY(W[w].Crystal_y[i*Nspin] += R[j][1] + dr[1]);
        CaseZ(W[w].Crystal_z[i*Nspin] += R[j][2] + dr[2]);
      }
      CaseX(W[w].Crystal_x[i*Nspin] /= (DOUBLE)Nspin);
      CaseY(W[w].Crystal_y[i*Nspin] /= (DOUBLE)Nspin);
      CaseZ(W[w].Crystal_z[i*Nspin] /= (DOUBLE)Nspin);
      for (j = (i*Nspin) + 1; j < Nspin*(i + 1); j++) {
        CaseX(W[w].Crystal_x[j] = W[w].Crystal_x[i*Nspin]);
        CaseY(W[w].Crystal_y[j] = W[w].Crystal_y[i*Nspin]);
        CaseZ(W[w].Crystal_z[j] = W[w].Crystal_z[i*Nspin]);
      }
    }
}
#endif

/************************** Calculate All Chains Center Of Mass of One Walker***************************/
#ifdef CHAINS_PINNED_TO_CENTER_OF_MASS
void CalculateAllChainsCenterOfMassWalker(struct Walker *walker) {
  int i, j;
  DOUBLE dr[3] = { 0.,0.,0 },r2;
  for (i = 0; i < N / Nspin; i++) { // i = Number of chain = {0,1,2,...,N/Nspin - 1}
    CaseX(W[walker->w].Crystal_x[i*Nspin] = 0.);
    CaseY(W[walker->w].Crystal_y[i*Nspin] = 0.);
    CaseZ(W[walker->w].Crystal_z[i*Nspin] = 0.);

    for (j = i*Nspin; j < Nspin*(i + 1); j++) { // j  particle index
      CaseX(dr[0] = walker->x[i*Nspin] - walker->x[j]);
      CaseY(dr[1] = walker->y[i*Nspin] - walker->y[j]);
      CaseZ(dr[2] = walker->z[i*Nspin] - walker->z[j]);

      r2 = FindNearestImageCenterOfMass(&dr[0], &dr[1], &dr[2]);

      CaseX(W[walker->w].Crystal_x[i*Nspin] += walker->x[j] + dr[0]);
      CaseY(W[walker->w].Crystal_y[i*Nspin] += walker->y[j] + dr[1]);
      CaseZ(W[walker->w].Crystal_z[i*Nspin] += walker->z[j] + dr[2]);
    }
    CaseX(W[walker->w].Crystal_x[i*Nspin] /= (DOUBLE)Nspin);
    CaseY(W[walker->w].Crystal_y[i*Nspin] /= (DOUBLE)Nspin);
    CaseZ(W[walker->w].Crystal_z[i*Nspin] /= (DOUBLE)Nspin);
    for (j = (i*Nspin) + 1; j < Nspin*(i + 1); j++) {
      CaseX(W[walker->w].Crystal_x[j] = W[walker->w].Crystal_x[i*Nspin]);
      CaseY(W[walker->w].Crystal_y[j] = W[walker->w].Crystal_y[i*Nspin]);
      CaseZ(W[walker->w].Crystal_z[j] = W[walker->w].Crystal_z[i*Nspin]);
    }
  }
}
#endif

/************************** Calculate One Chain Center Of Mass ***************************/
#ifdef CHAINS_PINNED_TO_CENTER_OF_MASS
void CalculateOneChainCenterOfMass(DOUBLE x, DOUBLE y, DOUBLE z, int i, int w) {
  int Nchain, j;
  DOUBLE dr[3] = { 0.,0.,0 }, r2;
  Nchain = i / Nspin;  // Number of chain
  CaseX(W[w].Crystal_x[Nchain*Nspin] = 0.);
  CaseY(W[w].Crystal_y[Nchain*Nspin] = 0.);
  CaseZ(W[w].Crystal_z[Nchain*Nspin] = 0.);
  
  if(i != Nchain*Nspin){
 
  for (j = Nchain*Nspin; j < Nspin*(Nchain + 1); j++) { // j  Chain index
    if (j != i) {
    CaseX(dr[0] = W[w].x[Nchain*Nspin] - W[w].x[j]);
    CaseY(dr[1] = W[w].y[Nchain*Nspin] - W[w].y[j]);
    CaseZ(dr[2] = W[w].z[Nchain*Nspin] - W[w].z[j]);
    r2 = FindNearestImageCenterOfMass(&dr[0], &dr[1], &dr[2]);
      CaseX(W[w].Crystal_x[Nchain*Nspin] += W[w].x[j] + dr[0]);
      CaseY(W[w].Crystal_y[Nchain*Nspin] += W[w].y[j] + dr[1]);
      CaseZ(W[w].Crystal_z[Nchain*Nspin] += W[w].z[j] + dr[2]);
    }
    else {
    CaseX(dr[0] = W[w].x[Nchain*Nspin] - x);
    CaseY(dr[1] = W[w].y[Nchain*Nspin] - y);
    CaseZ(dr[2] = W[w].z[Nchain*Nspin] - z);
    r2 = FindNearestImageCenterOfMass(&dr[0], &dr[1], &dr[2]);
      CaseX(W[w].Crystal_x[Nchain*Nspin] += x + dr[0]);
      CaseY(W[w].Crystal_y[Nchain*Nspin] += y + dr[1]);
      CaseZ(W[w].Crystal_z[Nchain*Nspin] += z + dr[2]);
    }
  }
  }
  else{

for (j = Nchain*Nspin; j < Nspin*(Nchain + 1); j++) { // j  Chain index
    if (j != i) {
    CaseX(dr[0] = x - W[w].x[j]);
    CaseY(dr[1] = y - W[w].y[j]);
    CaseZ(dr[2] = z - W[w].z[j]);
    r2 = FindNearestImageCenterOfMass(&dr[0], &dr[1], &dr[2]);
      CaseX(W[w].Crystal_x[Nchain*Nspin] += W[w].x[j] + dr[0]);
      CaseY(W[w].Crystal_y[Nchain*Nspin] += W[w].y[j] + dr[1]);
      CaseZ(W[w].Crystal_z[Nchain*Nspin] += W[w].z[j] + dr[2]);
    }
    else {
      CaseX(W[w].Crystal_x[Nchain*Nspin] += x);
      CaseY(W[w].Crystal_y[Nchain*Nspin] += y);
      CaseZ(W[w].Crystal_z[Nchain*Nspin] += z);
    }
  }

  }
  CaseX(W[w].Crystal_x[Nchain*Nspin] /= (DOUBLE)Nspin);
  CaseY(W[w].Crystal_y[Nchain*Nspin] /= (DOUBLE)Nspin);
  CaseZ(W[w].Crystal_z[Nchain*Nspin] /= (DOUBLE)Nspin);
  for (j = Nchain*Nspin + 1; j < Nspin*(Nchain + 1); j++) {
    CaseX(W[w].Crystal_x[j] = W[w].Crystal_x[Nchain*Nspin]);
    CaseY(W[w].Crystal_y[j] = W[w].Crystal_y[Nchain*Nspin]);
    CaseZ(W[w].Crystal_z[j] = W[w].Crystal_z[Nchain*Nspin]);
  }
}
#endif

void ConstructGridCrystalOneBody(void) {

  Message("  Constructing a lattice used in one-body terms, size NCrystal = %i\n", NCrystal);

  alpha_R = alpha_latt;

  ConstructGridCrystal(&Crystal, NCrystal, n * (DOUBLE) NCrystal / (DOUBLE) Ndens);
}

#endif
