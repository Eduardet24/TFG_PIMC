/*gencoord.c*/
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <time.h>
#include "main.h"
#include "memory.h"
#include "gencoord.h"
#include "vmc.h"
#include "randnorm.h"
#include "utils.h"
#include "compatab.h"
#include MATHINCLUDE
#include "rw.h"

/******************************* Generate Coordinates ************************/
int GenerateCoordinates(void) {
  int w,i;
  int go;
  int iter;

#ifdef SPINFULL_GENERATE_COORDENATES_IN_CHAIN
  int spin;
#endif 

  Message("  Generating particle coordinates ... ");

  if(generate_new_coordinates_from_lattice) Message("\n  Cubic lattice ordering will be used\n");

  Randomize();

  for(w=0; w<Nwalkers; w++) {
    go = ON;
    iter = 0;
    while(go) {
      for(i=0; i<N; i++) GenerateOneParticleCoordinates(&W[w].x[i],&W[w].y[i],&W[w].z[i], i);

#ifdef SPINFULL_GENERATE_COORDENATES_IN_CHAIN
      for(spin=1; spin<Nspin; spin++) {
        for (i = 0; i < N / Nspin; i++) {
          CaseX(W[w].x[Nspin*spin + i] = W[w].x[i] + Random()*L / (DOUBLE)N *0.01);
          CaseY(W[w].y[Nspin*spin + i] = W[w].y[i] + Random()*L / (DOUBLE)N *0.01);
          CaseZ(W[w].z[Nspin*spin + i] = W[w].z[i] + Random()*L / (DOUBLE)N *0.01);
        }
      }
#endif

#ifdef FIRST_PARTICLE_IS_NOT_MOVED
      CaseX(W[w].x[0] = 0.);
      CaseY(W[w].y[0] = 0.);
      CaseZ(W[w].z[0] = 0.);
#endif

      if(MC == PIGS) CopyWalkerCoord(&W[w], &W[0]); // for PIGS start from a classical configuration
      if(MC == PIMC) CopyWalkerCoord(&W[w], &W[0]); // for PIMC start from a classical configuration

      if(!OverlappingWalker(&W[w])) go = OFF;

      if(iter++ > MAX_GEN_ITER) Error("Can't generate initial configuration\n");
    }
    //Message("  generated %i-th walker %i iterations\n", w+1,iter);
#ifdef CENTER_OF_MASS_IS_NOT_MOVED // put the center of mass to zero
  AdjustCenterOfMassWalker(&W[w]);
#endif
#ifdef CENTER_OF_MASS_Z_IS_NOT_MOVED // put the center of mass to zero
  AdjustCenterOfMassWalker(&W[w]);
#endif
#ifdef CENTER_OF_MASS_IS_NOT_MOVED_TWO_PARTS // put the center of mass to zero for the first half and the second half of the system
  AdjustCenterOfMassWalker(&W[w]);
#endif
  }
  Nwalkersw = (DOUBLE) Nwalkers;

  if(CheckOverlapping()) Error(" generated configuration is wrong");

  SaveCoordinates("in3Dprev.in");

  Message("  done\n");
  return 0;
}

/******************************* Generate Coordinates ************************/
void GenerateOneParticleCoordinates(DOUBLE *x, DOUBLE *y, DOUBLE *z, int i) {
#ifdef VEXT_COS2
  int nx,ny,nz;
#endif
  int ix, iy, iz, index = 0, Nmax; // crystal indices
  DOUBLE l; // crystal spacing
  int search = ON;

  *x = *y = *z = 0.;

#ifdef TRIAL_3D //  i.e. 3D simulation
  if(boundary == NO_BOUNDARY_CONDITIONS) {
    if(alpha_x == 0 || alpha_y == 0 || alpha_z == 0) {
      Warning("  set alpha different from zero, initializing with 1\n");
      alpha_x = alpha_y = alpha_z = 1;
    }
    *x = (1-2.*RandomSys())*0.5/alpha_x;
    *y = (1-2.*RandomSys())*0.5/alpha_y;
    *z = (1-2.*RandomSys())*0.5/alpha_z;
  }
  else if(boundary == ONE_BOUNDARY_CONDITION) {
    if(alpha_x == 0 || alpha_y == 0) {
      Warning("  set alpha_x and alpha_y different from zero, initializing with 1\n");
      alpha_x = alpha_y = 1;
    }
    *x = (1-2.*RandomSys())*0.5/alpha_x;
    *y = (1-2.*RandomSys())*0.5/alpha_y;
    *z = L * RandomSys();
  }
  else if(boundary == TWO_BOUNDARY_CONDITIONS) {
    *x = L * RandomSys();
    *y = L * RandomSys();
    if(alpha_z == 0) {
      Warning("  set alpha_z different from zero, initializing with 1\n");
      alpha_z = 1;
    }
    *z = (1-2. * RandomSys())*0.5/alpha_z;
  }
  else if(boundary == THREE_BOUNDARY_CONDITIONS) {
    *x = L * RandomSys();
    *y = L * RandomSys();
    *z = L * RandomSys();
  }
#endif

#ifdef TRIAL_2D //  i.e. 2D simulation
  if(boundary == NO_BOUNDARY_CONDITIONS) {
    if(alpha_x == 0 || alpha_y == 0) {
      Warning("  set alpha_x and alpha_y different from zero, initializing with 1\n");
      alpha_x = alpha_y = 1;
    }
    *x = (1-2.*RandomSys())/alpha_x;
    *y = (1-2.*RandomSys())/alpha_y;
  }
  else if(boundary == ONE_BOUNDARY_CONDITION) {
    if(alpha_y == 0) Error("  set alpha_y different from zero\n");
    *x = L * RandomSys();
    *y = (1-2. * RandomSys())*0.5/alpha_y;
  }
  else if(boundary == TWO_BOUNDARY_CONDITIONS) {
    *x = L * RandomSys();
    *y = L * RandomSys();
  }
#endif

#ifdef TRIAL_1D //  i.e. 1D simulation
  if(boundary == NO_BOUNDARY_CONDITIONS) {
    if(alpha_z == 0) Error("  set alpha_z different from zero\n");
    *z = 1-2.*RandomSys()/alpha_z;
  }
  else if(boundary == ONE_BOUNDARY_CONDITION) { // there is a restriction in the z direction
    *z = L * RandomSys();
  }
#endif

#ifdef CRYSTAL
  if(boundary == THREE_BOUNDARY_CONDITIONS) {
    *x = Crystal.x[i] + (1-2.*RandomSys())*0.1*pow(L,1./DIMENSION);
    *y = Crystal.y[i] + (1-2.*RandomSys())*0.1*pow(L,1./DIMENSION);
    *z = Crystal.z[i] + (1-2.*RandomSys())*0.1*pow(L,1./DIMENSION);
  }
  else if(boundary == TWO_BOUNDARY_CONDITIONS) {
    *x = Crystal.x[i] + (1-2.*RandomSys())*L/(10*N);
    *y = Crystal.y[i] + (1-2.*RandomSys())*L/(10*N);
  }
  else if(boundary == ONE_BOUNDARY_CONDITION) {
    *z = Crystal.x[i] + (1-2.*RandomSys())*L/(10*N);
   }
#endif

  if(generate_new_coordinates_from_lattice) {
    if(boundary == THREE_BOUNDARY_CONDITIONS) {
      Nmax = ceil(pow((DOUBLE)N, 1./3.)); // find an integer such that Nmax^3 >= N
      l = pow(n, -1./3.); 
      for(ix=0; ix<Nmax && search; ix++)
        for(iy=0; iy<Nmax && search; iy++)
          for(iz=0; iz<Nmax && search; iz++)
            if(index == i)
              search = OFF;
            else
              index++;
      *x = l*((DOUBLE)ix-0.5) + (1-2.*RandomSys())*0.1*pow(L,1./DIMENSION);
      *y = l*((DOUBLE)iy-0.5) + (1-2.*RandomSys())*0.1*pow(L,1./DIMENSION);
      *z = l*((DOUBLE)iz-0.5) + (1-2.*RandomSys())*0.1*pow(L,1./DIMENSION);
    }
    else if(boundary == TWO_BOUNDARY_CONDITIONS) {
      Nmax = ceil(pow((DOUBLE)N, 1./2.)); // find an integer such that Nmax^2 >= N
      l = pow(n, -1./2.); 
      for(ix=0; ix<Nmax && search; ix++)
        for(iy=0; iy<Nmax && search; iy++)
          if(index == i)
            search = OFF;
          else
            index++;
      *x = l*((DOUBLE)ix-0.5) + (1-2.*RandomSys())*0.1*pow(L,1./DIMENSION);
      *y = l*((DOUBLE)iy-0.5) + (1-2.*RandomSys())*0.1*pow(L,1./DIMENSION);
    }
    else if(boundary == ONE_BOUNDARY_CONDITION) {
      *z = pow(n, -1.) * (0.5 + (DOUBLE) i)  + (1-2.*RandomSys())*L/(10*N);
    }
  }
 
#ifdef VEXT_COS2
#ifdef BC_3DPBC_CUBE
    //Message("%i  ", i);
   // 0<= i < Nlattice^3
    nx = i/(Nlattice*Nlattice);
    i -= nx*Nlattice*Nlattice;
    // 0<= i < Nlattice^2
    ny = i/Nlattice;
    i -= ny*Nlattice;
    // 0<= i < Nlattice
    nz = i;
    //Message("%i %i %i\n", nx, ny, nz);
    *x = PI*(0.51+nx)/kL;
    *y = PI*(0.51+ny)/kL;
    *z = PI*(0.51+nz)/kL;
#endif
#ifdef BC_2DPBC_SQUARE
    // 0<= i < Nlattice^2
    nx = i/(Nlattice*Nlattice);
    i -= ny*Nlattice*Nlattice;
    // 0<= i < Nlattice
    ny = i;

    *x = PI*(0.5+nx)/kL;
    *y = PI*(0.5+ny)/kL;
    *z = 0.;
#endif
#ifdef BC_1DPBC_Z
    *x = 0.;
    *y = 0.;
    *z = PI*(0.5+i)/kL;
#endif
#endif

  ReduceToTheBoxXYZ(x, y, z);
}

/***************************** Save Coordinates **************************/
void GenSaveParticles(void) {
  int w, i;
  FILE *out;
   
  out = fopen(file_particles, "w");
  fprintf(out, "%i\n", Nwalkers);
  for(w=0; w<Nwalkers; w++) {
    fprintf(out, "0.\n");
    for(i=0; i<N; i++) {
      fprintf(out, "%" LE "  %" LE "  %" LE " \n", W[w].x[i], W[w].y[i], W[w].z[i]);
    }
  }
  fclose(out);
}
