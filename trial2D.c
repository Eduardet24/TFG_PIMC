/*trial2D.c*/
/* best trial wavefunction Psi_T = Prod_{i<j} f(r_ij)
   U = ln f
   Fp = f' / f
3D Eloc = [-(f"+2f'/r)/f] + (f'/f)^2
2D Eloc = [-(f"/f + f'/r)/f] + (f'/f)^2
1D Eloc = [-f"/f] + (f'/f)^2 */

#include "trial2D.h"
#include "main.h"
#include "utils.h"
#include "memory.h"
#include <stdio.h>
#include "compatab.h"
#include MATHINCLUDE 

DOUBLE Uwell;

/************************** Load Trial Wave Function 2D ********************/
void LoadTrialWaveFunction2D(struct Grid2D *G) {
  FILE *in;
  DOUBLE R1, R2, a_check, dummy;
  int sizex1, sizex2;
  int i, j;

  Message("Loading trial wavefunction (2D)... ");

  in = fopen("wf2D.in", "r");

  if(in == NULL) {
    perror("\nError: can't load trial wave function from file,\nexiting ...");
    Error("\nError: can't load " INPATH "%s file,\nexiting ...", "wf2D.in");
  }

  fscanf(in, "a= %" LF " V= %" LF "\n", &a_check, &Uwell); // load a, Uwell
  if(a_check != a) Error("Can't continue: a = %" LG "  although in 'inwf.in' a = %" LG "\n", a, a_check);

  fscanf(in, "N1= %i R1= %" LF "\n", &sizex1, &R1); // load N1, R1
  fscanf(in, "N2= %i R2= %" LF "\n", &sizex2, &R2); // load N2, R2

  //AllocateWFGrid2D(G, sizex1, sizex2);
  G->sizex1 = sizex1;
  G->sizex2 = sizex2;

  fscanf(in, "x1=\n"); // load x1
  for(i=0; i<sizex1; i++) fscanf(in, "%" LF, &G->x1[i]);

  fscanf(in, "x2=\n"); // load x2
  for(i=0; i<sizex2; i++) fscanf(in, "%" LF, &G->x2[i]);

  fscanf(in, "f=\n"); // load lnf
  for(i=0; i<sizex1; i++) 
    for(j=0; j<sizex2; j++)
      fscanf(in, "%" LF, &G->lnf[i][j]);

  fscanf(in, "f=\n"); // skip f
  for(i=0; i<sizex1; i++) 
    for(j=0; j<sizex2; j++)
      fscanf(in, "%" LF, &dummy);
  
  /*fscanf(in, "fp_f=\n"); // load fp_f
  for(i=0; i<sizex1; i++) 
    for(j=0; j<sizex2; j++)
      fscanf(in, "%" LF, &G->fp[i][j]);

  fscanf(in, "Eloc=\n"); // load Eloc
  for(i=0; i<sizex1; i++) 
    for(j=0; j<sizex2; j++)
      fscanf(in, "%" LF, &G->E[i][j]);*/

  fclose(in);

  G->minx1 = R1 / (DOUBLE) sizex1;
  G->maxx1 = R1;
  G->stepx1 = G->maxx1 / (DOUBLE) sizex1;
  G->I_stepx1 = 1. / G->stepx1;
  G->maxx12 = G->maxx1*G->maxx1;

  G->minx2 = 0;
  G->maxx2 = R2;
  G->stepx2 =  G->maxx2 / (DOUBLE) (sizex2-1);
  G->I_stepx2 = 1. / G->stepx2;
  G->maxx22 = G->maxx2*G->maxx2;

  Message("done\n");
}

/***************************** 2D Interpolation ******************************/
DOUBLE Interpolate2DU(struct Grid2D *Grid2D, DOUBLE x1, DOUBLE x2) {
  DOUBLE I_dx1, I_dx2;
  DOUBLE t, u;
  int i, j;

  I_dx1 = Grid2D->I_stepx1;
  I_dx2 = Grid2D->I_stepx2;

  i = (int) ((x1-Grid2D->minx1) * I_dx1);
  j = (int) ((x2-Grid2D->minx2) * I_dx2);

  if(i >= Grid2D->sizex1-1) i = Grid2D->sizex1-2;
  if(j >= Grid2D->sizex2-1) j = Grid2D->sizex2-2;
  if(i < 0) i = 0;
  if(j < 0) j = 0;

  t = (x1-Grid2D->x1[i]) * I_dx1;
  u = (x2-Grid2D->x2[j]) * I_dx2;

  if(t<0) t = 0;
  if(u<0) u = 0;
  if(t>1) t = 1;
  if(u>1) u = 1;

  return (1-t)*(1-u)*Grid2D->lnf[i][j]
           + t*(1-u)*Grid2D->lnf[i+1][j]
           + t *  u *Grid2D->lnf[i+1][j+1]
           + (1-t)*u*Grid2D->lnf[i][j+1];
}
/*
  struct Grid2D G2D;

  LoadTrialWaveFunction2D(&G2D);
  for(i=0; i<500; i++) {
    for(w=0; w<500; w++) {
        Message("%" LF, Interpolate2DU(&G2D, (DOUBLE)i/200., (DOUBLE)w/200.));
    }
    Message("\n");
  }
  
  return 0;
*/
