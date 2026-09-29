/*speckle.c*/
#include <stdio.h>
#include <malloc.h>
#include "rw.h"
#include "main.h"
#include "memory.h"
#include "utils.h"
#include "vmc.h"
#include "quantities.h"
#include "trial.h"
#include "spline.h"
#include "speckles.h"
#include "randnorm.h"
#include "compatab.h"
#include "utils.h"
#include MATHINCLUDE

DOUBLE **intensity;
DOUBLE **wavefunction;
DOUBLE **wavefunction_t;
DOUBLE **V;
DOUBLE **wf;
DOUBLE *x_speckles;
DOUBLE speckles_dx, speckles_dx_inv, speckles_dx2, speckles_dx2_inv;
DOUBLE interp_dx, interp_dx_inv, interp_dx2, interp_dx2_inv;
DOUBLE **y2a,**y2ax2,**y2ax1;
int Npoints;

#ifdef POINT4
DOUBLE v[N_POINTS][N_POINTS];
DOUBLE vx[N_POINTS][N_POINTS];
DOUBLE vy[N_POINTS][N_POINTS];
DOUBLE vxx[N_POINTS][N_POINTS];
DOUBLE vxy[N_POINTS][N_POINTS];
DOUBLE vyy[N_POINTS][N_POINTS];
DOUBLE vxxx[N_POINTS][N_POINTS];
DOUBLE vxxy[N_POINTS][N_POINTS];
DOUBLE vxyy[N_POINTS][N_POINTS];
DOUBLE vyyy[N_POINTS][N_POINTS];
DOUBLE vxxxy[N_POINTS][N_POINTS];
DOUBLE vxyyy[N_POINTS][N_POINTS];

DOUBLE w[N_POINTS][N_POINTS];
DOUBLE wx[N_POINTS][N_POINTS];
DOUBLE wy[N_POINTS][N_POINTS];
DOUBLE wxx[N_POINTS][N_POINTS];
DOUBLE wxy[N_POINTS][N_POINTS];
DOUBLE wyy[N_POINTS][N_POINTS];
DOUBLE wxxx[N_POINTS][N_POINTS];
DOUBLE wxxy[N_POINTS][N_POINTS];
DOUBLE wxyy[N_POINTS][N_POINTS];
DOUBLE wyyy[N_POINTS][N_POINTS];
DOUBLE wxxxy[N_POINTS][N_POINTS];
DOUBLE wxyyy[N_POINTS][N_POINTS];
#else
DOUBLE **v,**vx,**vy,**vxx,**vxy,**vyy,**vxxx,**vxxy,**vxyy,**vyyy,**vxxxy,**vxyyy;
DOUBLE **w,**wx,**wy,**wxx,**wxy,**wyy,**wxxx,**wxxy,**wxyy,**wyyy,**wxxxy,**wxyyy;
#endif

/************************** Speckles Load *********************************/
void SpecklesLoad(void) {
  FILE *in, *out;
  DOUBLE dummy1, dummy2;
  int i,j;
  DOUBLE x1,x2,Fx,Fy,Fx2,Fy2;
  DOUBLE Lmax, Lmin;
  DOUBLE Integral;
  int Nmax;
  int factor=1;

  x_speckles = vector(0, N_POINTS+5);
  speckles_dx = Lx/(DOUBLE)(N_POINTS-1);
  speckles_dx_inv = 1./speckles_dx;
  speckles_dx2 = speckles_dx*speckles_dx;
  speckles_dx2_inv = 1./speckles_dx2;

  for(i=1; i<=N_POINTS; i++) x_speckles[i] = speckles_dx*(DOUBLE)(i-1);

  Message("  Size of the speckle matrix %ix%i\n", N_POINTS, N_POINTS);
  Message("  Loading speckle pattern ... ");
  in = fopen("speckle.in","r");
  if(in == NULL) Error("\nError: Cannot open file ``speckle.in''.");
  intensity = matrix(1, N_POINTS, 1, N_POINTS);
  for(i=1; i<=N_POINTS; i++) {
    for(j=1; j<=N_POINTS; j++) {
      fscanf(in, "%" LF " %" LF " %" LF  "\n", &dummy1, &dummy2, &intensity[i][j]);
      if(j==1 && dummy2>1e-5) 
        Error("  wrong speckle size. Expected size is %ix%i, check number of lines in 'speckle.in'\n", N_POINTS, N_POINTS);
       //intensity[i][j] = x_speckles[i]-Lhalf; // linear
       //intensity[i][j] = x_speckles[j]-Lhalf; // linear
       //intensity[i][j] = (x_speckles[i]-Lhalf)*(x_speckles[i]-Lhalf); // square
       //intensity[i][j] = (x_speckles[j]-Lhalf)*(x_speckles[j]-Lhalf); // square
       //intensity[i][j] = (x_speckles[i]-Lhalf)*(x_speckles[i]-Lhalf)*(x_speckles[i]-Lhalf); // cubic
       //intensity[i][j] = ((x_speckles[i]-Lhalf)*(x_speckles[i]-Lhalf)+(x_speckles[j]-Lhalf)*(x_speckles[j]-Lhalf)); // parabolic
    }
  }
  fclose(in);
  Message("done\n");

  Message("  Loading one-body wave function ... ");
  wavefunction_t = matrix(1, N_POINTS, 1, N_POINTS);
  wavefunction = matrix(1, N_POINTS, 1, N_POINTS);
  in = fopen("density.in","r");
  if(in == NULL) Error("\nError: Cannot open the file ``density.in''.");
  for(i=1; i<=N_POINTS; i++) {
    for(j=1; j<=N_POINTS; j++) {
      fscanf(in, "%" LF " %" LF " %" LF  "\n", &dummy1, &dummy2, &wavefunction[i][j]);
      //wavefunction[i][j] = Exp(-10*(((x_speckles[i]-Lhalf)*(x_speckles[i]-Lhalf)+(x_speckles[j]-Lhalf)*(x_speckles[j]-Lhalf)))); // HO
      //wavefunction[i][j] = Exp(x_speckles[i]-Lhalf); // linear
      //wavefunction[i][j] = Exp((x_speckles[i]-Lhalf)*(x_speckles[i]-Lhalf)); // square
      //wavefunction[i][j] = Exp(Sin(10*(x_speckles[i]-Lhalf))); // sine
      //wavefunction[i][j] = 1.;
      //wavefunction[i][j] = i;
    }
  }
  fclose(in);
  Message("done\n");

  for(i=1; i<=N_POINTS; i++) { // convert density to w.f.
    for(j=1; j<=N_POINTS; j++) {
      wavefunction[i][j] = Sqrt(wavefunction[i][j]);
    }
  }

#ifdef INTERPOLATE_LOG_ONE_BODY
  Message("  Log of one-body w.f. will be  used\n");
  if(inwf_power != 1) Warning("power of one-body w.f. is changed to %lf !\n", inwf_power);// adjust power
  for(i=0; i<N_POINTS; i++) {
    for(j=0; j<N_POINTS; j++) {
      if(wavefunction[i+1][j+1]<0) Error("  cannot take log of a negative speckle w.f.!\n");
      wavefunction[i+1][j+1] = Log(wavefunction[i+1][j+1]);
      wavefunction[i+1][j+1] *= inwf_power;
    }
  }
#endif

  Message("\n  checking w.f. normalization, matrix (%ix%i) dx = %lf, L = %lf\n", N_POINTS, N_POINTS, Lx/(DOUBLE) (N_POINTS-1), Lx);
  Integral = 0.;
  for(i=1; i<N_POINTS; i++) {
    x1 = ((DOUBLE) i)/(DOUBLE) (N_POINTS-1) * Lx;
    for(j=1; j<N_POINTS; j++) {
      x2 = ((DOUBLE) j)/(DOUBLE) (N_POINTS-1) * Lx;
      Integral += exp(2.*wavefunction[i][j]);
    }
  }
  Integral *= (Lx/(DOUBLE)(N_POINTS-1))*(Lx/(DOUBLE)(N_POINTS-1));
  Message("  Checking 2D normalization of the spline w.f.: %lf\n", Integral);
  for(i=1; i<=N_POINTS; i++) {
    for(j=1; j<=N_POINTS; j++) {
      wavefunction[i][j] -= 0.5*log(Integral);
    }
  }
  Integral = 0.;
  for(i=1; i<N_POINTS; i++) {
    x1 = ((DOUBLE) i)/(DOUBLE) (N_POINTS-1) * Lx;
    for(j=1; j<N_POINTS; j++) {
      x2 = ((DOUBLE) j)/(DOUBLE) (N_POINTS-1) * Lx;
      Integral += exp(2.*wavefunction[i][j]);
    }
  }
  Integral *= (Lx/(DOUBLE)(N_POINTS-1))*(Lx/(DOUBLE)(N_POINTS-1));
  Message("  corrected: Checking 2D normalization of the spline w.f.: %lf\n", Integral);

  for(i=1; i<=N_POINTS; i++) { // invert indexes
    for(j=1; j<=N_POINTS; j++) {
      wavefunction_t[i][j] = wavefunction[j][i];
    }
  }

  //preliminary steps for cubic-spline interpolation
  Message("  Constructing 2D cubic splines ... ");
  y2a=matrix(1, N_POINTS, 1, N_POINTS);
  splie2(x_speckles, x_speckles, intensity, N_POINTS, N_POINTS, y2a);
  y2ax2=matrix(1, N_POINTS, 1, N_POINTS);
  splie2(x_speckles, x_speckles, wavefunction, N_POINTS, N_POINTS, y2ax2);
  y2ax1 = matrix(1, N_POINTS, 1, N_POINTS);
  splie2(x_speckles, x_speckles, wavefunction_t, N_POINTS, N_POINTS, y2ax1);
  Message("done\n");

#ifndef SPLINES_2D
#ifdef POINT9
  Message("  Fast 2D interpolation using 9 nearest neighbors will be used\n");
#endif
#ifdef POINT16
  Message("  Fast 2D interpolation using 16 nearest neighbors will be used\n  Cubic polynom is used.\n");
#endif
#ifdef POINT16CONTINUOUS
  Message("  Fast 2D interpolation using 16 nearest neighbors will be used\n  Polynom of fifth power is used\n");
#endif
#ifdef POINT4
  Message("  Fast 2D interpolation using 4 nearest neighbors and 2nd derivatives\n  Cubic polynom is used.\n");
#endif

  /*factor = 1;
  Message("  Increasing accuracy of interpolation, resampling on a fine grid using splines\n");
  Message("  New spacing is different by a factor of 1/%i. Be patient ...\n", factor);
  Npoints = factor*N_POINTS;
  V = matrix(0, Npoints+5, 0, Npoints+5); // adding 5 for PBC
  wf = matrix(0, Npoints+5, 1, Npoints+5);
  Message("  ................\n  ");
  for(i=0; i<Npoints; i++) {
    x1 = ((DOUBLE) i)/(DOUBLE) (Npoints-1) * Lx;
    for(j=0; j<Npoints; j++) {
      x2 = ((DOUBLE) j)/(DOUBLE) (Npoints-1) * Lx;
      if(factor == 1) {
#ifndef INTERPOLATE_LOG_ONE_BODY
        if(wavefunction[i+1][j+1]<0) Warning("  negative w.f.\n");
        if(wavefunction[i+1][j+1]==0) Warning("  zero w.f.\n");
#endif
        wf[i+1][j+1] = wavefunction[i+1][j+1];
        V[i+1][j+1] = intensity[i+1][j+1];
      }
      else {
        wf[i+1][j+1] = SpecklesUSplines(x1, x2);
        V[i+1][j+1] = SpecklesEpotSplines(x1, x2);
      }
    }
    //Message("%i/%i\n", i, Npoints);
    if(i%(Npoints/16)== 0) Message(".");
  }*/

  Npoints = N_POINTS;
  V = matrix(0, Npoints+5, 0, Npoints+5); // adding 5 for PBC
  wf = matrix(0, Npoints+5, 1, Npoints+5);
  for(i=0; i<Npoints; i++) {
    for(j=0; j<Npoints; j++) {
#ifndef INTERPOLATE_LOG_ONE_BODY
        if(wavefunction[i+1][j+1]<0) Warning("  negative w.f.\n");
        if(wavefunction[i+1][j+1]==0) Warning("  zero w.f.\n");
#endif
      wf[i+1][j+1] = wavefunction[i+1][j+1];
      V[i+1][j+1] = intensity[i+1][j+1];
    }
  }

  interp_dx = Lx/(DOUBLE)(Npoints-1);
  interp_dx_inv = 1./interp_dx;
  interp_dx2 = interp_dx*interp_dx;
  interp_dx2_inv = 1./interp_dx2;

  // adding virtual points according to PBC
  // the last column is already repeated in the source file
  // i.e. columns [1,Npoints-1] are unique
  // column [0] is the same as [Npoints-1]
#ifdef POINT16 // next 5 points to the right are the same
  for(i=0; i<Npoints; i++) { // adding virtual points according to PBC
    for(j=0; j<5; j++) {
      V[Npoints+j+1][i+1] = V[j+1][i+1];
      V[i+1][Npoints+j+1] = V[i+1][j+1];
      if(i<5) V[Npoints+i+1][Npoints+j+1] = V[i+1][j+1];
    }
  }

  for(i=0; i<Npoints; i++) { // adding virtual points according to PBC
    for(j=0; j<5; j++) {
      wf[Npoints+j+1][i+1] = wf[j+1][i+1];
      wf[i+1][Npoints+j+1] = wf[i+1][j+1];
      if(i<5) wf[Npoints+i+1][Npoints+j+1] = wf[i+1][j+1];
    }
  }

  for(j=0; j<5; j++) x_speckles[Npoints+j+1] = x_speckles[j+1];
#else // duplicate for POINT4 case
  for(i=1; i<=Npoints; i++) {
    V[0][i] = V[Npoints][i];
    //V[Npoints+1][i] = V[1][i];
    V[i][0] = V[i][Npoints];
    //V[i][Npoints+1] = V[i][1];
  }
#endif

#ifdef POINT4
  SpecklesInitialize4();
#endif

  Message("\n  Saving speckle difference... ");
  out = fopen("speckle_diff.dat","w");
  fprintf(out, "x y Vpot lnf Fx Fy Fpp\n");
  if(out == NULL) Error("\nError: Cannot open the file.");
  //Lmax = Lx;
  Lmin = Lx/2;
  Lmax = Lmin + interp_dx*3.;
  Nmax = 31;
  for(i=0; i<Nmax; i++) {
    x1 = Lmin + ((DOUBLE) i+1.)/(DOUBLE) (Nmax) * (Lmax-Lmin) + (Random()-0.5)*interp_dx*0;
    for(j=0; j<Nmax; j++) {
      x2 = Lmin + ((DOUBLE) j+1.)/(DOUBLE) (Nmax) * (Lmax-Lmin) + (Random()-0.5)*interp_dx*0;

      Fx = Fy = 0.;
      SpecklesFpSplines(&Fx, &Fy, x1, x2);

      Fx2 = Fy2= 0.;
      SpecklesFp(&Fx2, &Fy2, x1, x2);
      // 1-7
      fprintf(out,"%.4" LF " %.4" LG " %.6" LG " %.4" LG " %.4" LG " %.4" LG " %.4" LG, 
        (x1-Lmin)/interp_dx, // 1
        (x2-Lmin)/interp_dx,  // 2
          SpecklesEpot(x1,x2)-SpecklesEpotSplines(x1,x2), // 3
          SpecklesU(x1,x2)-SpecklesUSplines(x1,x2), // 4
          Fx2-Fx, // 5
          Fy2-Fy, // 6
          SpecklesFpp(x1,x2)-SpecklesFppSplines(x1,x2)); // 7
      fprintf(out," %.6" LG " %.4" LG " %.4" LG " %.4" LG " %.4" LG, 
        SpecklesEpot(x1,x2), // 8
        SpecklesU(x1,x2), // 9 
        Fx2, // 10
        Fy2, // 11
        SpecklesFpp(x1,x2)); // 12
      fprintf(out,"\n");
    }
    fprintf(out,"\n");
  }
  fclose(out);
#endif

  Message("\n  Saving speckle w.f. and its derivatives of w.f. ... ");
  out = fopen("speckle.dat","w");
  fprintf(out, "x y Vpot lnf Fx Fy Fpp\n");
  if(out == NULL) Error("\nError: Cannot open the file.");
  Lmax = Lx;
  Nmax = 100;
  for(i=0; i<Nmax; i++) {
    x1 = ((DOUBLE) i+1.)/(DOUBLE) (Nmax) * Lmax+(Random()-0.5)*interp_dx;
    for(j=0; j<Nmax; j++) {
      x2 = ((DOUBLE) j+1.)/(DOUBLE) (Nmax) * Lmax+(Random()-0.5)*interp_dx;

      Fx = Fy = 0.;
      SpecklesFp(&Fx, &Fy, x1, x2);

      fprintf(out,"%.16" LE " %.16" LE " %.16" LE " %.4" LG " %.4" LG " %.4" LG " %.4" LG "\n", x1, x2, SpecklesEpot(x1,x2), SpecklesU(x1,x2), Fx, Fy, SpecklesFpp(x1,x2));
    }
    fprintf(out,"\n");
  }
  fclose(out);

//   for(factor=1; factor<4; factor++) {
//     Integral = 0.;
//     for(i=0; i<factor*(N_POINTS-1); i++) {
//       x1 = ((DOUBLE) i)/(DOUBLE) (factor*N_POINTS-1) * Lx;
//       for(j=0; j<factor*(N_POINTS-1); j++) {
//         x2 = ((DOUBLE) j)/(DOUBLE) (factor*N_POINTS-1) * Lx;
//         Integral += exp(2.*SpecklesUSplines(x1, x2));
//       }
//       if(i%(Npoints/16)== 0) Message(".");
//     }
//     Integral *= (Lx/(DOUBLE)(factor*N_POINTS-1))*(Lx/(DOUBLE)(factor*N_POINTS-1));
//     Message("%i %lf\n", factor, Integral);
//   }

  Message("done\n");
}

/************************** Speckles Epot *********************************/
// Speckle potential
DOUBLE SpecklesEpotSplines(DOUBLE x1, DOUBLE x2) {
  DOUBLE y;

  splin2(x_speckles, x_speckles, intensity, y2a, N_POINTS, N_POINTS, x1, x2, &y);

  return y;
}

/************************** Speckles U ************************************/
// logarithm of one-body w.f.
// ln f
DOUBLE SpecklesUSplines(DOUBLE x1, DOUBLE x2) {
  DOUBLE y;

  splin2(x_speckles, x_speckles, wavefunction, y2ax2, N_POINTS, N_POINTS, x1, x2, &y);

#ifdef INTERPOLATE_LOG_ONE_BODY
  return y;
#else
  return Log(y);
#endif
}

/************************** Speckles Fp ***********************************/
// logarithmic derivative 
// f'/f = (fx/f, fy/f)
void SpecklesFpSplines(DOUBLE *Fx, DOUBLE *Fy, DOUBLE x1, DOUBLE x2) {
  DOUBLE y,f1,f2,f11,f22;

  splin2x1(x_speckles, x_speckles, wavefunction, y2ax2, N_POINTS, N_POINTS, x1, x2, &y, &f1, &f11);
  splin2x1(x_speckles, x_speckles, wavefunction_t, y2ax1, N_POINTS, N_POINTS, x2, x1, &y, &f2, &f22);

#ifndef INTERPOLATE_LOG_ONE_BODY
  f1 /= y;
  f2 /= y;
#endif

  //first derivative with respect to x1 divided by the wavefunction
  *Fx += f1;
  //first derivative with respect to x2 divided by the wavefunction
  *Fy += f2;
}

/************************** Speckles Fpp **********************************/
// One-body contribution to the kinetic energy
// -f''/f + (f'/f)^2
DOUBLE SpecklesFppSplines(DOUBLE x1, DOUBLE x2) {
  DOUBLE y,f1,f11,f2,f22;

  splin2x1(x_speckles, x_speckles, wavefunction, y2ax2, N_POINTS, N_POINTS, x1, x2, &y, &f1, &f11);
  splin2x1(x_speckles, x_speckles, wavefunction_t, y2ax1, N_POINTS, N_POINTS, x2, x1, &y, &f2, &f22);

#ifndef INTERPOLATE_LOG_ONE_BODY
  f1 /= y;
  f2 /= y;
  f11 /= y;
  f22 /= y;
  return - f11 - f22 + f1*f1 + f2*f2;
#else
  return - f11 - f22;
#endif

}

/**************************** fast interpolation ***************************/
// Speckle potential
DOUBLE SpecklesEpot9point(DOUBLE x, DOUBLE y) {
  int x0,x1,xm1;
  int y0,y1,ym1;
  DOUBLE dx,dy,dx2,dy2,dx2half,dy2half,dxhalf,dyhalf;

  x0 = (int) (x*interp_dx_inv-0.5) + 1;
  if(x0 == Npoints) x0 = Npoints-1;
  x1 = x0+1;
  xm1 = x0-1;
  if(xm1<1) {
    xm1 = 1;
    x0 = 2;
    x1 = 3;
  }
  y0 = (int) (y*interp_dx_inv-0.5) + 1;
  if(y0 == Npoints) y0 = Npoints-1;
  y1 = y0+1;
  ym1 = y0-1;
  if(ym1<1) {
    ym1 = 1;
    y0 = 2;
    y1 = 3;
  }

#ifdef SECURE
  if(xm1<1) Warning("  splint call: xm1 out of range\n");
  if(x1>Npoints)  Warning("  splint call: x1 out of range (Epot)\n");
  if(ym1<1) Warning("  splint call: ym1 out of range (Epot)\n");
  if(y1>Npoints)  Warning("  splint call: y1 out of range\n");
#endif

  dx = x - interp_dx*x0;
  dy = y - interp_dx*y0;
  dxhalf = 0.5*dx;
  dx2 = dx*dx;
  dyhalf = 0.5*dy;
  dy2 = dy*dy;
  dx2half = 0.5*dx2;
  dy2half = 0.5*dy2;

  return
    V[ x0][ y0]*(1.-dx2-dy2)
  + V[ x1][ y0]*(dxhalf+dx2half)
  + V[xm1][ y0]*(-dxhalf+dx2half)
  + V[ x0][ y1]*(dyhalf+dy2half)
  + V[ x0][ym1]*(-dyhalf+dy2half)
  + dxhalf*dyhalf*(V[x1][y1]-V[x1][ym1]-V[xm1][y1]+V[xm1][ym1]);
}

/************************** Speckles U ************************************/
// logarithm of one-body w.f.
// ln f
DOUBLE SpecklesU9point(DOUBLE x, DOUBLE y) {
  int x0,x1,xm1;
  int y0,y1,ym1;
  DOUBLE dx,dy,dx2,dy2,dx2half,dy2half,dxhalf,dyhalf;


#ifndef INTERPOLATE_LOG_ONE_BODY
  Error("enable INTERPOLATE_LOG_ONE_BODY\n");
#endif

  x0 = (int) (x*interp_dx_inv-0.5) + 1;
  if(x0 == Npoints) x0 = Npoints-1;
  x1 = x0+1;
  xm1 = x0-1;
  if(xm1<1) {
    xm1 = 1;
    x0 = 2;
    x1 = 3;
  }
  y0 = (int) (y*interp_dx_inv-0.5) + 1;
  if(y0 == Npoints) y0 = Npoints-1;
  y1 = y0+1;
  ym1 = y0-1;
  if(ym1<1) {
    ym1 = 1;
    y0 = 2;
    y1 = 3;
  }
  dx = x - interp_dx*x0;
  dy = y - interp_dx*y0;

#ifdef SECURE
  if(xm1<1) Warning("  splint call: xm1 out of range\n");
  if(x1>Npoints) Warning("  splint call: x1 out of range (U)\n");
  if(ym1<1) Warning("  splint call: ym1 out of range (U)\n");
  if(y1>Npoints)  Warning("  splint call: y1 out of range\n");
  if(fabs(dx)>interp_dx && xm1 != 1 && x1 != Npoints) Warning("  dx>h (U)\n");
  if(fabs(dy)>interp_dx && ym1 != 1 && y1 != Npoints) Warning("  dy>h (U)\n");
#endif

  dxhalf = 0.5*dx;
  dx2 = dx*dx;
  dyhalf = 0.5*dy;
  dy2 = dy*dy;
  dx2half = 0.5*dx2;
  dy2half = 0.5*dy2;

  return
    wf[ x0][ y0]*(1.-dx2-dy2)
  + wf[ x1][ y0]*(dxhalf+dx2half)
  + wf[xm1][ y0]*(-dxhalf+dx2half)
  + wf[ x0][ y1]*(dyhalf+dy2half)
  + wf[ x0][ym1]*(-dyhalf+dy2half)
  + dxhalf*dyhalf*(wf[x1][y1]-wf[x1][ym1]-wf[xm1][y1]+wf[xm1][ym1]);
}

/************************** Speckles Fp ***********************************/
// logarithmic derivative 
// u' = (ux, uy)
void SpecklesFp9point(DOUBLE *Fx, DOUBLE *Fy, DOUBLE x, DOUBLE y) {
  int x0,x1,xm1;
  int y0,y1,ym1;
  DOUBLE dx,dy;

  x0 = (int) (x*interp_dx_inv-0.5) + 1;
  if(x0 == Npoints) x0 = Npoints-1;
  x1 = x0+1;
  xm1 = x0-1;
  if(xm1<1) {
    xm1 = 1;
    x0 = 2;
    x1 = 3;
  }
  y0 = (int) (y*interp_dx_inv-0.5) + 1;
  if(y0 == Npoints) y0 = Npoints-1;
  y1 = y0+1;
  ym1 = y0-1;
  if(ym1<1) {
    ym1 = 1;
    y0 = 2;
    y1 = 3;
  }
  dx = x - interp_dx*x0;
  dy = y - interp_dx*y0;

#ifdef SECURE
  if(xm1<1) Warning("  splint call: xm1 out of range\n");
  if(x1>Npoints) Warning("  splint call: x1 out of range (Fp)\n");
  if(ym1<1) Warning("  splint call: ym1 out of range (Fp)\n");
  if(y1>Npoints)  Warning("  splint call: y1 out of range\n");
  if(fabs(dx)>interp_dx && xm1 != 1 && x1 != Npoints) Warning("  dx>h (Fp)\n");
  if(fabs(dy)>interp_dx && ym1 != 1 && y1 != Npoints) Warning("  dy>h (Fp)\n");
#endif

#ifndef INTERPOLATE_LOG_ONE_BODY
  Error("enable INTERPOLATE_LOG_ONE_BODY\n");
#endif

  *Fx = (wf[x0][y0]*(-2.*dx) + wf[x1][y0]*(dx+0.5) + wf[xm1][y0]*(dx-0.5)
  + 0.25*dy*(wf[x1][y1]-wf[x1][ym1]-wf[xm1][y1]+wf[xm1][ym1]))*interp_dx_inv;
  *Fy = (wf[x0][y0]*(-2.*dy) + wf[x0][y1]*(dy+0.5) + wf[x0][ym1]*(dy-0.5)
  + 0.25*dx*(wf[x1][y1]-wf[xm1][y1]-wf[x1][ym1]+wf[xm1][ym1]))*interp_dx_inv;
}

/************************** Speckles Fpp **********************************/
// One-body contribution to the kinetic energy
// -f''/f + (f'/f)^2
// -u''+ (u')^2
DOUBLE SpecklesFpp9point(DOUBLE x, DOUBLE y) {
  int x0,x1,xm1;
  int y0,y1,ym1;
  DOUBLE dx,dy;
  DOUBLE f1,f11,f2,f22;

  x0 = (int) (x*interp_dx_inv-0.5) + 1;
  if(x0 == Npoints) x0 = Npoints-1;
  x1 = x0+1;
  xm1 = x0-1;
  if(xm1<1) {
    xm1 = 1;
    x0 = 2;
    x1 = 3;
  }
  y0 = (int) (y*interp_dx_inv-0.5) + 1;
  if(y0 == Npoints) y0 = Npoints-1;
  y1 = y0+1;
  ym1 = y0-1;
  if(ym1<1) {
    ym1 = 1;
    y0 = 2;
    y1 = 3;
  }
  dx = x - interp_dx*x0;
  dy = y - interp_dx*y0;

#ifdef SECURE
  if(xm1<1) Warning("  splint call: xm1 out of range\n");
  if(x1>Npoints) Warning("  splint call: x1 out of range (Fpp)\n");
  if(ym1<1) Warning("  splint call: ym1 out of range (Fpp)\n");
  if(y1>Npoints)  Warning("  splint call: y1 out of range\n");
  if(fabs(dx)>interp_dx && xm1 != 1 && x1 != Npoints) Warning("  dx>h (Fpp)\n");
  if(fabs(dy)>interp_dx && ym1 != 1 && y1 != Npoints) Warning("  dy>h (Fpp)\n");
#endif

  f1 = wf[x0][y0]*(-2.*dx) + wf[x1][y0]*(dx+0.5) + wf[xm1][y0]*(dx-0.5)
  + 0.25*dy*(wf[x1][y1]-wf[x1][ym1]-wf[xm1][y1]+wf[xm1][ym1]);
  f2 = wf[x0][y0]*(-2.*dy) + wf[x0][y1]*(dy+0.5) + wf[x0][ym1]*(dy-0.5)
  + 0.25*dx*(wf[x1][y1]-wf[xm1][y1]-wf[x1][ym1]+wf[xm1][ym1]);
  f11 = wf[x1][y0]+wf[xm1][y0]-2.*wf[x0][y0];
  f22 = wf[x0][y1]+wf[x0][ym1]-2.*wf[x0][y0];

  //return (-f11 -f22 + f1*f1 + f2*f2)*interp_dx2_inv;
  return (-f11 - f22)*interp_dx2_inv;
}

/*******************************************************************************
Given an m by n tabulated function ya[1..m][1..n],
and tabulated independent variables x2a[1..n],
this routine constructs one-dimensional natural cubic splines of the rows of ya
and returns the second-derivatives in the array y2a[1..m][1..n].
(The array x1a[1..m] is included in the argument list merely for consistency
 with routine splin2.)
*******************************************************************************/
void splie2(DOUBLE x1a[], DOUBLE x2a[], DOUBLE **ya, int m, int n, DOUBLE **y2a) {
  int j;

  for(j=1 ; j <= m ; j++ ) spline(x2a,ya[j],n,1.0e30,1.0e30,y2a[j]); //Values 1Ã-10^30 signal a natural spline
}

/*******************************************************************************
Given an m by n tabulated function ya[1..m][1..n],
and tabulated independent variables x2a[1..n],
this routine constructs one-dimensional natural cubic splines of the rows of ya
and returns the second-derivatives in the array y2a[1..m][1..n].
(The array x1a[1..m] is included in the argument list merely for consistency
 with routine splin2.)
*******************************************************************************/
void splie2x1(DOUBLE x1a[], DOUBLE x2a[], DOUBLE **ya, int m, int n, DOUBLE **y2a) {
  DOUBLE *yax1,*y2;
  int i, j;

  yax1=vector(1,m);
  y2=vector(1,m);

  for(j=1; j<=n ; j++) {
    for(i=1; i<=m ; i++) yax1[i]=ya[i][j];

    spline(x1a,yax1,m,1.0e30,1.0e30,y2); //Values 1Ã-10^30 signal a natural spline

    for(i=1; i<=m ; i++) y2a[i][j]=y2[i];
  }

  free_vector(yax1,1,m);
  free_vector(y2,1,m);
}

/*******************************************************************************
Given x1a, x2a, ya, m, n as described in splie2 and y2a
as produced by that routine;
and given a desired interpolating point x1,x2;
this routine returns an interpolated function value y
by bicubic spline interpolation.
*******************************************************************************/
void splin2(DOUBLE x1a[], DOUBLE x2a[], DOUBLE **ya, DOUBLE **y2a, int m, int n, DOUBLE x1, DOUBLE x2, DOUBLE *y) {
  int j;
  DOUBLE *ytmp,*yytmp;

  ytmp=vector(1,m);
  yytmp=vector(1,m);

  //Perform m evaluations of the row splines constructed by splie2,
  //using the one-dimensional spline evaluator splint.

  for(j=1; j<=m; j++) splint(x2a,ya[j],y2a[j],n,x2,&yytmp[j]);

  spline(x1a,yytmp,m,1.0e30,1.0e30,ytmp);
  //Construct the one-dimensional column spline and evaluate it.
  splint(x1a,yytmp,ytmp,m,x1,y);

  free_vector(yytmp,1,m);
  free_vector(ytmp,1,m);
}

/*******************************************************************************
Given arrays x[1..n] and y[1..n] containing a tabulated function, i.e.,
yi = f(xi), with x1 < x2 < .. . < xN,
and given values yp1 and ypn for the first derivative of the interpolating
function at points 1 and n, respectively, this routine returns
an array y2[1..n] that contains the second derivatives of
the interpolating function at the tabulated points xi.
If yp1 and/or ypn are equal to 1 Ã- 10^30 or larger,
the routine is signaled to set the corresponding boundary
condition for a natural spline, with zero second derivative on that boundary.
*******************************************************************************/
    void spline(DOUBLE x[], DOUBLE y[], int n, DOUBLE yp1, DOUBLE ypn, DOUBLE y2[]) {
      int i,k;
      DOUBLE p,qn,un,*u;

      u=vector(1,n-1);

      if(yp1>0.99e30) //The lower boundary condition is set either to be â__naturalâ__
        y2[1]=u[1]=0.0;
      else{ //or else to have a specified first derivative.
        y2[1] = -0.5;
        //u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
        u[1]=(3.*speckles_dx_inv)*((y[2]-y[1])*speckles_dx_inv-yp1);
      }

      for(i=2; i<=n-1; i++) {
        //This is the decomposition loop of the tridiagonal algorithm.
        //y2 and u are used for temporary storage of the decomposed factors.
        //sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
        //p=sig*y2[i-1]+2.0;
        //y2[i]=(sig-1.0)/p;
        //u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
        //u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;

        p = 0.5*y2[i-1] + 2.;
        y2[i] = -0.5/p;
        //u[i] = (y[i+1]-y[i])*speckles_dx_inv - (y[i]-y[i-1])*speckles_dx_inv;
        u[i] = (y[i+1]-2.*y[i]+y[i-1])*speckles_dx_inv;
        u[i] = (3.*u[i]*speckles_dx_inv-0.5*u[i-1])/p;
      }

      if(ypn>0.99e30) //The upper boundary condition is set either to be â__naturalâ__
        qn = un = 0.;
      else { //or else to have a specified first derivative.
        qn = 0.5;
        //un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
        un = (3.*speckles_dx_inv)*(ypn-(y[n]-y[n-1])*speckles_dx_inv);
      }

      y2[n] = (un-qn*u[n-1])/(qn*y2[n-1]+1.);

      //This is the back substitution loop of the tridiagonal algorithm.
      for(k=n-1; k>=1; k-- ) y2[k]=y2[k]*y2[k+1]+u[k];

      free_vector(u,1,n-1);
    }

/*******************************************************************************
Given the arrays xa[1..n] and ya[1..n],
which tabulate a function (with the xaiâ_Ts in order),
and given the array y2a[1..n], which is the output from spline above,
and given a value of x,
this routine returns a cubic-spline interpolated value y.
*******************************************************************************/
void splint(DOUBLE xa[], DOUBLE ya[], DOUBLE y2a[], int n, DOUBLE x, DOUBLE *y) {

  DOUBLE b,a;
  int klo,khi;

  /*We will find the right place in the table by means of bisection.
  This is optimal if sequential calls to this routine are at random values of x.
  If sequential calls are in order, and closely spaced, one would do better
  to store previous values of klo and khi and test if
  they remain appropriate on the next call.*/

  klo = (int) (x*speckles_dx_inv) + 1;
  if(klo == n) klo = n-1;
  khi = klo+1;

#ifdef SECURE
  if(klo<1) {
    klo = 1;
    Warning("  splint call: klo out of range\n");
  }
  if(khi>n) {
    klo = n-1;
    khi = n;
    Warning("  splint call: khi out of range\n");
  }
#endif

  a=(xa[khi]-x)*speckles_dx_inv;
  b=(x-xa[klo])*speckles_dx_inv; //Cubic spline polynomial is now evaluated.

  *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*speckles_dx2/6.;
}

/*******************************************************************************
Given the arrays xa[1..n] and ya[1..n],
which tabulate a function (with the xaiâ_Ts in order),
and given the array y2a[1..n], which is the output from spline above,
and given a value of x,
this routine returns the first derivative y1 obtained by a cubic-spline interpolation.
*******************************************************************************/
void splder1(DOUBLE xa[], DOUBLE ya[], DOUBLE y2a[], int n, DOUBLE x, DOUBLE *y1) {

  int klo,khi;
  DOUBLE b,a;

  klo = (int) (x*speckles_dx_inv) + 1;
  if(klo == n) klo = n-1;
  khi = klo+1;

  a=(xa[khi]-x)*speckles_dx_inv;
  b=(x-xa[klo])*speckles_dx_inv;

  //First derivative is now evaluated.
  *y1 = (ya[khi]-ya[klo])*speckles_dx_inv + ((-3.*a*a+1.)*y2a[klo]+(3.*b*b-1.)*y2a[khi])*speckles_dx/6.;
}

/*******************************************************************************
Given the arrays xa[1..n] and ya[1..n],
which tabulate a function (with the xaiâ_Ts in order),
and given the array y2a[1..n], which is the output from spline above,
and given a value of x,
this routine returns the second derivative y2 obtained by a cubic-spline interpolation.
*******************************************************************************/
void splder2(DOUBLE xa[], DOUBLE ya[], DOUBLE y2a[], int n, DOUBLE x, DOUBLE *y2) {

  int klo,khi;
  DOUBLE b,a;

  klo = (int) (x*speckles_dx_inv) + 1;
  if(klo == n) klo = n-1;
  khi = klo+1;

  a=(xa[khi]-x)*speckles_dx_inv;
  b=(x-xa[klo])*speckles_dx_inv;

  //Second derivative is now evaluated.
  *y2 = a*y2a[klo]+b*y2a[khi];
}

/* allocate a DOUBLE matrix with subscript range m[nrl..nrh][ncl..nch] */
DOUBLE **matrix(long nrl, long nrh, long ncl, long nch) {
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  DOUBLE **m;

  // allocate pointers to rows 
  m=(DOUBLE **) malloc((size_t)((nrow+NR_END)*sizeof(DOUBLE*)));
  if(!m) Error("allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;

  // allocate rows and set pointers to them
  m[nrl]=(DOUBLE *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(DOUBLE)));
  if(!m[nrl]) Error("allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  return m; //return pointer to array of pointers to rows
}

/* free a DOUBLE vector allocated with vector() */
void free_vector(DOUBLE *v, long nl, long nh) {
  free((FREE_ARG) (v+nl-NR_END));
}

/* allocate a DOUBLE vector with subscript range V[nl..nh] */
DOUBLE *vector(long nl, long nh) {
  DOUBLE *v;

  v=(DOUBLE *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(DOUBLE)));
  if(!v) Error("allocation failure in vector()");
  return v-nl+NR_END;
}

/* allocate a DOUBLE 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
DOUBLE ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh) {
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  DOUBLE ***t;

  /* allocate pointers to pointers to rows */
  t=(DOUBLE ***) malloc((size_t)((nrow+NR_END)*sizeof(DOUBLE**)));
  if(!t) Error("allocation failure 1 in f3tensor()");
  t += NR_END;
  t -= nrl;

  /* allocate pointers to rows and set pointers to them */
  t[nrl]=(DOUBLE **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(DOUBLE*)));
  if(!t[nrl]) Error("allocation failure 2 in f3tensor()");
  t[nrl] += NR_END;
  t[nrl] -= ncl;

  /* allocate rows and set pointers to them */
  t[nrl][ncl]=(DOUBLE *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(DOUBLE)));
  if(!t[nrl][ncl]) Error("allocation failure 3 in f3tensor()");
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;

  for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
  for(i=nrl+1;i<=nrh;i++) {
    t[i]=t[i-1]+ncol;
    t[i][ncl]=t[i-1][ncl]+ncol*ndep;
    for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
  }

  return t; // return pointer to array of pointers to rows
}


/*******************************************************************************
Given x1a, x2a, ya, m, n as described in splie2 and y2a
as produced by that routine;
and given a desired interpolating point x1,x2;
this routine returns an interpolated function value y
by bicubic spline interpolation
as well as
the first and second derivatives with respect to x1.
*******************************************************************************/
void splin2x1(DOUBLE x1a[], DOUBLE x2a[], DOUBLE **ya, DOUBLE **y2ax2, int m, int n, DOUBLE x1, DOUBLE x2, DOUBLE *y, DOUBLE *d11y, DOUBLE *d12y) {

  int j;
  DOUBLE *ytmp,*yytmp;

  ytmp=vector(1,m);
  yytmp=vector(1,m);

  //Perform m evaluations of the row splines constructed by splie2,
  //using the one-dimensional spline evaluator splint.

  for(j=1; j<=m; j++ ) splint(x2a,ya[j],y2ax2[j],n,x2,&yytmp[j]);

  spline(x1a,yytmp,m,1.0e30,1.0e30,ytmp);

  //Construct the one-dimensional column spline and evaluate it.
  //splint(x1a,yytmp,ytmp,m,x1,y);

  //Construct the one-dimensional column spline and evaluate it
  //together with the first and second derivatives with respect to x1.
  splintlong(x1a, yytmp, ytmp, m, x1, y, d11y, d12y);

  free_vector(yytmp,1,m);
  free_vector(ytmp,1,m);
}

/*******************************************************************************
Given the arrays xa[1..n] and ya[1..n],
which tabulate a function (with the xaiâ€™s in order),
and given the array y2a[1..n], which is the output from spline above,
and given a value of x,
this routine returns a cubic-spline interpolated value y
as well as
the first and second derivatives.
*******************************************************************************/
void splintlong(DOUBLE xa[], DOUBLE ya[], DOUBLE y2a[], int n, DOUBLE x, DOUBLE *y, DOUBLE *d1y, DOUBLE *d2y) {
  int klo,khi;
  DOUBLE b,a;

  /*We will find the right place in the table by means of bisection.
  This is optimal if sequential calls to this routine are at random values of x.
  If sequential calls are in order, and closely spaced, one would do better
  to store previous values of klo and khi and test if
  they remain appropriate on the next call.*/

  /*klo=1;

  khi=n;

  while(khi-klo > 1) {
   k=(khi+klo) >> 1;

   if( xa[k] > x ) khi=k;
   else klo=k;
  }

  //klo and khi now bracket the input value of x.

  h=xa[khi]-xa[klo];

  if( h == 0.0 ) Error("Bad xa input to routine splint"); //The xaâ€™s must be distinct.*/

  klo = (int) (x*speckles_dx_inv) + 1;
  if(klo == n) klo = n-1;
  khi = klo+1;

  a=(xa[khi]-x)*speckles_dx_inv;
  b=(x-xa[klo])*speckles_dx_inv;

  //Cubic spline polynomial is now evaluated.
  *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*speckles_dx2/6.;

  //First derivative is now evaluated.
  *d1y = (ya[khi]-ya[klo])*speckles_dx_inv + ((1.-3.*a*a)*y2a[klo]+(3.*b*b-1.)*y2a[khi])*speckles_dx/6.0;

  //Second derivative is now evaluated. Note that `a` and `b` are in units of [1/dx] 
  *d2y = a*y2a[klo]+b*y2a[khi];
}


/**************************** fast interpolation ***************************/
// Speckle potential
DOUBLE SpecklesEpot16point(DOUBLE x, DOUBLE y) {
  DOUBLE dx,dy;
  int x1,x2,x3,x4,y1,y2,y3,y4;

  x1 = (int) (x*interp_dx_inv);
  if(x1 == 0) {
    x2 = x1+1;
    dx = x*interp_dx_inv - (double)(x2-1);
    x1 = Npoints;
    x2 = x1+1;
  }
  else {
    x2 = x1+1;
    dx = x*interp_dx_inv - (double)(x2-1);
  }
  x3 = x1+2;
  x4 = x1+3;

  y1 = (int) (y*interp_dx_inv);
  if(y1 == 0) {
    y2 = y1+1;
    dy = y*interp_dx_inv - (double)(y2-1);
    y1 = Npoints;
    y2 = y1+1;
  }
  else {
    y2 = y1+1;
    dy = y*interp_dx_inv - (double)(y2-1);
  }
  y3 = y1+2;
  y4 = y1+3;

#ifdef SECURE
  if(dx<0 || dx>1) Warning("  dx out of range (0 < %lf < %lf), x = %lf\n", dx, speckles_dx, x);
  if((x<x_speckles[x2] || x>x_speckles[x3]) && x2 != N_POINTS) Warning("  dx out of range (2)\n");

  if(dy<0 || dy>1) Warning("  dy out of range (0 < %lf < %lf), y = %lf\n", dy, speckles_dx, y);
  if((y<x_speckles[y2] || y>x_speckles[y3]) && y2 != N_POINTS) Warning("  dy out of range (2)\n");
#endif
  // zeroth accuracy
  // return V[x2][y2];
  // linear accuracy
  //return (1.-dx) *(1.-dy)*V[x2][y2] + (1.-dx)*dy*V[x2][y3] + dx*(1.-dy)*V[x3][y2] + dx*dy*V[x3][y3];
  // cubic accuracy
  return 1./36.*((-2 + dx)*(-1 + dx)*dx*((-2 + dy)*((-1 + dy)*(dy*V[x1][y1] - 3*(1 + dy)*V[x1][y2]) + 3*dy*(1 + dy)*V[x1][y3]) - dy*(-1 + dy*dy)*V[x1][y4]) + 3*(-2 + dx)*(-1 + dx)*(1 + dx)*((-2 + dy)*(-((-1 + dy)*(dy*V[x2][y1] - 3*(1 + dy)*V[x2][y2])) - 3*dy*(1 + dy)*V[x2][y3]) + dy*(-1 + dy*dy)*V[x2][y4]) + 3*(-2 + dx)*dx*(1 + dx)*((-2 + dy)*((-1 + dy)*(dy*V[x3][y1] - 3*(1 + dy)*V[x3][y2]) + 3*dy*(1 + dy)*V[x3][y3]) - dy*(-1 + dy*dy)*V[x3][y4]) + dx*(-1 + dx*dx)*((-2 + dy)*(-((-1 + dy)*(dy*V[x4][y1] - 3*(1 + dy)*V[x4][y2])) - 3*dy*(1 + dy)*V[x4][y3]) + dy*(-1 + dy*dy)*V[x4][y4]));
}

/************************** Speckles U ************************************/
// logarithm of one-body w.f.
// ln f
DOUBLE SpecklesU16point(DOUBLE x, DOUBLE y) {
  DOUBLE dx,dy;
  int x1,x2,x3,x4,y1,y2,y3,y4;

#ifndef INTERPOLATE_LOG_ONE_BODY
  Error("enable INTERPOLATE_LOG_ONE_BODY\n");
#endif

  x1 = (int) (x*interp_dx_inv);
  if(x1 == 0) {
    x2 = x1+1;
    dx = x*interp_dx_inv - (double)(x2-1);
    x1 = Npoints;
    x2 = x1+1;
  }
  else {
    x2 = x1+1;
    dx = x*interp_dx_inv - (double)(x2-1);
  }
  x3 = x1+2;
  x4 = x1+3;

  y1 = (int) (y*interp_dx_inv);
  if(y1 == 0) {
    y2 = y1+1;
    dy = y*interp_dx_inv - (double)(y2-1);
    y1 = Npoints;
    y2 = y1+1;
  }
  else {
    y2 = y1+1;
    dy = y*interp_dx_inv - (double)(y2-1);
  }
  y3 = y1+2;
  y4 = y1+3;

  return 1./36.*((-2 + dx)*(-1 + dx)*dx*((-2 + dy)*((-1 + dy)*(dy*wf[x1][y1] - 3*(1 + dy)*wf[x1][y2]) + 3*dy*(1 + dy)*wf[x1][y3]) - dy*(-1 + dy*dy)*wf[x1][y4]) + 3*(-2 + dx)*(-1 + dx)*(1 + dx)*((-2 + dy)*(-((-1 + dy)*(dy*wf[x2][y1] - 3*(1 + dy)*wf[x2][y2])) - 3*dy*(1 + dy)*wf[x2][y3]) + dy*(-1 + dy*dy)*wf[x2][y4]) + 3*(-2 + dx)*dx*(1 + dx)*((-2 + dy)*((-1 + dy)*(dy*wf[x3][y1] - 3*(1 + dy)*wf[x3][y2]) + 3*dy*(1 + dy)*wf[x3][y3]) - dy*(-1 + dy*dy)*wf[x3][y4]) + dx*(-1 + dx*dx)*((-2 + dy)*(-((-1 + dy)*(dy*wf[x4][y1] - 3*(1 + dy)*wf[x4][y2])) - 3*dy*(1 + dy)*wf[x4][y3]) + dy*(-1 + dy*dy)*wf[x4][y4]));
}

/************************** Speckles Fp ***********************************/
// logarithmic derivative 
// u' = (ux, uy)
void SpecklesFp16point(DOUBLE *Fx, DOUBLE *Fy, DOUBLE x, DOUBLE y) {
  DOUBLE dx,dy;
  int x1,x2,x3,x4,y1,y2,y3,y4;

  x1 = (int) (x*interp_dx_inv);
  if(x1 == 0) {
    x2 = x1+1;
    dx = x*interp_dx_inv - (double)(x2-1);
    x1 = Npoints;
    x2 = x1+1;
  }
  else {
    x2 = x1+1;
    dx = x*interp_dx_inv - (double)(x2-1);
  }
  x3 = x1+2;
  x4 = x1+3;

  y1 = (int) (y*interp_dx_inv);
  if(y1 == 0) {
    y2 = y1+1;
    dy = y*interp_dx_inv - (double)(y2-1);
    y1 = Npoints;
    y2 = y1+1;
  }
  else {
    y2 = y1+1;
    dy = y*interp_dx_inv - (double)(y2-1);
  }
  y3 = y1+2;
  y4 = y1+3;

  *Fx = 1./36.*((2 + 3*(-2 + dx)*dx)*(-2 + dy)*(-1 + dy)*dy*wf[x1][y1] - 3*(2 + 3*(-2 + dx)*dx)*(-2 + dy)*(-1 + dy)*(1 + dy)*wf[x1][y2] + 3*(2 + 3*(-2 + dx)*dx)*(-2 + dy)*dy*(1 + dy)*wf[x1][y3] - (2 + 3*(-2 + dx)*dx)*dy*(-1 + dy*dy)*wf[x1][y4] - 3*(-1 + dx*(-4 + 3*dx))*(-2 + dy)*(-1 + dy)*dy*wf[x2][y1] + 9*(-1 + dx*(-4 + 3*dx))*(-2 + dy)*(-1 + dy)*(1 + dy)*wf[x2][y2] - 9*(-1 + dx*(-4 + 3*dx))*(-2 + dy)*dy*(1 + dy)*wf[x2][y3] + 3*(-1 + dx*(-4 + 3*dx))*dy*(-1 + dy*dy)*wf[x2][y4] + 3*(-2 + dx*(-2 + 3*dx))*(-2 + dy)*(-1 + dy)*dy*wf[x3][y1] - 9*(-2 + dx*(-2 + 3*dx))*(-2 + dy)*(-1 + dy)*(1 + dy)*wf[x3][y2] + 9*(-2 + dx*(-2 + 3*dx))*(-2 + dy)*dy*(1 + dy)*wf[x3][y3] - 3*(-2 + dx*(-2 + 3*dx))*dy*(-1 + dy*dy)*wf[x3][y4] - (-1 + 3*dx*dx)*(-2 + dy)*(-1 + dy)*dy*wf[x4][y1] + 3*(-1 + 3*dx*dx)*(-2 + dy)*(-1 + dy)*(1 + dy)*wf[x4][y2] + 3*(1 - 3*dx*dx)*(-2 + dy)*dy*(1 + dy)*wf[x4][y3] + (-1 + 3*dx*dx)*dy*(-1 + dy*dy)*wf[x4][y4]);
  *Fy = 1./36.*((-2 + dx)*(-1 + dx)*dx*(2 + 3*(-2 + dy)*dy)*wf[x1][y1] - 3*(-2 + dx)*(-1 + dx)*dx*(-1 + dy*(-4 + 3*dy))*wf[x1][y2] + 3*(-2 + dx)*(-1 + dx)*dx*(-2 + dy*(-2 + 3*dy))*wf[x1][y3] - (-2 + dx)*(-1 + dx)*dx*(-1 + 3*dy*dy)*wf[x1][y4] - 3*(-2 + dx)*(-1 + dx)*(1 + dx)*(2 + 3*(-2 + dy)*dy)*wf[x2][y1] + 9*(-2 + dx)*(-1 + dx)*(1 + dx)*(-1 + dy*(-4 + 3*dy))*wf[x2][y2] - 9*(-2 + dx)*(-1 + dx)*(1 + dx)*(-2 + dy*(-2 + 3*dy))*wf[x2][y3] + 3*(-2 + dx)*(-1 + dx)*(1 + dx)*(-1 + 3*dy*dy)*wf[x2][y4] + 3*(-2 + dx)*dx*(1 + dx)*(2 + 3*(-2 + dy)*dy)*wf[x3][y1] - 9*(-2 + dx)*dx*(1 + dx)*(-1 + dy*(-4 + 3*dy))*wf[x3][y2] + 9*(-2 + dx)*dx*(1 + dx)*(-2 + dy*(-2 + 3*dy))*wf[x3][y3] + 3*(-2 + dx)*dx*(1 + dx)*(1 - 3*dy*dy)*wf[x3][y4] - dx*(-1 + dx*dx)*(2 + 3*(-2 + dy)*dy)*wf[x4][y1] + 3*dx*(-1 + dx*dx)*(-1 + dy*(-4 + 3*dy))*wf[x4][y2] - 3*dx*(-1 + dx*dx)*(-2 + dy*(-2 + 3*dy))*wf[x4][y3] + dx*(-1 + dx*dx)*(-1 + 3*dy*dy)*wf[x4][y4]);

  *Fx *= interp_dx_inv;
  *Fy *= interp_dx_inv;
}

/************************** Speckles Fpp **********************************/
// One-body contribution to the kinetic energy
// -f''/f + (f'/f)^2
// -u''+ (u')^2
DOUBLE SpecklesFpp16point(DOUBLE x, DOUBLE y) {
  DOUBLE dx,dy;
  int x1,x2,x3,x4,y1,y2,y3,y4;
  DOUBLE Fxx;

  x1 = (int) (x*interp_dx_inv);
  if(x1 == 0) {
    x2 = x1+1;
    dx = x*interp_dx_inv - (double)(x2-1);
    x1 = Npoints;
    x2 = x1+1;
  }
  else {
    x2 = x1+1;
    dx = x*interp_dx_inv - (double)(x2-1);
  }
  x3 = x1+2;
  x4 = x1+3;

  y1 = (int) (y*interp_dx_inv);
  if(y1 == 0) {
    y2 = y1+1;
    dy = y*interp_dx_inv - (double)(y2-1);
    y1 = Npoints;
    y2 = y1+1;
  }
  else {
    y2 = y1+1;
    dy = y*interp_dx_inv - (double)(y2-1);
  }
  y3 = y1+2;
  y4 = y1+3;

  //Fx = 1./36.*((2 + 3*(-2 + dx)*dx)*(-2 + dy)*(-1 + dy)*dy*wf[x1][y1] - 3*(2 + 3*(-2 + dx)*dx)*(-2 + dy)*(-1 + dy)*(1 + dy)*wf[x1][y2] + 3*(2 + 3*(-2 + dx)*dx)*(-2 + dy)*dy*(1 + dy)*wf[x1][y3] - (2 + 3*(-2 + dx)*dx)*dy*(-1 + dy*dy)*wf[x1][y4] - 3*(-1 + dx*(-4 + 3*dx))*(-2 + dy)*(-1 + dy)*dy*wf[x2][y1] + 9*(-1 + dx*(-4 + 3*dx))*(-2 + dy)*(-1 + dy)*(1 + dy)*wf[x2][y2] - 9*(-1 + dx*(-4 + 3*dx))*(-2 + dy)*dy*(1 + dy)*wf[x2][y3] + 3*(-1 + dx*(-4 + 3*dx))*dy*(-1 + dy*dy)*wf[x2][y4] + 3*(-2 + dx*(-2 + 3*dx))*(-2 + dy)*(-1 + dy)*dy*wf[x3][y1] - 9*(-2 + dx*(-2 + 3*dx))*(-2 + dy)*(-1 + dy)*(1 + dy)*wf[x3][y2] + 9*(-2 + dx*(-2 + 3*dx))*(-2 + dy)*dy*(1 + dy)*wf[x3][y3] - 3*(-2 + dx*(-2 + 3*dx))*dy*(-1 + dy*dy)*wf[x3][y4] - (-1 + 3*dx*dx)*(-2 + dy)*(-1 + dy)*dy*wf[x4][y1] + 3*(-1 + 3*dx*dx)*(-2 + dy)*(-1 + dy)*(1 + dy)*wf[x4][y2] + 3*(1 - 3*dx*dx)*(-2 + dy)*dy*(1 + dy)*wf[x4][y3] + (-1 + 3*dx*dx)*dy*(-1 + dy*dy)*wf[x4][y4]);
  //Fy = 1./36.*((-2 + dx)*(-1 + dx)*dx*(2 + 3*(-2 + dy)*dy)*wf[x1][y1] - 3*(-2 + dx)*(-1 + dx)*dx*(-1 + dy*(-4 + 3*dy))*wf[x1][y2] + 3*(-2 + dx)*(-1 + dx)*dx*(-2 + dy*(-2 + 3*dy))*wf[x1][y3] - (-2 + dx)*(-1 + dx)*dx*(-1 + 3*dy*dy)*wf[x1][y4] - 3*(-2 + dx)*(-1 + dx)*(1 + dx)*(2 + 3*(-2 + dy)*dy)*wf[x2][y1] + 9*(-2 + dx)*(-1 + dx)*(1 + dx)*(-1 + dy*(-4 + 3*dy))*wf[x2][y2] - 9*(-2 + dx)*(-1 + dx)*(1 + dx)*(-2 + dy*(-2 + 3*dy))*wf[x2][y3] + 3*(-2 + dx)*(-1 + dx)*(1 + dx)*(-1 + 3*dy*dy)*wf[x2][y4] + 3*(-2 + dx)*dx*(1 + dx)*(2 + 3*(-2 + dy)*dy)*wf[x3][y1] - 9*(-2 + dx)*dx*(1 + dx)*(-1 + dy*(-4 + 3*dy))*wf[x3][y2] + 9*(-2 + dx)*dx*(1 + dx)*(-2 + dy*(-2 + 3*dy))*wf[x3][y3] + 3*(-2 + dx)*dx*(1 + dx)*(1 - 3*dy*dy)*wf[x3][y4] - dx*(-1 + dx*dx)*(2 + 3*(-2 + dy)*dy)*wf[x4][y1] + 3*dx*(-1 + dx*dx)*(-1 + dy*(-4 + 3*dy))*wf[x4][y2] - 3*dx*(-1 + dx*dx)*(-2 + dy*(-2 + 3*dy))*wf[x4][y3] + dx*(-1 + dx*dx)*(-1 + 3*dy*dy)*wf[x4][y4]);
  Fxx = 1./36.*(6*(-1 + dx)*(-1 + dy)*((-2 + dx)*dx + (-2 + dy)*dy)*wf[x1][y1] - 6*(-1 + dx)*(dx*(4 - 6*dy) + 3*(-2 + dy)*(-1 + dy)*(1 + dy) + dx*dx*(-2 + 3*dy))*wf[x1][y2] + 6*(-1 + dx)*(dx*(2 - 6*dy) + 3*(-2 + dy)*dy*(1 + dy) + dx*dx*(-1 + 3*dy))*wf[x1][y3] - 6*(-1 + dx)*dy*(-1 + (-2 + dx)*dx + dy*dy)*wf[x1][y4] - 6*(-1 + dy)*(6 + 3*dx*(-1 + (-2 + dx)*dx) + 4*dy - 6*dx*dy + (-2 + 3*dx)*dy*dy)*wf[x2][y1] + 18*(-8 + dx*dx*(4 - 6*dy) + dx*dx*dx*(-2 + 3*dy) - 2*dy*(-4 + (-2 + dy)*dy) + dx*(8 + 3*dy*(-2 + (-2 + dy)*dy)))*wf[x2][y2] - 18*(-2 + dx + 2*dx*dx - dx*dx*dx + 10*dy + 3*(-3 + dx)*dx*(1 + dx)*dy - (-2 + 3*dx)*dy*dy + (-2 + 3*dx)*dy*dy*dy)*wf[x2][y3] + 6*dy*(8 - 2*dy*dy + 3*dx*(-2 + (-2 + dx)*dx + dy*dy))*wf[x2][y4] + 6*(-1 + dy)*(3*(-2 + dx)*dx*(1 + dx) + 2*dy - 6*dx*dy + (-1 + 3*dx)*dy*dy)*wf[x3][y1] - 18*(-2 + dx*dx*(2 - 3*dy) + dy - (-2 + dy)*dy*dy + dx*dx*dx*(-2 + 3*dy) + dx*(10 + 3*(-3 + dy)*dy*(1 + dy)))*wf[x3][y2] + 18*(dx*dx*(1 - 3*dy) + dx*dx*dx*(-1 + 3*dy) + dy*(2 + dy - dy*dy) + dx*(2 + 3*dy*(-4 + (-1 + dy)*dy)))*wf[x3][y3] - 6*dy*(1 - dy*dy + 3*dx*(-3 + (-1 + dx)*dx + dy*dy))*wf[x3][y4] - 6*dx*(-1 + dy)*(-1 + dx*dx + (-2 + dy)*dy)*wf[x4][y1] + 6*dx*(8 + dx*dx*(-2 + 3*dy) + 3*dy*(-2 + (-2 + dy)*dy))*wf[x4][y2] - 6*dx*(1 + dx*dx*(-1 + 3*dy) + 3*dy*(-3 + (-1 + dy)*dy))*wf[x4][y3] + 6*dx*dy*(-2 + dx*dx + dy*dy)*wf[x4][y4]);

  //return (-Fxx + Fx*Fx + Fy*Fy)*interp_dx2_inv;
  return -Fxx*interp_dx2_inv;
}

/**************************** fast interpolation ***************************/
// Speckle potential
DOUBLE SpecklesEpot16pointContinuous(DOUBLE x, DOUBLE y) {
  DOUBLE dx,dy;
  int x1,x2,x3,x4,y1,y2,y3,y4;

  x1 = (int) (x*interp_dx_inv);
  if(x1 == 0) {
    x2 = x1+1;
    dx = x*interp_dx_inv - (double)(x2-1);
    x1 = Npoints;
    x2 = x1+1;
  }
  else {
    x2 = x1+1;
    dx = x*interp_dx_inv - (double)(x2-1);
  }
  x3 = x1+2;
  x4 = x1+3;

  y1 = (int) (y*interp_dx_inv);
  if(y1 == 0) {
    y2 = y1+1;
    dy = y*interp_dx_inv - (double)(y2-1);
    y1 = Npoints;
    y2 = y1+1;
  }
  else {
    y2 = y1+1;
    dy = y*interp_dx_inv - (double)(y2-1);
  }
  y3 = y1+2;
  y4 = y1+3;

#ifdef SECURE
  if(dx<0 || dx>speckles_dx) Warning(" dx out of range\n");
  if((x<x_speckles[x2] || x>x_speckles[x3]) && x2 != N_POINTS) Warning("  dx out of range (2)\n");

  if(dy<0 || dy>speckles_dx) Warning("  dy out of range\n");
  if((y<x_speckles[y2] || y>x_speckles[y3]) && y2 != N_POINTS) Warning("  dy out of range (2)\n");
#endif
  // zeroth accuracy
  // return V[x2][y2];
  // linear accuracy
  //return (1.-dx) *(1.-dy)*V[x2][y2] + (1.-dx)*dy*V[x2][y3] + dx*(1.-dy)*V[x3][y2] + dx*dy*V[x3][y3];
  // cubic accuracy
  //return 1./36.*((-2 + dx)*(-1 + dx)*dx*((-2 + dy)*((-1 + dy)*(dy*V[x1][y1] - 3*(1 + dy)*V[x1][y2]) + 3*dy*(1 + dy)*V[x1][y3]) - dy*(-1 + dy*dy)*V[x1][y4]) + 3*(-2 + dx)*(-1 + dx)*(1 + dx)*((-2 + dy)*(-((-1 + dy)*(dy*V[x2][y1] - 3*(1 + dy)*V[x2][y2])) - 3*dy*(1 + dy)*V[x2][y3]) + dy*(-1 + dy*dy)*V[x2][y4]) + 3*(-2 + dx)*dx*(1 + dx)*((-2 + dy)*((-1 + dy)*(dy*V[x3][y1] - 3*(1 + dy)*V[x3][y2]) + 3*dy*(1 + dy)*V[x3][y3]) - dy*(-1 + dy*dy)*V[x3][y4]) + dx*(-1 + dx*dx)*((-2 + dy)*(-((-1 + dy)*(dy*V[x4][y1] - 3*(1 + dy)*V[x4][y2])) - 3*dy*(1 + dy)*V[x4][y3]) + dy*(-1 + dy*dy)*V[x4][y4]));
#ifdef POINT16CONTINUOUS_2deriv
  return ((dx-1)*(dx-1)*(dx-1)*dx*(1 + 2*dx)*(2*V[x1][y2] + dy*((dy-1)*(dy-1)*(dy-1)*(1 + 2*dy)*V[x1][y1] + V[x1][y3] + dy*(-((2 + 3*(-1 + dy)*dy*(-3 + 2*dy))*V[x1][y2]) + V[x1][y3] + (-1 + dy)*dy*(-3 + 2*dy)*(3*V[x1][y3] - V[x1][y4])))) - (-1 + dx)*(2 + 2*dx - 9*dx*dx*dx + 6*dx*dx*dx*dx)*(2*V[x2][y2] + dy*((dy-1)*(dy-1)*(dy-1)*(1 + 2*dy)*V[x2][y1] + V[x2][y3] + dy*(-((2 + 3*(-1 + dy)*dy*(-3 + 2*dy))*V[x2][y2]) + V[x2][y3] + (-1 + dy)*dy*(-3 + 2*dy)*(3*V[x2][y3] - V[x2][y4])))) + dx*(1 + dx*(1 + 3*(-1 + dx)*dx*(-3 + 2*dx)))*(2*V[x3][y2] + dy*((dy-1)*(dy-1)*(dy-1)*(1 + 2*dy)*V[x3][y1] + V[x3][y3] + dy*(-((2 + 3*(-1 + dy)*dy*(-3 + 2*dy))*V[x3][y2]) + V[x3][y3] + (-1 + dy)*dy*(-3 + 2*dy)*(3*V[x3][y3] - V[x3][y4])))) - (-1 + dx)*dx*dx*dx*(-3 + 2*dx)*(2*V[x4][y2] + dy*((dy-1)*(dy-1)*(dy-1)*(1 + 2*dy)*V[x4][y1] + V[x4][y3] + dy*(-((2 + 3*(-1 + dy)*dy*(-3 + 2*dy))*V[x4][y2]) + V[x4][y3] + (-1 + dy)*dy*(-3 + 2*dy)*(3*V[x4][y3] - V[x4][y4])))))/4.;
#else
  return ((dx-3)*(dx-3)*(-2 + dx)*((dy-3)*(dy-3)*(-2 + dy)*V[x1][y1] - (-3 + dy)*(14 + dy*(-14 + 3*dy))*V[x1][y2] + (-2 + dy)*((19 + dy*(-16 + 3*dy))*V[x1][y3] - (-3 + dy)*(-2 + dy)*V[x1][y4])) - (-3 + dx)*(14 + dx*(-14 + 3*dx))*((dy-3)*(dy-3)*(-2 + dy)*V[x2][y1] - (-3 + dy)*(14 + dy*(-14 + 3*dy))*V[x2][y2] + (-2 + dy)*((19 + dy*(-16 + 3*dy))*V[x2][y3] - (-3 + dy)*(-2 + dy)*V[x2][y4])) + (-2 + dx)*(19 + dx*(-16 + 3*dx))*((dy-3)*(dy-3)*(-2 + dy)*V[x3][y1] - (-3 + dy)*(14 + dy*(-14 + 3*dy))*V[x3][y2] + (-2 + dy)*((19 + dy*(-16 + 3*dy))*V[x3][y3] - (-3 + dy)*(-2 + dy)*V[x3][y4])) - (-3 + dx)*(dx-2)*(dx-2)*((dy-3)*(dy-3)*(-2 + dy)*V[x4][y1] - (-3 + dy)*(14 + dy*(-14 + 3*dy))*V[x4][y2] + (-2 + dy)*((19 + dy*(-16 + 3*dy))*V[x4][y3] - (-3 + dy)*(-2 + dy)*V[x4][y4])))/4.;
#endif
}

/************************** Speckles U ************************************/
// logarithm of one-body w.f.
// ln f
DOUBLE SpecklesU16pointContinuous(DOUBLE x, DOUBLE y) {
  DOUBLE dx,dy;
  int x1,x2,x3,x4,y1,y2,y3,y4;

#ifndef INTERPOLATE_LOG_ONE_BODY
  Error("enable INTERPOLATE_LOG_ONE_BODY\n");
#endif

  x1 = (int) (x*interp_dx_inv);
  if(x1 == 0) {
    x2 = x1+1;
    dx = x*interp_dx_inv - (double)(x2-1);
    x1 = Npoints;
    x2 = x1+1;
  }
  else {
    x2 = x1+1;
    dx = x*interp_dx_inv - (double)(x2-1);
  }
  x3 = x1+2;
  x4 = x1+3;

  y1 = (int) (y*interp_dx_inv);
  if(y1 == 0) {
    y2 = y1+1;
    dy = y*interp_dx_inv - (double)(y2-1);
    y1 = Npoints;
    y2 = y1+1;
  }
  else {
    y2 = y1+1;
    dy = y*interp_dx_inv - (double)(y2-1);
  }
  y3 = y1+2;
  y4 = y1+3;

#ifdef POINT16CONTINUOUS_2deriv
  return ((dx-1)*(dx-1)*(dx-1)*dx*(1 + 2*dx)*(2*wf[x1][y2] + dy*((dy-1)*(dy-1)*(dy-1)*(1 + 2*dy)*wf[x1][y1] + wf[x1][y3] + dy*(-((2 + 3*(-1 + dy)*dy*(-3 + 2*dy))*wf[x1][y2]) + wf[x1][y3] + (-1 + dy)*dy*(-3 + 2*dy)*(3*wf[x1][y3] - wf[x1][y4])))) - (-1 + dx)*(2 + 2*dx - 9*dx*dx*dx + 6*dx*dx*dx*dx)*(2*wf[x2][y2] + dy*((dy-1)*(dy-1)*(dy-1)*(1 + 2*dy)*wf[x2][y1] + wf[x2][y3] + dy*(-((2 + 3*(-1 + dy)*dy*(-3 + 2*dy))*wf[x2][y2]) + wf[x2][y3] + (-1 + dy)*dy*(-3 + 2*dy)*(3*wf[x2][y3] - wf[x2][y4])))) + dx*(1 + dx*(1 + 3*(-1 + dx)*dx*(-3 + 2*dx)))*(2*wf[x3][y2] + dy*((dy-1)*(dy-1)*(dy-1)*(1 + 2*dy)*wf[x3][y1] + wf[x3][y3] + dy*(-((2 + 3*(-1 + dy)*dy*(-3 + 2*dy))*wf[x3][y2]) + wf[x3][y3] + (-1 + dy)*dy*(-3 + 2*dy)*(3*wf[x3][y3] - wf[x3][y4])))) - (-1 + dx)*dx*dx*dx*(-3 + 2*dx)*(2*wf[x4][y2] + dy*((dy-1)*(dy-1)*(dy-1)*(1 + 2*dy)*wf[x4][y1] + wf[x4][y3] + dy*(-((2 + 3*(-1 + dy)*dy*(-3 + 2*dy))*wf[x4][y2]) + wf[x4][y3] + (-1 + dy)*dy*(-3 + 2*dy)*(3*wf[x4][y3] - wf[x4][y4])))))/4.;
#else
  return ((dx-1)*(dx-1)*dx*(-2*wf[x1][y2] + dy*((dy-1)*(dy-1)*wf[x1][y1] - wf[x1][y3] + dy*((5 - 3*dy)*wf[x1][y2] + (-4 + 3*dy)*wf[x1][y3] - (-1 + dy)*wf[x1][y4]))) - (2 + dx*dx*(-5 + 3*dx))*(-2*wf[x2][y2] + dy*((dy-1)*(dy-1)*wf[x2][y1] - wf[x2][y3] + dy*((5 - 3*dy)*wf[x2][y2] + (-4 + 3*dy)*wf[x2][y3] - (-1 + dy)*wf[x2][y4]))) + dx*(-1 + dx*(-4 + 3*dx))*(-2*wf[x3][y2] + dy*((dy-1)*(dy-1)*wf[x3][y1] - wf[x3][y3] + dy*((5 - 3*dy)*wf[x3][y2] + (-4 + 3*dy)*wf[x3][y3] - (-1 + dy)*wf[x3][y4]))) - (-1 + dx)*dx*dx*(-2*wf[x4][y2] + dy*((dy-1)*(dy-1)*wf[x4][y1] - wf[x4][y3] + dy*((5 - 3*dy)*wf[x4][y2] + (-4 + 3*dy)*wf[x4][y3] - (-1 + dy)*wf[x4][y4]))))/4.;
#endif
}

/************************** Speckles Fp ***********************************/
// logarithmic derivative 
// u' = (ux, uy)
void SpecklesFp16pointContinuous(DOUBLE *Fx, DOUBLE *Fy, DOUBLE x, DOUBLE y) {
  DOUBLE dx,dy;
  int x1,x2,x3,x4,y1,y2,y3,y4;

  x1 = (int) (x*interp_dx_inv);
  if(x1 == 0) {
    x2 = x1+1;
    dx = x*interp_dx_inv - (double)(x2-1);
    x1 = Npoints;
    x2 = x1+1;
  }
  else {
    x2 = x1+1;
    dx = x*interp_dx_inv - (double)(x2-1);
  }
  x3 = x1+2;
  x4 = x1+3;

  y1 = (int) (y*interp_dx_inv);
  if(y1 == 0) {
    y2 = y1+1;
    dy = y*interp_dx_inv - (double)(y2-1);
    y1 = Npoints;
    y2 = y1+1;
  }
  else {
    y2 = y1+1;
    dy = y*interp_dx_inv - (double)(y2-1);
  }
  y3 = y1+2;
  y4 = y1+3;

#ifdef POINT16CONTINUOUS_2deriv
  *Fx = (2*(dx-1)*(dx-1)*(dx-1)*dx*(2*wf[x1][y2]+dy*((dy-1)*(dy-1)*(dy-1)*(1+2*dy)*wf[x1][y1]+wf[x1][y3]+dy*(-((2+3*(-1+dy)*dy*(-3+2*dy))*wf[x1][y2])+wf[x1][y3]+(-1+dy)*dy*(-3+2*dy)*(3*wf[x1][y3]-wf[x1][y4]))))+(dx-1)*(dx-1)*(dx-1)*(1+2*dx)*(2*wf[x1][y2]+dy*((dy-1)*(dy-1)*(dy-1)*(1+2*dy)*wf[x1][y1]+wf[x1][y3]+dy*(-((2+3*(-1+dy)*dy*(-3+2*dy))*wf[x1][y2])+wf[x1][y3]+(-1+dy)*dy*(-3+2*dy)*(3*wf[x1][y3]-wf[x1][y4]))))+3*(dx-1)*(dx-1)*dx*(1+2*dx)*(2*wf[x1][y2]+dy*((dy-1)*(dy-1)*(dy-1)*(1+2*dy)*wf[x1][y1]+wf[x1][y3]+dy*(-((2+3*(-1+dy)*dy*(-3+2*dy))*wf[x1][y2])+wf[x1][y3]+(-1+dy)*dy*(-3+2*dy)*(3*wf[x1][y3]-wf[x1][y4]))))-(2+2*dx-9*dx*dx*dx+6*dx*dx*dx*dx)*(2*wf[x2][y2]+dy*((dy-1)*(dy-1)*(dy-1)*(1+2*dy)*wf[x2][y1]+wf[x2][y3]+dy*(-((2+3*(-1+dy)*dy*(-3+2*dy))*wf[x2][y2])+wf[x2][y3]+(-1+dy)*dy*(-3+2*dy)*(3*wf[x2][y3]-wf[x2][y4]))))-(-1+dx)*(2+3*dx*dx*(-9+8*dx))*(2*wf[x2][y2]+dy*((dy-1)*(dy-1)*(dy-1)*(1+2*dy)*wf[x2][y1]+wf[x2][y3]+dy*(-((2+3*(-1+dy)*dy*(-3+2*dy))*wf[x2][y2])+wf[x2][y3]
+(-1+dy)*dy*(-3+2*dy)*(3*wf[x2][y3]-wf[x2][y4]))))+(1+dx*(1+3*(-1+dx)*dx*(-3+2*dx)))*(2*wf[x3][y2]+dy*((dy-1)*(dy-1)*(dy-1)*(1+2*dy)*wf[x3][y1]+wf[x3][y3]+dy*(-((2+3*(-1+dy)*dy*(-3+2*dy))*wf[x3][y2])+wf[x3][y3]+(-1+dy)*dy*(-3+2*dy)*(3*wf[x3][y3]-wf[x3][y4]))))+dx*(1+3*dx*(6+dx*(-15+8*dx)))*(2*wf[x3][y2]+dy*((dy-1)*(dy-1)*(dy-1)*(1+2*dy)*wf[x3][y1]+wf[x3][y3]+dy*(-((2+3*(-1+dy)*dy*(-3+2*dy))*wf[x3][y2])+wf[x3][y3]+(-1+dy)*dy*(-3+2*dy)*(3*wf[x3][y3]-wf[x3][y4]))))-2*(-1+dx)*dx*dx*dx*(2*wf[x4][y2]+dy*((dy-1)*(dy-1)*(dy-1)*(1+2*dy)*wf[x4][y1]+wf[x4][y3]+dy*(-((2+3*(-1+dy)*dy*(-3+2*dy))*wf[x4][y2])+wf[x4][y3]+(-1+dy)*dy*(-3+2*dy)*(3*wf[x4][y3]-wf[x4][y4]))))-3*(-1+dx)*dx*dx*(-3+2*dx)*(2*wf[x4][y2]+dy*((dy-1)*(dy-1)*(dy-1)*(1+2*dy)*wf[x4][y1]+wf[x4][y3]+dy*(-((2+3*(-1+dy)*dy*(-3+2*dy))*wf[x4][y2])+wf[x4][y3]+(-1+dy)*dy*(-3+2*dy)*(3*wf[x4][y3]-wf[x4][y4]))))-dx*dx*dx*(-3+2*dx)*(2*wf[x4][y2]+dy*((dy-1)*(dy-1)*(dy-1)*(1+2*dy)*wf[x4][y1]+wf[x4][y3]+dy*(-((2+3*(-1+dy)*dy*(-3+2*dy))*wf[x4][y2])+wf[x4][y3]+(-1+dy)*dy*(-3+2*dy)*(3*wf[x4][y3]-wf[x4][y4])))))/4.;
  *Fy = ((dx-1)*(dx-1)*(dx-1)*dx*(1+2*dx)*((dy-1)*(dy-1)*(-1+10*dy*dy)*wf[x1][y1]+wf[x1][y3]+dy*(-((4+3*dy*(9+10*(-2+dy)*dy))*wf[x1][y2])+2*wf[x1][y3]+dy*(9+10*(-2+dy)*dy)*(3*wf[x1][y3]-wf[x1][y4])))-(-1+dx)*(2+2*dx-9*dx*dx*dx+6*dx*dx*dx*dx)*((dy-1)*(dy-1)*(-1+10*dy*dy)*wf[x2][y1]+wf[x2][y3]+dy*(-((4+3*dy*(9+10*(-2+dy)*dy))*wf[x2][y2])+2*wf[x2][y3]+dy*(9+10*(-2+dy)*dy)*(3*wf[x2][y3]-wf[x2][y4])))+dx*(1+dx*(1+3*(-1+dx)*dx*(-3+2*dx)))*((dy-1)*(dy-1)*(-1+10*dy*dy)*wf[x3][y1]+wf[x3][y3]+dy*(-((4+3*dy*(9+10*(-2+dy)*dy))*wf[x3][y2])+2*wf[x3][y3]+dy*(9+10*(-2+dy)*dy)*(3*wf[x3][y3]-wf[x3][y4])))-(-1+dx)*dx*dx*dx*(-3+2*dx)*((dy-1)*(dy-1)*(-1+10*dy*dy)*wf[x4][y1]+wf[x4][y3]+dy*(-((4+3*dy*(9+10*(-2+dy)*dy))*wf[x4][y2])+2*wf[x4][y3]+dy*(9+10*(-2+dy)*dy)*(3*wf[x4][y3]-wf[x4][y4]))))/4.;
#else
  *Fx = ((dx-1)*(dx-1)*(-2*wf[x1][y2] + dy*((dy-1)*(dy-1)*wf[x1][y1] - wf[x1][y3] + dy*((5 - 3*dy)*wf[x1][y2] + (-4 + 3*dy)*wf[x1][y3] - (-1 + dy)*wf[x1][y4]))) + 2*(-1 + dx)*dx*(-2*wf[x1][y2] + dy*((dy-1)*(dy-1)*wf[x1][y1] - wf[x1][y3] + dy*((5 - 3*dy)*wf[x1][y2] + (-4 + 3*dy)*wf[x1][y3] - (-1 + dy)*wf[x1][y4]))) - dx*(-10 + 9*dx)*(-2*wf[x2][y2] + dy*((dy-1)*(dy-1)*wf[x2][y1] - wf[x2][y3] + dy*((5 - 3*dy)*wf[x2][y2] + (-4 + 3*dy)*wf[x2][y3] - (-1 + dy)*wf[x2][y4]))) + (-1 + dx)*(1 + 9*dx)*(-2*wf[x3][y2] + dy*((dy-1)*(dy-1)*wf[x3][y1] - wf[x3][y3] + dy*((5 - 3*dy)*wf[x3][y2] + (-4 + 3*dy)*wf[x3][y3] - (-1 + dy)*wf[x3][y4]))) - 2*(-1 + dx)*dx*(-2*wf[x4][y2] + dy*((dy-1)*(dy-1)*wf[x4][y1] - wf[x4][y3] + dy*((5 - 3*dy)*wf[x4][y2] + (-4 + 3*dy)*wf[x4][y3] - (-1 + dy)*wf[x4][y4]))) - dx*dx*(-2*wf[x4][y2] + dy*((dy-1)*(dy-1)*wf[x4][y1] - wf[x4][y3] + dy*((5 - 3*dy)*wf[x4][y2] + (-4 + 3*dy)*wf[x4][y3] - (-1 + dy)*wf[x4][y4]))))/4.;
  *Fy = ((dx-1)*(dx-1)*dx*((-1 + dy)*(-1 + 3*dy)*wf[x1][y1] - wf[x1][y3] + dy*((10 - 9*dy)*wf[x1][y2] + (-8 + 9*dy)*wf[x1][y3] + (2 - 3*dy)*wf[x1][y4])) - (2 + dx*dx*(-5 + 3*dx))*((-1 + dy)*(-1 + 3*dy)*wf[x2][y1] - wf[x2][y3] + dy*((10 - 9*dy)*wf[x2][y2] + (-8 + 9*dy)*wf[x2][y3] + (2 - 3*dy)*wf[x2][y4])) + dx*(-1 + dx*(-4 + 3*dx))*((-1 + dy)*(-1 + 3*dy)*wf[x3][y1] - wf[x3][y3] + dy*((10 - 9*dy)*wf[x3][y2] + (-8 + 9*dy)*wf[x3][y3] + (2 - 3*dy)*wf[x3][y4])) - (-1 + dx)*dx*dx*((-1 + dy)*(-1 + 3*dy)*wf[x4][y1] - wf[x4][y3] + dy*((10 - 9*dy)*wf[x4][y2] + (-8 + 9*dy)*wf[x4][y3] + (2 - 3*dy)*wf[x4][y4])))/4.;
#endif
  *Fx *= interp_dx_inv;
  *Fy *= interp_dx_inv;
}

/************************** Speckles Fpp **********************************/
// One-body contribution to the kinetic energy
// -f''/f + (f'/f)^2
// -u''+ (u')^2
DOUBLE SpecklesFpp16pointContinuous(DOUBLE x, DOUBLE y) {
  DOUBLE dx,dy;
  int x1,x2,x3,x4,y1,y2,y3,y4;
  DOUBLE Fxx_yy;

  x1 = (int) (x*interp_dx_inv);
  if(x1 == 0) {
    x2 = x1+1;
    dx = x*interp_dx_inv - (double)(x2-1);
    x1 = Npoints;
    x2 = x1+1;
  }
  else {
    x2 = x1+1;
    dx = x*interp_dx_inv - (double)(x2-1);
  }
  x3 = x1+2;
  x4 = x1+3;

  y1 = (int) (y*interp_dx_inv);
  if(y1 == 0) {
    y2 = y1+1;
    dy = y*interp_dx_inv - (double)(y2-1);
    y1 = Npoints;
    y2 = y1+1;
  }
  else {
    y2 = y1+1;
    dy = y*interp_dx_inv - (double)(y2-1);
  }
  y3 = y1+2;
  y4 = y1+3;

#ifdef POINT16CONTINUOUS_2deriv
  //Fx = (2*(dx-1)*(dx-1)*(dx-1)*dx*(2*wf[x1][y2]+dy*((dy-1)*(dy-1)*(dy-1)*(1+2*dy)*wf[x1][y1]+wf[x1][y3]+dy*(-((2+3*(-1+dy)*dy*(-3+2*dy))*wf[x1][y2])+wf[x1][y3]+(-1+dy)*dy*(-3+2*dy)*(3*wf[x1][y3]-wf[x1][y4]))))+(dx-1)*(dx-1)*(dx-1)*(1+2*dx)*(2*wf[x1][y2]+dy*((dy-1)*(dy-1)*(dy-1)*(1+2*dy)*wf[x1][y1]+wf[x1][y3]+dy*(-((2+3*(-1+dy)*dy*(-3+2*dy))*wf[x1][y2])+wf[x1][y3]+(-1+dy)*dy*(-3+2*dy)*(3*wf[x1][y3]-wf[x1][y4]))))+3*(dx-1)*(dx-1)*dx*(1+2*dx)*(2*wf[x1][y2]+dy*((dy-1)*(dy-1)*(dy-1)*(1+2*dy)*wf[x1][y1]+wf[x1][y3]+dy*(-((2+3*(-1+dy)*dy*(-3+2*dy))*wf[x1][y2])+wf[x1][y3]+(-1+dy)*dy*(-3+2*dy)*(3*wf[x1][y3]-wf[x1][y4]))))-(2+2*dx-9*dx*dx*dx+6*dx*dx*dx*dx)*(2*wf[x2][y2]+dy*((dy-1)*(dy-1)*(dy-1)*(1+2*dy)*wf[x2][y1]+wf[x2][y3]+dy*(-((2+3*(-1+dy)*dy*(-3+2*dy))*wf[x2][y2])+wf[x2][y3]+(-1+dy)*dy*(-3+2*dy)*(3*wf[x2][y3]-wf[x2][y4]))))-(-1+dx)*(2+3*dx*dx*(-9+8*dx))*(2*wf[x2][y2]+dy*((dy-1)*(dy-1)*(dy-1)*(1+2*dy)*wf[x2][y1]+wf[x2][y3]+dy*(-((2+3*(-1+dy)*dy*(-3+2*dy))*wf[x2][y2])+wf[x2][y3]
  //+(-1+dy)*dy*(-3+2*dy)*(3*wf[x2][y3]-wf[x2][y4]))))+(1+dx*(1+3*(-1+dx)*dx*(-3+2*dx)))*(2*wf[x3][y2]+dy*((dy-1)*(dy-1)*(dy-1)*(1+2*dy)*wf[x3][y1]+wf[x3][y3]+dy*(-((2+3*(-1+dy)*dy*(-3+2*dy))*wf[x3][y2])+wf[x3][y3]+(-1+dy)*dy*(-3+2*dy)*(3*wf[x3][y3]-wf[x3][y4]))))+dx*(1+3*dx*(6+dx*(-15+8*dx)))*(2*wf[x3][y2]+dy*((dy-1)*(dy-1)*(dy-1)*(1+2*dy)*wf[x3][y1]+wf[x3][y3]+dy*(-((2+3*(-1+dy)*dy*(-3+2*dy))*wf[x3][y2])+wf[x3][y3]+(-1+dy)*dy*(-3+2*dy)*(3*wf[x3][y3]-wf[x3][y4]))))-2*(-1+dx)*dx*dx*dx*(2*wf[x4][y2]+dy*((dy-1)*(dy-1)*(dy-1)*(1+2*dy)*wf[x4][y1]+wf[x4][y3]+dy*(-((2+3*(-1+dy)*dy*(-3+2*dy))*wf[x4][y2])+wf[x4][y3]+(-1+dy)*dy*(-3+2*dy)*(3*wf[x4][y3]-wf[x4][y4]))))-3*(-1+dx)*dx*dx*(-3+2*dx)*(2*wf[x4][y2]+dy*((dy-1)*(dy-1)*(dy-1)*(1+2*dy)*wf[x4][y1]+wf[x4][y3]+dy*(-((2+3*(-1+dy)*dy*(-3+2*dy))*wf[x4][y2])+wf[x4][y3]+(-1+dy)*dy*(-3+2*dy)*(3*wf[x4][y3]-wf[x4][y4]))))-dx*dx*dx*(-3+2*dx)*(2*wf[x4][y2]+dy*((dy-1)*(dy-1)*(dy-1)*(1+2*dy)*wf[x4][y1]+wf[x4][y3]+dy*(-((2+3*(-1+dy)*dy*(-3+2*dy))*wf[x4][y2])+wf[x4][y3]+(-1+dy)*dy*(-3+2*dy)*(3*wf[x4][y3]-wf[x4][y4])))))/4.;
  //Fy = ((dx-1)*(dx-1)*(dx-1)*dx*(1+2*dx)*((dy-1)*(dy-1)*(-1+10*dy*dy)*wf[x1][y1]+wf[x1][y3]+dy*(-((4+3*dy*(9+10*(-2+dy)*dy))*wf[x1][y2])+2*wf[x1][y3]+dy*(9+10*(-2+dy)*dy)*(3*wf[x1][y3]-wf[x1][y4])))-(-1+dx)*(2+2*dx-9*dx*dx*dx+6*dx*dx*dx*dx)*((dy-1)*(dy-1)*(-1+10*dy*dy)*wf[x2][y1]+wf[x2][y3]+dy*(-((4+3*dy*(9+10*(-2+dy)*dy))*wf[x2][y2])+2*wf[x2][y3]+dy*(9+10*(-2+dy)*dy)*(3*wf[x2][y3]-wf[x2][y4])))+dx*(1+dx*(1+3*(-1+dx)*dx*(-3+2*dx)))*((dy-1)*(dy-1)*(-1+10*dy*dy)*wf[x3][y1]+wf[x3][y3]+dy*(-((4+3*dy*(9+10*(-2+dy)*dy))*wf[x3][y2])+2*wf[x3][y3]+dy*(9+10*(-2+dy)*dy)*(3*wf[x3][y3]-wf[x3][y4])))-(-1+dx)*dx*dx*dx*(-3+2*dx)*((dy-1)*(dy-1)*(-1+10*dy*dy)*wf[x4][y1]+wf[x4][y3]+dy*(-((4+3*dy*(9+10*(-2+dy)*dy))*wf[x4][y2])+2*wf[x4][y3]+dy*(9+10*(-2+dy)*dy)*(3*wf[x4][y3]-wf[x4][y4]))))/4.;
  Fxx_yy = (2*(dx-1)*(dx-1)*(dx-1)*(2*wf[x1][y2]+dy*((dy-1)*(dy-1)*(dy-1)*(1+2*dy)*wf[x1][y1]+wf[x1][y3]+dy*(-((2+3*(-1+dy)*dy*(-3+2*dy))*wf[x1][y2])+wf[x1][y3]+(-1+dy)*dy*(-3+2*dy)*(3*wf[x1][y3]-wf[x1][y4]))))+6*(dx-1)*(dx-1)*dx*(2*wf[x1][y2]+dy*((dy-1)*(dy-1)*(dy-1)*(1+2*dy)*wf[x1][y1]+wf[x1][y3]+dy*(-((2+3*(-1+dy)*dy*(-3+2*dy))*wf[x1][y2])+wf[x1][y3]+(-1+dy)*dy*(-3+2*dy)*(3*wf[x1][y3]-wf[x1][y4]))))+3*(dx-1)*(dx-1)*(1+2*dx)*(2*wf[x1][y2]+dy*((dy-1)*(dy-1)*(dy-1)*(1+2*dy)*wf[x1][y1]+wf[x1][y3]+dy*(-((2+3*(-1+dy)*dy*(-3+2*dy))*wf[x1][y2])+wf[x1][y3]+(-1+dy)*dy*(-3+2*dy)*(3*wf[x1][y3]-wf[x1][y4]))))
+3*(-1+dx)*dx*(1+2*dx)*(2*wf[x1][y2]+dy*((dy-1)*(dy-1)*(dy-1)*(1+2*dy)*wf[x1][y1]+wf[x1][y3]+dy*(-((2+3*(-1+dy)*dy*(-3+2*dy))*wf[x1][y2])+wf[x1][y3]+(-1+dy)*dy*(-3+2*dy)*
(3*wf[x1][y3]-wf[x1][y4]))))+(dx-1)*(dx-1)*(dx-1)*dx*(1+2*dx)*((1+dy*(9+10*dy*(-3+2*dy)))*wf[x1][y1]-2*wf[x1][y2]+wf[x1][y3]-dy*(9+10*dy*(-3+2*dy))*(3*wf[x1][y2]-3*wf[x1]
[y3]+wf[x1][y4]))-9*(-1+dx)*dx*(-3+4*dx)*(2*wf[x2][y2]+dy*((dy-1)*(dy-1)*(dy-1)*(1+2*dy)*wf[x2][y1]+wf[x2][y3]+dy*(-((2+3*(-1+dy)*dy*(-3+2*dy))*wf[x2][y2])+wf[x2][y3]
+(-1+dy)*dy*(-3+2*dy)*(3*wf[x2][y3]-wf[x2][y4]))))-(2+3*dx*dx*(-9+8*dx))*(2*wf[x2][y2]+dy*((dy-1)*(dy-1)*(dy-1)*(1+2*dy)*wf[x2][y1]+wf[x2][y3]+dy*(-((2+3*(-1+dy)*dy*(-3+2*dy))*wf[x2][y2])+wf[x2][y3]+(-1+dy)*dy*(-3+2*dy)*(3*wf[x2][y3]-wf[x2][y4]))))-(-1+dx)*(2+2*dx-9*dx*dx*dx+6*dx*dx*dx*dx)*((1+dy*(9+10*dy*(-3+2*dy)))*wf[x2][y1]-2*wf[x2][y2]+wf[x2][y3]-dy*(9+10*dy*(-3+2*dy))*(3*wf[x2][y2]-3*wf[x2][y3]+wf[x2][y4]))+9*(-1+dx)*dx*(-1+4*dx)*(2*wf[x3][y2]+dy*((dy-1)*(dy-1)*(dy-1)*(1+2*dy)*wf[x3][y1]+wf[x3][y3]+dy*(-((2+3*(-1+dy)*dy*(-3+2*dy))*wf[x3][y2])+wf[x3][y3]+(-1+dy)*dy*(-3+2*dy)*(3*wf[x3][y3]-wf[x3][y4]))))+(1+3*dx*(6+dx*(-15+8*dx)))*(2*wf[x3][y2]+dy*((dy-1)*(dy-1)*(dy-1)*(1+2*dy)*wf[x3][y1]+wf[x3][y3]+dy*(-((2+3*(-1+dy)*dy*(-3+2*dy))*wf[x3][y2])+wf[x3][y3]+(-1+dy)*dy*(-3+2*dy)*(3*wf[x3][y3]-wf[x3][y4]))))+dx*(1+dx*(1+3*(-1+dx)*dx*(-3+2*dx)))*((1+dy*(9+10*dy*(-3+2*dy)))*wf[x3][y1]-2*wf[x3][y2]+wf[x3][y3]-dy*(9+10*dy*(-3+2*dy))*(3*wf[x3][y2]-3*wf[x3][y3]+wf[x3][y4]))-6*(-1+dx)*dx*dx*(2*wf[x4][y2]+dy*((dy-1)*(dy-1)*(dy-1)*(1+2*dy)*wf[x4][y1]+wf[x4][y3]+dy*(-((2+3*(-1+dy)*dy*(-3+2*dy))*wf[x4][y2])+wf[x4][y3]+(-1+dy)*dy*(-3+2*dy)*(3*wf[x4][y3]-wf[x4][y4]))))-2*dx*dx*dx*(2*wf[x4][y2]+dy*((dy-1)*(dy-1)*(dy-1)*(1+2*dy)*wf[x4][y1]+wf[x4][y3]+dy*(-((2+3*(-1+dy)*dy*(-3+2*dy))*wf[x4][y2])+wf[x4][y3]+(-1+dy)*dy*(-3+2*dy)*(3*wf[x4][y3]-wf[x4][y4]))))-3*(-1+dx)*dx*(-3+2*dx)*(2*wf[x4][y2]+dy*((dy-1)*(dy-1)*(dy-1)*(1+2*dy)*wf[x4][y1]+wf[x4][y3]+dy*(-((2+3*(-1+dy)*dy*(-3+2*dy))*wf[x4][y2])+wf[x4][y3]+(-1+dy)*dy*(-3+2*dy)*(3*wf[x4][y3]-wf[x4][y4]))))-3*dx*dx*(-3+2*dx)*(2*wf[x4][y2]+dy*((dy-1)*(dy-1)*(dy-1)*(1+2*dy)*wf[x4][y1]+wf[x4][y3]+dy*(-((2+3*(-1+dy)*dy*(-3+2*dy))*wf[x4][y2])+wf[x4][y3]+(-1+dy)*dy*(-3+2*dy)*(3*wf[x4][y3]-wf[x4][y4]))))-(-1+dx)*dx*dx*dx*(-3+2*dx)*((1+dy*(9+10*dy*(-3+2*dy)))*wf[x4][y1]-2*wf[x4][y2]+wf[x4][y3]-dy*(9+10*dy*(-3+2*dy))*(3*wf[x4][y2]-3*wf[x4][y3]+wf[x4][y4])))/2.;
#else
  //Fx = ((dx-1)*(dx-1)*(-2*wf[x1][y2] + dy*((dy-1)*(dy-1)*wf[x1][y1] - wf[x1][y3] + dy*((5 - 3*dy)*wf[x1][y2] + (-4 + 3*dy)*wf[x1][y3] - (-1 + dy)*wf[x1][y4]))) + 2*(-1 + dx)*dx*(-2*wf[x1][y2] + dy*((dy-1)*(dy-1)*wf[x1][y1] - wf[x1][y3] + dy*((5 - 3*dy)*wf[x1][y2] + (-4 + 3*dy)*wf[x1][y3] - (-1 + dy)*wf[x1][y4]))) - dx*(-10 + 9*dx)*(-2*wf[x2][y2] + dy*((dy-1)*(dy-1)*wf[x2][y1] - wf[x2][y3] + dy*((5 - 3*dy)*wf[x2][y2] + (-4 + 3*dy)*wf[x2][y3] - (-1 + dy)*wf[x2][y4]))) + (-1 + dx)*(1 + 9*dx)*(-2*wf[x3][y2] + dy*((dy-1)*(dy-1)*wf[x3][y1] - wf[x3][y3] + dy*((5 - 3*dy)*wf[x3][y2] + (-4 + 3*dy)*wf[x3][y3] - (-1 + dy)*wf[x3][y4]))) - 2*(-1 + dx)*dx*(-2*wf[x4][y2] + dy*((dy-1)*(dy-1)*wf[x4][y1] - wf[x4][y3] + dy*((5 - 3*dy)*wf[x4][y2] + (-4 + 3*dy)*wf[x4][y3] - (-1 + dy)*wf[x4][y4]))) - dx*dx*(-2*wf[x4][y2] + dy*((dy-1)*(dy-1)*wf[x4][y1] - wf[x4][y3] + dy*((5 - 3*dy)*wf[x4][y2] + (-4 + 3*dy)*wf[x4][y3] - (-1 + dy)*wf[x4][y4]))))/4.;
  //Fy = ((dx-1)*(dx-1)*dx*((-1 + dy)*(-1 + 3*dy)*wf[x1][y1] - wf[x1][y3] + dy*((10 - 9*dy)*wf[x1][y2] + (-8 + 9*dy)*wf[x1][y3] + (2 - 3*dy)*wf[x1][y4])) - (2 + dx*dx*(-5 + 3*dx))*((-1 + dy)*(-1 + 3*dy)*wf[x2][y1] - wf[x2][y3] + dy*((10 - 9*dy)*wf[x2][y2] + (-8 + 9*dy)*wf[x2][y3] + (2 - 3*dy)*wf[x2][y4])) + dx*(-1 + dx*(-4 + 3*dx))*((-1 + dy)*(-1 + 3*dy)*wf[x3][y1] - wf[x3][y3] + dy*((10 - 9*dy)*wf[x3][y2] + (-8 + 9*dy)*wf[x3][y3] + (2 - 3*dy)*wf[x3][y4])) - (-1 + dx)*dx*dx*((-1 + dy)*(-1 + 3*dy)*wf[x4][y1] - wf[x4][y3] + dy*((10 - 9*dy)*wf[x4][y2] + (-8 + 9*dy)*wf[x4][y3] + (2 - 3*dy)*wf[x4][y4])))/4.;
  Fxx_yy = ((dx-1)*(dx-1)*dx*((-2 + 3*dy)*wf[x1][y1] + (5 - 9*dy)*wf[x1][y2] + (-4 + 9*dy)*wf[x1][y3] + wf[x1][y4] - 3*dy*wf[x1][y4]) + 2*(-1 + dx)*(-2*wf[x1][y2] + dy*((dy-1)*(dy-1)*wf[x1][y1] - wf[x1][y3] + dy*((5 - 3*dy)*wf[x1][y2] + (-4 + 3*dy)*wf[x1][y3] - (-1 + dy)*wf[x1][y4]))) + dx*(-2*wf[x1][y2] + dy*((dy-1)*(dy-1)*wf[x1][y1] - wf[x1][y3] + dy*((5 - 3*dy)*wf[x1][y2] + (-4 + 3*dy)*wf[x1][y3] - (-1 + dy)*wf[x1][y4]))) - (2 + dx*dx*(-5 + 3*dx))*((-2 + 3*dy)*wf[x2][y1] + (5 - 9*dy)*wf[x2][y2] + (-4 + 9*dy)*wf[x2][y3] + wf[x2][y4] - 3*dy*wf[x2][y4]) - (-5 + 9*dx)*(-2*wf[x2][y2] + dy*((dy-1)*(dy-1)*wf[x2][y1] - wf[x2][y3] + dy*((5 - 3*dy)*wf[x2][y2] + (-4 + 3*dy)*wf[x2][y3] - (-1 + dy)*wf[x2][y4]))) + dx*(-1 + dx*(-4 + 3*dx))*((-2 + 3*dy)*wf[x3][y1] + (5 - 9*dy)*wf[x3][y2] + (-4 + 9*dy)*wf[x3][y3] + wf[x3][y4] - 3*dy*wf[x3][y4]) + (-4 + 9*dx)*(-2*wf[x3][y2] + dy*((dy-1)*(dy-1)*wf[x3][y1] - wf[x3][y3] + dy*((5 - 3*dy)*wf[x3][y2] + (-4 + 3*dy)*wf[x3][y3] - (-1 + dy)*wf[x3][y4]))) - (-1 + dx)*dx*dx*((-2 + 3*dy)*wf[x4][y1] + (5 - 9*dy)*wf[x4][y2] + (-4 + 9*dy)*wf[x4][y3] + wf[x4][y4] - 3*dy*wf[x4][y4]) - (-1 + dx)*(-2*wf[x4][y2] + dy*((dy-1)*(dy-1)*wf[x4][y1] - wf[x4][y3] + dy*((5 - 3*dy)*wf[x4][y2] + (-4 + 3*dy)*wf[x4][y3] - (-1 + dy)*wf[x4][y4]))) + 2*dx*(2*wf[x4][y2] + dy*(-((dy-1)*(dy-1)*wf[x4][y1]) + wf[x4][y3] + dy*((-5 + 3*dy)*wf[x4][y2] + (4 - 3*dy)*wf[x4][y3] + (-1 + dy)*wf[x4][y4]))))/2.;
#endif
  //return (-Fxx_yy + Fx*Fx + Fy*Fy)*interp_dx2_inv;
  return -Fxx_yy*interp_dx2_inv;
}

/**************************** fast interpolation ***************************/
// Speckle potential
#ifdef POINT4

void SpecklesInitialize4(void) {
  int i,j;
  DOUBLE f00,f01,f02,f03,f10,f11,f12,f13,f20,f21,f22,f23,f30,f31,f32,f33;
  DOUBLE f11xx,f21xx,f12xx,f22xx,f11yy,f12yy,f21yy,f22yy;

#ifndef INTERPOLATE_LOG_ONE_BODY
  Error("enable INTERPOLATE_LOG_ONE_BODY\n");
#endif

  // virtual points i = 0 and i == N are already defined
  for(i=1; i<Npoints; i++) {
    for(j=1; j<Npoints; j++) {
      // interpolate potential
      f00 = V[i-1][j-1];
      f01 = V[i-1][j];
      f02 = V[i-1][j+1];
      f03 = V[i-1][j+2];
      f10 = V[i][j-1];
      f11 = V[i][j];
      f12 = V[i][j+1];
      f13 = V[i][j+2];
      f20 = V[i+1][j-1];
      f21 = V[i+1][j];
      f22 = V[i+1][j+1];
      f23 = V[i+1][j+2];
      f30 = V[i+1][j-1];
      f31 = V[i+2][j];
      f32 = V[i+2][j+1];
      f33 = V[i+2][j+2];

      f11xx = f21 - 2.*f11 + f01;
      f21xx = f31 - 2.*f21 + f11;
      f12xx = f22 - 2.*f12 + f02;
      f22xx = f32 - 2.*f22 + f12;
      f11yy = f12 - 2.*f11 + f10;
      f12yy = f13 - 2.*f12 + f11;
      f21yy = f22 - 2.*f21 + f20;
      f22yy = f23 - 2.*f22 + f21;

      v[i-1][j-1] = f11;
      vx[i-1][j-1] = f21 - f11 - (2.*f11xx + f21xx)/6.;
      vy[i-1][j-1] = f12 - f11 - (2.*f11yy + f12yy)/6.;
      vxx[i-1][j-1] = f11xx/2.;
      vyy[i-1][j-1] = f11yy/2.;
      vxy[i-1][j-1] = f11 - f12 - f21 + f22 + (f11xx + f11yy - f12xx - f21yy)/3. + (f12yy + f21xx - f22xx - f22yy)/6.;
      vxxx[i-1][j-1] = (f21xx - f11xx)/6.;
      vyyy[i-1][j-1] = (f12yy - f11yy)/6.;
      vxxy[i-1][j-1] =  (f12xx - f11xx)/2.;
      vxyy[i-1][j-1] = (f21yy - f11yy)/2.;
      vxxxy[i-1][j-1] = (f11xx + f22xx - f12xx - f21xx)/6.;
      vxyyy[i-1][j-1] = (f11yy + f22yy - f12yy - f21yy)/6.;

      // interpolate w.f.
      f00 = wf[i-1][j-1];
      f01 = wf[i-1][j];
      f02 = wf[i-1][j+1];
      f03 = wf[i-1][j+2];
      f10 = wf[i][j-1];
      f11 = wf[i][j];
      f12 = wf[i][j+1];
      f13 = wf[i][j+2];
      f20 = wf[i+1][j-1];
      f21 = wf[i+1][j];
      f22 = wf[i+1][j+1];
      f23 = wf[i+1][j+2];
      f30 = wf[i+1][j-1];
      f31 = wf[i+2][j];
      f32 = wf[i+2][j+1];
      f33 = wf[i+2][j+2];

      f11xx = f21 - 2.*f11 + f01;
      f21xx = f31 - 2.*f21 + f11;
      f12xx = f22 - 2.*f12 + f02;
      f22xx = f32 - 2.*f22 + f12;
      f11yy = f12 - 2.*f11 + f10;
      f12yy = f13 - 2.*f12 + f11;
      f21yy = f22 - 2.*f21 + f20;
      f22yy = f23 - 2.*f22 + f21;

      w[i-1][j-1] = f11;
      wx[i-1][j-1] = f21 - f11 - (2.*f11xx + f21xx)/6.;
      wy[i-1][j-1] = f12 - f11 - (2.*f11yy + f12yy)/6.;
      wxx[i-1][j-1] = f11xx/2.;
      wyy[i-1][j-1] = f11yy/2.;
      wxy[i-1][j-1] = f11 - f12 - f21 + f22 + (f11xx + f11yy - f12xx - f21yy)/3. + (f12yy + f21xx - f22xx - f22yy)/6.;
      wxxx[i-1][j-1] = (f21xx - f11xx)/6.;
      wyyy[i-1][j-1] = (f12yy - f11yy)/6.;
      wxxy[i-1][j-1] =  (f12xx - f11xx)/2.;
      wxyy[i-1][j-1] = (f21yy - f11yy)/2.;
      wxxxy[i-1][j-1] = (f11xx + f22xx - f12xx - f21xx)/6.;
      wxyyy[i-1][j-1] = (f11yy + f22yy - f12yy - f21yy)/6.;
    }
  }
}
#endif

DOUBLE SpecklesEpot4point(DOUBLE x, DOUBLE y) {
  DOUBLE dx,dy;
  int i,j;

  // i goes from 1 to N
  i = (int) (x*interp_dx_inv);
  j = (int) (y*interp_dx_inv);
  // dx = (x - x[0]) / h
  dx = x*interp_dx_inv - (double)i;
  dy = y*interp_dx_inv - (double)j;

  return v[i][j]  + vx[i][j]*dx + vxx[i][j]*dx*dx + vxxx[i][j]*dx*dx*dx 
    + vy[i][j]*dy + vxy[i][j]*dx*dy + vxxy[i][j]*dx*dx*dy 
    + vxxxy[i][j]*dx*dx*dx*dy + vyy[i][j]*dy*dy + vxyy[i][j]*dx*dy*dy 
    + vyyy[i][j]*dy*dy*dy + vxyyy[i][j]*dx*dy*dy*dy;
}

/************************** Speckles U ************************************/
// logarithm of one-body w.f.
// ln f
DOUBLE SpecklesU4point(DOUBLE x, DOUBLE y) {
  DOUBLE dx,dy;
  int i,j;

  // i goes from 1 to N
  i = (int) (x*interp_dx_inv);
  j = (int) (y*interp_dx_inv);
  // dx = (x - x[0]) / h
  dx = x*interp_dx_inv - (double)i;
  dy = y*interp_dx_inv - (double)j;

  return w[i][j]  + wx[i][j]*dx + wxx[i][j]*dx*dx + wxxx[i][j]*dx*dx*dx 
    + wy[i][j]*dy + wxy[i][j]*dx*dy + wxxy[i][j]*dx*dx*dy 
    + wxxxy[i][j]*dx*dx*dx*dy + wyy[i][j]*dy*dy + wxyy[i][j]*dx*dy*dy 
    + wyyy[i][j]*dy*dy*dy + wxyyy[i][j]*dx*dy*dy*dy;
}

/************************** Speckles Fp ***********************************/
// logarithmic derivative 
// u' = (ux, uy)
void SpecklesFp4point(DOUBLE *Fx, DOUBLE *Fy, DOUBLE x, DOUBLE y) {
  DOUBLE dx,dy;
  int i,j;

  // i goes from 1 to N
  i = (int) (x*interp_dx_inv);
  j = (int) (y*interp_dx_inv);
  // dx = (x - x[0]) / h
  dx = x*interp_dx_inv - (double)i;
  dy = y*interp_dx_inv - (double)j;

  *Fx = wx[i][j] + 2.*wxx[i][j]*dx + 3.*wxxx[i][j]*dx*dx + wxy[i][j]*dy
+ 2.*wxxy[i][j]*dx*dy + 3.*wxxxy[i][j]*dx*dx*dy + wxyy[i][j]*dy*dy + wxyyy[i][j]*dy*dy*dy;

  *Fy = wy[i][j] + 2.*wyy[i][j]*dy + 3.*wyyy[i][j]*dy*dy + wxy[i][j]*dx
+ 2.*wxyy[i][j]*dx*dy + 3.*wxyyy[i][j]*dx*dy*dy + wxxy[i][j]*dx*dx + wxxxy[i][j]*dx*dx*dx;

  *Fx *= interp_dx_inv;
  *Fy *= interp_dx_inv;
}

/************************** Speckles Fpp **********************************/
// One-body contribution to the kinetic energy
// -f''/f + (f'/f)^2
// -u''+ (u')^2
DOUBLE SpecklesFpp4point(DOUBLE x, DOUBLE y) {
  DOUBLE dx,dy;
  int i,j;

  // i goes from 1 to N
  i = (int) (x*interp_dx_inv);
  j = (int) (y*interp_dx_inv);
  // dx = (x - x[0]) / h
  dx = x*interp_dx_inv - (double)i;
  dy = y*interp_dx_inv - (double)j;

  return -interp_dx2_inv*( 2.*(wxx[i][j] + wyy[i][j] + wxyy[i][j]*dx + wxxy[i][j]*dy) + 6.*(wxxx[i][j]*dx + wyyy[i][j]*dy + wxxxy[i][j]*dx*dy + wxyyy[i][j]*dx*dy) );
  //return -Fxx*interp_dx2_inv;
}
