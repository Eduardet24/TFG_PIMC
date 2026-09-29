/*spline.c*/
#include <stdlib.h>
#include <stdio.h>
#include "spline.h"
#include "main.h"
#include "utils.h"
#include "trial.h"
#include "memory.h"
#include "compatab.h"
#include "libigor.h"
#include MATHINCLUDE

#define SPLINE_EPSILON 1e-10

/****************************** SaveWaveFunctionSpline *****************************/
void SaveWaveFunctionSpline(struct Grid *G, char *file_out) {
  int i;
  int Nspline_export = 1000;
  FILE *out;
  DOUBLE x, yInter, ypInter, yppInter;
  static int first_time = ON;

  out = fopen(file_out, "w");

  if(first_time) Message("Saving wavefunction (spline).\n");

  fprintf(out, "x f fp fpp\n");

  if(Nspline_export>10*G->size) Nspline_export = G->size;

  for(i=0; i<Nspline_export; i++) {
    x = G->min + 0.1*G->step * (DOUBLE) i;
    InterpolateSpline(G, x, &yInter, &ypInter, &yppInter);
    //yInter = Exp(InterpolateSplineU(G, x));
    //ypInter = InterpolateSplineFp(G, x);
    //yppInter = InterpolateSplineE(G, x);

    //if(yInter<0) Message("  Bad value of Spline function (negative), distance %" LF " = %" LF " x dx, value %" LF "\n", x, 0.1*(DOUBLE) i, yInter);
    fprintf(out, "%.15" LF " %.15" LF " %.15" LF " %.15" LF "\n", x, yInter, ypInter, yppInter);
  }

  fclose(out);

  first_time = OFF;
}

/************************************ spline construction *******************************/
void SplineConstruct(struct Grid *G, DOUBLE yp1, DOUBLE ypn) {
  int i, k;
  DOUBLE p, qn, sig, un;
  int n;

  n = G->size;

  //G->fpp[0] = G->fp[0] = 0.; // natural spline
  if(yp1>1e30) // zero second derivative
    G->fpp[0] = G->fp[0]=0.;
  else { // first derivative equal to yp1
    G->fpp[0] = -0.5;
    G->fp[0] = (3./(G->x[1]-G->x[0]))*((G->f[1]-G->f[0])/(G->x[1]-G->x[0])-yp1);
  }

  for(i=1; i<n-1; i++) {
    sig=(G->x[i]-G->x[i-1])/(G->x[i+1]-G->x[i-1]);
    p=sig*G->fpp[i-1]+2.0;
    G->fpp[i]=(sig-1.0)/p;
    G->fp[i]=(G->f[i+1]-G->f[i])/(G->x[i+1]-G->x[i]) - (G->f[i]-G->f[i-1])/(G->x[i]-G->x[i-1]);
    G->fp[i]=(6.0*G->fp[i]/(G->x[i+1]-G->x[i-1])-sig*G->fp[i-1])/p;
  }

  if(ypn > 1e30)
    qn = un = 0.;
  else {
    qn=0.5;
    un=(3.0/(G->x[n-1]-G->x[n-2]))*(ypn-(G->f[n-1]-G->f[n-2])/(G->x[n-1]-G->x[n-2]));
  }
  G->fpp[n-1]=(un-qn*G->fp[n-2])/(qn*G->fpp[n-2]+1.);
  for(k=n-2; k>=0; k--) G->fpp[k]=G->fpp[k]*G->fpp[k+1]+G->fp[k];

  // u = (DOUBLE*) Calloc("u", G->size, sizeof(DOUBLE));
  // use lnf array as temporary
  /*n = G->size;

  if(yp1>1e30) // zero second derivative
    G->fpp[1] = G->lnf[1]=0.;
  else { // first derivative equal to yp1
    G->fpp[1] = -0.5;
    G->lnf[1] = (3./(G->x[2]-G->x[1]))*((G->f[2]-G->f[1])/(G->x[2]-G->x[1])-yp1);
  }

  for(i=2; i<=G->size-1; i++) {
    sig = (G->x[i]-G->x[i-1])/(G->x[i+1]-G->x[i-1]);
    p = sig*G->fpp[i-1]+2.;
    G->fpp[i] = (sig-1.)/p;
    G->lnf[i] = (G->f[i+1]-G->f[i])/(G->x[i+1]-G->x[i]) - (G->f[i]-G->f[i-1])/(G->x[i]-G->x[i-1]);
    G->lnf[i] = (6.*G->lnf[i]/(G->x[i+1]-G->x[i-1])-sig*G->lnf[i-1])/p;
  }

  if(ypn>1e30) // zero second derivative
    qn = un = 0.;
  else { // first derivative equal to ypn
    qn = 0.5;
    un = (3./(G->x[n]-G->x[n-1]))*(ypn-(G->f[n]-G->f[n-1])/(G->x[n]-G->x[n-1]));
  }

  G->fpp[n] = (un-qn*G->lnf[n-1])/(qn*G->fpp[n-1]+1.);

  // second derivative
  for(k=n-1;k>=1;k--) G->fpp[k] = G->fpp[k]*G->fpp[k+1] + G->lnf[k];*/
}

/************************************ interpolate spline  *******************************/
void InterpolateSpline(struct Grid *G, DOUBLE r, DOUBLE *y, DOUBLE *yp, DOUBLE *ypp) {
  DOUBLE a,b;
  int i;

  /*int klo,khi,k; // uncomment for non uniform grid
  DOUBLE h;
  
  klo=0;
  khi=G->size-1;
  while(khi-klo > 1) {
    //k=(khi+klo) >> 1;
    k = (int)(0.5*((DOUBLE)(khi+klo)));
    if(G->x[k]>r)
      khi=k;
    else 
      klo=k;
  }
  h=G->x[khi]-G->x[klo];
  i = klo;
  if(h == 0.) Warning("Bad input to InterpolateSpline");*/

  // uniform grid
  if(r<=G->min) r = G->min; // adjust lower bound
  //i = (int) ((r-G->min)*G->I_step + 0.5);
  //i = (int) ((r-G->min)*G->I_step);
  i = (int) ((r-G->min)*G->I_step+SPLINE_EPSILON); // shift by a small amount for the check at points G->x[i]

//#ifdef SECURE
  if(i<G->size-1 && (r+SPLINE_EPSILON<G->x[i] || r-SPLINE_EPSILON>G->x[i+1])) 
    Warning("Interp. spline check failed i=%i (%" LE " <%" LE " <%" LE " )\n", i, G->x[i], r, G->x[i+1]);
//#endif
  if(i >= G->size-1) i = G->size-2; // adjust upper bound

  a = (G->x[i+1]-r)*G->I_step;
  b = (r-G->x[i])*G->I_step;

  // (3.3.3)
  *y   = a*G->f[i] + b*G->f[i+1] + ((a*a*a-a)*G->fpp[i] + (b*b*b-b)*G->fpp[i+1])*(G->step*G->step)/6.;
  *yp  = (G->f[i+1]-G->f[i])*G->I_step + ((1.-3.*a*a)*G->fpp[i]+(3.*b*b-1.)*G->fpp[i+1])*G->step/6.;
  *ypp = a*G->fpp[i] + b*G->fpp[i+1];
}

/************************ Interpolate U ************************************/
// U = ln(f)
DOUBLE InterpolateSplineU(struct Grid *G, DOUBLE r) {
  DOUBLE a,b,f;
  int i;

  if(r<=G->min) r = G->min; // adjust lower bound
  i = (int) ((r-G->min)*G->I_step+SPLINE_EPSILON); // shift by a small amount for the check at points G->x[i]
#ifdef SECURE
  if(i<G->size-1 && (r+SPLINE_EPSILON<G->x[i] || r-SPLINE_EPSILON>G->x[i+1])) Warning("  Interpolate spline U check failed, i=%i (%" LE " < %" LE "  <%" LE " )\n", i, G->x[i], r, G->x[i+1]);
#endif
  //if(i >= G->size-1) i = G->size-2; // adjust upper bound
  if(i >= G->size-1) 
#ifdef INTERPOLATE_LOG // u(r)
    //return G->f[G->size-1];
    return 0.;
#else
    //return Log(G->f[G->size-1]);
    return 0.;
#endif

  a = (G->x[i+1]-r)*G->I_step;
  b = (r-G->x[i])*G->I_step;

  f = a*G->f[i] + b*G->f[i+1] + ((a*a*a-a)*G->fpp[i] + (b*b*b-b)*G->fpp[i+1])*(G->step*G->step)/6.;

#ifdef INTERPOLATE_LOG // u(r)
  return f;
#else // f(r)
#ifdef SECURE
  if(f<0) Warning("Interpolate spline, ln() of negative argument\n");
#endif
  if(f<0) return -1e8;

  return Log(f);
#endif
}

/************************ Interpolate F ************************************/
// same as U, but allows negative values (i.e. interpolation of Vint(r))
DOUBLE InterpolateSplineF(struct Grid *G, DOUBLE r) {
  DOUBLE a,b,f;
  int i;

  if(r<=G->min) r = G->min; // adjust lower bound
  i = (int) ((r-G->min)*G->I_step+SPLINE_EPSILON); // shift by a small amount for the check at points G->x[i]
  //if(i >= G->size-1) i = G->size-2; // adjust upper bound
  if(i >= G->size-1) return 0.;

  a = (G->x[i+1]-r)*G->I_step;
  b = (r-G->x[i])*G->I_step;

  f = a*G->f[i] + b*G->f[i+1] + ((a*a*a-a)*G->fpp[i] + (b*b*b-b)*G->fpp[i+1])*(G->step*G->step)/6.;

  return f;
}

/************************ InterpolateSplineFp ******************************/
// fp = f'/f
DOUBLE InterpolateSplineFp(struct Grid *G, DOUBLE r) {
  DOUBLE a,b;
  int i;
#ifndef INTERPOLATE_LOG
  DOUBLE fp, f;
#endif

  if(r<=G->min) r = G->min; // adjust lower bound
  i = (int) ((r-G->min)*G->I_step+SPLINE_EPSILON); // shift by a small amount for the check at points G->x[i]
#ifdef SECURE
  if(i<G->size-1 && (r+SPLINE_EPSILON<G->x[i] || r-SPLINE_EPSILON>G->x[i+1])) Warning("  Interpolate spline Fp check failed, i=%i (%" LE " < %" LE "  <%" LE " )\n", i, G->x[i], r, G->x[i+1]);
#endif
  //if(i >= G->size-1) i = G->size-2; // adjust upper bound
  if(i >= G->size-1) return 0.;

  a = (G->x[i+1]-r)*G->I_step;
  b = (r-G->x[i])*G->I_step;

#ifdef INTERPOLATE_LOG // f'/f = u'
  return (G->f[i+1]-G->f[i])*G->I_step + ((1.-3.*a*a)*G->fpp[i]+(3.*b*b-1.)*G->fpp[i+1])*G->step/6.;
#else // f'/f
  f = a*G->f[i] + b*G->f[i+1] + ((a*a*a-a)*G->fpp[i] + (b*b*b-b)*G->fpp[i+1])*(G->step*G->step)/6.;
  fp  = (G->f[i+1]-G->f[i])*G->I_step + ((1.-3.*a*a)*G->fpp[i]+(3.*b*b-1.)*G->fpp[i+1])*G->step/6.;
#ifdef SECURE
  if(f<0) Warning("Interpolate spline Fp of negative argument\n");
#endif
  if(f<0) return 0;

  return fp/f;
#endif
}

/************************ InterpolateSpline E*******************************/
// Eloc = -(f" +2/r f')/f+(f'/f)^2
DOUBLE InterpolateSplineE(struct Grid *G, DOUBLE r) {
  DOUBLE a,b;
  int i;
#ifdef INTERPOLATE_LOG
  DOUBLE up,upp;
#else
  DOUBLE f,fp,fpp;
#endif

  if(r == 0) return 0.;

  if(r<=G->min) r = G->min; // adjust lower bound
  i = (int) ((r-G->min)*G->I_step+SPLINE_EPSILON); // shift by a small amount for the check at points G->x[i]
#ifdef SECURE
  if(i<G->size-1 && (r+SPLINE_EPSILON<G->x[i] || r-SPLINE_EPSILON>G->x[i+1])) Warning("  Interpolate spline E check failed, i=%i (%" LE " < %" LE "  <%" LE " )\n", i, G->x[i], r, G->x[i+1]);
#endif
  //if(i >= G->size-1) i = G->size-2; // adjust upper bound
  if(i >= G->size-1) return 0.;

  a = (G->x[i+1]-r)*G->I_step;
  b = (r-G->x[i])*G->I_step;

#ifdef INTERPOLATE_LOG // (-f" + 2/r f')/f + (f'/f)^2 = - u" - u' 2/r
  up  = (G->f[i+1]-G->f[i])*G->I_step + ((1.-3.*a*a)*G->fpp[i]+(3.*b*b-1.)*G->fpp[i+1])*G->step/6.;
  upp = a*G->fpp[i] + b*G->fpp[i+1];

#ifdef TRIAL_1D
  return - upp;
#endif
#ifdef TRIAL_2D
  return - upp - up/r;
#endif
#ifdef TRIAL_3D
  return - upp - 2.*up/r;
#endif
#else // (-f" + 2/r f')/f + (f'/f)^2
  f = a*G->f[i] + b*G->f[i+1] + ((a*a*a-a)*G->fpp[i] + (b*b*b-b)*G->fpp[i+1])*(G->step*G->step)/6.;
  fp  = (G->f[i+1]-G->f[i])*G->I_step + ((1.-3.*a*a)*G->fpp[i]+(3.*b*b-1.)*G->fpp[i+1])*G->step/6.;
  fpp = a*G->fpp[i] + b*G->fpp[i+1];

  fp /= f;
  fpp  /= f;

  if(r==0.) return 0.;
#ifdef SECURE
  if(f<0) Warning("Interpolate spline E of negative argument\n");
#endif
  if(f<0) return 0;

#ifdef TRIAL_1D
  return -fpp + fp*fp;
#endif
#ifdef TRIAL_2D
  return -fpp - fp/r + fp*fp;
#endif
#ifdef TRIAL_3D
  return -fpp - 2.*fp/r + fp*fp;
#endif

#endif
}

/************************** Linear Interpolation Construct ******************************/
// initialize Fp and E
void LinearInterpolationConstruct(struct Grid *G) {
  int i, k;
  DOUBLE p, qn, sig, un;
  DOUBLE Fp, Fpp;
  DOUBLE dr,dr2;
  int n;

  n = G->size;

  /*if(yp1>1e30) // zero second derivative
    G->fpp[0] = G->fp[0]=0.;
  else { // first derivative equal to yp1
    G->fpp[0] = -0.5;
    G->fp[0] = (3./(G->x[1]-G->x[0]))*((G->f[1]-G->f[0])/(G->x[1]-G->x[0])-yp1);
  }*/

  /*Message("  damping two-body solution\n");
  Message("x f(x) Escat\n");
  dr = G->step;
  dr2 = dr*dr;

  for(i=1; i<n-1; i++) {
    Fp = (G->f[i+1]-G->f[i-1])*0.5/dr;
    Fpp = (G->f[i+1]-2.*G->f[i]+G->f[i-1])/dr2;

    Fp /= G->f[i];
    Fpp /= G->f[i];

    //Message("%lf %e %e\n", G->x[i], G->f[i], -(Fpp +(DIMENSION-1.)*Fp/G->x[i]) + InteractionEnergy(G->x[i]));
  }*/

  for(i=1; i<n-1; i++) {
    Fp = (G->f[i+1]-G->f[i-1])/(G->x[i+1]-G->x[i-1]);
    Fpp = (G->f[i+1]-2.*G->f[i]+G->f[i-1])/(G->x[i+1]-G->x[i])/(G->x[i]-G->x[i-1]);

    Fp /= G->f[i];
    Fpp /= G->f[i];

    G->lnf[i] = log(G->f[i]);
    G->fp[i] = Fp;
//3D Eloc = [-(f"+2f'/r)/f] + (f'/f)^2
//2D Eloc = [-(f" + f'/r)/f] + (f'/f)^2
//1D Eloc = [-f"/f] + (f'/f)^2
    G->E[i] = -(Fpp +(DIMENSION-1.)*Fp/G->x[i]) + Fp*Fp;
  }

  G->lnf[0] = log(G->f[0]);
  G->fp[0] = G->fp[1];
  G->E[0] = G->E[1];

  G->lnf[n-1] = log(G->f[n-1]);
  G->fp[n-1] = G->fp[n-2];
  G->E[n-1] = G->E[n-2];

  //Warning("  Pair interaction potential is added in the Jastrow spline\n");
  //for(i=0; i<n; i++) {
  //  G->E[i] += InteractionEnergy(G->x[i]);
  //}
  //Warning("  Setting: gaussian_alpha = gaussian_beta = 0\n");
  //gaussian_alpha = gaussian_beta = 0;
}

/************************** Kurbakov spline *********************************************/
/*#define TRIAL_KURBAKOV_NUMERIC// uses Rpar, Epar, Dpar, Mpar, a, aa, Kurbakov_BC
#define INTERACTION_KURBAKOV_GENERAL// uses D, Apar, Ipar
DOUBLE Rpar,Epar,Dpar,Mpar,Apar,Ipar,D,aa,a,n,Lwf,Lhalfwf;//  external
int Kurbakov_BC,N;//                                parameters
DOUBLE InteractionEnergy(DOUBLE r){
#ifdef INTERACTION_KURBAKOV
  return((2./r-2./sqrt(r*r+Apar*Apar))/(Apar*Apar));
#endif
#ifdef INTERACTION_KURBAKOV_GENERAL
#endif
}*/
/************************** Construct Grid Kurbakov Numeric **********************/
// // // // // // // // // // // // // // // // // // // // // // // // // //
//                                 structures
//
struct param{// structure for auxiliary parameters and arrays
  int memX,memmin,memmax,memadd;// memory volumes
  DOUBLE Doubl,doubl;// cutoffs for DOUBLE-type quantities
  int BC;// boundary condition for u(t) at small r: first (BC=1) or second (BC=2) kind only
  DOUBLE ac;// hard-core radius using in numeric scattering problem
  DOUBLE aa;// s-wave scattering length of pseudopotential
  DOUBLE deriv;// (dot u(a))/u(a), where a=asinh(log(Rpar)) and u(t)=asinh(psi(exp(sinh(t))))
  DOUBLE dpsiR;// (dot psi(Rpar))/psi(Rpar)
  int jm,jj;// lengths of arrays for u(t)
  DOUBLE xi;// parameter xi of trial WF
            // Kurbakov A-FESHBACH
  DOUBLE n0,g,Nxi;// Bogoliubov's condensate fraction, dressed coupling constant, and n/(2*g*n)^(dim/2)
  DOUBLE *yj;// auxiliary array for interpolation of pair interaction
};
struct wfnum{// structure for trial WF
  DOUBLE rmin;// parameters of grid: x(r)=((r-rmin)/(rmax-rmin))^(1/4),
  DOUBLE rmax;//                     r_j=rmin+((j-j0+.5)*Dx)^4*(rmax-rmin),
  DOUBLE Dx;//                       x_j=(j-j0+.5)/(J-2*j0),
  int J;//                           J=1/Dx+2*j0, 0<=j<J<=memJ
  int m,j0,nmax;// m is order of scheme, j0=2*m+1, nmax=2*j0/3
  int memJ;// maximal memory for lnf[]
  DOUBLE mulEpar,DmulEpar;// Epar --> Epar*mulEpar; mulEpar/=DmulEpar, DmulEpar*=DmulEpar
  DOUBLE acut;// minimal value of r (hard-core radius using in Monte Carlo)
  DOUBLE rMAX;// maximal value of r
  DOUBLE *lnf;// array for log(psi(x(r)))
};
struct Ugrid{// structure for pair interaction
  DOUBLE rmin,rmax;// parameters of grid: x(r)=((r-rmin)/(rmax-rmin))^(1/K),
                   //                                        r_j=rmin+(x_j)^K*(rmax-rmin),
                   //                                        x_j=(j-j0)*Dx, Dx=1./jmax,
                   //                                        0<=j<J<=memJ, K=32/(1+DIMENSION/2),
  DOUBLE a1,a2;//                         m=4, j0=2*m+1, J=jmax+2*j0+1,
  int m,j0,jmax,J,K;//                    a1=1/(rmax-rmin), a2=1./K
  int memJ;// maximal memory for array Urj[]
  DOUBLE U0;// U0=4*pi*n is characteristic energy for pair interaction
  DOUBLE err;// error of interpolation
  DOUBLE *Urj;// array for U(x(r))
};
struct Runge scat;
struct param tnk;
struct wfnum twf;
struct Ugrid ugr;
struct spline Ur;
// // // // // // // // // // // // // // // // // // // // // // // // // //
//            block for primary calculating dot u(asinh(log(Rpar)))
//
// calculating of square root for tunnelling
DOUBLE tunnsqrt(DOUBLE t){
  DOUBLE x,y;
  x=exp(sinh(t));
  y=(Dpar-1)*(Dpar-3)/4+x*x*(InteractionEnergy(x)*Mpar-Epar);
  if(y<0)y=0;else y=sqrt(y)*cosh(t);
  return(y);
}
// calculating of tunnelling exponent: primary correcting ac
void tunnexp(void){
  int j,jj;
  DOUBLE a,b,dt,t,y,s;
  if((DIMENSION==3||DIMENSION==1)&&tnk.ac<tnk.doubl*tnk.doubl)tnk.ac=tnk.doubl*tnk.doubl;
  if((DIMENSION==3||DIMENSION==1)&&tnk.BC==3&&tnk.ac<tnk.doubl)tnk.ac=tnk.doubl;
  if(DIMENSION==2&&tnk.ac<1/tnk.Doubl)tnk.ac=1/tnk.Doubl;
  jj=tnk.memmax;a=asinh(log(Rpar));b=asinh(log(tnk.ac));dt=(b-a)/jj;
  for(s=j=0;j<=jj;j++){
    t=a+j*dt;y=tunnsqrt(t);s+=-dt*y;if(s>log(tnk.Doubl)){
      if(j){j=jj;tnk.ac=exp(sinh(t));}else{
        Error("\nFATAL ERROR: ac and Rpar are quite coincided");
      } } } }
// dot u(t)
DOUBLE dut(int _l,DOUBLE *X,DOUBLE *ul){
  DOUBLE f,t;
  t=X[_l];f=ul[2];
  return(f);
}
// dot dot u(t)
DOUBLE ddut(int _l,DOUBLE *X,DOUBLE *ul){
  DOUBLE f,t;
  t=X[_l];
  f=-ul[2]*ul[2]*tanh(ul[1])+(tanh(t)+(2-Dpar)*cosh(t))*ul[2]
    +(InteractionEnergy(exp(sinh(t)))*Mpar-Epar)*exp(2*sinh(t))*cosh(t)*cosh(t)*tanh(ul[1]);
  return(f);
}
// array of functions
void scatut(int _l,DOUBLE *X,DOUBLE *ul,DOUBLE *fk){
  fk[1]=dut(_l,X,ul);fk[2]=ddut(_l,X,ul);
}
// runge integration at r<=Rpar
DOUBLE utd(DOUBLE *X){
  DOUBLE x;
  x=X[0];
  scat.a=asinh(log(Rpar));scat.b=asinh(log(tnk.ac));
  scat.s=1;scat.nn=(int)(fabs(scat.b-scat.a)*tnk.memmin);
  if(scat.nn<tnk.memmin-4)scat.nn=tnk.memmin-4;
  if(scat.nn>tnk.memmax-4-tnk.memadd)scat.nn=tnk.memmax-4-tnk.memadd;
  scat.yk0[1]=asinh((DOUBLE)1.);scat.yk0[2]=x;scat._l=1;
  _runge(X,&scat,scatut);
  if(tnk.BC==1||tnk.BC==2)return(scat.ykxj[tnk.BC]);
  if(tnk.BC==3){// pseudopotential
    if(tnk.aa<=tnk.ac){
      Error("\nFATAL ERROR: amplitude of pseudopotential is too small");
    }if(DIMENSION==3)return((tnk.ac/tnk.aa+1)*sinh(scat.ykxj[1])+
      scat.ykxj[2]*cosh(scat.ykxj[1])/cosh(scat.b));// 3D
    if(DIMENSION==2)return(sinh(scat.ykxj[1])*cosh(scat.b)+(log(2*tnk.aa)-gammaE-
      sinh(scat.b))*scat.ykxj[2]*cosh(scat.ykxj[1]));// 2D
    if(DIMENSION==1)return(tnk.ac/tnk.aa*sinh(scat.ykxj[1])*cosh(scat.b)+
      scat.ykxj[2]*cosh(scat.ykxj[1]));// 1D
  }Error("\nutd: FATAL ERROR in code");return(0);
}
// primary calculating dot u(asinh(log(Rpar)))
void dua(DOUBLE *X){
  int j,_l;
  DOUBLE a1,a2,j1,j2;
  _l=0;if(tnk.BC==3)tnk.deriv=_zero(-1000,1000,_l,X,utd);else{
    X[0]=j1=-10;a1=utd(X);X[_l]=j2=2;a2=utd(X);
    for(j2=j=2;j<32;j++)if(a1*a2>0){
      a1=a2;j1=j2;j2*=j;X[_l]=j2;a2=utd(X);
    }tnk.deriv=_zero(j1,j2,_l,X,utd);
  } }
// // // // // // // // // // // // // // // // // // // // // // // // // //
//        block for calculating ac, psi'(Rpar) and numeric data array
//
// correcting ac
void corrac(void){
  int j;
  DOUBLE t;
  if(0)// from r=Rpar to r=0
    for(j=0;j<scat.nn;j++){
      t=scat.a+j*(scat.b-scat.a)/scat.nn;
      if(fabs(scat.arr[j])<fabs(scat.arr[j+1])){
        tnk.ac=exp(sinh(t));j=scat.nn;
      } }else// from r=0 to r=Rpar
      for(j=scat.nn;j>0;j--){
        t=scat.a+(j-1)*(scat.b-scat.a)/scat.nn;
        if(fabs(scat.arr[j])<fabs(scat.arr[j-1])){
          tnk.ac=exp(sinh(t));j=0;
        } } }
// tetsing positivity of trial WF
int testpositivity(void){
  int j;
  for(j=0;j<scat.nn;j++)if(scat.arr[j]<=0){
    Warning("\nWarning: trial wave function changes a sign\n  scaling: Epar/=%g",twf.DmulEpar);return(1);
  }return(0);
}
// runge integration with 3 points at r<ac
void utdwrac(DOUBLE *X){
  scat.a=asinh(log(Rpar));scat.b=asinh(log(tnk.ac));
  scat.b+=3*(scat.b-scat.a)/tnk.jj;scat.nn=tnk.jj+3;// with 3 points at r<ac
  scat.s=1;scat.yk0[1]=asinh((DOUBLE)1.);scat.yk0[2]=tnk.deriv;scat._l=1;
  _runge(X,&scat,scatut);
}
// d psi(Rpar) / d r (one mind psi(Rpar)=1)
void dpsiRpardr(void){
  tnk.dpsiR=tnk.deriv*cosh(asinh((DOUBLE)1.))/Rpar/cosh(asinh(log(Rpar)));
}
// calculating jm
void nnadd(void){
  int j,jm;
  for(jm=j=0;j<tnk.memadd;j++)if(sinh(scat.arr[j])>sinh(-0.5))jm++;// !! bad
  tnk.jm=jm;
}
// runge integration at r>=Rpar
void utdadd(DOUBLE *X){
  int j;
  DOUBLE a1;
  for(j=tnk.jj+3;j>=0;j--)scat.arr[j+tnk.jm]=scat.arr[j];// shifting data by jm
  scat.a=asinh(log(Rpar));scat.b=scat.a+tnk.jm*(scat.a-asinh(log(tnk.ac)))/tnk.jj;
  scat.yk0[1]=asinh((DOUBLE)1.);scat.yk0[2]=tnk.deriv;scat.s=1;scat.nn=tnk.jm;scat._l=1;
  _runge(X,&scat,scatut);
  for(j=0;j<tnk.jm/2;j++){// reflecting data
    a1=scat.arr[j];scat.arr[j]=scat.arr[tnk.jm-j];scat.arr[tnk.jm-j]=a1;
  }scat.nn=tnk.jj;scat.b=asinh(log(tnk.ac));
}
// final correcting ac and deriv; primary calculating dot psi(Rpar); calculating j0
int acderivdpsiRj0(DOUBLE *X){
  int j;
  DOUBLE _deriv;
  X[0]=tnk.deriv;utd(X);_deriv=tnk.deriv;
  if(tnk.BC!=3)for(j=0;j<8;j++){corrac();X[0]=tnk.deriv;utd(X);}
  if(testpositivity())return(1);
  if(tnk.BC!=3)dua(X);
  if(tnk.BC!=3)for(j=0;j<8;j++){corrac();X[0]=tnk.deriv;utd(X);}
  if(testpositivity())return(1);
  tnk.jj=scat.nn;utdwrac(X);
  dpsiRpardr();
  nnadd();utdadd(X);
  if(tnk.dpsiR<=0){
    Error("\nFATAL ERROR: hydrodynamic part of trial wave function cannot be constructed");
  }if(fabs(_deriv-tnk.deriv)>fabs(_deriv)*0.1){
    Error("\nFATAL ERROR: bad precision of calculating trial wave function");
  }if(fabs(_deriv-tnk.deriv)>fabs(_deriv)*1E-5){
    Warning("\nWarning: the precision of calculating trial wave function is not good");
  }return(0);
}
// // // // // // // // // // // // // // // // // // // // // // // // // //
//                    block for a joint with hydrodynamics
//
// calculating xi
void calcxi(void){
  tnk.xi=0;
  if(DIMENSION==1)tnk.xi=tnk.dpsiR/(1/Rpar-1/(Lwf-Rpar));
  if(DIMENSION==2)tnk.xi=tnk.dpsiR/(Lwf/(Rpar*Rpar)-Lwf/((Lwf-Rpar)*(Lwf-Rpar)));
  if(DIMENSION==3)tnk.xi=tnk.dpsiR/(Lwf*Lwf/(Rpar*Rpar*Rpar)-Lwf*Lwf/((Lwf-Rpar)*(Lwf-Rpar)*(Lwf-Rpar)));
}
// ln psi(r) for HD part of trial WF
DOUBLE lnpsihd(DOUBLE r){
  DOUBLE y;
  y=0;
  if(r>=Lwf*.5)return(0);
  if(DIMENSION==1)y=-tnk.xi*log(Lwf*Lwf/(4*r*(Lwf-r)));
  if(DIMENSION==2)y=-tnk.xi*(Lwf/r+Lwf/(Lwf-r)-4);
  if(DIMENSION==3)y=-tnk.xi*(Lwf*Lwf*.5/(r*r)+Lwf*Lwf*.5/((Lwf-r)*(Lwf-r))-4);
  return(y);
}
// dot ln psi(r) for HD part of trial WF
DOUBLE dlnpsihd(DOUBLE r){
  DOUBLE y;
  y=0;
  if(r>=Lwf*.5)return(0);
  if(DIMENSION==1)y=tnk.xi*(1/r-1/(Lwf-r));
  if(DIMENSION==2)y=tnk.xi*(Lwf/(r*r)-Lwf/((Lwf-r)*(Lwf-r)));
  if(DIMENSION==3)y=tnk.xi*(Lwf*Lwf/(r*r*r)-Lwf*Lwf/((Lwf-r)*(Lwf-r)*(Lwf-r)));
  return(y);
}
// dot dot ln psi(r) for HD part of trial WF
DOUBLE ddlnpsihd(DOUBLE r){
  DOUBLE y;
  y=0;
  if(r>=Lwf*.5)return(0);
  if(DIMENSION==1)y=-1*tnk.xi*(1/(r*r)+1/((Lwf-r)*(Lwf-r)));
  if(DIMENSION==2)y=-2*tnk.xi*(Lwf/(r*r*r)+Lwf/((Lwf-r)*(Lwf-r)*(Lwf-r)));
  if(DIMENSION==3)y=-3*tnk.xi*(Lwf*Lwf/(r*r*r*r)+Lwf*Lwf/((Lwf-r)*(Lwf-r)*(Lwf-r)*(Lwf-r)));
  return(y);
}
// normalization of numerically calculated trial WF
void mulpsi(void){
  int j;
  DOUBLE mul;
  mul=exp(lnpsihd(Rpar));
  for(j=-tnk.jm;j<=tnk.jj+3;j++)scat.arr[j+tnk.jm]=asinh(mul*sinh(scat.arr[j+tnk.jm]));
}
// 3-th order interpolation of numerically calculated trial WF
DOUBLE lnpsinum(DOUBLE r){
  int j;
  DOUBLE y,t,Dt,dt,*a;
  a=scat.arr;t=asinh(log(r));Dt=(scat.b-scat.a)/tnk.jj;
  j=(int)(tnk.jm+(t-scat.a)/Dt+.5);dt=(t-scat.a)/Dt+tnk.jm-j;
  if(j<3){j=3;dt=-.5;}if(j>tnk.jm+tnk.jj){j=tnk.jm+tnk.jj;dt=.5;}
  y=(((-.0625*(a[j+1]-a[j-1])+1/48.*(a[j+3]-a[j-3]))*dt
    -.25*a[j]+.125*(a[j+2]+a[j-2]))*dt
    +.5625*(a[j+1]-a[j-1])-1/48.*(a[j+3]-a[j-3]))*dt+a[j];
  y=sinh(y);if(y<tnk.doubl)y=tnk.doubl*exp(y/tnk.doubl-1);y=log(y);
  return(y);
}
// constructing grid for ln psi(r)
void constructgridwf(int J,int m){
  int j;
  DOUBLE Dt,x,r;
  if(Rpar>=Lwf/2){
    Error("\nFATAL ERROR: hydrodynamic part of trial wave function is failed");
  }twf.J=J;twf.m=m;twf.j0=2*m+1;twf.nmax=2*twf.j0/3;twf.rmin=tnk.ac;twf.rmax=Rpar;
  if(Rpar<=twf.rmin){
    Error("\nFATAL ERROR: value of Rpar is inside a hard core");
  }twf.Dx=1./(twf.J-2*twf.j0);Dt=(asinh(log(Rpar))-asinh(log(tnk.ac)))/tnk.jj;
  x=(twf.J-twf.j0-.5)/(twf.J-2*twf.j0);twf.rMAX=twf.rmin+x*x*x*x*(twf.rmax-twf.rmin);
  if((int)(tnk.jm-(asinh(log(twf.rMAX))-asinh(log(Rpar)))/Dt+.5)<3){
    Error("\nFATAL ERROR: interpolation of trial wave function at r<Rpar is impossible");
  }for(j=0;j<twf.J;j++){
    x=(j-twf.j0+.5)/(twf.J-2*twf.j0);
    r=twf.rmin+x*x*x*x*(twf.rmax-twf.rmin);
    twf.lnf[j]=lnpsinum(r);
  } }
// function for caclulating acut
DOUBLE fcalcacut(DOUBLE *X){
  int j,step;
  DOUBLE r,x,y;
  r=X[0];
  x=sqrt(sqrt((r-twf.rmin)/(twf.rmax-twf.rmin)));
  j=(int)(twf.j0+x*(twf.J-2*twf.j0))-twf.j0;step=1;
  y=_fxjkx(x-(j+.5)*twf.Dx,twf.Dx,0,twf.m,twf.nmax,twf.j0,step,twf.lnf+j,M8);
  y=exp(y)-tnk.doubl*4;
  return(y);
}
// caclulating acut
void calcacut(void){
  int _l;
  DOUBLE _x;
  _l=0;
  if(fcalcacut(&twf.rmin)*fcalcacut(&twf.rmax)>=0)twf.acut=tnk.ac;
  else twf.acut=_zero(twf.rmin,twf.rmax,_l,&_x,fcalcacut);
  if(Rpar<=twf.acut){
    Error("\nFATAL ERROR: value of Rpar is inside a hard core");
  } }
// // // // // // // // // // // // // // // // // // // // // // // // // //
//                   block for interpolating pair interaction
//
// U(r(x))
DOUBLE fUx(DOUBLE *x){
  DOUBLE r,U;
  r=ugr.rmin+exp(ugr.K*log(*x**x)/2)*(ugr.rmax-ugr.rmin);
  U=InteractionEnergy(r)+1/tnk.Doubl+ugr.U0;
  return(U);
}
// constructing grid for pair interaction U(r)
void constructgridU(void){
  int j,J,Jm,m,j0,jmax,K,_l,stma;
  DOUBLE a,b,U0,err,rmin,rmax,_a,_x,*yj,*ydj;
  m=4;j0=2*m+1;K=32/(1+DIMENSION/2);a=0;b=1;stma=1;_l=0;ugr.m=m;ugr.j0=j0;ugr.K=K;
  _a=exp(1./DIMENSION*log(tnk.doubl*tnk.doubl/n));rmin=tnk.ac;if(rmin<_a)rmin=_a;
  rmax=16*Lwf;U0=4*pi*n;Jm=ugr.memJ;yj=tnk.yj;ydj=yj+Jm;
  ugr.a1=1/(rmax-rmin);ugr.a2=1./K;ugr.rmin=rmin;ugr.rmax=rmax;ugr.U0=U0;
  jmax=200;err=_calcstep(a,b,m,Jm,&jmax,_l,&_x,yj,ydj,M8,fUx);
  // Kurbakov 9-CORE
  if(err>1E-4)
#ifdef KURBAKOV_SOFT_CORE
    Warning("\nWarning: precision for interpolation of pair interaction is BAD!");
#else
    Error("\nFATAL ERROR: interpolation of pair interaction is failed");
#endif
  if(jmax)J=jmax+2*j0+1;else{
    if(err>1E-14)Warning("\nWarning: precision for interpolation of pair interaction is not achieved");
    J=Jm;jmax=J-2*j0-1;
  }Ur.y=ugr.Urj=(DOUBLE*)calloc(J,sizeof(DOUBLE));
  _splineprep(a,b,jmax,m,stma,J,&Ur);_splinecalc(_l,&_x,Ur,fUx);
  ugr.jmax=jmax;ugr.J=J;ugr.err=err;for(j=0;j<J;j++)Ur.y[j]-=U0;
}
// numeric interpolation of pair interaction potential
DOUBLE U0rnum(DOUBLE r){
  DOUBLE y,x;
  if(r<ugr.rmin)r=ugr.rmin;if(r>ugr.rmax)r=ugr.rmax;
  x=exp(ugr.a2*log((r-ugr.rmin)*ugr.a1+1./tnk.Doubl));
  y=_spline9(x,Ur);
  return(y);
}
// // // // // // // // // // // // // // // // // // // // // // // // // //
//                               final block
//
//#ifdef TRIAL_KURBAKOV_NUMERIC
// executed only if TRIAL_KURBAKOV_NUMERIC is defined
DOUBLE InterpolateKurbakovU(DOUBLE r) {
  // U=ln f
  int j,step;
  DOUBLE x,y;
  if(twf.acut<r&&r<twf.rmax){
    x=sqrt(sqrt((r-twf.rmin)/(twf.rmax-twf.rmin)));
    j=(int)(twf.j0+x*(twf.J-2*twf.j0))-twf.j0;step=1;
    y=_fxjkx(x-(j+.5)*twf.Dx,twf.Dx,0,twf.m,twf.nmax,twf.j0,step,twf.lnf+j,M8);
    return(y);
  }if(r<=twf.acut)return(-709);
  if(r>Lhalfwf)return(0);
  return(lnpsihd(r));
}
DOUBLE InterpolateKurbakovFp(DOUBLE r) {
  // Fp=f'/f=U'=U'(x)/r'(x)
  int j,step;
  DOUBLE x,y;
  if(twf.acut<r&&r<twf.rmax){
    x=sqrt(sqrt((r-twf.rmin)/(twf.rmax-twf.rmin)));
    j=(int)(twf.j0+x*(twf.J-2*twf.j0))-twf.j0;step=1;
    y=_fxjkx(x-(j+.5)*twf.Dx,twf.Dx,1,twf.m,twf.nmax,twf.j0,step,twf.lnf+j,M8);
    y/=4*x*x*x*(twf.rmax-twf.rmin);
    return(y);
  }if(r<=twf.acut||r>Lhalfwf)return(0);
  return(dlnpsihd(r));
}
DOUBLE InterpolateKurbakovE(DOUBLE r) {
  // Eloc=Fp^2-f"/f+(1-D)*Fp/r=(1-D)*U'/r-U", D=1,2,3,
  // Fp^2-f"/f=(r"(x)*Fp(x)-U"(x))/(r'(x))^2, Fp(x)=U'(x)/r'(x)
  int j,step;
  DOUBLE x,y,r1,r2,y1,y2;
  if(twf.acut<r&&r<twf.rmax){
    x=sqrt(sqrt((r-twf.rmin)/(twf.rmax-twf.rmin)));
    j=(int)(twf.j0+x*(twf.J-2*twf.j0))-twf.j0;step=1;
    r1=4*x*x*x*(twf.rmax-twf.rmin);
    r2=12*x*x*(twf.rmax-twf.rmin);
    y1=_fxjkx(x-(j+.5)*twf.Dx,twf.Dx,1,twf.m,twf.nmax,twf.j0,step,twf.lnf+j,M8);
    y2=_fxjkx(x-(j+.5)*twf.Dx,twf.Dx,2,twf.m,twf.nmax,twf.j0,step,twf.lnf+j,M8);
    y=(r2*y1/r1-y2)/r1/r1+(1-DIMENSION)*y1/r1/r;
    return(y);
  }if(r<=twf.acut||r>Lhalfwf)return(0);
  y=(1-DIMENSION)*dlnpsihd(r)/r-ddlnpsihd(r);
  return(y);
}
// Kurbakov A-FESHBACH
// calculation of dressed interaction potential
DOUBLE DressedInteractionPotential(DOUBLE p){
  int j,jm,dec,l,lm;
  DOUBLE r,x,der,a1,U,J,y;
  jm=16384;dec=256;lm=32;for(y=j=0;j<jm;j++){
    r=Lwf/2*exp(dec*(j-jm+.5)/jm);der=2*pi*r*dec*r/jm;x=p*r;if(r>twf.acut){
      U=InteractionEnergy(r)*exp(1*InterpolateKurbakovU(r));
      J=0;if(DIMENSION==1)J=cos(x);if(DIMENSION==3){if(x)J=sin(x)/x;else J=1;}if(DIMENSION==2){
        if(x<=3./64){a1=.25*x*x;J=1-a1*(1-a1*(1./4-a1*(1./36)));}else{
          if(x<=1){a1=.25*x*x;J=1-a1*(1-a1*(1./4-a1*(1./36-a1*(1./576-a1*(1./14400-a1*(1./518400-a1*(1./5040/5040-a1*(1./40320/40320))))))));}else{
            if(x<=50)for(J=l=0;l<lm;l++)J+=cos(x*cos(pi/2*(l+.5)/lm))/lm;else{
              a1=1./(x*x);J=sqrt(2./(pi*x))*(
                cos(x-pi*.25)*(1.+(-9/128.+(105*105/24./4096+(-10395*10395/720./262144.+10395*10395.*195*195/40320./4096./4096.*a1)*a1)*a1)*a1)+
                sin(x-pi*.25)*.125*a1*x*(1+(-75/128.+(945*945/120./4096+(-945*945.*143*143/262144./5040.+10395*10395.*195*195*289/362880./4096./4096.*a1)*a1)*a1)*a1));
            } } } }y+=U*J*der;
    } }return(y);
}

void SavingSplineArrays(void){
  FILE *ip;
  int j;
  DOUBLE t,x,r;
// saving u(t)
  ip=fopen("scat.dat","w");
  fprintf(ip,"#x*1E4\t\tnumer*1E4\tHD*1E4\t\tpsi(x)*1E4\n");
  for(j=scat.nn;j>=-tnk.jm;j--){
    t=scat.a+j*(scat.b-scat.a)/scat.nn;x=exp(sinh(t));
    fprintf(ip,"%lf\t%lf\t%lf\t%lf\n",x*1E4, sinh(scat.arr[j+tnk.jm])*1E4,exp(lnpsihd(x))*1E4,exp(lnpsinum(x))*1E4);
  }fclose(ip);
// saving lnf[]
  ip=fopen("lnf.dat","w");
  fprintf(ip,"#j\tr_j\t\tlnf[j]*1E4\n");
  for(j=twf.j0;j<twf.J-twf.j0;j++){
    x=(j-twf.j0+.5)/(twf.J-2*twf.j0);
    r=twf.rmin+x*x*x*x*(twf.rmax-twf.rmin);
    fprintf(ip,"%d\t%le\t%lf\n",j,r,twf.lnf[j]*1E4);
  }fclose(ip);
// saving Urj[]
  ip=fopen("urj.dat","w");
  fprintf(ip,"#j\tr_j\t\tUrj[j]\t\terr\n");for(j=ugr.j0;j<ugr.J-ugr.j0;j++){
    x=(j-ugr.j0)*Ur.dx;r=ugr.rmin+exp(ugr.K*log(x*x)/2)*(ugr.rmax-ugr.rmin);
    fprintf(ip,"%d\t%le\t%le\t%le\n",j,r,ugr.Urj[j],InteractionEnergy(r)-ugr.Urj[j]);
  }fclose(ip);
// saving trial WF
  ip=fopen("wfdat.dat","w");
  fprintf(ip,"#r\t\tU\t\tFp\t\tEloc\t\tEloc-Fp^2+U(r)\tU(r)-Unumeric(r)\n");
  for(j=1;j<=10000;j++){
    r=Lwf/2*j/10000.;fprintf(ip,"%le\t%le\t%le\t%le\t%le\t%le\n",r,
      InterpolateKurbakovU(r),InterpolateKurbakovFp(r),InterpolateKurbakovE(r),
      InterpolateKurbakovE(r)-InterpolateKurbakovFp(r)*InterpolateKurbakovFp(r)+InteractionEnergy(r),
      InteractionEnergy(r)-U0rnum(r));
  }fclose(ip);
}
//#endif

void ConstructGridKurbakovNumeric(struct Grid *G){
  int er;
  AllocateWFGrid(G, 3);
  Message("Kurbakov Numeric trial wavefunction\n");
  Message("  Rpar %" LG "  [Lwf/2] -> %" LG "  [R]\n", Rpar, Rpar*Lhalfwf);
  Message("  Epar %" LG " -> %" LG "\n", Epar, Epar*4*pi*n);
  Rpar *= Lhalfwf;Epar*=4*pi*n;
  DOUBLE *X;
  //_teacher(2,5);getchar();
  _extrmatr(8,&M8);_rungematr(2,&scat);
  tnk.memX=2;tnk.memmin=8192;tnk.memmax=131072-tnk.memX;tnk.memadd=2048;
  tnk.Doubl=1E250;tnk.doubl=1E-9;twf.memJ=512;ugr.memJ=8192/(1+DIMENSION/2);
  twf.mulEpar=1;twf.DmulEpar=1.001;
  scat.arr=(DOUBLE*)calloc(tnk.memmax+tnk.memX,sizeof(DOUBLE));
  X=scat.arr;scat.arr+=tnk.memX;
  twf.lnf=(DOUBLE*)calloc(twf.memJ,sizeof(DOUBLE));
  tnk.yj=(DOUBLE*)calloc(ugr.memJ*2,sizeof(DOUBLE));
  if(a<0.) { // attractive pseudopotential
    tnk.aa=-a; // s-wave scattering length f'(r)/f(r) = -1/a
    Message("  Pseudopotential with s-wave scattering length a = %lf\n", a);
    tnk.ac=0.; // zero hard-core size
  }
  else {
    // Kurbakov 9-CORE
#ifdef KURBAKOV_SOFT_CORE
    Kurbakov_softcore=a;
    Message("  Kurbakov soft core with diameter a = %lf\n", a);
    tnk.ac=a=0;
#else
    Message("  Hard-core potential with diameter a = %lf\n", a);
    //tnk.ac=hard_core_diameter; // hard-core size
    tnk.ac=a; // hard-core diameter equal to a
#endif
  }
  if(Kurbakov_BC) { // if not defined in MC.cfg use main.h
    tnk.BC=Kurbakov_BC;
  }
  else { // 1 - f(0)=0, 2 - u(t->-inf)/dt = 0, f(r) = sinh u (asinh (ln (r)))
#ifdef JASTROW_LEFT_BOUNDARY_ZERO // f(0) = 0
    tnk.BC=1;
#endif
#ifdef JASTROW_LEFT_BOUNDARY_U_MINUS_INF_ZERO  // u(t->-inf)/dt = 0, f(r) = sinh u (asinh (ln (r)))
    tnk.BC=2;
#endif
    //#ifdef JASTROW_LEFT_BOUNDARY_ZERO_DERIVATIVE // f'(0) = 0
    //    tnk.BC=0;
    //#endif
#ifdef JASTROW_LEFT_BOUNDARY_PSEUDOPOTENTIAL // (1D) f'(0)/f(0) = -1/a
    tnk.BC=3;
#endif
  }

  do{
    tunnexp();
    Message("\nac=%g, ",tnk.ac);
    dua(X);
    Message("deriv=%g",tnk.deriv);
    er=acderivdpsiRj0(X);if(er){
      Epar/=twf.DmulEpar;twf.mulEpar/=twf.DmulEpar;twf.DmulEpar*=twf.DmulEpar;
    } }while(twf.mulEpar>0.99&&er);
    Message("\nac=%g,deriv=%g,dpsiR=%g",tnk.ac,tnk.deriv,tnk.dpsiR);
    calcxi();
    mulpsi();
    Message("\nxi=%g; ",tnk.xi);
    constructgridwf(twf.memJ,5);
    calcacut();
    Message("rmin_wf=%g,rMAX_wf=%g,acut=%g\n",twf.rmin,twf.rMAX,twf.acut);
    Message("  Minimal distance in w.f. will be %lf\n", twf.acut);
    G->min = twf.acut;
    G->max = Lhalfwf;
    G->max2 = G->max*G->max;
    a=twf.acut;
    a2=a*a;
    constructgridU();
    Message("J=%d,rmin_U=%g,err=%g\n",ugr.J,ugr.rmin,ugr.err);
    SavingSplineArrays();
    // Kurbakov A-FESHBACH
    Message("dressed interaction: g=%g,n0/n=%g,n*xi^dim=%g\n",tnk.g,tnk.n0,tnk.Nxi);
}
