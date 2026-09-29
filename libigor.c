#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"
#include "libigor.h"
#include "utils.h"
#include "compatab.h"

//DOUBLE pi=3.14159265358979323846264338327950288
//DOUBLE gammaE=0.577215664901532860606512090082402431042
//DOUBLE pi=3.14159265358979323846;
//DOUBLE gammaE=0.57721566490153286;
struct exint M8;

/************************** Construct Kurbakov's library **********************/
// // // // // // // // // // // // // // // // // // // // // // // // // //
//                                block for matrices
//
// binominal coefficients n!/m!/(n-m)!
DOUBLE _Cnm(int n,int m){
  DOUBLE y;
  int k;
  y=0;if(m>=0&&m<=n)for(y=1,k=1;k<=m;k++){y*=(n-k+1.);y/=k;}
  return(y);
}
// n^m
DOUBLE _po(DOUBLE n,int m){
  int k;
  DOUBLE p;
  if(m<0){n=1./n;m=-m;}
  for(p=1,k=0;k<m;k++)p*=n;
  return(p);
}
// n! (n>=0)
DOUBLE _fa(int n){
  int k;
  DOUBLE f;
  for(f=1,k=1;k<=n;k++)f*=k;
  return(f);
}
// calculating D_n^l
DOUBLE _Dnlfunc(int n,int j_,struct exint e){
  int j,k,l,sig;
  DOUBLE y;
  j=j_-1;
  for(sig=-1,y=0,k=0;k<=n+2*j;k++){
    sig=-sig;for(l=0;l<=j;l++){
      y+=sig*e.Dnl[n*e.M+l]*_Cnm(n+2*l,k-j+l)*_po(n/2.+j-k,n+2*j+2)/_fa(n+2*j+2);
  } }return(y);
}
// calculating A_kn^m
DOUBLE _Aknmfunc(int k,int n,int m,struct exint e){
  int l;
  DOUBLE y;
  if(n){
    for(y=0,l=0;l<=m;l++)y+=(DOUBLE)(e.Dnl[n*e.M+l]*_Cnm(n+2*l,k-m+l));
    for(l=0;l<k+m;l++)y=-y;
  }else{
    if(k==m)y=1;else y=0;
  }return(y);
}
// calculation of W_k^m
DOUBLE _Wkmfunc(int k,int m,struct exint e){
  int n;
  DOUBLE c0,c1;
  if(k>m)k=2*m-k;
  for(c0=e.Aknm[k+2*0*(e.M+1)+m*(e.M+1)*2*e.M],c1=1,n=1;n<=m;n++){
    c1/=8.*n*(2*n+1);c0+=e.Aknm[k+2*n*(e.M+1)+(m-n)*(e.M+1)*2*e.M]*c1;
  }return(c0);
}
// calculation of I_l^m
DOUBLE _Ilmfunc(int l,int m,struct exint e){
  int k;
  DOUBLE c0;
  for(c0=0,k=0;k<=l;k++)if(k<=2*m)c0+=e.Wkm[(2*m-k)*e.M+m];
  return(c0);
}
// calculation of B_l^m
DOUBLE _Blmfunc(int l,int m,struct exint e){
  int n,sig;
  DOUBLE y;
  for(sig=1,y=0,n=0;n<=m;n++){
    sig=-sig;if(n>=l){
      y+=sig*_Cnm(m,n)*_Cnm(2*n,n-l)*_po(1/4.,n);
  } }return(y);
}
// calculation of c_l for dual grid
void _cldu(struct exint e){
  int k,n;
  DOUBLE ss;
  e.cl_du[0]=1;for(n=1;n<e.M;n++){
    for(ss=k=0;k<n;k++)ss+=e.cl_du[k]*_Cnm(2*n,2*k);e.cl_du[n]=-(DOUBLE)(ss);
} }
// calculation of \tilde D_k^m
DOUBLE _Dkmfunc(int k,int m,struct exint e){
  int l;
  DOUBLE mul,y;
  if(k>m)return(0);if(!m)return(0.5);
  for(mul=0.5,y=l=0;l<=m;l++){
    y+=e.cl_du[l]*e.Aknm[m-k+(e.M+1)*2*l+(e.M+1)*2*e.M*(m-l)]*mul;
    if(k<m)y+=e.cl_du[l]*e.Aknm[m-k-1+(e.M+1)*2*l+(e.M+1)*2*e.M*(m-l)]*mul;
    mul/=4*(2*l+1)*(2*l+2);
  }return(y);
}
// calculation of matrices for D_n^l and A_kn^m and array of values 1./n
void _extrmatr(int M,struct exint *e){
  int k,l,m,n;
  e->M=M;
  e->Dnl=(DOUBLE*)calloc(2*M*M,sizeof(DOUBLE));
  e->Aknm=(DOUBLE*)calloc(2*M*(M+1)*M,sizeof(DOUBLE));
  e->Wkm=(DOUBLE*)calloc(2*M*M,sizeof(DOUBLE));
  e->Ilm=(DOUBLE*)calloc(2*M*M,sizeof(DOUBLE));
  e->Blm=(DOUBLE*)calloc(M*M,sizeof(DOUBLE));
  e->Dkm=(DOUBLE*)calloc(M*M,sizeof(DOUBLE));
  e->cl_du=(DOUBLE*)calloc(M,sizeof(DOUBLE));
  e->v1n=(DOUBLE*)calloc(64,sizeof(DOUBLE));
  e->yxj=(DOUBLE*)calloc(64,sizeof(DOUBLE));
  for(n=1;n<2*M;n++)for(e->Dnl[n*M+0]=1,l=1;l<M;l++)e->Dnl[n*M+l]=_Dnlfunc(n,l,*e);
  for(k=0;k<=M;k++)for(n=0;n<2*M;n++)for(m=0;m<M;m++)e->Aknm[k+n*(M+1)+m*(M+1)*2*M]=_Aknmfunc(k,n,m,*e);
  for(k=0;k<2*M;k++)for(m=0;m<M;m++)e->Wkm[k*M+m]=_Wkmfunc(k,m,*e);
  for(l=0;l<2*M;l++)for(m=0;m<M;m++)e->Ilm[l*M+m]=_Ilmfunc(l,m,*e);
  for(l=0;l<M;l++)for(m=0;m<M;m++)e->Blm[l*M+m]=_Blmfunc(l,m,*e);
  _cldu(*e);for(m=0;m<M;m++)for(k=0;k<M;k++)e->Dkm[k*M+m]=_Dkmfunc(k,m,*e);
  for(n=1;n<64;n++)e->v1n[n]=1./n;
}
// // // // // // // // // // // // // // // // // // // // // // // // // //
//               block for differentiation and interpolation
//
// quick direct calculation of n-th derivative of f(x) in (m-n/2)-th order
// array f(x-j0*step*dx),f(x-(j0*step-1)*dx),...f(x+j0*step*dx) is given by
//   fxj[0],...fxj[j0*step*2] (hence, fxj[j0*step]=f(x));
// one must be 0<=n/2<=m<=7 and 1+2m<=j0<=15
DOUBLE _fxjnmx(DOUBLE dx,int n,int m,int j0,int step,DOUBLE *fxj,struct exint e){
  int k,i,jj,si;
  DOUBLE ynx,c1;
  m-=n>>1;if(m<0)m=0;// order of differentiation is m-n/2
  if(n+2*m>j0||n<0||m>7||j0>15){
    Error("\nFATAL ERROR in fxjnmx: n=%d, m=%d, j0=%d, step=%d",n,m+n/2,j0,step);getchar();
  }j0*=step;if(n==0)return(fxj[j0]);
  if(n==2*(n>>1))si=1;else si=-1;// si=(-1)^n
  i=(e.M+1)*(n+2*e.M*m);jj=(n+2*m)*step;
  for(ynx=k=0;k<=(n+2*m-1)>>1;k++,jj-=step+step)// lateral points
    ynx+=e.Aknm[k+i]*(fxj[j0+jj]+si*fxj[j0-jj]);
  if(si==1)ynx+=e.Aknm[k+i]*fxj[j0];// central point (if there is)
  dx*=2*step;for(c1=dx,k=1;k<n;k++)c1*=dx;// division by (2*dx)^n
  return(ynx/c1);
}
// quick Tailor series for k-th derivative of f(x) in m-th order at x near to x0
// grid f(x0-j0*step*dx),f(x0-(j0*step-1)*dx),...f(x0+j0*step*dx) is given by
//   fxj[0],...fxj[j0*step*2] (hence, fxj[j0*step]=f(x0), where x0 is the nearest site to x)
// one must be 0<=nmax/2<=m<=7 and 2*m+1<=j0<=15
DOUBLE _fxjkx(DOUBLE xminusx0,DOUBLE dx,int k,int m,int nmax,int j0,int step,DOUBLE *fxj,struct exint e){
  int n;
  DOUBLE f,c1;
  if(2*m+1>j0||nmax/2>m||nmax<0||m>7||j0>15){
    Error("\nFATAL ERROR in fxjkx: nmax=%d, m=%d, j0=%d, step=%d",nmax,m,j0,step);getchar();
  }f=_fxjnmx(dx,k,m,j0,step,fxj,e);
  if(xminusx0)for(c1=1,n=1;n<=nmax-k;n++){
    c1*=xminusx0*e.v1n[n];f+=_fxjnmx(dx,n+k,m,j0,step,fxj,e)*c1;
  }return(f);
}
// ultraquick calculation of 5-th order spline for function f(x) (21 multiplications)
// f(x0)=y[0],f(x0+dx)=y[1],f(x0+2*dx)=y[2],...f(x0+jm*dx)=y[jm];   _dx=1/dx
// one must be 4.5*dx<x-x0<(jm-4.5)*dx and jm>=10
DOUBLE _spline5(DOUBLE x,struct spline s){
  int j;
  DOUBLE xx,*y;
  xx=(x-s.x0)*s._dx;j=(int)(xx+.5);xx-=j;y=s.y;
  if(5<=j&&j<=s.jm-5)
    return(y[j]+
      xx*(3/1280.*(y[j+5]-y[j-5])-25/768.*(y[j+3]-y[j-3])+75/128.*(y[j+1]-y[j-1])+
          xx*(-1/96.*(y[j+4]+y[j-4])+1/6.*(y[j+2]+y[j-2])-5/16.*y[j]+
              xx*(-1/384.*(y[j+5]-y[j-5])+13/384.*(y[j+3]-y[j-3])-17/192.*(y[j+1]-y[j-1])+
                  xx*(1/384.*(y[j+4]+y[j-4])-1/96.*(y[j+2]+y[j-2])+1/64.*y[j]+
                      xx*(1/3840.*(y[j+5]-y[j-5])-1/768.*(y[j+3]-y[j-3])+1/384.*(y[j+1]-y[j-1])))))));
  Error("\nFATAL ERROR in spline5: x=%g is out of the range of array y[] (j=%d, jm=%d)",x,j,s.jm);
  getchar();return(0);
}
DOUBLE _deriv5(DOUBLE x,struct spline s){
  int j;
  DOUBLE xx,*y;
  xx=(x-s.x0)*s._dx;j=(int)(xx+.5);xx-=j;y=s.y;
  if(5<=j&&j<=s.jm-5)return(
      s._dx*(3/1280.*(y[j+5]-y[j-5])-25/768.*(y[j+3]-y[j-3])+75/128.*(y[j+1]-y[j-1])+
          xx*(-1/48.*(y[j+4]+y[j-4])+1/3.*(y[j+2]+y[j-2])-5/8.*y[j]+
              xx*(-1/128.*(y[j+5]-y[j-5])+13/128.*(y[j+3]-y[j-3])-17/64.*(y[j+1]-y[j-1])+
                  xx*(1/96.*(y[j+4]+y[j-4])-1/24.*(y[j+2]+y[j-2])+1/16.*y[j]+
                      xx*(1/768.*(y[j+5]-y[j-5])-5/768.*(y[j+3]-y[j-3])+5/384.*(y[j+1]-y[j-1])))))));
  Error("\nFATAL ERROR in deriv5: x=%g is out of the range of array y[] (j=%d, jm=%d)",x,j,s.jm);
  getchar();return(0);
}
// ultraquick calculation of 7-th order spline for function f(x) (31 multiplications)
// f(x0)=y[0],f(x0+dx)=y[1],f(x0+2*dx)=y[2],...f(x0+jm*dx)=y[jm];   _dx=1/dx
// one must be 6.5*dx<x-x0<(jm-6.5)*dx and jm>=14
DOUBLE _spline7(DOUBLE x,struct spline s){
  int j;
  DOUBLE xx,*y;
  xx=(x-s.x0)*s._dx;j=(int)(xx+.5);xx-=j;y=s.y;
  if(7<=j&&j<=s.jm-7)
    return(y[j]+
      xx*(-5/14336.*(y[j+7]-y[j-7])+49/10240.*(y[j+5]-y[j-5])-245/6144.*(y[j+3]-y[j-3])+1225/2048.*(y[j+1]-y[j-1])+
          xx*(1/720.*(y[j+6]+y[j-6])-3/160.*(y[j+4]+y[j-4])+3/16.*(y[j+2]+y[j-2])-49/144.*y[j]+
              xx*(37/92160.*(y[j+7]-y[j-7])-499/92160.*(y[j+5]-y[j-5])+1299/30720.*(y[j+3]-y[j-3])-1891/18432.*(y[j+1]-y[j-1])+
                  xx*(-1/2304.*(y[j+6]+y[j-6])+1/192.*(y[j+4]+y[j-4])-13/768.*(y[j+2]+y[j-2])+7/288.*y[j]+
                      xx*(-1/18432.*(y[j+7]-y[j-7])+59/92160.*(y[j+5]-y[j-5])-5/2048.*(y[j+3]-y[j-3])+83/18432.*(y[j+1]-y[j-1])+
                          xx*(1/46080.*(y[j+6]+y[j-6])-1/7680.*(y[j+4]+y[j-4])+1/3072.*(y[j+2]+y[j-2])-1/2304.*y[j])))))));
  Error("\nFATAL ERROR in spline7: x=%g is out of the range of array y[] (j=%d, jm=%d)",x,j,s.jm);
  getchar();return(0);
}
DOUBLE _deriv7(DOUBLE x,struct spline s){
  int j;
  DOUBLE xx,*y;
  xx=(x-s.x0)*s._dx;j=(int)(xx+.5);xx-=j;y=s.y;
  if(7<=j&&j<=s.jm-7)return(
      s._dx*(-5/14336.*(y[j+7]-y[j-7])+49/10240.*(y[j+5]-y[j-5])-245/6144.*(y[j+3]-y[j-3])+1225/2048.*(y[j+1]-y[j-1])+
          xx*(1/360.*(y[j+6]+y[j-6])-3/80.*(y[j+4]+y[j-4])+3/8.*(y[j+2]+y[j-2])-49/72.*y[j]+
              xx*(37/30720.*(y[j+7]-y[j-7])-499/30720.*(y[j+5]-y[j-5])+1299/10240.*(y[j+3]-y[j-3])-1891/6144.*(y[j+1]-y[j-1])+
                  xx*(-1/576.*(y[j+6]+y[j-6])+1/48.*(y[j+4]+y[j-4])-13/192.*(y[j+2]+y[j-2])+7/72.*y[j]+
                      xx*(-5/18432.*(y[j+7]-y[j-7])+59/18432.*(y[j+5]-y[j-5])-25/2048.*(y[j+3]-y[j-3])+415/18432.*(y[j+1]-y[j-1])+
                          xx*(1/7680.*(y[j+6]+y[j-6])-1/1280.*(y[j+4]+y[j-4])+1/512.*(y[j+2]+y[j-2])-1/384.*y[j])))))));
  Error("\nFATAL ERROR in deriv7: x=%g is out of the range of array y[] (j=%d, jm=%d)",x,j,s.jm);
  getchar();return(0);
}
// ultraquick calculation of 9-th order spline for function f(x) (37 multiplications)
// f(x0)=y[0],f(x0+dx)=y[1],f(x0+2*dx)=y[2],...f(x0+jm*dx)=y[jm];   _dx=1/dx
// one must be 8.5*dx<x-x0<(jm-8.5)*dx and jm>=18
DOUBLE _spline9(DOUBLE x,struct spline s){
  int j;
  DOUBLE xx,*y;
  xx=(x-s.x0)*s._dx;j=(int)(xx+.5);xx-=j;y=s.y;
  if(9<=j&&j<=s.jm-9)
    return(y[j]+
      xx*(35/589824.*(y[j+9]-y[j-9])-405/458752.*(y[j+7]-y[j-7])+567/81920.*(y[j+5]-y[j-5])-735/16384.*(y[j+3]-y[j-3])+19845/32768.*(y[j+1]-y[j-1])+
          xx*(-1/4480.*(y[j+8]+y[j-8])+1/315.*(y[j+6]+y[j-6])-1/40.*(y[j+4]+y[j-4])+1/5.*(y[j+2]+y[j-2])-205/576.*y[j]+
              xx*(-3229/46448640.*(y[j+9]-y[j-9])+589/573440.*(y[j+7]-y[j-7])-227/28672.*(y[j+5]-y[j-5])+26611/552960.*(y[j+3]-y[j-3])-4561/40960.*(y[j+1]-y[j-1])+
                  xx*(7/92160.*(y[j+8]+y[j-8])-1/960.*(y[j+6]+y[j-6])+169/23040.*(y[j+4]+y[j-4])-61/2880.*(y[j+2]+y[j-2])+91/3072.*y[j]+
                      xx*(47/4423680.*(y[j+9]-y[j-9])-221/1474560.*(y[j+7]-y[j-7])+377/368640.*(y[j+5]-y[j-5])-1229/368640.*(y[j+3]-y[j-3])+4307/737280.*(y[j+1]-y[j-1])+
                          xx*(-1/184320.*(y[j+8]+y[j-8])+1/15360.*(y[j+6]+y[j-6])-13./46080.*(y[j+4]+y[j-4])+29./46080.*(y[j+2]+y[j-2])-5/6144.*y[j])))))));
  Error("\nFATAL ERROR in spline9: x=%g is out of the range of array y[] (j=%d, jm=%d)",x,j,s.jm);
  getchar();return(0);
}
DOUBLE _deriv9(DOUBLE x,struct spline s){
  int j;
  DOUBLE xx,*y;
  xx=(x-s.x0)*s._dx;j=(int)(xx+.5);xx-=j;y=s.y;
  if(9<=j&&j<=s.jm-9)return(
      s._dx*(35/589824.*(y[j+9]-y[j-9])-405/458752.*(y[j+7]-y[j-7])+567/81920.*(y[j+5]-y[j-5])-735/16384.*(y[j+3]-y[j-3])+19845/32768.*(y[j+1]-y[j-1])+
          xx*(-1/2240.*(y[j+8]+y[j-8])+2/315.*(y[j+6]+y[j-6])-1/20.*(y[j+4]+y[j-4])+2/5.*(y[j+2]+y[j-2])-205/288.*y[j]+
              xx*(-3229/15482880.*(y[j+9]-y[j-9])+1767/573440.*(y[j+7]-y[j-7])-681/28672.*(y[j+5]-y[j-5])+26611/184320.*(y[j+3]-y[j-3])-13683/40960.*(y[j+1]-y[j-1])+
                  xx*(7/23040.*(y[j+8]+y[j-8])-1/240.*(y[j+6]+y[j-6])+169/5760.*(y[j+4]+y[j-4])-61/720.*(y[j+2]+y[j-2])+91/768.*y[j]+
                      xx*(47/884736.*(y[j+9]-y[j-9])-221/294912.*(y[j+7]-y[j-7])+377/73728.*(y[j+5]-y[j-5])-1229/73728.*(y[j+3]-y[j-3])+4307/147456.*(y[j+1]-y[j-1])+
                          xx*(-1/30720.*(y[j+8]+y[j-8])+1/2560.*(y[j+6]+y[j-6])-13./7680.*(y[j+4]+y[j-4])+29./7680.*(y[j+2]+y[j-2])-5/1024.*y[j])))))));
  Error("\nFATAL ERROR in deriv9: x=%g is out of the range of array y[] (j=%d, jm=%d)",x,j,s.jm);
  getchar();return(0);
}
// analytic prodolzhenie: Tailor series for k-th derivative of f(z) in m-th order
// f(x0)=y[0],f(x0+dx)=y[1],f(x0+2*dx)=y[2],...f(x0+jm*dx)=y[jm];
// one must be j0*step*dx<re-x0<(jm-j0*step)*dx and jm>=2*j0*step,
// where step=1,2,3,...stma, j0=nmax=2*m+1, z=re+i*im, f^(k)(z)=Re+i*Im,
void _spline_anal(DOUBLE re,DOUBLE im,DOUBLE *Re,DOUBLE *Im,int k,
    struct spline s,struct exint e){
  int jm,step,m,j,nmax,n,mm,l,si,il,ij1,ij2;
  DOUBLE x0,dx,_dx,xj,a1,b1,a2,b2,c,C,yjn,*y;
  x0=s.x0;dx=s.dx;_dx=s._dx;jm=s.jm;step=s.step;m=s.m;y=s.y;nmax=2*m+1;
  c=.5*_dx*e.v1n[step];j=(int)((re-x0)*_dx+.5);xj=x0+j*dx;
  if((2*m+1)*step>j||j>jm-(2*m+1)*step||k<0||k>nmax||m<0||m>7||step<1||dx<=0){
    Error("\nFATAL ERROR in spline_anal: one must be 0<=m<=7, 0<=k<=2*m+1, step>=1 and jm>=(4*m+2)*step\n(m=%d, k=%d, step=%d, jm=%d, dx=%g)",m,k,step,jm,dx);getchar();
  }for(*Re=*Im=0,n=k;n<=nmax;n++){
    if(n==2*(n>>1))si=1;else si=-1;// si=(-1)^n is sign
    if(n==k){// first time
      a1=1;b1=0;for(C=1,l=0;l<k;l++)C*=c;
    }else{// a1+i*b1=(z-xj)^(n-k)/(n-k)!, C=1/(2*dx*step)^n
      a2=(re-xj)*a1-im*b1;b2=(re-xj)*b1+im*a1;
      a1=a2*e.v1n[n-k];b1=b2*e.v1n[n-k];C*=c;
    }mm=m-(n>>1);// mm=m-n/2 is order of differentiation (mm>=0)
    if(n==0)yjn=y[j];else{// yjn*C=f^(n)(xj)
      il=(e.M+1)*(n+2*e.M*mm);ij1=j+(n+2*mm)*step;ij2=j-(n+2*mm)*step;
      for(yjn=l=0;l<=(n+2*mm-1)>>1;l++,ij1-=step+step,ij2+=step+step)
        yjn+=e.Aknm[l+il]*(y[ij1]+si*y[ij2]);// sum of lateral points
      if(si==1)yjn+=e.Aknm[l+il]*y[ij1];// plus central point (if there is)
    }*Re+=yjn*C*a1;*Im+=yjn*C*b1;// Tailor series for f^(k)(z)
} }
// // // // // // // // // // // // // // // // // // // // // // // // // //
//                      block for treatment of splines
//
// calculation of step (b-a)/_jmax of function f(x)=f(X[0],X[1],X[2],X[3],...),
//   where x=X[_l], at [a;b] in m-th order
// memory for y[Jmax] and yd[Jmax] must be allocated; error is returned
DOUBLE _calcstep(DOUBLE a,DOUBLE b,int m,int Jmax,int *_jmax,int _l,DOUBLE *X,
    DOUBLE *y,DOUBLE *yd,struct exint e,DOUBLE (f)(DOUBLE *X)){
  int j,k,l,j0,jmax,er;
  DOUBLE inc,rat,x,dx,S,err,err0;
  inc=sqrt(3/2.);rat=1.4;// increment and ratio
  jmax=*_jmax;if(jmax<3)jmax=3;j0=2*m+1;dx=(b-a)/jmax;err=0;do{
    for(j=0;j<=2*j0+jmax;j++){X[_l]=x=a+(j-j0)*dx;y[j]=f(X);}// spline
    for(j=m+1;j<=2*j0-m+jmax;j++)for(yd[j]=k=0;k<=m;k++)//   dual
      yd[j]+=e.Dkm[k*e.M+m]*(y[j-1-k]+y[j+k]);//             spline
    for(err0=err,l=0,j=j0;j<=j0+jmax;j++){// calculation of error
      if(y[j]*y[j+1]>0&&y[j]*y[j-1]>0)er=1;else er=0;// elimination of statistical glitch
      for(S=k=0;k<=m;k++)S+=e.Dkm[k*e.M+m]*(yd[j-k]+yd[j+1+k]);// reciprocal dual transformation
      if(er)err=sqrt((err*err*l+(S/y[j]-1)*(S/y[j]-1))/(l+1));l+=er;// error
    }if(!err0)err0=1E308;// first time
    jmax=(int)(jmax*inc+.5);dx=(b-a)/jmax;// increment for jmax and dx
    if(jmax<=Jmax-2*j0-1)er=0;else er=1;// not enough memory
  }while(err0>err*rat&&!er);
  if(er)*_jmax=0;else *_jmax=jmax;
  return(err);
}
// preparing of a structure for a spline
void _splineprep(DOUBLE a,DOUBLE b,int jmax,int m,int stma,int Jm,struct spline *s){
  int j0,j0st,jm;
  j0=2*m+1;j0st=j0*stma;jm=jmax+2*j0st;
  if(a>=b||jmax<1||m<0||m>7||stma<1||jm+1>Jm){
    Error("\nFATAL ERROR in preparing a spline. One must be:\na<b, jmax>=1, 0<=m<=7, stma>=1 and jm+1<=Jm\n(a=%g, b=%g, jmax=%d, m=%d, stma=%d, j0st=%d, jm=%d, Jm=%d)",a,b,jmax,m,stma,j0st,jm,Jm);getchar();
  }s->m=m;s->j0=j0;s->step=1;s->stma=stma;
  s->j0st=j0st;s->jmax=jmax;s->jm=jm;s->Jm=Jm;
  s->a=a;s->b=b;s->dx=(b-a)/jmax;s->_dx=1./s->dx;
  s->x0=a-j0st*s->dx;
}
// calculating a spline y[j]=f(X[0],X[1],X[2],X[3],...) with X[_l]=x0+j*dx
void _splinecalc(int _l,DOUBLE *X,struct spline s,DOUBLE (f)(DOUBLE *X)){
  int j;
  for(j=0;j<=s.jm;j++){X[_l]=s.x0+j*s.dx;s.y[j]=f(X);}
}
// dual transformation in m-th order: yd[j] is between y[j] and y[j+1]
// y[0],...y[jm] is at input; yd[m+1],...yd[jm-m] is at output
// the corresponding DOUBLE grid is yd[m+1],y[m+1],yd[m+2],...y[jm-m-1],yd[jm-m]
// memory for yd[0],...yd[jm] must be allocated
void _dual(int m,int jm,DOUBLE *y,DOUBLE *yd,struct exint e){
  int j,k;
  for(j=m+1;j<=jm-m;j++)for(yd[j]=k=0;k<=m;k++)
    yd[j]+=e.Dkm[k*e.M+m]*(y[j-1-k]+y[j+k]);
}
// calculation of dual spline yd[] for a spline y[]
void _dualspline(struct spline s,struct exint e){
  _dual(s.m,s.jm,s.y,s.yd,e);
}
// converting a dual spline in s into a new spline S
void _convertdual(struct spline s,struct spline *S){
  if(s.stma<2){Error("\nFATAL ERROR in convertdual: one must be stma>=2 (stma=%d)",s.stma);getchar();}
  S->m=s.m;S->j0=s.j0;S->step=1;S->Jm=s.Jm;S->dx=s.dx;S->_dx=s._dx;
  S->jmax=s.jmax+s.j0;S->a=s.a-s.j0*s.dx/2;S->b=s.b+s.j0*s.dx/2;
  S->stma=s.stma-1;S->j0st=s.j0st-s.j0;S->jm=s.jm-s.j0;
  S->x0=s.x0+s.j0*s.dx/2;S->y=s.yd+s.m+1;
}
// doubling the spline in m-th order:
// 1) yd[m+1]-->y[0],y[m+1]-->y[1],...y[j0]-->y[j0],...
//    y[jm-j0]-->y[2*jm-3*j0],...yd[jm-m]-->y[2*(jm-j0)]
// 2) jmax_new=2*jmax, stma_new=2*stma-1, j0st_new=2*j0st-j0, jm_new=2*(jm-j0),
//    dx_new=dx/2, x0_new=a-j0st_new*dx_new
// y[0],...y[jm] is at input; y[0],...y[jm_new] is at output
// memory for dual spline yd[] must be allocated
void _doublingspline(struct spline *s,struct exint e){
  int j;
  if(2*(s->jm-s->j0)+1>s->Jm){Error("\nFATAL ERROR: doubling a spline is impossible, not enough memory\n(jm=%d, j0=%d, Jm=%d)",s->jm,s->j0,s->Jm);getchar();}
  _dualspline(*s,e);// calculating dual spline
  for(j=0;j<s->jm-s->j0;j++)s->y[1+j]=s->y[s->m+1+j];// shifting a spline by m
  for(j=s->jm-s->j0-1;j>0;j--)s->y[1+2*j]=s->y[1+j];// 2-scattering of spline
  for(j=0;j<=s->jm-s->j0;j++)s->y[2*j]=s->yd[s->m+1+j];// inserting the dual spline
  s->jmax*=2;s->stma=2*s->stma-1;s->j0st=2*s->j0st-s->j0;
  s->jm=2*(s->jm-s->j0);s->dx/=2;s->_dx*=2;s->x0=s->a-s->j0st*s->dx;
}
// caclulating k-th derivative of spline in m-th order
// differentiation is performed with a step of dx*step (step=1,2,3,...stma-1)
// y[0],y[1],...y[jm] is at input; y1[j0*step],y1[j0*step+1],...y1[jm-j0*step] is at output
void _splinederivative(int k,struct spline s,struct exint e){
  int j;
  if(s.step<1||s.step>s.stma-1){Error("\nFATAL ERROR in splinederivative: one must be 1<=step<=stma-1 (step=%d, stma=%d)",s.step,s.stma);getchar();}
  for(j=s.j0*s.step;j<=s.jm-s.j0*s.step;j++)
    s.y1[j]=_fxjnmx(s.dx,k,s.m,s.j0,s.step,s.y+j-s.j0*s.step,e);
}
// converting k-th derivative of spline s into a new spline S
void _convertderivative(struct spline s,struct spline *S){
  S->m=s.m;S->j0=s.j0;S->step=1;S->stma=1;S->j0st=s.j0;S->jmax=s.jmax;
  S->Jm=s.Jm;S->a=s.a;S->b=s.b;S->dx=s.dx;S->_dx=s._dx;
  S->jm=s.jm-2*s.j0*s.step;S->x0=s.x0+s.j0*s.step*s.dx;S->y=s.y1+s.j0*s.step;
}
// calculating errors in spline (in m-th order with averaging over 2*Jer+1 values)
// y[0],y[1],...y[jm] is at input; er[Jer+m],er[Jer+m+1],...er[jm-Jer-m] is at output
// total error is returned;
// Jer>=0 and 2*(Jer+m)<=jm; memory for er[0],...er[jm] must be allocated
DOUBLE _splineerror(int Jer,struct spline s,struct exint e){
  int j,l;
  DOUBLE err,ysm,a1,a2;
  if(Jer<0||2*(Jer+s.m)>s.jm){Error("\nFATAL ERROR in splineerror: one must be Jer>=0 and 2*(Jer+m)<=jm (Jer=%d,m=%d,jm=%d)",Jer,s.m,s.jm);getchar();}
  for(err=0,j=s.m;j<=s.jm-s.m;j++){
    ysm=s.y[j]*(1+e.Blm[0*e.M+s.m]);//                                  m-th order
    for(l=1;l<=s.m;l++)ysm+=(s.y[j+l]+s.y[j-l])*e.Blm[l*e.M+s.m];//   smoothing
    err+=(s.y[j]-ysm)*(s.y[j]-ysm);s.er[j-s.m]=s.y[j]-ysm;// preparing errors
  }err=sqrt(err/(s.jm+1-2*s.m)*sqrt(2*s.m*pi));// total error
  a2=sqrt(2*s.m*pi)/(2*Jer+1);for(j=0;j<=s.jm-2*(Jer+s.m);j++){//           averaging errors
    for(a1=l=0;l<=2*Jer;l++)a1+=s.er[j+l]*s.er[j+l];s.er[j]=sqrt(a1*a2);//  over 2*Jer+1 values
  }for(j=s.jm-2*(Jer+s.m);j>=0;j--)s.er[j+Jer+s.m]=s.er[j];// shifting error array by Jer+m
  return(err);
}
// predicting J (-J) values of function f(x) followed (preceded) by its known 2*m values
// grid f(x),f(x+dx),...f(x+2*m*dx) is given by fxj[0],...fxj[2*m]
// at J>0 grid f(x+(2*m+1)*dx)...f(x+(2*m+J)*dx) is returned as fxj[2*m+1]...fxj[2*m+J]
// at J<0 grid f(x-dx),f(x-2*dx)...f(x+J*dx) is returned as fxj[-1],fxj[-2]...fxj[J]
// memory for fxj[2*m+1]...fxj[2*m+J] or fxj[-1],fxj[-2]...fxj[J] must be reserved
void _predict(int J,int m,DOUBLE *fxj,struct exint e){
  int k,j;
  DOUBLE C;
  if(m<1||m>7){Error("\nFATAL ERROR in predict: m=%d must be form 1 to 7",m);getchar();}
  if(J>0)for(j=1;j<=J;j++){
    fxj[2*m+j]=fxj[j-1];C=-2*m-1;
    for(k=1;k<=m;k++,C=-C*(2*m-k+2)*e.v1n[k])
      fxj[2*m+j]+=C*(fxj[j+k-1]-fxj[2*m+j-k]);
  }else for(j=-1;j>=J;j--){
    fxj[j]=fxj[j+2*m+1];C=-2*m-1;
    for(k=1;k<=m;k++,C=-C*(2*m-k+2)*e.v1n[k])fxj[j]+=C*(fxj[j+2*m+1-k]-fxj[j+k]);
} }
// // // // // // // // // // // // // // // // // // // // // // // // // //
//      block for Runge integration, random number generator and other
//
// calculation of coefficients m0,m1,m2,m3 in Runge integration
void _coeffm(DOUBLE *X,struct Runge *r,void (funcs)(int _l,DOUBLE *X,DOUBLE *yl,DOUBLE *fk)){
  int l,k,n,_l;
  DOUBLE dx,xj;
  _l=r->_l;n=r->n;dx=(r->b-r->a)/r->nn;xj=X[_l];
  for(l=1;l<=n;l++)r->yl[l]=r->ykxj[l];
  X[_l]=xj;funcs(_l,X,r->yl,r->fk);for(k=1;k<=n;k++)r->mik0[k]=r->fk[k];
  for(l=1;l<=n;l++)r->yl[l]=r->ykxj[l]+r->mik0[l]*dx/2;
  X[_l]=xj+dx/2;funcs(_l,X,r->yl,r->fk);for(k=1;k<=n;k++)r->mik1[k]=r->fk[k];
  for(l=1;l<=n;l++)r->yl[l]=r->ykxj[l]+r->mik1[l]*dx/2;
  X[_l]=xj+dx/2;funcs(_l,X,r->yl,r->fk);for(k=1;k<=n;k++)r->mik2[k]=r->fk[k];
  for(l=1;l<=n;l++)r->yl[l]=r->ykxj[l]+r->mik2[l]*dx;
  X[_l]=xj+dx;funcs(_l,X,r->yl,r->fk);for(k=1;k<=n;k++)r->mik3[k]=r->fk[k];
}
// Runge integration of n equations from x=a to x=b with step dx=(b-a)/nn
void _runge(DOUBLE *X,struct Runge *r,
    void (funcs)(int _l,DOUBLE *X,DOUBLE *yl,DOUBLE *fk)){
  int j,k,s,n,nn;
  DOUBLE a,b,xj,dx;
  a=r->a;b=r->b;n=r->n;nn=r->nn;s=r->s;dx=(b-a)/nn;
  for(k=1;k<=n;k++){
    if(k<=s)r->arr[(k-1)*(nn+1)]=r->ykxj[k]=r->yk0[k];
    else r->ykxj[k]=r->yk0[k];
  }for(j=1,xj=a;j<=nn;xj=a+j*dx,j++){
    X[r->_l]=xj;_coeffm(X,r,funcs);for(k=1;k<=n;k++){
      r->ykxj[k]+=dx/6*(r->mik0[k]+2*r->mik1[k]+2*r->mik2[k]+r->mik3[k]);
      if(k<=s)r->arr[j+(k-1)*(nn+1)]=r->ykxj[k];
} } }
// memory for Runge matrices
void _rungematr(int n,struct Runge *r){
  r->n=n;r->s=0;
  r->yl=(DOUBLE*)calloc(n+1,sizeof(DOUBLE));
  r->fk=(DOUBLE*)calloc(n+1,sizeof(DOUBLE));
  r->mik0=(DOUBLE*)calloc(n+1,sizeof(DOUBLE));
  r->mik1=(DOUBLE*)calloc(n+1,sizeof(DOUBLE));
  r->mik2=(DOUBLE*)calloc(n+1,sizeof(DOUBLE));
  r->mik3=(DOUBLE*)calloc(n+1,sizeof(DOUBLE));
  r->ykxj=(DOUBLE*)calloc(n+1,sizeof(DOUBLE));
  r->yk0=(DOUBLE*)calloc(n+1,sizeof(DOUBLE));
}
// generation of rectangularly distributed random numbers 0<=arr[0],...arr[15]<1 of DOUBLE type (53 bits)
void _randcalc(struct random *s){
  int i,j,k,l,m;
  DOUBLE a,a31,a22,*arr;
  a31=65536.*32768;a22=65536*64;// 2^31 and 2^22 (22+31=53 bits in DOUBLE type)
  arr=s->arr;m=(int)(arr[(int)(16*arr[0])]+.5);// m=0 or m=1
  for(i=j=k=l=0;i<16;i++){// generation
    j+=i+10+m;if(j>=26)j-=26;k+=i+j;if(k>=40)k-=40;l+=i+j+k+m;if(l>=80)l-=80;// "random" indices
    a=s->c1+arr[i]*arr[(int)(k*.4)]*s->c0+arr[(int)(j*.6)]*s->c2+
      arr[(int)(16*arr[i])]*s->c4+.071*(i*3+j*2+k)*arr[(int)(l*.2)]*s->c3;// 0.23<a<4
    a*=a31/4;arr[i]=1./a22*(int)((a-(int)(a))*a22);// high 22 bits
    a=s->c1+arr[i]*arr[(int)(k*.4)]*s->c0+arr[(int)(j*.6)]*s->c2+
      arr[(int)(16*arr[i])]*s->c4+.071*(i*3+j*2+k)*arr[(int)(l*.2)]*s->c3;
    a*=a22/4;arr[i]+=1./a31/a22*(int)((a-(int)(a))*a31);// plus lower 31 bits
} }
// initialization of generator (-2^31<=init<2^31)
void _randinit(int init,struct random *s){
  int i;
  s->c0=0.9242163423728435;
  s->c1=0.2301457072317184;
  s->c2=0.9840354236571643;
  s->c3=0.1052134963331603;
  s->c4=0.8213468021468023;
  for(i=0;i<16;i++){
    s->arr[i]=(i+90.4)/log(2.97+fabs(init*(init+13.4*i)+19.1*i));
    s->arr[i]-=(int)(s->arr[i]);
  }_randcalc(s);s->index=16;s->Xi=s->arr[0];s->xi=s->arr[1];s->aa=s->arr[2];
}
// rectangular distribution (53 bits, 50 tacts)
DOUBLE _randrect(struct random *s){
  if(s->index==16){s->index=0;_randcalc(s);}
  s->Xi+=s->arr[s->index];s->Xi-=(int)(s->Xi);s->index+=1;
  return(s->Xi);
}
// quick rectangular distribution (53 bits, 25 tacts)
DOUBLE _randrect_lowquality(struct random *s){
  DOUBLE a;
  a=1048576.*s->c0*s->xi+524288.*s->c2*s->Xi+524288.*s->c4*s->aa;a-=(int)(a);// formula
  s->aa=s->xi;s->xi+=a+1./65536./32768.*s->aa;s->xi-=(int)(s->xi);// random shift
  if(a<1./64)_randrect(s);// updating Xi every 64 calls
  return(s->xi);
}
// normal distribution (exact formula, 170 tacts)
DOUBLE _randnorm(struct random *s){
  DOUBLE a;
  a=sqrt(-2*log(1-(_randrect(s))));
  a*=cos(acos(-1.)*(_randrect(s)));
  return(a);
}
// searching a zero of f(X) (or 1/f(X)), where X=X[0],X[1],X[2],X[3]... and
//   x=X[_l], at segment [x1;x2] (dividing a segment by two)
DOUBLE _zero(DOUBLE x1,DOUBLE x2,int _l,DOUBLE *X,DOUBLE (f)(DOUBLE *X)){
  DOUBLE x,x3,y1,y2,y3;
  X[_l]=x1;y1=f(X);if(!y1)return(x1);y1/=fabs(y1);
  X[_l]=x2;y2=f(X);if(!y2)return(x2);y2/=fabs(y2);x3=x1;
  if(fabs(y1)!=1||fabs(y2)!=1){Error("\nFATAL ERROR in zero: incorrect f(x1) or f(x2)");getchar();}
  if(y1*y2>0){Error("\nFATAL ERROR in zero: search of the root of f(x) at [x1;x2] is impossible");getchar();}
  do{
    x=x3;x3=(x1+x2)*.5;X[_l]=x3;y3=f(X);
    if(!(y3&&1/y3))return(x3);else y3/=fabs(y3);
    if(fabs(y3)!=1){Error("\nFATAL ERROR in zero: incorrect f(x)");getchar();}
    if(y1*y3<0){x2=x3;y2=y3;}if(y2*y3<0){x1=x3;y1=y3;}
  }while(x!=x3);return(x);
}
// accumulating average value and error at n-th iteration
// n=0,1,2,3,...
void _accumulate(int n,DOUBLE a,DOUBLE *av,DOUBLE *er){
  if(n)*er=(*er**er*(n-1)*(n+1)*(n+1)+(a-*av)*(a-*av))/((n+1.)*(n+1.)*(n+1.));
  *av=(n**av+a)/(n+1.);if(n)*er=sqrt(*er+(a-*av)*(a-*av)/(n*(n+1.)));
}
// // // // // // // // // // // // // // // // // // // // // // // // // //
//                block for summation, integration and limits
//
// quick integration of f(x) on dx over [a;b] in m-th order;
// array f(a+(.5-m)*dx),f(a+(1.5-m)*dx),...f(b+(m-.5)*dx)
//   is given by fxj[0],fxj[1],...fxj[jmax+2*m-1], with b-a=dx*jmax and jmax>0
DOUBLE _integ_arr(DOUBLE a,DOUBLE b,int jmax,int m,DOUBLE *fxj,struct exint e){
  DOUBLE I0,A;
  int j;
  if(a==b)return(0);
  if(m<0||m>7||jmax<1){Error("\nFATAL ERROR in integ_arr: one must be 0<=m<=7 and jmax>=1 (m=%d,jmax=%d)",m,jmax);getchar();}
  if(jmax>=2*m){
    for(I0=0,j=m;j<jmax-m;j++)I0+=fxj[j+m];
    for(j=-m;j<m;j++)I0+=e.Ilm[(j+m)*e.M+m]*fxj[j+m];
    for(j=-m;j<m;j++)I0+=(1-e.Ilm[(j+m)*e.M+m])*fxj[j+jmax+m];
  }else{
    for(I0=0,j=-m;j<jmax+m;j++){
      if(j>=m)A=1;else A=e.Ilm[(j+m)*e.M+m];
      if(j>=jmax-m)A-=e.Ilm[(j+m-jmax)*e.M+m];
      I0+=A*fxj[j+m];
  } }return(I0*(b-a)/jmax);
}
// quick integration of f(X[0],X[1],X[2],X[3],...) on dX[_l] over [a;b] in m-th order
DOUBLE _integ_quick(DOUBLE a,DOUBLE b,int jmax,int m,int _l,DOUBLE *X,
    struct exint e,DOUBLE (f)(DOUBLE *X)){
  DOUBLE I0,dx,A;
  int j;
  if(a==b)return(0);
  if(m<0||m>7||jmax<1){Error("\nFATAL ERROR in integ_quick: one must be jmax>=1 and 0<=m<=7 (jmax=%d, m=%d)",jmax,m);getchar();}
  I0=0;dx=(b-a)/jmax;if(jmax>=2*m){
    for(j=m;j<jmax-m;j++){X[_l]=a+(j+.5)*dx;I0+=f(X);}
    for(j=-m;j<m;j++){X[_l]=a+(j+.5)*dx;I0+=e.Ilm[(j+m)*e.M+m]*f(X);}
    for(j=-m;j<m;j++){X[_l]=a+(j+jmax+.5)*dx;I0+=(1-e.Ilm[(j+m)*e.M+m])*f(X);}
  }else{
    for(j=-m;j<jmax+m;j++){
      if(j>=m)A=1;else A=e.Ilm[(j+m)*e.M+m];
      if(j>=jmax-m)A-=e.Ilm[(j+m-jmax)*e.M+m];
      X[_l]=a+(j+.5)*dx;I0+=A*f(X);
  } }return(I0*dx);
}
// summation of aj(X[0],X[1],X[2],X[3],...) with X[_l]=j from j=1 to j=+infty:
// direct joint with the tail at j=jm in m-th order
// \int_{jm+.5}^{\infty}a_jdj must be added; one must be jm>=m and 0<=m<=7
DOUBLE _sumaj(int jm,int m,int _l,DOUBLE *X,struct exint e,DOUBLE (aj)(DOUBLE *X)){
  DOUBLE s0;
  int j;
  if(m<0||m>=8||jm<m){Error("\nFATAL ERROR in sumaj: one must be 0<=m<=7 and jm>=m (m=%d, jm=%d)",m,jm);getchar();}
  s0=0;for(j=1;j<=jm-m;j++){X[_l]=j;s0+=aj(X);}
  for(j=jm-m+1;j<=jm+m;j++){X[_l]=j;s0+=(1-e.Ilm[(j+m-jm-1)*e.M+m])*aj(X);}
  return(s0);
}
// direct taking of a limit: f(x0)=f(x0-dx)=fxj[0], f(x0+dx)=f(x0-2*dx)=fxj[1],
//   ..., f(x+(j0-1)*dx)=f(x-j0*dx)=fxj[j0-1], j0=3,4,5,...16,
// limit_{x-->(x0-dx/2)}f(x)=f(x0-dx/2)
DOUBLE _limit_dir(int j0,struct limit l,struct exint e){
  int j,m,n,step;
  DOUBLE y;
  m=(j0-1)/2;if(m>7)m=7;n=j0-1;if(n>10)n=10;step=1;// orders
  for(j=0;j<j0;j++)// simmetrizing and errors
    l.a[j+j0]=l.a[j0-1-j]=l.fxj[j]+l.erj[j]*l.raj[j];
  y=_fxjkx(-.5,1.,0,m,n,j0-1,step,l.a+1,e);// Taylor series for f(x0-dx/2)
  return(y);
}
// advanced taking a limit: first point fxj[0] is skipped, j0-1 are given:
// f(x0+dx)=fxj[1],f(x0+2*dx)=fxj[2],...f(x0+(j0-1)*dx)=fxj[j0-1];  j0=3,4,5,...16
// limit_{x-->(x0-dx/2)}f(x)=f(x0-dx/2);
// recovered point f(x0) is loaded to fxj[0]
DOUBLE _limit_skip(int j0,struct limit l,struct exint e){
  int j,k,m,n,step;
  DOUBLE y,dy,a1,B,C;
  m=(j0-1)/2;if(m>7)m=7;k=2*m;if(k>14)k=14;n=j0-1;if(n>10)n=10;step=1;// orders
  // simmetrizing (first point is skipped) and errors
    for(j=1;j<j0;j++)l.a[j+j0]=l.a[j0-1-j]=l.fxj[j]+l.erj[j]*l.raj[j];
  // calculating dy (its characteristic value for the function)
    for(dy=0,j=1;j<j0;j++)dy+=fabs(l.fxj[j]-l.fxj[1]);
    if(!dy){l.fxj[0]=l.fxj[1];return(l.fxj[1]);}// at f(x)=const
  // calculating the quantity B
    // calculating the quantity C
      a1=l.fxj[1];l.a[j0]=l.a[j0-1]=a1;// calculate without increment
      C=_fxjkx(-.5,1.,k,m,1,j0-1,step,l.a+1,e);
      a1=l.fxj[1]-dy;l.a[j0]=l.a[j0-1]=a1;// calculate with increment
      C-=_fxjkx(-.5,1.,k,m,1,j0-1,step,l.a+1,e);// and subtract
      C/=dy;// and divide by the increment
    a1=l.fxj[1];l.a[j0]=l.a[j0-1]=a1;// set the unknown value of the function to be equal to "flat bottom"
    B=-_fxjkx(-.5,1.,k,m,1,j0-1,step,l.a+1,e);// and calculate Taylor series with it
    a1=C;B+=a1*l.fxj[1];// then, add the "flat bottom" explicitly
  // solving an equation: fxj[0]=B/C
    a1=B/C;// solve and find the skipped value of the function
    l.fxj[0]=l.a[j0]=l.a[j0-1]=a1;// и заносим ее в массив
  // final Taylor series with the RECOVERED function value
    y=_fxjkx(-.5,1.,0,m,n,j0-1,step,l.a+1,e);
  return(y);
}
// calculating average value and error in taking limit, j0=3,4,5,...16
// N is statistics (N=2,3,4,...); at N=0 errors erj[j] not use
// the calculated error is ONLY due to an error in input array fxj[]
// (i.e., an error of interpolation do not be taken into account)
DOUBLE _limit_err(int j0,int N,DOUBLE *err,struct limit l,struct exint e,
    struct random *r,DOUBLE (lim)(int j0,struct limit l,struct exint e)){
  int n,j;
  DOUBLE a,av,er;
  if(!N){
    for(j=0;j<j0;j++)l.raj[j]=0;av=lim(j0,l,e);er=0;
  }else for(av=er=n=0;n<N;n++){
    for(j=0;j<j0;j++)l.raj[j]=_randnorm(r);// generating random numbers
    a=lim(j0,l,e);// n-th taking a limit
    _accumulate(n,a,&av,&er);// average value and error
  }*err=er;return(av);
}
// memory for arrays used for taking limit
void _limit_mem(struct limit *l){
  l->fxj=(DOUBLE*)calloc(80,sizeof(DOUBLE));
  l->erj=l->fxj+16;l->raj=l->erj+16;l->a=l->raj+16;
}

/************************** Construct Kurbakov's teacher **********************/
// // // // // // // // // // // // // // // // // // // // // // // // // //
//                             functions for examples
//
// simplest examples of a functions
DOUBLE _f1(DOUBLE x){
  return(cosh(x)+sin(x));
}
DOUBLE _f2(DOUBLE x){
  return(exp(-x*log(2.5))*log(2.5));
}
DOUBLE _f3(DOUBLE x){
  return(exp((x-1)*(x-1)/2)+sin(x-1));
}
DOUBLE _f3_(DOUBLE *X){
  return(_f3(X[0]));
}
DOUBLE _f4(DOUBLE x){
  return(sin(x*cos(x))/(exp(sin(x))-1));
}
// another examples of functions
DOUBLE _fvp(DOUBLE *x){// 1/(2-x-x^2)
  return(1/(2-*x-*x**x));
}
DOUBLE _aj(DOUBLE *j){// a_j=2/[(2j+1)(2j-1)]
  return(2/((*j*2+1)*(*j*2-1)));
}
DOUBLE _Saj(DOUBLE *X){// S=1-2+3-4+5-6+7-8+9-...
  int j,Jm;
  DOUBLE eps,b,aj,S;
  eps=X[0];b=10;Jm=(int)(b/eps);// j*exp(-eps^2*j^2/2)=0 at j>Jm
  for(S=0,j=1;j<=Jm;j++){
     if(j==2*(j/2))aj=-j;else aj=j;// sequence
     S+=aj*exp(-eps*eps*j*j/2.);// slow cutoff, eps-->0
  }return(S);
}
DOUBLE _fFt(DOUBLE *X){
  DOUBLE eps,u,x,a,f;
  eps=X[0];u=X[1];x=X[2];a=X[3];
  f=a*2/pi/(a*a+x*x)*cos(u*x);// integrand
  f*=exp(-eps*eps/2*x*x/a/a);// cutoff, eps-->0
  return(f);
}
// example of integrating a multi argument function
DOUBLE _fXu(DOUBLE *X){
  DOUBLE x,y,z,t,u,v,w,o,f;
  x=X[0];y=X[1];z=X[2];t=X[3];u=X[4];v=X[5];w=X[6];o=X[7];
  f=exp(-u*u/2-v*v/2-w*w/2-o*o/2)*cos(u*x)*cos(v*y)*cos(w*z)*cos(o*t)*4/pi/pi;
  return(f);
}
DOUBLE _fXv(DOUBLE *X){
  int jmax;
  DOUBLE f,b;
  b=X[8];jmax=(int)(X[9]);f=_integ_quick(0,b,jmax,7,4,X,M8,_fXu);
  return(f);
}
DOUBLE _fXw(DOUBLE *X){
  int jmax;
  DOUBLE f,b;
  b=X[8];jmax=(int)(X[10]);f=_integ_quick(0,b,jmax,7,5,X,M8,_fXv);
  return(f);
}
DOUBLE _fXo(DOUBLE *X){
  int jmax;
  DOUBLE f,b;
  b=X[8];jmax=(int)(X[11]);f=_integ_quick(0,b,jmax,7,6,X,M8,_fXw);
  return(f);
}
DOUBLE _fX(DOUBLE *X){
  int jmax;
  DOUBLE f,b;
  b=X[8];jmax=(int)(X[12]);f=_integ_quick(0,b,jmax,7,7,X,M8,_fXo);
  return(f);
}
// examples of functions for Runge intergation
DOUBLE _func1(int _l,DOUBLE *X,DOUBLE *yl){
  DOUBLE x;
  x=X[_l];
  return((yl[1]+exp(x*(1+yl[2]*yl[2]+yl[3]*yl[3])/2))/2);
}
DOUBLE _func2(int _l,DOUBLE *X,DOUBLE *yl){
  DOUBLE x;
  x=X[_l];
  return((1-sin(x)*sin(x))/yl[3]);
}
DOUBLE _func3(int _l,DOUBLE *X,DOUBLE *yl){
  DOUBLE x;
  x=X[_l];
  return(-yl[2]*sin(x)*sin(x)-sin(2*x)*yl[3]/2);
}
void _funcs(int _l,DOUBLE *X,DOUBLE *yl,DOUBLE *fk){// array of functions
  fk[1]=_func1(_l,X,yl);fk[2]=_func2(_l,X,yl);fk[3]=_func3(_l,X,yl);
}
// // // // // // // // // // // // // // // // // // // // // // // // // //
//                             I. Basic theory.
//
// I. Basic theory. 1. Differentiation.
void main_I_1(void){
  int j,j0,m,n,step;
  DOUBLE x,dx,A,A1,A2;
  Message("I. Basic theory. 1. Differentiation\n  ");
  _extrmatr(8,&M8);// initializing
  j0=15;m=j0/2;n=3;step=1;x=0;// parameters
  dx=0.14;for(j=0;j<=j0*2;j++)M8.yxj[j]=_f1(x+(j-j0)*dx);// calculating array
  A1=_fxjnmx(dx,n,m,j0,step,M8.yxj,M8);// differentiation
  dx=0.10;for(j=0;j<=j0*2;j++)M8.yxj[j]=_f1(x+(j-j0)*dx);//   the same with
  A2=_fxjnmx(dx,n,m,j0,step,M8.yxj,M8);//                     smaller step
  A=-1;Message("f'''(%g)=%g, err=%g=%g",x,A1,fabs(A1/A-1),fabs(A2/A-1));
}
// I. Basic theory. 2. Integration.
void main_I_2(void){
  int j,m,jmax,Jmax;
  DOUBLE dx,a,b,A,A1,A2,*yj;
  Message("I. Basic theory. 2. Integration.\n  ");
  _extrmatr(8,&M8);// initializing
  Jmax=1024;yj=(DOUBLE*)calloc(Jmax,sizeof(DOUBLE));// memory
  a=-2;b=2;m=7;// segment of intergation [a;b] and order
  jmax=16;dx=(b-a)/jmax;for(j=0;j<jmax+2*m;j++)yj[j]=_f1(a+(j-m+.5)*dx);// array
  A1=_integ_arr(a,b,jmax,m,yj,M8);// intergation
  jmax=18;dx=(b-a)/jmax;for(j=0;j<jmax+2*m;j++)yj[j]=_f1(a+(j-m+.5)*dx);//  the same with
  A2=_integ_arr(a,b,jmax,m,yj,M8);//                                        smaller step
  A=exp(2.)-exp(-2.);Message("I=%g, err=%g=%g",A1,fabs(A1/A-1),fabs(A2/A-1));
}
// I. Basic theory. 3. Interpolation.
void main_I_3(void){
  int j,m,jmax,nmax,j0,jx,k,step,Jmax;
  DOUBLE x,x0,dx,a,A,A1,A2,*yj;
  Message("I. Basic theory. 3. Interpolation\n  ");
  _extrmatr(8,&M8);// initializing
  Jmax=1024;yj=(DOUBLE*)calloc(Jmax,sizeof(DOUBLE));// memory
  a=-0.549;jmax=96;nmax=10;j0=15;m=j0/2;k=2;step=2;x=0;// parameters
  dx=0.12/step;for(j=0;j<=j0*step*2+jmax;j++)yj[j]=_f1(a+(j-j0*step)*dx);// calculating array
  jx=(int)((x-a)/dx+.5);x0=a+dx*jx;// calculating the nearest point
  A1=_fxjkx(x-x0,dx,k,m,nmax,j0,step,yj+jx,M8);// interpolation
  dx=0.10/step;for(j=0;j<=j0*step*2+jmax;j++)yj[j]=_f1(a+(j-j0*step)*dx);//  the same
  jx=(int)((x-a)/dx+.5);x0=a+dx*jx;//                                        with the 
  A2=_fxjkx(x-x0,dx,k,m,nmax,j0,step,yj+jx,M8);//                            smaller step
  A=1;Message("f''(%g)=%g, err=%g=%g",x,A1,fabs(A1/A-1),fabs(A2/A-1));
}
// I. Basic theory. 4. Runge integration.
void main_I_4(void){
  int n;
  DOUBLE x,err,*X;
  struct Runge ru;
  Message("I. Basic theory. 4. Runge integration.\n  ");
  n=3;x=1;// number of functions and segment of integration [0;x]
  _rungematr(n,&ru);// initializing
  ru.arr=(DOUBLE*)calloc(6144,sizeof(DOUBLE));// memory for saving functions
  X=(DOUBLE*)calloc(16,sizeof(DOUBLE));// array of variables
  ru.yk0[1]=1;ru.yk0[2]=0;ru.yk0[3]=1;ru._l=0;// initial conditions and number of a variable
  ru.a=0;ru.b=x;ru.nn=2000;ru.s=ru.n;// segment of integration, number of points, saving functions
  _runge(X,&ru,_funcs);// Runge integration
  err=sqrt(// error
    (ru.ykxj[1]/exp(x)-1)*(ru.ykxj[1]/exp(x)-1)//     y1=exp(x)
    +(ru.ykxj[2]/sin(x)-1)*(ru.ykxj[2]/sin(x)-1)//    y2=sin(x)
    +(ru.ykxj[3]/cos(x)-1)*(ru.ykxj[3]/cos(x)-1));//  y3=cos(x)
  Message("exp(%g)=%g, sin(%g)=%g, cos(%g)=%g;  err=%g",x,ru.ykxj[1],x,ru.ykxj[2],x,ru.ykxj[3],err);
}
// I. Basic theory. 5. Generator of random numbers.
void main_I_5(void){
  int init;
  struct random ra;
  Message("I. Basic theory. 5. Generator of random numbers.\n  ");
  init=2;// number of random path
  _randinit(init,&ra);// initialization
  Message("rect=%lf, ",_randrect(&ra));// rectangular distribution
  Message("norm=%lf",_randnorm(&ra));// normal distribution
}
// // // // // // // // // // // // // // // // // // // // // // // // // //
//                        II. Ultraquick calculations.
//
// II. Ultraquick calculations. 1. First derivative.
void main_II_1(void){
  int j,j0,jmax,Jmax;
  DOUBLE a,x,dx,a1,A,*yj;
  Message("II. Ultraquick calculations. 1. First derivative.\n  ");
  Jmax=1024;yj=(DOUBLE*)calloc(Jmax,sizeof(DOUBLE));// memory
  a=-1;dx=0.1;a1=1/dx;j0=7;jmax=100;// parameters
  for(j=0;j<=2*j0+jmax;j++)yj[j]=_f3(a+(j-j0)*dx);// array
  j=j0+20;x=a+(j-j0)*dx;// a choice of a point (x=1)
  // 2-th order
  A=a1*(.5*yj[j+1]-.5*yj[j-1]);Message("f'(%g)*1E4=%lf",x,A*1E4);
  // 4-th order
  A=a1*(27/48.*(yj[j+1]-yj[j-1])-1/48.*(yj[j+3]-yj[j-3]));Message("=%lf",A*1E4);
  // 6-th order
  A=a1*(2250/3840.*(yj[j+1]-yj[j-1])-125/3840.*(yj[j+3]-yj[j-3])+
    9/3840.*(yj[j+5]-yj[j-5]));Message("=%lf",A*1E4);
  // 8-th order
  A=a1*(128625/420./512.*(yj[j+1]-yj[j-1])-
    8575/420./512.*(yj[j+3]-yj[j-3])+1029/420./512.*(yj[j+5]-yj[j-5])-
    75/420./512.*(yj[j+7]-yj[j-7]));Message("=%lf",A*1E4);
  // exact
  Message("=1E4");
}
// II. Ultraquick calculations. 2. Second derivative.
void main_II_2(void){
  int j,j0,jmax,Jmax;
  DOUBLE a,x,dx,a1,A,*yj;
  Message("II. Ultraquick calculations. 2. Second derivative.\n  ");
  Jmax=1024;yj=(DOUBLE*)calloc(Jmax,sizeof(DOUBLE));// memory
  a=-1;dx=0.04;a1=1/dx/dx;j0=8;jmax=100;// parameters
  for(j=0;j<=2*j0+jmax;j++)yj[j]=_f3(a+(j-j0)*dx);// array
  j=j0+50;x=a+(j-j0)*dx;// a choice of a point (x=1)
  // 2-nd order
  A=a1*(-.5*yj[j]+.25*(yj[j+2]+yj[j-2]));Message("f''(%g)*1E4=%lf",x,A*1E4);
  // 4-th order
  A=a1*(-30/48.*yj[j]+16/48.*(yj[j+2]+yj[j-2])-
    1/48.*(yj[j+4]+yj[j-4]));Message("=%lf",A*1E4);
  // 6-th order
  A=a1*(-490/720.*yj[j]+270/720.*(yj[j+2]+yj[j-2])-
    27/720.*(yj[j+4]+yj[j-4])+2/720.*(yj[j+6]+yj[j-6]));Message("=%lf",A*1E4);
  // 8-th order
  A=a1*(-14350/120./168.*yj[j]+8064/120./168.*(yj[j+2]+yj[j-2])-
    1008/120./168.*(yj[j+4]+yj[j-4])+128/120./168.*(yj[j+6]+yj[j-6])-
    9/120./168.*(yj[j+8]+yj[j-8]));Message("=%lf",A*1E4);
  // exact
  Message("=1E4");
}
// II. Ultraquick calculations. 3. Integral.
void main_II_3(void){
  int j,jmax,m,Jmax;
  DOUBLE a,b,dx,A,*yj;
  Message("II. Ultraquick calculations. 3. Integral.\n  ");
  Jmax=1024;yj=(DOUBLE*)calloc(Jmax,sizeof(DOUBLE));// memory
  a=0;b=1;jmax=4;dx=(b-a)/jmax;// diapazone and step
  // 2-nd order
  for(A=j=0;j<jmax;j++)A+=_f2(a+(j+.5)*dx)*dx;// array
  Message("I*1E4=%lf",A*1E4);
  // 4-th order
  m=1;for(j=0;j<jmax+2*m;j++)yj[j]=_f2(a+(j-m+.5)*dx);// array
  for(A=0,j=m;j<jmax+m;j++)A+=yj[j]*dx;// sum inside [a;b]
  A+=1/24.*(yj[0]-yj[1]-yj[jmax+0]+yj[jmax+1])*dx;// difference scheme
  Message("=%lf",A*1E4);
  // 6-th order
  m=2;for(j=0;j<jmax+2*m;j++)yj[j]=_f2(a+(j-m+.5)*dx);
  for(A=0,j=m;j<jmax+m;j++)A+=yj[j]*dx;
  A+=1/5760.*(-17*(yj[0]-yj[3]-yj[jmax+0]+yj[jmax+3])
             +291*(yj[1]-yj[2]-yj[jmax+1]+yj[jmax+2]))*dx;
  Message("=%lf",A*1E4);
  // 8-th order
  m=3;for(j=0;j<jmax+2*m;j++)yj[j]=_f2(a+(j-m+.5)*dx);
  for(A=0,j=m;j<jmax+m;j++)A+=yj[j]*dx;
  A+=1/967680.*(  367*(yj[0]-yj[5]-yj[jmax+0]+yj[jmax+5])
                -4691*(yj[1]-yj[4]-yj[jmax+1]+yj[jmax+4])
               +52558*(yj[2]-yj[3]-yj[jmax+2]+yj[jmax+3]))*dx;
  Message("=%lf",A*1E4);
  // exact
  Message("=0.6E4");
}
// II. Ultraquick calculations. 4. High-order spline.
void main_II_4(void){
  struct spline sp;
  int _l,m,Jm,Jmax,jmax,stma;
  DOUBLE x,a,b,A,A1,A2,A3,*X,*yj;
  Message("II. Ultraquick calculations. 4. High-order spline.\n  ");
  X=(DOUBLE*)calloc(16,sizeof(DOUBLE));// array of variables
  Jmax=1024;yj=(DOUBLE*)calloc(Jmax,sizeof(DOUBLE));// memory for spline
  a=-0.526;b=1.514;m=7;stma=1;_l=0;Jm=Jmax;x=1;sp.y=yj;// parameters
  jmax=(int)((b-a)*400);_splineprep(a,b,jmax,m,stma,Jm,&sp);// making a spline structure
  _splinecalc(_l,X,sp,_f3_);// calculating a spline array
  A1=_spline5(x,sp);// interpolation with error of O(\Delta x^6)
  jmax=(int)((b-a)*125);_splineprep(a,b,jmax,m,stma,Jm,&sp);//   the same for
  _splinecalc(_l,X,sp,_f3_);A2=_spline7(x,sp);//                 O(\Delta x^8)
  jmax=(int)((b-a)*50);_splineprep(a,b,jmax,m,stma,Jm,&sp);//  the same for
  _splinecalc(_l,X,sp,_f3_);A3=_spline9(x,sp);//               O(\Delta x^10)
  A=1;Message("f(%g)=%g=%g=%g, err=%g=%g=%g",x,A1,A2,A3,fabs(A1/A-1),fabs(A2/A-1),fabs(A3/A-1));
}
// II. Ultraquick calculations. 5. First derivative of high-order spline.
void main_II_5(void){
  struct spline sp;
  int _l,m,Jm,Jmax,jmax,stma;
  DOUBLE x,a,b,A,A1,A2,A3,*X,*yj;
  Message("II. Ultraquick calculations. 5. First derivative of high-order spline.\n  ");
  X=(DOUBLE*)calloc(16,sizeof(DOUBLE));// array of variables
  Jmax=1024;yj=(DOUBLE*)calloc(Jmax,sizeof(DOUBLE));// memory for spline
  a=-0.516;b=1.524;m=7;stma=1;_l=0;Jm=Jmax;x=1;sp.y=yj;// parameters
  jmax=(int)((b-a)*400);_splineprep(a,b,jmax,m,stma,Jm,&sp);// making a spline structure
  _splinecalc(_l,X,sp,_f3_);// calculating a spline array
  A1=_deriv5(x,sp);// interpolation with error of O(\Delta x^6)
  jmax=(int)((b-a)*125);_splineprep(a,b,jmax,m,stma,Jm,&sp);//   the same for
  _splinecalc(_l,X,sp,_f3_);A2=_deriv7(x,sp);//                  O(\Delta x^8)
  jmax=(int)((b-a)*50);_splineprep(a,b,jmax,m,stma,Jm,&sp);//  the same for
  _splinecalc(_l,X,sp,_f3_);A3=_deriv9(x,sp);//                O(\Delta x^10)
  A=1;Message("f'(%g)=%g=%g=%g, err=%g=%g=%g",x,A1,A2,A3,fabs(A1/A-1),fabs(A2/A-1),fabs(A3/A-1));
}
// // // // // // // // // // // // // // // // // // // // // // // // // //
//                         III. Treatment of splines.
//
// III. Treatment of splines. 1. Calculation of step and making a spline.
void main_III_1(void){
  struct spline sp;
  int _l,m,Jm,Jmax,jmax,stma;
  DOUBLE a,b,err,*X,*yj,*ydj;
  Message("III. Treatment of splines. 1. Calculation of step and making a spline.\n  ");
  _extrmatr(8,&M8);// initializing
  X=(DOUBLE*)calloc(16,sizeof(DOUBLE));// array of variables
  Jmax=1024;yj=(DOUBLE*)calloc(2*Jmax,sizeof(DOUBLE));// memories for a spline and a dual spline
  a=-1;b=3;m=4;jmax=3;Jm=Jmax;_l=0;stma=1;ydj=yj+Jm;// parameters
  err=_calcstep(a,b,m,Jm,&jmax,_l,X,yj,ydj,M8,_f3_);// calculating a step
  _splineprep(a,b,jmax,m,stma,Jm,&sp);// making a spline structure
  sp.y=yj;_splinecalc(_l,X,sp,_f3_);// calculating a spline array
  Message("jmax=%d, dx=%g, err=%g",sp.jmax,sp.dx,err);
}
// III. Treatment of splines. 2. Doubling a spline.
void main_III_2(void){
  struct spline sp;
  int _l,m,Jm,Jmax,jmax,stma;
  DOUBLE x,a,b,A,A1,A2,A3,err,*X,*yj;
  Message("III. Treatment of splines. 2. Doubling a spline.\n  ");
  _extrmatr(8,&M8);// initializing
  X=(DOUBLE*)calloc(16,sizeof(DOUBLE));// array of variables
  Jmax=1024;yj=(DOUBLE*)calloc(2*Jmax,sizeof(DOUBLE));// memory for a spline and a dual spline
  a=-1;b=3;x=2.8;m=7;jmax=4;Jm=Jmax;_l=0;stma=1;sp.y=yj;sp.yd=yj+Jm;// parameters
  err=_calcstep(a,b,m,Jm,&jmax,_l,X,sp.y,sp.yd,M8,_f3_);// calculating a step
  _splineprep(a,b,jmax,m,stma,Jm,&sp);_splinecalc(_l,X,sp,_f3_);// making a spline
  A1=_spline9(x,sp);// interpolation with large dx
  _doublingspline(&sp,M8);A2=_spline9(x,sp);// doubling a spline and interpolation with smaller dx
  _doublingspline(&sp,M8);A3=_spline9(x,sp);// once more
  A=_f3_(&x);Message("  f(%g)=%g, err=%g-->%g-->%g, dx=%lf",x,A3,fabs(A1/A-1),fabs(A2/A-1),fabs(A3/A-1),sp.dx);
}
// III. Treatment of splines. 3. Dual transformation of a spline.
void main_III_3(void){
  struct spline sp,sp_;
  int _l,m,Jm,Jmax,jmax,stma;
  DOUBLE x,a,b,A,A1,err,*X,*yj;
  Message("III. Treatment of splines. 3. Dual transformation of a spline.\n  ");
  _extrmatr(8,&M8);// initializing
  X=(DOUBLE*)calloc(16,sizeof(DOUBLE));// array of variables
  Jmax=1024;yj=(DOUBLE*)calloc(2*Jmax,sizeof(DOUBLE));// memory for a spline and a dual spline
  a=-1;b=3;x=1.8;m=7;jmax=4;Jm=Jmax;_l=0;stma=2;sp.y=yj;sp.yd=yj+Jm;// parameters
  err=_calcstep(a,b,m,Jm,&jmax,_l,X,sp.y,sp.yd,M8,_f3_);// calculating a step
  _splineprep(a,b,jmax,m,stma,Jm,&sp);_splinecalc(_l,X,sp,_f3_);// making a spline
  _doublingspline(&sp,M8);_doublingspline(&sp,M8);// doubling a spline 2 times
  _dualspline(sp,M8);// calculating a dual spline
  _convertdual(sp,&sp_);// making dual spline as a new spline sp_
  A=_f3_(&x);A1=_spline9(x,sp_);// exact value and interpolation
  Message("f(%g)=%g, err=%g",x,A,fabs(A1/A-1));
}
// III. Treatment of splines. 4. Differentiation of a spline.
void main_III_4(void){
  struct spline sp,sp_;
  int _l,m,k,Jm,Jmax,jmax,stma;
  DOUBLE x,a,b,A,A1,err,*X,*yj;
  Message("III. Treatment of splines. 4. Differentiation of a spline.\n  ");
  _extrmatr(8,&M8);// initializing
  X=(DOUBLE*)calloc(16,sizeof(DOUBLE));// array of variables
  Jmax=1024;yj=(DOUBLE*)calloc(2*Jmax,sizeof(DOUBLE));// memory for a spline and its k-th derivative
  a=-1;b=3;x=1;m=7;k=3;jmax=4;Jm=Jmax;_l=0;stma=2;sp.y=yj;sp.yd=yj+Jm;// parameters
  err=_calcstep(a,b,m,Jm,&jmax,_l,X,sp.y,sp.yd,M8,_f3_);// calculating a step
  _splineprep(a,b,jmax,m,stma,Jm,&sp);_splinecalc(_l,X,sp,_f3_);// making a spline
  _doublingspline(&sp,M8);_doublingspline(&sp,M8);// doubling a spline 2 times
  sp.y1=yj+Jm;sp.step=2*2;_splinederivative(k,sp,M8);// calculating k-th derivative of a spline
  _convertderivative(sp,&sp_);// making k-th derivative as a new spline
  A1=_spline9(x,sp_);A=-1;// interpolation of k-th derivative and its exact value
  Message("f'''(%g)=%g, err=%g",x,A,fabs(A1/A-1));
}
// III. Treatment of splines. 5. Calculation of errors in a spline data.
void main_III_5(void){
  struct random ra;
  struct spline sp;
  int _l,m,Jm,Jmax,jmax,j,Jer,stma,init;
  DOUBLE a,b,err,sig,*X,*yj;
  Message("III. Treatment of splines. 5. Calculation of errors in a spline data.\n  ");
  _extrmatr(8,&M8);init=2;_randinit(init,&ra);// initialization
  X=(DOUBLE*)calloc(16,sizeof(DOUBLE));// array of variables
  Jmax=1024;yj=(DOUBLE*)calloc(2*Jmax,sizeof(DOUBLE));// memory for spline
  a=-1;b=3;m=7;jmax=30;Jm=Jmax;_l=0;stma=1;sp.y=yj;sp.yd=yj+Jm;// parameters
  err=_calcstep(a,b,m,Jm,&jmax,_l,X,sp.y,sp.yd,M8,_f3_);// calculating a step
  _splineprep(a,b,jmax,m,stma,Jm,&sp);_splinecalc(_l,X,sp,_f3_);// making a spline
  sig=0E-16;for(j=0;j<=sp.jm;j++)sp.y[j]+=sig*_randnorm(&ra);// adding a noise to a spline
  _doublingspline(&sp,M8);// doubling a spline (for a quality of interpolation)
  sp.er=yj+Jm;Jer=2;err=_splineerror(Jer,sp,M8);// calculating errors
  Message("err=%g, er[%d]=%g,er[%d]=%g,er[%d]=%g",err,
    Jer+sp.m,sp.er[Jer+sp.m],sp.jm/2,sp.er[sp.jm/2],sp.jm-Jer-sp.m,sp.er[sp.jm-Jer-sp.m]);
}
// III. Treatment of splines. 6. Analiticheskoe prodolzhenie of a spline.
void main_III_6(void){
  struct spline sp;
  int _l,m,k,Jm,Jmax,jmax,stma;
  DOUBLE rez,imz,Ref,Imf,a,b,A1,A2,err,*X,*yj;
  Message("III. Treatment of splines. 6. Analiticheskoe prodolzhenie of a spline.\n  ");
  _extrmatr(8,&M8);// initializing
  X=(DOUBLE*)calloc(16,sizeof(DOUBLE));// array of variables
  Jmax=1024;yj=(DOUBLE*)calloc(2*Jmax,sizeof(DOUBLE));// memory for spline
  a=-2;b=3;k=0;rez=1.5;imz=.04;// diapazone, order of differentiation, rez+i*imz
  m=7;jmax=4;Jm=Jmax;_l=0;stma=2;sp.y=yj;sp.yd=yj+Jm;// other parameters
  err=_calcstep(a,b,m,Jm,&jmax,_l,X,sp.y,sp.yd,M8,_f3_);// calculating a step
  _splineprep(a,b,jmax,m,stma,Jm,&sp);_splinecalc(_l,X,sp,_f3_);// making a spline
  _doublingspline(&sp,M8);_doublingspline(&sp,M8);// doubling a spline 2 times
  sp.step=2*2;_spline_anal(rez,imz,&Ref,&Imf,k,sp,M8);
  jmax+=jmax/4;_splineprep(a,b,jmax,m,stma,Jm,&sp);_splinecalc(_l,X,sp,_f3_);
  _doublingspline(&sp,M8);_doublingspline(&sp,M8);
  sp.step=2*2+1;_spline_anal(rez,imz,&A1,&A2,k,sp,M8);// the same with smaller dx
  err=sqrt(((Ref-A1)*(Ref-A1)+(Imf-A2)*(Imf-A2))/(A1*A1+A2*A2));// "error"
  Message("f(%g+i*%g)=%g+i*%g, err=%g",rez,imz,Ref,Imf,err);
}
// // // // // // // // // // // // // // // // // // // // // // // // // //
//                  IV. Advanced theory and applications.
//
// IV. Advanced theory and applications. 1. Calculation of infinite sums.
void main_IV_1(void){
  int _l,m,k,jm,jmax,Jmax;
  DOUBLE x,dx,A,A1,der,*X,*yj;
  Message("IV. Advanced theory and applications. 1. Calculation of infinite sums.\n  ");
  _extrmatr(8,&M8);// initializing
  X=(DOUBLE*)calloc(16,sizeof(DOUBLE));// array of variables
  Jmax=1024;yj=(DOUBLE*)calloc(Jmax,sizeof(DOUBLE));// memory
  jm=20;m=7;_l=0;A1=_sumaj(jm,m,_l,X,M8,_aj);// summation of a_j
  jmax=14;dx=1./jmax;for(k=-m;k<jmax+m;k++){// \int_{jm+.5}^{\infty}a_{j(x)}(dj(x)/dx)dx
    x=(k+.5)*dx;der=(jm+.5)/(1-x)/(1-x);X[_l]=(jm+.5)/(1-x);
    yj[k+m]=_aj(X)*der;// j(x)=(jm+.5)/(1-x), j'(x)=(jm+.5)/(1-x)/(1-x)
  }A1+=_integ_arr(0,1,jmax,m,yj,M8);// plus integral of the tail
  A=1;Message("S=%g, err=%g",A1,fabs(A1/A-1));
}
// IV. Advanced theory and applications. 2. Taking limits.
void main_IV_2(void){
  struct random ra;
  struct limit li;
  int j,j0,N,init;
  DOUBLE x,dx,A,A1,A2,A3,err;
  Message("IV. Advanced theory and applications. 2. Taking limits.\n  ");
  _extrmatr(8,&M8);_limit_mem(&li);init=1;_randinit(init,&ra);// initialization
  j0=16;N=1000;for(j=0;j<j0;j++)li.erj[j]=1E-15;// random error in an initial data
  dx=0.047;for(j=0;j<j0;j++){x=(j+.5)*dx;li.fxj[j]=(_f4(x)+_f4(-x))/2;}// array
  A1=_limit_err(j0,N,&err,li,M8,&ra,_limit_dir);// direct taking a limit
  dx=0.016;for(j=1;j<j0;j++){x=(j+.5)*dx;li.fxj[j]=(_f4(x)+_f4(-x))/2;}// array
  A2=_limit_err(j0,N,&err,li,M8,&ra,_limit_skip);// advanced taking a limit: the first point is skipped
  A3=_limit_skip(j0,li,M8);// advanced taking limit without calculating error
  A=1;Message("lim=%g, err=%g, err_dir=%g, err_skip=%g",A3,fabs(A3/A-1),fabs(A1/A-1),fabs(A2/A-1));
}
// IV. Advanced theory and applications. 3. Prediction. Integration over segment.
void main_IV_3(void){
  int j,jmax,m,Jmax;
  DOUBLE a,b,dx,A,A1,*yj,*ydj;
  Message("IV. Advanced theory and applications. 3. Prediction. Integration over segment.\n  ");
  _extrmatr(8,&M8);// initialization
  Jmax=1024;yj=(DOUBLE*)calloc(2*Jmax,sizeof(DOUBLE));// memory
  a=-2;b=2;m=7;jmax=24;dx=(b-a)/jmax;// diapazone and step
  for(j=0;j<=jmax;j++)yj[j+2*m]=_f1(a+j*dx);// calculating array
  _predict(-2*m,m,yj+2*m,M8);// predicting left 2*m points
  _predict(2*m,m,yj+2*m+(jmax-2*m),M8);// predicting right 2*m points
  ydj=yj+Jmax;_dual(m,jmax+4*m,yj,ydj,M8);// calculating dual array ydj[]
  A1=_integ_arr(a,b,jmax,m,ydj+m+1,M8);// integration with using dual array
  A=exp(2.)-exp(-2.);Message("I=%g, err=%g",A1,fabs(A1/A-1));
}
// IV. Advanced theory and applications. 4. Principal value of integrals.
void main_IV_4(void){
  int j,jmax,m,_l,Jmax;
  DOUBLE a,b,x,_x,x0,dx,A,A1,*X,*yj;
  Message("IV. Advanced theory and applications. 4. Principal value of integrals.\n  ");
  _extrmatr(8,&M8);// initialization
  X=(DOUBLE*)calloc(16,sizeof(DOUBLE));// array of variables
  Jmax=1024;yj=(DOUBLE*)calloc(Jmax,sizeof(DOUBLE));// memory
  a=0;b=4;_l=0;x0=_zero(a,b,_l,X,_fvp);// searching a singularity
  jmax=14;m=7;dx=(x0-a)/jmax;for(j=-m;j<jmax+m;j++){// array
    x=a+(j+.5)*dx;_x=2*x0-x;//    symmetrizing a singular integrand
    yj[j+m]=_fvp(&x)+_fvp(&_x);//   around x0: f(x)-->f(x)+f(2*x0-x)
  }A1=_integ_arr(a,x0,jmax,m,yj,M8);// intergating a singular part
  jmax=40;dx=(b-2*x0)/jmax;for(j=-m;j<jmax+m;j++){
    x=2*x0+(j+.5)*dx;yj[j+m]=_fvp(&x);
  }A1+=_integ_arr(2*x0,b,jmax,m,yj,M8);// adding regular integral of tail
  A=0;Message("v.p. int=%g, err=%g (x0=%g)",A1,fabs(A1-A),x0);
}
// IV. Advanced theory and applications. 5. Integrating a multi argument function.
void main_IV_5(void){
  int jmax;
  DOUBLE x,y,z,t,b,A,A1,*X;
  Message("IV. Advanced theory and applications. 5. Integrating a multi argument function.\n  ");
  _extrmatr(8,&M8);// initialization
  X=(DOUBLE*)calloc(16,sizeof(DOUBLE));// array of variables
  X[0]=x=.8;X[1]=y=1.3;X[2]=z=.9;X[3]=t=1.1;X[8]=b=8.5;// x,y,z,t\lesssim1 and exp(-b^2/2)=0
  X[9]=X[10]=X[11]=X[12]=jmax=26;A1=_fX(X);// integral
  A=exp(-(x*x+y*y+z*z+t*t)/2);// exact value
  Message("I(%g,%g,%g,%g)=%g, err=%g",x,y,z,t,A1,fabs(A1-A));
}
// IV. Advanced theory and applications. 6. Fourier transform.
void main_IV_6(void){
  struct limit li;
  int j,j0,jmax,m,_l;
  DOUBLE a,u,b,eps,deps,A,A1,*X;
  Message("IV. Advanced theory and applications. 6. Fourier transform.\n  ");
  _extrmatr(8,&M8);_limit_mem(&li);// initialization
  X=(DOUBLE*)calloc(16,sizeof(DOUBLE));// array of variables
  j0=16;m=7;_l=2;X[3]=a=1;X[1]=u=1./a;u=fabs(u);// parameters (scaled on a)
  for(j=1;j<j0;j++){// array
    deps=0.016;if(0<u*a||u*a<1)deps*=sqrt(u*a);// fit "by hand" (морока)
    X[0]=eps=(j+.5)*(j+.5)*deps*deps;// eps-->0
    b=8.5*a/eps;// int_0^infty=int_0^b with exp(-8.5^2/2)=0
    jmax=(int)(5*b/a+4*b*u+32);// sum of jmax over each function
    li.fxj[j]=_integ_quick(0,b,jmax,m,_l,X,M8,_fFt);li.erj[j]=0;//  calculating Fourier
  }A1=_limit_skip(j0,li,M8);//                                      transform, eps-->0
  A=exp(-u*a);Message("F(%g)=%g, err=%g",u,A1,fabs(A1-A));
}
// IV. Advanced theory and applications. 7. Divergent sums.
void main_IV_7(void){
  struct random ra;
  struct limit li;
  int j,j0,init;
  DOUBLE eps,deps,A,A1;
  Message("IV. Advanced theory and applications. 7. Divergent sums.\n  ");
  _extrmatr(8,&M8);_limit_mem(&li);init=1;_randinit(init,&ra);// initialization
  j0=16;;for(j=0;j<j0;j++)li.erj[j]=1E-14;// random error in an initial data
  deps=0.01;for(j=0;j<j0;j++){eps=(j+.5)*deps;li.fxj[j]=_Saj(&eps);}// array
  A1=_limit_skip(j0,li,M8);// taking the limit eps-->0
  A=.25;Message("lim=%g, err=%g",A1,fabs(A1/A-1));
}
// // // // // // // // // // // // // // // // // // // // // // // // // //
//                               teacher
//
void _teacher(int Roman,int arabic){
  if(Roman==1){// I. Basic theory
    if(arabic==1)main_I_1();// 1. Differentiation
    if(arabic==2)main_I_2();// 2. Integration
    if(arabic==3)main_I_3();// 3. Interpolation
    if(arabic==4)main_I_4();// 4. Runge integration
    if(arabic==5)main_I_5();// 5. Generator of random numbers
  }if(Roman==2){// II. Ultraquick calculations
    if(arabic==1)main_II_1();// 1. First derivative
    if(arabic==2)main_II_2();// 2. Second derivative
    if(arabic==3)main_II_3();// 3. Integral
    if(arabic==4)main_II_4();// 4. High-order spline
    if(arabic==5)main_II_5();// 5. First derivative of high-order spline
  }if(Roman==3){// III. Treatment of splines.
    if(arabic==1)main_III_1();// 1. Calculation of step and making a spline.
    if(arabic==2)main_III_2();// 2. Doubling a spline.
    if(arabic==3)main_III_3();// 3. Dual transformation of a spline.
    if(arabic==4)main_III_4();// 4. Differentiation of a spline.
    if(arabic==5)main_III_5();// 5. Calculation of errors in a spline data.
    if(arabic==6)main_III_6();// 6. Analiticheskoe prodolzhenie of a spline.
  }if(Roman==4){// IV. Advanced theory and applications.
    if(arabic==1)main_IV_1();// 1. Calculation of infinite sums.
    if(arabic==2)main_IV_2();// 2. Taking limits.
    if(arabic==3)main_IV_3();// 3. Prediction. Integration over segment.
    if(arabic==4)main_IV_4();// 4. Principal value of integrals.
    if(arabic==5)main_IV_5();// 5. Integrating a multi argument function.
    if(arabic==6)main_IV_6();// 6. Fourier transform.
    if(arabic==7)main_IV_7();// 7. Divergent sums.
} }