// // // // // // // // // // // // // // // // // // // // // // // // // //
//                                 structures
//
struct random{// structure for random number generator
  DOUBLE c0,c1,c2,c3,c4;// constants
  DOUBLE arr[16];// array of 16 rectangularly distributed random numbers (0<=arr[j]<1; 53 bits)
  DOUBLE xi,Xi,aa;// rectangularly distributed random numbers (0<=xi,Xi,aa<1; 53 bits)
  int index;// index
};

struct Runge{// structure for Runge integration
// d/dx y1(X)=f1(X,y1(X),... yn(X)),
// ...
// d/dx yn(X)=fn(X,y1(X),... yn(X)),
// x=X[_l], y1(a)=y01,...yn(a)=y0n
  DOUBLE a;// initial point
  DOUBLE b;// final point
  int _l;// number of the variable (x=X[_l])
  int nn;// number of steps (nn=(b-a)/dx)
  int n;// number of functions
  int s;// save values of functions y1(xj),...ys(xj) (1<=j<=nn, xj=a+dx*j)
  DOUBLE *arr;// array for saving y1(xj),...ys(xj) (memory s*(nn+1) MUST BE allocated)
  DOUBLE *yl;// reserve array
  DOUBLE *fk;// array of functions
  DOUBLE *mik0,*mik1,*mik2,*mik3;// coefficients m0,m1,m2,m3
  DOUBLE *yk0;// initial values
  DOUBLE *ykxj;// final values
};

struct exint{// structure for extrapolation and integration
  DOUBLE *Dnl;// coefficients D_n^l (2M*M)
  DOUBLE *Aknm;// coefficients A_kn^m ((M+1)*2M*M)
  DOUBLE *Wkm;// coefficients W_k^m (2M*M)
  DOUBLE *Ilm;// coefficients I_l^m (2M*M)
  DOUBLE *Blm;// coefficients B_l^m (M*M)
  DOUBLE *Dkm;// coefficients D_k^m (M*M)
  DOUBLE *cl_du;// coefficients cl for dual grid  (M)
  DOUBLE *v1n;// values of 1./n (64)
  DOUBLE *yxj;// auxiliary array (64)
  int M;// maximal order of schemes
};

struct limit{// structure for taking limiting value numerically
  DOUBLE *fxj;// array of function values (16)
  DOUBLE *erj;// array of their errors (16)
  DOUBLE *raj;// array of random numbers (16)
  DOUBLE *a;// symmetrized array (32)
};

struct spline{// structure for a spline y[j]=f(x0+j*dx,x1,x2,x3,...), j=0,1,2,3,...jm
  int m,j0;// 0<=m<=7 is order of scheme, j0=2*m+1
  int step;// step*dx is step for differentiation (1<=step<=stma),
  int stma;// stma is maximal value of step,
  int j0st;// j0st=(2*m+1)*stma is number of points to the left and to the right of the interpolated region
  int jmax;// jmax>=1 (dx*jmax is lenght of the interpolated region)
  int jm;// jm+1=jmax+2*j0st+1 is total number of points in the spline
  int Jm;// maximal number of points in the spline: jm+1<=Jm
  DOUBLE a,b;// diapazon of interpolated region (f(a)=y[j0st], f(b)=y[j0st+jmax])
  DOUBLE dx,_dx;// dx=(b-a)/jmax>0 is step of the spline, _dx=1./dx
  DOUBLE x0;// x0=a-j0st*dx is first point of the spline
  DOUBLE *y;// array for the spline (Jm)                        memory for these
  DOUBLE *yd;// array for the dual spline (Jm)                  arrays (if they
  DOUBLE *er;// array for errors in spline points (Jm)          are used) must be
  DOUBLE *y1;// array for k-th derivative of the spline (Jm)    allocated apart
};

DOUBLE _Cnm(int n,int m); // binominal coefficients n!/m!/(n-m)!
DOUBLE _po(DOUBLE n,int m); 
DOUBLE _fa(int n); // n! (n>=0)
DOUBLE _Dnlfunc(int n,int j_,struct exint e); // calculating D_n^l
DOUBLE _Aknmfunc(int k,int n,int m,struct exint e); // calculating A_kn^m
DOUBLE _Wkmfunc(int k,int m,struct exint e); // calculation of W_k^m
DOUBLE _Ilmfunc(int l,int m,struct exint e); // calculation of I_l^m
DOUBLE _Blmfunc(int l,int m,struct exint e); // calculation of B_l^m
void _cldu(struct exint e); // calculation of c_l for dual grid
DOUBLE _Dkmfunc(int k,int m,struct exint e); // calculation of \tilde D_k^m
void _extrmatr(int M,struct exint *e); // calculation of matrices for D_n^l and A_kn^m and array of values 1./n
DOUBLE _fxjnmx(DOUBLE dx,int n,int m,int j0,int step,DOUBLE *fxj,struct exint e); // quick direct calculation of n-th derivative of f(x) in (m-n/2)-th order
DOUBLE _fxjkx(DOUBLE xminusx0,DOUBLE dx,int k,int m,int nmax,int j0,int step,DOUBLE *fxj,struct exint e); //quick Tailor series for k-th derivative of f(x) in m-th order at x near to x0
DOUBLE _spline5(DOUBLE x,struct spline s); // ultraquick calculation of 5-th order spline for function f(x) (21 multiplications)
DOUBLE _deriv5(DOUBLE x,struct spline s); // spline derivative
DOUBLE _spline7(DOUBLE x,struct spline s); // ultraquick calculation of 7-th order spline for function f(x) (31 multiplications)
DOUBLE _deriv7(DOUBLE x,struct spline s); //
DOUBLE _spline9(DOUBLE x,struct spline s); // ultraquick calculation of 9-th order spline for function f(x) (37 multiplications)
DOUBLE _deriv9(DOUBLE x,struct spline s);
void _spline_anal(DOUBLE re,DOUBLE im,DOUBLE *Re,DOUBLE *Im,int k,struct spline s,struct exint e); // analytic continuation: Tailor series for k-th derivative of f(z) in m-th order
DOUBLE _calcstep(DOUBLE a,DOUBLE b,int m,int Jmax,int *_jmax,int _l,DOUBLE *X, DOUBLE *y,DOUBLE *yd,struct exint e,DOUBLE (f)(DOUBLE *X)); // calculation of step (b-a)/_jmax of function f(x)=f(X[0],X[1],X[2],X[3],...),
void _splineprep(DOUBLE a,DOUBLE b,int jmax,int m,int stma,int Jm,struct spline *s); // preparing of a structure for a spline
void _splinecalc(int _l,DOUBLE *X,struct spline s,DOUBLE (f)(DOUBLE *X)); // calculating a spline y[j]=f(X[0],X[1],X[2],X[3],...) with X[_l]=x0+j*dx
void _dual(int m,int jm,DOUBLE *y,DOUBLE *yd,struct exint e); // dual transformation in m-th order: yd[j] is between y[j] and y[j+1]
void _dualspline(struct spline s,struct exint e); // calculation of dual spline yd[] for a spline y[]
void _convertdual(struct spline s,struct spline *S); // converting a dual spline in s into a new spline S
void _doublingspline(struct spline *s,struct exint e); // doubling the spline in m-th order:
void _splinederivative(int k,struct spline s,struct exint e); // caclulating k-th derivative of spline in m-th order
void _convertderivative(struct spline s,struct spline *S); // converting k-th derivative of spline s into a new spline S
DOUBLE _splineerror(int Jer,struct spline s,struct exint e); // calculating errors in spline (in m-th order with averaging over 2*Jer+1 values)
void _predict(int J,int m,DOUBLE *fxj,struct exint e); // predicting J (-J) values of function f(x) followed (preceded) by its known 2*m values
void _coeffm(DOUBLE *X,struct Runge *r,void (funcs)(int _l,DOUBLE *X,DOUBLE *yl,DOUBLE *fk)); // calculation of coefficients m0,m1,m2,m3 in Runge integration
void _runge(DOUBLE *X,struct Runge *r, void (funcs)(int _l,DOUBLE *X,DOUBLE *yl,DOUBLE *fk)); // Runge integration of n equations from x=a to x=b with step dx=(b-a)/nn
void _rungematr(int n,struct Runge *r); // memory for Runge matrices
void _randcalc(struct random *s); // generation of rectangularly distributed random numbers 0<=arr[0],...arr[15]<1 of DOUBLE type (53 bits)
void _randinit(int init,struct random *s); // initialization of generator (-2^31<=init<2^31)
DOUBLE _randrect(struct random *s); // generates rectangular distribution (53 bits, 50 tacts)
DOUBLE _randrect_lowquality(struct random *s); // quick rectangular distribution (53 bits, 25 tacts)
DOUBLE _randnorm(struct random *s); // normal distribution (exact formula, 170 tacts)
DOUBLE _zero(DOUBLE x1,DOUBLE x2,int _l,DOUBLE *X,DOUBLE (f)(DOUBLE *X)); // searching a zero of f(X) (or 1/f(X)), where X=X[0],X[1],X[2],X[3]... and
void _accumulate(int n,DOUBLE a,DOUBLE *av,DOUBLE *er); // accumulating average value and error at n-th iteration
DOUBLE _integ_arr(DOUBLE a,DOUBLE b,int jmax,int m,DOUBLE *fxj,struct exint e); // quick integration of f(x) on dx over [a;b] in m-th order;
DOUBLE _integ_quick(DOUBLE a,DOUBLE b,int jmax,int m,int _l,DOUBLE *X, struct exint e,DOUBLE (f)(DOUBLE *X)); // quick integration of f(X[0],X[1],X[2],X[3],...) on dX[_l] over [a;b] in m-th order
DOUBLE _sumaj(int jm,int m,int _l,DOUBLE *X,struct exint e,DOUBLE (aj)(DOUBLE *X)); // summation of aj(X[0],X[1],X[2],X[3],...) with X[_l]=j from j=1 to j=+infty:
DOUBLE _limit_dir(int j0,struct limit l,struct exint e); // direct taking of the limit
DOUBLE _limit_skip(int j0,struct limit l,struct exint e); // advanced taking the limit
DOUBLE _limit_err(int j0,int N,DOUBLE *err,struct limit l,struct exint e, struct random *r,DOUBLE (lim)(int j0,struct limit l,struct exint e)); // calculating average value and error in taking limit, j0=3,4,5,...16
void _limit_mem(struct limit *l); // memory allocation for arrays used for taking limit

extern struct exint M8;