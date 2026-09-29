/*speckles.h*/
#ifndef _SPECKLES_H_
#define _SPECKLES_H_

#define N_POINTS 513
//#define N_POINTS 129
//#define N_POINTS 257

#define NR_END 1
#define FREE_ARG char*

void SpecklesLoad(void);
// usingg splines
DOUBLE SpecklesEpotSplines(DOUBLE x1, DOUBLE x2); // Speckle potential
DOUBLE SpecklesUSplines(DOUBLE x1, DOUBLE x2); // logarithm of one-body w.f.
DOUBLE SpecklesFppSplines(DOUBLE x1, DOUBLE x2); // second derivative
void SpecklesFpSplines(DOUBLE *Fx, DOUBLE *Fy, DOUBLE x1, DOUBLE x2); // logarithmic derivative

DOUBLE SpecklesEpot9point(DOUBLE x1, DOUBLE x2); // Speckle potential
DOUBLE SpecklesU9point(DOUBLE x1, DOUBLE x2); // logarithm of one-body w.f.
DOUBLE SpecklesFpp9point(DOUBLE x1, DOUBLE x2); // second derivative
void SpecklesFp9point(DOUBLE *Fx, DOUBLE *Fy, DOUBLE x1, DOUBLE x2); // logarithmic derivative

DOUBLE SpecklesEpot16point(DOUBLE x1, DOUBLE x2); // Speckle potential
DOUBLE SpecklesU16point(DOUBLE x1, DOUBLE x2); // logarithm of one-body w.f.
DOUBLE SpecklesFpp16point(DOUBLE x1, DOUBLE x2); // second derivative
void SpecklesFp16point(DOUBLE *Fx, DOUBLE *Fy, DOUBLE x1, DOUBLE x2); // logarithmic derivative

DOUBLE SpecklesEpot16pointContinuous(DOUBLE x1, DOUBLE x2); // Speckle potential
DOUBLE SpecklesU16pointContinuous(DOUBLE x1, DOUBLE x2); // logarithm of one-body w.f.
DOUBLE SpecklesFpp16pointContinuous(DOUBLE x1, DOUBLE x2); // second derivative
void SpecklesFp16pointContinuous(DOUBLE *Fx, DOUBLE *Fy, DOUBLE x1, DOUBLE x2); // logarithmic derivative

void SpecklesInitialize4(void);
DOUBLE SpecklesEpot4point(DOUBLE x, DOUBLE y);
DOUBLE SpecklesU4point(DOUBLE x1, DOUBLE x2); // logarithm of one-body w.f.
DOUBLE SpecklesFpp4point(DOUBLE x1, DOUBLE x2); // second derivative
void SpecklesFp4point(DOUBLE *Fx, DOUBLE *Fy, DOUBLE x1, DOUBLE x2); // logarithmic derivative

DOUBLE *vector(long nl, long nh); //arrays.c
void free_vector(DOUBLE *v, long nl, long nh); //arrays.c
void spline(DOUBLE x[], DOUBLE y[], int n, DOUBLE yp1, DOUBLE ypn, DOUBLE y2[]);
void splie2(DOUBLE x1a[], DOUBLE x2a[], DOUBLE **ya, int m, int n, DOUBLE **y2a);
void splie2x1(DOUBLE x1a[], DOUBLE x2a[], DOUBLE **ya, int m, int n, DOUBLE **y2a);
void splint(DOUBLE xa[], DOUBLE ya[], DOUBLE y2a[], int n, DOUBLE x, DOUBLE *y);
void splin2(DOUBLE x1a[], DOUBLE x2a[], DOUBLE **ya, DOUBLE **y2a, int m, int n, DOUBLE x1, DOUBLE x2, DOUBLE *y);
void nrerror(char error_text[]);
void splder1(DOUBLE xa[], DOUBLE ya[], DOUBLE y2a[], int n, DOUBLE x, DOUBLE *y1);
void splder2(DOUBLE xa[], DOUBLE ya[], DOUBLE y2a[], int n, DOUBLE x, DOUBLE *y2);
void nrerror(char error_text[]);

DOUBLE **matrix(long nrl, long nrh, long ncl, long nch);
DOUBLE *vector(long nl, long nh);
void free_vector(DOUBLE *v, long nl, long nh);
DOUBLE ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void splin2x1(DOUBLE x1a[], DOUBLE x2a[], DOUBLE **ya, DOUBLE **y2ax2, int m, int n, DOUBLE x1, DOUBLE x2, DOUBLE *y, DOUBLE *d11y, DOUBLE *d12y);
void splintlong(DOUBLE xa[], DOUBLE ya[], DOUBLE y2a[], int n, DOUBLE x, DOUBLE *y, DOUBLE *d1y, DOUBLE *d2y);

#ifdef SPLINES_2D
#  define SpecklesEpot SpecklesEpotSplines
#  define SpecklesU SpecklesUSplines
#  define SpecklesF SpecklesFSplines
#  define SpecklesFp SpecklesFpSplines
#  define SpecklesFpp SpecklesFppSplines
#endif

#ifdef POINT9
#  define SpecklesEpot SpecklesEpot9point
#  define SpecklesU SpecklesU9point
#  define SpecklesF SpecklesF9point
#  define SpecklesFp SpecklesFp9point
#  define SpecklesFpp SpecklesFpp9point
#endif

#ifdef POINT16
#  define SpecklesEpot SpecklesEpot16point
#  define SpecklesU SpecklesU16point
#  define SpecklesF SpecklesF16point
#  define SpecklesFp SpecklesFp16point
#  define SpecklesFpp SpecklesFpp16point
#endif

#ifdef POINT16CONTINUOUS
#  define SpecklesEpot SpecklesEpot16pointContinuous
#  define SpecklesU SpecklesU16pointContinuous
#  define SpecklesF SpecklesF16pointContinuous
#  define SpecklesFp SpecklesFp16pointContinuous
#  define SpecklesFpp SpecklesFpp16pointContinuous
#endif

#ifdef POINT4
#  define SpecklesEpot SpecklesEpot4point
#  define SpecklesU SpecklesU4point
#  define SpecklesF SpecklesF4point
#  define SpecklesFp SpecklesFp4point
#  define SpecklesFpp SpecklesFpp4point
#endif

#ifdef POINT4
extern DOUBLE v[N_POINTS][N_POINTS];
extern DOUBLE vx[N_POINTS][N_POINTS];
extern DOUBLE vy[N_POINTS][N_POINTS];
extern DOUBLE vxx[N_POINTS][N_POINTS];
extern DOUBLE vxy[N_POINTS][N_POINTS];
extern DOUBLE vyy[N_POINTS][N_POINTS];
extern DOUBLE vxxx[N_POINTS][N_POINTS];
extern DOUBLE vxxy[N_POINTS][N_POINTS];
extern DOUBLE vxyy[N_POINTS][N_POINTS];
extern DOUBLE vyyy[N_POINTS][N_POINTS];
extern DOUBLE vxxxy[N_POINTS][N_POINTS];
extern DOUBLE vxyyy[N_POINTS][N_POINTS];

extern DOUBLE w[N_POINTS][N_POINTS];
extern DOUBLE wx[N_POINTS][N_POINTS];
extern DOUBLE wy[N_POINTS][N_POINTS];
extern DOUBLE wxx[N_POINTS][N_POINTS];
extern DOUBLE wxy[N_POINTS][N_POINTS];
extern DOUBLE wyy[N_POINTS][N_POINTS];
extern DOUBLE wxxx[N_POINTS][N_POINTS];
extern DOUBLE wxxy[N_POINTS][N_POINTS];
extern DOUBLE wxyy[N_POINTS][N_POINTS];
extern DOUBLE wyyy[N_POINTS][N_POINTS];
extern DOUBLE wxxxy[N_POINTS][N_POINTS];
extern DOUBLE wxyyy[N_POINTS][N_POINTS];
#else
extern DOUBLE **v,**vx,**vy,**vxx,**vxy,**vyy,**vxxx,**vxxy,**vxyy,**vyyy,**vxxxy,**vxyyy;
extern DOUBLE **w,**wx,**wy,**wxx,**wxy,**wyy,**wxxx,**wxxy,**wxyy,**wyyy,**wxxxy,**wxyyy;
#endif

#endif
