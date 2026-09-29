/*vmc.h*/

#ifndef _VMC_H
#define _VMC_H

#include "main.h"
#include "utils.h"
#include <math.h>

DOUBLE Psi(struct Walker Walker);

DOUBLE U(struct Walker Walker);
DOUBLE U_sum_of_two_wavefunctions(struct Walker *W);
DOUBLE Uvector(DOUBLE **R, int w);
DOUBLE OneBodyUWalker(struct Walker Walker);

DOUBLE OneBodyU(DOUBLE x, DOUBLE y, DOUBLE z, int i, int w);
void OneBodyFp(DOUBLE *Fx, DOUBLE *Fy, DOUBLE *Fz, DOUBLE x, DOUBLE y, DOUBLE z, int i, int w);
DOUBLE OneBodyE(DOUBLE x, DOUBLE y, DOUBLE z, int i, int w);

DOUBLE WalkerEnergy2(struct Walker *W);

void VMCMoveOneByOne(int w);
void VMCMoveAll(int w);
void VMCMoveDrift(int w);
void VMCMoveDriftOneByOne(int w);
void VMCMoveSumOfTwoWavefunctions(int w);

void ClassicalMoveOneByOne(int w);

int CheckInteractionConditionWF(const DOUBLE x, const DOUBLE y, const DOUBLE z, const DOUBLE r2);
int CheckInteractionConditionPotential(const DOUBLE x, const DOUBLE y, const DOUBLE z, const DOUBLE r2);
int CheckParticlesInTheBox(const struct Walker W);
int CheckWalkersAlive(void);
int CheckWalkersInTheBox(void);

void UpdateScalableTable(struct Walker *W);
int FindCell(DOUBLE x, DOUBLE y, DOUBLE z);
int FindNearestCells(DOUBLE x, DOUBLE y, DOUBLE z, int *s);
void AdjustCenterOfMassWalker(struct Walker *Walker);

/**************************** INLINE *************************************/
#ifdef __linux__
#  define INLINE static inline
#else
#  define INLINE static __inline
#endif

/********************** Reduce To The Box XYZ *****************************/
INLINE void ReduceToTheBoxXYZ(DOUBLE *x, DOUBLE *y, DOUBLE *z) {

//x=x-L*int(x/L)
//x=x-L*int(2*x/L)

#ifdef BC_3DPBC_CUBE // infinite 3D system
#ifdef IMAGES_FLOOR
  *x -= (floor)(*x*L_inv_x)*Lx;
  *y -= (floor)(*y*L_inv_y)*Ly;
  *z -= (floor)(*z*L_inv_z)*Lz;
#else
  if(*x < 0) *x += Lx;
  if(*x > Lx) *x -= Lx;
  if(*y < 0) *y += Ly;
  if(*y > Ly) *y -= Ly;
  if(*z < 0) *z += Lz;
  if(*z > Lz) *z -= Lz;

  if(*x>Lx || *x<0 || *y>Ly || *y<0 || *z>Lz || *z<0) { // in 3D PBC (x,y,z) must be inside of the box
    if(verbosity > 1) Warning("Reduce to the box : (%" LF ", %" LF ", %" LF "); L = %" LF "\n", *x, *y, *z, Lx);
    ReduceToTheBoxXYZ(x, y, z);
  }
#endif
#endif

#ifdef BC_3DPBC_TRUNCATED_OCTAHEDRON
  // reduce to the cube
  *x -= (floor)(*x*L_inv_x)*Lx;
  *y -= (floor)(*y*L_inv_y)*Ly;
  *z -= (floor)(*z*L_inv_z)*Lz;

  // check the truncated regions
  if(fabs(*x-Lhalf)+fabs(*y-Lhalf)+fabs(*z-Lhalf)>0.75*L) {
    *x += (*x>Lhalf)?(-Lhalf):(Lhalf);
    *y += (*y>Lhalf)?(-Lhalf):(Lhalf);
    *z += (*z>Lhalf)?(-Lhalf):(Lhalf);
  }

#ifdef SECURE
  if(*x>Lx && *y>Ly && *z>Lz) Error("bad particle coordinates (%" LF ", %" LF ", %" LF "); L = %" LF "\n", *x, *y, *z, L);
  if(*x<0 && *y<0 && *z<0) Error("bad particle coordinates (%" LF ", %" LF ", %" LF ")" LF "\n", *x, *y, *z);
#endif
#endif

#ifdef BC_SPHERE // put make to the sphere
  Rnow = sqrt((*x)*(*x) + (*y)*(*y) + (*z)*(*z)); // calculate distance from the center
  Rnow /= Rsphere; //  normalize to the sphere radius
  *x /= Rnow;
  *y /= Rnow;
  *z /= Rnow;
#endif

#ifdef BC_2DPBC_SQUARE // infinite 2D system
#ifdef IMAGES_FLOOR
  *x -= (floor)(*x*L_inv_x)*Lx;
  *y -= (floor)(*y*L_inv_y)*Ly;
#else
  if(*x < 0) *x += Lx;
  if(*x > Lx) *x -= Lx;
  if(*y < 0) *y += Ly;
  if(*y > Ly) *y -= Ly;

  if(*x>Lx || *x<0 || *y>Ly || *y<0 ) { // in 2D PBC (x,y) must be inside of the box
    if(verbosity > 1) Warning("Reduce to the box : (%" LF ", %" LF ", %" LF "); L = %" LF "\n", *x, *y, *z, Lx);
    ReduceToTheBoxXYZ(x, y, z);
  }
#endif
#endif

#ifdef BC_2DPBC_HEXAGON // infinite 2D system
  *x -= (floor)(*x*L_inv_x)*Lx;
  if(fabs(*x-L_half_x)+sqrt(3.)*fabs(*y-L_half_y)>Lx) { // check angles
    *x -= (*x>L_half_x)?(L_half_x):(-L_half_x);
    *y -= (*y>L_half_y)?(0.75*Ly):(-0.75*Ly);
    //*y -= (*y>L_half_y)?(L_half_x*sqrt(3.)):(-L_half_x*sqrt(3.));
  }
  if(*x>Lx || *x<0 || fabs(*x-L_half_x)+sqrt(3.)*fabs(*y-L_half_y)>Lx) {
    if(verbosity > 1) Warning("Reduce to the box : (%" LF ", %" LF ", %" LF "); L = %" LF "\n", *x, *y, *z, Lx);
    ReduceToTheBoxXYZ(x, y, z);
  }
#endif

#ifdef BC_1DPBC_Z // tube
#ifdef IMAGES_FLOOR
  *z -= (floor)(*z*L_inv_z)*Lz;
#else
  if(*z < 0.) *z += Lz;
  if(*z > L)  *z -= Lz;
  if(*z>Lz || *z<0) { // in 1D PBC (z) coordinate must be inside of the box
    if(verbosity > 1) Warning("Reduce to the box : (%" LF ", %" LF ", %" LF "); L = %" LF "\n", *x, *y, *z, Lz);
    ReduceToTheBoxXYZ(x, y, z);
  }
#endif
#endif

#ifdef BC_1DPBC_X // tube
#ifdef IMAGES_FLOOR
  *x -= (floor)(*x*L_inv_x)*Lx;
#else
  if(*x < 0.) *x += Lx;
  if(*x > L)  *x -= Lx;
  if(*x>Lx || *x<0) { // in 1D PBC (x) coordinate must be inside of the box
    if(verbosity > 1) Warning("Reduce to the box : (%" LF ", %" LF ", %" LF "); L = %" LF "\n", *x, *y, *z, Lz);
    ReduceToTheBoxXYZ(x, y, z);
  }
#endif
#endif
}

/********************** Reduce To The Box *********************************/
INLINE void ReduceToTheBox(DOUBLE *x) {
  ReduceToTheBoxXYZ(&x[0], &x[1], &x[2]);
}

/********************** Reduce Walker To The Box **************************/
INLINE void ReduceWalkerToTheBox(struct Walker *Walker) {
  int i;

#ifdef BC_ABSENT
  return;
#endif

  for(i=0; i<N; i++) ReduceToTheBoxXYZ(&Walker->x[i], &Walker->y[i], &Walker->z[i]);
}

/********************** Reduce Vector To The Box **************************/
INLINE void ReduceVectorToTheBox(DOUBLE **R) {
  int i;

#ifdef BC_ABSENT
  return;
#endif

  for(i=0; i<N; i++) ReduceToTheBoxXYZ(&R[i][0], &R[i][1], &R[i][2]);
}

/*************************** Find Image **************************************/
/* find for the particle i the image of the pair particle j within L/2 */
INLINE DOUBLE FindNearestImageCenterOfMass(DOUBLE *x, DOUBLE *y, DOUBLE *z) {
#ifdef BC_ABSENT
  CaseX(*x = (floor)(*x*L_inv_x + 0.5)*Lx);
  CaseY(*y = (floor)(*y*L_inv_y + 0.5)*Ly);
  CaseZ(*z = (floor)(*z*L_inv_z + 0.5)*Lz);
#endif
  return *x**x + *y**y + *z**z;
}
/*************************** Find Image **************************************/
/* find for the particle i the image of the pair particle j within L/2 */
INLINE DOUBLE FindNearestImage(DOUBLE *x, DOUBLE *y, DOUBLE *z) {

#ifdef BC_3DPBC_CUBE // homogeneous 3D system
#ifdef IMAGES_FLOOR
  *x -= (floor)(*x*L_inv_x + 0.5)*Lx;
  *y -= (floor)(*y*L_inv_y + 0.5)*Ly;
  *z -= (floor)(*z*L_inv_z + 0.5)*Lz;
#endif

#ifdef IMAGES_IF_IF
  if(*x > L_half_x) {
    *x -= Lx;
  }
  else if(*x < -L_half_x) {
    *x += Lx;
  }
  if(*y > L_half_y) {
    *y -= Ly;
  }
  else if(*y < -L_half_y) {
    *y += Ly;
  }
  if(*z > L_half_z) {
    *z -= Lz;
  }
  else if(*z < -L_half_z) {
    *z += Lz;
  }
#endif

#ifdef IMAGES_ABS_IF
  if(fabs(*x) > L_half_x) {
    if(*x > 0.) {
      *x -= Lx;
    }
    else {
      *x += Lx;
    }
  }
  if(fabs(*y) > L_half_y) {
    if(*y > 0.) {
      *y -= Ly;
    }
    else {
      *y += Ly;
    }
  }
  if(fabs(*z) > L_half_z) {
    if(*z > 0.) {
      *z -= Lz;
    }
    else {
      *z += Lz;
    }
  }
#endif

#ifdef SECURE
  if(fabs(*x)>Lx && fabs(*y)>Ly && fabs(*z)>Lz) Error("bad particle coordinates (%" LF ", %" LF ", %" LF "); L = %" LF "\n", *x, *y, *z, Lx);
#endif
#endif

#ifdef  BC_3DPBC_TRUNCATED_OCTAHEDRON
  // reduce to the cube
  *x -= (floor)(*x*L_inv_x + 0.5)*Lx;
  *y -= (floor)(*y*L_inv_y + 0.5)*Ly;
  *z -= (floor)(*z*L_inv_z + 0.5)*Lz;

  // check the truncated regions
  if(fabs(*x)+fabs(*y)+fabs(*z)>0.75*L) {
    *x += (*x>0)?(-Lhalf):(Lhalf);
    *y += (*y>0)?(-Lhalf):(Lhalf);
    *z += (*z>0)?(-Lhalf):(Lhalf);
  }

#ifdef SECURE
  if(fabs(*x)>Lhalf && fabs(*y)>Lhalf && fabs(*z)>Lhalf) 
    Error("bad particle coordinates (%" LF ", %" LF ", %" LF "); L = %" LF "\n", *x, *y, *z, Lx);
#endif
#endif

#ifdef BC_2DPBC_NON_ORTHOGONAL
  *x -= (int)(*x*Crystal.e1x*Crystal.e1abs_inv2)*Crystal.e1x + (int)(*x*Crystal.e2x*Crystal.e2abs_inv2)*Crystal.e2x;
  *y -= (int)(*y*Crystal.e1y*Crystal.e1abs_inv2)*Crystal.e1y + (int)(*y*Crystal.e2y*Crystal.e2abs_inv2)*Crystal.e2y;
#endif

#ifdef BC_2DPBC_SQUARE // homogeneous 2D system (x,y) plane
#ifdef IMAGES_FLOOR
  *x -= (floor)(*x*L_inv_x + 0.5)*Lx;
  *y -= (floor)(*y*L_inv_y + 0.5)*Ly;
#endif

#ifdef IMAGES_IF_IF
  if(*x > L_half_x) {
    //if((*x -Lx) != (*x-(floor)(*x*L_inv_x + 0.5)*Lx)) Message("aaa x");
    *x -= Lx;
  }
  else if(*x < -L_half_x) {
    //if((*x +Lx) != *x-(floor)(*x*L_inv_x + 0.5)*Lx) Message("aaa x2");
    *x += Lx;
  }
  if(*y > L_half_y) {
    //if(*y -Ly != *y-(floor)(*y*L_inv_y + 0.5)*Ly) Message("aaa y");
    *y -= Ly;
  }
  else if(*y < -L_half_y) {
    //if(*y +Ly != *y-(floor)(*y*L_inv_y + 0.5)*Ly) Message("aaa y2");
    *y += Ly;
  }
#endif

#ifdef IMAGES_ABS_IF
  if(fabs(*x) > L_half_x) {
    if(*x > 0.) {
      *x -= Lx;
    }
    else {
      *x += Lx;
    }
  }
  if(fabs(*y) > L_half_y) {
    if(*y > 0.) {
      *y -= Ly;
    }
    else {
      *y += Ly;
    }
  }
#ifdef SECURE
  if(fabs(*x)>Lx && fabs(*y)>Ly) Error("bad particle coordinates (%" LF ", %" LF ", %" LF "); L = %" LF "\n", *x, *y, *z, Lx);
#endif
#endif
#endif

#ifdef BC_2DPBC_HEXAGON // infinite 2D system
  *x -= (floor)(*x*L_inv_x + 0.5)*Lx;
  if(fabs(*x)+sqrt(3.)*fabs(*y)>Lx) { // check angles
    *x -= (*x>0)?(L_half_x):(-L_half_x);
    *y -= (*y>0)?(0.75*Ly):(-0.75*Ly);
  }
#endif

#ifdef BC_1DPBC_Z // tube (z direction)
#ifdef IMAGES_FLOOR
  *z -= (floor)(*z*L_inv_z + 0.5)*Lz;
#endif

#ifdef IMAGES_IF_IF
  if(*z > L_half_z) {
    *z -= Lz;
  }
  else if(*z < -L_half_z) {
    *z += Lz;
  }
#endif

#ifdef IMAGES_ABS_IF
  if(fabs(*z) > L_half_z) {
    if(*z > 0.) {
      *z -= Lz;
    }
    else {
      *z += Lz;
    }
  }
#endif

#ifdef SECURE
  if(fabs(*z)>Lz) Error("bad particle coordinates (%" LF ", %" LF ", %" LF "); L = %" LF "\n", *x, *y, *z, Lz);
#endif
#endif

#ifdef BC_1DPBC_X // tube (x direction)
  if(fabs(*x) > L_half_x) {
    if(*x > 0.) {
      *x -= Lx;
    }
    else {
      *x += Lx;
    }
  }
#ifdef SECURE
  if(fabs(*x)>Lx) Error("bad particle coordinates (%" LF ", %" LF ", %" LF "); L = %" LF "\n", *x, *y, *z, Lz);
#endif
#endif

#ifdef TRIAL_1D
#ifdef BC_1DPBC_X // in 1D trap return z*z
  return *x**x;
#else
  return *z**z;
#endif
#endif

#ifdef BC_SPHERE
  return *x**x + *y**y + *z**z;
#endif

#ifdef TRIAL_2D
  return *x**x + *y**y;
#endif
#ifdef TRIAL_3D
  return *x**x + *y**y + *z**z;
#endif
}

/********************************** Overlapping *******************************/
#ifndef HARD_SPHERE
#  define Overlapping(R) 0
#  define OverlappingWalker(W) 0
#  define CheckOverlapping() 0
#  define CheckWalkerOverlapping(W) 0
#else
int Overlapping(DOUBLE **R);
int OverlappingWalker(const struct Walker* W);
int CheckOverlapping(void);
int CheckWalkerOverlapping(const struct Walker W);
#endif

int CheckWalkerHyperradiusOverlapping(const struct Walker W);
int CheckWalkerHyperradiusOverlappingCount(struct Walker *W);
int CheckVectorHyperradiusOverlapping(DOUBLE **R);
void U_psiT_sigma_inv(struct Walker Walker);

/***************************** 3D image ***************************************/
#ifdef BC_ABSENT
#define FindNearestImage3D(r2, x, y, z) (r2 = x*x + y*y + z*z)
#else //!BC_ABSENT
#ifdef BC_1DPBC_Z // z direction
#define FindNearestImage3D(r2, x, y, z) {\
  z -= (floor)(z*L_inv_z + 0.5)*Lz;\
  r2 = z*z;}
#else
#ifdef BC_1DPBC_X // x direction
#define FindNearestImage3D(r2, x, y, z) {\
  x -= (floor)(x*L_inv_x + 0.5)*Lx;\
  r2 = x*x;}
#else
#define FindNearestImage3D(r2, x, y, z) (r2 = FindNearestImage(&x,&y,&z))
#endif // BC_1DPBC_X
#endif // BC_1DPBC_Z
#endif // BC_ABSENT

void ExcludedVolumeCoordinatesWrap(struct Walker Walker);
void ExcludedVolumeCoordinatesUnwrap(struct Walker Walker);

#endif
