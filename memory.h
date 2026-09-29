/*memory.h*/

#ifndef _MEMORY_H_
#define _MEMORY_H_
#include <malloc.h>
#include "trial.h"
#include "main.h"
#include "utils.h"

void AllocateWFGrid(struct Grid *G, unsigned long int size);
void AllocateSingleWalker(struct Walker* W);
void AllocateWalkers(void);
void AllocateGrids(void);
void CopyWalker(struct Walker *out, const struct Walker *in);
void CopyWalkerCoord(struct Walker *out, const struct Walker *in);
void CopyRawVectorToWalkerDisplace_i_j(DOUBLE *r, struct Walker *W, int index1, DOUBLE dx1, int index2, DOUBLE dx2);
int WalkerCompare(const struct Walker x, const struct Walker y);
int VectorToWalkerCompare(const DOUBLE **x, const struct Walker y);

extern void* Calloc(const char* name, unsigned length, size_t size);
extern void* CallocContiguous2D(const char* name, unsigned dim1, unsigned dim2, char *type);
extern void* CallocContiguous3D(const char* name, unsigned dim1, unsigned dim2, unsigned dim3, char *type);

#ifdef MEMORY_CONTIGUOUS
#define ArrayCalloc2D(p2D, name, i, dim1, dim2, type_static, type_dynamic) \
   p2D = (type_static**) CallocContiguous2D(name, dim1, dim2, type_dynamic)
#else 
#define ArrayCalloc2D(p2D, name, i, dim1, dim2, type_static, type_dynamic) { \
  p2D = (type_static**) Calloc(name, dim1, sizeof(type_static*)); \
  for(i=0; i<dim1; i++) p2D[i]  = (type_static*) Calloc(name, dim2, sizeof(type_static)); \
}
#define ArrayCalloc3D(p3D, name, i, dim1, j, dim2, dim3, type_static, type_dynamic) { \
  p3D = (type_static***) Calloc(name, dim1, sizeof(type_static**)); \
  for(i=0; i<dim1; i++) { \
    p3D[i] = (type_static**) Calloc(name, dim2, sizeof(type_static*)); \
    for(j=0; j<dim2; j++) { \
      p3D[i][j] = (type_static*) Calloc(name, dim3, sizeof(type_static)); \
    } \
  } \
}
#endif

#ifdef MEMORY_CONTIGUOUS
#define ArrayCopy1D(pointer_out, pointer_in, i, dim) \
  memcpy(pointer_out, pointer_in, sizeof(*pointer_in)*dim);
#else 
#define ArrayCopy1D(pointer_out, pointer_in, i, dim) \
  for(i=0; i<dim; i++) pointer_out[i] = pointer_in[i];
#endif

#ifdef MEMORY_CONTIGUOUS
#define ArrayCopy2D(pointer_out, pointer_in, i, dim1, j, dim2) \
  memcpy(pointer_out[0], pointer_in[0], sizeof(**pointer_in)*dim1*dim2);
#else 
#define ArrayCopy2D(pointer_out, pointer_in, i, dim1, j, dim2) \
  for(i=0; i<dim1 ; i++) { \
    for(j=0; j<dim2; j++) { \
      pointer_out[i][j] = pointer_in[i][j]; \
    } \
  }
#endif
#define ArrayCopy3D(pointer_out, pointer_in, i, dim1, j, dim2, k, dim3) \
  for(i=0; i<dim1 ; i++) { \
    for(j=0; j<dim2; j++) { \
      for(k=0; k<dim3; k++) { \
      pointer_out[i][j][k] = pointer_in[i][j][k]; \
    } \
   } \
  }

#ifdef MEMORY_CONTIGUOUS
#define ArrayEmpty1D(pointer_array, i, dim) \
  memset(pointer_array, 0, dim*sizeof(*pointer_array));
#else
#define ArrayEmpty1D(pointer_array, i, dim) \
  for(i=0; i<dim; i++) { \
    pointer_array[i] = 0; \
  }
#endif

//#ifdef MEMORY_CONTIGUOUS
//#define ArrayEmpty2D(pointer_array, i, dim1, j, dim2, type_static) \
//  memset(pointer_array[0], 0, dim1*dim2*sizeof(type_static)); // does not work!
//#else
#define ArrayEmpty2D(pointer_array, i, dim1, j, dim2, type_static) \
  for(i=0; i<dim1; i++) { \
    for(j=0; j<dim2; j++) { \
      pointer_array[i][j] = 0; \
    } \
  }
//#endif

#endif
