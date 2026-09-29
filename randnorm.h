/*randnorm.h*/

#ifndef __RANDOM_NORMAL_H_
#define __RANDOM_NORMAL_H_

#include <stdlib.h>
#include "main.h"

#define USE_RAN2_GENERATOR // else use system random generator

/*#ifdef __alpha
# define MAXR  2147483647
#else
# define MAXR RAND_MAX
#endif*/

#ifndef RAN2_SEED
# define RAN2_SEED 1023
#endif

//#define Random() ran2(RAN2_SEED)

#define Randomize() srand((unsigned)time(NULL))
#define RandomSys() ((DOUBLE) rand() / (DOUBLE) RAND_MAX)

void RandomNormal(DOUBLE *x1, DOUBLE *x2, const DOUBLE mu, const DOUBLE sigma);
void RandomNormal3(DOUBLE *x1, DOUBLE *x2, DOUBLE *x3, const DOUBLE mu, const DOUBLE sigma);
void RandomNormal3slow(DOUBLE *x1, DOUBLE *x2, DOUBLE *x3, const DOUBLE mu, const DOUBLE sigma);
DOUBLE ran2(long *idum);
extern long int rand_seed;

#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))
#define ON 1
#define OFF 0

void RandInit(void);
#ifdef _OPENMP
#include "parallel.h"

extern int rand_seed_parallel_length;
extern long int rand_seed_parallel[];
//long int rand_seed_parallel[OPENMP_MAX_NPROC];
DOUBLE ran2parallel(long *idum);
#endif


// define Random()
#ifdef _OPENMP
#  define Random() ran2parallel(&rand_seed_parallel[omp_get_thread_num()])
//#  define Random() RandomCritical()
//#  define Random() ran2(&rand_seed)
#else // not defined _OPENMP
#ifdef USE_RAN2_GENERATOR
#  define Random() ran2(&rand_seed)
  DOUBLE ran2(long *idum);
#else
//#  define Random() RandomSys()
  double genrand(void);
#  define Random() genrand()
#endif
#endif

#endif

