/*gencoord.h*/
#ifndef _GENERATE_COORDINATES_H_
#define _GENERATE_COORDINATES_H_

#define ON 1
#define OFF 0
#define MAX 1000
#define MAX_GEN_ITER 100000

int GenerateCoordinates(void);
void GenerateOneParticleCoordinates(DOUBLE *x, DOUBLE *y, DOUBLE *z, int i);
void GenSaveParticles(void);

#endif
