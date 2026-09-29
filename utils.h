/*utils.h*/

#ifndef _ERROR_H
#define _ERROR_H

#include "main.h" //define DOUBLE

void Error(const char * format, ...);
void Warning(const char * format, ...);
void Message(const char * format, ...);
void Exit(int status, const char * format, ...);
void CheckMantissa(void);
int CheckMemcpy(void);
int CheckContigiousArray(int **a, int dim1, int dim2);
int CheckEnergyAndDerivatives(void);
void SaveMatrix(char *fileout, DOUBLE **matrix, int size1, int size2);
void SaveGrid(char *fileout, struct Grid *G);
void MessageTimeElapsed(int sec);

#define ON 1
#define OFF 0

#endif
