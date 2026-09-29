/*utils.c*/

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include "compatab.h"
#include "utils.h"
#include "main.h"
#include "vmc.h"
#include "memory.h"
#include MATHINCLUDE

#ifdef MPI
#  include "mpi.h"
#  include "parallel.h"
#endif

/************************************ Error ********************************/
void Error(const char * format, ...) {
  va_list arg;
  FILE *out;
  char text[100];

#ifndef MPI
  out = fopen(OUTPATH "MC.log","a");
#else
  char num[2]=" ";
  num[0] = '0'+myid;
  strcpy(text, OUTPATH "MC");
  strcat(text, num);
  strcat(text, ".log");
  out = fopen(text, "a");
#endif
  
  va_start(arg, format);
  vsprintf(text, format, arg);
  va_end(arg);

  fprintf(stderr,"\nFATAL ERROR:\n");
  fprintf(stderr, "%s\n", text);

  fprintf(out, "\nFATAL ERROR:\n %s", text);
  fclose(out);

#ifdef MPI
// MPI_Finalize();
#endif

  exit(0);
}

/************************************** Warning ****************************/
void Warning(const char *format, ...) {
  va_list arg;
  static int CR = ON, x=20;
  FILE *out;
  char text[200];

#ifndef MPI
  out = fopen(OUTPATH "MC.log", "a");
#else
  char num[2]=" ";
  num[0] = '0'+myid;
  strcpy(text, OUTPATH "MC");
  strcat(text, num);
  strcat(text, ".log");
  out = fopen(text, "a");
#endif

  va_start(arg, format);
  vsprintf(text, format, arg);

  if(video == OFF) {
    printf("WARNING: %s", text);
  }
  else{ 
    if(CR) {
      setcolor(BLACK);
      bar(0,20, 1200, 40);
    }
    setcolor(LIGHTRED);
    outtextxy(x,21,text);
    if(text[strlen(text)-1] != '\n') {
      CR = OFF;
      x += strlen(text)*9;
    }
    else {
      CR = ON;
      x = 20;
    }
  }

  if(out) {
    fprintf(out, "WARNING: %s", text);
    fclose(out);
  }
}

/************************************** Warning ****************************/
void Message(const char * format, ...) {
  va_list arg;
  static int initialized = OFF;
  static int CR = ON, x=20;
  FILE *out;
  char text[200];

#ifndef MPI
  out = fopen(OUTPATH "MC.log", initialized++?"a":"w");
#else
  char num[2]=" ";
  num[0] = '0'+myid;
  strcpy(text, OUTPATH "MC");
  strcat(text, num);
  strcat(text, ".log");
  out = fopen(text, initialized++?"a":"w");
#endif

  va_start(arg, format);
  vsprintf(text, format, arg);

  if(video == OFF) {
    printf("%s", text);
  }
  else{
    if(CR) {
      setcolor(BLACK);
      bar(0,0, 1200, 20);
    }
    setcolor(WHITE);
    outtextxy(x,1,text);
    if(text[strlen(text)-1] != '\n') {
      CR = OFF;
      x += strlen(text)*9;
    }
    else {
      CR = ON;
      x = 20;
    }
  }

  if(out) {
    fprintf(out, "%s", text);
    fclose(out);
  }
}

/*************************************** Exit *******************************/
void Exit(int status, const char * format, ...) {

  va_list arg;

  Warning(format);
  fprintf(stderr,"\nExit message: ");
  va_start(arg, format);
  fprintf(stderr, format, arg);
  va_end(arg);

  exit(status);
}

/********************************** Check Mantissa **************************/
/*

mantissa
        linux     windows/debug  windows/release
        gcc icc   icc visual     icc visual
float   6   18    6   6          14  14
double  14  18    14  14         14  14
long    18  18    14  14         14  14 */
void CheckMantissa(void) {
  int i;
  float a1 = 1.1f;
  double a2 = 1.1;
  long double a3 = 1.1;
  DOUBLE a4 = 1.1;
  int l1 = 0;
  int l2 = 0;
  int l3 = 0;
  int l4 = 0;
  double d1,d2,d3,d4;
  for(i=0; i<20; i++) {
    d1 = (a1-1.)*pow(10.,(double)(i+1));
    d2 = (a2-1.)*pow(10.,(double)(i+1));
    d3 = (a3-1.)*pow(10.,(double)(i+1));
    d4 = (a4-1.)*pow(10.,(double)(i+1));
    //printf("%i %" LF " %i \n", i, d1, (d1>0.)?(1):(0));
    a1 = 0.1f*(a1-1.f) + 1.f;
    a2 = 0.1*(a2-1.) + 1.;
    a3 = 0.1*(a3-1.) + 1.;
    a4 = 0.1*(a4-1.) + 1.;
    if(d1>0) l1 = i+1;
    if(d2>0) l2 = i+1;
    if(d3>0) l3 = i+1;
    if(d4>0) l4 = i+1;
  }
  Message("  Checking datatypes, number of digits in mantissa ...\n");
  Message("    float       %i \n", l1);
  Message("    double      %i \n", l2);
  Message("    long double %i \n", l3);
  Message("    DOUBLE      %i \n\n", l4);
}

/******************* Check Memory Copy **************************************/
int CheckMemcpy(void) {
  int i;
  int t, repeat_times=1000000;
  clock_t time_memcpy;
  DOUBLE time1, time2;
  int check_passed = ON;
  int index_copy, memory_size;

  if(Nwalkers == NwalkersMax) return -1; // comparison cannot be done

  Message("  Checking memory copy speed ... \n");
  // Copy first walker to the last walker in the memory using different methods
  index_copy = NwalkersMax-1;

  CopyWalker(&W[index_copy], &W[0]);

  time_memcpy = clock();
  for(t=0; t<repeat_times; t++) {
    for(i=0; i<N; i++) {
      W[index_copy].x[i] = W[0].x[i];
      W[index_copy].y[i] = W[0].y[i];
      W[index_copy].z[i] = W[0].z[i];
    }
  }
  time1 = (DOUBLE)(clock()-time_memcpy)/(DOUBLE)CLOCKS_PER_SEC;
  Message("    loop copy: %e sec\n", time1);
  if(WalkerCompare(W[index_copy], W[0])) check_passed = OFF;

  for(i=0; i<N; i++) {
    W[index_copy].x[i] = W[0].x[i]+1;
    W[index_copy].y[i] = W[0].y[i]+1;
    W[index_copy].z[i] = W[0].z[i]+1;
  }

  //memory_size = sizeof(DOUBLE)*N;
  memory_size = sizeof(*W[0].x)*N;
  time_memcpy = clock();
  for(t=0; t<repeat_times; t++) {
    //memcpy(W[index_copy].x, W[0].x, sizeof(DOUBLE)*3*N);
    memcpy(W[index_copy].x, W[0].x, memory_size);
    memcpy(W[index_copy].y, W[0].y, memory_size);
    memcpy(W[index_copy].z, W[0].z, memory_size);
  }
  time2 = (DOUBLE)(clock()-time_memcpy)/(DOUBLE)CLOCKS_PER_SEC;
  if(WalkerCompare(W[index_copy], W[0])) check_passed = OFF;
  if(check_passed == OFF) Error("  walkers compare: check not passed.");
  Message("    memcopy: %e sec\n", time2);
  Message("    speed gain: %g times\n", time1/time2);

  /*time_memcpy = clock();
  for(t=0; t<repeat_times; t++) {
    W[index_copy] = W[0];
  }
  time2 = (DOUBLE)(clock()-time_memcpy)/(DOUBLE)CLOCKS_PER_SEC;
  Message("    struct copy: %e sec\n", time2);
  Message("    speed gain: %g times\n", time1/time2);*/

  Message("  done\n");

  return check_passed;
}

/******************* Check Contigious Array *********************************/
int CheckContigiousArray(int **a, int dim1, int dim2) {
// checks order in memory of elements
  int i,j;

  for(i=0; i<dim1; i++) {
    for(j=0; j<dim2; j++) {
    //Message("%i %i a[i][j]=%i &a[i][j]=%i &a[i][j]-&a[i]=%i &a[i][j]-&a=%i\n", i, j, a[i][j], &a[i][j], &a[i][j]-a[i], &a[i][j]-&a[0][0]);
      if(&a[i][j] - &a[0][0] != i*dim2+j) {
        Warning("  memory allocation of this array is not contiguous!\n");
        return 1;
      }
    }
  }
  return 0;
}

/***************************** Save Matrix **********************************/
void SaveMatrix(char *fileout, DOUBLE **matrix, int size1, int size2) {
  FILE *out;
  int i,j;

  out = fopen(fileout, "w");

  for(i=0; i<size1; i++) {
    for(j=0; j<size2; j++) {
      fprintf(out, "%.15e ", matrix[i][j]);
    }
    fprintf(out, "\n");
  }

  fclose(out);
}

void SaveGrid(char *fileout, struct Grid *G) {
  FILE *out;
  int i;

  out = fopen(fileout, "w");

  for(i=0; i<G->size; i++) fprintf(out, "%.15e %.15e\n", G->x[i], G->f[i]);

  fclose(out);
}

/******************* Message Time Elapsed ***********************************/
void MessageTimeElapsed(int sec) {
  int day, min, hour;

  if(sec<60) {
    Message("%2dsec ", sec);
  }
  else if(sec<3600) {
    min = sec/60;
    sec = sec%60;
    Message("%2d:%2d ", min, sec);
  }
  else if(sec<86400) {
    min = sec/60;
    sec = sec%60;
    hour = min/60;
    min = min%60;
    Message("%2d:%2d:%2d ", hour, min, sec);
  }
  else {
    min = sec/60;
    sec = sec%60;
    hour = min/60;
    day = hour / 24;
    min = min%60;
    hour = hour%24;
    Message("%2dd %2dh %2dm %2ds ", day, hour, min, sec);
  }
}

/****************** Check Energy And Derivatives ****************************/
int CheckEnergyAndDerivatives(void) {
  int w;
  DOUBLE Ekin;
  int Nmax = 10;

  if(Nmax > Nwalkers) Nmax = Nwalkers;

#ifdef CALCULATE_DERIVATIVES_NUMERICALLY
  return 0; // cannot check
#endif

    Message("\n  Checking calculation of the kinetic energy by discrete differences\n");

#pragma omp parallel for
    for(w=0; w<Nmax; w++) {
      W[w].E = WalkerEnergy(&W[w], &W[w].Epot, &W[w].Ekin, &W[w].EFF, &W[w].Edamping, &W[w].Eint, &W[w].Eext);
      Ekin = MeasureKineticEnergyNumerically(&W[w]);
      if(fabs(W[w].Ekin)<1e-100)
        Warning("  extremely small (zero?) value of the kinetic energy is found\n");
      else
        Message("    %i Analytic %.3"LE"  discrete %.3"LE" Err: rel %.2"LF" %%  abs %.2"LE"\n", 
        w, W[w].Ekin/(DOUBLE)N, Ekin/(DOUBLE)N, fabs((W[w].Ekin-Ekin)/W[w].Ekin)*100., fabs(W[w].Ekin-Ekin)/(DOUBLE)N);
    }
    Message("  done\n");

    return 0;
}