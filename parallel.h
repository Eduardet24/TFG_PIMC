/*parallel.h*/

#ifndef _PARALLEL_H_
#define _PARALLEL_H_

void ParallelRandInit(void);

#ifdef MPI

#include <mpi.h>

#define MASTER 0
extern int myid; // check if this runs is a mastere or a slave
extern int Nslaves;

#define MSG_NWALKER 1
#define MSG_WALKER 2
#define MSG_REPLY 3
#define MSG_RADIAL_DISTRIBUTION 4
#define MSG_OBDM 5

extern int NwalkersAll; // total number of walkers on all machines
extern char *buffer;
extern int message_size;
extern int *NwalkersSl; // number of walkers at a given machine
extern int *walker_indices; // original indices of walkers distributed to a given machine
extern double time_done_in_parallel;

extern MPI_Status mpi_status; // error status

void ParallelInitialize(int argn, char **argv);
void ParallelInitialize2(void);
void ParallelStop(void);
void CollectAllWalkers(void);
void DistributeAllWalkers(void);

void ParallelSaveEnergyVMC(void);
void ParallelSaveEnergyDMC(DOUBLE E);
void ParallelSaveMomentumDistribution(void);

void CollectAllWalkersData(void);
int SendWalkers(int dest, int wmin, int Nwalk, int tag, MPI_Request *request);
void ReceiveWalkersData(int src, int wmin, int tag);
void ReceiveWalkersUnpack(int src, int wmin);

void ReceiveWalkers(int src, int wmin, int tag, MPI_Request *request);
int SendWalkersData(int dest, int wmin, int Nwalk, int tag);

void BranchingMPI(void);

#endif

/**************************** Open MP ****************************************/
#ifdef _OPENMP

#include <omp.h>
#include "main.h"

#define MASTER 0

extern int myid; // check if this runs is a mastere or a slave
extern int Nslaves;

void ParallelInitialize(int argn, char **argv);
DOUBLE RandomCritical(void);

#endif

#endif


