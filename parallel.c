/*parallel.c*/
#include "main.h"

#ifdef MPI

#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "parallel.h"
#include "utils.h"
#include "memory.h"
#include "vmc.h"
#include "rw.h"
#include "randnorm.h"

int myid;    // process identification
int Nslaves; // number of jobs
MPI_Status mpi_status; // error status
char *buffer; // buffer for collecting global data
char **bufferSl; // array of buffers for non blocking send
char *buffer_local; // buffer for local data
int *buffer_int;
double *buffer_double;
int *buffer_displacement; // needed for SCATTER
int *buffer_elements; // needed for SCATTER
int buffer_size; // size of all packed walker in memory
int message_size; // size of the local packet
int measurements_in_block; // number of measurements of energy in one block
MPI_Request *requestSl;
MPI_Status *mpi_status_Sl;

double time_parallel_total = 0.;    // time of the calculation
double time_parallel_transfer = 0.; // time lost for transferring data
double time_parallel_transfer_pack = 0.; // lost due to packing / unpacking
double time_parallel_transfer_send = 0.;  //due to collect / gather
double time_done_in_parallel = 0.;   // time of distributed calculation

int NwalkersAll; // total number of walkers, now Nwalkers corresponds to a given machine
int *walker_indices, *NwalkersSl; // number of walkers on each of the processors

// storing measurements, i.e saves measurements in the block
double *storeE, *storeEFF, *storeNwalkersw;
int *storeNwalkers;
int storeEN; // number of stored measurements
double *bufferEall; // buffer for MPI
double *bufferE, *bufferEFF, *bufferNwalkersw;
int *bufferNwalkers;

/******************************* Parallel Initialize *************************/
// MPI initialization, number of particles is not yet known
void ParallelInitialize(int argn, char **argv) {
  int namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];

  Message("\nStarting MPI parallelization ...\n");

  MPI_Init(&argn, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &Nslaves);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  Message("parameters are %s\n",argv[0]);
  Message("total number of processors is %i\n", Nslaves);

  time_parallel_total = MPI_Wtime();
  Message("  time count started.\n");

  MPI_Get_processor_name(processor_name, &namelen);
  Message("  Machine name %s, process id number %i\n", processor_name, myid);

  Message("done\n");
}

/******************************* Parallel Initialize 2 ***********************/
// now number of particles is known, allocate buffers
void ParallelInitialize2(void) {
  int i, w;

  // Allocate arrays
  Message("Parallel Initialize2: Allocating memory for the buffer\n");
  walker_indices = (int*) Calloc("walker_indices", Nslaves, sizeof(int));
  NwalkersSl = (int*) Calloc("NwalkersSl", Nslaves, sizeof(int));
  requestSl =  (MPI_Request*) Calloc("requestSl", Nslaves, sizeof(MPI_Request));
  mpi_status_Sl =  (MPI_Status*) Calloc("mpi_status_Sl", Nslaves, sizeof(MPI_Status));
  NwalkersAll = Nwalkers;

  //Allocate pack buffer 
  //  Check number of passed elements!!!
  buffer_size = NwalkersMax*(3*N*sizeof(double)+3*sizeof(int)); // walkers coordinates
  buffer_size += Niter*2; //E
  if(measure_OBDM) buffer_size += NwalkersMax*2*OBDM.size*sizeof(double);
  if(measure_PairDistr) buffer_size += NwalkersMax*((PD.size)*sizeof(double)+sizeof(int));
  if(measure_PairDistr && MC == DIFFUSION) buffer_size += NwalkersMax*((PD.size)*sizeof(double)+sizeof(int));
  if(measure_SD && MC == DIFFUSION) {
    buffer_size += NwalkersMax*3*SD.size*sizeof(double); // W[].CM
    buffer_size += NwalkersMax*N*3*sizeof(double); // W[].rreal
    buffer_size += NwalkersMax*SD.size*sizeof(double); // SD.CM2
    buffer_size += NwalkersMax*SD.size*sizeof(int); // SD.N
    buffer_size += NwalkersMax*N*3*sizeof(double); // rreal
  }
  buffer_size += 10; // additional 10 spaces are added

  if(measure_energy) { // allocate arrays for storing energy
    storeEN = 0;
    measurements_in_block = Niter/Nmeasure;
    storeE = (double*) Calloc("storeE", measurements_in_block, sizeof(double)); // store all values of E in one block
    storeEFF = (double*) Calloc("storeEFF", measurements_in_block, sizeof(double)); // store all values of EFF in one block

    bufferEall = (double*) Calloc("bufferE", NwalkersMax*measurements_in_block, sizeof(double)); // store E Nwalkers
    bufferE = (double*) Calloc("bufferE", measurements_in_block, sizeof(double));
    bufferEFF = (double*) Calloc("bufferEFF", measurements_in_block, sizeof(double));
    buffer_size += 2*NwalkersMax*measurements_in_block*sizeof(double); // E, EFF
    if(MC == DIFFUSION) {
      storeNwalkersw = (double*) Calloc("storeNwalkersw", measurements_in_block, sizeof(double));
      storeNwalkers = (int*) Calloc("storeEFF", measurements_in_block, sizeof(int));
      bufferNwalkersw = (double*) Calloc("bufferNwalkersw", measurements_in_block, sizeof(double));
      bufferNwalkers = (int*) Calloc("bufferNwalkers", measurements_in_block, sizeof(int));
      buffer_size += NwalkersMax*measurements_in_block*(sizeof(double)+sizeof(int)); // Nwalkers, Nwalkersw
    }
  }

  Message("  MPI: buffer size is %i\n", buffer_size);
  buffer = (char*) Calloc("buffer", buffer_size, sizeof(char)); // buffer used for SCATTER and GATHER

  buffer_local = (char*) Calloc("buffer_local", buffer_size/Nslaves+10, sizeof(char));  // buffer used locally on each node
  buffer_int = (int*) Calloc("buffer_int", buffer_size/Nslaves+10, sizeof(int));
  buffer_double = (double*) Calloc("buffer_int", buffer_size/Nslaves+10, sizeof(double));
  buffer_displacement = (int*) Calloc("buffer_displacement", Nslaves, sizeof(int)); // needed for SCATTER
  buffer_elements = (int*) Calloc("buffer_elements", Nslaves, sizeof(int));     // needed for SCATTER

  bufferSl = (char**) Calloc("buffer_local", Nslaves, sizeof(char*));
  for(i=0; i<Nslaves; i++) {
    bufferSl[i] = (char*) Calloc("buffer_local", buffer_size/Nslaves+10, sizeof(char));
  }

  // Check the network
  Message("\nChecking that network spreading works correctly ... \n");
  for(w=0; w<NwalkersAll; w++) {

    /*SaveCoordinates("in3Dprev.in");
    for(w=0; w<Nwalkers; w++) {
      for(i=0; i<N; i++) {
        W[w].x[i] = w*100+i;
        W[w].y[i] = w*100+i+0.01;
        W[w].z[i] = w*100+i+0.02;
      }
    }
    SaveCoordinates("in3Dprev0.in");
    */
    if(myid == MASTER) CopyWalker(&Wp[w], &W[w]);//w
    if(verbosity) Message(" distribution of walkers started\n");

    DistributeAllWalkers();

    if(verbosity) Message(" distribution of walkers finished\n");

    if(myid == MASTER && w>NwalkersSl[0]) {
      for(i=0; i<N; i++) W[w].x[i] = W[w].y[i] = W[w].z[i] = 0.;
      W[w].E  = W[w].U = 0.;
    }
    //Message("! walker %i, Nwalkers %i, NwalkersSl[0] %i\n", w, Nwalkers, NwalkersSl[0]);//!!!

    CollectAllWalkers();

    if(myid == MASTER) {
      if(WalkerCompare(wp, W[w])) {
        Error(" failed for walker W[%i] !\n", w);
      }
    }
  }
  Message("done\n\n");
}

/**************************** ParallelRandInit ********************************/
void ParallelRandInit(void) {
#ifndef USE_RAN2_GENERATOR
  Error("define USE_RAN2_GENERATOR");
#endif
  rand_seed = (long int)(pow(2.,myid+3.)-1.);

  Message(" Ran2() random number generator started with seed %i\n", rand_seed);

  if(RandomSys()<0) Error("bad implementation of RandomSys(). It should not produce negative numbers");
  if(RandomSys()>1) Error("bad implementation of RandomSys(). It should not produce numbers larger than 1");
}

/******************************* Parallel stop *******************************/
void ParallelStop(void) {
  double effective_single;

  //time_done_in_parallel += MPI_Wtime();
  time_parallel_total = MPI_Wtime() - time_parallel_total;

  Message("\nMPI time statistics, time:\n");
  Message("    total:   : %.1f sec (100%%)\n", time_parallel_total);
  Message("    serial   : %.1f sec (%2.lf%%)\n", time_parallel_total-time_done_in_parallel, (time_parallel_total-time_done_in_parallel)/time_parallel_total*100);
  Message("    parallel : %.1f sec (%2.lf%%), effective %.1f sec of a single CPU\n", time_done_in_parallel, time_done_in_parallel/time_parallel_total*100, time_done_in_parallel*Nslaves);
  Message("    transfer : %.1f sec (%2.lf%%)\n", time_parallel_transfer, time_parallel_transfer/time_parallel_total*100);
  Message("    transfer (pack) : %.1f sec (%2.lf%%)\n", time_parallel_transfer_pack, time_parallel_transfer_pack/time_parallel_total*100);
  Message("    transfer (send) : %.1f sec (%2.lf%%)\n", time_parallel_transfer_send, time_parallel_transfer_send/time_parallel_total*100);
  effective_single = time_parallel_total - time_parallel_transfer + (Nslaves-1)*time_done_in_parallel;
  Message("    speed gain : %.1f times\n", effective_single/time_parallel_total);

  MPI_Finalize();
}

/**************************** Distribute All Walkers ****************************/
void DistributeAllWalkers(void) {
  int w,i;
  static int tag = 0;

  tag++;

  if(myid == MASTER) { // the master code
    if(verbosity) Message("Distribution of %i walkers to %i nodes\n", Nwalkers, Nslaves);
    NwalkersAll = Nwalkers;
    NwalkersSl[0] = (int)ceil((double)NwalkersAll / (double) Nslaves); // distribute walkers in a homogeneous way
    walker_indices[0] = 0; // initialize

    for(i=1; i<Nslaves; i++) {
      walker_indices[i] = walker_indices[i-1] + NwalkersSl[i-1]; // first index of node i
      if(i<Nslaves-1)
        NwalkersSl[i] = (int)ceil((double)NwalkersAll / (double) Nslaves); // number of walkers on node i
      else
        NwalkersSl[i] = NwalkersAll - walker_indices[i]; // last node takes the rest
    }
    //walker_indices[Nslaves] = NwalkersAll;// outside of array???

    if(verbosity) {
      Message("\nSending configuration coordinates to slaves...\n");
      Message("  => MASTER has   %i walkers :", NwalkersSl[0]);
      Message("  [%i,%i]\n", walker_indices[0] + 1, walker_indices[0]+NwalkersSl[0]);
    }
    Nwalkers = NwalkersSl[0];

    for(i=1; i<Nslaves; i++) {
      if(verbosity) Message("  => SL %i sending %i walkers :  [%i,%i]\n", i, NwalkersSl[i], walker_indices[i]+1, walker_indices[i]+NwalkersSl[i]);
      SendWalkers(i, walker_indices[i], NwalkersSl[i], tag, &requestSl[i]);
      if(verbosity) Message("  => SL %i sending %i walkers :  [%i,%i] (done)\n", i, NwalkersSl[i], walker_indices[i]+1, walker_indices[i]+NwalkersSl[i]);
    }
    MPI_Waitall(Nslaves-1, &requestSl[1], mpi_status_Sl); // wait until non-blocking send is completed

    if(verbosity) Message("done\n");
  }
  else { // the slave code
    if(verbosity) Message("Receiving walkers ...\n");
    ReceiveWalkers(MASTER, 0, tag, &requestSl[myid]);
    MPI_Wait(&requestSl[myid], &mpi_status);
    ReceiveWalkersUnpack(MASTER, 0);
    for(w=0; w<Nwalkers; w++) W[w].w = w; // update local index
    if(verbosity) Message("done\n");
  }

  if(MC == DIFFUSION) {
    for(w=0; w<Nwalkers; w++) W[w].status = ALIVE;
    for(w=Nwalkers; w<NwalkersAll; w++) W[w].status = DEAD;
  }
}

/**************************** Collect All Walkers ****************************/
void CollectAllWalkers(void) {
  int i;
  static int tag=1000;

  tag++;
  if(verbosity) Message("Collect all walkers: packing local %i walkers ... ", Nwalkers);

  if(myid == MASTER) {
    for(i=1; i<Nslaves; i++) ReceiveWalkers(i, walker_indices[i], tag, &requestSl[i]);
    MPI_Waitall(Nslaves-1, &requestSl[1], mpi_status_Sl); // wait until non-blocking receive is completed
    for(i=1; i<Nslaves; i++) ReceiveWalkersUnpack(i, walker_indices[i]);
    Nwalkers = NwalkersAll;
  }
  else {
    SendWalkers(MASTER, 0, Nwalkers, tag, &requestSl[myid]);
  }
}

/***************************** Send Walkers ***********************************/
// send walkers W[wmin],...,W[wmin+Nwalk] to machine No dest
int SendWalkers(int dest, int wmin, int Nwalk, int tag, MPI_Request *request) {
  int w;
  int buffer_position = 0;

  if(verbosity) Message("\npacking walkers dest %i, wmin %i, Nwalk %i, tag %i...", dest, wmin, Nwalk, tag);
#ifdef SECURE
  if(wmin<0) Error("SendWalkers called with negative w\n");
#endif

  time_parallel_transfer -= MPI_Wtime();
  time_parallel_transfer_pack -= MPI_Wtime();
  MPI_Pack(&Nwalk, 1, MPI_INT, bufferSl[dest], buffer_size, &buffer_position, MPI_COMM_WORLD);
  for(w=wmin; w<wmin+Nwalk; w++) {
    MPI_Pack(W[w].x, N, MPI_DOUBLE, bufferSl[dest], buffer_size, &buffer_position, MPI_COMM_WORLD);
    MPI_Pack(W[w].y, N, MPI_DOUBLE, bufferSl[dest], buffer_size, &buffer_position, MPI_COMM_WORLD);
    MPI_Pack(W[w].z, N, MPI_DOUBLE, bufferSl[dest], buffer_size, &buffer_position, MPI_COMM_WORLD);
    MPI_Pack(&W[w].U, 1, MPI_DOUBLE, bufferSl[dest], buffer_size, &buffer_position, MPI_COMM_WORLD);
    MPI_Pack(&W[w].E, 1, MPI_DOUBLE, bufferSl[dest], buffer_size, &buffer_position, MPI_COMM_WORLD);
    if(measure_SD && MC == DIFFUSION) {
      MPI_Pack(W[w].CM[0], 3*SD.size, MPI_DOUBLE, bufferSl[dest], buffer_size, &buffer_position, MPI_COMM_WORLD);
      MPI_Pack(W[w].rreal[0], N*3, MPI_DOUBLE, bufferSl[dest], buffer_size, &buffer_position, MPI_COMM_WORLD);
    }
  }
  time_parallel_transfer_pack += MPI_Wtime();

  if(verbosity) Message(" sending walkers ...");

  time_parallel_transfer_send -= MPI_Wtime();

  //int MPI_Send(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
  //int MPI_Isend(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request) 

  //MPI_Send(bufferSl[dest], buffer_position, MPI_PACKED, dest, tag, MPI_COMM_WORLD);
  MPI_Isend(bufferSl[dest], buffer_position, MPI_PACKED, dest, tag, MPI_COMM_WORLD, request);
  //MPI_Wait(request,&mpi_status);
  time_parallel_transfer_send += MPI_Wtime();
  time_parallel_transfer += MPI_Wtime();

  if(verbosity) Message(" done\n");
  return buffer_position;
}

/**************************** Receive Walkers *********************************/
void ReceiveWalkers(int src, int wmin, int tag, MPI_Request *request) {
  if(verbosity) Message("receiving walkers ... ");
  time_parallel_transfer -= MPI_Wtime();
  time_parallel_transfer_send -= MPI_Wtime(); // store moment when the function was called

  //int MPI_Recv( void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status)
  //int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request)
  //MPI_Recv(buffer, buffer_size, MPI_PACKED, src, tag, MPI_COMM_WORLD, &mpi_status);
  MPI_Irecv(bufferSl[src], buffer_size, MPI_PACKED, src, tag, MPI_COMM_WORLD, request);
  time_parallel_transfer_send += MPI_Wtime();
  time_parallel_transfer += MPI_Wtime();
  if(verbosity) Message("done\n");
}

void ReceiveWalkersUnpack(int src, int wmin) {
  int buffer_position = 0;
  int nwalkers; // nwalkers to be received
  int w;

  time_parallel_transfer -= MPI_Wtime();
  time_parallel_transfer_pack -= MPI_Wtime();
  MPI_Unpack(bufferSl[src], buffer_size, &buffer_position, &nwalkers, 1, MPI_INT, MPI_COMM_WORLD);
  if(verbosity) Message("unpacking %i walkers...", nwalkers);
  for(w=wmin; w<wmin+nwalkers; w++) {
    MPI_Unpack(bufferSl[src], buffer_size, &buffer_position, W[w].x, N, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Unpack(bufferSl[src], buffer_size, &buffer_position, W[w].y, N, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Unpack(bufferSl[src], buffer_size, &buffer_position, W[w].z, N, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Unpack(bufferSl[src], buffer_size, &buffer_position, &W[w].U, 1, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Unpack(bufferSl[src], buffer_size, &buffer_position, &W[w].E, 1, MPI_DOUBLE, MPI_COMM_WORLD);
    W[w].w = w; // update local index
    W[w].status = ALIVE;
    if(measure_SD && MC == DIFFUSION) {
      MPI_Unpack(bufferSl[src], buffer_size, &buffer_position, W[w].CM[0], 3*SD.size, MPI_DOUBLE, MPI_COMM_WORLD);
      MPI_Unpack(bufferSl[src], buffer_size, &buffer_position, W[w].rreal[0], N*3, MPI_DOUBLE, MPI_COMM_WORLD);
    }
  }
  time_parallel_transfer_pack += MPI_Wtime();
  time_parallel_transfer += MPI_Wtime();

  if(wmin == 0) Nwalkers = nwalkers; // update Nwalkers only if filled from zero

  if(verbosity) Message("done\n");
}

/******************************** Save Energy VMC ************************/
void ParallelSaveEnergyVMC(void) {
  static int file_open = OFF;
  int i;
  int buffer_position = 0;
  FILE *out;

  // accumulate E in array
  storeE[storeEN] = E;
  storeEFF[storeEN] = EFF;

  if(storeEN == measurements_in_block) {
    storeEN = -1;
    // transfer data
    time_parallel_transfer -= MPI_Wtime(); // store moment when the function was called
    time_parallel_transfer_send -= MPI_Wtime();

//   int MPI_Gather(void *sendbuf, int sendcnt, MPI_Datatype sendtype, 
//                 void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm )
//       sendbuf 	starting address of send buffer (choice)
//       sendcount 	number of elements in send buffer (integer)
//       sendtype 	data type of send buffer elements (handle)
//       recvbuf    address of receive buffer (choice, significant only at root) 
//       recvcount 	number of elements for any single receive (integer, significant only at root)
//       recvtype 	data type of recv buffer elements (significant only at root) (handle)
//       root 	rank of receiving process (integer)
//       comm 	communicator (handle) 
    //MPI_Gather(&E, 1, MPI_DOUBLE, bufferEall, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); // one number
    //MPI_Gather(storeE, measurements_in_block, MPI_DOUBLE, bufferEall, measurements_in_block, MPI_DOUBLE, 0, MPI_COMM_WORLD); // array of double

    //void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm)
    MPI_Reduce(storeE, bufferE, measurements_in_block, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(storeEFF, bufferEFF, measurements_in_block, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

//     int MPI_Pack(void *inbuf, int incount, MPI_Datatype datatype, void *outbuf, int outcount, int *position, MPI_Comm comm )
//       inbuf 	input buffer start (choice)
//       incount 	number of input data items (integer)
//       datatype 	datatype of each input data item (handle)
//       outbuf -output buffer start (choice)
//       outcount 	output buffer size, in bytes (integer)
//       position 	current position in buffer, in bytes (integer)
//       comm 	communicator for packed message (handle)

    //MPI_Pack(storeE, measurements_in_block, MPI_DOUBLE, buffer_local, buffer_size, &buffer_position, MPI_COMM_WORLD);
    //MPI_Gather(buffer_local, buffer_position, MPI_PACKED, buffer, buffer_position, MPI_PACKED, 0, MPI_COMM_WORLD); // array of doubles

    //MPI_Unpack(buffer_local, message_size, &buffer_position, storeE, measurements_in_block, MPI_DOUBLE, MPI_COMM_WORLD);
    //MPI_Gather(storeE, buffer_position, MPI_PACKED, bufferEall, buffer_position, MPI_PACKED, 0, MPI_COMM_WORLD);

    time_parallel_transfer_send += MPI_Wtime();

    if(myid == MASTER) {
      time_parallel_transfer_pack -= MPI_Wtime();
      out = fopen(file_energy, file_open?"a":"w");
      if(file_open == OFF) file_open = ON; 
      if(out == NULL) {
        perror("\nError:");
        Warning("can't write to energy file %s\n", file_energy);
        return;
      }
      time_parallel_transfer_pack += MPI_Wtime();
      for(i=0; i<measurements_in_block; i++) {
        fprintf(out, "%.15e %.15e\n", bufferE[i]/(double)Nslaves, bufferEFF[i]/(double)Nslaves);
      }
      fclose(out);
    }
    time_parallel_transfer += MPI_Wtime(); //update transfer time*/
  }

  storeEN++;
}

/******************************** Save Energy VMC ************************/
void ParallelSaveEnergyDMC(DOUBLE E) {
  static int file_open = OFF;
  int i;
  int buffer_position = 0;
  FILE *out;

  // accumulate E in array
  storeE[storeEN] = E;
  storeEFF[storeEN] = EFF;
  storeNwalkers[storeEN] = Nwalkers;
  storeNwalkersw[storeEN] = Nwalkersw;

  if(storeEN == measurements_in_block-1) {
    storeEN = -1;
    // transfer data
    time_parallel_transfer -= MPI_Wtime(); // store moment when the function was called
    time_parallel_transfer_send -= MPI_Wtime();
    MPI_Reduce(storeE, bufferE, measurements_in_block, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(storeEFF, bufferEFF, measurements_in_block, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(storeNwalkers, bufferNwalkers, measurements_in_block, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(storeNwalkersw, bufferNwalkersw, measurements_in_block, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    time_parallel_transfer_send += MPI_Wtime();
    if(myid == MASTER) {
      time_parallel_transfer_pack -= MPI_Wtime();
      out = fopen(file_energy, file_open?"a":"w");
      if(file_open == OFF) file_open = ON; 
      if(out == NULL) {
        perror("\nError:");
        Warning("can't write to energy file %s\n", file_energy);
        return;
      }
      time_parallel_transfer_pack += MPI_Wtime();
      for(i=0; i<measurements_in_block; i++) {
        fprintf(out, "%.15" LE " %.15" LE " %i %.2" LF "\n", bufferE[i]/(double)Nslaves, bufferEFF[i]/(double)Nslaves, bufferNwalkers[i], bufferNwalkersw[i]);
      }
      fclose(out);
    }
    time_parallel_transfer += MPI_Wtime(); //update transfer time
  }

  storeEN++;
}

/**************************** Collect All Walkers Data ***********************/
void CollectAllWalkersData(void) {
  int i;
  static int tag = 2000;

  time_parallel_transfer -= MPI_Wtime(); // store moment when the function was called

  if(myid == MASTER) {
    for(i=1; i<Nslaves; i++) {
      ReceiveWalkersData(i, walker_indices[i], tag);
    }
  }
  else {
    SendWalkersData(MASTER, 0, Nwalkers, tag);
  }
  time_parallel_transfer += MPI_Wtime(); //update transfer time
}

/***************************** Send Walkers Data ******************************/
// send walkers W[wmin],...,W[wmin+Nwalk] and arrays to machine No dest
int SendWalkersData(int dest, int wmin, int Nwalk, int tag) {
  //int w;
  int buffer_position = 0;

  if(verbosity) Message("Send walker data..");
  time_parallel_transfer -= MPI_Wtime();
  time_parallel_transfer_pack -= MPI_Wtime();
//   MPI_Pack(&Nwalk, 1, MPI_INT, buffer, buffer_size, &buffer_position, MPI_COMM_WORLD);
//   for(w=wmin; w<wmin+Nwalk; w++) {
//     MPI_Pack(W[w].x, N, MPI_DOUBLE, buffer, buffer_size, &buffer_position, MPI_COMM_WORLD);
//     MPI_Pack(W[w].y, N, MPI_DOUBLE, buffer, buffer_size, &buffer_position, MPI_COMM_WORLD);
//     MPI_Pack(W[w].z, N, MPI_DOUBLE, buffer, buffer_size, &buffer_position, MPI_COMM_WORLD);
//     MPI_Pack(&W[w].U, 1, MPI_DOUBLE, buffer, buffer_size, &buffer_position, MPI_COMM_WORLD);
//     MPI_Pack(&W[w].E, 1, MPI_DOUBLE, buffer, buffer_size, &buffer_position, MPI_COMM_WORLD);
//     if(measure_SD && MC == DIFFUSION) {
//       MPI_Pack(W[w].CM[0], 3*SD.size, MPI_DOUBLE, buffer, buffer_size, &buffer_position, MPI_COMM_WORLD);
//       MPI_Pack(W[w].rreal[0], N*3, MPI_DOUBLE, buffer, buffer_size, &buffer_position, MPI_COMM_WORLD);
//     }
//   }

  /*if(measure_energy) {
    MPI_Pack(storeE,   storeEN, MPI_DOUBLE, buffer, buffer_size, &buffer_position, MPI_COMM_WORLD);
    MPI_Pack(storeEFF, storeEN, MPI_DOUBLE, buffer, buffer_size, &buffer_position, MPI_COMM_WORLD);
  }
  // empty energy arrays
  storeEN = 0;*/

  if(measure_OBDM) { // SaveOBDM()
    MPI_Pack(OBDM.f, OBDM.size, MPI_DOUBLE, buffer, buffer_size, &buffer_position, MPI_COMM_WORLD);
    MPI_Pack(OBDM.N, OBDM.size, MPI_DOUBLE, buffer, buffer_size, &buffer_position, MPI_COMM_WORLD);
  }
  if(measure_PairDistr) { //SavePairDistribution()
    MPI_Pack(&PD.times_measured, 1, MPI_DOUBLE, buffer, buffer_size, &buffer_position, MPI_COMM_WORLD);
    MPI_Pack(PD.N, PD.size, MPI_INT, buffer, buffer_size, &buffer_position, MPI_COMM_WORLD);
    if(MC == DIFFUSION) {
      MPI_Pack(&PD_pure.times_measured, 1, MPI_DOUBLE, buffer, buffer_size, &buffer_position, MPI_COMM_WORLD);
      MPI_Pack(PD_pure.N, PD.size, MPI_INT, buffer, buffer_size, &buffer_position, MPI_COMM_WORLD);
    }
  }
  if(measure_SD && MC == DIFFUSION) {
    MPI_Pack(SD.CM2, SD.size, MPI_DOUBLE, buffer, buffer_size, &buffer_position, MPI_COMM_WORLD);
    MPI_Pack(SD.N, SD.size, MPI_INT, buffer, buffer_size, &buffer_position, MPI_COMM_WORLD);
    MPI_Pack(rreal[0], N*3, MPI_DOUBLE, buffer, buffer_size, &buffer_position, MPI_COMM_WORLD);
  }

  /*if(measure_TBDM) SaveTBDM();
  if(measure_TBDM_MATRIX) SaveTBDMMatrix();
  if(measure_OBDM_MATRIX) SaveOBDMMatrix();
  if(measure_RadDistr && boundary != ONE_BOUNDARY_CONDITION) SaveRadialDistribution();
  if(measure_Sk) {
    SaveStaticStructureFactor();
    SaveStaticStructureFactorAngularPart();
  }
  if(measure_OP) SaveOrderParameter();
  if(measure_Lind) SaveLindemannRatio();
  if(var_par_array) SaveVarParWeights();*/
  time_parallel_transfer_pack += MPI_Wtime();

  time_parallel_transfer_send -= MPI_Wtime();
  MPI_Send(buffer, buffer_position, MPI_PACKED, dest, tag, MPI_COMM_WORLD);
  time_parallel_transfer_send += MPI_Wtime();
  time_parallel_transfer += MPI_Wtime();
  if(verbosity) Message(".done\n");

  return buffer_position;
}

/**************************** Receive Walkers Data ****************************/
void ReceiveWalkersData(int src, int wmin, int tag) {
  //int w;
  int buffer_position = 0;
  int i,j;

  if(verbosity) Message("receiving walkers data, wmin= %i ...", wmin);
  time_parallel_transfer -= MPI_Wtime(); // store moment when the function was called
  time_parallel_transfer_send -= MPI_Wtime();
  MPI_Recv(buffer, buffer_size, MPI_PACKED, src, tag, MPI_COMM_WORLD, &mpi_status);
  time_parallel_transfer_send += MPI_Wtime();

  time_parallel_transfer_pack -= MPI_Wtime();
//   MPI_Unpack(buffer, buffer_size, &buffer_position, &Nwalkers, 1, MPI_INT, MPI_COMM_WORLD);
//   for(w=wmin; w<wmin+Nwalkers; w++) {
//     MPI_Unpack(buffer, buffer_size, &buffer_position, W[w].x, N, MPI_DOUBLE, MPI_COMM_WORLD);
//     MPI_Unpack(buffer, buffer_size, &buffer_position, W[w].y, N, MPI_DOUBLE, MPI_COMM_WORLD);
//     MPI_Unpack(buffer, buffer_size, &buffer_position, W[w].z, N, MPI_DOUBLE, MPI_COMM_WORLD);
//     MPI_Unpack(buffer, buffer_size, &buffer_position, &W[w].U, 1, MPI_DOUBLE, MPI_COMM_WORLD);
//     MPI_Unpack(buffer, buffer_size, &buffer_position, &W[w].E, 1, MPI_DOUBLE, MPI_COMM_WORLD);
//     W[w].w = w; // update local index
//     if(measure_SD && MC == DIFFUSION) {
//       MPI_Unpack(buffer, buffer_size, &buffer_position, W[w].CM[0], 3*SD.size, MPI_DOUBLE, MPI_COMM_WORLD);
//       MPI_Unpack(buffer, buffer_size, &buffer_position, W[w].rreal[0], N*3, MPI_DOUBLE, MPI_COMM_WORLD);
//     }
//   }

  /*if(measure_energy) {
  MPI_Unpack(buffer, buffer_size, &buffer_position, buffer_double, storeEN, MPI_DOUBLE, MPI_COMM_WORLD);
  for(i=0; i<storeEN; i++) storeE[i] += buffer_double[i];
    MPI_Unpack(buffer, buffer_size, &buffer_position, buffer_double, storeEN, MPI_DOUBLE, MPI_COMM_WORLD);
    for(i=0; i<storeEN; i++) storeEFF[i] += buffer_double[i];
#ifdef LATTICE
    MPI_Unpack(buffer, buffer_size, &buffer_position, buffer_double, storeEN, MPI_DOUBLE, MPI_COMM_WORLD);
    for(i=0; i<storeEN; i++) bufferErelease[i] += buffer_double[i];
    MPI_Unpack(buffer, buffer_size, &buffer_position, buffer_double, storeEN, MPI_DOUBLE, MPI_COMM_WORLD);
    for(i=0; i<storeEN; i++) storeElat[i] += buffer_double[i];
#endif
  }*/

  if(measure_OBDM) { // SaveOBDM()
    MPI_Unpack(buffer, buffer_size, &buffer_position, buffer_double, OBDM.size, MPI_DOUBLE, MPI_COMM_WORLD);
    for(i=0; i<OBDM.size; i++) OBDM.f[i] += buffer_double[i];
    MPI_Unpack(buffer, buffer_size, &buffer_position, buffer_double, OBDM.size, MPI_DOUBLE, MPI_COMM_WORLD);
    for(i=0; i<OBDM.size; i++) OBDM.N[i] += buffer_double[i];
  }

  if(measure_PairDistr) { //SavePairDistribution();
    MPI_Unpack(buffer, buffer_size, &buffer_position, buffer_double, 1, MPI_DOUBLE, MPI_COMM_WORLD);
    PD.times_measured += buffer_double[0];
    MPI_Unpack(buffer, buffer_size, &buffer_position, buffer_int, PD.size, MPI_INT, MPI_COMM_WORLD);
    for(i=0; i<PD.size; i++) PD.N[i] += buffer_int[i];

    if(MC == DIFFUSION) {
      MPI_Unpack(buffer, buffer_size, &buffer_position, buffer_double, 1, MPI_DOUBLE, MPI_COMM_WORLD);
      PD_pure.times_measured += buffer_double[0];
      MPI_Unpack(buffer, buffer_size, &buffer_position, buffer_int, PD.size, MPI_INT, MPI_COMM_WORLD);
      for(i=0; i<PD.size; i++) PD_pure.N[i] += buffer_int[i];
    }
  }

  if(measure_SD && MC == DIFFUSION) {
    MPI_Unpack(buffer, buffer_size, &buffer_position, buffer_double, SD.size, MPI_DOUBLE, MPI_COMM_WORLD);
    for(i=0; i<SD.size; i++) SD.CM2[i] += buffer_double[i];
    MPI_Unpack(buffer, buffer_size, &buffer_position, buffer_double, SD.size, MPI_INT, MPI_COMM_WORLD);
    for(i=0; i<SD.size; i++) SD.N[i] += buffer_int[i];
    MPI_Unpack(buffer, buffer_size, &buffer_position, buffer_int, N*3, MPI_INT, MPI_COMM_WORLD);
    for(i=0; i<3; i++) {
      for(j=0; j<N; i++) {
        rreal[i][j] = buffer[i*N+j];
      }
    }
  }
  time_parallel_transfer_pack += MPI_Wtime();
  time_parallel_transfer += MPI_Wtime();
  if(verbosity) Message("done\n");
}

/******************************* Branching Walker ****************************/
void BranchingMPI(void) {
  int Nwalkers_pure_check;
  int i, j, w, wp, swap_index;
  DOUBLE dE, weight, log_weight;
  int w_branching, Nbranching; // not all walkers have to be branched
  int Npop_new;
  DOUBLE wmax = 0.;
  int src;
  int buffer_position = 0;
  static int tag = 123456;
  MPI_Request request;
  int protocol_length_sl[1000]; //!!! Nslaves
  int protocol_sl[100][201];// Nsl Nsl*2+1
  int wmin;
  double dummy;
  int sort_Nwalkers_max;
  int slaves_ordered_by_Nwalkers[1000]; //!!! ordered from large Nwalkers to smaller
  int Nwalkers_diff[1000];//!!!
  int NwalkersMean_branching;
  int go;
  int dst, tag_send;
  int flag_adjust_multiplicity; // 0 - no change, -1 reduce, +1 increase population

  //if(CheckWalkersAlive()) Message("before branching)");//???

  // Measure pure estimators
  if((iteration_global-1) % Nmeasure == 0) {
    //if(measure_R2) SaveMeanR2DMC();
    if(measure_OBDM) MeasureOBDM();
    if(measure_TBDM) MeasureTBDM();
  }

  if(!branchng_present) return;

  // calculate weights
  w_branching = 0;
  for(w=0; w<Nwalkers; w++) {
   if(W[w].status == ALIVE) {
      if(SmartMC == DMC_MOVE_ONE || SmartMC == DMC_QUARTIC) {
        weight = W[w].weight;
      }
      else {
        dE = Eo - 0.5*(W[w].E + W[w].Eold);
#ifdef UNITS_SET_2M_TO_1
        log_weight = dt*dE;
#else
        log_weight = dt_2x*dE;
#endif
        weight = Exp(log_weight);
      }

      branching_weight[w_branching] = weight; // assign weights array
      W[w].weight = weight;
      if(weight>wmax) wmax = weight; // find maximal weight
      w_branching++;
    }
  }
  Nbranching = w_branching;

  Npop_new = 0;
  for(w=0; w<Nbranching; w++) {
    branching_xi[w] = Random(); // throw random numbers
    branching_multiplicity[w] = (int)(branching_weight[w]+branching_xi[w]); // define number of sons
    if(branching_multiplicity[w]<0) {
      Warning("Branching: huge multiplicity (larger than INT_MAX)\n");
      branching_multiplicity[w] = NwalkersMax;
    }
    Npop_new += branching_multiplicity[w]; // check if the new population size is within the proper limits
  }

  // do the branching
  w_branching = 0;
  for(w=0; w<NwalkersMax; w++) {
    if(W[w].status == VIRTUAL) { // do not branch virtual (OBDM pure) walker
      W[w].weight += log(branching_weight[w_branching]); // virtual walker stores logarithmic weight
      if((iteration_global-1) % Nmeasure == 0) { // in case of OBDM pure measurement
        W[w].OBDMweight[W[w].OBDM_position] = W[w].weight; // store total weight
        //if(W[w].OBDM_position) W[w].OBDMweight[W[w].OBDM_position] /= W[w].OBDMweight[W[w].OBDM_position-1]; // convert to relative weight
        W[w].OBDM_position++;
        if(W[w].OBDM_position>grid_pure_block) Warning("OBDM_position (%i) of Walker %i is out of range (%i), iter %i!\n", W[w].OBDM_position, w, grid_pure_block, iteration_global);
      }
    }
    else if(W[w].status == ALIVE){ // not virtual walker
      if(branching_multiplicity[w_branching] == 0) {
        W[w].status = KILLED;
        W[w].weight = branching_weight[w_branching];
      }
      else {
        W[w].weight = branching_weight[w_branching];// / (DOUBLE) multiplicity;
        for(i=1; i<branching_multiplicity[w_branching]; i++) { // make replicae
          wp = dead_walkers_storage[Ndead_walkers-1];
          CopyWalker(&W[wp], &W[w]);
          W[wp].status = REINCARNATION;
          W[wp].weight = 0.;
          Ndead_walkers--;
        }
      }
    }
    w_branching++;
  }

  //Eo = EnergyO();

  Nwalkers = 0;
  Nwalkersw = 0;
  Nwalkers_weight_killed = 0.;
  for(w=0; w<NwalkersMax; w++) {
    if(W[w].status == REINCARNATION)
      W[w].status = ALIVE;
    else if(W[w].status == KILLED) {
      W[w].status = DEAD;
      Nwalkers_weight_killed += W[w].weight;
    }
  }

  for(w=0; w<NwalkersMax; w++) { // order all alive walkers
    if(W[w].status == ALIVE) {
      if(w != Nwalkers) { // walker should be moved to Nwalkers position
        if(W[Nwalkers].status == DEAD) { // overwrite
          CopyWalker(&W[Nwalkers], &W[w]);
          W[w].status = DEAD;
        }
        else { // exchange
          CopyWalker(&W[dead_walkers_storage[Ndead_walkers-1]], &W[Nwalkers]);
          W[dead_walkers_storage[Ndead_walkers-1]].status = W[Nwalkers].status;
          CopyWalker(&W[Nwalkers], &W[w]);
          W[Nwalkers].status = ALIVE;
          CopyWalker(&W[w], &W[dead_walkers_storage[Ndead_walkers-1]]);
          W[w].status = W[dead_walkers_storage[Ndead_walkers-1]].status;
          W[dead_walkers_storage[Ndead_walkers-1]].status = DEAD;
        }
      }
      Nwalkers++;
      Nwalkersw += W[w].weight;
    }
  }

  Nwalkers_pure_check = 0;
  for(w=0; w<NwalkersMax; w++) { // order all virtual walkers
    if(W[w].status == VIRTUAL) {
      if(w != Nwalkers+Nwalkers_pure_check) {
        CopyWalker(&W[Nwalkers+Nwalkers_pure_check], &W[w]);
        W[Nwalkers+Nwalkers_pure_check].status = VIRTUAL;
        W[w].status = DEAD;
      }
      Nwalkers_pure_check++;
    }
  }
  if(Nwalkers_pure_check != Nwalkers_pure) {
    Warning("  Branching: Nwalkers_pure is not self consistent (real %i, expected %i)\n", Nwalkers_pure_check, Nwalkers_pure);
    Nwalkers_pure = Nwalkers_pure_check;
  }

  Ndead_walkers = 0;
  for(w=NwalkersMax-1; w>=Nwalkers+Nwalkers_pure; w--) {
    dead_walkers_storage[Ndead_walkers++] = w;
  }
  // end of serial part

  tag++;
  if(myid == MASTER) { // the master code
    time_parallel_transfer -= MPI_Wtime();
    time_parallel_transfer_send -= MPI_Wtime();
    //for(src=1; src<Nslaves; src++) MPI_Irecv(bufferSl[src], buffer_size, MPI_PACKED, src, tag, MPI_COMM_WORLD, &request);
    for(src=1; src<Nslaves; src++) MPI_Recv(bufferSl[src], buffer_size, MPI_PACKED, src, tag, MPI_COMM_WORLD, &mpi_status);
    time_parallel_transfer += MPI_Wtime();
    time_parallel_transfer_send += MPI_Wtime();
    //if(verbosity) Message("wait2 E\n");
    //MPI_Waitall(Nslaves-1, &requestSl[1], mpi_status_Sl); // wait until non-blocking receive is completed
    //if(verbosity) Message("wait2 %i done\n");
    tag++;

    // initialize with zeros
    Eo = 0.; // new energy shift in the branching
    NwalkersAll = 0;

    // fill data from master into arrays
    src = 0;
    wmin = 0;
    NwalkersSl[src] = Nwalkers;

    for(w=wmin; w<wmin+NwalkersSl[src]; w++) {
      Eo += W[w].E;
    }
    wmin += NwalkersSl[src];
    NwalkersAll += NwalkersSl[src];

    // fill data from slaves into arrays
    time_parallel_transfer_pack -= MPI_Wtime();
    for(src=1; src<Nslaves; src++) {
      MPI_Unpack(bufferSl[src], buffer_size, &buffer_position, &NwalkersSl[src], 1, MPI_INT, MPI_COMM_WORLD);
      for(w=wmin; w<wmin+NwalkersSl[src]; w++) {
        MPI_Unpack(bufferSl[src], buffer_size, &buffer_position, &dummy, 1, MPI_DOUBLE, MPI_COMM_WORLD);
        Eo += W[w].E;
      }
      wmin += NwalkersSl[src];
      NwalkersAll += NwalkersSl[src];
    }
    time_parallel_transfer_pack += MPI_Wtime();

    // define new energy shift
    Eo /= (double) NwalkersAll;

    // check limits of the population and adjust if necessary
    flag_adjust_multiplicity = 0; // so far no change
    if(dt>0) { // amplify > 1
        //Message("%lf + (%lf)-> ", Eo, -500.*log((double)NwalkersAll_branching/(double)Npop));
        //Eo -= 500.*log((double)NwalkersAll_branching/(double)Npop);

      if(NwalkersAll > Npop_max) {
        flag_adjust_multiplicity = -1;
        //multiplicity = (int) ((DOUBLE) multiplicity*reduce + Random());
        //Message("  Nw > Npop_max (%i > %i) iter = %i, Eo = %lf corrected to ", NwalkersAll_branching, Npop_max, iteration_global, Eo);
        //Eo -= log(amplify)/dt;
        //Eo -= log(NwalkersAll_branching/(double)Npop)/dt;
        //Eo *= 1.05;
        //Message("%lf\n", Eo);
        //if(NwalkersAll_branching > NwalkersMax) Error("branching: population size is too large\n");
      }
      else if(NwalkersAll < Npop_min) {
        flag_adjust_multiplicity = +1;
        //multiplicity = (int) ((DOUBLE) multiplicity*amplify + Random());
        //Message("  Nw < Npop_min (%i < %i) iter = %i, Eo = %lf corrected to \n", NwalkersAll_branching, Npop_max, iteration_global, Eo);
        //Eo += log(amplify)/dt;
        //Eo -= log(NwalkersAll_branching/(double)Npop)/dt;
        //Eo /= 1.05;
        //Message("%lf\n", Eo);
      }
    }

    // sort and distribute walkers
    for(i=0; i<Nslaves; i++) slaves_ordered_by_Nwalkers[i] = i; // initialize for sorting
    NwalkersMean_branching = (int)ceil((double)NwalkersAll /(double) Nslaves);
    if(verbosity) Message("Nwalkers all %i, mean %i\n", NwalkersAll, NwalkersMean_branching);

    if(verbosity) Message("\nsorting\nNw sl exc\n");
    for(i=0; i<Nslaves; i++) { // bubble sort
      sort_Nwalkers_max = 0;
      for(j=i; j<Nslaves; j++) {
        if(NwalkersSl[slaves_ordered_by_Nwalkers[j]]>=sort_Nwalkers_max) {
          sort_Nwalkers_max = NwalkersSl[slaves_ordered_by_Nwalkers[j]];
          swap_index = slaves_ordered_by_Nwalkers[j];
          slaves_ordered_by_Nwalkers[j] = slaves_ordered_by_Nwalkers[i];
          slaves_ordered_by_Nwalkers[i] = swap_index;
        }
      }
      Nwalkers_diff[slaves_ordered_by_Nwalkers[i]] = NwalkersSl[slaves_ordered_by_Nwalkers[i]] - NwalkersMean_branching;
      if(verbosity) Message("%i %i %i\n",NwalkersSl[slaves_ordered_by_Nwalkers[i]],slaves_ordered_by_Nwalkers[i],Nwalkers_diff[slaves_ordered_by_Nwalkers[i]]);
    }

    j=0;
    for(i=0; i<Nslaves; i++) j += Nwalkers_diff[slaves_ordered_by_Nwalkers[i]];
    if(verbosity) Message("Mean disbalance %i\n", j);

    // decide how to distribute walkers
    for(i=0; i<Nslaves; i++) {
      protocol_length_sl[i] = 1; // for a moment zero message length
      protocol_sl[i][0] = 0; // [total number of walkers to pass / message length] [Nw1] [sl1] [Nw2] [sl2] ...
    }
    j = Nslaves - 1; // take walkers from [i] and send them to [j]
    for(i=0; i<Nslaves; i++) {
      if(Nwalkers_diff[slaves_ordered_by_Nwalkers[i]]>0) {
        if(verbosity) Message("sl %i has to distribute %i walkers: ", slaves_ordered_by_Nwalkers[i], Nwalkers_diff[slaves_ordered_by_Nwalkers[i]]);
        go = ON;
        while(go) {
          if(Nwalkers_diff[slaves_ordered_by_Nwalkers[i]]+Nwalkers_diff[slaves_ordered_by_Nwalkers[j]]<=0) { // i.e. all walkers from [i] can fit to [j]
            go = OFF;
            if(verbosity) Message("%i to %i\n", Nwalkers_diff[slaves_ordered_by_Nwalkers[i]], slaves_ordered_by_Nwalkers[j]);

            protocol_sl[slaves_ordered_by_Nwalkers[i]][0] += Nwalkers_diff[slaves_ordered_by_Nwalkers[i]];
            protocol_sl[slaves_ordered_by_Nwalkers[j]][0] -= Nwalkers_diff[slaves_ordered_by_Nwalkers[i]];
            protocol_sl[slaves_ordered_by_Nwalkers[i]][protocol_length_sl[slaves_ordered_by_Nwalkers[i]]++] = Nwalkers_diff[slaves_ordered_by_Nwalkers[i]];
            protocol_sl[slaves_ordered_by_Nwalkers[i]][protocol_length_sl[slaves_ordered_by_Nwalkers[i]]++] = slaves_ordered_by_Nwalkers[j];
            protocol_sl[slaves_ordered_by_Nwalkers[j]][protocol_length_sl[slaves_ordered_by_Nwalkers[j]]++] = -Nwalkers_diff[slaves_ordered_by_Nwalkers[i]];
            protocol_sl[slaves_ordered_by_Nwalkers[j]][protocol_length_sl[slaves_ordered_by_Nwalkers[j]]++] = slaves_ordered_by_Nwalkers[i];

            Nwalkers_diff[slaves_ordered_by_Nwalkers[j]] += Nwalkers_diff[slaves_ordered_by_Nwalkers[i]];
            Nwalkers_diff[slaves_ordered_by_Nwalkers[i]] = 0;
          }
          else { // i.e. walkers fit only partially
            if(verbosity) Message("%i to %i, ", -Nwalkers_diff[slaves_ordered_by_Nwalkers[j]], slaves_ordered_by_Nwalkers[j]);

            protocol_sl[slaves_ordered_by_Nwalkers[i]][0] -= Nwalkers_diff[slaves_ordered_by_Nwalkers[j]];
            protocol_sl[slaves_ordered_by_Nwalkers[j]][0] += Nwalkers_diff[slaves_ordered_by_Nwalkers[j]];
            protocol_sl[slaves_ordered_by_Nwalkers[i]][protocol_length_sl[slaves_ordered_by_Nwalkers[i]]++] = -Nwalkers_diff[slaves_ordered_by_Nwalkers[j]];
            protocol_sl[slaves_ordered_by_Nwalkers[i]][protocol_length_sl[slaves_ordered_by_Nwalkers[i]]++] = slaves_ordered_by_Nwalkers[j];
            protocol_sl[slaves_ordered_by_Nwalkers[j]][protocol_length_sl[slaves_ordered_by_Nwalkers[j]]++] = Nwalkers_diff[slaves_ordered_by_Nwalkers[j]];
            protocol_sl[slaves_ordered_by_Nwalkers[j]][protocol_length_sl[slaves_ordered_by_Nwalkers[j]]++] = slaves_ordered_by_Nwalkers[i];
            Nwalkers_diff[slaves_ordered_by_Nwalkers[i]] += Nwalkers_diff[slaves_ordered_by_Nwalkers[j]];
            Nwalkers_diff[slaves_ordered_by_Nwalkers[j]] = 0;
            j--;
          }
        }
      }
    }
    if(verbosity) Message("after distribution\n");
    for(i=0; i<Nslaves; i++) 
      if(verbosity) Message("%i %i\n", slaves_ordered_by_Nwalkers[i], Nwalkers_diff[slaves_ordered_by_Nwalkers[i]]);

    j=0;
    for(i=0; i<Nslaves; i++) j += Nwalkers_diff[slaves_ordered_by_Nwalkers[i]];
    if(verbosity) Message("Mean disbalance %i\n", j);

    if(verbosity) Message("protocol\n");
    for(i=0; i<Nslaves; i++) {
      if(verbosity) Message("\n Sl[%i]", i);
      for(j=0; j<protocol_length_sl[i]; j++) {
        if(verbosity) Message(" %i", protocol_sl[i][j]);
      }
    }

    // send new energy and protocol to all slaves
    if(verbosity) Message("Eo = %lf\n", Eo);
    time_parallel_transfer -= MPI_Wtime();
    for(i=1; i<Nslaves; i++) {
      time_parallel_transfer_pack -= MPI_Wtime();
      buffer_position = 0;
      MPI_Pack(&Eo, 1, MPI_DOUBLE, bufferSl[i], buffer_size, &buffer_position, MPI_COMM_WORLD); // pack Eo
      MPI_Pack(&flag_adjust_multiplicity, 1, MPI_INT, bufferSl[i], buffer_size, &buffer_position, MPI_COMM_WORLD); // pack flag
      protocol_sl[i][0] = protocol_length_sl[i]; // change Nwalkers to length of the protocol message
      MPI_Pack(&protocol_sl[i], protocol_length_sl[i], MPI_INT, bufferSl[i], buffer_size, &buffer_position, MPI_COMM_WORLD); // pack protocol
      time_parallel_transfer_pack += MPI_Wtime();
      time_parallel_transfer_send -= MPI_Wtime();
      //MPI_Isend(bufferSl[i], buffer_position, MPI_PACKED, i, tag, MPI_COMM_WORLD, &requestSl[i]);
      MPI_Send(bufferSl[i], buffer_position, MPI_PACKED, i, tag, MPI_COMM_WORLD);
      time_parallel_transfer_send += MPI_Wtime();
    }
    if(verbosity) Message("wait3 protocol\n");
    //MPI_Waitall(Nslaves-1, &requestSl[1], mpi_status_Sl); // wait until non-blocking receive is completed
    if(verbosity) Message("wait3 %i done\n");
    time_parallel_transfer += MPI_Wtime();
  }
  else { //slave code
    // pass walker energy to master
    time_parallel_transfer -= MPI_Wtime();
    time_parallel_transfer_pack -= MPI_Wtime();
    buffer_position = 0;

    MPI_Pack(&Nwalkers, 1, MPI_INT, bufferSl[myid], buffer_size, &buffer_position, MPI_COMM_WORLD);
    for(w=0; w<Nwalkers; w++) MPI_Pack(&W[w].E, 1, MPI_DOUBLE, bufferSl[myid], buffer_size, &buffer_position, MPI_COMM_WORLD);
    time_parallel_transfer_pack += MPI_Wtime();
    time_parallel_transfer_send -= MPI_Wtime();
    MPI_Send(bufferSl[myid], buffer_position, MPI_PACKED, MASTER, tag, MPI_COMM_WORLD);
    //MPI_Isend(bufferSl[myid], buffer_position, MPI_PACKED, MASTER, tag, MPI_COMM_WORLD, &request);
    //MPI_Wait(&request,&mpi_status);
    time_parallel_transfer_send += MPI_Wtime();
    time_parallel_transfer += MPI_Wtime();
    tag++;

    // receive protocol
    time_parallel_transfer -= MPI_Wtime();
    time_parallel_transfer_send -= MPI_Wtime();
    MPI_Recv(bufferSl[myid], buffer_size, MPI_PACKED, 0, tag, MPI_COMM_WORLD, &mpi_status);
    //MPI_Irecv(bufferSl[myid], buffer_size, MPI_PACKED, 0, tag, MPI_COMM_WORLD, &request);
    //MPI_Wait(&request, &mpi_status);
    time_parallel_transfer_send += MPI_Wtime();
    buffer_position = 0;
    MPI_Unpack(bufferSl[myid], buffer_size, &buffer_position, &Eo, 1, MPI_DOUBLE, MPI_COMM_WORLD); // get Eo
    if(verbosity) Message("Eo = %lf\n", Eo);
    MPI_Unpack(bufferSl[myid], buffer_size, &buffer_position, &flag_adjust_multiplicity, 1, MPI_INT, MPI_COMM_WORLD); // get flag
    MPI_Unpack(bufferSl[myid], buffer_size, &buffer_position, &protocol_sl[myid][0], 1, MPI_INT, MPI_COMM_WORLD); // message length
    protocol_length_sl[myid] = protocol_sl[myid][0];
    if(verbosity) Message("protocol: length %i", protocol_length_sl[myid]);
    MPI_Unpack(bufferSl[myid], buffer_size, &buffer_position, &protocol_sl[myid][1], protocol_length_sl[myid]-1, MPI_INT, MPI_COMM_WORLD);
    time_parallel_transfer += MPI_Wtime();
  }

  for(i=1; i<protocol_length_sl[myid]; i++) if(verbosity) Message(" %i", protocol_sl[myid][i]);

  // send or receive data [M/S]
  j=0; // counts number of messages
  for(i=1; i<protocol_length_sl[myid]; i+=2) {
    tag_send = iteration_global*abs(protocol_sl[myid][i]) % 16000;
    if(protocol_sl[myid][i]>0) { // i.e. send walkers
      dst = protocol_sl[myid][i+1];
      if(verbosity) Message("  => SL %i sending %i walkers to %i\n", myid, protocol_sl[myid][i], dst);
      Nwalkers -= protocol_sl[myid][i]; // delete walkers which were sent
      SendWalkers(dst, Nwalkers-1, protocol_sl[myid][i], tag_send, &requestSl[j]);
    }
    else { // i.e. recieve walkers
      src = protocol_sl[myid][i+1];
      if(verbosity) Message("  => SL %i receiving %i walkers from %i\n", myid, -protocol_sl[myid][i], src);
      ReceiveWalkers(src, Nwalkers-1, tag_send, &requestSl[j]);
      MPI_Wait(&requestSl[j], &mpi_status);
      ReceiveWalkersUnpack(src, Nwalkers-1);
      Nwalkers -= protocol_sl[myid][i]; // delete walkers which were sent
    }
    j++;
    NwalkersSl[myid] -= protocol_sl[myid][i];
    //NwalkersSl[dst] += protocol_sl[myid][i];???
  }
  // wait only if walkers were transferred
  //f(j>0) MPI_Waitall(j, &requestSl[0], mpi_status_Sl); // wait until non-blocking send is completed
  // walkers are not still unpacked!

  // adjust population limits if requested
  if(flag_adjust_multiplicity != 0) {
    Message(flag_adjust_multiplicity>0?"(+)":"(-)");
    Message("%i (%i)-> ", Nwalkers, Nbranching);//???
    if(Nwalkers == 0) Error("done");

    Npop_new = 0;
    for(w=0; w<Nbranching; w++) {
      //branching_xi[w] = Random(); // throw random numbers
      //branching_multiplicity[w] = (int)(branching_weight[w]+branching_xi[w]); // define number of sons
      if(flag_adjust_multiplicity =- -1) {
        //multiplicity = (int) ((DOUBLE) multiplicity*reduce + Random());
        branching_multiplicity[w] = (int) ((DOUBLE) branching_multiplicity[w]*reduce + Random());
      }
      else {
        branching_multiplicity[w] = (int) ((DOUBLE) branching_multiplicity[w]*amplify + Random());
      }
      Npop_new += branching_multiplicity[w]; // check if the new population size is within the proper limits
    }

    // do the branching
    w_branching = 0;
    for(w=0; w<NwalkersMax; w++) {
      if(W[w].status == VIRTUAL) { // do not branch virtual (OBDM pure) walker
        W[w].weight += log(branching_weight[w_branching]); // virtual walker stores logarithmic weight
        if((iteration_global-1) % Nmeasure == 0) { // in case of OBDM pure measurement
          W[w].OBDMweight[W[w].OBDM_position] = W[w].weight; // store total weight
          W[w].OBDM_position++;
          if(W[w].OBDM_position>grid_pure_block) Warning("OBDM_position (%i) of Walker %i is out of range (%i), iter %i!\n", W[w].OBDM_position, w, grid_pure_block, iteration_global);
        }
      }
      else if(W[w].status == ALIVE){ // not virtual walker
        if(branching_multiplicity[w_branching] == 0) {
          W[w].status = DEAD;
          W[w].weight = branching_weight[w_branching];
        }
        else {
          W[w].weight = branching_weight[w_branching];// / (DOUBLE) multiplicity;
          for(i=1; i<branching_multiplicity[w_branching]; i++) { // make replicae
            wp = dead_walkers_storage[Ndead_walkers-1];
            CopyWalker(&W[wp], &W[w]);
            W[wp].status = ALIVE;
            W[wp].weight = 0.;
            Ndead_walkers--;
          }
        }
      }
      w_branching++;
    }

    Nwalkers = 0;
    Nwalkersw = 0;
    Nwalkers_weight_killed = 0.;
    for(w=0; w<NwalkersMax; w++) { // order all alive walkers
      if(W[w].status == ALIVE) {
        if(w != Nwalkers) { // walker should be moved to Nwalkers position
          if(W[Nwalkers].status == DEAD) { // overwrite
            CopyWalker(&W[Nwalkers], &W[w]);
            W[w].status = DEAD;
          }
          else { // exchange
            CopyWalker(&W[dead_walkers_storage[Ndead_walkers-1]], &W[Nwalkers]);
            W[dead_walkers_storage[Ndead_walkers-1]].status = W[Nwalkers].status;
            CopyWalker(&W[Nwalkers], &W[w]);
            W[Nwalkers].status = ALIVE;
            CopyWalker(&W[w], &W[dead_walkers_storage[Ndead_walkers-1]]);
            W[w].status = W[dead_walkers_storage[Ndead_walkers-1]].status;
            W[dead_walkers_storage[Ndead_walkers-1]].status = DEAD;
          }
        }
        Nwalkers++;
        Nwalkersw += W[w].weight;
      }
    }

    Nwalkers_pure_check = 0;
    for(w=0; w<NwalkersMax; w++) { // order all virtual walkers
      if(W[w].status == VIRTUAL) {
        if(w != Nwalkers+Nwalkers_pure_check) {
          CopyWalker(&W[Nwalkers+Nwalkers_pure_check], &W[w]);
          W[Nwalkers+Nwalkers_pure_check].status = VIRTUAL;
          W[w].status = DEAD;
        }
        Nwalkers_pure_check++;
      }
    }
    if(Nwalkers_pure_check != Nwalkers_pure) {
      Warning("  Branching: Nwalkers_pure is not self consistent (real %i, expected %i)\n", Nwalkers_pure_check, Nwalkers_pure);
      Nwalkers_pure = Nwalkers_pure_check;
    }

    Ndead_walkers = 0;
    for(w=NwalkersMax-1; w>=Nwalkers+Nwalkers_pure; w--) {
      dead_walkers_storage[Ndead_walkers++] = w;
    }

    // update number of walkers on master
    tag++;
    if(myid == MASTER) { // the master code
      time_parallel_transfer -= MPI_Wtime();
      time_parallel_transfer_send -= MPI_Wtime();

      NwalkersSl[0] = Nwalkers; // number of walkers on master
      // receive number of walkers on slaves
      for(src=1; src<Nslaves; src++) MPI_Recv(&NwalkersSl[src], buffer_size, MPI_PACKED, src, tag, MPI_COMM_WORLD, &mpi_status);
      time_parallel_transfer += MPI_Wtime();
      time_parallel_transfer_send += MPI_Wtime();

      NwalkersAll = 0;
      for(i=0; i<Nslaves; i++) NwalkersAll += NwalkersSl[i]; // update total number of walkers
      walker_indices[0] = 0; // update walker index
      for(i=1; i<Nslaves; i++) walker_indices[i] = walker_indices[i-1] + NwalkersSl[i-1]; // first index of node i
    }
    else { //slave code
      // pass walker energy to master
      time_parallel_transfer -= MPI_Wtime();
      time_parallel_transfer_send -= MPI_Wtime();
      buffer_position = 0;
      MPI_Send(&Nwalkers, buffer_position, MPI_DOUBLE, MASTER, tag, MPI_COMM_WORLD);
      time_parallel_transfer_send += MPI_Wtime();
      time_parallel_transfer += MPI_Wtime();
    }
  }

  //Message("\nNwalkers M %i S %i (%i < %i <%i) ", NwalkersSl[0], NwalkersSl[1], Npop_min, NwalkersSl[0]+NwalkersSl[1], Npop_max);!!!

  if(verbosity) Message("MPI branching done\n");
}

#endif

/**************************** Open MP ****************************************/
#ifdef _OPENMP

#include <stdio.h>
#include <math.h>
#include "omp.h"
#include "parallel.h"
#include "utils.h"
#include "memory.h"
#include "vmc.h"
#include "rw.h"
#include "randnorm.h"

int myid;    // process identification
int Nslaves; // number of jobs

void ParallelInitialize(int argn, char **argv) {

  Message("\nStarting OpenMP parallelization ...\n");
#pragma omp parallel
  {
    Nslaves = omp_get_num_threads();
  }
  if(Nslaves>OPENMP_MAX_NPROC) {
    Error(" Cannot use all nodes: available %i, used %i\n  Change OPENMP_MAX_NPROC and recompile\n", Nslaves, OPENMP_MAX_NPROC);
    //omp_set_num_threads(OPENMP_MAX_NPROC);
  }

#pragma omp parallel
  {
    myid = omp_get_thread_num();  // Obtain thread number
    Message("  Hello World from thread = %d\n", myid);

    if(myid == 0) { // Only master thread does this
      Message("  Hello from Master\n");
      Nslaves = omp_get_num_threads();
    }
  } // All threads join master thread and disband

  Message("  total number of processors is %i\n", Nslaves);
  Message("done\n\n");
}

/**************************** Open MP: Random Critical ************************/
DOUBLE RandomCritical(void) {
  DOUBLE xi;
#pragma omp critical
  {
  xi = ran2(&rand_seed);
  }
  return xi;
}

/**************************** ParallelRandInit ********************************/
void ParallelRandInit(void) {
#ifndef USE_RAN2_GENERATOR
  Error("define USE_RAN2_GENERATOR");
#endif

#pragma omp parallel
  {
  myid = omp_get_thread_num();
  //rand_seed_parallel[myid] = (long int)(pow(2.,myid+3.)-1.);
  if(myid > rand_seed_parallel_length) Error("  Not enough seeds for parallel random number generator, %i > %i, recompile\n", myid, rand_seed_parallel_length);
  Message(" node %i, Ran2() random number generator started with seed %i\n", myid, rand_seed_parallel[myid]);

  if(RandomSys()<0) Error("bad implementation of RandomSys(). It should not produce negative numbers");
  if(RandomSys()>1) Error("bad implementation of RandomSys(). It should not produce numbers larger than 1");
}
}

#endif
