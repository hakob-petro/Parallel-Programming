#include <stdio.h>
#include <mpi.h>
#include <sys/time.h>
int main(int argc, char *argv[])
{
  int myrank, nprocs, n = 0;
  struct timeval t1, t2, dt;
  struct timezone tz;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  if (myrank == 0) {
    gettimeofday(&t1, &tz);
    MPI_Send(&t1, 1, MPI_INT, 1, 99, MPI_COMM_WORLD);
    printf("I'm the node %d, and I sent the message.\n", myrank);
  }
  else if (myrank == 1) {
    MPI_Recv(&t1, 1, MPI_INT, 0, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    gettimeofday(&t2, &tz);
    dt.tv_sec = t2.tv_sec - t1.tv_sec;
    dt.tv_usec = t2.tv_usec - t1.tv_usec;
    if(dtv.tv_usec<0) dtv.tv_sec--; dtv.tv_usec+=1000000; 
    printf("I'm the node %d, and I receive the message. Time passed: %f .\n", myrank, dt.tv_sec*1000+dt_tv_usec/1000);
  }

  MPI_Finalize();
  return 0;
}