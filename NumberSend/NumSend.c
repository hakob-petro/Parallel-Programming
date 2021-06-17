#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
  int myrank, nprocs, n = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  if (myrank == 0) {
    printf("I am node %d of %d, and the number = %d\n", myrank, nprocs, n);
    n += 1;
    MPI_Send(&n, 1, MPI_INT, myrank + 1, 99, MPI_COMM_WORLD);
  }
  else if (myrank != nprocs - 1) {
  	const double t_1 = MPI_Wtime();
    MPI_Recv(&n, 1, MPI_INT, myrank - 1, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printf("I am node %d of %d, and the number = %d\n", myrank, nprocs, n);
    n += 1;
    MPI_Send(&n, 1, MPI_INT, myrank + 1, 99, MPI_COMM_WORLD);
    const double t_2 = MPI_Wtime();
    printf("%f\n", &(t_2-t_1));
  }
  else {
    MPI_Recv(&n, 1, MPI_INT, myrank - 1, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printf("I am node %d of %d, and the number = %d\n", myrank, nprocs, n);
  }


  MPI_Finalize();
  return 0;
}