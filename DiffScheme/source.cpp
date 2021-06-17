/*
-------------MIPT-------------
--Petrosyan Hakob, group 814--
--------No copyright.---------
*/

#include <iostream>
#include <cmath>
#include <vector>
#include <mpi.h>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;

const static double pi = 3.14159265359;
const static double X = 1.0;
const static double T = 1.0;
const static int K = 20000;
const static int M = 10000;
const static double tau = (T / K);
const static double h = (X / M);

/*
====================================================================
= The equation: delta(u)/delta(t) + 2 * delta(u)/delta(x) = x + t  =
=                    0 < x <= 1, 0 < t <= 1                                =
= Initial conditions:  u(t)|x=0 == e^(-t), u(x)|t=0 == cos(pi*x)   =
=                                                                  =
= Grid:                                                            =
= u[K][0] ... u[K][m] ... u[K][M]                                  =
= .                                                                =
= .                                                                =
= u[k][0] ... u[k][m] ... u[k][M]                                  =
= .                                                                =
= .                                                                =
= u[0][0] ... u[0][m] ... u[0][M]                                  = 
=                                                                  =
= The scheme is stable when 2*tau > h                              =
====================================================================
*/

// Initial conditions
double u_t(const double &t) {
  return exp(-t);
}
double u_x(const double &x) {
  return cos(pi * x);
}

// Non-linear part
double f_x_t(const double &t, const double &x) {
  return t + x;
}

// This function pushs the result into files
void push_result(vector<vector<double> > &u, string &path,
                 int &myrank, double &parallel_time) {
  ofstream output(path.c_str(), ios::app);
  output << "I'm the node with rank: " << myrank << " and I completed my part in " << parallel_time << " seconds." << endl;
  if(myrank == 0){
    for (int m = 0; m <= M; m++) {
      output << u[K][m] << " ";
    }
    output << endl;
  }
}

int main(int argc, char *argv[]) {
  int myrank, nprocs;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  // Result 2D-vector
  vector<vector<double> > u(K+1, vector<double>(M+1));

  // Mesurement start time
  const double t_1 = MPI_Wtime();

  // Calculate initial conditions
  if(myrank == 0){
    for (int k = 0; k <= K; k++) {
      u[k][0] = u_t(k*tau);
    }
    for (int m = 0; m <= M/nprocs; m++) {
      u[0][m] = u_x(m*h);
    }
  } 
  else{
    for (int m = myrank*M/nprocs - 1; m <= (myrank+1)*M/nprocs; m++) {
      u[0][m] = u_x(m*h);
    }
  }

  // Calculate inner values
  for (int k = 0; k <= K-1; k++) {
    if (myrank == 0) {
      for (int m = 1; m <= M/nprocs - 1; m++) {
        u[k+1][m] = tau*f_x_t(k*tau, m*h) - (u[k][m+1] - u[k][m-1])*tau/h + (u[k][m+1] + u[k][m-1])/2;
      }
    }
    else {
      for (int m = myrank*M/nprocs; m <= (myrank+1)*M/nprocs - 1; m++) {
        u[k+1][m] = tau*f_x_t(k*tau, m*h) - (u[k][m+1] - u[k][m-1])*tau/h + (u[k][m+1] + u[k][m-1])/2;
      }
    }

    // Calculating the value for the rightmost column
    if (myrank == nprocs - 1) {
      u[k+1][M] = tau*f_x_t(k*tau, M*h) - 2*(u[k][M] - u[k][M-1])*tau/h + u[k][M];
    }
    
    // For odd numbers
    if (myrank % 2) {
      if (myrank > 0) MPI_Send(&u[k+1][myrank*M/nprocs], 1, MPI_DOUBLE, myrank - 1, 99, MPI_COMM_WORLD);
      if (myrank > 0) MPI_Recv(&u[k+1][myrank*M/nprocs - 1], 1, MPI_DOUBLE, myrank - 1, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      if (myrank < nprocs - 1) MPI_Recv(&u[k + 1][(myrank+1)*M/nprocs], 1, MPI_DOUBLE, myrank + 1, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      if (myrank < nprocs - 1) MPI_Send(&u[k + 1][(myrank+1)*M/nprocs - 1], 1, MPI_DOUBLE, myrank + 1, 99, MPI_COMM_WORLD);
    }
    // For even numbers
    else {
      if (myrank < nprocs - 1) MPI_Recv(&u[k+1][(myrank+1)*M/nprocs], 1, MPI_DOUBLE, myrank + 1, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      if (myrank < nprocs - 1) MPI_Send(&u[k+1][(myrank+1)*M/nprocs - 1], 1, MPI_DOUBLE, myrank + 1, 99, MPI_COMM_WORLD);
      if (myrank > 0) MPI_Send(&u[k+1][myrank*M/nprocs], 1, MPI_DOUBLE, myrank - 1, 99, MPI_COMM_WORLD);
      if (myrank > 0) MPI_Recv(&u[k+1][myrank*M/nprocs - 1], 1, MPI_DOUBLE, myrank - 1, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
  }

  // Collect all values with MPI_Gather
  MPI_Gather(&u[K][myrank*M/nprocs], M/nprocs, MPI_DOUBLE, &u[K][0], M/nprocs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  // The last element we can ignore if you don't need that;
  if (myrank == nprocs - 1) MPI_Send(&u[K][M], 1, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD);
  if (myrank == 0) MPI_Recv(&u[K][M], 1, MPI_DOUBLE, nprocs - 1, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


  const double t_2 = MPI_Wtime();
  double parallel_time = t_1 - t_2;

  ostringstream spath;
  spath << nprocs << "nodes/parallel.txt";
  string path = spath.str();
  push_result(u, path, nprocs, myrank, parallel_time);

  MPI_Finalize();  
  return 0;
}