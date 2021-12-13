#include "solver.h"

// Incrementally improve X as a solution to A X = F by performing a jacobi
// sweep, storing the results in Xnew
void Solver::_jacobi_update(GridFunction &X, GridFunction &Xnew) {
  double total_start = MPI_Wtime();

  double halo_start_start = MPI_Wtime();
  X.halo_start();
  double halo_start_end = MPI_Wtime();
  _time_halo_start += halo_start_end - halo_start_start;

  double halo_end_start = halo_start_end;
  X.halo_end();
  double halo_end_end = MPI_Wtime();
  _time_halo_end += halo_end_end - halo_end_start;

  double sweep_start = halo_end_end;

  auto local_dims = X.local_dims();
  for (ssize_t i = 0; i < (ssize_t) local_dims[0]; i++) {
    for (ssize_t j = 0; j < (ssize_t) local_dims[1]; j++) {
      #pragma omp simd
      for (ssize_t k = 0; k < (ssize_t) local_dims[2]; k++) {
        Xnew[{i,j,k}] = (1./6.) *
                        (
                         _F[{i,j,k}]
                         + X[{i-1,j,k}]
                         + X[{i+1,j,k}]
                         + X[{i,j-1,k}]
                         + X[{i,j+1,k}]
                         + X[{i,j,k-1}]
                         + X[{i,j,k+1}]
                        );
      }
    }
  }
  double sweep_end = MPI_Wtime();
  _sweep_time += sweep_end - sweep_start;

  double total_end = sweep_end;
  _sweep_count++;
  _total_time += total_end - total_start;
}

// Measure the residual: how far X is from solving the equations A X = F
double Solver::_residual_norm(const GridFunction &X) {
  double res2 = 0.;
  auto local_dims = X.local_dims();
  for (ssize_t i = 0; i < (ssize_t) local_dims[0]; i++) {
    for (ssize_t j = 0; j < (ssize_t) local_dims[1]; j++) {
      #pragma omp simd reduction(+:res2)
      for (ssize_t k = 0; k < (ssize_t) local_dims[2]; k++) {
        double diff = _F[{i,j,k}]
                      -6*X[{i,j,k}]
                      +  X[{i-1,j,k}]
                      +  X[{i+1,j,k}]
                      +  X[{i,j-1,k}]
                      +  X[{i,j+1,k}]
                      +  X[{i,j,k-1}]
                      +  X[{i,j,k+1}];
        res2 += diff*diff;
      }
    }
  }
  double global_res2;
  int mpierr = MPI_Allreduce(&res2, &global_res2, 1, MPI_DOUBLE, MPI_SUM, X.comm()); MPI_ECHK(mpierr);

  return sqrt(global_res2);
}

// Measure the error: how far X is from a known exact solution to A X = F
double Solver::_error_norm(const GridFunction &X) {
  double err2 = 0.;

  auto local_dims = X.local_dims();
  for (ssize_t i = 0; i < (ssize_t) local_dims[0]; i++) {
    for (ssize_t j = 0; j < (ssize_t) local_dims[1]; j++) {
      #pragma omp simd reduction(+:err2)
      for (ssize_t k = 0; k < (ssize_t) local_dims[2]; k++) {
        double diff = X[{i,j,k}] - _Xexact[{i,j,k}];
        err2 += diff * diff;
      }
    }
  }

  double global_err2;
  int mpierr = MPI_Allreduce(&err2, &global_err2, 1, MPI_DOUBLE, MPI_SUM, X.comm()); MPI_ECHK(mpierr);

  return sqrt(global_err2);
}
