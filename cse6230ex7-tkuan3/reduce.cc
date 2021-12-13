
#include <mpi.h>
#include "util.h"
#include <vector>
#include <iostream>

int allreduce_vec(std::vector<double> &vec, MPI_Op op, double &result, MPI_Comm comm)
{
  int err;

  double local_val = vec[0];
  for (int i = 1; i < vec.size(); i++) {
    err = MPI_Reduce_local(&vec[i], &local_val, 1, MPI_DOUBLE, op); MPI_CHK(err);
  }

  double val;
  err = MPI_Allreduce(&local_val, &val, 1, MPI_DOUBLE, op, comm); MPI_CHK(err);

  result = val;
  return 0;
}

int allreduce_vec_reference(std::vector<double> &vec, MPI_Op op, double &result, MPI_Comm comm)
{
  int err;

  std::vector<double> gathered_vec;
  err = gather_vec(vec, gathered_vec, 0, comm); MPI_CHK(err);

  int rank;
  err = MPI_Comm_rank(comm, &rank); MPI_CHK(err);

  double val = 0.;
  if (rank == 0) {
    int n = gathered_vec.size();
    if (n) {
      val = gathered_vec[0];
      for (int i = 1; i < n; i++) {
        err = MPI_Reduce_local(&gathered_vec[i], &val, 1, MPI_DOUBLE, op); MPI_CHK(err);
      }
    }
  }

  err = MPI_Bcast(&val, 1, MPI_DOUBLE, 0, comm); MPI_CHK(err);

  result = val;

  return 0;
}
