
#include <mpi.h>
#include <vector>
#include "util.h"

int inclusive_prefix_sum(std::vector<double> &vec, std::vector<double> &result, MPI_Comm comm)
{
  int err;

  int n = result.size();
  result[0] = vec[0];
  for (int i = 1; i < n; i++) {
    result[i] = vec[i] + result[i-1];
  }

  double last_sum = result[n-1];
  double prev_sum;
  err = MPI_Exscan(&last_sum, &prev_sum, 1, MPI_DOUBLE, MPI_SUM, comm);

  int rank;
  err = MPI_Comm_rank(comm, &rank); MPI_CHK(err);

  if (rank != 0) {
    for (int i = 0; i < n; i++) {
      result[i] = prev_sum + result[i];
    }
  }

  return 0;
}

int inclusive_prefix_sum_reference(std::vector<double> &vec, std::vector<double> &result, MPI_Comm comm)
{
  int err;

  std::vector<double> gathered_vec;
  err = gather_vec(vec, gathered_vec, 0, comm); MPI_CHK(err);

  int rank;
  err = MPI_Comm_rank(comm, &rank); MPI_CHK(err);

  std::vector<double> gathered_result(rank == 0 ? gathered_vec.size() : 0);
  if (rank == 0) {
    double sum = 0;
    int n = gathered_vec.size();
    for (int i = 0; i < n; i++) {
      sum += gathered_vec[i];
      gathered_result[i] = sum;
    }
  }

  err = scatter_vec(gathered_result, result, 0, comm); MPI_CHK(err);

  return 0;
}

