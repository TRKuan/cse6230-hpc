#include <mpi.h>
#include <vector>
#include "util.h"

// Transpose a matrix that has a 1D distribution of the columns
//
// input A: contains the local columns for this process: it has size [nrows_global x ncols_local]
// and is stored in column major order
//
// output AT: contains the local columns for this process: it has size [ncols_global x nrows_local]
// and is stored in column major order
//
// You may assume that ncols_local is the same for each process and
// nrows_local is the same for each process, so each process has
// exactly the same number of input and output values
int transpose_1d_col(int nrows_global, int ncols_global, std::vector<double> &A, std::vector<double> &AT, MPI_Comm comm)
{
  int err;
  int size;

  err = MPI_Comm_size(comm, &size); MPI_CHK(err);

  int block_size = nrows_global * ncols_global / size;
  std::vector<double> A_row_major(block_size);
  int pos = 0;
  for (int i = 0; i < nrows_global; i++) {
    for (int j = 0; j < ncols_global / size; j++) {
      A_row_major[pos] = A[i + j*nrows_global];
      pos += 1;
    }
  }

  std::vector<double> AT_stack(block_size);
  err = MPI_Alltoall(&A_row_major[0], block_size / size, MPI_DOUBLE, &AT_stack[0], block_size / size, MPI_DOUBLE, comm); MPI_CHK(err);

  int t_ncols_global = nrows_global;
  int t_nrows_global = ncols_global;
  pos = 0;
  for (int k = 0; k < size; k++) {
    for (int j = 0; j < t_ncols_global / size; j++) {
      for (int i = k*t_nrows_global/size; i < (k+1)*t_nrows_global/size; i++) {
        AT[i + j*t_nrows_global] = AT_stack[pos];
        pos += 1;
      }
    }
  }

  return 0;
}

int transpose_1d_col_reference(int nrows_global, int ncols_global, std::vector<double> &A, std::vector<double> &AT, MPI_Comm comm)
{
  int err;

  std::vector<double> A_gathered;
  err = gather_vec(A, A_gathered, 0, comm); MPI_CHK(err);

  int rank;
  err = MPI_Comm_rank(comm, &rank); MPI_CHK(err);

  std::vector<double> AT_gathered(rank == 0 ? A_gathered.size() : 0);
  if (rank == 0) {
    for (int j = 0; j < nrows_global; j++) {
      for (int i = 0; i < ncols_global; i++) {
        AT_gathered[i + j*ncols_global] = A_gathered[j + i*nrows_global];
      }
    }
  }

  err = scatter_vec(AT_gathered, AT, 0, comm); MPI_CHK(err);
  return 0;
}
