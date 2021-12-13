#include <mpi.h>
#include <vector>
#include "util.h"
#include <iostream>

// Transpose a matrix that has a 2D distribution of the values
//
// input A: contains the local columns for this process: it has size
// [nrows_local x ncols_local] and is stored in column major order
//
// output AT: the transpose of A. contains the local columns for this process:
// it has size [ncols_local x nrows_local] and is stored in column major order
//
// You may assume that ncols_local is the same for each process and
// nrows_local is the same for each process, so each process has
// exactly the same number of input and output values
//
// You may assume comm is a Cartesian communicator, so that
// MPI_Cart_coords and MPI_Cart_rank can be used to convert ranks to/from the
// [block row, block column] index of the process's values in A and AT.
//
// You may assume the the Cartesian communicator is square:
// # of processes per row == # processses per column = sqrt(size)
int transpose_2d(int nrows_local, int ncols_local, std::vector<double> &A, std::vector<double> &AT, MPI_Comm comm)
{
  int err;

  int rank;
  int t_rank;
  int coords[2];
  err = MPI_Comm_rank(comm, &rank); MPI_CHK(err);
  err = MPI_Cart_coords(comm, rank, 2, coords); MPI_CHK(err);
  int t_coords[2] = {coords[1], coords[0]};
  err = MPI_Cart_rank(comm, t_coords, &t_rank); MPI_CHK(err);

  int block_size = nrows_local * ncols_local;
  std::vector<double> AT_col_major(block_size);
  if (coords[0] != coords[1]) {
    err = MPI_Send(&A[0], block_size, MPI_DOUBLE, t_rank, 0, comm); MPI_CHK(err);
    err = MPI_Recv(&AT_col_major[0], block_size, MPI_DOUBLE, t_rank, 0, comm, MPI_STATUS_IGNORE); MPI_CHK(err);
  } else {
    AT_col_major = A;
  }

  int t_nrows_local = ncols_local;
  int t_ncols_local = nrows_local;
  int pos = 0;
  for (int i = 0; i < t_nrows_local; i++) {
    for (int j = 0; j < t_ncols_local; j++) {
      AT[i + j*t_nrows_local] = AT_col_major[pos];
      pos += 1;
    }
  }

  return 0;
}

int transpose_2d_reference(int nrows_local, int ncols_local, std::vector<double> &A, std::vector<double> &AT, MPI_Comm comm)
{
  int err;
  int coords[2];

  int rank;
  err = MPI_Comm_rank(comm, &rank); MPI_CHK(err);

  err = MPI_Cart_coords(comm, rank, 2, coords); MPI_CHK(err);

  MPI_Comm row_comm; // color by row
  err = MPI_Comm_split(comm, coords[0], rank, &row_comm); MPI_CHK(err);

  MPI_Comm col_comm; // color by column
  err = MPI_Comm_split(comm, coords[1], rank, &col_comm); MPI_CHK(err);

  int row_size;
  err = MPI_Comm_size(row_comm, &row_size); MPI_CHK(err);
  int col_size;
  err = MPI_Comm_size(col_comm, &col_size); MPI_CHK(err);
  int row_rank;
  err = MPI_Comm_rank(row_comm, &row_rank); MPI_CHK(err);
  int col_rank;
  err = MPI_Comm_rank(col_comm, &col_rank); MPI_CHK(err);

  int ncols_global = ncols_local * row_size;
  int nrows_global = nrows_local * col_size;

  std::vector<double> A_row_gathered;
  err = gather_vec(A, A_row_gathered, 0, row_comm); MPI_CHK(err);

  std::vector<double> A_rowT_gathered(row_rank == 0 ? ncols_global * nrows_local : 0);
  std::vector<double> AT_row_gathered(row_rank == 0 ? nrows_global * ncols_local : 0);

  if (row_rank == 0) {

    for (int j = 0; j < nrows_local; j++) {
      for (int i = 0; i < ncols_global; i++) {
        A_rowT_gathered[i + j * ncols_global] = A_row_gathered[j + i * nrows_local];
      }
    }

    std::vector<double> A_gathered;
    err = gather_vec(A_rowT_gathered, A_gathered, 0, col_comm); MPI_CHK(err);

    std::vector<double> AT_gathered(col_rank == 0 ? A_gathered.size() : 0);
    if (col_rank == 0) {
      for (int j = 0; j < ncols_global; j++) {
        for (int i = 0; i < nrows_global; i++) {
          AT_gathered[i + j * nrows_global] = A_gathered[j + i * ncols_global];
        }
      }
    }

    std::vector<double> AT_rowT_gathered(nrows_global * ncols_local);
    scatter_vec(AT_gathered, AT_rowT_gathered, 0, col_comm); MPI_CHK(err);
    for (int j = 0; j < nrows_global; j++) {
      for (int i = 0; i < ncols_local; i++) {
        AT_row_gathered[i + j * ncols_local] = AT_rowT_gathered[j + i * nrows_global];
      }
    }
  }
  scatter_vec(AT_row_gathered, AT, 0, row_comm); MPI_CHK(err);

  return 0;
}
