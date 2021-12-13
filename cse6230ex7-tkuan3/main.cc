#include <iostream>
#include <mpi.h>
#include <vector>
#include <cmath>
#include "util.h"

int test_reduce(MPI_Comm comm, int nvalues_per_proc, int &nfailed)
{
  int rank;
  int size;
  int err;

  err = MPI_Comm_size(comm, &size); MPI_CHK(err);
  err = MPI_Comm_rank(comm, &rank); MPI_CHK(err);
  std::vector<double> X(nvalues_per_proc + rank);

  int nvalues_prev = nvalues_per_proc * rank + (rank * (rank - 1)) / 2;
  for (int i = 0; i < nvalues_per_proc + rank; i++) {
    X[i] = (double) nvalues_prev + i + 1;
  }

  int nvalues_total = nvalues_per_proc * size + (size * (size - 1)) / 2;

  double out_expected = nvalues_total * (nvalues_total + 1) / 2.;
  double out;
  if (rank == 0) {
    std::cout << "Testing allreduce_vec_reference: sum" << std::endl;
  }
  err = allreduce_vec_reference(X, MPI_SUM, out, comm); MPI_CHK(err);
  if (rank == 0) {
    std::cout << "Expected: " << out_expected << ", got: " << out << std::endl;
  }
  if (rank == 0) {
    std::cout << "Testing allreduce_vec: sum" << std::endl;
  }
  out = 0;
  err = allreduce_vec(X, MPI_SUM, out, comm); MPI_CHK(err);
  if (rank == 0) {
    std::cout << "Expected: " << out_expected << ", got: " << out << std::endl;
  }
  if (fabs(out - out_expected) > 1.e-7 * out_expected) {
    if (rank == 0) {
      std::cout << "not ok" << std::endl;
    }
    nfailed++;
  } else {
    if (rank == 0) {
      std::cout << "ok" << std::endl;
    }
  }

  out_expected = nvalues_total;
  out = 0.;
  if (rank == 0) {
    std::cout << "Testing allreduce_vec_reference: max" << std::endl;
  }
  err = allreduce_vec_reference(X, MPI_MAX, out, comm); MPI_CHK(err);
  if (rank == 0) {
    std::cout << "Expected: " << out_expected << ", got: " << out << std::endl;
  }
  if (rank == 0) {
    std::cout << "Testing allreduce_vec: max" << std::endl;
  }
  out = 0;
  err = allreduce_vec(X, MPI_MAX, out, comm); MPI_CHK(err);
  if (rank == 0) {
    std::cout << "Expected: " << out_expected << ", got: " << out << std::endl;
  }
  if (out != out_expected) {
    if (rank == 0) {
      std::cout << "not ok" << std::endl;
    }
    nfailed++;
  } else {
    if (rank == 0) {
      std::cout << "ok" << std::endl;
    }
  }

  out_expected = 1.;
  out = 0.;
  if (rank == 0) {
    std::cout << "Testing allreduce_vec_reference: min" << std::endl;
  }
  err = allreduce_vec_reference(X, MPI_MIN, out, comm); MPI_CHK(err);
  if (rank == 0) {
    std::cout << "Expected: " << out_expected << ", got: " << out << std::endl;
  }
  if (rank == 0) {
    std::cout << "Testing allreduce_vec: min" << std::endl;
  }
  out = 0;
  err = allreduce_vec(X, MPI_MIN, out, comm); MPI_CHK(err);
  if (rank == 0) {
    std::cout << "Expected: " << out_expected << ", got: " << out << std::endl;
  }
  if (out != out_expected) {
    if (rank == 0) {
      std::cout << "not ok" << std::endl;
    }
    nfailed++;
  } else {
    if (rank == 0) {
      std::cout << "ok" << std::endl;
    }
  }
  return 0;
}

int test_scan(MPI_Comm comm, int nvalues_per_proc, int &nfailed)
{
  int rank;
  int size;
  int err;

  err = MPI_Comm_size(comm, &size); MPI_CHK(err);
  err = MPI_Comm_rank(comm, &rank); MPI_CHK(err);
  std::vector<double> X(nvalues_per_proc + rank);
  std::vector<double> Y(nvalues_per_proc + rank);
  std::vector<double> Yref(nvalues_per_proc + rank);
  std::vector<double> Yexpected(nvalues_per_proc + rank);

  for (int i = 0; i < nvalues_per_proc + rank; i++) {
    X[i] = (double) 1;
  }

  int nvalues_prev = nvalues_per_proc * rank + (rank * (rank - 1)) / 2;
  for (int i = 0; i < nvalues_per_proc + rank; i++) {
    Yexpected[i] = (double) nvalues_prev + i + 1;
  }

  if (rank == 0) {
    std::cout << "Testing inclusive_prefix_sum_reference" << std::endl;
  }
  err = inclusive_prefix_sum_reference(X, Yref, comm); MPI_CHK(err);
  double error = 0.;
  for (int i = 0; i < nvalues_per_proc + rank; i++) {
    error = std::max(fabs(Yref[i] - Yexpected[i]),error);
  }
  double global_error;
  err = MPI_Allreduce(&error, &global_error, 1, MPI_DOUBLE, MPI_MAX, comm); MPI_CHK(err);
  if (rank == 0) {
    std::cout << "Error: " << global_error << std::endl;
  }
  if (rank == 0) {
    std::cout << "Testing inclusive_prefix_sum" << std::endl;
  }
  err = inclusive_prefix_sum(X, Y, comm); MPI_CHK(err);
  error = 0.;
  for (int i = 0; i < nvalues_per_proc + rank; i++) {
    error = std::max(fabs(Y[i] - Yexpected[i]),error);
  }
  err = MPI_Allreduce(&error, &global_error, 1, MPI_DOUBLE, MPI_MAX, comm); MPI_CHK(err);
  if (rank == 0) {
    std::cout << "Error: " << global_error << std::endl;
  }
  if (error > 1.e-7) {
    if (rank == 0) {
      std::cout << "not ok" << std::endl;
    }
    nfailed++;
  } else {
    if (rank == 0) {
      std::cout << "ok" << std::endl;
    }
  }
  return 0;
}

int test_transpose_1d(MPI_Comm comm, int nrows_per_proc, int ncols_per_proc, int &nfailed)
{
  int rank;
  int size;
  int err;

  err = MPI_Comm_size(comm, &size); MPI_CHK(err);
  err = MPI_Comm_rank(comm, &rank); MPI_CHK(err);
  std::vector<double> A(ncols_per_proc * nrows_per_proc * size);
  std::vector<double> AT(ncols_per_proc * nrows_per_proc * size);
  std::vector<double> ATref(ncols_per_proc * nrows_per_proc * size);
  std::vector<double> ATexpected(ncols_per_proc * nrows_per_proc * size);
  int joffset = ncols_per_proc * rank;
  for (int j = 0; j < ncols_per_proc; j++) {
    for (int i = 0; i < nrows_per_proc * size; i++) {
      A[i + j*(nrows_per_proc*size)] = i + (j+joffset)*(nrows_per_proc*size);
    }
  }
  joffset = nrows_per_proc * rank;
  for (int j = 0; j < nrows_per_proc; j++) {
    for (int i = 0; i < ncols_per_proc * size; i++) {
      ATexpected[i + j*(ncols_per_proc*size)] = (j+joffset) + i*(nrows_per_proc*size);
    }
  }

  if (rank == 0) {
    std::cout << "Testing transpose_1d_col_reference" << std::endl;
  }
  err = transpose_1d_col_reference(nrows_per_proc*size, ncols_per_proc*size, A, ATref, comm); MPI_CHK(err);
  double error = 0.;
  for (int i = 0; i < ncols_per_proc*nrows_per_proc*size; i++) {
    error = std::max(fabs(ATref[i] - ATexpected[i]),error);
  }
  double global_error;
  err = MPI_Allreduce(&error, &global_error, 1, MPI_DOUBLE, MPI_MAX, comm); MPI_CHK(err);
  if (rank == 0) {
    std::cout << "Error: " << global_error << std::endl;
  }

  if (rank == 0) {
    std::cout << "Testing transpose_1d_col" << std::endl;
  }
  err = transpose_1d_col(nrows_per_proc*size, ncols_per_proc*size, A, AT, comm); MPI_CHK(err);
  error = 0.;
  for (int i = 0; i < ncols_per_proc*nrows_per_proc*size; i++) {
    error = std::max(fabs(AT[i] - ATexpected[i]),error);
  }
  err = MPI_Allreduce(&error, &global_error, 1, MPI_DOUBLE, MPI_MAX, comm); MPI_CHK(err);
  if (rank == 0) {
    std::cout << "Error: " << global_error << std::endl;
  }
  if (error > 0.) {
    if (rank == 0) {
      std::cout << "not ok" << std::endl;
    }
    nfailed++;
  } else {
    if (rank == 0) {
      std::cout << "ok" << std::endl;
    }
  }
  return 0;
}

int test_transpose_2d(MPI_Comm comm, int nrows_per_proc, int ncols_per_proc, int &nfailed)
{
  int rank;
  int size;
  int err;
  int dims[2] = {0, 0};
  int periods[2] = {0, 0};
  int coords[2];

  err = MPI_Comm_size(comm, &size); MPI_CHK(err);

  MPI_Comm cart_comm;
  err = MPI_Dims_create(size, 2, dims); MPI_CHK(err);
  err = MPI_Cart_create(comm, 2, dims, periods, 0, &cart_comm); MPI_CHK(err);

  err = MPI_Comm_rank(cart_comm, &rank); MPI_CHK(err);
  err = MPI_Cart_coords(cart_comm, rank, 2, coords); MPI_CHK(err);

  std::vector<double> A(ncols_per_proc * nrows_per_proc);
  std::vector<double> AT(ncols_per_proc * nrows_per_proc);
  std::vector<double> ATref(ncols_per_proc * nrows_per_proc);
  std::vector<double> ATexpected(ncols_per_proc * nrows_per_proc);
  int ioffset = nrows_per_proc * coords[0];
  int joffset = ncols_per_proc * coords[1];
  for (int j = 0; j < ncols_per_proc; j++) {
    for (int i = 0; i < nrows_per_proc; i++) {
      A[i + j*nrows_per_proc] = i+ioffset + (j+joffset)*(nrows_per_proc*dims[0]);
    }
  }
  ioffset = ncols_per_proc * coords[0];
  joffset = nrows_per_proc * coords[1];
  for (int j = 0; j < nrows_per_proc; j++) {
    for (int i = 0; i < ncols_per_proc; i++) {
      ATexpected[i + j*ncols_per_proc] = (j+joffset) + (i+ioffset)*(nrows_per_proc*dims[0]);
    }
  }

  if (rank == 0) {
    std::cout << "Testing transpose_2d_reference" << std::endl;
  }
  err = transpose_2d_reference(nrows_per_proc, ncols_per_proc, A, ATref, cart_comm); MPI_CHK(err);
  double error = 0.;
  for (int i = 0; i < ncols_per_proc*nrows_per_proc; i++) {
    error = std::max(fabs(ATref[i] - ATexpected[i]),error);
  }
  double global_error;
  err = MPI_Allreduce(&error, &global_error, 1, MPI_DOUBLE, MPI_MAX, comm); MPI_CHK(err);
  if (rank == 0) {
    std::cout << "Error: " << global_error << std::endl;
  }

  if (rank == 0) {
    std::cout << "Testing transpose_2d" << std::endl;
  }
  err = transpose_2d(nrows_per_proc, ncols_per_proc, A, AT, cart_comm); MPI_CHK(err);
  error = 0.;
  for (int i = 0; i < ncols_per_proc*nrows_per_proc; i++) {
    error = std::max(fabs(AT[i] - ATexpected[i]),error);
  }
  err = MPI_Allreduce(&error, &global_error, 1, MPI_DOUBLE, MPI_MAX, comm); MPI_CHK(err);
  if (rank == 0) {
    std::cout << "Error: " << global_error << std::endl;
  }
  if (error > 0.) {
    if (rank == 0) {
      std::cout << "not ok" << std::endl;
    }
    nfailed++;
  } else {
    if (rank == 0) {
      std::cout << "ok" << std::endl;
    }
  }
  err = MPI_Comm_free(&cart_comm); MPI_CHK(err);
  return 0;
}

int main(int argc, char **argv)
{
  int err;
  int nfailed = 0;
  err = MPI_Init(&argc, &argv); MPI_CHK(err);

  err = test_reduce(MPI_COMM_WORLD, 13, nfailed); MPI_CHK(err);
  err = test_scan(MPI_COMM_WORLD, 17, nfailed); MPI_CHK(err);
  err = test_transpose_1d(MPI_COMM_WORLD, 19, 23, nfailed); MPI_CHK(err);
  err = test_transpose_2d(MPI_COMM_WORLD, 29, 31, nfailed); MPI_CHK(err);

  if (nfailed > 0) {
    int rank;

    err = MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_CHK(err);
    if (!rank) {
      std::cout << "Total failed tests: " << nfailed << std::endl;
    }
  }

  err = MPI_Finalize();
  return err;
}
