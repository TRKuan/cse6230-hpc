#if !defined(UTIL_H)
#define UTIL_H

#include <vector>
#include <mpi.h>

#define MPI_CHK(err) if (err != MPI_SUCCESS) return err

int allreduce_vec(std::vector<double> &vec, MPI_Op op, double &result, MPI_Comm comm);
int allreduce_vec_reference(std::vector<double> &vec, MPI_Op op, double &result, MPI_Comm comm);
int inclusive_prefix_sum(std::vector<double> &vec, std::vector<double> &result, MPI_Comm comm);
int inclusive_prefix_sum_reference(std::vector<double> &vec, std::vector<double> &result, MPI_Comm comm);
int transpose_1d_col(int nrows_global, int ncols_global, std::vector<double> &A, std::vector<double> &AT, MPI_Comm comm);
int transpose_1d_col_reference(int nrows_global, int ncols_global, std::vector<double> &A, std::vector<double> &AT, MPI_Comm comm);
int transpose_2d(int nrows_local, int ncols_local, std::vector<double> &A, std::vector<double> &AT, MPI_Comm comm);
int transpose_2d_reference(int nrows_local, int ncols_local, std::vector<double> &A, std::vector<double> &AT, MPI_Comm comm);

int gather_vec(std::vector<double> &vec, std::vector<double> &gathered_vec, int root,  MPI_Comm comm);

int scatter_vec(std::vector<double> &gathered_vec, std::vector<double> &vec, int root,  MPI_Comm comm);
#endif
