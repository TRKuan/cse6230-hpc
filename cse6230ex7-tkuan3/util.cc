#include <mpi.h>
#include <vector>
#include "util.h"

int gather_vec(std::vector<double> &vec, std::vector<double> &gathered_vec, int root,  MPI_Comm comm)
{
  int err;
  int local_size = vec.size();

  int size;
  err = MPI_Comm_size(comm, &size); MPI_CHK(err);

  int rank;
  err = MPI_Comm_rank(comm, &rank); MPI_CHK(err);

  std::vector<int> recvcounts(rank == root ? size : 0);
  err = MPI_Gather(&local_size, 1, MPI_INT,
                   recvcounts.data(), 1, MPI_INT,
                   root, comm); MPI_CHK(err);

  std::vector<int> displs(rank == root ? (size + 1) : 0);
  if (rank == root) {
    displs[0] = 0;
    for (int i = 0; i < size; i++) {
      displs[i + 1] = displs[i] + recvcounts[i];
    }
  }

  gathered_vec = std::vector<double>(rank == root ? displs[size] : 0);

  err = MPI_Gatherv(vec.data(), local_size, MPI_DOUBLE,
                    gathered_vec.data(), recvcounts.data(), displs.data(), MPI_DOUBLE,
                    root, comm); MPI_CHK(err);
  return 0;
}

int scatter_vec(std::vector<double> &gathered_vec, std::vector<double> &vec, int root,  MPI_Comm comm)
{
  int err;
  int local_size = vec.size();

  int size;
  err = MPI_Comm_size(comm, &size); MPI_CHK(err);

  int rank;
  err = MPI_Comm_rank(comm, &rank); MPI_CHK(err);

  std::vector<int> recvcounts(rank == root ? size : 0);
  err = MPI_Gather(&local_size, 1, MPI_INT,
                   recvcounts.data(), 1, MPI_INT,
                   root, comm); MPI_CHK(err);

  std::vector<int> displs(rank == root ? (size + 1) : 0);
  if (rank == root) {
    displs[0] = 0;
    for (int i = 0; i < size; i++) {
      displs[i + 1] = displs[i] + recvcounts[i];
    }
  }

  err = MPI_Scatterv(gathered_vec.data(), recvcounts.data(), displs.data(), MPI_DOUBLE,
                     vec.data(), local_size, MPI_DOUBLE,
                     root, comm); MPI_CHK(err);
  return 0;
}
