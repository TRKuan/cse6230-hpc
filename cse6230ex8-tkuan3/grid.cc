#include "grid.h"
#include "util.h"
#include <iostream>

GridFunction::GridFunction(
  MPI_Comm non_cart_comm,
  size_t n_global_x,
  size_t n_global_y,
  size_t n_global_z,
  bool var)
{
  _var = var;
  // 1. construct a Cartesian communicator (the grid is not periodic)
  int err, size, rank, coords[3] = {0, 0, 0}, dims[3] = {0, 0, 0}, periods[3] = {0, 0, 0};

  err = MPI_Comm_size(non_cart_comm, &size); MPI_ECHK(err);
  err = MPI_Dims_create(size, 3, dims); MPI_ECHK(err);
  err = MPI_Cart_create(non_cart_comm, 3, dims, periods, 0, &_comm); MPI_ECHK(err);
  err = MPI_Comm_rank(_comm, &rank); MPI_ECHK(err);
  err = MPI_Cart_coords(_comm, rank, 3, coords); MPI_ECHK(err);

  _global_dims = std::array<size_t, 3> ({n_global_x, n_global_y, n_global_z});
  _grid_dims = std::array<size_t, 3> ({dims[0], dims[1], dims[2]});

  // 2. determine the local size and the local offsets
  for (size_t d = 0; d < 3; d++) {
    size_t n_global = _global_dims[d];
    size_t coord = coords[d];
    size_t dim_size = dims[d];

    size_t rem = n_global % dim_size;
    size_t div = n_global / dim_size;
    auto offset = [div, rem, dim_size] (size_t c) -> size_t {return c * div + (rem * c) / dim_size;};
    size_t this_offset = offset(coord);
    size_t next_offset = offset(coord+1);
    _local_offsets[d] = this_offset;
    _local_dims[d] = (next_offset - this_offset);
  }

  // 3. allocate an array
  size_t local_grid_size = 1;
  for (auto d: _local_dims) local_grid_size *= (d + 2);
  _data = std::vector<double>(local_grid_size);
  for (auto &v: _data) v = 0.;

  // 4. compute the strides
  _local_strides[2] = 1;
  _local_strides[1] = _local_strides[2] * (_local_dims[2] + 2);
  _local_strides[0] = _local_strides[1] * (_local_dims[1] + 2);

  // 5. Buffer
  if (!var) {
    _buffer = std::vector<std::vector<std::vector<double>>>(3);
    for (int coord_dim = 0; coord_dim < 3; coord_dim++) {
      int buffer_size = 1;
      for (int d = 0; d < 3; d++) {
        if (d != coord_dim) {
          buffer_size *= _local_dims[d];
        }
      }

      _buffer[coord_dim] = std::vector<std::vector<double>>(2, std::vector<double>(buffer_size));
    }
  }

  //6. Vector
  if (_var) {
    MPI_Type_vector(_local_dims[1], _local_dims[2], _local_strides[1], MPI_DOUBLE, &_yz_vector);
    MPI_Type_commit(&_yz_vector);

    MPI_Type_vector(_local_dims[0], _local_dims[2], _local_strides[0], MPI_DOUBLE, &_xz_vector);
    MPI_Type_commit(&_xz_vector);

    MPI_Type_vector(_local_dims[1], 1, _local_strides[1], MPI_DOUBLE, &_y_vector);
    MPI_Type_commit(&_y_vector);
    MPI_Type_hvector(_local_dims[0], 1, _local_strides[0] * sizeof(double), _y_vector, &_xy_vector);
    MPI_Type_commit(&_xy_vector);
  }
}

// Begin the exchange of values
// this function should be non-blocking
// the data values in the GridFunction should not be modified
// and the halo values are not guaranteed to be current when this function
// returns
void GridFunction::halo_start() {
  int err;
  int rank;
  int rank_dest;
  MPI_Request request;
  int coords[3];

  err = MPI_Comm_rank(_comm, &rank); MPI_ECHK(err);
  err = MPI_Cart_coords(_comm, rank, 3, coords); MPI_ECHK(err);

  for (int coord_dim = 0; coord_dim < 3; coord_dim++) {
    for (int dir: {1, -1}) {
      if (coords[coord_dim] + dir < 0 || coords[coord_dim] + dir >= _grid_dims[coord_dim]) {
        continue;
      }

      coords[coord_dim] += dir;
      err = MPI_Cart_rank(_comm, coords, &rank_dest); MPI_ECHK(err);
      coords[coord_dim] -= dir;

      int buffer_size = 1;
      for (int d = 0; d < 3; d++) {
        if (d != coord_dim) {
          buffer_size *= _local_dims[d];
        }
      }

      if (_var) {
        if (coord_dim == 0) {
          int start_idx;
          if (dir == -1) {
            start_idx = _local_strides[0] + _local_strides[1] + _local_strides[2];
          } else {
            start_idx = _local_dims[0] * _local_strides[0] + _local_strides[1] + _local_strides[2];
          }
          err = MPI_Isend(&_data[start_idx], 1, _yz_vector, rank_dest, 0, _comm, &request); MPI_ECHK(err);
        } else if (coord_dim == 1) {
          int start_idx;
          if (dir == -1) {
            start_idx = _local_strides[0] + _local_strides[1] + _local_strides[2];
          } else {
            start_idx = _local_strides[0] + _local_dims[1] * _local_strides[1] + _local_strides[2];
          }
          err = MPI_Isend(&_data[start_idx], 1, _xz_vector, rank_dest, 0, _comm, &request); MPI_ECHK(err);
        } else {
          int start_idx;
          if (dir == -1) {
            start_idx = _local_strides[0] + _local_strides[1] + _local_strides[2];
          } else {
            start_idx = _local_strides[0] + _local_strides[1] + _local_dims[2] * _local_strides[2];
          }
          err = MPI_Isend(&_data[start_idx], 1, _xy_vector, rank_dest, 0, _comm, &request); MPI_ECHK(err);
        }
      } else {
        int dim_major;
        int dim_minor;
        if (coord_dim == 0) {
          dim_major = 1;
          dim_minor = 2;
        } else if (coord_dim == 1) {
          dim_major = 0;
          dim_minor = 2;
        } else {
          dim_major = 0;
          dim_minor = 1;
        }

        ssize_t ijk[3] = {0, 0, 0};
        ijk[coord_dim] = dir == -1 ? 0 : _local_dims[coord_dim] - 1;

        for (int n = 0; n < buffer_size; n++) {
          _buffer[coord_dim][dir == -1?0:1][n] = (*this)[ijk];

          ijk[dim_minor] += 1;
          if (ijk[dim_minor] >= _local_dims[dim_minor]) {
            ijk[dim_minor] = 0;
            ijk[dim_major] += 1;
            if (ijk[dim_major] >= _local_dims[dim_major]) {
              break;
            }
          }
        }

        err = MPI_Isend(&_buffer[coord_dim][dir == -1?0:1][0], buffer_size, MPI_DOUBLE, rank_dest, 0, _comm, &request); MPI_ECHK(err);
      }
    }
  }
}

// Complete the exchange of values
// this function should be blocking until:
// - it is safe to modify the values in the GridFunction data without
//   affecting the halo values for other processes
// - the halo values can be read: they contain to up-to-date values
//   from other processes
void GridFunction::halo_end() {
  int err;
  int rank;
  int rank_source;
  int coords[3];

  err = MPI_Comm_rank(_comm, &rank); MPI_ECHK(err);
  err = MPI_Cart_coords(_comm, rank, 3, coords); MPI_ECHK(err);

  for (int coord_dim = 0; coord_dim < 3; coord_dim++) {
    for (int dir: {1, -1}) {
      if (coords[coord_dim] + dir < 0 || coords[coord_dim] + dir >= _grid_dims[coord_dim]) {
        continue;
      }

      coords[coord_dim] += dir;
      err = MPI_Cart_rank(_comm, coords, &rank_source); MPI_ECHK(err);
      coords[coord_dim] -= dir;

      int buffer_size = 1;
      for (int d = 0; d < 3; d++) {
        if (d != coord_dim) {
          buffer_size *= _local_dims[d];
        }
      }

      if (_var) {
        if (coord_dim == 0) {
          int start_idx;
          if (dir == -1) {
            start_idx = _local_strides[1] + _local_strides[2];
          } else {
            start_idx = (_local_dims[0] + 1) * _local_strides[0] + _local_strides[1] + _local_strides[2];
          }
          err = MPI_Recv(&_data[start_idx], 1, _yz_vector, rank_source, 0, _comm, MPI_STATUS_IGNORE); MPI_ECHK(err);
        } else if (coord_dim == 1) {
          int start_idx;
          if (dir == -1) {
            start_idx = _local_strides[0] + _local_strides[2];
          } else {
            start_idx = _local_strides[0] + (_local_dims[1] + 1) * _local_strides[1] + _local_strides[2];
          }
          err = MPI_Recv(&_data[start_idx], 1, _xz_vector, rank_source, 0, _comm, MPI_STATUS_IGNORE); MPI_ECHK(err);
        } else {
          int start_idx;
          if (dir == -1) {
            start_idx = _local_strides[0] + _local_strides[1];
          } else {
            start_idx = _local_strides[0] + _local_strides[1] + (_local_dims[2] + 1) * _local_strides[2];
          }
          err = MPI_Recv(&_data[start_idx], 1, _xy_vector, rank_source, 0, _comm, MPI_STATUS_IGNORE); MPI_ECHK(err);
        }
      } else {
        err = MPI_Recv(&_buffer[coord_dim][dir == -1?0:1][0], buffer_size, MPI_DOUBLE, rank_source, 0, _comm, MPI_STATUS_IGNORE); MPI_ECHK(err);

        int dim_major;
        int dim_minor;
        if (coord_dim == 0) {
          dim_major = 1;
          dim_minor = 2;
        } else if (coord_dim == 1) {
          dim_major = 0;
          dim_minor = 2;
        } else {
          dim_major = 0;
          dim_minor = 1;
        }

        ssize_t ijk[3] = {0, 0, 0};
        ijk[coord_dim] = dir == -1 ? -1 : _local_dims[coord_dim];

        for (int n = 0; n < buffer_size; n++) {
          (*this)[ijk] = _buffer[coord_dim][dir == -1?0:1][n];

          ijk[dim_minor] += 1;
          if (ijk[dim_minor] >= _local_dims[dim_minor]) {
            ijk[dim_minor] = 0;
            ijk[dim_major] += 1;
            if (ijk[dim_major] >= _local_dims[dim_major]) {
              break;
            }
          }
        }
      }
    }
  }
}
