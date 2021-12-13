#if !defined(GRID_H)
#define GRID_H

#include <vector>
#include <mpi.h>
#include <cstddef>
#include <array>
#include <cmath>
#include "util.h"

class GridFunction {
  private:
    MPI_Comm _comm; // the Cartesian communicator for this grid

    std::array<size_t, 3> _global_dims;   // the total number of grid points
                                          // in each dimension: each grid
                                          // point is logically indexed by a
                                          // triplet in
                                          // [0, _global_dims[0]) x [0, _global_dims[1]) x [0, _global_dims[2])
                                          //
                                          // NOTE: the halo of boundary points
                                          // around the grid is not counted:
                                          // grid functions are always zero at
                                          // all of the boundary points, so
                                          // they aren't updated by any
                                          // process.

    std::array<size_t, 3> _local_dims;    // the number of **owned** grid points in each
                                          // dimension on this process's portion of the grid: each local grid
                                          // point is logically indexed by a
                                          // triplet in
                                          // [0, _local_dims[0]) x [0, _local_dims[1]) x [0, _local_dims[2])
                                          //
                                          // NOTE: each process will have a
                                          // grid with a halo of one grid
                                          // point in each direction,
                                          // so the number of grid points
                                          // *present* on each processor will
                                          // be _local_dims[d] + 2,
                                          // but the points in the halo region
                                          // are not owned by the local
                                          // process: they are either
                                          // duplicates of points on another
                                          // process or boundary points
                                          // (which aren't counted)

    std::array<size_t, 3> _local_offsets; // The grid point with local index {i, j, k} corresponds to the global
                                          // grid point with index
                                          // {i + _local_offsets[0], j + _local_offsets[1], k + _local_offsets[2]}

    std::vector<double> _data;            // one contiguous vector that
                                          // contains the portion of the 3D
                                          // grid function that is on this
                                          // process.  It should have size
                                          // { _local_dims[0] + 2, _local_dims[1] + 2, _local_dims[2] + 2 }
                                          // because it represent both the
                                          // grid values owned by the current
                                          // process and the grid values in the
                                          // halo region
    std::array<size_t, 3> _local_strides;

    bool _var; // whether the halo_start() and halo_end() methods should use the
               // variant algorithm or not
               // true: in place buffers with MPI vector types
               // false: packing and unpacking send and receive buffers
    
    std::array<size_t, 3> _grid_dims;
    std::vector<std::vector<std::vector<double>>> _buffer;
    MPI_Datatype _yz_vector;
    MPI_Datatype _xz_vector;
    MPI_Datatype _xy_vector;
    MPI_Datatype _y_vector;

  public:
    // accessors
    const MPI_Comm& comm() const { return _comm; }
    const std::array<size_t,3>& global_dims() const { return _global_dims; }
    const std::array<size_t,3>& local_dims() const { return _local_dims; }
    const std::array<size_t,3>& local_offsets() const { return _local_offsets; }

    // default constructor: from non-cartesian communicator and global grid size
    // creates a cartesian communicator and the data for the function
    GridFunction(MPI_Comm non_cart_comm, size_t n_global_x, size_t n_global_y, size_t n_global_z, bool var);

    // copy constructor: duplicate the cartesian communicator
    GridFunction(const GridFunction &grid):
      _global_dims(grid._global_dims),
      _local_dims(grid._local_dims),
      _local_offsets(grid._local_offsets),
      _data(grid._data),
      _local_strides(grid._local_strides),
      _var(grid._var),
      _grid_dims(grid._grid_dims),
      _buffer(grid._buffer),
      _yz_vector(grid._yz_vector),
      _xz_vector(grid._xz_vector),
      _xy_vector(grid._xy_vector),
      _y_vector(grid._y_vector)
    {
      int err;

      err = MPI_Comm_dup(grid._comm, &_comm); MPI_ECHK(err);
    };

    // assignment copies the data
    void operator = (const GridFunction &grid) {
      _data = grid._data;
    };

    // free the cartesian communicator
    ~GridFunction() {
      MPI_Comm_free(&_comm);
    }

    double& operator[](const ssize_t(& ijk)[3]) {
      return _data[(ijk[0] + 1) * _local_strides[0] + (ijk[1] + 1) *_local_strides[1] + (ijk[2] + 1)];
    }

    const double& operator[](const ssize_t(& ijk)[3]) const {
      return _data[(ijk[0] + 1) * _local_strides[0] + (ijk[1] + 1) *_local_strides[1] + (ijk[2] + 1)];
    }

    double norm() const {
      double local_norm2 = 0.;
      auto ldims = local_dims();
      for (ssize_t i = 0; i < (ssize_t) ldims[0]; i++) {
        for (ssize_t j = 0; j < (ssize_t) ldims[1]; j++) {
          for (ssize_t k = 0; k < (ssize_t) ldims[1]; k++) {
            double val = (*this)[{i,j,k}];
            local_norm2 += val * val;
          }
        }
      }
      double global_norm2;
      int err = MPI_Allreduce(&local_norm2, &global_norm2, 1, MPI_DOUBLE, MPI_SUM, _comm); MPI_ECHK(err);
      return sqrt(global_norm2);
    }

    // Begin the exchange of values
    // this function should be non-blocking
    // the data values in the GridFunction should not be modified
    // and the halo values are not guaranteed to be current when this function
    // returns
    void halo_start();

    // Complete the exchange of values
    // this functin should be blocking until:
    // - it is safe to modify the values in the GridFunction data without
    //   affecting the halo values for other processes
    // - the halo values can be read: they contain to up-to-date values
    //   from other processes
    void halo_end();
};



#endif
