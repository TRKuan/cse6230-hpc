#if !defined(SOLVER_H)
#define SOLVER_H

#define _USE_MATH_DEFINES
#include <cmath>
#include "grid.h"
#include "params.h"

class Solver {
  private:
    const GridFunction& _F;
    const GridFunction& _Xexact;
    GridFunction _Xtemp;
    double _T;
    double _S;
    double _time_halo_start;
    double _time_halo_end;
    double _sweep_time;
    double _total_time;
    double _sweep_count;
    Json& _outlog;

    double _residual_norm(const GridFunction &X);
    double _error_norm(const GridFunction &X);

    // not const on X because the halo values are updated in X
    // during this call
    void _jacobi_update(GridFunction &X, GridFunction &Xnew);

  public:
    Solver(Params &params, Json& outlog, const GridFunction& F, const GridFunction& Xexact):
      _F(F),
      _Xexact(Xexact),
      _Xtemp(F),
      _T(params.T),
      _S(params.S),
      _time_halo_start(0.),
      _time_halo_end(0.),
      _sweep_time(0.),
      _total_time(0.),
      _sweep_count(0.),
      _outlog(outlog) {
        double exact_rnorm = _residual_norm(_Xexact);
        _outlog.write("Exact solution residual", exact_rnorm);
      };

    ~Solver() {
      int size;

      MPI_Comm_size(_F.comm(), &size);
      double time_min[4] = {_time_halo_start, _time_halo_end, _sweep_time, _total_time};
      double time_max[4] = {_time_halo_start, _time_halo_end, _sweep_time, _total_time};
      double time_avg[4] = {_time_halo_start / size, _time_halo_end / size, _sweep_time / size, _total_time / size};
      double global_time_min[4];
      double global_time_max[4];
      double global_time_avg[4];

      MPI_Allreduce(&time_min[0], &global_time_min[0], 4, MPI_DOUBLE, MPI_MIN, _F.comm());
      MPI_Allreduce(&time_max[0], &global_time_max[0], 4, MPI_DOUBLE, MPI_MAX, _F.comm());
      MPI_Allreduce(&time_avg[0], &global_time_avg[0], 4, MPI_DOUBLE, MPI_SUM, _F.comm());
      _outlog.write("Minimum halo start time (seconds)", global_time_min[0]);
      _outlog.write("Minimum halo end time (seconds)", global_time_min[1]);
      _outlog.write("Minimum sweep time (seconds)", global_time_min[2]);
      _outlog.write("Minimum total jacobi time (seconds)", global_time_min[3]);
      _outlog.write("Maximum halo start time (seconds)", global_time_max[0]);
      _outlog.write("Maximum halo end time (seconds)", global_time_max[1]);
      _outlog.write("Maximum sweep time (seconds)", global_time_max[2]);
      _outlog.write("Maximum total jacobi time (seconds)", global_time_max[3]);
      _outlog.write("Mean halo start time (seconds)", global_time_avg[0]);
      _outlog.write("Mean halo end time (seconds)", global_time_avg[1]);
      _outlog.write("Mean sweep time (seconds)", global_time_avg[2]);
      _outlog.write("Mean total jacobi time (seconds)", global_time_avg[3]);
      auto global_dims = _F.global_dims();
      double total_lattice_updates = _sweep_count * ((double) global_dims[0] * (double) global_dims[1] * (double) global_dims[2]);
      double rate = total_lattice_updates / global_time_avg[3];
      _outlog.write("Mean lattice updates per second", rate);
    };

    GridFunction solve() {
      GridFunction Xout(_F);
      GridFunction *X = &Xout, *Xnew = &_Xtemp;

      std::vector<double> rs;
      std::vector<double> es;
      std::vector<double> cs;
      double prev_r = _residual_norm(*X);
      double prev_e = _error_norm(*X);
      rs.push_back(prev_r);
      es.push_back(prev_e);
      for (size_t i = 0, t = 0; t < _T; t++) {
        for (size_t s = 0; s < _S; s++, i++) {
          _jacobi_update(*X, *Xnew);
          std::swap(X, Xnew);
        }
        double r = _residual_norm(*X);
        double e = _error_norm(*X);
        double c = log(e / prev_e) / _S;
        rs.push_back(r);
        es.push_back(e);
        cs.push_back(c);
        prev_r = r;
        prev_e = e;
      }
      Xout = *X;
      _outlog.write_vec("residual history", rs);
      _outlog.write_vec("error history", es);
      _outlog.write_vec("convergence factor history", cs);
      return Xout;
    }
};

#endif
