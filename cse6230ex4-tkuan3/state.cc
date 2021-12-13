#include "state.h"
#include "rand.h"

/// Create unitialized state
State::State(size_t Np) {
  _Np = Np;
  _x = new double[Np];
  _y = new double[Np];
  _z = new double[Np];
}

/// Copy existing state
State::State(State &s) {
  _Np = s._Np;
  _x = new double[_Np];
  _y = new double[_Np];
  _z = new double[_Np];

  OMP_PARALLEL_FOR(simd schedule(static))
  for (size_t i = 0; i < _Np; i++) {
    _x[i] = s._x[i];
    _y[i] = s._y[i];
    _z[i] = s._z[i];
  }
}

/// Initialize state with random values in [lo, hi]
State::State(size_t Np, randgen::Rand &gen, double lo, double hi) {
  _Np = Np;
  _x = new double[Np];
  _y = new double[Np];
  _z = new double[Np];
  auto tag = gen.get_tag();

  size_t Np_div_4 = Np / 4;
  size_t Np_mod_4 = Np % 4;
  size_t Np_even = Np - Np_mod_4;

  // the random number generator generates 4 random numbers at a time,
  // so work in batches of 4 particles
  OMP_PARALLEL_FOR(schedule(static))
  for (size_t ib = 0; ib < Np_div_4; ib ++) {
    size_t i = 4 * ib;
    double rx[4], ry[4], rz[4];

    gen.urand_hash_4(tag, i, 0, 0, rx, lo, hi);
    gen.urand_hash_4(tag, i, 1, 0, ry, lo, hi);
    gen.urand_hash_4(tag, i, 2, 0, rz, lo, hi);
    OMP_SIMD()
    for (size_t j = 0; j < 4; j++) {
      _x[i+j] = rx[j];
      _y[i+j] = ry[j];
      _z[i+j] = rz[j];
    }
  }
  if (Np_mod_4) {
    size_t i = Np_even;
    double rx[4], ry[4], rz[4];

    gen.urand_hash_4(tag, i, 0, 0, rx, lo, hi);
    gen.urand_hash_4(tag, i, 1, 0, ry, lo, hi);
    gen.urand_hash_4(tag, i, 2, 0, rz, lo, hi);
    for (size_t j = 0; j < Np - i; j++) {
      _x[i+j] = rx[j];
      _y[i+j] = ry[j];
      _z[i+j] = rz[j];
    }
  }

}

/// Copy without creating
void State::copy(State &s)
{
  OMP_PARALLEL_FOR(simd schedule(static))
  for (size_t i = 0; i < _Np; i++) {
    _x[i] = s._x[i];
    _y[i] = s._y[i];
    _z[i] = s._z[i];
  }

}

/// Average the square distance of each particle from location in state s
double State::mean_square_distance(State &s) {
  double sum = 0.;

  OMP_PARALLEL_FOR(schedule(static) reduction(+:sum))
  for (size_t p = 0; p < _Np; p++) {
    double dx = _x[p] - s._x[p];
    double dy = _y[p] - s._y[p];
    double dz = _z[p] - s._z[p];

    double R2 = dx*dx + dy*dy + dz*dz;
    sum += R2;
  }
  double mean = sum / _Np;
  return mean;
}

/// Get for routines that need raw arrays
void State::get_arrays(size_t &Np, double *&x, double *&y, double *&z)
{
  Np = _Np;
  x = _x;
  y = _y;
  z = _z;
}

/// Destroy
State::~State()
{
  delete [] _x;
  delete [] _y;
  delete [] _z;
}
