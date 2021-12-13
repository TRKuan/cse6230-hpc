#include <math.h>
#include <tictoc.h>
#include "cloud_util.h"
#include "timestep.h"

using namespace randgen;

Timestepper::Timestepper(Params &params, Rand &gen):
  _accel(params)
{
  _gen = &gen;
  _params = &params;
  _total_particle_steps = 0.;
  _total_time = 0.;
}

/// Take Nt timesteps from initialize position X
void Timestepper::step(size_t Nt, double dt, State &X) {
  // combine the Brownian noise coefficient and the timestep into a
  // noise timestep
  double dt_noise = sqrt(2. * _params->d * dt);

  // create a temporary vector for hold the velocities
  State U(X.Np());

  TicTocTimer timer = tic();
  for (size_t t = 0; t < Nt; t++) {
    // compute the velocities from the positions
    _accel.accelerate(X, U);
    // compute the positions from the velocities
    _stream_and_noise (dt, dt_noise, X, U);
  }
  _total_time += toc(&timer);
  _total_particle_steps += X.Np() * Nt;
}

void Timestepper::write(Json &log) {
  log.write("total timestepping time", _total_time);
  double total_rate = _total_particle_steps / _total_time;
  log.write("particle timesteps per second", total_rate);
  _accel.write(log);
}

/// Given the velocities and a Brownian noise time step, update the position
//  of each particle.  The Brownian noise multiplies a normally distributed
//  random update.
void Timestepper::_stream_and_noise(double dt_stream, double dt_noise, State &X, State &U) {
  size_t Np;
  double *_x, *_y, *_z, *_u, *_v, *_w;

  X.get_arrays(Np, _x, _y, _z);
  U.get_arrays(Np, _u, _v, _w);

  double  *x = _x;
  double  *y = _y;
  double  *z = _z;
  const double  *u = _u;
  const double  *v = _v;
  const double  *w = _w;

  auto tag = _gen->get_tag();

  // the random number generator generates 4 random numbers at a time,
  // so work in batches of 4 particles
  OMP_PARALLEL_FOR(schedule(static))
  for (size_t ib = 0; ib < (Np + 3) / 4; ib++) {
    size_t i = 4 * ib;
    double rval[3][4];

    for (size_t d = 0; d < 3; d++) {
      _gen->nrand_hash_4(tag, i, d, 0, &rval[d][0]);
    }

    if (i + 4 <= Np) {
      OMP_SIMD()
      for (size_t j = 0; j < 4; j++) {
        x[i + j] += dt_stream * u[i + j] + dt_noise * rval[0][j];
        y[i + j] += dt_stream * v[i + j] + dt_noise * rval[1][j];
        z[i + j] += dt_stream * w[i + j] + dt_noise * rval[2][j];
      }
    }
    else {
      for (size_t j = 0; j < Np - i; j++) {
        x[i + j] += dt_stream * u[i + j] + dt_noise * rval[0][j];
        y[i + j] += dt_stream * v[i + j] + dt_noise * rval[1][j];
        z[i + j] += dt_stream * w[i + j] + dt_noise * rval[2][j];
      }
    }
  }
}

