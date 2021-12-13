#include <iostream>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include "elements.h"
#include "state.h"
#include "parse.h"
#include "timestep.h"

/// Convert elements to state vector arrays
template<typename Real>
  static void elems_to_state_arrays(size_t n_asteroids, double mu, double time, std::vector<double> elems,
                                    Real *__restrict__ x,
                                    Real *__restrict__ y,
                                    Real *__restrict__ z,
                                    Real *__restrict__ u,
                                    Real *__restrict__ v,
                                    Real *__restrict__ w)
{
  using namespace elements;
  size_t n_asteroids_true = elems.size() / NUM_EL;
  #pragma omp parallel for
  for (size_t i = 0; i < n_asteroids; i++) {
    // if more asteroids are requested than are in the file, duplicate some
    double *elem = &elems[NUM_EL * (i % n_asteroids_true)];
    double statev[NUM_EL];
    Real statev_r[NUM_EL];

    // Compute the state vectors in double precision
    elements_to_state_vector<double>(mu, time, elem, statev);

    // Convert double precision to the requested precision
    for (size_t j = 0; j < NUM_EL; j++) {
      statev_r[j] = (Real) statev[j];
    }

    // copy into the arrays
    state::unpack(statev_r,
                  x[i], y[i], z[i],
                  u[i], v[i], w[i]);
  }
}

/// Compute the error: the maximum root least squares error between the
// prediction and the analytically computed value
//
// An error that is on the order of 1 means that either the final position is
// off by 1 AU (the distance from the Earth to the sun), the velocity is
// off by 1 AU (about 16% the orbital velocity of the earth
template<typename Real>
  static double compute_error(size_t n_asteroids,
                              Real *__restrict__ x,
                              Real *__restrict__ y,
                              Real *__restrict__ z,
                              Real *__restrict__ u,
                              Real *__restrict__ v,
                              Real *__restrict__ w,
                              double *__restrict__ x_true,
                              double *__restrict__ y_true,
                              double *__restrict__ z_true,
                              double *__restrict__ u_true,
                              double *__restrict__ v_true,
                              double *__restrict__ w_true)
{
  double error = 0.;
  // The error is the largest square error between the computed final state
  // and the true final state
  #pragma omp parallel for reduction(max: error)
  for (size_t i = 0; i < n_asteroids; i++) {
    double err_x = (double) x[i] - x_true[i];
    double err_y = (double) y[i] - y_true[i];
    double err_z = (double) z[i] - z_true[i];
    double err_u = (double) u[i] - u_true[i];
    double err_v = (double) v[i] - v_true[i];
    double err_w = (double) w[i] - w_true[i];
    double diff = std::sqrt
                  (
                    err_x*err_x +
                    err_y*err_y +
                    err_z*err_z +
                    err_u*err_u +
                    err_v*err_v +
                    err_w*err_w
                  );
    error = std::max(error,diff);
  }
  return error;
}

/// Read asteroid data from file and simulate one year using classical mechanics on the state vector.
//
// Time the simulation, and compare the final state vector to one derived analytically from the
// orbital elements.
template<typename Real>
  static void simulate(size_t &n_asteroids, size_t n_steps, double &seconds, double &error)
{
  using namespace elements;
  const double mu = 39.47841760435743;; // 4 pi^2 [AU^3 / y^2], sun's gravitational constant

  std::string filename("elements.txt");
  // load orbital elements at time T=0
  auto elems = parse::asteroid_file<double>(filename);
  size_t n_asteroids_true = elems.size() / NUM_EL;
  if (n_asteroids == 0) {
    n_asteroids = n_asteroids_true;
  }

  // Allocate arrays: not using automatically managed memory
  // to make performance simpler to understand
  Real *x = new Real[n_asteroids];
  Real *y = new Real[n_asteroids];
  Real *z = new Real[n_asteroids];
  Real *u = new Real[n_asteroids];
  Real *v = new Real[n_asteroids];
  Real *w = new Real[n_asteroids];
  double *x_final = new double[n_asteroids];
  double *y_final = new double[n_asteroids];
  double *z_final = new double[n_asteroids];
  double *u_final = new double[n_asteroids];
  double *v_final = new double[n_asteroids];
  double *w_final = new double[n_asteroids];

  // convert from orbital elements to state vectors
  elems_to_state_arrays(n_asteroids, mu, 0., elems,
                        x, y, z, u, v, w);
  elems_to_state_arrays(n_asteroids, mu, n_steps > 0 ? 1. : 0., elems,
                        x_final, y_final, z_final,
                        u_final, v_final, w_final);

  // == simulation ==
  seconds = 0.;
  if (n_steps > 0) {
    Real dt = ((Real) 1.) / n_steps;
    const Real mu_r = (Real) mu;

    auto start_time = std::chrono::high_resolution_clock::now();
    #pragma omp parallel
    {
      // Move the asteroids half a time step by their initial velocity
      timestep::position_timestep(n_asteroids, (Real) 0.5 * dt,
                                  &x[0], &y[0], &z[0],
                                  &u[0], &v[0], &w[0]);
      // Leap-frog time stepping
      for (size_t i = 0; i < n_steps; i++) {
        timestep::timestep(n_asteroids, mu_r, dt,
                          &x[0], &y[0], &z[0],
                          &u[0], &v[0], &w[0]);
      }
      // Leap frog time stepping ends with the asteroids at time 1 + 0.5 * dt:
      // correct this by moving them backwards half a time step
      timestep::position_timestep(n_asteroids, (Real) -0.5 * dt,
                                  &x[0], &y[0], &z[0],
                                  &u[0], &v[0], &w[0]);
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> dur = end_time - start_time;

    seconds = dur.count();
  }

  error = compute_error(n_asteroids,
                        x, y, z, u, v, w,
                        x_final, y_final, z_final,
                        u_final, v_final, w_final);

  // clean up
  delete [] w_final;
  delete [] v_final;
  delete [] u_final;
  delete [] z_final;
  delete [] y_final;
  delete [] x_final;
  delete [] w;
  delete [] v;
  delete [] u;
  delete [] z;
  delete [] y;
  delete [] x;
}

static void parse_args(int argc, char **argv, size_t &n_asteroids, size_t &n_steps, std::string &precision)
{
  if (argc < 3 || argc > 4 ) {
    std::cout << "Usage: " << argv[0] << " {float,double} <T> [<N>]" << std::endl;
    std::cout << "     :" << std::endl;
    std::cout << "     : Simulate a year of small asteroid orbits from" << std::endl;
    std::cout << "     : <https://ssd.jpl.nasa.gov/?sb_elem> and print statistics" << std::endl;
    std::cout << "     : about the simulation runtime and its final accuracy" << std::endl;
    std::cout << "     :" << std::endl;
    std::cout << "     : {float,double} precision for arithmetic" << std::endl;
    std::cout << "     : <T> number of time steps to simulate 1 year" << std::endl;
    std::cout << "     : <N> (optional) number of asteroids to include in the simulation" << std::endl;
    std::exit(-1);
  }

  precision = std::string(argv[1]);
  if (precision != "double" && precision != "float") {
    std::cout << "unrecognized precision: " << precision << std::endl;
    std::exit(-1);
  }
  n_steps = std::strtoull(argv[2], NULL, 10);
  n_asteroids = 0;
  if (argc == 4) {
    n_asteroids = std::strtoull(argv[3], NULL, 10);
  }
}

int main(int argc, char **argv)
{
  size_t n_asteroids, n_steps;
  std::string precision;
  parse_args(argc, argv, n_asteroids, n_steps, precision);

  double seconds, error;
  if (precision == "double") {
    simulate<double>(n_asteroids, n_steps, seconds, error);
  } else {
    simulate<float>(n_asteroids, n_steps, seconds, error);
  }

  std::cout << "{" << std::endl;
  std::cout << "  \"precision\" : \"" << precision << "\"," << std::endl;
  std::cout << "  \"n_asteroids\" : \"" << n_asteroids << "\"," << std::endl;
  std::cout << "  \"n_steps\" : \"" << n_steps << "\"," << std::endl;
  std::cout << "  \"runtime (seconds)\" : \"" << seconds << "\"," << std::endl;
  std::cout << "  \"asteroid time steps per second\" : \"" << (double) n_steps * (double) n_asteroids / seconds << "\"," << std::endl;
  std::cout << "  \"simulated years per second\" : \"" << 1.0 / seconds << "\"," << std::endl;
  std::cout << "  \"error\" : \"" << error << "\"" << std::endl;
  std::cout << "}" << std::endl;

  return 0;
}
