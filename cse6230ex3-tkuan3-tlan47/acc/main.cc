#include <iostream>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include "state.h"
#include "timestep.h"

using namespace state;
using namespace timestep;

static void parse_args(int argc, char **argv, size_t &n_charges_per_dir, size_t &n_steps, std::string &precision)
{
  if (argc < 3 || argc > 4 ) {
    std::cout << "Usage: " << argv[0] << " {float,double} <T> <N>" << std::endl;
    std::cout << "     :" << std::endl;
    std::cout << "     : Perform an n-body simulation of charged particles trapped in a box" << std::endl;
    std::cout << "     : The particles are initially arranged in a regular N^3 grid" << std::endl;
    std::cout << "     :" << std::endl;
    std::cout << "     : {float,double} precision for arithmetic" << std::endl;
    std::cout << "     : <T> number of time steps to simulate" << std::endl;
    std::cout << "     : <N> number of charges per direction" << std::endl;
    std::exit(-1);
  }

  precision = std::string(argv[1]);
  if (precision != "double" && precision != "float") {
    std::cout << "unrecognized precision: " << precision << std::endl;
    std::exit(-1);
  }
  n_steps = std::strtoull(argv[2], NULL, 10);
  n_charges_per_dir = 0;
  if (argc == 4) {
    n_charges_per_dir = std::strtoull(argv[3], NULL, 10);
  }
}

template <typename Real>
  static void simulate(size_t n_charges_per_dir, size_t n_steps)
{
    State<Real> state(n_charges_per_dir);
    Real *x, *y, *z, *u, *v, *w;

    state.arrays(&x, &y, &z, &u, &v, &w);

    verlet_timesteps<Real>(n_steps, state.n_charges(), 1.e-7,
                           x, y, z,
                           u, v, w);
}

static void print_statistics(size_t n_charges_per_dir, size_t n_steps, std::string &precision, double seconds)
{
  size_t n_c = n_charges_per_dir;
  size_t n_charges = n_c * n_c * n_c;

  std::cout << "{" << std::endl;
  std::cout << "  \"precision\" : \"" << precision << "\"," << std::endl;
  std::cout << "  \"num charges per direction\" : \"" << n_c << "\"," << std::endl;
  std::cout << "  \"num charges\" : \"" << n_charges << "\"," << std::endl;
  std::cout << "  \"n_steps\" : \"" << n_steps << "\"," << std::endl;
  std::cout << "  \"runtime (seconds)\" : \"" << seconds << "\"," << std::endl;
  std::cout << "  \"charge interactions per second\" : \"" << (double) n_charges * (double) n_charges * (double) n_steps / seconds << "\"," << std::endl;
  std::cout << "}" << std::endl;
}

int main (int argc, char **argv)
{
  size_t n_charges_per_dir, n_steps;
  std::string precision;

  parse_args(argc, argv, n_charges_per_dir, n_steps, precision);

  auto start_time = std::chrono::high_resolution_clock::now();
  if (precision == "float") {
    simulate<float>(n_charges_per_dir, n_steps);
  } else {
    simulate<double>(n_charges_per_dir, n_steps);
  }
  auto end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> dur = end_time - start_time;
  double seconds = dur.count();

  print_statistics(n_charges_per_dir, n_steps, precision, seconds);

  return 0;
}

