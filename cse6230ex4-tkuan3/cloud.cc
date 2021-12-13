
#include <numeric>
#include "rand.h"
#include "timestep.h"
#include "state.h"
#include "json.h"
#include "params.h"

using namespace randgen;

int main (int argc, char **argv)
{
  Params params(argc, argv);
  Json log = Json();
  params.write(log);

  Rand gen(params.seed);
  auto L = params.L;
  State X0(params.Np, gen, -L/2, L/2);
  State X = X0; // copy constructor
  Timestepper ts(params, gen);

  // for each chunk of time steps,
  // evolve the state and measure the square mean
  // displacement
  //
  // this should grow proportionally to the time
  // and the ratio is the estimated diffusivity
  // of the gas
  std::vector<double> est_coeffs(0);
  for (size_t t = 0; t < params.Nt; t += params.Nc) {
    size_t n_steps = MIN(params.Nc, params.Nt - t);
    double T = n_steps * params.dt;
    ts.step(n_steps, params.dt, X);
    double mean_dist2 = X0.mean_square_distance(X);
    double est_coeff = mean_dist2 / T / 6.0;
    est_coeffs.push_back(est_coeff);
    X0.copy(X);
  }
  log.write_vec("estimated diffusivities", est_coeffs);
  double mean_est_coeff = std::accumulate(est_coeffs.begin(), est_coeffs.end(), 0.0) / est_coeffs.size();
  log.write("mean estimated diffusivity", mean_est_coeff);
  ts.write(log);

  return 0;

}
