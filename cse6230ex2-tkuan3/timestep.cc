
#include <cstddef>
#include <cmath>

namespace timestep {

  /// Update the state by one time step
  template<typename Real>
    void timestep(size_t n_asteroids,
                  Real mu, Real dt,
                  Real *__restrict__ x,
                  Real *__restrict__ y,
                  Real *__restrict__ z,
                  Real *__restrict__ u,
                  Real *__restrict__ v,
                  Real *__restrict__ w)
    {
      #pragma omp for nowait
      for (size_t i = 0; i < n_asteroids; i++) {
        Real r2 = x[i]*x[i] + y[i]*y[i] + z[i]*z[i];
        Real r = std::sqrt(r2);
        Real factor = (-mu * dt) / (r2 * r);

        u[i] += factor * x[i];
        v[i] += factor * y[i];
        w[i] += factor * z[i];
        x[i] += dt * u[i];
        y[i] += dt * v[i];
        z[i] += dt * w[i];
      }
    }

  /// Only update position by velocity
  template<typename Real>
    void position_timestep(size_t n_asteroids,
                           Real dt,
                           Real *__restrict__ x,
                           Real *__restrict__ y,
                           Real *__restrict__ z,
                           Real *__restrict__ u,
                           Real *__restrict__ v,
                           Real *__restrict__ w)
    {
      #pragma omp for nowait
      for (size_t i = 0; i < n_asteroids; i++) {
        x[i] += dt * u[i];
        y[i] += dt * v[i];
        z[i] += dt * w[i];
      }
    }

  template
    void timestep<double>(size_t n_asteroids,
                  double mu, double dt,
                  double *__restrict__ x,
                  double *__restrict__ y,
                  double *__restrict__ z,
                  double *__restrict__ u,
                  double *__restrict__ v,
                  double *__restrict__ w);

  template
    void timestep<float>(size_t n_asteroids,
                  float mu, float dt,
                  float *__restrict__ x,
                  float *__restrict__ y,
                  float *__restrict__ z,
                  float *__restrict__ u,
                  float *__restrict__ v,
                  float *__restrict__ w);

  template
    void position_timestep<double>(size_t n_asteroids,
                                   double dt,
                                   double *__restrict__ x,
                                   double *__restrict__ y,
                                   double *__restrict__ z,
                                   double *__restrict__ u,
                                   double *__restrict__ v,
                                   double *__restrict__ w);

  template
    void position_timestep<float>(size_t n_asteroids,
                                  float dt,
                                  float *__restrict__ x,
                                  float *__restrict__ y,
                                  float *__restrict__ z,
                                  float *__restrict__ u,
                                  float *__restrict__ v,
                                  float *__restrict__ w);
}
