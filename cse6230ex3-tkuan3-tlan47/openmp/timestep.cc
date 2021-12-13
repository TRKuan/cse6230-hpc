
#include <cstddef>
#include <cmath>
#include <algorithm>

namespace timestep {

  /// Only update velocity by n-body potential from other charges
  template<typename Real>
    static void velocity_timestep(size_t n_charges,
                                  Real dt,
                                  const Real *__restrict__ x,
                                  const Real *__restrict__ y,
                                  const Real *__restrict__ z,
                                  Real *__restrict__ u,
                                  Real *__restrict__ v,
                                  Real *__restrict__ w)
    {
      #pragma omp parallel for schedule (static)
      for (size_t i = 0; i < n_charges; i++) {
        Real d_u = 0;
        Real d_v = 0;
        Real d_w = 0;
        #pragma omp parallel for reduction(+:d_u, d_v, d_w) schedule (static)
        for (size_t j = 0; j < n_charges; j++) {
          if (i == j) {
            continue;
          }
          Real diff_x = x[i] - x[j];
          Real diff_y = y[i] - y[j];
          Real diff_z = z[i] - z[j];
          Real r2 = diff_x * diff_x + diff_y * diff_y + diff_z * diff_z;
          Real r = std::sqrt(r2);
          Real factor = dt / (r2 * r);

          d_u += factor * diff_x;
          d_v += factor * diff_y;
          d_w += factor * diff_z;
        }
        u[i] += d_u;
        v[i] += d_v;
        w[i] += d_w;
      }
      
    }

#define MAX(a,b) (a < b ? b: a)
#define MIN(a,b) (a < b ? a: b)

  /// Only update position by velocity
  template<typename Real>
    static void position_timestep(size_t n_charges,
                                  Real dt,
                                  Real *__restrict__ x,
                                  Real *__restrict__ y,
                                  Real *__restrict__ z,
                                  const Real *__restrict__ u,
                                  const Real *__restrict__ v,
                                  const Real *__restrict__ w)
    {
      #pragma omp for schedule (static)
      for (size_t i = 0; i < n_charges; i++) {
        x[i] += dt * u[i];
        x[i] = MAX(x[i], (Real) -1);
        x[i] = MIN(x[i], (Real) 1);
        y[i] += dt * v[i];
        y[i] = MAX(y[i], (Real) -1);
        y[i] = MIN(y[i], (Real) 1);
        z[i] += dt * w[i];
        z[i] = MAX(z[i], (Real) -1);
        z[i] = MIN(z[i], (Real) 1);
      }
    }

  /// Leap frog verlet timesteps
  template<typename Real>
    void verlet_timesteps(size_t n_steps,
                          size_t n_charges,
                          Real dt,
                          Real *__restrict__ x,
                          Real *__restrict__ y,
                          Real *__restrict__ z,
                          Real *__restrict__ u,
                          Real *__restrict__ v,
                          Real *__restrict__ w)
    {
      if (n_steps == 0) {
        return;
      }
      position_timestep<Real>(n_charges, 0.5 * dt,
                              x, y, z,
                              u, v, w);
      velocity_timestep<Real>(n_charges, dt,
                              x, y, z,
                              u, v, w);
      for (size_t i = 0; i < n_steps - 1; i++) {
        position_timestep<Real>(n_charges, dt,
                                x, y, z,
                                u, v, w);
        velocity_timestep<Real>(n_charges, dt,
                                x, y, z,
                                u, v, w);
      }
      position_timestep<Real>(n_charges, 0.5 * dt,
                              x, y, z,
                              u, v, w);
    }

  template
    void verlet_timesteps<double>(size_t n_steps,
                                  size_t n_charges,
                                  double dt,
                                  double *__restrict__ x,
                                  double *__restrict__ y,
                                  double *__restrict__ z,
                                  double *__restrict__ u,
                                  double *__restrict__ v,
                                  double *__restrict__ w);

  template
    void verlet_timesteps<float>(size_t n_steps,
                                  size_t n_charges,
                                  float dt,
                                  float *__restrict__ x,
                                  float *__restrict__ y,
                                  float *__restrict__ z,
                                  float *__restrict__ u,
                                  float *__restrict__ v,
                                  float *__restrict__ w);
}
