#include <cmath>
#include <cstddef>
#include <limits>
#include "elements.h"
#include "state.h"

namespace elements {

  static const double pi_d = 3.141592653589793;

  /// Convert orbital elements to state vector: see
  //
  // https://downloads.rene-schwarz.com/download/M001-Keplerian_Orbit_Elements_to_Cartesian_State_Vectors.pdf
  template<typename Real>
    void elements_to_state_vector(Real mu,        // mass x gravity of origin (the sum)
                                  Real t,             // time since the epoch
                                  const Real elems[],
                                  Real state[])
    {
      const Real pi = (Real) pi_d;
      Real a, e, i, w, Node, M0;
      Real x, y, z, v_x, v_y, v_z;

      unpack(elems, a, e, i, w, Node, M0);
      // anomaly at epoch to anomaly at t (years)
      Real M = M0 + t * std::sqrt(mu / (a * a * a));
      // normalize to [0, 2 * pi)
      M -= std::floor(M / (2. * pi))* (2. * pi);
      // Newton's method to solve for eccentric anomaly
      Real E = M; Real sinE = std::sin(E); Real cosE = std::cos(E); Real g = (1. - e * cosE);
      for (size_t it = 0; it < 7; it++) {
        Real E_new = E - (E - e * sinE - M) / g;
        E = E_new; sinE = std::sin(E); cosE = std::cos(E); g = (1. - e * cosE);
      }
      // true anomaly
      Real nu = 2. * std::atan2(std::sqrt(1. + e) * std::sin(E / 2.), std::sqrt(1. - e) * std::cos(E / 2.));
      Real sin_nu = std::sin(nu); Real cos_nu = std::cos(nu);
      // distance from the sun
      Real r = a * g;
      // velocity scale
      Real scale = std::sqrt(mu * a) / r;
      // position and velocity within the plane with semimajor axis in the x
      // direction
      Real o_x = r * cos_nu; Real o_u = - (scale * sinE);
      Real o_y = r * sin_nu; Real o_v = scale * cosE * std::sqrt(1 - e * e);

      // rotate the semimajor axis
      Real sin_w = std::sin(w); Real cos_w = std::cos(w);
      Real w_x = cos_w * o_x - sin_w * o_y; Real w_u = cos_w * o_u - sin_w * o_v;
      Real w_y = sin_w * o_x + cos_w * o_y; Real w_v = sin_w * o_u + cos_w * o_v;

      // incline the plane
      Real sin_i = std::sin(i); Real cos_i = std::cos(i);
      Real i_x = w_x;         Real i_u = w_u;
      Real i_y = cos_i * w_y; Real i_v = cos_i * w_v;
      z        = sin_i * w_y;      v_z = sin_i * w_v;

      // rotate the ascending node
      Real sin_O = std::sin(Node); Real cos_O = std::cos(Node);
      x = cos_O * i_x - sin_O * i_y; v_x = cos_O * i_u - sin_O * i_v;
      y = sin_O * i_x + cos_O * i_y; v_y = sin_O * i_u + cos_O * i_v;
      state::pack(x, y, z, v_x, v_y, v_z, state);
    }

  template void elements_to_state_vector<double>(double, double, const double[], double[]);
  template void elements_to_state_vector<float>(float, float, const float[], float[]);

  /// Convert state vector to orbital elements: see
  //
  // https://downloads.rene-schwarz.com/download/M002-Cartesian_State_Vectors_to_Keplerian_Orbit_Elements.pdf
  template<typename Real>
    void state_vector_to_elements(Real mu_sun,
                                  Real t,
                                  const Real state[],
                                  Real elems[])
    {
      const Real pi = (Real) pi_d;
      const Real small_eps = 10. * std::numeric_limits<Real>::epsilon();
      Real a, e, i, w, Node, M0;
      Real x, y, z, v_x, v_y, v_z;

      state::unpack(state, x, y, z, v_x, v_y, v_z);
      // orbital momentum
      Real h_x = y * v_z - z * v_y;
      Real h_y = z * v_x - x * v_z;
      Real h_z = x * v_y - y * v_x;
      Real h_norm = std::sqrt(h_x * h_x + h_y * h_y + h_z * h_z);

      Real v2 = v_x * v_x + v_y * v_y + v_z * v_z;
      Real r = std::sqrt(x*x + y*y + z*z);
      Real v_r_dot = (x * v_x + y * v_y + z * v_z);

      // eccentricity vector
      Real e_c_r = v2 / mu_sun - 1. / r;
      Real e_c_v = v_r_dot / mu_sun;
      Real e_x = e_c_r * x - e_c_v * v_x;
      Real e_y = e_c_r * y - e_c_v * v_y;
      Real e_z = e_c_r * z - e_c_v * v_z;
      e = std::sqrt(e_x * e_x + e_y * e_y + e_z * e_z);

      // towards ascending node
      Real n_x = -h_y;
      Real n_y = h_x;
      Real n_z = 0.;
      Real n_norm = std::sqrt(n_x * n_x + n_y * n_y + n_z * n_z);

      // true anomaly
      Real nu;
      if (std::abs(r * e) <= small_eps * (e_x * x + e_y * y + e_z * z)) {
        // perfect circle, arbitrarily decide current position is 0 anomaly
        nu = 0.;
      } else {
        Real e_r_angle = (e_x * x + e_y * y + e_z * z) / (r * e);
        if (v_r_dot >= 0.) {
          nu = std::acos(e_r_angle);
        } else {
          nu = 2. * pi - std::acos(e_r_angle);
        }
      }

      // inclination
      if (h_norm == 0.) {
        // orbital plane is the ecliptic
        i = 0.;
      } else {
        i = std::acos(h_z / h_norm);
      }

      // eccentric anomaly
      Real E = 2. * std::atan(std::tan(nu / 2.) / sqrt((1. + e) / (1. - e)));

      // ascending node
      if (n_norm == 0.) {
        // by convention
        Node = 0.;
      } else {
        if (n_y >= 0.) {
          Node = std::acos(n_x / n_norm);
        } else {
          Node = 2. * pi - std::acos(n_x / n_norm);
        }
      }

      // periapsis
      if (std::abs(e * n_norm) < small_eps * (e_x * n_x + e_y * n_y + e_z * n_z)) {
        if (h_z >= 0.) {
          w = std::atan2(e_y, e_x);
        } else {
          w = 2. * pi - std::atan2(e_y, e_x);
        }
      } else {
        Real e_n_angle = (e_x * n_x + e_y * n_y + e_z * n_z) / (e * n_norm);
        if (e_z >= 0.) {
          w = std::acos(e_n_angle);
        } else {
          w = 2. * pi - std::acos(e_n_angle);
        }
      }

      // current mean anomaly
      Real M = E - e * std::sin(E);

      // semi-major axis
      a = 1. / ((2./r) - (v2 / mu_sun));

      // mean anomaly at epoc
      M0 = M - t * std::sqrt(mu_sun / (a * a * a));
      M0 -= std::floor(M0 / (2. * pi)) * (2. * pi);
      pack(a, e, i, w, Node, M0, elems);
    }

  template void state_vector_to_elements<double>(double, double, const double[], double[]);
  template void state_vector_to_elements<float>(float, float, const float[], float[]);
};
