#include <cstddef>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string_view>
#include <cmath>

using std::size_t;

/// Orbital elements
enum {
  AX, /// [AU ] semi-major axis
  EC, /// [1  ] eccentricity
  IN, /// [rad] inclination
  PE, /// [rad] argument of periapsis
  NO, /// [rad] longitude of ascending node
  ME, /// [rad] mean anomaly
  NUM_EL};

#if !defined(BATCH)
#define BATCH 1
#endif

enum{NAME_LEN=17,NAME_STRIDE};

const double pi = 3.141592653589793;
const double deg_to_rad = 3.141592653589793 / 180.;
const double mu_sun = 39.47841760435743;; // 4 pi^2 [AU^3 / y^2], sun's gravitational constant
const double small_eps = 1.e-15;

static inline void elements_to_cartesian(double a, double e, double i, double w, double Node, double M0, double t, double &x, double &y, double &z, double &v_x, double &v_y, double &v_z)
{
  // anomaly at epoch to anomaly at t (years)
  double M = M0 + t * std::sqrt(mu_sun / (a * a * a));
  // normalize to [0, 2 * pi)
  M -= std::floor(M / (2. * pi))* (2. * pi);
  // Newton's method to solve for eccentric anomaly
  double E = M; double sinE = std::sin(E); double cosE = std::cos(E); double g = (1. - e * cosE);
  for (size_t it = 0; it < 7; it++) {
    double E_new = E - (E - e * sinE - M) / g;
    E = E_new; sinE = std::sin(E); cosE = std::cos(E); g = (1. - e * cosE);
  }
  // true anomaly
  double nu = 2. * std::atan2(std::sqrt(1. + e) * std::sin(E / 2.), std::sqrt(1. - e) * std::cos(E / 2.));
  double sin_nu = std::sin(nu); double cos_nu = std::cos(nu);
  // distance from the sun
  double r = a * g;
  // velocity scale
  double scale = std::sqrt(mu_sun * a) / r;
  // position and velocity within the plane with semimajor axis in the x
  // direction
  double o_x = r * cos_nu; double o_u = - (scale * sinE);
  double o_y = r * sin_nu; double o_v = scale * cosE * std::sqrt(1 - e * e);

  // rotate the semimajor axis
  double sin_w = std::sin(w); double cos_w = std::cos(w);
  double w_x = cos_w * o_x - sin_w * o_y; double w_u = cos_w * o_u - sin_w * o_v;
  double w_y = sin_w * o_x + cos_w * o_y; double w_v = sin_w * o_u + cos_w * o_v;

  // incline the plane
  double sin_i = std::sin(i); double cos_i = std::cos(i);
  double i_x = w_x;         double i_u = w_u;
  double i_y = cos_i * w_y; double i_v = cos_i * w_v;
  z          = sin_i * w_y;        v_z = sin_i * w_v;

  // rotate the ascending node
  double sin_O = std::sin(Node); double cos_O = std::cos(Node);
  x = cos_O * i_x - sin_O * i_y; v_x = cos_O * i_u - sin_O * i_v;
  y = sin_O * i_x + cos_O * i_y; v_y = sin_O * i_u + cos_O * i_v;
}

static inline void cartesian_to_elements(double x, double y, double z, double v_x, double v_y, double v_z, double t, double &a, double &e, double &i, double &w, double &Node, double &M0)
{
  // orbital momentum
  double h_x = y * v_z - z * v_y;
  double h_y = z * v_x - x * v_z;
  double h_z = x * v_y - y * v_x;
  double h_norm = std::sqrt(h_x * h_x + h_y * h_y + h_z * h_z);

  double v2 = v_x * v_x + v_y * v_y + v_z * v_z;
  double r = std::sqrt(x*x + y*y + z*z);
  double v_r_dot = (x * v_x + y * v_y + z * v_z);

  // eccentricity vector
  double e_c_r = v2 / mu_sun - 1. / r;
  double e_c_v = v_r_dot / mu_sun;
  double e_x = e_c_r * x - e_c_v * v_x;
  double e_y = e_c_r * y - e_c_v * v_y;
  double e_z = e_c_r * z - e_c_v * v_z;
  e = std::sqrt(e_x * e_x + e_y * e_y + e_z * e_z);

  // towards ascending node
  double n_x = -h_y;
  double n_y = h_x;
  double n_z = 0.;
  double n_norm = std::sqrt(n_x * n_x + n_y * n_y + n_z * n_z);

  // true anomaly
  double nu;
  if (std::abs(r * e) <= small_eps * (e_x * x + e_y * y + e_z * z)) {
    // perfect circle, arbitrarily decide current position is 0 anomaly
    nu = 0.;
  } else {
    double e_r_angle = (e_x * x + e_y * y + e_z * z) / (r * e);
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
  double E = 2. * std::atan(std::tan(nu / 2.) / sqrt((1. + e) / (1. - e)));

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
    double e_n_angle = (e_x * n_x + e_y * n_y + e_z * n_z) / (e * n_norm);
    if (e_z >= 0.) {
      w = std::acos(e_n_angle);
    } else {
      w = 2. * pi - std::acos(e_n_angle);
    }
  }

  // current mean anomaly
  double M = E - e * std::sin(E);

  // semi-major axis
  a = 1. / ((2./r) - (v2 / mu_sun));

  // mean anomaly at epoc
  M0 = M - t * std::sqrt(mu_sun / (a * a * a));
  M0 -= std::floor(M0 / (2. * pi)) * (2. * pi);
}

void verlet_step(float dt,
                 float &__restrict__ x, float &__restrict__ y, float &__restrict__ z,
                 float &__restrict__ u, float &__restrict__ v, float &__restrict__ w)
{
  x += 0.5f * dt * u;
  y += 0.5f * dt * v;
  z += 0.5f * dt * w;
  float r2 = x*x + y*y + z*z;
  float ir = 1.0f / sqrtf(r2);
  float ir2 = 1.0f / r2;
  float factor = (- dt * ((float) mu_sun)) * ir * ir2;
  u += x * factor;
  v += y * factor;
  w += z * factor;
  x += 0.5f * dt * u;
  y += 0.5f * dt * v;
  z += 0.5f * dt * w;
}

void verlet(double dt, size_t n_steps,
                   double x, double y, double z,
                   double u, double v, double w,
                   double &x_final, double &y_final, double &z_final,
                   double &u_final, double &v_final, double &w_final)
{
  for (size_t k = 0; k < n_steps; k++) {
    x += 0.5 * dt * u;
    y += 0.5 * dt * v;
    z += 0.5 * dt * w;
    double r2 = x*x + y*y + z*z;
    double ir = 1. / sqrtf(r2);
    double ir2 = 1. / r2;
    double factor = (- dt * mu_sun) * ir * ir2;
    u += x * factor;
    v += y * factor;
    w += z * factor;
    x += 0.5 * dt * u;
    y += 0.5 * dt * v;
    z += 0.5 * dt * w;
  }
  x_final = x;
  y_final = y;
  z_final = z;
  u_final = u;
  v_final = v;
  w_final = w;
}

void verlet_step_array(float dt, size_t n_asteroids,
                       float *__restrict__ x,
                       float *__restrict__ y,
                       float *__restrict__ z,
                       float *__restrict__ u,
                       float *__restrict__ v,
                       float *__restrict__ w)
{
  for (size_t k = 0; k < n_asteroids; k++) {
    verlet_step(dt,
                x[k], y[k], z[k],
                u[k], v[k], w[k]);
  }
}

template <size_t B = BATCH>
class Elements {
  public:
    size_t num_bodies;

    ~Elements() {
      delete[] _elements;
      delete[] _name_buf;
      delete[] _names;
    }

    static Elements<B> from_numbered_asteroid_format(std::ifstream &ist, bool test = false)
    {
      std::string name_buf;
      std::string line;
      std::vector<double> elems;
      size_t n_bodies = 0;

      std::getline(ist, line); // skip field names
      std::getline(ist, line); // skip ----- -----
      while (std::getline(ist, line)) {
        char name[NAME_STRIDE] = {'\0'};
        std::istringstream str(line);
        n_bodies++;

        size_t id;
        str >> id;

        // " namept1 namept2  "
        //   |---------------| (17 chars)
        str.ignore(1);
        str.read(name, NAME_LEN);
        size_t this_name_len = NAME_LEN - 1;
        while(this_name_len > 0 && name[this_name_len] == ' ') name[this_name_len--] = '\0';
        name_buf.append(name, NAME_STRIDE);

        str.ignore(7); // skip Epoch

        double these_elems[NUM_EL];
        str >> these_elems[AX];
        str >> these_elems[EC];
        str >> these_elems[IN];
        str >> these_elems[PE];
        str >> these_elems[NO];
        str >> these_elems[ME];

        these_elems[IN] *= deg_to_rad;
        these_elems[PE] *= deg_to_rad;
        these_elems[NO] *= deg_to_rad;
        these_elems[ME] *= deg_to_rad;

        elems.insert(elems.end(), std::begin(these_elems), std::end(these_elems));
      }

      Elements<B> asteroid_elems(n_bodies);

      // setup names
      std::copy(name_buf.begin(), name_buf.end(), asteroid_elems._name_buf);
      for (size_t k = 0; k < n_bodies; k++) {
        asteroid_elems._names[k] = std::string_view(&(asteroid_elems._name_buf[k * NAME_STRIDE]), NAME_STRIDE);
      }

      // test conversion
      if (test) {
        double max_rel_err = 0.;
        double max_vel = 0.;
        for (size_t k = 0; k < n_bodies; k++) {
          double *elms = &elems[NUM_EL * k];
          double x, y, z, v_x, v_y, v_z;
          double x_f, y_f, z_f, v_x_f, v_y_f, v_z_f;
          double a, e, i, w, Node, M0;
          double n_steps = 10000;
          double dt = 1. / n_steps;
          double t = dt * n_steps;
          max_vel = std::max(max_vel, std::sqrt(mu_sun * (1. + elms[EC]) / (elms[AX] * (1. - elms[EC]))));
          elements_to_cartesian(elms[AX], elms[EC], elms[IN], elms[PE], elms[NO], elms[ME], 0., x, y, z, v_x, v_y, v_z);
          verlet(dt, n_steps,
                 x, y, z, v_x, v_y, v_z,
                 x_f, y_f, z_f, v_x_f, v_y_f, v_z_f);
          cartesian_to_elements(x_f, y_f, z_f, v_x_f, v_y_f, v_z_f, t, a, e, i, w, Node, M0);
          double diff[NUM_EL] = {
                                elms[AX] - a,
                                elms[EC] - e,
                                elms[IN] - i,
                                elms[PE] - w,
                                elms[NO] - Node,
                                elms[ME] - M0};
          double diff_norm = std::sqrt(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2] + diff[3] * diff[3] + diff[4] * diff[4] + diff[5] * diff[5]);
          double rel_err = diff_norm / std::sqrt(a*a + e*e + i*i + w*w + Node*Node + M0*M0);
          if (rel_err > max_rel_err) {
            std::cout << rel_err << std::endl;
            for (size_t j = 0; j < NUM_EL; j++) std::cout << diff[j] << " ";
            std::cout << std::endl;
          }
          max_rel_err = std::max(max_rel_err, rel_err);
        }
        std::cout << "Maximum relative error of conversion: " << max_rel_err << std::endl;
        std::cout << "Maximum velocity: " << max_vel << std::endl;
      }

      // copy element data in batch form
      for (size_t batch = 0; batch < asteroid_elems._num_batches; batch++) {
        double * src = &asteroid_elems._elements[batch * B * NUM_EL];
        double * dest = &asteroid_elems._elements[batch * B * NUM_EL];
        if ((batch + 1) * B <= n_bodies) {
          for (size_t el = 0; el < NUM_EL; el++) {
            for (size_t i = 0; i < B; i++) {
              dest[el * B + i] = src[i * NUM_EL + el];
            }
          }
        } else {
          size_t lim = n_bodies - batch * B;
          for (size_t el = 0; el < NUM_EL; el++) {
            for (size_t i = 0; i < lim; i++) {
              dest[el * B + i] = src[i * NUM_EL + el];
            }
          }
        }
      }
      return asteroid_elems;
    }

  private:
    size_t _num_batches;
    double * _elements;
    char * _name_buf;
    std::string_view * _names;

    Elements(size_t n_bodies)
      : num_bodies(n_bodies)
    {
      _num_batches = num_bodies / B;
      if (num_bodies % B) _num_batches++;
      _elements = new double[NUM_EL * _num_batches];
      _name_buf = new char[num_bodies * NAME_STRIDE];
      _names = new std::string_view[num_bodies];
    }
};


int main() {
  std::ifstream asteroid_file("elements.txt", std::ifstream::in);
  if (!asteroid_file.is_open()) {
    std::exit(EXIT_FAILURE);
  }
  auto asteroid_elems = Elements<>::from_numbered_asteroid_format(asteroid_file, true);
}
