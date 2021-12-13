
#include <cmath>
#include <vector>
#include <functional>
#include <tictoc.h>
#include <getopt.h>
#include "cloud_util.h"
#include "accelerate.h"

/// Particle interaction kernel: compute the effect of the in particles on the
//  out particles by adding their acceleration to the velocities of the out
//  particles
//
//  Do not put a parallel for in this function: use this kernel to update
//  independent output particles in each thread
//
//  L: size of the periodic domain
//  k: strength of the particle field
//  r: the particle radius
void accelerate_batch (double L, double k, double r,
                       size_t n_in,  const double *x_in, const double *y_in, const double *z_in,
                       size_t n_out, const double *x_out, const double *y_out, const double *z_out,
                       double *u_out, double *v_out, double *w_out)
{
  // this is our n-body O(n^2) loop, like from the last exercise
  for (size_t i = 0; i < n_out; i++) {
    double u = 0;
    double v = 0;
    double w = 0;
    double x1 = x_out[i];
    double y1 = y_out[i];
    double z1 = z_out[i];

    for (size_t j = 0; j < n_in; j++) {

      double du, dv, dw;
      {
        double x2 = x_in[j];
        double y2 = y_in[j];
        double z2 = z_in[j];
        double dx, dy, dz;
        double R2;

        // get the _periodic_ distance and displacement between particles
        R2 = dist_and_disp (x1, y1, z1, x2, y2, z2, L, &dx, &dy, &dz);
        double same_particle = (R2 == 0.);
        double R = sqrt(R2);
        {
          // The interaction strength starts at 0 when they are just touching,
          // becoming infinite as the distance becomes zero
          //
          // If the particles are not close together, most interactions will
          // have zero effect!
          double denom = R + same_particle;
          double strength = (2. * r - R) / denom;
          strength = strength < 0. ? 0. : strength;

          du = k * strength * dx;
          dv = k * strength * dy;
          dw = k * strength * dz;
        }
      }
      u += du;
      v += dv;
      w += dw;
    }
    u_out[i] += u;
    v_out[i] += v;
    w_out[i] += w;
  }
}

/// Implementation details for reducing the number of computed interactions
//  by binning particles: placing them into a spatial grid of boxes that are
//  larger than twice the radius of the particles, so that all non-zero
//  interactions are between neighboring boxes
class Accel::Impl {
  private:

    size_t _boxes_per_dimension; // the number of boxes in each direction: cube this to get the number of boxes
    double _box_width;           // width of each box
    double _min;                 // most negative coordinate of the domain
    double _L;                   // periodic domain width
    double _k;                   // strength of the particle field coefficient
    double _r;                   // radius of the particles
    std::vector<double> _xb; // work vectors for holding the sorting the particles into their boxes
    std::vector<double> _yb;
    std::vector<double> _zb;
    std::vector<double> _ub;
    std::vector<double> _vb;
    std::vector<double> _wb;
    State *_X = nullptr;
    State *_X_last = nullptr;
    double _target_ppb; // a desired number of particles per box: more
                        // particles per box means more zero interactions will be computed, but
                        // the computation will also become more regular, more
                        // vectorizable, and their will be fewer boxes to
                        // manage, so there is a tradeoff here
                        //
                        // control this value with -ppb X on the command line
    double _accel_time;
    double _accel_interactions;

    bool _use_direct; // if true, use the full O(n^2) n-body calculation
                      // instead of binning: useful for debugging and judging
                      // performance
                      //
                      // conrol with with the -use_direct flag on the command
                      // line
    std::vector<size_t> _box_offsets;
    std::vector<size_t> _particle_box;
    std::vector<size_t> _particle_perm;

    void process_options(int argc, char **argv) {
      static struct option long_options[] =
      {
        /* 0 */ {"ppb",        required_argument, NULL, 0},
        /* 1 */ {"use_direct", no_argument,       NULL, 0},
                {NULL,         no_argument,       NULL, 0},
      };

      optind = 1;

      char **argv_copy = new char*[argc];
      for (int a = 0; a < argc; a++) argv_copy[a] = argv[a];

      while (true) {
        int option_index;
        char c = getopt_long_only(argc, argv_copy, "", long_options, &option_index);

        if (c == -1) {
          break;
        }

        if (c == 0) {
          switch (option_index) {
          case 0: _target_ppb = atof(optarg); break;
          case 1: _use_direct = true; break;
          default:
            break;
          }
        }
      }

      delete [] argv_copy;
    }

    /// Take particles in their original order and
    // reorder consecutively into boxes
    //
    // x, y, z: input coordinates of particles
    //
    // particle_box: work array for determining which box each particle is
    //               contained in
    //
    // particle_perm: output array, the location of each particle in box
    //                ordering (have to pass this to unbin_particles)
    // box_offsets: output array (boxes_per_dimension)**3 + 1
    //
    //   box_offsets[b] is the start of box b
    //   box_offsets[b+1] is the start of the next box
    //
    // xb, yb, zb: the particle coordinates after they have been sorted
    // into boxes
    void bin_particles(size_t Np,
                       const double *x,
                       const double *y,
                       const double *z,
                       size_t *particle_box,
                       size_t *particle_perm,
                       size_t *box_offsets,
                       double *xb,
                       double *yb,
                       double *zb)
    {
      size_t bpd = _boxes_per_dimension;
      size_t Nb = bpd * bpd * bpd;

      // start with offsets holding the counts
      #pragma omp for schedule(static)
      for (size_t b = 0; b < Nb + 1; b++) {
        box_offsets[b] = 0;
      }

      #pragma omp for schedule(static)
      for (size_t p = 0; p < Np; p++) {
        // periodic logic for figuring out the x,y,z coordinates of the box
        size_t x_box = (size_t) floor((_remainder(x[p], _L) - _min) / _box_width);
        x_box = x_box >= bpd ? bpd - 1 : x_box;
        size_t y_box = (size_t) floor((_remainder(y[p], _L) - _min) / _box_width);
        y_box = y_box >= bpd ? bpd - 1 : y_box;
        size_t z_box = (size_t) floor((_remainder(z[p], _L) - _min) / _box_width);
        z_box = z_box >= bpd ? bpd - 1 : z_box;

        // box coordinates to index
        size_t box = (x_box * bpd + y_box) * bpd + z_box;

        // note which box it is in
        particle_box[p] = box;
        // increment the counter for the box, which is the location of this
        // particle within this box
        #pragma omp atomic capture
        particle_perm[p] = box_offsets[box+1]++;
      }

      #pragma omp single
      {
        // Do a prefix sum (scan) on the counts to get the final box_offsets
        for (size_t b = 0; b < Nb; b++) {
          box_offsets[b + 1] += box_offsets[b];
        }
      }

      #pragma omp for schedule(static)
      for (size_t p = 0; p < Np; p++) {
        // adding the offsets to the position within the box gives
        // each particle's position with the reordered vector
        particle_perm[p] += box_offsets[particle_box[p]];
        size_t pp = particle_perm[p];
        xb[pp] = x[p];
        yb[pp] = y[p];
        zb[pp] = z[p];
      }
    }

    bool need_rebin(const double *x,
                    const double *y,
                    const double *z)
    {
      if(_X_last == nullptr) {
        _X_last = new State(*_X);
        return true;
      }

      double *x_last, *y_last, *z_last;
      size_t Np_last;
      _X_last->get_arrays(Np_last, x_last, y_last, z_last);
      
      static double max_dist = 0;
      #pragma omp for reduction(max:max_dist)
      for (size_t p = 0; p < Np_last; p++) {
        double dx, dy, dz;
        double dist = sqrt(dist_and_disp(x[p], y[p], z[p], x_last[p], y_last[p], z_last[p], _L, &dx, &dy, &dz));
        max_dist = std::max(dist, max_dist);
      }
      
      if(_box_width - 2 * max_dist > 2 * _r) {
        return false;
      }
      else {
        _X_last->copy(*_X);
        return true;
      }
    }

    // Using the particle permutation computed in bin_particles(), move the
    // velocities back to their original order
    void unbin_particles(size_t Np,
                         const double *ub,
                         const double *vb,
                         const double *wb,
                         size_t *particle_perm,
                         double *u,
                         double *v,
                         double *w)
    {
      #pragma omp for schedule(static)
      for (size_t p = 0; p < Np; p++) {
        size_t pp = particle_perm[p];
        u[p] = ub[pp];
        v[p] = vb[pp];
        w[p] = wb[pp];
      }
    }

    /// Permute the particles into boxes,
    //  compute the velocities from the interactions in neighboring
    //  boxes, and then permute the velocities back to the original order
    void accelerate_indirect(double L, double K, double r, size_t Np,
                             const double *x,
                             const double *y,
                             const double *z,
                             double *u,
                             double *v,
                             double *w)
    {
      // make sure the work variables are long enough

      #pragma omp single
      {
        _xb.resize(Np);
        _yb.resize(Np);
        _zb.resize(Np);

        _ub.resize(Np);
        _vb.resize(Np);
        _wb.resize(Np);

        _particle_box.resize(Np);
        _particle_perm.resize(Np);
      }

      size_t bpd = _boxes_per_dimension;

      double *xb = &_xb[0];
      double *yb = &_yb[0];
      double *zb = &_zb[0];

      double *ub = &_ub[0];
      double *vb = &_vb[0];
      double *wb = &_wb[0];

      // we are going to accumulate in the boxed velocities in ub, vb, wb,
      // so zero them to start
      #pragma omp for schedule(static) nowait
      for (size_t p = 0; p < Np; p++) {
        ub[p] = 0.;
        vb[p] = 0.;
        wb[p] = 0.;
      }

      if(need_rebin(x, y, z)) {
        // reorder the particles
        bin_particles(Np, x, y, z, &_particle_box[0], &_particle_perm[0], &_box_offsets[0], xb, yb, zb);
      }

      // triple loop over boxes: by design, each instance of the main loop
      // body can be done independently
      #pragma omp for collapse(3) schedule(static)
      for (size_t i = 0; i < bpd; i++) {
        for (size_t j = 0; j < bpd; j++) {
          for (size_t k = 0; k < bpd; k++) {
            // start main loop body

            // coordinates to index
            size_t id = (i * bpd + j) * bpd + k;
            // index to the subset of particles in this box
            size_t n_t = _box_offsets[id + 1] - _box_offsets[id];
            const double *xt = &xb[_box_offsets[id]];
            const double *yt = &yb[_box_offsets[id]];
            const double *zt = &zb[_box_offsets[id]];
            double *ut = &ub[_box_offsets[id]];
            double *vt = &vb[_box_offsets[id]];
            double *wt = &wb[_box_offsets[id]];

            // optimization: determine if all the z neighbors are consecutive
            // (not at the periodic boundary)
            size_t nkm = (k + bpd - 1) % bpd;
            size_t nkp = (k + bpd + 1) % bpd;
            size_t nkmin = MIN(nkm, nkp);
            nkmin = MIN(nkmin, k);
            size_t nkmax = MAX(nkm, nkp);
            nkmax = MAX(nkmax, k);

            // for all neighboring boxes
            for (int di = -1; di <= 1; di++) {
              size_t ni = (i + bpd + di) % bpd;
              for (int dj = -1; dj <= 1; dj++) {
                size_t nj = (j + bpd + dj) % bpd;


                // if the z neighbors are consecutive work on three at a
                // time
                if (nkmin + 2 == nkmax) {
                  size_t nidmin = (ni * bpd + nj) * bpd + nkmin;
                  size_t nidmax = (ni * bpd + nj) * bpd + nkmax;
                  size_t n_in = _box_offsets[nidmax + 1] - _box_offsets[nidmin];

                  const double *xn = &xb[_box_offsets[nidmin]];
                  const double *yn = &yb[_box_offsets[nidmin]];
                  const double *zn = &zb[_box_offsets[nidmin]];

                  accelerate_batch(L, K, r, n_in, xn, yn, zn, n_t, xt, yt, zt, ut, vt, wt);
                  _accel_interactions += n_in * n_t;
                } else {
                  // otherwise grab each neighboring box separately
                  for (int dk = -1; dk <= 1; dk++) {
                    size_t nk = (k + bpd + dk) % bpd;

                    // coordinates to index
                    size_t nid = (ni * bpd + nj) * bpd + nk;

                    // get the particles in this neighboring box
                    size_t n_in = _box_offsets[nid + 1] - _box_offsets[nid];
                    const double *xn = &xb[_box_offsets[nid]];
                    const double *yn = &yb[_box_offsets[nid]];
                    const double *zn = &zb[_box_offsets[nid]];

                    // compute the effect of the neighboring (ni, nj, nk) particles on the
                    // (i,j,k) particles
                    accelerate_batch(L, K, r, n_in, xn, yn, zn, n_t, xt, yt, zt, ut, vt, wt);
                    _accel_interactions += n_in * n_t;
                  }
                }
              }
            }
            // end main loop body
          }
        }
      }
      // copy the permuted velocities back to the original order
      unbin_particles(Np, ub, vb, wb, &_particle_perm[0], u, v, w);
    }


  public:

    /// Setup the necessary data for binning the periodic cube of length L
    //  into boxes that have a width at least l
    Impl(size_t Np, double L, double l, double k, double r, int argc, char **argv):
      _k(k),
      _r(r),
      _target_ppb(16),
      _accel_time(0.),
      _accel_interactions(0.),
      _use_direct(false)
    {

      _L = L;
      _min = -L / 2;

      process_options(argc, argv);

      // compute the boxes per dimension (bpd) and box width:
      // - the minimum width places a maximum on bpd
      size_t max_bpd = MAX((size_t) floor(L / l), (size_t) 1);
      // - the target particles per box suggests a target boxes per dimension
      //   floor => err on the side of too many particles per box
      size_t target_bpd = (size_t) floor(pow((double) Np / _target_ppb, 1./3.));
      _boxes_per_dimension = MIN(max_bpd, target_bpd);
      _box_width = L / _boxes_per_dimension;
      size_t bpd = _boxes_per_dimension;

      // construct work vectors
      _xb = std::vector<double>();
      _yb = std::vector<double>();
      _zb = std::vector<double>();
      _ub = std::vector<double>();
      _vb = std::vector<double>();
      _wb = std::vector<double>();
      _box_offsets = std::vector<size_t>(bpd*bpd*bpd+1);
      _particle_box = std::vector<size_t>();
      _particle_perm = std::vector<size_t>();
    }

    void write(Json &log)
    {
      log.write("use direct interactions", _use_direct);
      log.write("target particles per box", _target_ppb);
      log.write("boxes per dimension", _boxes_per_dimension);
      log.write("total acceleration time", _accel_time);
      double rate = _accel_interactions / _accel_time;
      log.write("particle interactions per second", rate);
    }

    void accelerate(State &X, State &U) {
      TicTocTimer timer = tic();

      double *x, *y, *z;
      double *u, *v, *w;
      size_t Np;
      X.get_arrays(Np, x, y, z);
      U.get_arrays(Np, u, v, w);
      if (_use_direct || _boxes_per_dimension < 4) {
        for (size_t i = 0; i < Np; i++) {
          u[i] = 0.;
          v[i] = 0.;
          w[i] = 0.;
        }
        accelerate_batch(_L, _k, _r, Np, x, y, z, Np, x, y, z, u, v, w);
        _accel_interactions += Np * Np;
      } else {
        _X = &X;
        #pragma omp parallel
        {
          accelerate_indirect(_L, _k, _r, Np, x, y, z, u, v, w);
        }
      }
      _accel_time += toc(&timer);
    }

};

Accel::~Accel() {
  delete _impl;
}

Accel::Accel(Params &params) {
  _impl = new Impl(params.Np, params.L, 2*params.r, params.k, params.r, params._argc, params._argv);
}

void Accel::accelerate(State &X, State &U) {
  _impl->accelerate(X, U);
}

void Accel::write(Json &log) {
  _impl->write(log);
}


