#if !defined(STATE_H)
#define      STATE_H

#include <cstddef>
#include <memory>
#include <vector>
#include "rand.h"

/// The positions of the particles
class State {
  private:
    size_t _Np;
    double *_x;
    double *_y;
    double *_z;

  public:
    State(State &);

    State(size_t);

    State(size_t, randgen::Rand &, double lo=-1.0, double hi=1.0);

    void get_arrays(size_t &Np, double *&x, double *&y, double *&z);

    size_t Np() {return _Np;}

    void copy(State &);

    double mean_square_distance(State &);

    ~State();
};
#endif
