#if !defined(TIMESTEP_H)
#define      TIMESTEP_H

#include "rand.h"
#include "accelerate.h"
#include "state.h"
#include "params.h"
#include "json.h"

class Timestepper {
  public:
    Timestepper(Params &, randgen::Rand &);

    void step(size_t, double, State &);
    void write(Json &);

  private:

    void _stream_and_noise(double, double, State &, State &);
    Params *_params;
    randgen::Rand *_gen;
    Accel _accel;
    double _total_particle_steps;
    double _total_time;
};

#endif
