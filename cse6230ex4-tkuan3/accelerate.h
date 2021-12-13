#if !defined(ACCELERATE_H)
#define      ACCELERATE_H

#include "params.h"
#include "state.h"

class Accel {
  public:
    Accel(Params &);
    ~Accel();
    void accelerate(State &, State &);
    void write(Json &);
  private:
    class Impl;
    Impl *_impl;
};

#endif
