#if !defined(TIMESTEP_H)
#define TIMESTEP_H

#include <cstddef>

namespace timestep {

  template<typename Real>
    void verlet_timesteps(size_t n_steps,
                          size_t n_charges,
                          Real dt,
                          Real *__restrict__ x,
                          Real *__restrict__ y,
                          Real *__restrict__ z,
                          Real *__restrict__ u,
                          Real *__restrict__ v,
                          Real *__restrict__ w);
}

#endif
