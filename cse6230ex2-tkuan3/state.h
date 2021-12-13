#if !defined(STATE_H)
#define      STATE_H

namespace state {

  // State vector components
  enum {
    X, // position components
    Y,
    Z,
    U, // velocity components
    V,
    W
  };

  template<typename Real>
    static inline void pack(Real x, Real y, Real z, Real u, Real v, Real w, Real state[])
    {
      state[X] = x;
      state[Y] = y;
      state[Z] = z;
      state[U] = u;
      state[V] = v;
      state[W] = w;
    }

  template<typename Real>
    static inline void unpack(const Real state[], Real &x, Real &y, Real &z, Real &u, Real &v, Real &w)
    {
      x = state[X];
      y = state[Y];
      z = state[Z];
      u = state[U];
      v = state[V];
      w = state[W];
    }
}
#endif
