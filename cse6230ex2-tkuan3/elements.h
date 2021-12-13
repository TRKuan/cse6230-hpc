#if !defined(ELEMENTS_H)
#define ELEMENTS_H

namespace elements
{
  /// Orbital elements
  enum {
    AX, /// [AU ] semi-major axis
    EC, /// [1  ] eccentricity
    IN, /// [rad] inclination
    PE, /// [rad] argument of periapsis
    NO, /// [rad] longitude of ascending node
    ME, /// [rad] mean anomaly
    NUM_EL
  };

  template<typename Real>
    static inline void pack(Real a, Real e, Real i, Real w, Real Node, Real M, Real elements[])
    {
      elements[AX] = a;
      elements[EC] = e;
      elements[IN] = i;
      elements[PE] = w;
      elements[NO] = Node;
      elements[ME] = M;
    }

  template<typename Real>
    static inline void unpack(const Real elements[], Real &a, Real &e, Real &i, Real &w, Real &Node, Real &M)
    {
      a = elements[AX];
      e = elements[EC];
      i = elements[IN];
      w = elements[PE];
      Node = elements[NO];
      M = elements[ME];
    }

  template<typename Real> void elements_to_state_vector(Real mu, Real t, const Real elems[], Real state[]);

  template<typename Real> void state_vector_to_elements(Real mu_sun, Real t, const Real state[], Real elems[]);
}

#endif /* define ELEMENTS_H */
