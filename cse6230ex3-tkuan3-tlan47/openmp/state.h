#if !defined(STATE_H)
#define      STATE_H

#include <cstddef>

namespace state {

  template<typename Real>
  class State {

    public:

      State(size_t n_charges_per_dir);

      size_t n_charges();

      void arrays(Real *x[], Real *y[], Real *z[], Real *u[], Real *v[], Real *w[]);

      ~State();

    private:

      size_t _n_charges;
      Real *_x;
      Real *_y;
      Real *_z;
      Real *_u;
      Real *_v;
      Real *_w;

  };

  template class State<float>;
  template class State<double>;
}
#endif
