#define _USE_MATH_DEFINES
#include <cmath>
#include "exact.h"

GridFunction create_rhs(const Params &params) {
  size_t n = params.N;
  size_t np1 = n+1;

  GridFunction F(MPI_COMM_WORLD, n, n, n, params.var);

  auto o = F.local_offsets();
  auto l = F.local_dims();

  for (ssize_t i = -1; i <= (ssize_t) l[0]; i++) {
    for (ssize_t j = -1; j <= (ssize_t) l[1]; j++) {
      for (ssize_t k = -1; k <= (ssize_t) l[2]; k++) {
        F[{i,j,k}] = sin(M_PI * 1. * ((double) (i+o[0]+1) / (double) np1)) *
                     sin(M_PI * 2. * ((double) (j+o[1]+1) / (double) np1)) *
                     sin(M_PI * 3. * ((double) (k+o[2]+1) / (double) np1));
      }
    }
  }
  return F;
}

GridFunction create_sol(const Params &params, const GridFunction& F) {
  size_t n = params.N;
  size_t np1 = n+1;
  auto l = F.local_dims();

  GridFunction Xexact(F);

  double factor = (2.*(1. - cos(M_PI * 1. / (double) np1)) +
                   2.*(1. - cos(M_PI * 2. / (double) np1)) +
                   2.*(1. - cos(M_PI * 3. / (double) np1)));

  for (ssize_t i = -1; i <= (ssize_t) l[0]; i++) {
    for (ssize_t j = -1; j <= (ssize_t) l[1]; j++) {
      for (ssize_t k = -1; k <= (ssize_t) l[2]; k++) {
        Xexact[{i,j,k}] = F[{i,j,k}] / factor;
      }
    }
  }

  return Xexact;
}

double convergence_factor(const Params &params) {
  size_t n = params.N;
  size_t np1 = n+1;
  double factor = (2.*(1. - cos(M_PI * 1. / (double) np1)) +
                   2.*(1. - cos(M_PI * 2. / (double) np1)) +
                   2.*(1. - cos(M_PI * 3. / (double) np1)));

  return log(1. - (1./6.)*factor);
}
