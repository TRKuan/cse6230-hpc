#if !defined(EXACT_H)
#define EXACT_H

#include "params.h"
#include "grid.h"

GridFunction create_rhs(const Params &params);
GridFunction create_sol(const Params &params, const GridFunction& F);
double convergence_factor(const Params &params);

#endif
