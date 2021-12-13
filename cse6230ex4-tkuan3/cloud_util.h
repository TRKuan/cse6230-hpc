#if !defined(CLOUD_UTIL_H)
#define      CLOUD_UTIL_H

#include <cstddef>
#include <cstdlib>
#include <cstdio>
#include <cstdint>
#include <cmath>

#if defined(M_PI)
#define CSE6230_PI M_PI
#else
#define CSE6230_PI 3.1415926535897932384626433832795L
#endif

#define PRAGMA(x) _Pragma(#x)

#if defined(_OPENMP)
#define OMP_ATOMIC(x) PRAGMA(omp atomic x)
#define OMP_CRITICAL(name) PRAGMA(omp critical(name))
#define OMP_PARALLEL_FOR(x) PRAGMA(omp parallel for x)
#define OMP_SIMD(x) PRAGMA(omp simd x)
#else
#define OMP_ATOMIC(x)
#define OMP_CRITICAL(name)
#define OMP_PARALLEL_FOR(x)
#define OMP_SIMD(x) PRAGMA(omp simd x)
#endif

#define MIN(a,b) ((a) <= (b)) ? (a) : (b)
#define MAX(a,b) ((a) > (b)) ? (a) : (b)

#pragma omp declare simd
static inline double _remainder(double x, double y)
{
  double n = round(x / y);
  return x - n * y;
}

/* get the square distance and displacement between two particles under periodic
 * conditions */
#pragma omp declare simd
static inline double
dist_and_disp (double x1, double y1, double z1, /* The center of the first particle */
               double x2, double y2, double z2, /* The center of the second particle */
               double L, /* The width of the periodic domain */
               double *Dx, double *Dy, double *Dz)
{
  double dx, dy, dz;
  *Dx = dx = _remainder(x1 - x2, L);
  *Dy = dy = _remainder(y1 - y2, L);
  *Dz = dz = _remainder(z1 - z2, L);

  return dx*dx + dy*dy + dz*dz;
}

#endif
