#if !defined(CSE6230RAND_H)
#define      CSE6230RAND_H

#include "cloud_util.h"

/* Wrappers to Random123 */

#include <cstdint>
#include <cmath>

#include <Random123/include/Random123/threefry.h>

#if defined(RAND_ARS)
#include <Random123/include/Random123/ars.h>
#define RANDINIT ars4x32keyinit
#define RANDKEY  ars4x32_key_t
#define RANDCTR  ars4x32_ctr_t
#define RANDGEN  ars4x32
#define RANDMAX  UINT32_MAX
#define RAND_t   uint32_t
#else
#include <Random123/include/Random123/threefry.h>
#define RANDINIT threefry4x64keyinit
#define RANDKEY  threefry4x64_key_t
#define RANDCTR  threefry4x64_ctr_t
#define RANDGEN  threefry4x64
#define RANDMAX  UINT64_MAX
#define RAND_t   uint64_t
#endif

namespace randgen {
  class Rand {
    public:
      /// create a pseudo-random number generator with a given seed
      //
      // the PRNG has a mutating state for calls to urand() but can also
      // be used with hash-based random number generation in urand_hash_4()
      // and nrand_hash_4()
      Rand(int seed) {
        RANDKEY uk = {{0}};
        uk.v[0] = seed;
        _k = RANDINIT(uk);
        for (int i = 0; i < 4; i++) {
          _c.v[i] = 0;
        }
        _r = RANDGEN(_c,_k);
        _count = 0;
        _tag = 1; // tag 0 is reserved for rand(), other tags for rand_hash_4()
      }

      /// extract a double from the PRNG
      //
      // this modifies the st
      double urand(double lo=0.0, double hi=1.0) {
        const double scale = (hi - lo) / ((double) RANDMAX + 1.);
        const double shift = lo + scale / 2.;
        double       ret;

        OMP_CRITICAL(urand)
        {
          size_t       mod   = (_count++) % 4;

          ret = _r.v[mod] * scale + shift;

          if (mod == 3) {
            for (int i = 3; i > 0; i--) {
              if (_c.v[i] == RANDMAX) {
                _c.v[i] = 0;
              } else {
                _c.v[i]++;
              }
              break;
            }
            _r = RANDGEN(_c,_k);
          }
        }

        return ret;
      }

      /// get a new unique tag for use with urand_hash_4() and nrand_hash_4()
      size_t get_tag()
      {
        size_t this_tag;
        OMP_ATOMIC(capture)
          this_tag = _tag++;
        return this_tag;
      }

      void urand_hash_4(size_t tag, size_t index1, size_t index2, size_t index3, double out[], double lo=0.0, double hi=1.0)
      {
        const double scale = (hi - lo) / ((double) RANDMAX + 1.);
        const double shift = lo + scale / 2.;
        RANDCTR cc;
        RANDCTR rr;

        cc.v[0] = tag;
        cc.v[1] = index1;
        cc.v[2] = index2;
        cc.v[3] = index3;

        rr = RANDGEN(cc, _k);
        OMP_SIMD()
        for (int i = 0; i < 4; i++) {
          out[i] = rr.v[i] * scale + shift;
        }
      }

      void nrand_hash_4(size_t tag, size_t index1, size_t index2, size_t index3, double out[], double mu=0.0, double sigma=1.0)
      {
        const double uscale = 1. / ((double) RANDMAX + 1.);
        const double ushift = uscale / 2.;
        RANDCTR cc;
        RANDCTR rr;

        cc.v[0] = tag;
        cc.v[1] = index1;
        cc.v[2] = index2;
        cc.v[3] = index3;

        rr = RANDGEN(cc, _k);
        for (int i = 0; i < 4; i += 2) {
          double u1 = rr.v[i] * uscale + ushift;
          double u2 = rr.v[i+1] * uscale + ushift;
          double z1 = sqrt(-2. * log(u1)) * cos(CSE6230_PI * 2. * u2);
          double z2 = sqrt(-2. * log(u1)) * sin(CSE6230_PI * 2. * u2);
          out[i] = z1 * sigma + mu;
          out[i+1] = z2 * sigma + mu;
        }
      }

    private:
      RANDCTR _c;
      RANDKEY _k;
      RANDCTR _r;
      size_t _count;
      size_t _tag;
  };
}

typedef struct _cse6230rand
{
  threefry4x64_ctr_t c;
  threefry4x64_key_t k;
  threefry4x64_ctr_t r;
  size_t count;
  size_t tag;
}
cse6230rand_t;

static inline void cse6230rand_seed(int seed, cse6230rand_t *__restrict__ rand)
{
  threefry4x64_key_t uk = {{0}};

  uk.v[0] = seed;
  rand->k = threefry4x64keyinit(uk);
  rand->c.v[0] = 0;
  rand->c.v[1] = 1;
  rand->c.v[2] = 2;
  rand->c.v[3] = 3;
  rand->r = threefry4x64(rand->c,rand->k);
  rand->count = 0;
  rand->tag = 0;
}

static inline double cse6230rand(cse6230rand_t *__restrict__ rand)
{
  const double scale = 1. / ((double) UINT64_MAX + 1.);
  const double shift = scale / 2.;
  size_t       mod   = (rand->count++) % 4;
  double       ret;

  ret = rand->r.v[mod] * scale + shift;

  if (mod == 3) {
    for (int i = 0; i < 4; i++) {rand->c.v[i] += 4;}
    rand->r = threefry4x64(rand->c,rand->k);
  }

  return ret;
}

static inline size_t cse6230rand_get_tag(cse6230rand_t *__restrict__ rand)
{
  return rand->tag++;
}

/* gives four random numbers for four indices */
static inline void cse6230rand_hash(cse6230rand_t *__restrict__ rand,
                                    size_t tag,
                                    size_t index1,
                                    size_t index2,
                                    size_t index3,
                                    double rand_out[])
{
  const double scale = 1. / ((double) UINT64_MAX + 1.);
  const double shift = scale / 2.;
  threefry4x64_ctr_t c;
  threefry4x64_ctr_t r;
  int i;

  c.v[0] = tag;
  c.v[1] = index1;
  c.v[2] = index2;
  c.v[3] = index3;

  r = threefry4x64(c, rand->k);
  for (i = 0; i < 4; i++) {
    rand_out[i] = r.v[i] * scale + shift;
  }
}

/* gives three normally distributed random numbers from four indices */
static inline void cse6230rand_normal_hash(cse6230rand_t *__restrict__ rand,
                                           size_t tag,
                                           size_t index1,
                                           size_t index2,
                                           size_t index3,
                                           double rand_out[])
{
  const double scale = 1. / (UINT64_MAX + 1.);
  const double shift = scale / 2.;
  threefry4x64_ctr_t c;
  threefry4x64_ctr_t r;

  c.v[0] = tag;
  c.v[1] = index1;
  c.v[2] = index2;
  c.v[3] = index3;

  r = threefry4x64(c, rand->k);
  for (int i = 0; i < 4; i += 2) {
    double u1 = r.v[i] * scale + shift;
    double u2 = r.v[i+1] * scale + shift;
    double z1 = sqrt(-2. * log(u1)) * cos(CSE6230_PI * 2. * u2);
    double z2 = sqrt(-2. * log(u1)) * sin(CSE6230_PI * 2. * u2);
    rand_out[i] = z1;
    rand_out[i+1] = z2;
  }
}

typedef struct _cse6230nrand
{
  cse6230rand_t urand;
  double        z[4];
  size_t        count;
}
cse6230nrand_t;

static inline void _cse6230nrand_batch(cse6230nrand_t *__restrict__ nrand)
{
  const double scale = 1. / (UINT64_MAX + 1.);
  const double shift = scale / 2.;

  for (int i = 0; i < 4; i++) {nrand->urand.c.v[i] += 4;}
  nrand->urand.r = threefry4x64(nrand->urand.c,nrand->urand.k);
  for (int i = 0; i < 4; i++) {nrand->z[i] = nrand->urand.r.v[i] * scale + shift;}
  nrand->urand.count += 4;

  for (int i = 0; i < 4; i += 2) {
    double u1 = nrand->z[i];
    double u2 = nrand->z[i+1];
    double z1 = sqrt(-2. * log(u1)) * cos(CSE6230_PI * 2. * u2);
    double z2 = sqrt(-2. * log(u1)) * sin(CSE6230_PI * 2. * u2);
    nrand->z[i]   = z1;
    nrand->z[i+1] = z2;
  }
}

static inline void cse6230nrand_seed(int seed, cse6230nrand_t *__restrict__ nrand)
{
  cse6230rand_seed(seed,&(nrand->urand));
  _cse6230nrand_batch(nrand);
  nrand->count = 0;
}

static inline double cse6230nrand(cse6230nrand_t *__restrict__ nrand)
{
  size_t       mod   = (nrand->count++) % 4;
  double       ret;

  ret = nrand->z[mod];

  if (mod == 3) {
    _cse6230nrand_batch(nrand);
  }

  return ret;
}

#endif
