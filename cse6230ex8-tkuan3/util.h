#if !defined(UTIL_H)
#define UTIL_H
#include <iostream>
#include <vector>
#include <exception>

// throws an exception (so this can be called in constructors)
#define MPI_ECHK(err)         \
  do {                        \
    if (err != MPI_SUCCESS) { \
      throw err;              \
    }                         \
  } while(0)

#endif
