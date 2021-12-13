#if !defined(PARAMS_H)
#define PARAMS_H

#include <cstddef>
#include <string>
#include "rand.h"
#include "json.h"

class Params {
  public:
    size_t Np;
    size_t Nt;
    size_t Nc;
    double dt;
    double k;
    double d;
    double L;
    double r;
    RAND_t seed;
    std::string program_name;
    int    _argc;
    char **_argv;

    Params(int argc, char **argv);
    ~Params();
    void write(Json &);
};

#endif
