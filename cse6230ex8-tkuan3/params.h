#if !defined(PARAMS_H)
#define PARAMS_H

#include <cstddef>
#include <string>
#include "json.h"

class Params {
  public:
    size_t N;
    size_t T;
    size_t S;
    bool   var;
    std::string program_name;
    int    _argc;
    char **_argv;

    Params(int argc, char **argv);
    ~Params();
    void write(Json &);
};

#endif
