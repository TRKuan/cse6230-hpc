#if !defined(PARSE_H)
#define PARSE_H

#include <string>
#include <vector>
#include "elements.h"

namespace parse {
  template<typename Real> std::vector<Real> asteroid_file(std::string &filename);
}

#endif // PARSE_H
