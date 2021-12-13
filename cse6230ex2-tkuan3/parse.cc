#include <fstream>
#include <sstream>
#include "parse.h"

namespace parse {

  template<typename Real>
    std::vector<Real> asteroid_file(std::string &filename)
    {
      using namespace elements;
      const Real deg_to_rad = 3.141592653589793 / 180.;

      std::ifstream ist(filename, std::ifstream::in);
      std::string name_buf;
      std::string line;
      std::vector<Real> elems;

      std::getline(ist, line); // skip field names
      std::getline(ist, line); // skip ----- -----
      while (std::getline(ist, line)) {
        std::istringstream str(line);
        str.ignore(31); // skip name and epoch

        Real these_elems[NUM_EL];
        str >> these_elems[AX];
        str >> these_elems[EC];
        str >> these_elems[IN];
        str >> these_elems[PE];
        str >> these_elems[NO];
        str >> these_elems[ME];

        these_elems[IN] *= deg_to_rad;
        these_elems[PE] *= deg_to_rad;
        these_elems[NO] *= deg_to_rad;
        these_elems[ME] *= deg_to_rad;

        for (size_t i = 0; i < NUM_EL; i++) {
          elems.push_back(these_elems[i]);
        }
      }

      return elems;
    }

    template std::vector<double> asteroid_file<double>(std::string &);
    template std::vector<float> asteroid_file<float>(std::string &);
}
