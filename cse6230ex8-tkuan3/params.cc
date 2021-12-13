
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <getopt.h>
#include "params.h"

const size_t N_d = 10;
const size_t T_d = 10;
const size_t S_d = 10;

void print_options(char * prog_name) {
      std::cerr << "Usage: " << prog_name  << " [-N N] [-T T] [-S S] [-var]" << std::endl;

      std::cerr << "  -N N : number of grid points in each direction ("    << N_d << ")" << std::endl;
      std::cerr << "  -T T : number of outer steps in the Jacobi solver (" << T_d << ")" << std::endl;
      std::cerr << "  -S S : number of inner steps in the Jacobi solver (" << S_d << ")" << std::endl;
      std::cerr << "  -var : choose the variant halo exchange algorithm (false)"         << std::endl;
}

Params::Params(int argc, char **argv):
  N(N_d),
  T(T_d),
  S(S_d),
  var(false)
{
  program_name = std::string(argv[0]);
  static struct option long_options[] =
  {
    /*  0 */ {"N",   required_argument, NULL, 0},
    /*  1 */ {"T",   required_argument, NULL, 0},
    /*  2 */ {"S",   required_argument, NULL, 0},
    /*  3 */ {"var", no_argument,       NULL, 0},
             {NULL,  no_argument,       NULL, 0},
  };

  char **argv_copy = new char*[argc];
  for (int a = 0; a < argc; a++) argv_copy[a] = argv[a];

  while (true) {

    int option_index;
    opterr = 0;
    char c = getopt_long_only(argc, argv_copy, "h", long_options, &option_index);

    if (c == -1) {
      break;
    }

    switch (c) {
    case 0:
      switch (option_index) {
      case 0: N = strtoul(optarg, NULL, 10); break;
      case 1: T = strtoul(optarg, NULL, 10); break;
      case 2: S = strtoul(optarg, NULL, 10); break;
      case 3: var = true; break;
      default: break;
      }
      break;
    case 'h':
      print_options(argv[0]);
      exit(0);
    case '?':
      continue;
    default:
      break;
    }
  }

  delete [] argv_copy;

  _argc = argc;
  _argv = argv;
}

Params::~Params() {
}

void Params::write(Json &log) {
  log.write("program name", program_name);
  log.write("number of grid points per direction", N);
  log.write("number of outer solver steps", T);
  log.write("number of inner solver steps", S);
  std::string b(var ? "true" : "false");
  log.write("variant algorithm", b);
}
