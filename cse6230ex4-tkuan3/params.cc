
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <getopt.h>
#include "params.h"
#include "json.h"
#include "rand.h"

const size_t Np_d      = 256;
const size_t Nt_d      = 1000;
const double dt_d      = 1.e-4;
const double k_d       = 100.;
const double d_d       = 1.;
const double L_d       = 20.;
const double r_d       = 1.;
const bool   direct_d  = false;
const size_t Nc_d      = SIZE_MAX;
const RAND_t seed_d    = 6230;


void print_options(char * prog_name) {
      std::cerr << "Usage: " << prog_name  << " [-Np NP] [-Nt NT] [-dt DT] [-k K] [-d D] [-L L] [-r R] [-Nc NC] [-seed S]" << std::endl;

      std::cerr << "  -Np NP                     : number of particles      (default " << Np_d      << ")" << std::endl;
      std::cerr << "  -Nt NT                     : number of timesteps      (default " << Nt_d      << ")" << std::endl;
      std::cerr << "  -dt DT                     : timestep length          (default " << dt_d      << ")" << std::endl;
      std::cerr << "  -k  K                      : interaction strength     (default " << k_d       << ")" << std::endl;
      std::cerr << "  -d  D                      : diffusive strength       (default " << d_d       << ")" << std::endl;
      std::cerr << "  -L  L                      : periodic box length      (default " << L_d       << ")" << std::endl;
      std::cerr << "  -r  R                      : particle radius          (default " << r_d       << ")" << std::endl;
      std::cerr << "  -Nc NC                     : timesteps per checkpoint (default " << Nc_d      << ")" << std::endl;
      std::cerr << "  -seed                      : PRNG seed                (default " << seed_d    << ")" << std::endl;
}

Params::Params(int argc, char **argv):
  Np(Np_d),
  Nt(Nt_d),
  Nc(Nc_d),
  dt(dt_d),
  k(k_d),
  d(d_d),
  L(L_d),
  r(r_d),
  seed(seed_d)
{
  program_name = std::string(argv[0]);
  static struct option long_options[] =
  {
    /*  0 */ {"Np",             required_argument, NULL, 0},
    /*  1 */ {"Nt",             required_argument, NULL, 0},
    /*  2 */ {"Nc",             required_argument, NULL, 0},
    /*  3 */ {"dt",             required_argument, NULL, 0},
    /*  4 */ {"d",              required_argument, NULL, 0},
    /*  5 */ {"k",              required_argument, NULL, 0},
    /*  6 */ {"L",              required_argument, NULL, 0},
    /*  7 */ {"r",              required_argument, NULL, 0},
    /*  8 */ {"seed",           required_argument, NULL, 0},
             {NULL,             no_argument,       NULL, 0},
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
      case 0: Np = strtoul(optarg, NULL, 10); break;
      case 1: Nt = strtoul(optarg, NULL, 10); break;
      case 2: Nc = strtoul(optarg, NULL, 10); break;
      case 3: dt = atof(optarg); break;
      case 4: d  = atof(optarg); break;
      case 5: k  = atof(optarg); break;
      case 6: L  = atof(optarg); break;
      case 7: r  = atof(optarg); break;
      case 8: seed = (RAND_t) strtoul(optarg, NULL, 10); break;
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

  Nc = MIN(Nc, Nt);
  _argc = argc;
  _argv = argv;
}

Params::~Params() {
}

void Params::write(Json &log) {
  log.write("program name", program_name);
  log.write("number of particles", Np);
  log.write("number of time steps", Nt);
  log.write("time steps per checkpoint", Nc);
  log.write("time step dt", dt);
  log.write("interaction constant k", k);
  log.write("diffusive constant d", d);
  log.write("domain width L", L);
  log.write("particle radius r", r);
  log.write("seed", seed);
}
