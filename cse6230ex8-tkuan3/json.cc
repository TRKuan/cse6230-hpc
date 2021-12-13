
#include <iostream>
#include <mpi.h>
#include "json.h"
#include "util.h"

Json::~Json() {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    *_fs << "\n}" << std::endl;
    if (_fs != &std::cout) {
      delete _fs;
    }
  }
}

Json::Json(std::string &name): _delim("") {
  int rank;
  int err = MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_ECHK(err);
  if (rank == 0) {
    _fs = new std::ofstream(name);
    *_fs << "{";
  } else {
    _fs = NULL;
  }
}

Json::Json(): _delim("") {
  int rank;
  int err = MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_ECHK(err);
  if (rank == 0) {
    _fs = &std::cout;
    *_fs << "{";
  } else {
    _fs = NULL;
  }
}

template void Json::write<double>(const char *, double &);
template void Json::write<size_t>(const char *, size_t &);
template void Json::write<bool>(const char *, bool &);
template<typename T>
  void Json::write(const char *name, T& t) {
    int rank;
    int err = MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_ECHK(err);
    if (rank == 0) {
      *_fs << _delim << std::endl << "  \"" << name << "\": " << t;
      _delim = ",";
    }
  }

template<> void Json::write<std::string>(const char *name, std::string &t) {
    int rank;
    int err = MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_ECHK(err);
    if (rank == 0) {
      *_fs << _delim << std::endl << "  \"" << name << "\": \"" << t << "\"";
      _delim = ",";
    }
}

template<typename T>
  void Json::write_vec(const char *name, std::vector<T> &v) {
    int rank;
    int err = MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_ECHK(err);
    if (rank == 0) {
      *_fs << _delim << std::endl << "  \"" << name << "\": [";
      const char *_vdelim = "";
      for (T &t: v) {
        *_fs << _vdelim << " " << t;
        _vdelim = ",";
      }
      *_fs << "]";
      _delim = ",";
    }
  }

template void Json::write_vec<double>(const char *, std::vector<double> &);
