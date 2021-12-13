
#include <iostream>
#include "json.h"

Json::~Json() {
  *_fs << "\n}" << std::endl;
  if (_fs != &std::cout) {
    delete _fs;
  }
}

Json::Json(std::string &name): _delim("") {
  _fs = new std::ofstream(name);
  *_fs << "{";
}

Json::Json(): _delim("") {
  _fs = &std::cout;
  *_fs << "{";
}

template void Json::write<double>(const char *, double &);
template void Json::write<size_t>(const char *, size_t &);
template void Json::write<bool>(const char *, bool &);
template<typename T>
  void Json::write(const char *name, T& t) {
    *_fs << _delim << std::endl << "  \"" << name << "\": " << t;
    _delim = ",";
  }

template<> void Json::write<std::string>(const char *name, std::string &t) {
    *_fs << _delim << std::endl << "  \"" << name << "\": \"" << t << "\"";
    _delim = ",";
}

template<typename T>
  void Json::write_vec(const char *name, std::vector<T> &v) {
    *_fs << _delim << std::endl << "  \"" << name << "\": [";
    const char *_vdelim = "";
    for (T &t: v) {
      *_fs << _vdelim << " " << t;
      _vdelim = ",";
    }
    *_fs << "]";
    _delim = ",";
  }

template void Json::write_vec<double>(const char *, std::vector<double> &);
