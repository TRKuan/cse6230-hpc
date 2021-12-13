#if !defined(JSON_H)
#define JSON_H
#include <string>
#include <fstream>
#include <vector>
#include "cloud_util.h"

class Json {
  public:
    Json();
    Json(std::string &name);
    ~Json();
    template<typename T> void write(const char *, T&);
    template<typename T> void write_vec(const char *, std::vector<T> &);
  private:
    std::ostream *_fs;
    const char *_delim;
};

#endif
