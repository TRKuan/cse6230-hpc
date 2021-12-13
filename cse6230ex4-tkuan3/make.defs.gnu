CC=gcc-10
CVERFLAGS=-std=c11
CINFOFLAGS=-g3 -Wall -fopt-info-vec -fopt-info-omp
COPTFLAGS=-O3 -march=native -ffast-math -fopenmp-simd
OMPFLAGS=-fopenmp
CXX=g++-10
CXXVERFLAGS=-std=c++17
CXXINFOFLAGS=$(CINFOFLAGS)
CXXOPTFLAGS=$(COPTFLAGS)
