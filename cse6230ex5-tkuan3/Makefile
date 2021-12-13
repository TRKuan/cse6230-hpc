include make.defs.nvhpc
CXXVERFLAGS=-std=c++17
CXXINCFLAGS=-I/usr/local/pace-apps/manual/packages/nvhpc/20.7/Linux_x86_64/20.7/compilers/include-gcc70 -I/usr/local/pace-apps/manual/packages/nvhpc/20.7/Linux_x86_64/20.7/compilers/include-stdpar
CXX:=$(CXX) -stdpar=gpu
CXXINFOFLAGS=-g -Wall


all: qsort

qsort: qsort.cu
	$(CUDAC) $(CUDAVERFLAGS) $(CUDAOPTFLAGS) $(CUDAINFOFLAGS) -L${NVHPC_ROOT}/Linux_x86_64/20.7/math_libs/11.0/lib64 -lcurand -o $@ $<

clean:
	rm -f *.o qsort

.PHONY: clean
