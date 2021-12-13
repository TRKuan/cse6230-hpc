
include make.defs

all: main

%.o: %.cc
	$(CXX) $(OMPFLAGS) $(CXXVERFLAGS) $(CXXFLAGS) $(CXXINFOFLAGS) $(CXXOPTFLAGS) -c $< -o $@

main: main.o state.o timestep.o
	$(CXX) $(OMPFLAGS) $^ -o $@

clean:
	rm -rf *.o main

.PHONY: clean

