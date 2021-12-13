
include make.defs

all: simulation

timestep.o: timestep.cc
	$(CXX) $(OMPFLAGS) $(CXXVERFLAGS) $(CXXFLAGS) $(CXXINFOFLAGS) $(CXXOPTFLAGS) $(TIMESTEPFLAGS) -c $< -o $@

%.o: %.cc
	$(CXX) $(OMPFLAGS) $(CXXVERFLAGS) $(CXXFLAGS) $(CXXINFOFLAGS) $(CXXOPTFLAGS) -c $< -o $@

simulation: simulation.o parse.o elements.o timestep.o
	$(CXX) $(OMPFLAGS) $^ -o $@

clean:
	rm -rf *.o simulation

.PHONY: clean

