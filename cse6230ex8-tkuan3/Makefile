include make.defs

all: main

MPICXX=mpicxx

%.o: %.cc
	$(MPICXX) -c $(CXXVERFLAGS) $(CXXOPTFLAGS) $(CXXINFOFLAGS) -o $@ $<

exact.o: exact.h params.h json.h

solver.o: solver.h grid.h params.h json.h

main.o: grid.h solver.h params.h exact.h json.h

json.o: json.h

params.o: params.h json.h

main: main.o grid.o json.o params.o exact.o solver.o
	$(MPICXX) $(CXXVERFLAGS) $(CXXOPTFLAGS) $(CXXINFOFLAGS) -o $@ $^

clean:
	rm -f main *.o *.optrpt

.PHONY: clean

