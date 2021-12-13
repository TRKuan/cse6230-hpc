include make.defs

MPICXX=mpicxx

%.o: %.cc
	$(MPICXX) -c $(CXXVERFLAGS) $(CXXOPTFLAGS) $(CXXINFOFLAGS) -o $@ $<

main: main.o reduce.o scan.o transpose_1d.o transpose_2d.o util.o
	$(MPICXX) $(CXXVERFLAGS) $(CXXOPTFLAGS) $(CXXINFOFLAGS) -o $@ $^

clean:
	rm -f main *.o *.optrpt

.PHONY: clean

