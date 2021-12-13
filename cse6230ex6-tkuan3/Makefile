include make.defs.intel

MPICXX=mpicxx

ring-test: ring-test.cc
	$(MPICXX) $(CXXVERFLAGS) $(CXXOPTFLAGS) $(CXXINFOFLAGS) -o $@ $<

clean:
	rm -f ring-test

.PHONY: clean

