
MAKERULES?=intel
include make.defs.$(MAKERULES)
CPPFLAGS = -I./ -I./tictoc
DEFINES =
RM = rm -f
LIBS = -lm -lrt
LINKER = $(CXX) $(OMPFLAGS)

all: cloud

%.o: %.cc
	$(CXX) $(CXXVERFLAGS) $(CXXINFOFLAGS) $(CXXOPTFLAGS) $(CPPFLAGS) $(DEFINES) $(OMPFLAGS) -c -o $@ $<

cloud: timestep.o accelerate.o cloud.o state.o json.o params.o
	$(LINKER) -o $@ $^ $(LIBS)

PERF =
STREAM_N = 1000000

NP=256
DT=1.e-4
K=100.
D=1.
L=20.
R=1.
NT=1000
NC=$(NT)
SEED=6230

runcloud:
	$(PERF) ./cloud -Np $(NP) -Nt $(NT) -dt $(DT) -k $(K) -d $(D) -L $(L) -r $(R) -Nc $(NC) $(DIRECT) -seed $(SEED)

reproducibility_check: cloud
	if [ `(for i in {1..24}; do OMP_NUM_THREADS=$$i ./cloud -Np 16384 -Nt 1 | grep "mean estimated diffusivity"; done) | uniq | wc -l` == "1" ]; then \
		echo "reproducibility check passed"; \
	else \
		echo "reproducibiilty check failed: different coefficients computed with different thread counts" && exit -1; \
	fi

timings.json: cloud
	rm -f $@ && touch $@ && delim="[" && \
	for N_p in 1 2 4 8 16 32 64; do \
	  echo $$delim >> $@; \
		this_L=`echo "$$N_p 0.3333 20." | awk '{ print ($$3 * $$1^$$2); }'`; \
		this_T=`echo "$$N_p 100000" | awk '{ print ($$2 / ($$1 * $$1)); }'`; \
		OMP_NUM_THREADS=$$PBS_NP OMP_PROC_BIND=spread ./cloud -Np $$(( 256 * $$N_p )) -L $$this_L -Nt $$this_T >> $@; \
		delim=","; \
	done && echo "]" >> $@

timings.png: timings.json
	python3 timings.py

coefficient_test.json: cloud
	OMP_NUM_THREADS=$$PBS_NP OMP_PROC_BIND=spread ./cloud -Np 16384 -Nt 20000 -Nc 10000 -L 80. > $@

coefficient_test: coefficient_test.json
	python3 coefficient_test.py


clean:
	$(RM) *.o *.optrpt cloud timings.json coefficient_test.json timings.png

.PHONY: clean runcloud reproducibility_check
