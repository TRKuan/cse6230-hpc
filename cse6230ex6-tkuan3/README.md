
# Exercise 6: MPI ring test

**Due:** October 20, before class.

**Collaboration:** you may work in teams of two.  Your team may consult online
resources, others in the class, and the instructor and TA, but your work must
be your own.  Copying other students' code will be considered an honor code
violation.

## About the code

The goal of this exercise to get your first practice using the MPI interface,
and to time message latencies in different settings using a ring test.

MPI processes are given ranks [0,size) in the world communicator, where `size`
is requested number of processes in `mpiexec -n size`.

A ring test is when rank 0 sends a message to rank 1, which sends it to rank 2, and so on.
The process with rank (size-1) sends it back to rank 0. to complete the ring.

Repeating trips around the ring can give you a good idea of the latency of sending small messages
from one process to the next.

**Your task:** implement the ring test and helper functions labeled `TODO` in `ring-test.cc`.

The program `ring-test` is run using `mpiexec`: for example, when run with 12 processes, your
correcty implemented program should print something like the following:

```bash
$ mpiexec -f $PBS_NODEFILE -n 12 ./ring-test 2 v
MVAPICH2 Version      : 2.3.2
MVAPICH2 Release date : Fri August 9 22:00:00 EST 2019
MVAPICH2 Device       : ch3:mrail
MVAPICH2 configure    : --prefix=/usr/local/pace-apps/spack/packages/0.13/linux-rhel7-cascadelake/intel-19.0.5/mvapich2-2.3.2-hpgbkqoytbjh35qn2t63rdorepxcezek --enable-shared --enable-romio --disable-silent-rules --disable-new-dtags --enable-fortran=all --enable-threads=multiple --with-ch3-rank-bits=32 --enable-wrapper-rpath=yes --disable-alloca --enable-fast=all --disable-cuda --enable-registration-cache --with-pbs=/opt/torque/current --with-device=ch3:mrail --with-rdma=gen2 --disable-mcast --with-file-system=nfs+ufs
MVAPICH2 CC           : /var/pace-apps/spack/root/0.13/458c9a2/lib/spack/env/intel/icc    -DNDEBUG -DNVALGRIND -O2
MVAPICH2 CXX          : /var/pace-apps/spack/root/0.13/458c9a2/lib/spack/env/intel/icpc   -DNDEBUG -DNVALGRIND -O2
MVAPICH2 F77          : /var/pace-apps/spack/root/0.13/458c9a2/lib/spack/env/intel/ifort -L/lib -L/lib   -O2
MVAPICH2 FC           : /var/pace-apps/spack/root/0.13/458c9a2/lib/spack/env/intel/ifort   -O2

MPI # Procs: 12
MPI Wtime 1.63414e+09, precision 1e-09
MPI Wtime is global
MPI rank 0 of 12 on host: atl1-1-02-003-20-l.pace.gatech.edu
MPI rank 1 of 12 on host: atl1-1-02-003-20-l.pace.gatech.edu
MPI rank 2 of 12 on host: atl1-1-02-003-20-l.pace.gatech.edu
MPI rank 3 of 12 on host: atl1-1-02-003-20-l.pace.gatech.edu
MPI rank 4 of 12 on host: atl1-1-02-003-20-l.pace.gatech.edu
MPI rank 5 of 12 on host: atl1-1-02-003-20-l.pace.gatech.edu
MPI rank 6 of 12 on host: atl1-1-02-003-20-l.pace.gatech.edu
MPI rank 7 of 12 on host: atl1-1-02-003-20-l.pace.gatech.edu
MPI rank 8 of 12 on host: atl1-1-02-003-20-l.pace.gatech.edu
MPI rank 9 of 12 on host: atl1-1-02-003-20-l.pace.gatech.edu
MPI rank 10 of 12 on host: atl1-1-02-003-20-l.pace.gatech.edu
MPI rank 11 of 12 on host: atl1-1-02-003-20-l.pace.gatech.edu
MPI rank 0 of 12, message: 0
MPI rank 1 of 12, message: 1
MPI rank 2 of 12, message: 2
MPI rank 3 of 12, message: 3
MPI rank 4 of 12, message: 4
MPI rank 5 of 12, message: 5
MPI rank 6 of 12, message: 6
MPI rank 7 of 12, message: 7
MPI rank 8 of 12, message: 8
MPI rank 9 of 12, message: 9
MPI rank 10 of 12, message: 10
MPI rank 11 of 12, message: 11
MPI rank 0 of 12, message: 12
MPI rank 1 of 12, message: 13
MPI rank 2 of 12, message: 14
MPI rank 3 of 12, message: 15
MPI rank 4 of 12, message: 16
MPI rank 5 of 12, message: 17
MPI rank 6 of 12, message: 18
MPI rank 7 of 12, message: 19
MPI rank 8 of 12, message: 20
MPI rank 9 of 12, message: 21
MPI rank 10 of 12, message: 22
MPI rank 11 of 12, message: 23
MPI rank 0 of 12, message: 24
MPI rank 1 of 12, message: 25
MPI rank 2 of 12, message: 26
MPI rank 3 of 12, message: 27
MPI rank 4 of 12, message: 28
MPI rank 5 of 12, message: 29
MPI rank 6 of 12, message: 30
MPI rank 7 of 12, message: 31
MPI rank 8 of 12, message: 32
MPI rank 9 of 12, message: 33
MPI rank 10 of 12, message: 34
MPI rank 11 of 12, message: 35
{
 "size": 12,
 "number of trips": 3,
 "time (seconds)": 0.00248766,
 "seconds per message": 6.91017e-05
}
```

## Grading

- Correctness (70%):
  - Your code must send the message around the ring the requested number of times, incrementing the counter with each message. (40%)
  - Your code must, in verbose mode, correctly send the message value to the root processor to print, so that your output looks like the output above (30%)

- Figure and report (30%):
  - Run `test.pbs` on your code to generate timings of the average message time for different number of processes in `ring_test_X.json`.
    Generate a figure of the seconds per message for each number of processes (20%)
  - The script also generates `scatter_ring_test_X.json`, which runs the ring test with four processes on four separate coc-ice nodes.
    Compare the average seconds per message from that test to the values in `ring_test_X.json` and answer the question: why are
    these times different? (10%)
