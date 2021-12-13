
# Exercise 7: MPI algorithm construction

**Due:** October 27, before class.

**Collaboration:** you may work in teams of two.  Your team may consult online
resources, others in the class, and the instructor and TA, but your work must
be your own.  Copying other students' code will be considered an honor code
violation.

## About the code

The goal of this exercise to see how the MPI collective operations correspond
to the std::algorithms we have already encountered and to learn how to build
higher-level interfaces with those collectives.

You must implement four algorithms.

1. Distributed array reduce in `allreduce_vec()`: Compute the reduction of all
   of the elements of a distributed array and store the result on all processes.
   Edit this function in `reduce.cc`.

2. Distributed array prefix sum in `inclusive_prefix_sum()`: Given a distributed
   array `vec`, compute the distributed prefix sum of its values and store the result
   in `result`, which has the same layout.
   Edit this function in `scan.cc`.

3. Distributed matrix transpose in `transpose_1D_col()`: the matrix `A` has
   `nrows_global` rows and `ncols_global` columns.  It is distributed in a 1D fashion
   by distributing the columns, in order, to each process.
   The local portion of the matrix is stored in column-major order.
   
   You must store the transpose of `A` in the distributed matrix `AT` which has
   `ncols_global` rows and `nrows_global` columns, and which is
   column-distributed in the same way as `AT`.

   The number of processes divides `ncols_global` and `nrows_global`, so each process has the same
   number of columns of `A` and of `AT`.

   Edit this function in `transpose_1d.cc`.

4. 2D distributed matrix transpose in `transpose_2D()`: the matrix `A` is distributed
   on a square 2D Cartesian grid of processes, created with `MPI_Cart_create()`.
   Each process has a slice of `nrows_local` rows and `ncols_local` columns of `A`, stored
   in column-major order.  The rank with Cartesian coordinate `(i,j)` has a slice of rows
   `i*nrows_local:(i+1)*nrows_local` and columns `j*ncols_local:(j+1)*ncols_local`.

   You must store the transpose of `A` in the distributed matrix `AT` which is also
   distributed on a square 2D Cartesian grid of processes. Each process has `ncols_local` columns
   and `nrows_local` rows.

   Edit this function in `transpose_2d.cc`.

Your functions are tested in `main.cc`, and compared against reference implementations that only
use scatters, gathers and broadcasts for communication.

**Your goal** is that use only one MPI communication routine in each of your functions.

## Grading

- Correctness (12.5% per function): your functions should pass the `test.pbs`
  without errors.  Points may be deducted for correctness problems found by
  inspecting the code that did not cause the test to fail.

- Minimal communication (12.5% per function): In each of your functions, each
  process should participate in at most one collective call or at most one send
  and one receive.
