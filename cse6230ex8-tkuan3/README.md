
# Exercise 7: MPI stencil performance

**Due:** November 12, end of day.

**Collaboration:** you may work in teams of two.  Your team may consult online
resources, others in the class, and the instructor and TA, but your work must
be your own.  Copying other students' code will be considered an honor code
violation.

## About the code

The goals of this assignment are two-fold:

1. Practice correctly implementing a non-trivial MPI communication task -- the halo exchange for a 3D stencil operation.
2. Practice gathering data to compare competing implementations of the same operation

The code implements a simple partial differential equation solver for [Poisson's equation](https://en.wikipedia.org/wiki/Poisson's_equation)
on a 3D box domain.  The numerical details are not important: we only need to know that it works by iteratively applying
a Jacobi smoother (or "Jacobi sweep" as we have called it previously in the class).

The program runs the following way

```
mpiexec -n P ./main -N N -T T -S S [-var]
```

- `P`: it should be valid for any number of processors
- `N`: the number of grid points in each of the three directions
- `T`: the number of batches of Jacobi sweeps that are applied (three convergence calculations -- the [residual](https://en.wikipedia.org/wiki/Residual_(numerical_analysis)),
       the error (the difference between the computed and exact solution), and the [convergence factor](https://en.wikipedia.org/wiki/Rate_of_convergence).
- `S`: the number of Jacobi sweeps in a batch.
- `var`: a switch to choose a _variant_ algorithm (more on this below)

The data that is distributed on the grid is represented by a `GridFunction` class.

The code you have received is a correct serial program.  It is also a partially correct parallel program: in parallel it can

1. Divide up the 3D grid into a logically Cartesian layout of processes, each with its own slice of the distributed grid.
2. Verify that the exact solution is correct.

The solver is not correct, because the **halo** exchange is not implemented yet.

## Your task, part 1

You must implement the methods `GridFunction::halo_start()` and `GridFunction::halo_end()`.

- `halo_start()`: Calling this function means that the values in the interior of the local grid are up-to-date.
   This function should be non-blocking. The data values in the `GridFunction` should not be modified
   and the halo values are not guaranteed to be current when this function returns.

- `halo_end()`: Calling this functions means that you are ready to read the halo values from other processes.
   This function should be blocking until it is safe to modify the values in the `GridFunction` without affecting the
   halo values send to other processes, and until the halo values received from other processes are up-to-date.

How will you know your implementation is correct?

The data for this problem is set up in such a way that the _convergence factor_ -- the rate at which the computed
solution converges to the true solution -- is a fixed constant that can be predicted exactly.  If the `convergence factor history`
of your implementation matches the `Predicted convergence factor`, then the halo exchange is working properly.

### Guidance

You should construct the data structures that you will need for halo exchange
--- send / receive buffers, `MPI_Win`s, lists of neighboring processes -- by
adding additional fields to `GridFunction` and initializing them when the
`GridFunction` is created (either directly in `GridFunction::GridFunction()` or
in the copy constructor `GridFunction::GridFunction(const GridFunction& F)`.
That way the halo exchange includes only the actual communication cost.

## Your task, part 2

As discussed in class, there are many potential optimizations for stencil operations which could plausibly improve the performance.
You will **implement two variants** of halo exchange and gather measurements to compare them.

You can choose variants from any of the following pairs:

1. **packing** and **unpacking** send and receive buffers vs. **in place buffers with MPI vector types**
2. **alternating communication and computation** (calling `halo_start()` and `halo_end()` without any intermediate grid calculations) vs.
   **overlapping communication and computation** (performing the Jacobi sweep operations that do not depend on the halo values
   between the calls to `halo_start()` and `halo_end()`)
3. **two-sided communication** (using symmetric sends and receives) vs. **one-sided communication** (using remote memory access operations)

For `P = 64`, you will create a job script to collect data for varying sizes of `N` (from as small as `N = 8` to as large as you can run in a reasonable
amount of time).  You will plot the data that you have collected and see if you detect any trends.

## Grading

1. Correctness (50%): the existing test script `test.pbs` runs to completion, and the `convergence factor history` matches the `predicted convergence factor`
   for each of your two variants (25% each).  Points may be deducted for other MPI usage errors in your code that did not affect the convergence factor.

2. Reasonableness (10%): as the Guidance above says, the implementations of your two variants should have their setup take place during construction,
   not during the calls to the halo methods.

3. Reproducibility (20%): your repository should contain the file `figures.pbs` that recreates the data for the figures in your report.

4. Report (20%): Indicate which variants you tested.  Present figures comparing the variants as discussed above (10%).  If there is a robust 
   trend that shows the variants are different, do your best to explain why;
   otherwise, suggest a potential future test that could better differentiate
   them (10%).
