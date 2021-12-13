
# Exercise 3: n-body charged particle simulation

**Due:** September 20, before class.

**Collaboration:** you may work in teams of two.  Your team may consult online
resources, others in the class, and the instructor and TA, but your work must
be your own.  Copying other students' code will be considered an honor code
violation.

## Step 1 (5%): Getting the code

This project is hosted in a repository on <https://github.gatech.edu/cse6230fa21/cse6230ex3>.

Your team should clone the repository into a local copy that you will use to develop.

## About the code

The math of this simulation is very similar to the asteroid simulation from the last assignment, but now
we have charged particles that repel each other, and they are trapped in a box.


The program `./main PRECISION N_STEPS N_CHARGES_PER_DIM` simulates `N_STEPS` of `N_CHARGES_PER_DIM^3` particles.

As before, each asteroid is represented by the six components of its state vector:

- its position `x`, `y`, `z`
- its velocity `u`, `v`, `w`

This time, however, each particle acts on every other particle: this is the starter code
for the computation of the acceleration that updates the velocity components:

```c
for (size_t i = 0; i < n_charges; i++) {
  Real d_u = 0;
  Real d_v = 0;
  Real d_w = 0;
  for (size_t j = 0; j < n_charges; j++) {
    if (i == j) {
      continue;
    }
    Real diff_x = x[i] - x[j];
    Real diff_y = y[i] - y[j];
    Real diff_z = z[i] - z[j];
    Real r2 = diff_x * diff_x + diff_y * diff_y + diff_z * diff_z;
    Real r = std::sqrt(r2);
    Real factor = dt / (r2 * r);

    d_u += factor * diff_x;
    d_v += factor * diff_y;
    d_w += factor * diff_z;
  }
  u[i] += d_u;
  v[i] += d_v;
  w[i] += d_w;
}
```

The performance of an implementation is measured in **charged particle interactions per second**.

## Twin directories

This assignment wll be run on the `coc-ice-gpu` nodes, where we have the opportunity to compare OpenMP accleration on the host cores and OpenACC acceleration on an NVIDA V100 GPU.
You start with two identical versions of the code: one for each method of acceleration.

## Step 2 (30%): Parallelize the code with OpenMP

In the `openmp/` directory, you should optimize the performance of the calculation using OpenMP.

Your code should:

1. Be correct.
2. Use optimization flags that allow for floating-point reorganization into SIMD instructions (make these changes in `openmp/make.defs`)
3. Parallelize the loops in `timestep.cc`.
4. Parallelize the initialization in `state.cc` (You should use the `collapse()` OpenMP directive where it would make sense to do so.)
5. Minimize the overhead of thread synchronization at the start and end of a parallel region.
6. Explicitly control the `schedule()` of all loops.
7. Control the affinity of threads using the `OMP_PROC_BIND` environment variable (make this change in `openmp/omp_data.sh`)

You are welcome to make additional optimizations changes.

## Step 3 (30%): Parallelize the code with OpenACC

In the `acc/` directory, you should optimize the performance of the calculation using OpenACC.

Your code should:

1. Be correct.
2. Parallelize the loops in `timestep.cc`.
3. Minimize the overhead of data movement to and from the GPU.

You are welcome to make additional optimizations changes.

## Running calculations on coc-ice compute nodes

You will run your simulations on the `coc-ice-gpu` queue compute nodes.  You
should submit the `data.pbs` job script from the top directory of this project to generate data for your figure.

Note that for un-optimized serial code the script may run out of time, but once you optimize and parallel your code,
this should not be an issue.

**Do not run your code on the login nodes.**  If you want to develop interactively, get an interactive session.

## Step 4 (20%): make figures

Make two figures from your data: one for `float`, and one for `double`.

Your figure will compare the work rate (charged particle interactions per second) of the OpenMP and OpenACC implementations as a function of the number of charged particles.

Use the plotting software of your choice to create the figures.

## Step 5 (15%): make a report

Make a report `report.pdf` that includes your two figures.  You should also comment on the performance of your OpenMP and OpenACC code: why is one faster than the other?  A good answer will refer to both the capabilities of the GPU vs. the capabilities of the host, as well as the nature of this n-body problem.

## submit your work

- Check in to your repository your `report.pdf`, as well as any additional files you may have added that are necessary to build your program.

- Create one **private** repository for your team on `github.gatech.edu` and push your final version there.  The name of your repository should look like `cse6230ex2-gtusername1` or `cse6230ex2-gtusername1-gtusername2` for a team.

- When you are ready to submit, [transfer ownership](https://docs.github.com/en/github/administering-a-repository/managing-repository-settings/transferring-a-repository) of your repository to the group `cse6230fa21`.

- On canvas, your "submission" will just be the name of your repository.

- The TA will submit your `data.pbs` script to verify your code.


