
# Exercise 5: median of medians

**Due:** October 13, before class.

**Collaboration:** you may work in teams of two.  Your team may consult online
resources, others in the class, and the instructor and TA, but your work must
be your own.  Copying other students' code will be considered an honor code
violation.

## About the code

The goal of this project is to familiarize yourself with using parallel algorithm libraries
like [thrust](https://thrust.github.io).

This code implements quicksort using median of medians to find the pivot.  But, the median finding algorithm
is currently implemented on the host CPUs.

**Your task:** implement `select_device` by adapting `select_host` using
thrust.  You adaptation should work so that the elements in the array to be
sorted never leave the device.

The logic related to how big the subarrays are and which subarrays to use for
recursive calls will still be run in host code:

**Note:** you should `module purge && module load nvhpc gcc/9.2.0.lua` to get the right compilers for this assignment.

## Working with `thrust::device_vector<int>`

Our container on the host that manipulates the array on the device is `thrust::device_vector<int> A`.

Because thrust tries not to copy data back and forth from the host to device unneccessarily, operations
on a device vector are a little different. Here are some tips for working on the assignment:

- `A[i]` is a reference to value on the device: it has type `thrust::device_reference<int>`.
- Note: if you want, for example to swap `A[0]` and `A[pivot]`, you can use `thrust::swap(A[0],A[pivot])` to avoid copying either value to the host
- `&A[i]` is a `thrust::device_ptr<int>`.
- You can make a zip iterator from device pointers as follows: `thrust::make_zip_iterator(thrust::make_tuple(A[i0], A[i1], ...))`.
- Notice that I have defined several `iter_x_ref` types, which would be the iteration type for that kind of zip iterator.
- In addition to the algorithms that we have already discussed in class, you may want to use [`for_each`](https://thrust.github.io/doc/group__modifying_ga14cba62489aee67ffa6348eb74997b57.html) or [`for_each_n`](https://thrust.github.io/doc/group__modifying_ga00ad46c06f41dd29601d44d3c6c25819.html#ga00ad46c06f41dd29601d44d3c6c25819)
- Remember: lambdas that run on the device look like `[] __device__ (args) -> output { body }`.

## Grading

- Correctness (50%): Your code must properly sort all arrays that our TA tests it on.
- Use of the device (25%): If all thrust algorithms in `select_device` use `thrust::device` execution policy.
- No transfers (25%): If no values out of the array (including pivot values!) are copied from the device to the host during sorting.
