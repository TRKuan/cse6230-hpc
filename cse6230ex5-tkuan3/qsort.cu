#include <iostream>
#include <thrust/device_malloc.h>
#include <thrust/device_free.h>
#include <thrust/device_vector.h>
#include <thrust/transform_reduce.h>
#include <thrust/random.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/tuple.h>
#include <thrust/execution_policy.h>
#include <curand.h>

#define DEBUG 0

typedef thrust::tuple<int, int> iter_pair;
typedef thrust::tuple<int &, int &> iter_pair_ref;
typedef thrust::tuple<int &, int &, int &> iter_trip_ref;
typedef thrust::tuple<int &, int &, int &, int &> iter_quart_ref;
typedef thrust::tuple<int &, int &, int &, int &, int &> iter_quint_ref;

/// fill the array with random numbers
void generate(thrust::device_vector<int> &array) {
  curandGenerator_t gen;
  curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);
  curandSetPseudoRandomGeneratorSeed(gen, 0);
  int * A_start_raw = (&array[0]).get();
  curandGenerate(gen, (unsigned int *) A_start_raw, array.size());
} 

/// fill the array with 0..N
void generate_sorted(thrust::device_vector<int> &array) {
  thrust::counting_iterator<int> first(0);
  thrust::copy(first, first + array.size(), array.begin());
}

/// test if sorted
bool is_sorted(thrust::device_vector<int> &array) {
  size_t N = array.size();
  if (N < 2) {
    return true;
  }
  // get references to the first and second elements
  auto first = array.begin();
  auto second = next(first);
  // get references to the second-to-last and last elements
  auto last = array.end();
  auto beforelast = prev(last);

  // loop over consecutive pairs by making a zip iterator
  // that looks at pairs (array[0], array[1]), (array[1], array[2]),
  //
  // each pair (array[i], array[i+1]) evaluates to true
  // if array[i] <= array[i+1]
  // 
  // the list is sorted if this is true for all pairs, so we
  // reduce the {true,false} values with a logical and
  bool sorted = thrust::transform_reduce(
      thrust::device,
      thrust::make_zip_iterator(thrust::make_tuple(first, second)),
      thrust::make_zip_iterator(thrust::make_tuple(beforelast, last)),
      [] __device__ (iter_pair pair) -> bool {
        return thrust::get<0>(pair) <= thrust::get<1>(pair);
      },
      true,
      thrust::logical_and<bool>());
  return sorted;
}

// comparator swap
template<typename T>
void __host__ __device__ cswp(T &a1, T &a2) {
  T min = a1 <= a2 ? a1 : a2;
  T max = a1 <= a2 ? a2 : a1;
  a1 = min;
  a2 = max;
}

// see <https://stackoverflow.com/questions/3903086/standard-sorting-networks-for-small-values-of-n>
template<typename T>
void __host__ __device__ sort_3(T &a1, T &a2, T &a3) {
  cswp<T>(a1, a2);
  cswp<T>(a1, a3);
  cswp<T>(a2, a3);
}

template<typename T>
void __host__ __device__ sort_4(T &a1, T &a2, T &a3, T &a4) {
  cswp<T>(a1, a2); cswp<T>(a3, a4);
  cswp<T>(a1, a3); cswp<T>(a2, a4);
  cswp<T>(a2, a3);
}

template<typename T>
void __host__ __device__ sort_5(T &a1, T &a2, T &a3, T &a4, T &a5) {
  cswp<T>(a1, a2); cswp<T>(a3, a4);
  cswp<T>(a1, a3); cswp<T>(a2, a5);
  cswp<T>(a1, a2); cswp<T>(a3, a4);
  cswp<T>(a2, a3); cswp<T>(a4, a5);
  cswp<T>(a3, a4);
}

/// partition so that the kth element (by order) is in position k
void select_host(int *A, int *A_end, size_t k) {
  size_t n = A_end - A;
  size_t n_div_5 = n / 5;
  size_t n_rem_5 = n % 5;

  int *subA[5];
  size_t offsets[6];

  offsets[0] = 0;

  // break up A into 5 sub-array of equal length
  for (size_t i = 0; i < 5; i++) {
    subA[i] = &A[i * n_div_5 + (i * n_rem_5) / 5];
    offsets[i+1] = (i+1) * n_div_5 + ((i+1) * n_rem_5) / 5;
  }

#if DEBUG
  std::cout << "Pre " << std::endl;
  for (size_t j = 0; j < 5; j++) {
    for (size_t i = offsets[j]; i < offsets[j+1]; i++) {
      std::cout << A[i] << " ";
    }
    std::cout << std::endl;
  }
#endif

  // sort matching elements of the five arrays
  for (size_t j = 0; j < n_div_5; j++) {
    sort_5<int>(subA[0][j], subA[1][j], subA[2][j], subA[3][j], subA[4][j]);
  }

  // handle the remainders when some of the arrays are longer than others
  switch (n_rem_5) {
  case 4:
    sort_4<int>(subA[1][n_div_5], subA[2][n_div_5], subA[3][n_div_5], subA[4][n_div_5]);
    break;
  case 3:
    sort_3<int>(subA[1][n_div_5], subA[3][n_div_5], subA[4][n_div_5]);
    break;
  case 2:
    cswp<int>(subA[2][n_div_5], subA[4][n_div_5]);
    break;
  default:
    break;
  }

#if DEBUG
  std::cout << "After remainder " << std::endl;
  for (size_t j = 0; j < 5; j++) {
    for (size_t i = offsets[j]; i < offsets[j+1]; i++) {
      std::cout << A[i] << " ";
    }
    std::cout << std::endl;
  }
#endif

  if (n <= 5) {
    return;
  }

  // which array should have element k?
  size_t subindex = 0;
  if (k >= offsets[1]) {subindex++;}
  if (k >= offsets[2]) {subindex++;}
  if (k >= offsets[3]) {subindex++;}
  if (k >= offsets[4]) {subindex++;}

  // what should its location be?
  size_t subn = offsets[subindex+1] - offsets[subindex];
  size_t subk = k - offsets[subindex];
  select_host(subA[subindex], subA[subindex] + subn, subk);
  int pivot = subA[subindex][subk];

  // divide the partition into less than pivot and greater or equal to the pivot
  int *A_l = A;
  int *A_eq = thrust::partition(
      thrust::host,
      A, A_end,
      [pivot] __host__ (int v) -> bool {
        return v < pivot;
      }
      );
  int *A_g = thrust::partition(
      thrust::host,
      A_eq, A_end,
      [pivot] __host__ (int v) -> bool {
        return v == pivot;
      }
      );

  size_t n_lo = A_eq - A_l;
  size_t n_eq = A_g - A_eq;
  size_t n_hi = A_end - A_g;

  if (n_lo <= k && k < n_lo + n_eq) {
    return;
  } else if (k < n_lo) {
    select_host(A_l, A_l + n_lo, k);
  } else {
    select_host(A_g, A_g + n_hi, k - (n_lo + n_eq));
  }
#if DEBUG
  std::cout << "Final " << std::endl;
  for (size_t i = 0; i < k; i++) {
    std::cout << A[i] << " ";
  }
  std::cout << std::endl;
  std::cout << A[k] << std::endl;
  for (size_t i = k+1; i < n; i++) {
    std::cout << A[i] << " ";
  }
  std::cout << std::endl;
#endif
}

void select_device(thrust::device_ptr<int> A, thrust::device_ptr<int> A_end, size_t k) {
  size_t n = A_end - A;
  size_t n_div_5 = n / 5;
  size_t n_rem_5 = n % 5;
  
  thrust::device_ptr<int> subA[5];
  size_t offsets[6];

  offsets[0] = 0;

  // break up A into 5 sub-array of equal length
  for (size_t i = 0; i < 5; i++) {
    subA[i] = &A[i * n_div_5 + (i * n_rem_5) / 5];
    offsets[i+1] = (i+1) * n_div_5 + ((i+1) * n_rem_5) / 5;
  }

  // sort matching elements of the five arrays
  thrust::for_each_n(
    thrust::device,
    thrust::make_zip_iterator(thrust::make_tuple(&subA[0][0], &subA[1][0], &subA[2][0], &subA[3][0], &subA[4][0])),
    n_div_5,
    [] __device__ (auto pair) -> void {
      sort_5<int>(thrust::get<0>(pair), thrust::get<1>(pair), thrust::get<2>(pair), thrust::get<3>(pair), thrust::get<4>(pair));
    }
  );

  // handle the remainders when some of the arrays are longer than others
  thrust::for_each_n(
    thrust::device,
    thrust::make_zip_iterator(thrust::make_tuple(&subA[0][n_div_5], &subA[1][n_div_5], &subA[2][n_div_5], &subA[3][n_div_5], &subA[4][n_div_5])),
    1,
    [n_rem_5] __device__ (auto pair) -> void {
      switch (n_rem_5) {
        case 4:
          sort_4<int>(thrust::get<1>(pair), thrust::get<2>(pair), thrust::get<3>(pair), thrust::get<4>(pair));
          break;
        case 3:
          sort_3<int>(thrust::get<1>(pair), thrust::get<3>(pair), thrust::get<4>(pair));
          break;
        case 2:
          cswp<int>(thrust::get<2>(pair), thrust::get<4>(pair));
          break;
        default:
          break;
      }
    }
  );

  if (n <= 5) {
    return;
  }

  // which array should have element k?
  size_t subindex = 0;
  if (k >= offsets[1]) {subindex++;}
  if (k >= offsets[2]) {subindex++;}
  if (k >= offsets[3]) {subindex++;}
  if (k >= offsets[4]) {subindex++;}

  // what should its location be?
  size_t subn = offsets[subindex+1] - offsets[subindex];
  size_t subk = k - offsets[subindex];
  select_device(subA[subindex], subA[subindex] + subn, subk);
  //int pivot = subA[subindex][subk];
  thrust::device_ptr<int> pivot = thrust::device_malloc<int>(1);
  *pivot = subA[subindex][subk];

  // divide the partition into less than pivot and greater or equal to the pivot
  thrust::device_ptr<int> A_l = A;
  thrust::device_ptr<int> A_eq = thrust::partition(
      thrust::device,
      A, A_end,
      [pivot] __device__ (int v) -> bool {
        return v < *pivot;
      }
      );
  thrust::device_ptr<int> A_g = thrust::partition(
      thrust::device,
      A_eq, A_end,
      [pivot] __device__ (int v) -> bool {
        return v == *pivot;
      }
      );
  
  thrust::device_free(pivot);

  size_t n_lo = A_eq - A_l;
  size_t n_eq = A_g - A_eq;
  size_t n_hi = A_end - A_g;

  if (n_lo <= k && k < n_lo + n_eq) {
    return;
  } else if (k < n_lo) {
    select_device(A_l, A_l + n_lo, k);
  } else {
    select_device(A_g, A_g + n_hi, k - (n_lo + n_eq));
  }
}

void select(thrust::device_ptr<int> A, thrust::device_ptr<int> A_end, size_t k) {
#if 0
  size_t n = A_end - A;

  thrust::host_vector<int> A_host(n);

  thrust::copy(A, A_end, &A_host[0]);

  select_host(&A_host[0], &A_host[0] + n, k);

  thrust::copy(&A_host[0], &A_host[0] + n, A);
#else
  select_device(A, A_end, k);
#endif
}

void qsort(thrust::device_ptr<int> A, thrust::device_ptr<int> A_end) {
  size_t n = A_end - A;

  if (n <= 2) {
    return;
  }
  select(A, A_end, n / 2);
  qsort(A, &A[n/2]);
  qsort(&A[n/2], A_end);
}

int main(int argc, char **argv) {
  size_t N = 1000000;

  if (argc > 1) {
    N = (size_t) strtoull(argv[1], NULL, 10);
  }
  std::cout << "Sorting " << N << " ints" << std::endl;

  thrust::device_vector<int> array(N);
  generate(array);

  bool input_sorted = is_sorted(array);
  std::cout << "Input is " << (input_sorted ? "" : "not ") << "sorted" << std::endl;

  qsort(&array[0], &array[N]);

  bool output_sorted = is_sorted(array);
  std::cout << "Output is " << (output_sorted ? "" : "not ") << "sorted" << std::endl;
  return 0;
}
