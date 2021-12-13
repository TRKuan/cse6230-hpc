
#include <cstdio>
#include <cstdlib>
#include <mpi.h>

#define MPI_CHK(err) if (err != MPI_SUCCESS) return err

const int nametag = 1;

int print_host(void) {
  int err;
  char procname[MPI_MAX_PROCESSOR_NAME] = {0};
  int rank, size;
  int namelen;
  MPI_Comm comm = MPI_COMM_WORLD;

  err = MPI_Get_processor_name(procname, &namelen); MPI_CHK(err);
  err = MPI_Comm_size(comm, &size); MPI_CHK(err);
  err = MPI_Comm_rank(comm, &rank); MPI_CHK(err);
  printf("MPI rank %d of %d on host: %s\n", rank, size, procname);
  return 0;
}

int send_host_to_root(void) {
  int err;
  char procname[MPI_MAX_PROCESSOR_NAME] = {0};
  int rank, size;
  int namelen;
  MPI_Comm comm = MPI_COMM_WORLD;

  err = MPI_Get_processor_name(procname, &namelen); MPI_CHK(err);
  err = MPI_Send(procname, namelen, MPI_CHAR, 0, nametag, comm); MPI_CHK(err);
  return 0;
}

int print_host_from_other(void) {
  int err;
  char procname[MPI_MAX_PROCESSOR_NAME] = {0};
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Status status;

  int size;
  err = MPI_Comm_size(MPI_COMM_WORLD, &size); MPI_CHK(err);
  err = MPI_Recv(procname, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, MPI_ANY_SOURCE, nametag, comm, &status); MPI_CHK(err);
  printf("MPI rank %d of %d on host: %s\n", status.MPI_SOURCE, size, procname);
  return 0;
}

int print_mpi_info(void) {
  char library[MPI_MAX_LIBRARY_VERSION_STRING] = {0};
  int version, subversion, liblen, err;
  double time, timeprec;

  err = MPI_Get_version(&version, &subversion); MPI_CHK(err);
  err = MPI_Get_library_version(library, &liblen); MPI_CHK(err);

  int rank;
  err = MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_CHK(err);
  int size;
  err = MPI_Comm_size(MPI_COMM_WORLD, &size); MPI_CHK(err);


  if (rank == 0) {
    printf("MPI Version: %d.%d\n", version, subversion);
    printf("%s\n", library);
    printf("MPI # Procs: %d\n", size);

    time = MPI_Wtime();
    timeprec = MPI_Wtick();

    printf("MPI Wtime %g, precision %g\n", time, timeprec);
    if (MPI_WTIME_IS_GLOBAL) {
      printf("MPI Wtime is global\n");
    } else {
      printf("MPI Wtime is not global\n");
    }
  }

  if (rank == 0) {
    print_host();
    for (int i = 1; i < size; i++) {
      print_host_from_other();
    }
  } else {
    send_host_to_root();
  }
  return 0;
}

const int msgtag = 2;

int print_message(int message) {
  int err;
  int rank, size;
  MPI_Comm comm = MPI_COMM_WORLD;

  err = MPI_Comm_size(comm, &size); MPI_CHK(err);
  err = MPI_Comm_rank(comm, &rank); MPI_CHK(err);
  printf("MPI rank %d of %d, message: %d\n", rank, size, message);
  return 0;
}

// TODO: instead of having each process print their message for themselves (which
// can lead printing at the same time), have each process send its message to the
// host to print
int send_message_to_root(int message) {
  int err;
  MPI_Comm comm = MPI_COMM_WORLD;

  err = MPI_Send(&message, 1, MPI_INT, 0, msgtag, comm); MPI_CHK(err);
  return 0;
}

// TODO: receive a message sent with send_message_to_root() and print it
// can lead printing at the same time), have each process send its message to the
// host to print
int print_message_from_other(void) {
  int err;
  int size = -1;
  int message = -1;
  int rank = -1;
  MPI_Comm comm = MPI_COMM_WORLD;

  err = MPI_Comm_size(comm, &size); MPI_CHK(err);
  for (rank = 1; rank < size; rank++) {
    err = MPI_Recv(&message, 1, MPI_INT, MPI_ANY_SOURCE, msgtag, comm, MPI_STATUS_IGNORE); MPI_CHK(err);
    printf("MPI rank %d of %d, message: %d\n", rank, size, message);
  }
  return 0;
}

// TODO: pass messages around the world communicator in a ring num_trips times
// 
// if rank == 0:
//   message = initial_message
//   for each trip:
//     if verbose:
//       print the message
//     send_message = message + 1
//     send (send_message) to the next process in the ring
//     if verbose:
//       help the other processes print their messages
//     recv (message) from the previous process in the ring
// else:
//   for each trip:
//     recv (message) from the previous process in the ring
//     if verbose:
//       send the message to the root to print
//     send_message = message + 1
//     send (send_message) to the next process in the ring
const int tag = 0;
int ring_test(int num_trips, int verbose, int initial_message, int &final_message)
{
  int err;
  int size = -1;;
  int rank = -1;;
  int message = -1;;
  int send_message = -1;;
  MPI_Comm comm = MPI_COMM_WORLD;

  err = MPI_Comm_size(MPI_COMM_WORLD, &size); MPI_CHK(err);
  err = MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_CHK(err);

  if (rank == 0) {
    message = initial_message;
    for (int trip = 0; trip < num_trips; trip++) {
      if (verbose) {
        print_message(message);
      }
      send_message = message + 1;
      err = MPI_Send(&send_message, 1, MPI_INT, (rank + 1) % size, tag, comm); MPI_CHK(err);
      if (verbose) {
        print_message_from_other();
      }
      err = MPI_Recv(&message, 1, MPI_INT, MPI_ANY_SOURCE, tag, comm, MPI_STATUS_IGNORE); MPI_CHK(err);
    }
  } else {
    for (int trip = 0; trip < num_trips; trip++) {
      err = MPI_Recv(&message, 1, MPI_INT, MPI_ANY_SOURCE, tag, comm, MPI_STATUS_IGNORE); MPI_CHK(err);
      if (verbose) {
        send_message_to_root(message);
      }
      send_message = message + 1;
      err = MPI_Send(&send_message, 1, MPI_INT, (rank + 1) % size, tag, comm); MPI_CHK(err);
    }
  }
  final_message = message;
  return 0;
}

int main(int argc, char **argv)
{
  int      num_trips = 100;
  if (argc > 1) {
    num_trips = atoi(argv[1]);
  }
  bool verbose = (argc > 2);

  int      err;
  err = MPI_Init(&argc, &argv); MPI_CHK(err);


  if (verbose) {
    err = print_mpi_info(); MPI_CHK(err);
  }
  
  MPI_Comm comm = MPI_COMM_WORLD;
  int size;
  err = MPI_Comm_size(comm, &size); MPI_CHK(err);
  if (size == 1) {
    printf("More than one rank is needed for the ring test\n");
    return -1;
  }

  int rank;
  err = MPI_Comm_rank(comm, &rank); MPI_CHK(err);
  int initial_message = 0;
  int final_message = -1;
  int expected_final_message = rank > 0 ? ((num_trips - 1) * size + rank) : num_trips * size;

  double start = MPI_Wtime();
  err = MPI_Barrier(comm); MPI_CHK(err);
  err = ring_test(num_trips, verbose, initial_message, final_message); MPI_CHK(err);
  err = MPI_Barrier(comm); MPI_CHK(err);
  double end = MPI_Wtime();
  double time = end - start;

  int local_success = (final_message == expected_final_message);
  if (!local_success) {
    printf("%d: %d %d\n", rank, final_message, expected_final_message);
  }
  int global_success = 0;
  err = MPI_Reduce(&local_success, &global_success, 1, MPI_INT, MPI_LAND, 0, comm); MPI_CHK(err);
  if (rank == 0) {
    if (!global_success) {
      printf("Incorrect final message state\n");
      err = MPI_Abort(comm, -1);
      return -1;
    }
    printf("{\n \"size\": %d,\n \"number of trips\": %d,\n \"time (seconds)\": %g,\n \"seconds per message\": %g\n}\n",
        size, num_trips, time, time / (num_trips * size));
  }

  err = MPI_Finalize();
  return err;
}
