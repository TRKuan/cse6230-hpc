#include "params.h"
#include "exact.h"
#include "grid.h"
#include "solver.h"

int main(int argc, char **argv) {
  int err = MPI_Init(&argc, &argv); MPI_ECHK(err);
  {
    Json outlog = Json();

    Params params(argc, argv);
    params.write(outlog);

    GridFunction F = create_rhs(params);
    GridFunction Xexact = create_sol(params, F);

    double predicted_cf = convergence_factor(params);
    outlog.write("Predicted convergence factor", predicted_cf);

    Solver solver(params, outlog, F, Xexact);

    GridFunction X = solver.solve();
  }
  err = MPI_Finalize();
  return err;
}
