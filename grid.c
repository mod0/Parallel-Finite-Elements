#include <stdlib.h>
#include "grid.h"
#include "error.h"

double* uniform_partition(double lb, double ub, int N)
{
  int i;
  double h;
  double* partition;

  // Compute the step size
  h = (ub - lb)/N;

  // Allocate space for the parition
  partition = (double*) calloc(N+1, sizeof(double));
  for(i = 1, partition[0] = lb; i < N; i++)
  {
    partition[i] = partition[i - 1] + h;
  }
  partition[N] = ub;

  return partition;
}


grid* build_cartesian_grid(double lb_x, double ub_x, double lb_y, double ub_y, int N)
{
  grid* cartesian_grid;
  if(ub_x - lb_x != ub_y - lb_y)
  {
    error("The domain is not square.");
  }

  cartesian_grid = (grid*) malloc(sizeof(grid));
  cartesian_grid->N = N;
  cartesian_grid->lb_x = lb_x;
  cartesian_grid->ub_x = ub_x;
  cartesian_grid->lb_y = lb_y;
  cartesian_grid->ub_y = ub_y;
  cartesian_grid->h = (ub_x - lb_x)/N;
  cartesian_grid->grid_nodes_x = uniform_partition(lb_x, ub_x, N);  // Partition in x direction
  cartesian_grid->grid_nodes_y = uniform_partition(lb_y, ub_y, N);  // Partition in y direction
  return cartesian_grid;
}
