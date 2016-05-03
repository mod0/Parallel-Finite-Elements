#include <stdlib.h>
#include "grid.h"

double* uniform_partition(double lb, double ub, int N)
{
  int i;
  double h;
  double location;
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


grid* build_cartesian_grid(grid_prop gProperties)
{
  grid* cartesian_grid;
  cartesian_grid = (grid*) malloc(sizeof(grid));
  cartesian_grid->gProperties = gProperties;        // Grid property object
  cartesian_grid->grid_nodes_x = uniform_partition(gProperties.lb_x, gProperties.ub_x, gProperties.N);  // Partition in x direction
  cartesian_grid->grid_nodes_y = uniform_partition(gProperties.lb_y, gProperties.ub_y, gProperties.N);  // Partition in y direction
  return cartesian_grid;
}
