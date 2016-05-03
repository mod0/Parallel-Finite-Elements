#ifndef __GRID_H__
#define __GRID_H__

// Properties of cartesian grid - extent and resolution
typedef struct
{
  double lb_x;
  double ub_x;
  double lb_y;
  double ub_y;
  int N;
} grid_prop;

// Describes the cartesian grid used in the problem
typedef struct
{
  grid_prop gProperties;            // Grid property object
  double* grid_nodes_x;             // Partition in x direction
  double* grid_nodes_y;             // Partition in y direction
} grid;

// Partitions the X/Y axis into partitions pieces
double* uniform_partition(double lb, double ub, int N);

// Builds a cartesian grid with the given grid properties
grid* build_cartesian_grid(grid_prop gProperties);

#endif
