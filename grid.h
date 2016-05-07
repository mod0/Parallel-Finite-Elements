#ifndef __GRID_H__
#define __GRID_H__

// Describes the cartesian grid used in the problem
typedef struct
{
  int N;                            // Number of segments along each direction
  double lb_x;
  double ub_x;
  double lb_y;
  double ub_y;
  double h;                         // Max grid resolution
  double* grid_nodes_x;             // Partition in x direction
  double* grid_nodes_y;             // Partition in y direction
} grid;

// Partitions the X/Y axis into partitions pieces
double* uniform_partition(double lb, double ub, int N);

// Builds a cartesian grid with the given grid properties
grid* build_cartesian_grid(double lb_x, double ub_x, double lb_y, double ub_y, int N);

// Cleanup the grid
int cleanup_grid(grid* cartesian_grid);


#endif
