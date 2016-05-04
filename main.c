//
//  Finite Element Analysis of Poisson and Heat Equation
//  on 2D cartesian domain using triangular elements.
//

#include <time.h>
#include <stdlib.h>
#include <stdio.h>

#include "grid.h"
#include "domain.h"

int main(int argc, char** argv)
{
  // TODO: Either read the following data from the file or ask the user
  double lb_x = 0;    // The lower bound of the domain in the x direction
  double ub_x = 1;    // The upper bound of the domain in the y direction
  double lb_y = 0;    // The lower bound of the domain in the y direction`
  double ub_y = 1;    // The upper bound of the domain in the y direction
  int N = 100;        // The number of grid segments
  int overlap = 1;    // Amount of overlap in each direction
  // Create a grid with the grid properties object.
  grid* cartesian_grid = build_cartesian_grid(lb_x, ub_x, lb_y, ub_y, N);

  // Create a domain
  // domain* cartesian_domain = build_cartesian_domain(cartesian_grid, overlap);

  return 0;
}
