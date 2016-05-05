//
//  Finite Element Analysis of Poisson and Heat Equation
//  on 2D cartesian domain using triangular elements.
//

#include <time.h>
#include <stdlib.h>
#include <stdio.h>

#include "grid.h"
#include "domain.h"
#include "elements.h"
#include "assemble.h"
#include "fep.h"

int main(int argc, char** argv)
{
  int i;

  // TODO: Either read the following data from the file or ask the user
  double lb_x = 0;    // The lower bound of the domain in the x direction
  double ub_x = 1;    // The upper bound of the domain in the y direction
  double lb_y = 0;    // The lower bound of the domain in the y direction`
  double ub_y = 1;    // The upper bound of the domain in the y direction
  int N = 100;        // The number of grid segments in each direction
  int subdomains = 1; // The number of threads in the system
  int overlap_in_each_direction = 0; // Amount of overlap in each direction -> 2 times will be the amount of overlap

  // Create a grid with the grid properties
  grid* cartesian_grid = build_cartesian_grid(lb_x, ub_x, lb_y, ub_y, N);

  // Create a domain
  domain* cartesian_domain = build_cartesian_domain(cartesian_grid, subdomains);

  // Create subdomains with appropriate overlaps
  build_subdomains_in_domain(cartesian_domain, overlap_in_each_direction);

  // Create all vertices in the domain
  create_vertices_for_domain(cartesian_domain);

  // Create all triangular elements in the domain
  create_triangular_elements_for_cartesian_domain(cartesian_domain);

  // Add vertices to subdomains and create a vertex subdomain mapping
  for(i = 0; i < subdomains; i++)
  {
    create_vertex_subdomain_mapping(cartesian_domain, i);
  }

  // Add elements to subdomains
  for(i = 0; i < subdomains; i++)
  {
    add_triangular_elements_to_subdomains(cartesian_domain, i);
  }




  return 0;
}
