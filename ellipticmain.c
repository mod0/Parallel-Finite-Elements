//
//  Finite Element Analysis of Poisson and Heat Equation
//  on 2D cartesian domain using triangular elements.
//

#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#include "grid.h"
#include "domain.h"
#include "elements.h"
#include "assemble.h"
#include "parameters.h"
#include "solver.h"
#include "output.h"

int main(int argc, char** argv)
{
  int i;

  // TODO: Either read the following data from the file or ask the user
  double lb_x = 0;    // The lower bound of the domain in the x direction
  double ub_x = 1;    // The upper bound of the domain in the y direction
  double lb_y = 0;    // The lower bound of the domain in the y direction`
  double ub_y = 1;    // The upper bound of the domain in the y direction
  int N = 99;        // The number of grid segments in each direction
  int subdomains = 2; // The number of threads in the system
  int overlap_in_each_direction = 4; // Amount of overlap in each direction -> 2 times will be the amount of overlap

  if (argc > 1) {
      subdomains = atoi(argv[1]);
  }
  if (argc > 2) {
      N = atoi(argv[2]);
  }
  if (argc > 3) {
      overlap_in_each_direction = atoi(argv[3]);
  }



  #ifdef _OPENMP
    omp_set_num_threads(subdomains);
  #endif

  mgmres_parameters linear_solve_parameters = {.outerItr = 2, .innerItr = N, .absTol = 1e-8, .relTol = 1e-8};
  elliptic_solver_parameters solver_parameters = {.mgmresParameters = linear_solve_parameters, .outputProcessor = file_output_processor,
                                                  .solverRelTol = 1e-3, .maxItr = 99999};

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
  #pragma omp parallel for private(i)
  for(i = 0; i < subdomains; i++)
  {
    create_vertex_subdomain_mapping(cartesian_domain, i);
  }

  // Add elements to subdomains
  #pragma omp parallel for private(i)
  for(i = 0; i < subdomains; i++)
  {
    add_triangular_elements_to_subdomains(cartesian_domain, i);
  }

  // Call the solver
  ellipticsolver(cartesian_domain,  solver_parameters);

  // Clean up
  cleanup_triangular_elements();
  cleanup_subdomains(cartesian_domain);
  cleanup_domain(cartesian_domain);
  cleanup_grid(cartesian_grid);

  return 0;
}
