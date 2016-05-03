#include "grid.h"
#include "domain.h"
#include <stdlib.h>

// Create one big domain
domain* build_cartesian_domain(grid* cartesian_grid, int subdomain_count_x)
{
  domain* cartesian_domain;

  cartesian_domain = (domain*) calloc(1, sizeof(cartesian_domain));
  cartesian_domain->cartesian_grid = cartesian_grid;
  cartesian_domain->subdomain_count_x = subdomain_count_x;
  cartesian_domain->vertex_numbering_dir = vertex_numbering_dir;

  return cartesian_domain;
}


// Create subdomains and add them to the domain object
int build_subdomains_in_domain(domain* cartesian_domain, int overlap)
{
  int i;
  int subdomain_count_x;
  subdomain** all_subdomains;
  subdomain* cartesian_subdomain;

  subdomain_count_x = cartesian_domain->subdomain_count_x;

  all_subdomains = (subdomain**) calloc(subdomain_count_x, sizeof(subdomain*));

  for(i = 0; i < subdomain_count_x; i++)
  {
    cartesian_subdomain = (subdomain*) calloc(1, sizeof(subdomain));
    cartesian_subdomain->cartesian_grid = cartesian_domain->cartesian_grid;   // Global grid
    cartesian_subdomain->dimX;                                                // May include overlap if there are overlapping subdomains
    cartesian_subdomain->dimY;                                                // May include overlap if there are overlapping subdomains
    cartesian_subdomain->overlap;                                             // Number of overlapping nodes with adjacent subdomains
    cartesian_subdomain->bottom_left_x;                                       // Bottom left corner in global grid - x
    cartesian_subdomain->bottom_left_y;                                       // Bottom left corner in global grid - y
    cartesian_subdomain->top_right_x;                                         // Bottom left corner in global grid - x
    cartesian_subdomain->top_right_y;                                         // Bottom left corner in global grid - y
    cartesian_subdomain->left;                                                // Pointer to left subdomain
    cartesian_subdomain->right;                                               // Pointer to right subdomain
    cartesian_subdomain->ghost_subdomain_left;                                // Ghost cells into which neighboring threads will write info
    cartesian_subdomain->ghost_subdomain_right;                               // Ghost cells into which neighboring threads will write info
    cartesian_subdomain->subdomain_solution;                                  // Solution in the subdomain in the global grid vertex order
    cartesian_subdomain->subdomain_vertices;                                  // Vertices belonging to the subdomain in the global grid vertex order
  }
}
