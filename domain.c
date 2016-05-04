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
  subdomain* cartesian_subdomain;

  subdomain_count_x = cartesian_domain->subdomain_count_x;

  cartesian_domain->subdomains = (subdomain**) calloc(subdomain_count_x, sizeof(subdomain*));

  #define CDN (cartesian_domain->N)
  #define PROC (subdomain_count_x)

  for(i = 0; i < subdomain_count_x; i++)
  {
    cartesian_subdomain = (subdomain*) calloc(1, sizeof(subdomain));
    cartesian_subdomain->id = i + 1;                                                                                      // Start with 0 go to p-1
    cartesian_subdomain->cartesian_grid = cartesian_domain->cartesian_grid;                                               // Global grid
    cartesian_subdomain->overlap = overlap;                                                                               // Number of overlapping nodes with adjacent subdomains

    if(i == 0)
    {
      cartesian_subdomain->bottom_left_x = 1;                                                                             // Bottom left corner in global grid - x
      cartesian_subdomain->bottom_left_y = 1;                                                                             // Bottom left corner in global grid - y
      cartesian_subdomain->top_right_x = CDN/PROC + overlap;                                                              // Bottom left corner in global grid - x
      cartesian_subdomain->top_right_y = N;                                                                               // Bottom left corner in global grid - y
    }
    else if(i == subdomain_count_x - 1)
    {
      cartesian_subdomain->bottom_left_x = CDN - CDN/PROC + 1 - overlap;                                                  // Bottom left corner in global grid - x
      cartesian_subdomain->bottom_left_y = 1;                                                                             // Bottom left corner in global grid - y
      cartesian_subdomain->top_right_x = N;                                                                               // Bottom left corner in global grid - x
      cartesian_subdomain->top_right_y = N;                                                                               // Bottom left corner in global grid - y
      cartesian_subdomain->left = cartesian_domain->subdomains[i - 1];                                                    // Store the left subdomain
      cartesian_domain->subdomains[i-1]->right = cartesian_subdomain;                                                     // Store the right subdomain
    }
    else
    {
      cartesian_subdomain->bottom_left_x = (CDN/PROC) * i  + 1 - overlap;                                                 // Bottom left corner in global grid - x
      cartesian_subdomain->bottom_left_y = 1;                                                                             // Bottom left corner in global grid - y
      cartesian_subdomain->top_right_x = (CDN/PROC) * (i + 1) + overlap;                                                  // Bottom left corner in global grid - x
      cartesian_subdomain->top_right_y = N;                                                                               // Bottom left corner in global grid - y
      cartesian_subdomain->left = cartesian_domain->subdomains[i - 1];                                                    // Store the left subdomain
      cartesian_domain->subdomains[i-1]->right = cartesian_subdomain;                                                     // Store the right subdomain
    }

    cartesian_domain->subdomains[i] = cartesian_subdomain;

    // cartesian_subdomain->dimX;                                                // May include overlap if there are overlapping subdomains
    // cartesian_subdomain->dimY;                                                // May include overlap if there are overlapping subdomains
    // cartesian_subdomain->ghost_subdomain_left;                                // Ghost cells into which neighboring threads will write info
    // cartesian_subdomain->ghost_subdomain_right;                               // Ghost cells into which neighboring threads will write info
    // cartesian_subdomain->subdomain_solution;                                  // Solution in the subdomain in the global grid vertex order
    // cartesian_subdomain->subdomain_vertices;                                  // Vertices belonging to the subdomain in the global grid vertex order
  }

  #undef CDN
  #undef PROC
}
