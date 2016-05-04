#include "grid.h"
#include "domain.h"
#include <stdlib.h>

// Convention in ids etc
// Ids for vertices, triangles follow the matlab id convention - starting at 1

// Create one big domain
domain* build_cartesian_domain(grid* cartesian_grid, int subdomain_count_x)
{
  domain* cartesian_domain;

  cartesian_domain = (domain*) calloc(1, sizeof(cartesian_domain));
  cartesian_domain->cartesian_grid = cartesian_grid;
  cartesian_domain->subdomain_count_x = subdomain_count_x;
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

  #define CDN (cartesian_domain->N + 1)
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
      cartesian_subdomain->top_right_y = CDN;                                                                               // Bottom left corner in global grid - y
    }
    else if(i == subdomain_count_x - 1)
    {
      cartesian_subdomain->bottom_left_x = CDN - CDN/PROC + 1 - overlap;                                                  // Bottom left corner in global grid - x
      cartesian_subdomain->bottom_left_y = 1;                                                                             // Bottom left corner in global grid - y
      cartesian_subdomain->top_right_x = CDN;                                                                               // Bottom left corner in global grid - x
      cartesian_subdomain->top_right_y = CDN;                                                                               // Bottom left corner in global grid - y
      cartesian_subdomain->left = cartesian_domain->subdomains[i - 1];                                                    // Store the left subdomain
      cartesian_domain->subdomains[i-1]->right = cartesian_subdomain;                                                     // Store the right subdomain
    }
    else
    {
      cartesian_subdomain->bottom_left_x = (CDN/PROC) * i  + 1 - overlap;                                                 // Bottom left corner in global grid - x
      cartesian_subdomain->bottom_left_y = 1;                                                                             // Bottom left corner in global grid - y
      cartesian_subdomain->top_right_x = (CDN/PROC) * (i + 1) + overlap;                                                  // Bottom left corner in global grid - x
      cartesian_subdomain->top_right_y = CDN;                                                                               // Bottom left corner in global grid - y
      cartesian_subdomain->left = cartesian_domain->subdomains[i - 1];                                                    // Store the left subdomain
      cartesian_domain->subdomains[i-1]->right = cartesian_subdomain;                                                     // Store the right subdomain
    }

    cartesian_domain->subdomains[i] = cartesian_subdomain;
  }

  #undef CDN
  #undef PROC

  return 0;
}

// Create vertices in domain
int create_vertices_for_domain(domain* cartesian_domain)
{
  int i, j;
  vertex* grid_vertex;

  #define CG (cartesian_domain->cartesian_grid)
  #define CDN (cartesian_domain->cartesian_grid->N + 1)

  cartesian_domain->vertices = (vertices**) calloc(CDN*CDN, sizeof(vertex*));

  for(i = 1; i <= CDN; i++)                             // Go one row after another
  {
    for(j = 1; j <= CDN; j++)                           // Loop over columns
    {
      grid_vertex = (vertex*) calloc(1, sizeof(vertex));
      grid_vertex->global_vertex_id = (i - 1)*N + j;
      grid_vertex->x = CG->grid_nodes_x[j-1];
      grid_vertex->y = CG->grid_nodes_y[i-1];
      cartesian_domain->vertices[(i-1)*CDN + (j-1)] = grid_vertex;
    }
  }
  #undef CDN
  #undef CG
}

// Go over the list of vertices, add them to the subdomain vertex list.
int create_vertex_subdomain_mapping(domain* cartesian_domain, int subdomain_index)
{
  subdomain* cartesian_subdomain;

  cartesian_subdomain = cartesian_domain->subdomains[subdomain_index];
  cartesian_subdomain->dimX = cartesian_subdomain->top_right_x - cartesian_subdomain->bottom_left_x + 1;                                                // May include overlap if there are overlapping subdomains
  cartesian_subdomain->dimY = cartesian_subdomain->top_right_x - cartesian_subdomain->bottom_left_x + 1;                                                // May include overlap if there are overlapping subdomains
  cartesian_subdomain->ghost_subdomain_left;                                // Ghost cells into which neighboring threads will write info
  cartesian_subdomain->ghost_subdomain_right;                               // Ghost cells into which neighboring threads will write info
  cartesian_subdomain->subdomain_solution;                                  // Solution in the subdomain in the global grid vertex order
  cartesian_subdomain->subdomain_vertices;                                  // Vertices belonging to the subdomain in the global grid vertex order

}
