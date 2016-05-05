#include "grid.h"
#include "domain.h"
#include "vector.h"
#include <string.h>
#include <stdlib.h>

// Convention in ids etc
// Ids for vertices, triangles follow the matlab id convention - starting at 1

// Create one big domain
domain* build_cartesian_domain(grid* cartesian_grid, int subdomain_count_x)
{
  int i;
  domain* cartesian_domain;
  #define CGN (cartesian_grid->N + 1)
  cartesian_domain = (domain*) calloc(1, sizeof(cartesian_domain));
  cartesian_domain->cartesian_grid = cartesian_grid;
  cartesian_domain->subdomain_count_x = subdomain_count_x;
  cartesian_domain->subdomain_vertex_map = (int**) calloc(subdomain_count_x, sizeof(int*));
  #define CDM (cartesian_domain->subdomain_vertex_map)
  for(i = 0; i < subdomain_count_x; i++)
  {
    CDM[i] = (int*) calloc(CGN*CGN, sizeof(int));
    memset(CDM[i], -1, CGN*CGN);
  }
  return cartesian_domain;
}

// Create subdomains and add them to the domain object
int build_subdomains_in_domain(domain* cartesian_domain, int overlap)
{
  int i;

  #define CG (cartesian_domain->cartesian_grid)
  #define CDN (cartesian_domain->cartesian_grid->N + 1)
  #define CSD (cartesian_domain->subdomains)
  #define CSDN (cartesian_domain->subdomain_count_x)

  cartesian_domain->subdomains = (subdomain*) calloc(CSDN, sizeof(subdomain));

  for(i = 0; i < CSDN; i++)
  {
    CSD[i].id = i;                                                                                  // Start with 0 go to p-1
    CSD[i].cartesian_grid = CG;                                                                     // Global grid
    CSD[i].overlap = 2*overlap;                                                                     // Number of overlapping nodes with adjacent subdomains

    if(i == 0)
    {
      CSD[i].bottom_left_x = 0;                                                                     // Bottom left corner in global grid - x
      CSD[i].bottom_left_y = 0;                                                                     // Bottom left corner in global grid - y
      CSD[i].top_right_x = CDN/CSDN + overlap - 1;                                                  // Bottom left corner in global grid - x
      CSD[i].top_right_y = CDN - 1;                                                                 // Bottom left corner in global grid - y
    }
    else if(i == CSDN - 1)
    {
      CSD[i].bottom_left_x = CDN - CDN/CSDN - overlap;                                              // Bottom left corner in global grid - x
      CSD[i].bottom_left_y = 0;                                                                     // Bottom left corner in global grid - y
      CSD[i].top_right_x = CDN - 1;                                                                 // Bottom left corner in global grid - x
      CSD[i].top_right_y = CDN - 1;                                                                 // Bottom left corner in global grid - y
    }
    else
    {
      CSD[i].bottom_left_x = (CDN/CSDN) * i  - overlap;                                             // Bottom left corner in global grid - x
      CSD[i].bottom_left_y = 0;                                                                     // Bottom left corner in global grid - y
      CSD[i].top_right_x = (CDN/CSDN) * (i + 1) + overlap - 1;                                      // Bottom left corner in global grid - x
      CSD[i].top_right_y = CDN - 1;                                                                 // Bottom left corner in global grid - y
    }

    CSD[i].dimX = CSD[i].top_right_x - CSD[i].bottom_left_x + 1;                                    // May include overlap if there are overlapping subdomains
    CSD[i].dimY = CSD[i].top_right_x - CSD[i].bottom_left_x + 1;                                    // May include overlap if there are overlapping subdomains

    CSD[i].subdomain_solution.size = CSD[i].dimX * CSD[i].dimY;                                     // Size of the solution vector
    CSD[i].subdomain_solution.elements = (double*) calloc(CSD[i].dimX*CSD[i].dimY, sizeof(double)); // Solution in the subdomain in the global grid vertex order

    CSD[i].subdomain_vertices = (vertex**) calloc(CSD[i].dimX*CSD[i].dimY, sizeof(vertex*));        // Vertices belonging to the subdomain in the global grid vertex order

    if(i != 0)
    {
      CSD[i].ghost_subdomain_left.size = 2 * overlap * CSD[i].dimY;
      CSD[i].ghost_subdomain_left.elements = (double*) calloc(2*overlap*CSD[i].dimY, sizeof(double));  // Ghost cells into which neighboring threads will write info
    }

    if(i != CSDN - 1)
    {
      CSD[i].ghost_subdomain_right.size = 2 * overlap * CSD[i].dimY;
      CSD[i].ghost_subdomain_right.elements = (double*) calloc(2*overlap*CSD[i].dimY, sizeof(double)); // Ghost cells into which neighboring threads will write info
    }
  }

  #undef CDN
  #undef CSD
  #undef CSDN
  #undef CG

  return 0;
}

// Create vertices in domain
int create_vertices_for_domain(domain* cartesian_domain)
{
  int i, j;

  #define CG (cartesian_domain->cartesian_grid)
  #define CDN (cartesian_domain->cartesian_grid->N + 1)

  cartesian_domain->vertices = (vertex*) calloc(CDN*CDN, sizeof(vertex));

  #define CDV (cartesian_domain->vertices)

  for(i = 0; i < CDN; i++)                             // Go one row after another
  {
    for(j = 0; j < CDN; j++)                           // Loop over columns
    {
      CDV[i].id = i * CDN + j;
      CDV[i].x = CG->grid_nodes_x[j];
      CDV[i].y = CG->grid_nodes_y[i];
    }
  }
  #undef CDN
  #undef CG
  #undef CDV

  return 0;
}

// Go over the list of vertices, add them to the subdomain vertex list.
int create_vertex_subdomain_mapping(domain* cartesian_domain, int idx)
{
  int i, j, count;
  #define CSD (cartesian_domain->subdomains)
  #define CDV (cartesian_domain->vertices)
  #define CDM (cartesian_domain->subdomain_vertex_map)
  #define CDN (cartesian_domain->cartesian_grid->N + 1)

  count = 0;
  for(i = CSD[idx].bottom_left_y; i <= CSD[idx].top_right_y; i++)
  {
    for(j = CSD[idx].bottom_left_x; i <= CSD[idx].top_right_x; j++)
    {
      CSD[idx].subdomain_vertices[count] = &(CDV[i * CDN + j]);
      CDM[idx][i * CDN + j] = count;
      count++;
    }
  }

  #undef CSD
  #undef CDV
  #undef CDN
  #undef CDM
  return 0;
}

// Copy the solution to the right and left subdomains
int copy_ghost_overlap(domain* cartesian_domain, int idx)
{
  int i, j;
  int count;

  #define CSD (cartesian_domain->subdomains)
  #define CSDN (cartesian_domain->subdomain_count_x)

  if(idx > 0)  // copy left
  {
    count = 0;
    for(i = 0; i < CSD[idx].dimY; i++)
    {
      for(j = 0; j < CSD[idx].overlap; j++)
      {
        CSD[idx - 1].ghost_subdomain_right.elements[count++] = CSD[idx].subdomain_solution.elements[i*CSD[idx].dimX + j];
      }
    }
  }

  if(idx < CSDN - 1) // copy right
  {
    count = 0;
    for(i = 0; i < CSD[idx].dimY; i++)
    {
      for(j = CSD[idx].dimX - CSD[idx].overlap; j < CSD[idx].dimX; j++)
      {
        CSD[idx + 1].ghost_subdomain_left.elements[count++] = CSD[idx].subdomain_solution.elements[i*CSD[idx].dimX + j];
      }
    }
  }

  #undef CSD
  #undef CSDN

  return 0;
}
