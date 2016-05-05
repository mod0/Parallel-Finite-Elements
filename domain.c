#include <stdlib.h>

#include "grid.h"
#include "domain.h"
#include "vector.h"

// Convention in ids etc
// Ids for vertices, triangles follow the matlab id convention - starting at 1

// Create one big domain
domain* build_cartesian_domain(grid* cartesian_grid, int subdomain_count_x)
{
  int i,j;
  domain* cartesian_domain;
  #define CGN (cartesian_grid->N + 1)
  cartesian_domain = calloc(1, sizeof(domain));
  cartesian_domain->cartesian_grid = cartesian_grid;
  cartesian_domain->subdomain_count_x = subdomain_count_x;
  cartesian_domain->subdomain_vertex_map = calloc(2, sizeof(int*));
  #define CDM (cartesian_domain->subdomain_vertex_map)
  for(i = 0; i < subdomain_count_x; i++)
  {
    cartesian_domain->subdomain_vertex_map[i] = calloc(CGN*CGN, sizeof(int));
    for(j = 0; j < CGN*CGN; j++)
    {
      CDM[i][j] = -1;
    }
  }
  #undef CGN
  #undef CDM

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
    CSD[i].dimY = CSD[i].top_right_y - CSD[i].bottom_left_y + 1;                                    // May include overlap if there are overlapping subdomains

    CSD[i].subdomain_vertices = calloc(CSD[i].dimX*CSD[i].dimY, sizeof(vertex*));                   // Vertices belonging to the subdomain in the global grid vertex order

    vector_init(&(CSD[i].subdomain_solution), CSD[i].dimX*CSD[i].dimY);                             // Solution in the subdomain in the global grid vertex order
    if(i != 0)
    {
      vector_init(&(CSD[i].ghost_subdomain_left), 2*overlap*CSD[i].dimY);                            // Ghost cells into which neighboring threads will write info
    }

    if(i != CSDN - 1)
    {
      vector_init(&(CSD[i].ghost_subdomain_right), 2*overlap*CSD[i].dimY);                            // Ghost cells into which neighboring threads will write info
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
  int i, j, count;

  #define CG (cartesian_domain->cartesian_grid)
  #define CDN (cartesian_domain->cartesian_grid->N + 1)

  cartesian_domain->vertices = (vertex*) calloc(CDN*CDN, sizeof(vertex));

  #define CDV (cartesian_domain->vertices)

  count = 0;
  for(i = 0; i < CDN; i++)                             // Go one row after another
  {
    for(j = 0; j < CDN; j++)                           // Loop over columns
    {
      CDV[count].id = i * CDN + j;
      CDV[count].x = CG->grid_nodes_x[j];
      CDV[count].y = CG->grid_nodes_y[i];
      count++;
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
    for(j = CSD[idx].bottom_left_x; j <= CSD[idx].top_right_x; j++)
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

// Copy the solution to the right, left, or both subdomains
int copy_ghost_overlap(domain* cartesian_domain, int idx, int direction)
{
  int i, j;
  int count;

  #define CSD (cartesian_domain->subdomains)
  #define CSDN (cartesian_domain->subdomain_count_x)

  if(direction <= 0)
  {
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
  }

  if(direction >= 0)
  {
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
  }

  #undef CSD
  #undef CSDN

  return 0;
}
