#include <stdlib.h>

#include "elements.h"

// Creates elements and stores it in the triangular elements list
int create_triangular_elements_for_cartesian_domain(domain* cartesian_domain)
{
  int i, j;
  int nT, sqid, tid;
  vertex* vertex1;
  vertex* vertex2;
  vertex* vertex3;
  vertex* vertex4;

  #define CDN (cartesian_domain->cartesian_grid->N + 1)
  #define CDV (cartesian_domain->vertices)

  nT = 2 * (CDN - 1) * (CDN - 1);
  _total_grid_elements = 0;
  _triangular_elements = calloc(nT, sizeof(triangular_element));

  for(i = 0; i < CDN - 1; i++)
  {
    for(j = 0; j < CDN - 1; j++)
    {
      sqid = i * (CDN - 1) + j;
      tid = 2 * sqid;

      vertex1 = &CDV[i * CDN + j];
      vertex2 = &CDV[i * CDN + j + 1];
      vertex3 = &CDV[(i + 1) * CDN + j];
      vertex4 = &CDV[(i + 1) * CDN + j + 1];

      _triangular_elements[tid].id = tid;
      _triangular_elements[tid].cartesian_domain = cartesian_domain;
      _triangular_elements[tid].grid_vertex[0] = vertex3;
      _triangular_elements[tid].grid_vertex[1] = vertex1;
      _triangular_elements[tid].grid_vertex[2] = vertex4;
      _total_grid_elements++;

      _triangular_elements[tid + 1].id = tid + 1;
      _triangular_elements[tid + 1].cartesian_domain = cartesian_domain;
      _triangular_elements[tid + 1].grid_vertex[0] = vertex2;
      _triangular_elements[tid + 1].grid_vertex[1] = vertex4;
      _triangular_elements[tid + 1].grid_vertex[2] = vertex1;
      _total_grid_elements++;
    }
  }

  #undef CDN
  #undef CDV

  return 0;
}

// Add element indices to subdomains
int add_triangular_elements_to_subdomains(domain* cartesian_domain, int idx)
{
  int i, j, count;
  int nT, sqid, tid;

  #define CSD (cartesian_domain->subdomains)
  #define CDN (cartesian_domain->cartesian_grid->N + 1)

  nT = 2 * (CSD[idx].dimX - 1) * (CSD[idx].dimY - 1);
  CSD[idx].elements = calloc(nT, sizeof(int));

  count = 0;
  for(i = CSD[idx].bottom_left_y; i < CSD[idx].top_right_y; i++)
  {
    for(j = CSD[idx].bottom_left_x; j < CSD[idx].top_right_x; j++)
    {
      sqid = i * (CDN - 1) + j;
      tid = 2 * sqid;
      CSD[idx].elements[count++] = tid++;
      CSD[idx].elements[count++] = tid;
    }
  }

  #undef CSD
  #undef CDN

  return 0;
}


// Remove all triangular elements from heap
int cleanup_triangular_elements()
{
  if(_triangular_elements != NULL)
  {
    free(_triangular_elements);
  }

  return 0;
}
