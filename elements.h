#ifndef __ELEMENTS_H__
#define __ELEMENTS_H__

#include "domain.h"
#include "grid.h"


// Define triangular_element
// TODO: Make these generic elements
typedef struct
{
  int id;
  domain* cartesian_domain;
  const int element_vertex_count = 3;
  vertex* grid_vertex[3];
} triangular_element;

// Count of number of elements;
int _total_grid_elements;

// Create a global array of elements
triangular_element* _triangular_elements;

// Creates elements and stores it in the triangular elements list
int create_triangular_elements_for_cartesian_domain(domain* cartesian_domain);

// Add element indices to subdomains
int add_triangular_elements_to_subdomains(domain* cartesian_domain, int idx);

// Free memory allocated on heap
int cleanup_triangular_elements();

#endif
