#ifndef __ELEMENTS_H__
#define __ELEMENTS_H__

#include "domain.h"
#include "grid.h"


// Define triangular_element
typedef struct
{
  int id;
  domain* cartesian_domain;
  vertex* grid_vertex[3];
} triangular_element;

// Create a global array of elements
triangular_element* _triangular_elements;

// Creates elements and stores it in the triangular elements list
int create_triangular_elements_for_cartesian_domain(domain* cartesian_domain);

// Add element indices to subdomains
int add_triangular_elements_to_subdomains(domain* cartesian_domain, int idx);

#endif
