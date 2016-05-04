#ifndef __ELEMENTS_H__
#define __ELEMENTS_H__

#include "domain.h"
#include "grid.h"


// Define triangular_element
typedef struct
{
  int id;
  grid* cartesian_grid;
  vertex* grid_vertex[3];
} triangular_element;

// Create a global array of elements
triangular_element* _triangular_elements;

int create_triangular_elements_in_cartesian_domain();

#endif
