#ifndef __ELEMENTS_H__
#define __ELEMENTS_H__

#include "grid.h"
#include "domain.h"

triangular_element* _triangular_elements;

typedef struct
{
  int id;
  grid* cartesian_grid;
  vertex* grid_vertex_1;
  vertex* grid_vertex_2;
  vertex* grid_vertex_3;
} triangular_element;

int create_triangular_elements_in_cartesian_domain()

#endif
