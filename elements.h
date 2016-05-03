#ifndef __ELEMENTS_H__
#define __ELEMENTS_H__

#include "grid.h"
#include "domain.h"

triangular_element* elements;

typedef struct
{
  int id;
  grid* region;
  vertex* grid_vertex_1;
  vertex* grid_vertex_2;
  vertex* grid_vertex_3;
} triangular_element;

#endif
