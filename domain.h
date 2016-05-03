#ifndef __DOMAIN_H__
#define __DOMAIN_H__

#include "grid.h"
#include "elements.h"

// Information about domains and subdomains

// Grid vertices numbered in x direction
// 133 134 135 136 137 138 139 140 141 142 143 144
// 121 122 123 124 125 126 127 128 129 130 131 132
// 109 110 111 112 113 114 115 116 117 118 119 120
//  97  98  99 100 101 102 103 104 105 106 107 108
//  85  86  87  88  89  90  91  92  93  94  95  96
//  73  74  75  76  77  78  79  80  81  82  83  84
//  61  62  63  64  65  66  67  68  69  70  71  72
//  49  50  51  52  53  54  55  56  57  58  59  60
//  37  38  39  40  41  42  43  44  45  46  47  48
//  25  26  27  28  29  30  31  32  33  34  35  36
//  13  14  15  16  17  18  19  20  21  22  23  24
//   1   2   3   4   5   6   7   8   9  10  11  12

// Direction specifier for numbering vertices
typedef enum {dirx, diry} vertex_numbering_dir;

// Vertex object
typedef struct
{
  int global_vertex_id;     // The id based on the grid and direction
  double x;                 // The x location
  double y;                 // The y location
} vertex;

// Direction specifier for grid partition
typedef enum {dirx, diry, dirxy} grid_partition_direction;

// Subdomain object
typedef struct
{
  grid* region;
  int dimX;             // May include overlap if there are overlapping subdomains
  int dimY;             // May include overlap if there are overlapping subdomains
  int overlap;          // Number of overlapping nodes with adjacent subdomains
  int bottom_left_x;    // Bottom left corner in global grid - x
  int bottom_left_y;    // Bottom left corner in global grid - y
  subdomain* left;      // Pointer to left subdomain
  subdomain* right;     // Pointer to right subdomain
  subdomain* top;       // Pointer to top subdomain
  subdomain* bottom;    // Pointer to bottom subdomain
  double* ghost_subdomain_left;   // Ghost cells into which neighboring threads will write info
  double* ghost_subdomain_right;  // Ghost cells into which neighboring threads will write info
  double* ghost_subdomain_top;    // Ghost cells into which neighboring threads will write info
  double* ghost_subdomain_bottom; // Ghost cells into which neighboring threads will write info
  double* subdomain_solution;     // Solution in the subdomain in the global grid vertex order
  vertex* subdomain_vertices;     // Vertices belonging to the subdomain in the global grid vertex order
} subdomain;

// Domain object
typedef struct
{
  grid* region;
  int subdomain_count_x;
  int subdomain_count_y;
  subdomain* subdomains;
  vertex_numbering_dir direction;   // Direction in which grid vertices be numbered
} domain;

#endif
