#ifndef __DOMAIN_H__
#define __DOMAIN_H__

#include "grid.h"
#include "vector.h"

// Information about domains and subdomains

// Vertex object
typedef struct
{
  int id;     // The id based on the grid and direction
  double x;                 // The x location
  double y;                 // The y location
} vertex;

// Subdomain object
typedef struct subdomain
{
  int id;  // Index of the subdomain along the x direction
  grid* cartesian_grid; // Global grid
  int dimX;             // May include overlap if there are overlapping subdomains
  int dimY;             // May include overlap if there are overlapping subdomains
  int overlap;          // Number of overlapping nodes with adjacent subdomains
  int bottom_left_x;    // Bottom left corner in global grid - x
  int bottom_left_y;    // Bottom left corner in global grid - y
  int top_right_x;      // Top right corner in global grid - x
  int top_right_y;      // Top right corner in global grid - y
  vector ghost_subdomain_left;   // Ghost cells into which neighboring threads will write info
  vector ghost_subdomain_right;  // Ghost cells into which neighboring threads will write info
  vector subdomain_solution;     // Solution in the subdomain in the global grid vertex order
  vertex** subdomain_vertices;    // Vertices belonging to the subdomain in the global grid vertex order
  int* elements;                  // List of Elements Indices
} subdomain;

// Domain object
typedef struct
{
  grid* cartesian_grid;              // Reference to the grid object
  vertex* vertices;                  // Array of all vertexes in the domain
  int subdomain_count_x;             // Subdomains in X direction
  subdomain* subdomains;             // List of all subdomains
  int** subdomain_vertex_map;        // Map of subdomain to vertices, gives the local vertex order for assembly
} domain;

// Create a domain with as many number of
domain* build_cartesian_domain(grid* cartesian_grid, int subdomain_count_x);

// Create subdomains and add them to the domain object
int build_subdomains_in_domain(domain* cartesian_domain, int overlap);

// Create vertices in domain
int create_vertices_for_domain(domain* cartesian_domain);

// Create vertex subdomain mapping
int create_vertex_subdomain_mapping(domain* cartesian_domain, int subdomain_index);

// Copy overlapping solution to adjacent subdomains
int copy_ghost_overlap(domain* cartesian_domain, int idx);

#endif
