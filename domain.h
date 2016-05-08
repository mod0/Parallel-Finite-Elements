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
  int elements_count;              // The number of elements in the subdomain
  int* elements;                  // List of Elements Indices
  int converged;                  // Boolean flag to check whether the subdomain has converged
} subdomain;

// Domain object
typedef struct domain
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
int copy_overlap_to_adjacent_neighbours_ghost(domain* cartesian_domain, int idx, int direction);

// Copy the solution from the left or right ghost cell of the same subdomain
int copy_from_my_ghost_cell(domain* cartesian_domain, int idx, int direction);

// Boolean method returns true or false to write output for vertex. Each subdomain only writes the right overlap.
int write_output_for_vertex(domain* d, int subdomainIdx, vertex* v);

// Remove all vertex objects from heap
int cleanup_vertices(domain* cartesian_domain);

// Remove all subdomain objects from heap
int cleanup_subdomains(domain* cartesian_domain);

// Remove domain object from heap
int cleanup_domain(domain* cartesian_domain);

#endif
