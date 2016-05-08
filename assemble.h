#ifndef __ASSEMBLE_H__
#define __ASSEMBLE_H__

#include "matrix.h"

enum subdomain_wall{top, bottom, left, right};

double get_boundary_value(domain* D, int subdomain_idx, int vertex_id, int wall);

void assemble_global_load_vector(domain* D, int subdomain_idx, vector* global_load_vector);

void apply_boundary_operator_on_vector(domain* D, int subdomain_idx, vector* F);

void assemble_global_matrix(domain* D, int subdomain_idx, double local_matrix[][3], sparse_matrix* global_matrix,
                           int vector_sizes[], int diagonal_offsets[], int diag_count, double boundary_diagonal_value);

void apply_boundary_operator_on_matrix(domain* D,  int subdomain_idx, vector* bands, int vector_sizes[],
                                      int diagonal_offsets[], int diag_count, double boundary_diagonal_value);

double triangular_element_one_point_quadrature(domain* cartesian_domain, int subdomain_idx, int element_id);

#endif
