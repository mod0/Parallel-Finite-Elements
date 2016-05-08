#ifndef __ASSEMBLE_H__
#define __ASSEMBLE_H__

#include "matrix.h"

enum subdomain_wall{top, bottom, left, right};

void assemble_local_KF(sparse_matrix* K, vector* F, domain* D,
                       int subdomain_idx);

void boundary_op_local_F(vector* F, domain* D, int subdomain_idx);

double get_boundary_value(domain* D, int subdomain_idx, int vertex_id, int wall);

void boundary_op_K(vector* bands,
                   vector* F, domain* D, int subdomain_idx);



#endif
