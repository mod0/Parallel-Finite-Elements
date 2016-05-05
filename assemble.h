#ifndef __ASSEMBLE_H__
#define __ASSEMBLE_H__

#include "matrix.h"
enum bearing{top, bottom, left, right};

void assemble_local_KF(sparse_matrix* K, vector* F, domain* D,
                       int subdomain_idx);

void boundary_op_local_F(vector* F, domain* D, int subdomain_idx);

double bc_choice(domain* D,double u0,
                 int subdomain_idx, int local_coor, int where);

void boundary_op_KF(vector* bands,
                    vector* F, domain* D, int subdomain_idx);

void boundary_op_K(vector* bands,
                   vector* F, domain* D, int subdomain_idx);

#endif
