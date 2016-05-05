#ifndef __MATRIX_H__
#define __MATRIX_H__

#include "vector.h"

typedef struct
{
    int size;
    vector elements;
    int* rows;
    int* cols;
} sparse_matrix;

void sparse_matrix_init(sparse_matrix* m, size_t size, size_t numNonzero);

void sparse_symmetric_banded_init(sparse_matrix* m, size_t size, vector* bands, size_t numBands);

double sparse_matrix_get(sparse_matrix* m, int row, int col);

void sparse_matrix_vector_multiply(sparse_matrix* m, vector* x, vector* w);

double sparse_matrix_frobenius_norm(sparse_matrix* m);

void mgmres(sparse_matrix* matrix, vector* x, vector* rhs, int itr_max, int mr, double tol_abs, double tol_rel);

void sparse_matrix_print(sparse_matrix* m, size_t ellipsisThreshold);

void sparse_matrix_free(sparse_matrix* m);

#endif
