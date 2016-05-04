#ifndef __VECTOR_H__
#define __VECTOR_H__

#include <stdlib.h>

typedef struct
{
    double* elements;
    size_t size;
} vector;

void vector_init(vector* v, size_t size);

void vector_fill(vector* v, double value);

void vector_add(vector* v1, vector* v2, vector* sum);

void vector_subtract(vector* v1, vector* v2, vector* difference);

void vector_scalar_multiply(double scalar, vector* v, vector* product);

double vector_2_norm(vector* v);

void vector_free(vector* v);

#endif
