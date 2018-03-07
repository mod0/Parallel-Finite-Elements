#include <math.h>
#include "vector.h"
#include "error.h"

void vector_init(vector* v, size_t size) {
  v->elements = (double*) calloc(size, sizeof(double));
  v->size = size;
}


void vector_fill(vector* v, double value) {
  int i;
  for (i = 0; i < v->size; i++) {
    v->elements[i] = value;
  }
}


void vector_add(vector* v1, vector* v2, vector* sum) {
  int i;
  for (i = 0; i < sum->size; i++) {
    sum->elements[i] = v1->elements[i] + v2->elements[i];
  }
}


void vector_subtract(vector* v1, vector* v2, vector* difference) {
  int i;
  for (i = 0; i < difference->size; i++) {
    difference->elements[i] = v1->elements[i] - v2->elements[i];
  }
}


void vector_scalar_multiply(double scalar, vector* v, vector* product) {
  int i;
  for (i = 0; i < product->size; i++) {
    product->elements[i] = scalar * v->elements[i];
  }
}


double vector_2_norm(vector* v) {
  double ssq = 0;
  size_t i;
  for (i = 0; i < v->size; i++) {
    ssq += pow(v->elements[i], 2);
  }

  return sqrt(ssq);
}

void vector_print(vector* v) {
  int i;
  printf("[");
  for (i = 0; i < v->size; i++) {
    printf("%f%s", v->elements[i], i == v->size - 1 ? "" : " ");
  }
  printf("]\n");
}

void vector_free(vector* v) {
  free(v->elements);
}
