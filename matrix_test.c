#include <stdio.h>
#include <math.h>
#include "vector.h"
#include "matrix.h"
#include "parameters.h"

static void test_mgmres() {
# define N 20
# define NZ_NUM 3 * N - 2

  double a[NZ_NUM];
  int i;
  int ia[NZ_NUM];
  int j;
  int ja[NZ_NUM];
  int k;
  int n = N;
  int nz_num = NZ_NUM;
  vector rhs; vector_init(&rhs, N);
  int test;
  double x_error;
  double x_exact[N];
  mgmres_parameters params;
  params.relTol = 1.0E-08;
  params.absTol = 1.0E-08;

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Test MGMRES_ST on the simple -1,2-1 matrix.\n" );
  /*
    Set the matrix.
    Note that we use zero based index valuesin IA and JA.
  */
  k = 0;

  for ( i = 0; i < n; i++ )
    {
      if ( 0 < i )
        {
          ia[k] = i;
          ja[k] = i-1;
          a[k] = -1.0;
          k = k + 1;
        }

      ia[k] = i;
      ja[k] = i;
      a[k] = 2.0;
      k = k + 1;

      if ( i < n - 1 )
        {
          ia[k] = i;
          ja[k] = i+1;
          a[k] = -1.0;
          k = k + 1;
        }

    }
  /*
    Set the right hand side:
  */
  rhs.elements[n-1] = ( double ) ( n + 1 );
  /*
    Set the exact solution.
  */
  for ( i = 0; i < n; i++ )
    {
      x_exact[i] = ( double ) ( i + 1 );
    }

  for ( test = 1; test <= 3; test++ )
    {
      /*
        Set the initial solution estimate.
      */
      x_error = 0.0;
      for ( i = 0; i < n; i++ )
        {
          x_error = x_error + pow ( x_exact[i], 2 );
        }
      x_error = sqrt ( x_error );

      if ( test == 1 )
        {
          params.outerItr = 1;
          params.innerItr = 20;
        }
      else if ( test == 2 )
        {
          params.outerItr = 2;
          params.innerItr = 10;
        }
      else if ( test == 3 )
        {
          params.outerItr = 5;
          params.innerItr = 4;
        }

      printf ( "\n" );
      printf ( "  Test %d\n", test );
      printf ( "  Matrix order N = %d\n", n );
      printf ( "  Inner iteration limit = %d\n", params.innerItr );
      printf ( "  Outer iteration limit = %d\n", params.outerItr );
      printf ( "  Initial X_ERROR = %g\n", x_error );
      sparse_matrix m;
      sparse_matrix_init(&m, n, nz_num);
      m.elements.elements = a;
      m.rows = ia;
      m.cols = ja;

      vector x_estimate;
      vector_init(&x_estimate, N);

      mgmres(&m, &x_estimate, &rhs, params);

      x_error = 0.0;
      for ( i = 0; i < n; i++ )
        {
          x_error = x_error + pow ( x_exact[i] - x_estimate.elements[i], 2 );
        }
      x_error = sqrt ( x_error );

      printf ( "  Final X_ERROR = %g\n", x_error );
    }

  return;
# undef N
# undef NZ_NUM
}

void test_matrix_print() {
  int n = 16;
  int offsets[2] = {2, -1};
  vector v[2];
  vector_init(&v[0], n-offsets[0]); vector_init(&v[1], n-offsets[1]);
  vector_fill(&v[0], 5); vector_fill(&v[1], -3);

  sparse_matrix m;
  sparse_matrix_banded_init(&m, n, v, offsets, 2);

  sparse_matrix_print(&m, n);
}

int main() {
  test_mgmres();
  test_matrix_print();
}
