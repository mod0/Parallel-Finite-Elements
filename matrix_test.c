#include <stdio.h>
#include <math.h>
#include "vector.h"
#include "matrix.h"

void test01() {
# define N 20
# define NZ_NUM 3 * N - 2

  double a[NZ_NUM];
  int i;
  int ia[NZ_NUM];
  int itr_max;
  int j;
  int ja[NZ_NUM];
  int k;
  int mr;
  int n = N;
  int nz_num = NZ_NUM;
  vector rhs; vector_init(&rhs, N);
  int test;
  double tol_abs;
  double tol_rel;
  double x_error;
  double x_exact[N];

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
      itr_max = 1;
      mr = 20;
    }
    else if ( test == 2 )
    {
      itr_max = 2;
      mr = 10;
    }
    else if ( test == 3 )
    {
      itr_max = 5;
      mr = 4;
    }
    tol_abs = 1.0E-08;
    tol_rel = 1.0E-08;

    printf ( "\n" );
    printf ( "  Test %d\n", test );
    printf ( "  Matrix order N = %d\n", n );
    printf ( "  Inner iteration limit = %d\n", mr );
    printf ( "  Outer iteration limit = %d\n", itr_max );
    printf ( "  Initial X_ERROR = %g\n", x_error );
    sparse_matrix m;
    sparse_matrix_init(&m, n, nz_num);
    m.elements.elements = a;
    m.rows = ia;
    m.cols = ja;

    vector x_estimate;
    vector_init(&x_estimate, N);

    mgmres(&m, &x_estimate, &rhs, itr_max, mr, tol_abs, tol_rel);

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

int main() {
    test01();
}
