#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "matrix.h"
#include "error.h"

// Initializes the members of a sparse_matrix struct
void sparse_matrix_init(sparse_matrix* m, size_t size, size_t numNonzero) {
    m->size = size;
    vector_init(&m->elements, numNonzero);
    m->rows = (int*) calloc(numNonzero, sizeof(int));
    m->cols = (int*) calloc(numNonzero, sizeof(int));
}

// Helper function that writes a band to a sparse_matrix.  If the band is not on
// the main diagonal, it will write it twice, since the matrix is symmetric
static void fill_band(sparse_matrix* m, vector* band, int* elementIdx) {
    int bandOffset = m->size - band->size;

    int i;
    for(i = 0; i < band->size; i++) {
        double val = band->elements[i];

        m->rows[*elementIdx] = i;
        m->cols[*elementIdx] = i + bandOffset;
        m->elements.elements[*elementIdx] = val;

        (*elementIdx)++;

        if (bandOffset != 0) { // Write the band in the lower-left if the band is not on the diagonal
            m->rows[*elementIdx] = i + bandOffset;
            m->cols[*elementIdx] = i;
            m->elements.elements[*elementIdx] = val;

            (*elementIdx)++;
        }
    }
}


void sparse_symmetric_banded_init(sparse_matrix* m, size_t size, vector* bands, size_t numBands) {
    int i;
    int numNonzero = 0;
    for (i = 0; i < numBands; i++) {
        int bandLength = bands[i].size;
        numNonzero += (bandLength == size ? 1 : 2) * bandLength;
    }
    sparse_matrix_init(m, size, numNonzero);

    int curIdx = 0;
    for (i = 0; i < numBands; i++) {
        fill_band(m, &bands[i], &curIdx);
    }
}


double sparse_matrix_get(sparse_matrix* m, int row, int col) {
    int i;
    for (i = 0; i < m->elements.size; i++) {
        if (m->rows[i] == row && m->cols[i] == col) {
            return m->elements.elements[i];
        }
    }
    return 0;
}


static void sparse_matrix_array_multiply(sparse_matrix* m, double* x, double* w) {
    int i;
    int j;
    int k;

    for (i = 0; i < m->size; i++) {
        w[i] = 0.0;
    }

    for (k = 0; k < m->elements.size; k++) {
        i = m->rows[k];
        j = m->cols[k];
        w[i] += m->elements.elements[k] * x[j];
    }
}


void sparse_matrix_vector_multiply(sparse_matrix* m, vector* x, vector* w) {
    sparse_matrix_array_multiply(m, x->elements, w->elements);
}


double sparse_matrix_frobenius_norm(sparse_matrix* m) {
    return vector_2_norm(&m->elements);
}


static double r8vec_dot(int n, double a1[], double a2[]) {
    int i;
    double value;

    value = 0.0;
    for (i = 0; i < n; i++) {
        value = value + a1[i] * a2[i];
    }
    return value;
}


static void mult_givens(double c, double s, int k, vector* g) {
    double g1;
    double g2;

    g1 = c * g->elements[k] - s * g->elements[k+1];
    g2 = s * g->elements[k] + c * g->elements[k+1];

    g->elements[k]   = g1;
    g->elements[k+1] = g2;
}

static double** create_dense_matrix(int nrl, int nrh, int ncl, int nch) {
    int i;
    double **m;
    int nrow = nrh - nrl + 1;
    int ncol = nch - ncl + 1;
    /*
    Allocate pointers to the rows.
    */
    m = (double**) malloc((size_t) ((nrow + 1) * sizeof(double*)));

    m = m + 1;
    m = m - nrl;
    /*
    Allocate each row and set pointers to them.
    */
    m[nrl] = (double*) malloc((size_t) ((nrow * ncol + 1) * sizeof(double)));

    m[nrl] = m[nrl] + 1;
    m[nrl] = m[nrl] - ncl;

    for (i = nrl + 1; i <= nrh; i++) {
        m[i] = m[i-1] + ncol;
    }

    return m;
}


static void free_dense_matrix(double **m, int nrl, int nrh, int ncl, int nch) {
    free((char*) (m[nrl] + ncl - 1));
    free((char*) (m + nrl - 1));
}

// Restarted GMRES adapted from http://people.sc.fsu.edu/~jburkardt/c_src/mgmres/mgmres.html
void mgmres(sparse_matrix* m, vector* x, vector* rhs, int itr_max, int mr, double tol_abs, double tol_rel) {
    double av;
    vector c;
    double delta = 1.0e-03;
    vector g;
    double **h;
    double htmp;
    int i;
    int itr;
    int itr_used;
    int j;
    int k;
    int k_copy;
    double mu;
    vector r;
    double rho;
    double rho_tol;
    vector s;
    double **v;
    vector y;

    int n = m->size;
    if (n < mr) {
        error("MGMRES N < MR");
    }

    itr_used = 0;

    vector_init(&c, mr);
    vector_init(&g, mr + 1);
    vector_init(&r, n);
    vector_init(&s, mr);
    vector_init(&y, mr + 1);
    h = create_dense_matrix(0, mr, 0, mr - 1);
    v = create_dense_matrix(0, mr, 0, n - 1);

    for (itr = 0; itr < itr_max; itr++) {
        sparse_matrix_vector_multiply(m, x, &r);

        vector_subtract(rhs, &r, &r);
        rho = vector_2_norm(&r);

        if (itr == 0) {
            rho_tol = rho * tol_rel;
        }

        for (i = 0; i < n; i++) {
            v[0][i] = r.elements[i] / rho;
        }

        vector_fill(&g, 0);
        g.elements[0] = rho;

        for (i = 0; i < mr + 1; i++) {
            for (j = 0; j < mr; j++) {
                h[i][j] = 0.0;
            }
        }

        for (k = 0; k < mr; k++) {
            k_copy = k;

            sparse_matrix_array_multiply(m, v[k], v[k+1]);

            av = sqrt(r8vec_dot(n, v[k+1], v[k+1]));

            for (j = 0; j < k+1; j++) {
                h[j][k] = r8vec_dot (n, v[k+1], v[j]);
                for (i = 0; i < n; i++) {
                    v[k+1][i] -= h[j][k] * v[j][i];
                }
            }

            h[k+1][k] = sqrt(r8vec_dot(n, v[k+1], v[k+1]));

            if ((av + delta * h[k+1][k]) == av) {
                for (j = 0; j < k+1; j++) {
                    htmp = r8vec_dot(n, v[k+1], v[j]);
                    h[j][k] += htmp;
                    for (i = 0; i < n; i++) {
                        v[k+1][i] -= htmp * v[j][i];
                    }
                }
                h[k+1][k] = sqrt(r8vec_dot(n, v[k+1], v[k+1]));
            }

            if (h[k+1][k] != 0.0) {
                for (i = 0; i < n; i++) {
                    v[k+1][i] /= h[k+1][k];
                }
            }

            if (0 < k) {
                for (i = 0; i < k + 2; i++) {
                    y.elements[i] = h[i][k];
                }
                for (j = 0; j < k; j++) {
                    mult_givens(c.elements[j], s.elements[j], j, &y);
                }
                for (i = 0; i < k + 2; i++) {
                    h[i][k] = y.elements[i];
                }
            }

            mu = sqrt(h[k][k] * h[k][k] + h[k+1][k] * h[k+1][k]);
            c.elements[k] = h[k][k] / mu;
            s.elements[k] = -h[k+1][k] / mu;
            h[k][k] = c.elements[k] * h[k][k] - s.elements[k] * h[k+1][k];
            h[k+1][k] = 0.0;
            mult_givens(c.elements[k], s.elements[k], k, &g);

            rho = fabs(g.elements[k+1]);

            itr_used = itr_used + 1;

            if (rho <= rho_tol && rho <= tol_abs) {
                break;
            }
        }

        k = k_copy;

        y.elements[k] = g.elements[k] / h[k][k];
        for (i = k - 1; 0 <= i; i--) {
            y.elements[i] = g.elements[i];
            for (j = i+1; j < k + 1; j++) {
                y.elements[i] -= h[i][j] * y.elements[j];
            }
            y.elements[i] /= h[i][i];
        }

        for (i = 0; i < n; i++) {
            for (j = 0; j < k + 1; j++) {
                x->elements[i] += v[j][i] * y.elements[j];
            }
        }

        if (rho <= rho_tol && rho <= tol_abs) {
            break;
        }
    }

    vector_free(&c);
    vector_free(&g);
    vector_free(&r);
    vector_free(&s);
    vector_free(&y);
    free_dense_matrix(v, 0, mr, 0, n - 1);
    free_dense_matrix(h, 0, mr, 0, mr - 1);
}

typedef void (*PrinterFunction)(sparse_matrix*, int, int);
#define PRINT_NUM_DIGITS 8

static void printEllipsis(sparse_matrix* m, int row, int col) {
	printf("%*s ", PRINT_NUM_DIGITS + 5, "|");
}

static void printElement(sparse_matrix* m, int row, int col) {
	printf("%+*e ", PRINT_NUM_DIGITS, sparse_matrix_get(m, row, col));
}

static void printRow(sparse_matrix* m, size_t ellipsisThreshold, int row, PrinterFunction printer) {
    int numCols = m->size;

    if (numCols > ellipsisThreshold) {
        int numBeginningCols = ellipsisThreshold / 2;
        int numEndingCols = ellipsisThreshold - numBeginningCols - 1;
        int c;
        for (c = 0; c < numBeginningCols; c++) {
            printer(m, row, c);
        }

        printf("… ");

        for (c = numCols - numEndingCols; c < numCols; c++) {
            printer(m, row, c);
        }
    } else {
        int c;
        for (c = 0; c < numCols; c++) {
            printer(m, row, c);
        }
    }
    printf("\n");
}

void sparse_matrix_print(sparse_matrix* m, size_t ellipsisThreshold) {
    int numRows = m->size;

    if (numRows > ellipsisThreshold) {
        int numBeginningRows = ellipsisThreshold / 2;
        int numEndingRows = ellipsisThreshold - numBeginningRows - 1;
        int r;
        for (r = 0; r < numBeginningRows; r++) {
            printRow(m, ellipsisThreshold, r, &printElement);
        }

        printRow(m, ellipsisThreshold, r, &printEllipsis);

        for (r = numRows - numEndingRows; r < numRows; r++) {
            printRow(m, ellipsisThreshold, r, &printElement);
        }
    } else {
        int r;
        for (r = 0; r < numRows; r++) {
            printRow(m, ellipsisThreshold, r, &printElement);
        }
    }
}


void sparse_matrix_free(sparse_matrix* m) {
    vector_free(&m->elements);
    free(m->rows);
    free(m->cols);
}
