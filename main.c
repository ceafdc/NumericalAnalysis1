#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "matrix/matrix.h"

#define MAX(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define SINGULAR_MATRIX_ERROR (1)
#define MAX_IT_ERROR (2)

int
gauss_seidel(
    const matrix *A,
    const matrix *x,
    const matrix *b,
    double e,
    int itmax) {

    int n = A->m;
    double new_value;

    for (int it = 0; it < itmax; it++) {
        double max_delta = 0;
        for (int i = 0; i < n; i++) {
            new_value = b->data[i][0];
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    continue;
                }
                new_value -= A->data[i][j] * x->data[j][0];
            }
            if (fabs(A->data[i][i]) < DBL_EPSILON) {
                return SINGULAR_MATRIX_ERROR;
            }
            new_value /= A->data[i][i];
            max_delta = MAX(max_delta, fabs(x->data[i][0] - new_value));
            x->data[i][0] = new_value;
        }
        if (max_delta < e) {
            return EXIT_SUCCESS;
        }
    }
    return MAX_IT_ERROR;
}

double
norm_diff(const matrix *A, const matrix *B) {
    matrix *C = matrix_sub(A, B);
    double norm = matrix_norm_inf(C);
    matrix_free(C);
    return norm;
}

matrix *
create_matrix(int n) {
    matrix *A = matrix_create_zeros(n, n);
    for (int i = 0; i < n; i++) {
        A->data[i][i] = 4;
        if (i + 1 < n) {
            A->data[i][i + 1] = -1;
            A->data[i + 1][i] = -1;
        }
        if (i + 3 < n) {
            A->data[i][i + 3] = -1;
            A->data[i + 3][i] = -1;
        }
    }

    return A;
}

matrix *
create_vector_1(matrix *A) {
    matrix *b = matrix_create_zeros(A->m, 1);
    for (int i = 0; i < b->m; i++) {
        for (int j = 0; j < A->n; j++) {
            b->data[i][0] += A->data[i][j];
        }
    }
    return b;
}

void
create_data_1(
    int n,
    matrix **A,
    matrix **x,
    matrix **b) {
    *A = create_matrix(n);
    *x = matrix_create_zeros(n, 1);
    *b = create_vector_1(*A);
}

void
test_case_1(int n, double e, int itmax) {
    matrix *A = NULL;
    matrix *x = NULL;
    matrix *b = NULL;
    create_data_1(n, &A, &x, &b);

    int gs_error;
    if ((gs_error = gauss_seidel(A, x, b, e, itmax))) {
        switch (gs_error) {
            case 1:fprintf(stderr, "Division by zero\n");break;
            case 2:fprintf(stderr, "Max iterations\n");break;
        }
    }
    matrix *Ax = matrix_mult(A, x);

    matrix *ones = matrix_create_ones(n, 1);
    printf("||x - 1||∞ = %e\n", norm_diff(x, ones));
    printf("||Ax - b||∞ = %e\n", norm_diff(Ax, b));

    matrix_free(A);
    matrix_free(x);
    matrix_free(b);
    matrix_free(Ax);
    matrix_free(ones);
}

matrix *
create_vector_2(int n) {
    matrix *b = matrix_create_ones(n, 1);
    for (int i = 0; i < n; i++) {
        b->data[i][0] /= i + 1;
    }
    return b;
}

void
create_data_2(
    int n,
    matrix **A,
    matrix **x,
    matrix **b) {
    *A = create_matrix(n);
    *x = matrix_create_zeros(n, 1);
    *b = create_vector_2(n);
}

void
test_case_2(int n, double e, int itmax) {
    matrix *A = NULL;
    matrix *x = NULL;
    matrix *b = NULL;
    create_data_2(n, &A, &x, &b);

    int gs_error;
    if ((gs_error = gauss_seidel(A, x, b, e, itmax))) {
        switch (gs_error) {
            case 1:fprintf(stderr, "Division by zero\n");break;
            case 2:fprintf(stderr, "Max iterations\n");break;
        }
    }
    matrix *Ax = matrix_mult(A, x);

    printf("||Ax - b||∞ = %e\n", norm_diff(Ax, b));

    matrix_free(A);
    matrix_free(x);
    matrix_free(b);
    matrix_free(Ax);
}


int
main(int argc, char *argv[]) {
    double e = 1e-10;
    int itmax = 100000;
    puts("1º caso, n=50");
    test_case_1(50, e, itmax);
    puts("1º caso, n=100");
    test_case_1(100, e, itmax);
    puts("2º caso, n=100");
    test_case_2(100, e, itmax);

    return EXIT_SUCCESS;
}
