#include <stdio.h>
#include <stdlib.h>

#include "matrix/matrix.h"

int
gauss_seidel(
    matrix *A,
    matrix *b,
    double e,
    int itmax) {

    return 0;
}

matrix *
create_matrix(n) {
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
create_vector(matrix *A) {
    matrix *b = matrix_create_zeros(A->n, 1);
    for (int i = 0; i < A->n; i++) {
        for (int j = 0; j < A->m; j++) {
            b->data[i][0] += A->data[i][j];
        }
    }
    return b;
}

void
create_data(
    int n,
    matrix **A,
    matrix **b) {
    *A = create_matrix(n);
    *b = create_vector(*A);
}


int
main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s n\n", *argv);
        return EXIT_FAILURE;
    }
    int n = atoi(argv[1]);
    if (n <= 0) {
        fprintf(stderr, "n must be grater than 0\n");
    }

    matrix *A = NULL;
    matrix *b = NULL;
    create_data(n, &A, &b);

    matrix_print(A);
    puts("");
    matrix_print(b);

    matrix_free(A);
    matrix_free(b);

    return EXIT_SUCCESS;
}
