#include "mkl_spblas.h"
#include "mm-utils.h"
#include "omp.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, const char** argv) {
    if(argc != 3){
        fprintf(stderr, "Argument count doesn't match.\nCorrect usage:\t%s <matrix-filepath> <number-of-multiplications>\n",argv[0]);
        exit(EXIT_FAILURE);
    }
    COO coomatrix = ReadSparseMatrix((char*)(argv[1]));
    printf("%s has %d non-zeros, %d rows, %d columns\n",argv[1],coomatrix.nnz,coomatrix.m,coomatrix.n);
    
    sparse_matrix_t cooA ,csrA;
    sparse_status_t status;

    status = mkl_sparse_d_create_coo(&cooA, SPARSE_INDEX_BASE_ZERO, coomatrix.m, coomatrix.n, coomatrix.nnz, coomatrix.I, coomatrix.J, coomatrix.val);

    if(status != SPARSE_STATUS_SUCCESS){
        fprintf(stderr, "ERROR: failure to create mkl coo from internal coo.\n");
        exit(EXIT_FAILURE);
    }

    status = mkl_sparse_convert_csr(cooA, SPARSE_OPERATION_NON_TRANSPOSE, &csrA);
    if(status != SPARSE_STATUS_SUCCESS){
        fprintf(stderr, "ERROR: failure to convert mkl coo to mkl csr.\n");
        exit(EXIT_FAILURE);
    }

    double *x = malloc(sizeof(double)*coomatrix.n);
    for (int i = 0; i < coomatrix.n; i++)
        x[i] = i;
    double *y = malloc(sizeof(double)*coomatrix.m);

    struct matrix_descr descr = {
        SPARSE_MATRIX_TYPE_GENERAL,
        SPARSE_FILL_MODE_FULL,
        SPARSE_DIAG_NON_UNIT
    };

    int num_ops = atoi(argv[2]);

    double start = omp_get_wtime();
    for(int i = 0; i < num_ops; i++)
        mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, csrA, descr, x, 1.0, y);
    double end = omp_get_wtime();

    printf("%s Time taken: %lf\n",argv[1],end-start);

    // mkl_sparse_s_mv(const sparse_operation_t operation, const float alpha, const sparse_matrix_t A, const struct matrix_descr descr, const float *x, const float beta, float *y)
    return 0;
}