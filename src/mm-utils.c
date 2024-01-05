#include "include/mm-utils.h"
#include <stdio.h>
#include <stdlib.h>
#include "mmio.h"

//Note that this function allocates memory.
//It can be freed using freeSparseMatrix.
CSR CSC_to_CSR(CSC *in){
    CSR out = {
        .I = (int *)malloc(sizeof(int) * (in->m + 1)),
        .J = (int *)malloc(sizeof(int) * in->nnz),
        .val = (double *)malloc(sizeof(double) * in->nnz),
        .nnz = in->nnz,
        .m = in->m,
        .n = in->n
    };
    if (out.I == NULL || out.J == NULL || out.val == NULL){
        return (CSR){0};
    }
    out.I[0] = 0;
    out.I[out.m] = out.nnz;
    for(int i = 0, row = 0, col = 0; i < out.nnz; col++){
        if(col >= in ->n){
            row++;
            out.I[row] = i;
            col = 0;
        }
        for(int j[2] = {in->J[col], in->J[col+1]}; j[0] < j[1] & in->I[j[0]] <= row; j[0]++){
            if(in->I[j[0]]== row){
                out.J[i] = col;
                i++;
                break;
            }
        }
    }
    return out;
}

// returns an array containing the rank of each row.
// NEEDS TO BE FREED.
int* CalcWeights(CSR* in){
    int *out = (int*)malloc(in->m*sizeof(int));
    for(int i = 0; i < in->m; i++){
        out[i]= in->I[i+1] - in->I[i];
    }
    return out;
}

COO SparseTranspose(COO in){
    int *temp = in.I;
    in.I = in.J;
    in.J = temp;
    int swap = in.m;
    in.m = in.n;
    in.n = swap;
    return in;
}

void freeSparseMatrix(COO *in){
    free(in->I);
    free(in->J);
    free(in->val);
}


CSC COO_to_CSC(COO *in){
    CSC out = {
        .I = (int *)malloc(sizeof(int) * in->nnz),
        .J = (int *)malloc(sizeof(int) * (in->n + 1)),
        .val = (double *)malloc(sizeof(double) * in->nnz),
        .nnz=in->nnz,
        .m=in->m,
        .n=in->n
    };
    if (out.I == NULL || out.J == NULL || out.val == NULL){
        return (CSC){0};
    }
    out.J[0] = 0;
    out.J[out.n] = out.nnz;
    int i, col;
    for(i = 0, col = 0; i < out.nnz; i++){
        if(in->J[i] != col)
            out.J[col++ + 1] = i;
        out.I[i] = in->I[i];
        out.val[i] = in->val[i];
    }
    int last = out.J[out.n];
    if (last != out.nnz){
      out.J[out.n] = out.nnz;
    }
    return out;
}

COO ReadSparseMatrix(char *fname) {
  int ret_code;
  MM_typecode matcode;
  FILE *f;
  int M, N, nz;

  if ((f = fopen(fname, "r")) == NULL) {
    fprintf(stderr, "Could not open matrix file\n");
    exit(1);
  }

  if (mm_read_banner(f, &matcode) != 0) {
    printf("Could not process Matrix Market banner.\n");
    exit(1);
  }

  if (mm_is_complex(matcode) | !mm_is_matrix(matcode) | !mm_is_sparse(matcode) |
      mm_is_symmetric(matcode)) {
    printf("Sorry, this application does not support ");
    printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
    exit(1);
  }

  /* find out size of sparse matrix .... */

  if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) != 0) {
    exit(1);
  }

  /* reseve memory for matrices */

  COO smatrix = {.I = (int *)malloc(nz * sizeof(int)),
                 .J = (int *)malloc(nz * sizeof(int)),
                 .val = (double *)malloc(nz * sizeof(double)),
                 .nnz = nz,
                 .m = M,
                 .n = N};

  /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
  /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
  /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

  for (int i = 0; i < nz; i++) {
    if (fscanf(f, "%d %d %lg\n", &smatrix.I[i], &smatrix.J[i],
               &smatrix.val[i]) != 3) {
      fprintf(stderr, "Error reading matrix\n");
      exit(EXIT_FAILURE);
    }
    smatrix.I[i]--; /* adjust from 1-based to 0-based */
    smatrix.J[i]--;
  }
  fclose(f);

  // CSC cscmatrix = {0};
  // if ((cscmatrix = COO_to_CSC(&smatrix)).nnz == 0)
  //   exit(EXIT_FAILURE);
  // freeSparseMatrix(&smatrix);

  // CSR csrmatrix = {0};
  // if((csrmatrix = CSC_to_CSR(&cscmatrix)).nnz == 0)
  //   exit(EXIT_FAILURE);
  // freeSparseMatrix(&cscmatrix);

  return smatrix;
}