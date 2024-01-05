#pragma once

// Column sorted Coordinate.
typedef struct {
    int* I; // Row index
    int* J; // Column index
    double* val; // value
    int nnz; // number of non-zeros
    int m; // number of rows
    int n; // number of columns
} COO;

// Compressed Sparse Row.
// .I contains the compressed index.
typedef COO CSR;

// Compressed Sparse Column.
// .J contains the compressed index.
typedef COO CSC;


//Note that this function allocates memory.
//It can be freed using freeSparseMatrix.
CSR CSC_to_CSR(CSC *in);

// returns an array containing the rank of each row.
// NEEDS TO BE FREED.
int* CalcWeights(CSR* in);

COO SparseTranspose(COO in);
COO ReadSparseMatrix(char *fname);

void freeSparseMatrix(COO *in);
