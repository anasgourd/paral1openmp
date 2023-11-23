#ifndef MATRIX_OPERATIONS_H
#define MATRIX_OPERATIONS_H

// Function declarations related to matrix operations

void coo_to_csr(int *I, int *J, double *val, int nrows, int ncols, int nnz,
                int **row_ptr, int **col_ind, double **values);

int read_matrix_market_to_csr(const char *filename, int *nrows, int *ncols, int *nnz,
                              int **row_ptr, int **col_ind, double **values);

// Add more function declarations as needed

#endif // MATRIX_OPERATIONS_H
