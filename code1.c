#include <stdio.h>
#include <stdlib.h>
#include "mmio.h"
#include "matrix_operations.h"


void coo_to_csr(int *I, int *J, double *val, int nrows, int ncols, int nnz,
                int **row_ptr, int **col_ind, double **values) {
    // Allocate memory for CSR arrays
    *row_ptr = (int *)malloc((nrows + 1) * sizeof(int));
    *col_ind = (int *)malloc(nnz * sizeof(int));
    *values = (double *)malloc(nnz * sizeof(double));

    // Initialize row_ptr with zeros
    for (int i = 0; i <= nrows; i++) {
        (*row_ptr)[i] = 0;
    }

    // Count the number of non-zero elements in each row
    for (int i = 0; i < nnz; i++) {
        (*row_ptr)[I[i]]++;
    }

    // Cumulative sum to get the starting index of each row
    int sum = 0;
    for (int i = 0; i < nrows; i++) {
        int temp = (*row_ptr)[i];
        (*row_ptr)[i] = sum;
        sum += temp;
    }
    (*row_ptr)[nrows] = nnz;

    // Fill col_ind and values arrays
    for (int i = 0; i < nnz; i++) {
        int row = I[i];
        int index = (*row_ptr)[row];
        (*col_ind)[index] = J[i];
        (*values)[index] = val[i];
        (*row_ptr)[row]++;
    }

    // Reset row_ptr
    for (int i = nrows; i > 0; i--) {
        (*row_ptr)[i] = (*row_ptr)[i - 1];
    }
    (*row_ptr)[0] = 0;
}

int read_matrix_market_to_csr(const char *filename, int *nrows, int *ncols, int *nnz,
                              int **row_ptr, int **col_ind, double **values) {
    // Open the Matrix Market file for reading
    FILE *f = fopen(filename, "r");
    if (f == NULL) {
        fprintf(stderr, "Error opening file.\n");
        return 1;
    }

    // Read the Matrix Market header
    MM_typecode matcode;
    if (mm_read_banner(f, &matcode) != 0) {
        fprintf(stderr, "Error reading Matrix Market banner.\n");
        fclose(f);
        return 1;
    }

    // Check the matrix type
    if (!mm_is_matrix(matcode) || !mm_is_sparse(matcode) || !mm_is_real(matcode)) {
        fprintf(stderr, "This example only works with real sparse matrices.\n");
        fclose(f);
        return 1;
    }

    // Read matrix dimensions and number of non-zero elements
    if (mm_read_mtx_crd_size(f, nrows, ncols, nnz) != 0) {
        fprintf(stderr, "Error reading matrix dimensions.\n");
        fclose(f);
        return 1;
    }

    // Allocate memory for the sparse matrix in COO format
    int *I = (int *)malloc(*nnz * sizeof(int));
    int *J = (int *)malloc(*nnz * sizeof(int));
    double *val = (double *)malloc(*nnz * sizeof(double));

    // Read the sparse matrix data
    for (int i = 0; i < *nnz; i++) {
        fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
        I[i]--;  // Adjust indices to 0-based
        J[i]--;
    }

    // Close the file
    fclose(f);

    // Convert COO to CSR
    coo_to_csr(I, J, val, *nrows, *ncols, *nnz, row_ptr, col_ind, values);

    // Free allocated memory for COO format
    free(I);
    free(J);
    free(val);

    return 0;
}
