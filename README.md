# Sparse Matrix Multiplication (OpenMP Parallel)

## Overview

This C project implements efficient sparse matrix multiplication using CSR and CSC formats with OpenMP parallelization.
Intermediate matrices are computed indirectly via a vector that indicates the cluster assignment of each row. During multiplication, this vector is used to accumulate values directly into the intermediate sparse matrix, avoiding the creation of large sparse matrices, reducing memory usage, and improving cache efficiency.

This project is the parallel extension of the sequential implementation and was developed for the Parallel and Distributed Systems course at ECE AUTH (2023–2024).

Input matrices are in Matrix Market (.mtx) format (for example, this project uses the matrix: `af23560.mtx`).

**Compile :**

```bash
gcc main_code.c code1.c mmio.c -o omp_sparse_mult
```
**Run:**

```bash
./omp_sparse_mult
```
