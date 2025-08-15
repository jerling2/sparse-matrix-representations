# csr_to_ell.c

This program demonstrates the conversion between two common sparse matrix storage formats: Compressed Sparse Row (CSR) and ELLPACK (ELL).

## Overview

Sparse matrices are matrices in which most elements are zero. Efficient storage and computation with sparse matrices require specialized formats to save space and computation time. Two popular formats are:

- **CSR (Compressed Sparse Row)**: Stores nonzero values in a 1D array, along with two arrays indicating column indices and row boundaries.
- **ELL (ELLPACK)**: Stores nonzero values row-wise in a dense rectangular format, padding rows with fewer nonzeros to ensure uniform length.

## What the Program Does

1. **Sparse Matrix to CSR Conversion**  
   The program starts with a fixed 4x4 dense matrix (with mostly zero values) and converts it to the CSR format.  
   - The nonzero values, their column indices, and pointers to the start of each row are extracted.

2. **CSR to ELL Conversion**  
   The CSR representation is then converted to the ELL format.  
   - Rows are padded with special values to ensure all have the same length as the row with the most non-zeros.
   - The ELL format is stored in column-major order for efficient access.

3. **Printing Functions**  
   The program includes utilities to print CSR and ELL representations to the console, making it easy to visualize how the sparse matrix is stored in each format. Special colored output highlights padding elements in the ELL format.

4. **Memory Management**  
   All dynamic memory allocations are properly freed before program termination.

## Sample Matrix

The conversion is demonstrated on the following matrix:

```
[ 1, 7, 0, 0 ]
[ 5, 0, 3, 9 ]
[ 0, 2, 8, 0 ]
[ 0, 0, 0, 6 ]
```

## Key Functions

- `convert_sparse_to_csr`: Converts a dense matrix to CSR format.
- `convert_csr_to_ell`: Converts CSR format to ELL format.
- `print_csr` / `print_ell`: Print the matrix data in CSR or ELL format, respectively.

## How to Compile and Run

To compile the code, use a C compiler such as `gcc`:

```sh
gcc -o csr_to_ell csr_to_ell.c
```

To run the compiled program:

```sh
./csr_to_ell
```

The output will display the CSR and ELL representations for the included test matrix.

---
