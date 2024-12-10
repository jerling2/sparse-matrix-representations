#include <stdlib.h>
#include <stdio.h>

#define PADDING_ELEMENT ~0


void print_csr(int m, int nnz, unsigned int* csr_row_ptr, unsigned int* csr_col_ind, double* csr_vals)
{
    for (int i = 0; i < m + 1; i++)
    {
        if (i == 0) printf("CSR row: ");
        if (i != 0) printf(", ");
        printf("%u", csr_row_ptr[i]);
        if (i == m) printf("\n");
    }

    for (int i = 0; i < nnz; i++)
    {
        if (i == 0) printf("CSR col: ");
        if (i != 0) printf(", ");
        printf("%u", csr_col_ind[i]);
        if (i == nnz - 1) printf("\n");
    }

    for (int i = 0; i < nnz; i++)
    {
        if (i == 0) printf("CSR nnz: ");
        if (i != 0) printf(", ");
        printf("%.0f", csr_vals[i]);
        if (i == nnz - 1) printf("\n");
    }
}


void print_ell(int m, int n, unsigned int* ell_col_ind, double* ell_vals)
{
    for (int i = 0; i < m * n; i++)
    {
        if (i == 0) printf("ELL col: ");
        if (i != 0) printf(", ");
        if (ell_col_ind[i] == PADDING_ELEMENT) printf("\x1b[31;1m*\x1b[0m");
        else printf("%u", ell_col_ind[i]);
        if (i == m * n - 1) printf("\n");
    }
    for (int i = 0; i < m * n; i++)
    {
        if (i == 0) printf("ELL nnz: ");
        if (i != 0) printf(", ");
        if (ell_col_ind[i] == PADDING_ELEMENT) printf("\x1b[31;1m*\x1b[0m");
        else printf("%.0f", ell_vals[i]);
        if (i == m * n - 1) printf("\n");
    }

}


void convert_sparse_to_csr(double MATRIX[4][4], int m, int n, int *nnz,
                           unsigned int** csr_row_ptr, 
                           unsigned int** csr_col_ind, double** csr_vals)
{
    // First pass, count the nnz elements    
    for (int i = 0; i < m; i++)
        for (int l = 0; l < n; l++)
            if (MATRIX[i][l] != 0) (*nnz)++;

    *csr_row_ptr = (unsigned int*) malloc(sizeof(unsigned int) * (m + 1));
    *csr_col_ind = (unsigned int*) malloc(sizeof(unsigned int) * (*nnz));
    *csr_vals = (double *) malloc(sizeof(double) * (*nnz));

    **csr_row_ptr = 0;
    int row_ptr = 0;
    for (int i = 0; i < m; i ++) {
        for (int l = 0; l < n; l++) {
            if (MATRIX[i][l] != 0) {
                (*csr_col_ind)[row_ptr] = l;
                (*csr_vals)[row_ptr] = MATRIX[i][l];
                row_ptr ++;
            }
        }
        (*csr_row_ptr)[i+1] = row_ptr;
    }
}


void convert_csr_to_ell(unsigned int* csr_row_ptr, unsigned int* csr_col_ind,
                        double* csr_vals, int m, int n, int nnz, 
                        unsigned int** ell_col_ind, double** ell_vals, 
                        int* n_new)
{
    // First pass, count the minimum length of each row;
    *n_new = 0;
    for (int i = 1; i < m + 1; i++) {
        int row_len = csr_row_ptr[i] - csr_row_ptr[i-1];
        *n_new = row_len > *n_new ? row_len : *n_new;
    }

    unsigned int** ind_big_buf = (unsigned int**) malloc(sizeof(unsigned int*) * m);
    double** val_big_buf = (double**) malloc(sizeof(double*) * m);
    int count = 0;

    // Second pass to construct ell in row-major order.
    for (int i = 1; i < m + 1; i++) {
        unsigned int* ell_ind_buf = (unsigned int*) malloc(sizeof(unsigned int) * *n_new);
        double* ell_val_buf = (double*) malloc(sizeof(double) * *n_new);
        int row_len = csr_row_ptr[i] - csr_row_ptr[i-1];
        for (int l = 0; l < row_len; l++) {
            ell_ind_buf[l] = csr_col_ind[count + l];
            ell_val_buf[l] = csr_vals[count + l];
        }
        for (int k = row_len; k < *n_new; k++) {
            ell_ind_buf[k] += PADDING_ELEMENT;
            ell_val_buf[k] += PADDING_ELEMENT;
        }
        count += row_len;
        ind_big_buf[i-1] = ell_ind_buf;
        val_big_buf[i-1] = ell_val_buf;
    }

    *ell_col_ind = (unsigned int*) malloc(sizeof(unsigned int) * (*n_new * m));
    *ell_vals = (double*) malloc(sizeof(double) * (*n_new * m));
    count = 0;

    // Third pass to construct ell in col-major order.
    for (int j = 0; j < *n_new; j++) {
        for (int i = 0; i < m; i++) {
            (*ell_col_ind)[count] = ind_big_buf[i][j];
            (*ell_vals)[count] = val_big_buf[i][j];
            count++;
        }
    }
    
    // Free buffers
    for (int i = 0; i < m; i++) {
        free(ind_big_buf[i]);
        free(val_big_buf[i]);
    }
    free(ind_big_buf);
    free(val_big_buf);
}


int main(int argc, char **argv)
{
    int NUM_ROWS = 4;
    int NUM_COLS = 4;
    
    double MATRIX[4][4] = {
        {1,7,0,0},
        {5,0,3,9},
        {0,2,8,0},
        {0,0,0,6}
    };

    int nnz;
    unsigned int* csr_row_ptr;
    unsigned int* csr_col_ind;
    double* csr_vals;

    convert_sparse_to_csr(
        MATRIX, NUM_ROWS, NUM_COLS, 
        &nnz, &csr_row_ptr, &csr_col_ind, &csr_vals
    );

    print_csr(NUM_ROWS, nnz, csr_row_ptr, csr_col_ind, csr_vals);

    unsigned int* ell_col_ind;
    double* ell_vals;
    int new_n;

    convert_csr_to_ell(
        csr_row_ptr, csr_col_ind, csr_vals, NUM_ROWS, NUM_COLS, nnz, 
        &ell_col_ind, &ell_vals, &new_n
    );

    print_ell(NUM_ROWS, new_n, ell_col_ind, ell_vals);

    free(csr_row_ptr); free(csr_col_ind); free(csr_vals);
    free(ell_col_ind); free(ell_vals);

    return 0;
}