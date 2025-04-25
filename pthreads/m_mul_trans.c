
#include<stdio.h>
#include<stdlib.h>
#include<time.h>

typedef struct matrix {
    double * data;
    unsigned rows, cols;
} matrix;

double * transpose(matrix *M) {
    int i, j;

    double *temp = malloc(M->rows*M->cols*sizeof(double)); 

    for(i = 0; i < M->rows; i++) {
        for(j = 0; j < M->cols; j++) {
            temp[i*M->rows + j] = M->data[j*M->rows + i];
        }
    }

    return temp;
}

matrix* mmul(matrix A, matrix B) {
    matrix* r_matrix = malloc(sizeof(matrix));
    
    int i, j, k, l;
    double acumulator = 0;
    
    r_matrix->data = malloc(A.rows * B.cols * sizeof(double));
    r_matrix->rows = A.rows;
    r_matrix->cols = B.cols;

    B.data = transpose(&B);

    for(i = 0; i < A.rows; i++) {
        for(j = 0; j < A.cols; j++) {
            acumulator = 0.0;
            for(k = 0; k < B.rows; k++) {
                acumulator += A.data[i*A.rows + k] * B.data[i*A.rows + k];
            }
            r_matrix->data[i*A.rows + j] = acumulator;
        }
    }

    return r_matrix;
}

void initialize_matrix(matrix *M) {
    int i, j;

    for(i = 0; i < M->rows; i++)
        for(int j = 0; j < M->cols; j++) {
            M->data[i*M->rows+j] = i*M->rows+j * 1.0f; 
        }
}

int main(int argc, char const *argv[]) {

    int n = 2000, i;

    clock_t begin,
            end;

    double time_spent;

    matrix* result;
    matrix A = {malloc(n*n*sizeof(double)), n, n};
    matrix B = {malloc(n*n*sizeof(double)), n, n};

    initialize_matrix(&A);
    initialize_matrix(&B);

    begin = clock();

    result = mmul(A, B);

    end = clock();

    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    printf("Execution time: %.4lfs for A(%d, %d) and B(%d, %d)\n", time_spent, A.rows, A.rows, A.rows, A.rows);

    free(A.data);           
    free(B.data);           
    free(result->data);     
    free(result);           

    return 0;
}
