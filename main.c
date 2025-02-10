/*

Projekt mnożenia macierzy Filip Hojan & Igor Swiatek

*/

#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <time.h>
//#include <windows.h>
#include "omp.h"
#include <cmath>
#include <math.h>

//u mnie cache ma 12288 kb
#define N 5000
#define PP (3 * 1024 * 1024) // 3 MB w bajtach

FILE *result_file;

void multiply_matrices_IKJ(float matrix_a[N][N], float matrix_b[N][N], float matrix_c[N][N], int N)
{
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < N; k++) {
            float temp = matrix_a[i][k];
            for (int j = 0; j < N; j++) {
                matrix_c[i][j] += temp * matrix_b[k][j];
            }
        }
    }
}

void multiply_matrices_4_loop(float matrix_a[N][N], float matrix_b[N][N], float matrix_c[N][N], int N)
{
    size_t matrix_size = sizeof(float) * N * N;
    int k = (int)ceil((double)matrix_size / PP); 
    int r = (int)ceil((double)N / k); 


    

    for (int j = 0; j < N; j += r) {
        #pragma omp parallel for
        for (int i = 0; i < N; i++) {
            for (int k = 0; k < N; k++) {
                for (int jj = j; jj < j+r-1; jj++) {
                    matrix_c[i][jj] += matrix_a[i][k] * matrix_b[k][jj];
                }
            }
        }
    }
}

void multiply_matrices_6_loop(float matrix_a[N][N], float matrix_b[N][N], float matrix_c[N][N], int N)
{
    int r = floor(sqrt(PP/12.0));
    int ka = floor(n/r);


    for (int i = 0; i < N; i += r) {
        for (int j = 0; j < N; j += r) {
            for (int k = 0; k < N; k += r) {
                #pragma omp parallel for
                for (int ii = i; ii < i + r; ii++) {
                    for (int kk = k; kk < k + r; kk++) {
                        for (int jj = j; jj < j + r; jj++) {
                            matrix_c[ii][jj] += matrix_a[ii][kk] * matrix_b[kk][jj];
                        }
                    }
                }
            }
        }
    }
}

void initialize_matrices(float matrix_a[N][N], float matrix_b[N][N], float matrix_c[N][N], int N)
{
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matrix_a[i][j] = (float)rand() / RAND_MAX;
            matrix_b[i][j] = (float)rand() / RAND_MAX;
            matrix_c[i][j] = 0.0;
        }
    }
}

void clear_result_matrix(float matrix_c[N][N], N)
{
    #pragma omp parallel for 
    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLUMNS; j++) {
            matrix_c[i][j] = 0.0;
        }
    }
}

void print_elapsed_time(double start, char* description)
{
    double elapsed;
    elapsed = (double)clock() / CLK_TCK;
    printf("%s \t\t| czas: %8.4f sec \n", description, elapsed - start);
    fprintf(result_file, "%s \t\t| czas : %8.4f sec (%6.4f sec rozdzielczosc pomiaru)\n", description, elapsed-start, 1.0 / CLK_TCK);
}

void check_results(float matrix_c[N][N], float matrix_c_reference[N][N], int N)
{
    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLUMNS; j++) {
            if (abs(matrix_c_reference[i][j] - matrix_c[i][j]) > 0.001) {
                printf("Błąd: %f != %f w pozycji [%d][%d]\n", matrix_c_reference[i][j], matrix_c[i][j], i, j);
                return;
            }
        }
    }
    printf("Wyniki poprawne.\n");
}

void multiply_sequential(float matrix_a[N][N], float matrix_b[N][N], float matrix_c[N][N], int N)
{
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < N; k++) {
            for (int j = 0; j < N; j++) {
                matrix_c[i][j] += matrix_a[i][k] * matrix_b[k][j];
            }
        }
    }
}

int main(int argc, char* argv[])
{
    double start;

    float matrix_a[N][N];
    float matrix_b[N][N];
    float matrix_c[N][N];
    float matrix_c_reference[N][N];

    initialize_matrices(matrix_a, matrix_b, matrix_c, N);
    print_elapsed_time(clock() / CLK_TCK, "Inicjalizacja macierzy");
    multiply_sequential(matrix_a, matrix_b, matrix_c_reference, N);

    // Test dla funkcji IKJ
    clear_result_matrix(matrix_c, N);
    start = clock() / CLK_TCK;
    multiply_matrices_IKJ(matrix_a, matrix_b, matrix_c, N);
    print_elapsed_time(start, "Mnozenie macierzy IKJ");
    check_results(matrix_c, matrix_c_reference, N);

    // Test dla funkcji 4-pętlowej IKJ
    clear_result_matrix(matrix_c, N);
    start = clock() / CLK_TCK;
    multiply_matrices_4_loop(matrix_a, matrix_b, matrix_c, N);
    print_elapsed_time(start, "Mnozenie macierzy 4-petlowe");
    check_results(matrix_c, matrix_c_reference, N);

    // Test dla funkcji 6-pętlowej  IKJ
    clear_result_matrix(matrix_c, N);
    start = clock() / CLK_TCK;
    multiply_matrices_6_loop(matrix_a, matrix_b, matrix_c, N);
    print_elapsed_time(start, "Mnozenie macierzy 6-petlowe");
    check_results(matrix_c, matrix_c_reference, N);

    //fclose(result_file);
    return 0;
}
