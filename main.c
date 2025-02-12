/*

Projekt mnożenia macierzy Filip Hojan & Igor Swiatek

*/

#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <time.h>
#include "omp.h"
#include <math.h>
#include <stdlib.h>

#define N 5000   // lub 3072, 4096 w kolejnych eksperymentach

#define PP (3 * 1024 * 1024)  // 3 MB w bajtach


void multiply_matrices_IKJ(float matrix_a[N][N], float matrix_b[N][N], float matrix_c[N][N])
{
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < N; k++) {
            for (int j = 0; j < N; j++) {
                matrix_c[i][j] += matrix_a[i][k] * matrix_b[k][j];
            }
        }
    }
}

void multiply_matrices_4_loop(float matrix_a[N][N], float matrix_b[N][N], float matrix_c[N][N])
{
    size_t matrix_size = (size_t)(4 * N * N);
    size_t matrix_size_mb = matrix_size / (1024 * 1024); // Teraz poprawnie przeliczamy na MB

    printf("Wielkosc macierzy: %zu\n", matrix_size);
    printf("Wielkosc macierzy w MB: %zu\n", matrix_size_mb);

    int k = (int)ceil((double)matrix_size / PP); 
    int r = (int)ceil((double)N / k); 
    printf("Bloki: %d\n",k );
    printf("Ilosc slow w dlugosci wiersza na blok: %d\n",r );

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

void multiply_matrices_6_loop(float matrix_a[N][N], float matrix_b[N][N], float matrix_c[N][N])
{
    int r = (int)floor(sqrt(PP / (sizeof(float) * 3.0) ));
    int ka = (int)ceil(N/r);
    printf("Podbloki: %d\n",ka );
    printf("Ilosc slow na blok: %d\n",r );

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

void initialize_matrices(float matrix_a[N][N], float matrix_b[N][N], float matrix_c[N][N], float matrix_c_reference[N][N])
{
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matrix_a[i][j] = (float)rand() / RAND_MAX;
            matrix_b[i][j] = (float)rand() / RAND_MAX;
            matrix_c[i][j] = 0.0;
            matrix_c_reference[i][j] = 0.0;
        }
    }
}

void clear_result_matrix(float matrix_c[N][N])
{
    #pragma omp parallel for 
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matrix_c[i][j] = 0.0;
        }
    }
}

int check_results(float matrix_c[N][N], float matrix_c_reference[N][N])
{
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (abs(matrix_c_reference[i][j] - matrix_c[i][j]) > 0.001) {
                printf("Błąd: %f != %f w pozycji [%d][%d]\n", matrix_c_reference[i][j], matrix_c[i][j], i, j);
                return 0;
            }
        }
    }
    printf("Wyniki poprawne.\n");
}

void multiply_sequential(float matrix_a[N][N], float matrix_b[N][N], float matrix_c_reference[N][N])
{
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < N; k++) {
            for (int j = 0; j < N; j++) {
                matrix_c_reference[i][j] += matrix_a[i][k] * matrix_b[k][j];
            }
        }
    }
}

int main(int argc, char* argv[])
{

    clock_t start, end;
    double used_time;	
	
    float matrix_a[N][N];
    float matrix_b[N][N];
    float matrix_c[N][N];
    float matrix_c_reference[N][N];
    


    initialize_matrices(matrix_a, matrix_b, matrix_c, matrix_c_reference);

    multiply_sequential(matrix_a, matrix_b, matrix_c_reference);


    printf("*********************\n");
    printf("Metoda 3-petlowa\n");
    // Test dla funkcji IKJ
    clear_result_matrix(matrix_c);
    start = clock();
    multiply_matrices_IKJ(matrix_a, matrix_b, matrix_c);
    end= clock();
    used_time= ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Funkcja IJK zajela %f sekund\n", used_time);
    check_results(matrix_c, matrix_c_reference);
    printf("\n");


    printf("*********************\n");
    printf("Metoda 4-petlowa\n");
    // Test dla funkcji 4-pętlowej IKJ
    clear_result_matrix(matrix_c);
    start = clock();
    multiply_matrices_4_loop(matrix_a, matrix_b, matrix_c);
    end= clock();
    used_time= ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Funkcja 4_loops zajela %f sekund\n", used_time);
    check_results(matrix_c, matrix_c_reference);
    printf("\n");
    



    printf("*********************\n");
    printf("Metoda 6-petlowa\n");
    // Test dla funkcji 6-pętlowej  IKJ
    clear_result_matrix(matrix_c);
    start = clock();
    multiply_matrices_6_loop(matrix_a, matrix_b, matrix_c);
    end= clock();
    used_time= ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Funkcja 6_loops zajela %f sekund\n", used_time);
    check_results(matrix_c, matrix_c_reference);
    printf("\n");

    return 0;
}