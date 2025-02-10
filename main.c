/*

Projekt mnożenia macierzy Filip Hojan & Igor Swiatek

*/

#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <time.h>
//#include <windows.h>
#include "omp.h"
//#include <cmath>
#include <math.h>
#include <stdlib.h>

//u mnie cache ma 12288 kb
#define N 5000
#define PP 3145728 // 3 MB w bajtach

void multiply_matrices_IKJ(float matrix_a[N][N], float matrix_b[N][N], float matrix_c[N][N])
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

void multiply_matrices_4_loop(float matrix_a[N][N], float matrix_b[N][N], float matrix_c[N][N])
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

void multiply_matrices_6_loop(float matrix_a[N][N], float matrix_b[N][N], float matrix_c[N][N])
{
    int r = floor(sqrt(PP/12.0));
    int ka = floor(N/r);


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

void initialize_matrices(float matrix_a[N][N], float matrix_b[N][N], float matrix_c[N][N])
{
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matrix_a[i][j] = (float)rand() / RAND_MAX;
            matrix_b[i][j] = (float)rand() / RAND_MAX;
            matrix_c[i][j] = 0.0;
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

/*void print_elapsed_time(double start, char* description)
{
    double elapsed;
    elapsed = (double)clock() / CLK_TCK;
    printf("%s \t\t| czas: %8.4f sec \n", description, elapsed - start);
    fprintf(result_file, "%s \t\t| czas : %8.4f sec (%6.4f sec rozdzielczosc pomiaru)\n", description, elapsed-start, 1.0 / CLK_TCK);
}*/

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

void multiply_sequential(float matrix_a[N][N], float matrix_b[N][N], float matrix_c[N][N])
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
    clock_t start, end;
    double used_time;	
	
    float matrix_a[N][N];
    float matrix_b[N][N];
    float matrix_c[N][N];
    float matrix_c_reference[N][N];

    initialize_matrices(matrix_a, matrix_b, matrix_c);
    //print_elapsed_time(clock() / CLK_TCK, "Inicjalizacja macierzy");
    multiply_sequential(matrix_a, matrix_b, matrix_c_reference);

    // Test dla funkcji IKJ
    clear_result_matrix(matrix_c);
    start = clock();
    multiply_matrices_IKJ(matrix_a, matrix_b, matrix_c);
    //print_elapsed_time(start, "Mnozenie macierzy IKJ");
    end= clock();
    used_time= ((double) (end - start)) / CLOCKS_PER_SEC;
    check_results(matrix_c, matrix_c_reference);
    printf("Funkcja IJK zajela %f sekund\n", used_time);

    // Test dla funkcji 4-pętlowej IKJ
    clear_result_matrix(matrix_c);
    start = clock();
    multiply_matrices_4_loop(matrix_a, matrix_b, matrix_c);
    //print_elapsed_time(start, "Mnozenie macierzy 4-petlowe");
    end= clock();
    used_time= ((double) (end - start)) / CLOCKS_PER_SEC;
    check_results(matrix_c, matrix_c_reference);
    printf("Funkcja 4_loops zajela %f sekund\n", used_time);

    // Test dla funkcji 6-pętlowej  IKJ
    clear_result_matrix(matrix_c);
    start = clock();
    multiply_matrices_6_loop(matrix_a, matrix_b, matrix_c);
    //print_elapsed_time(start, "Mnozenie macierzy 6-petlowe");
    end= clock();
    used_time= ((double) (end - start)) / CLOCKS_PER_SEC;
    check_results(matrix_c, matrix_c_reference);
    printf("Funkcja 6_loops zajela %f sekund\n", used_time);

    //fclose(result_file);
    return 0;
}
