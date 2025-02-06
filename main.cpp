#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <time.h>
#include <windows.h>
#include "omp.h"
#include <cmath>

static const int ROWS = 1000;     // liczba wierszy macierzy  2048/3072/4096
static const int COLUMNS = 1000;  // liczba kolumn macierzy   2048/3072/4096
float matrix_a[ROWS][COLUMNS];    // lewy operand 
float matrix_b[ROWS][COLUMNS];    // prawy operand
float matrix_c[ROWS][COLUMNS];    // wynik
float matrix_c_reference[ROWS][COLUMNS];  // wynik referencyjny dla porównania

FILE *result_file;

void initialize_matrices();
void clear_result_matrix();
void print_elapsed_time(double start, char* description);
void check_results();

void multiply_matrices_IJK();
void multiply_matrices_IKJ();
void multiply_matrices_JIK();
void multiply_matrices_JKI();
void multiply_matrices_KIJ();
void multiply_matrices_KJI();
void multiply_matrices_4_loop();
void multiply_matrices_6_loop();

void multiply_matrices_IJK() 
{ 
    #pragma omp parallel for
    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLUMNS; j++) {
            float sum = 0.0;
            for (int k = 0; k < COLUMNS; k++) {
                sum += matrix_a[i][k] * matrix_b[k][j];
            }
            matrix_c[i][j] = sum;
        }
    }
}

void multiply_matrices_IKJ()
{
    #pragma omp parallel for
    for (int i = 0; i < ROWS; i++) {
        for (int k = 0; k < COLUMNS; k++) {
            float temp = matrix_a[i][k];
            for (int j = 0; j < COLUMNS; j++) {
                matrix_c[i][j] += temp * matrix_b[k][j];
            }
        }
    }
}

void multiply_matrices_JIK()
{
    #pragma omp parallel for
    for (int j = 0; j < COLUMNS; j++) {
        for (int i = 0; i < ROWS; i++) {
            float sum = 0.0;
            for (int k = 0; k < COLUMNS; k++) {
                sum += matrix_a[i][k] * matrix_b[k][j];
            }
            matrix_c[i][j] = sum;
        }
    }
}

void multiply_matrices_JKI()
{
    #pragma omp parallel for
    for (int j = 0; j < COLUMNS; j++) {
        for (int k = 0; k < COLUMNS; k++) {
            float temp = matrix_b[k][j];
            for (int i = 0; i < ROWS; i++) {
                matrix_c[i][j] += matrix_a[i][k] * temp;
            }
        }
    }
}

void multiply_matrices_KIJ()
{
    #pragma omp parallel for
    for (int k = 0; k < COLUMNS; k++) {
        for (int i = 0; i < ROWS; i++) {
            float temp = matrix_a[i][k];
            for (int j = 0; j < COLUMNS; j++) {
                matrix_c[i][j] += temp * matrix_b[k][j];
            }
        }
    }
}

void multiply_matrices_KJI()
{
    #pragma omp parallel for
    for (int k = 0; k < COLUMNS; k++) {
        for (int j = 0; j < COLUMNS; j++) {
            float temp = matrix_b[k][j];
            for (int i = 0; i < ROWS; i++) {
                matrix_c[i][j] += matrix_a[i][k] * temp;
            }
        }
    }
}

void multiply_matrices_4_loop()
{
    int r = 10;
    #pragma omp parallel for
    for (int k = 0; k < ROWS; k += r) {
        for (int ii = 0; ii < ROWS; ii++) {
            for (int kk = k; kk < k + r; kk++) {
                for (int jj = 0; jj < COLUMNS; jj++) {
                    matrix_c[ii][jj] += matrix_a[ii][kk] * matrix_b[kk][jj];
                }
            }
        }
    }
}

void multiply_matrices_6_loop()
{
    int r = 10;
    #pragma omp parallel for
    for (int i = 0; i < ROWS; i += r) {
        for (int j = 0; j < COLUMNS; j += r) {
            for (int k = 0; k < COLUMNS; k += r) {
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

void initialize_matrices()
{
    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLUMNS; j++) {
            matrix_a[i][j] = (float)rand() / RAND_MAX;
            matrix_b[i][j] = (float)rand() / RAND_MAX;
            matrix_c[i][j] = 0.0;
        }
    }
}

void clear_result_matrix()
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

void check_results()
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

int main(int argc, char* argv[])
{
    if ((result_file = fopen("classic.txt", "a")) == NULL) {
        fprintf(stderr, "Nie można otworzyć pliku wyniku \n");
        perror("classic");
        return(EXIT_FAILURE);
    }

    double start;

    // **Test dla funkcji IJK**
    initialize_matrices();
    print_elapsed_time(clock() / CLK_TCK, "Inicjalizacja macierzy");

    clear_result_matrix();
    start = clock() / CLK_TCK;
    multiply_matrices_IJK();
    print_elapsed_time(start, "Mnozenie macierzy IJK");

    // Skopiowanie wyniku referencyjnego do matrix_c_reference
    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLUMNS; j++) {
            matrix_c_reference[i][j] = matrix_c[i][j];
        }
    }

    // **Test dla funkcji IKJ**
    clear_result_matrix();
    start = clock() / CLK_TCK;
    multiply_matrices_IKJ();
    print_elapsed_time(start, "Mnozenie macierzy IKJ");
    check_results();

    // **Test dla funkcji JIK**
    clear_result_matrix();
    start = clock() / CLK_TCK;
    multiply_matrices_JIK();
    print_elapsed_time(start, "Mnozenie macierzy JIK");
    check_results();

    // **Test dla funkcji JKI**
    clear_result_matrix();
    start = clock() / CLK_TCK;
    multiply_matrices_JKI();
    print_elapsed_time(start, "Mnozenie macierzy JKI");
    check_results();

    // **Test dla funkcji KIJ**
    clear_result_matrix();
    start = clock() / CLK_TCK;
    multiply_matrices_KIJ();
    print_elapsed_time(start, "Mnozenie macierzy KIJ");
    check_results();

    // **Test dla funkcji KJI**
    clear_result_matrix();
    start = clock() / CLK_TCK;
    multiply_matrices_KJI();
    print_elapsed_time(start, "Mnozenie macierzy KJI");
    check_results();

    // **Test dla funkcji 4-pętlowej**
    clear_result_matrix();
    start = clock() / CLK_TCK;
    multiply_matrices_4_loop();
    print_elapsed_time(start, "Mnozenie macierzy 4-petlowe");
    check_results();

    // **Test dla funkcji 6-pętlowej**
    clear_result_matrix();
    start = clock() / CLK_TCK;
    multiply_matrices_6_loop();
    print_elapsed_time(start, "Mnozenie macierzy 6-petlowe");
    check_results();

    fclose(result_file);
    return 0;
}

