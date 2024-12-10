#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <filesystem>
using namespace std;

string testDir = "tests/";

void printArray(double *array, int size)
{
    for (int i = 0; i < size; ++i)
    {
        cout << array[i] << ' ';
    }
}


void printMatrix(int N, int L, double **A, double *f)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < L; j++)
        {
            cout << A[i][j] << '\t';
        }
        cout << "|\t" << f[i] << endl;
    }
    for (int i = 0; i < 40; i++)
        cout << '-';
    cout << endl;
}

void printMatrixWithoutF(int N, int L, double **A)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < L; j++)
        {
            cout << round(A[i][j] * 100.0) / 100.0 << '\t';
        }
        cout << endl;
    }
    for (int i = 0; i < 40; i++)
        cout << '-';
    cout << endl;
}

void printMatrixRounded(int N, int L, double **A, double *f)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < L; j++)
        {
            cout << round(A[i][j] * 100.0) / 100.0 << '\t';
        }
        cout << "|\t" << f[i] << endl;
    }
    for (int i = 0; i < 41; i++)
        cout << '-';
    cout << endl;
}


double* solveBottomBandMatrix(int N, int L, double** A, double* f)
{
    // Объявляем переменную для множителя и массив для хранения множителей текущей строки.
    double mult;
    double* mult_set = new double[L - 1];
    
    // Прямой ход метода Гаусса (снизу вверх)
    for (int i = N - 1; i >= 0; i--) // Проходим по строкам от последней к первой.
    {
        // Извлекаем диагональный элемент текущей строки.
        double diag = A[i][L - 1];

        // Нормализуем элемент в правой части уравнения.
        f[i] /= diag;
        A[i][L - 1] = 1.0; // Диагональный элемент приводим к 1.

        // Сохраняем множители для вычитания строк.
        for (int k = 1; k < L; k++) 
        {
            mult_set[k - 1] = A[i][L - k - 1]; // Сохраняем элементы слева от диагонали.
            A[i][L - k - 1] /= diag; // Делим элементы слева от диагонали на диагональный элемент.
        }

        // Вычитание текущей строки из строк выше, чтобы занулить соответствующие элементы.
        for (int j = 0; j < L - 1 && (i - j - 1) >= 0; j++) 
        {
            mult = mult_set[j]; // Берем множитель из сохраненного набора.
            for (int k = 1; k < L && (L - k - j - 1) >= 0; k++) 
            {
                // Обновляем элементы матрицы выше текущей строки.
                A[i - j - 1][L - k] -= mult * A[i][L - k - j - 1];
            }

            // Обновляем правую часть для строки выше.
            f[i - j - 1] -= f[i] * mult;
        }
    }

    // Обратный ход метода Гаусса.
    double* x = new double[N]; // Результирующий массив решений.
    for (int i = 0; i < N; i++)
    {
        x[i] = f[i]; // Инициализируем x текущим значением из правой части.

        // Учитываем вклад элементов, расположенных выше по диагонали.
        for (int k = 1; k < L && i - k >= 0; k++)
        {
            x[i] -= A[i][L - k - 1] * x[i - k];
        }
    }

    // Освобождаем память для массива множителей.
    delete[] mult_set;

    return x; // Возвращаем массив решений.
}


double **cloneMatrix(double **A, int N, int L)
{
    double **clone = new double *[N];
    for (int i = 0; i < N; i++)
        clone[i] = new double[L];

    for (int i = 0; i < N; i++)
        for (int j = 0; j < L; j++)
            clone[i][j] = A[i][j];

    return clone;
}