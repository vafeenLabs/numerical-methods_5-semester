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


// Функция для решения СЛАУ с нижнеленточной матрицей
double *solveBottomBandMatrix(int N, int L, double **A, double *f)
{
    double mult;
    double *cur_bottom = new double[L - 1];
    // Прямой ход метода Гаусса с единственным делением
    for (int i = 0; i < N; i++)
    {
        double diag = A[i][L - 1];
        A[i][L - 1] = 1.0;
        f[i] /= diag;

        for (int k = 1; k < L && i + k < N; k++) // сохраняем множители для операций со строками
        {
            cur_bottom[k - 1] = A[i + k][L - k - 1];
        }

        for (int j = 1; j < L && i + j < N; j++)
        {
            A[i + j][L - j - 1] /= diag;
        }

        for (int j = 1; j < L && i + j < N; j++) // вычитание текущей строки из следующих
        {
            mult = cur_bottom[j - 1];
            for (int k = j; k < L && i + k < N; k++)
            {
                A[i + k][L - (k - j) - 1] -= mult * A[i + k][L - k - 1];
                // printMatrixRounded(N, L, A, f);
            }
            mult = cur_bottom[j - 1];
            f[i + j] -= f[i] * mult;

            // printRectMatrixRounded(N, L, A, f);
        }
    }

    // Обратный ход
    double *x = new double[N];
    for (int i = N - 1; i >= 0; i--)
    {
        x[i] = f[i];

        for (int k = i + 1; k < N && k - i < L; k++)
        {
            x[i] -= A[k][L - (k - i) - 1] * x[k];
        }
        // x[i] /= A[i][L - 1];
    }

    delete[] cur_bottom;

    return x;
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