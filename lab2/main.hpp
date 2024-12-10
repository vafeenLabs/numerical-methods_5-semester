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
    double mult;
    double* mult_set = new double[L - 1];
    
    for (int i = N - 1; i >= 0; i--)
    {
        double diag = A[i][L - 1];

        f[i] /= diag;
        A[i][L - 1] = 1.0;

        for (int k = 1; k < L; k++) // сохраняем множители для операций со строками
        {
            mult_set[k - 1] = A[i][L - k - 1];
            A[i][L - k - 1] /= diag;
        }


        for (int j = 0; j < L - 1 && (i - j - 1) >= 0; j++) // вычитание текущей строки из следующих
        {
            mult = mult_set[j];
            for (int k = 1; k < L && (L - k - j - 1) >= 0; k++)
            {
                A[i - j - 1][L - k] -= mult * A[i][L - k - j - 1];
                //printRectMatrixRounded(N, L, A, f);
            }

            f[i - j - 1] -= f[i] * mult;

            //printRectMatrixRounded(N, L, A, f);
        }
        
        //printRectMatrixRounded(N, L, A, f);
    }

    // Обратный ход
    double* x = new double[N];
    for (int i = 0; i < N; i++)
    {
        x[i] = f[i];

        for (int k = 1; k < L && i - k >= 0; k++)
        {
            x[i] -= A[i][L - k - 1] * x[i - k];
        }
    }

    delete[] mult_set;

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