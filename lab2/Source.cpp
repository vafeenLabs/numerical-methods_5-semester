#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

void printRectMatrix(int N, int L, double** A, double* f)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < L; j++)
        {
            cout << A[i][j] << '\t';
        }
        cout << f[i] << endl;
    }
    for (int i = 0; i < 40; i++)
        cout << '-';
    cout << endl;
}

void printRectMatrixRounded(int N, int L, double** A, double* f)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < L; j++)
        {
            cout << round(A[i][j] * 100.0) / 100.0 << '\t';
        }
        cout << f[i] << endl;
    }
    for (int i = 0; i < 40; i++)
        cout << '-';
    cout << endl;
}

// решение с нижнеленточной матрицей
double* solveBottomBandMatrix(int N, int L, double** A, double* f)
{
    double mult;
    double* first_column_backup = new double[L - 1];
    // Прямой ход метода Гаусса с единственным делением
    for (int i = 0; i < N; i++)
    {
        double diag = A[i][L - 1];
        A[i][L - 1] = 1.0;
        f[i] /= diag;

        for (int k = 1; k < L && i + k < N; k++) // сохраняем множители для операций со строками
        {
            first_column_backup[k - 1] = A[i + k][L - k - 1];
        }

        for (int j = 1; j < L && i + j < N; j++)
        {
            A[i + j][L - j - 1] /= diag;
        }

        for (int j = 1; j < L && i + j < N; j++) // вычитание текущей строки из следующих
        {
            mult = first_column_backup[j - 1];
            for (int k = j; k < L && i + k < N; k++)
            {
                A[i + k][L - (k - j) - 1] -= mult * A[i + k][L - k - 1];
                //printRectMatrixRounded(N, L, A, f);
            }
            f[i + j] -= f[i] * mult;

            printRectMatrixRounded(N, L, A, f);
        }

    }

    // Обратный ход
    double* x = new double[N];
    for (int i = N - 1; i >= 0; i--)
    {
        x[i] = f[i];

        for (int k = i + 1; k < N && k - i < L; k++)
        {
            x[i] -= A[k][L - (k - i) - 1] * x[k];
        }
        //x[i] /= A[i][L - 1];
    }

    delete[] first_column_backup;

    return x;
}


int main()
{
    setlocale(LC_ALL, "rus");
    ifstream inputFile("test.txt");
    if (!inputFile)
    {
        cerr << "Ошибка: не удалось открыть файл input.txt" << endl;
        return 1;
    }

    // Считываем размерность и половину ширины ленты
    int N, L;
    inputFile >> N >> L;

    // Создаем динамический массив для ленты матрицы A и для вектора правой части f
    double* f = new double[N];
    double** A = new double* [N];
    for (int i = 0; i < N; i++)
    {
        A[i] = new double[L];
        for (int j = 0; j < L; j++)
        {
            inputFile >> A[i][j];
        }
        inputFile >> f[i];
    }

    // Закрываем файл
    inputFile.close();

    printRectMatrix(N, L, A, f);
    // Решение системы
    double* x = solveBottomBandMatrix(N, L, A, f);

    // Вывод результата
    cout << "Решение системы:" << endl;
    for (int i = 0; i < N; i++)
    {
        cout << "x[" << i << "] = " << x[i] << endl;
    }

    // Освобождение памяти
    for (int i = 0; i < N; i++)
    {
        delete[] A[i];
    }
    delete[] A;
    delete[] f;
    delete[] x;
    return 0;
}
