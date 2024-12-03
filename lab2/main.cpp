#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

using namespace std;

// вывод матрицы 
void printRectMatrixRounded(int N, int L, double **A, double *f)
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

// решение системы 
double *solveBottomMatrix(int N, int L, double **A, double *f)
{
    double mult;
    double *mult_set = new double[L - 1];

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
                // printRectMatrixRounded(N, L, A, f);
            }

            f[i - j - 1] -= f[i] * mult;

            // printRectMatrixRounded(N, L, A, f);
        }

        // printRectMatrixRounded(N, L, A, f);
    }

    // Обратный ход
    double *x = new double[N];
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

// генерация рандомной матрицы
double *generateRandMatrix(ofstream &file, int N, int L)
{
    double *genx = new double[N];
    double **newA = new double *[N];
    for (int i = 0; i < N; i++)
        newA[i] = new double[L];

    for (int i = 0; i < N; i++)
    {
        genx[i] = (double)rand() / (double)RAND_MAX * 20.0 - 10.0;
    }

    for (int i = 0; i < (L - 1); i++)
    {
        for (int j = 0; j < (L - i) - 1; j++)
            newA[i][j] = 0.0;

        for (int j = (L - i) - 1; j < L; j++)
            newA[i][j] = (double)rand() / (double)RAND_MAX * 20.0 - 10.0;
    }

    for (int i = L - 1; i < N; i++)
    {
        for (int j = 0; j < L; j++)
        {
            newA[i][j] = (double)rand() / (double)RAND_MAX * 20.0 - 10.0;
        }
    }

    double *f = new double[N];
    for (int i = N - 1; i >= 0; i--)
    {
        f[i] = 0.0;

        for (int k = i; k < N && k - i < L; k++)
            f[i] += newA[k][L - (k - i) - 1] * genx[k];

        for (int j = 1; j < L && (i - j) >= 0; j++)
            f[i] += newA[i][L - j - 1] * genx[i - j];
    }

    file << N << endl
         << L << endl;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < L; j++)
        {
            file << newA[i][j] << ' ';
        }
        file << f[i] << endl;
    }

    return genx;
}

// чтение из файла
void readMatrixFromFile(ifstream &inputFile, int N, int L, double *f, double **A)
{
    for (int i = 0; i < N; i++)
    {
        A[i] = new double[L];
        for (int j = 0; j < L; j++)
        {
            inputFile >> A[i][j];
        }
        inputFile >> f[i];
    }
}

// вывод динамического массива 
void printArray(double *array, int size)
{
    for (int i = 0; i < size; ++i)
    {
        cout << array[i] << ' ';
    }
}

// тесты 
void test(string testFileName, int N, int L, int numTests)
{

    double globDiff = 0.0, diff = 0.0, curDiff = 0.0;
    double *f = new double[N];
    double **A = new double *[N];
    for (int k = 0; k < numTests; k++)
    {
        ofstream fileWrite(testFileName + ".txt");
        double *sourceX = generateRandMatrix(fileWrite, N, L);
        if (N == 10)
        {
            cout << "Сгенерированные иксы:->";
            printArray(sourceX, N);
            cout << '\n';
        }
        fileWrite.close();
        ifstream fileRead(testFileName + ".txt");
        fileRead >> N >> L;
        readMatrixFromFile(fileRead, N, L, f, A);

        double *newX = solveBottomMatrix(N, L, A, f);
        if (N == 10)
        {
            cout << "Полученные иксы:->";
            printArray(newX, N);
            cout << "\n\n";
        }

        diff = 0.0;
        for (int i = 0; i < N; i++)
        {
            curDiff = abs(newX[i] - sourceX[i]);
            diff = max(diff, curDiff);
        }
        globDiff += diff;
    }
    cout << "N = " << N << "  L = " << L << endl;
    globDiff /= numTests;
    cout << "Порядок матрицы = " << N << "x" << N << ", диапазон эл-ов = [-10; 10]: " << globDiff << endl;
}

// бэкап матрицы 
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

int main()
{
    setlocale(LC_ALL, "rus");
    srand(time(NULL));
    int n;
    cout << "1: Матрица из файла\t2: Сгенерировать тесты\nВыбор: ";
    cin >> n;

    if (n == 1)
    {
        cout << "Введите имя файла: ";
        string fileName;
        cin >> fileName;
        fileName += ".txt";
        ifstream inputFile(fileName);
        if (!inputFile)
        {
            cerr << "Ошибка: не удалось открыть файл " << fileName << ".\n";
            return 1;
        }

        int N, L;
        inputFile >> N >> L;
        double *f = new double[N];
        double **A = new double *[N];

        readMatrixFromFile(inputFile, N, L, f, A);

        // Закрываем файл
        inputFile.close();

        printRectMatrixRounded(N, L, A, f);
        // Решение системы
        double *x = solveBottomMatrix(N, L, A, f);

        // Вывод результата
        cout << "Решение системы:" << endl;
        for (int i = 0; i < N; i++)
        {
            cout << "x[" << i << "] = " << x[i] << endl;
        }

        // измерение погрешности
        double **newA = cloneMatrix(A, N, L);
        double *genX = new double[N];
        double *newF = new double[N];
        for (int i = 0; i < N; i++)
        {
            genX[i] = (double)rand() / (double)RAND_MAX * 20.0 - 10.0;
        }

        for (int i = N - 1; i >= 0; i--)
        {
            newF[i] = 0.0;

            for (int k = i; k < N && k - i < L; k++)
                newF[i] += A[k][L - (k - i) - 1] * genX[k];

            for (int j = 1; j < L; j++)
                newF[i] += A[i][L - j - 1] * genX[i - j];
        }

        double *newX = solveBottomMatrix(N, L, newA, newF);

        double diff = 0.0, curDiff = 0.0;
        for (int i = 0; i < N; i++)
        {
            for (int i = 0; i < N; i++)
            {
                curDiff = abs(newX[i] - genX[i]);
                diff = max(diff, curDiff);
            }
        }
        cout << "Средняя относительная погрешность решения = " << diff << endl;

        // Освобождение памяти
        for (int i = 0; i < N; i++)
        {
            delete[] A[i];
        }
        delete[] A;
        for (int i = 0; i < N; i++)
            delete[] newA[i];
        delete[] newA;
        delete[] newF;
        delete[] f;
        delete[] x;

        return 0;
    }

    if (n == 2)
    {
        srand(time(NULL));

        int numTests = 20;
        string testFileName = "testMatrix";

        int N = 10;
        int L = 2;

        test(testFileName, N, L, numTests);

        N = 10;
        L = 5;
        test(testFileName, N, L, numTests);

        N = 50;
        L = 10;
        test(testFileName, N, L, numTests);

        N = 50;
        L = 20;
        test(testFileName, N, L, numTests);
    }

    return 0;
}
