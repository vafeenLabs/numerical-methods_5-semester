#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

using namespace std;

void printRectMatrix(int N, int L, double** A, double* f)
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

void printRectMatrixWithoutF(int N, int L, double** A)
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

void printRectMatrixRounded(int N, int L, double** A, double* f)
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
double* solveBottomBandMatrixIncorrect(int N, int L, double** A, double* f)
{
    double mult;
    double* cur_bottom = new double[L - 1];
    // Прямой ход метода Гаусса с единственным делением
    for (int i = 0; i < N; i++)
    {
        double diag = A[i][L - 1];
        /*if (diag == 0.0)
        {
            cerr << "Error: diagonal element [" << i << "][" << i << "] equals 0.";
            return nullptr;
        }*/
        f[i] /= diag;
        A[i][L - 1] = 1.0;

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
                //printRectMatrixRounded(N, L, A, f);
            }

            f[i + j] -= f[i] * mult;

            //printRectMatrixRounded(N, L, A, f);
        }
        printRectMatrixRounded(N, L, A, f);
    }

    // Обратный ход
    double* x = new double[N];
    for (int i = N - 1; i >= 0; i--)
    {
        x[i] = f[i];

        for (int k = L - 2; k >= 0; k--)
        {
            x[i] -= A[i][k] * x[k - i];
        }
    }

    delete[] cur_bottom;

    return x;
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

double* generateRandBand(ofstream& file, int N, int L)
{
    srand(time(NULL));
    double* genx = new double[N];
    double** newA = new double* [N];
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

    double* f = new double[N];
    for (int i = N - 1; i >= 0; i--)
    {
        f[i] = 0.0;

        for (int k = i; k < N && k - i < L; k++)
            f[i] += newA[k][L - (k - i) - 1] * genx[k];

        for (int j = 1; j < L && (i - j) >= 0; j++)
            f[i] += newA[i][L - j - 1] * genx[i - j];
    }

    file << N << endl << L << endl;
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

void MultiplyDiag(int N, int L, double** A, double* f, double* x, int k)
{
    double badK = pow(10, (-2) * k);
    for (int i = 0; i < N; i++)
    {
        A[i][L - 1] *= badK;
    }

    for (int i = N - 1; i >= 0; i--)
    {
        f[i] = 0.0;

        for (int k = i; k < N && k - i < L; k++)
            f[i] += A[k][L - (k - i) - 1] * x[k];

        for (int j = 1; j < L && (i - j) >= 0; j++)
            f[i] += A[i][L - j - 1] * x[i - j];
    }
}

void readMatrixFromFile(ifstream& inputFile, int N, int L, double* f, double** A)
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

void test1(string testFileName, int N, int L, int numTests)
{
    double globDiff = 0.0, diff = 0.0, curDiff = 0.0;
    double* f = new double[N];
    double** A = new double* [N];
    for (int k = 0; k < numTests; k++)
    {
        ofstream fileWrite(testFileName + ".txt");
        double* sourceX = generateRandBand(fileWrite, N, L);
        fileWrite.close();

        ifstream fileRead(testFileName + ".txt");
        fileRead >> N >> L;
        readMatrixFromFile(fileRead, N, L, f, A);

        double* newX = solveBottomBandMatrix(N, L, A, f);

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

void test2(string testFileName, int N, int L, int numTests)
{
    cout << "N = " << N << "  L = " << L << endl;

    double diff = 0.0, curDiff = 0.0;
    double* F = new double[N];
    double** A = new double* [N];
    double* globDiffers = new double[3];
    for (int i = 0; i < 3; i++)
        globDiffers[i] = 0.0;
    
    for (int j = 0; j < numTests; j++)
    {
        ofstream fileWrite(testFileName + ".txt");
        double* sourceX = generateRandBand(fileWrite, N, L);
        fileWrite.close();

        for (int k = 2; k < 8; k += 2)
        {
            ifstream fileRead(testFileName + ".txt");
            fileRead >> N >> L;
            readMatrixFromFile(fileRead, N, L, F, A);

            MultiplyDiag(N, L, A, F, sourceX, k);

            double* newX = solveBottomBandMatrix(N, L, A, F);

            diff = 0.0;
            for (int i = 0; i < N; i++)
            {
                curDiff = abs(newX[i] - sourceX[i]);
                diff = max(diff, curDiff);
            }
            globDiffers[k / 2 - 1] += diff;
        }
    }
    for (int i = 0; i < 3; i++)
    cout << "Порядок матрицы = " << N << "x" << N << ", диапазон эл-ов = [-10; 10], k = " << (i + 1) * 2 << ": " << globDiffers[i] / numTests << endl;
}

void testSolveFile(ifstream inputFile)
{
    int N, L;
    inputFile >> N >> L;
    double* f = new double[N];
    double** A = new double* [N];

    readMatrixFromFile(inputFile, N, L, f, A);
    inputFile.close();

    double* genX = new double[N];
    for (int i = 0; i < N; i++)
    {
        genX[i] = (double)rand() / (double)RAND_MAX * 20.0 - 10.0;
    }

    for (int i = N - 1; i >= 0; i--)
    {
        f[i] = 0.0;

        for (int k = i; k < N && k - i < L; k++)
            f[i] += A[k][L - (k - i) - 1] * genX[k];

        for (int j = 1; j < L; j++)
            f[i] += A[i][L - j - 1] * genX[i - j];
    }

    double* newX = solveBottomBandMatrix(N, L, A, f);

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
}

double** cloneMatrix(double** A, int N, int L)
{
    double** clone = new double* [N];
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
        double* f = new double[N];
        double** A = new double* [N];

        readMatrixFromFile(inputFile, N, L, f, A);

        // Закрываем файл
        inputFile.close();

        printRectMatrixRounded(N, L, A, f);
        // Решение системы
        double* x = solveBottomBandMatrix(N, L, A, f);

        // Вывод результата
        cout << "Решение системы:" << endl;
        for (int i = 0; i < N; i++)
        {
            cout << "x[" << i << "] = " << x[i] << endl;
        }

        
        double** newA = cloneMatrix(A, N, L);
        double* genX = new double[N];
        double* newF = new double[N];
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

        double* newX = solveBottomBandMatrix(N, L, newA, newF);

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
        int N = 10;
        int L = N / 10;
        int numTests = 20;
        string testFileName = "testRand";
        cout << "\nВычислительный эксперимент №1:\n\n";
        
        test1(testFileName, N, L, numTests);

        L = floor(sqrt(N));
        test1(testFileName, N, L, numTests);
        
        N = 500;
        L = 50;
        test1(testFileName, N, L, numTests);

        L = floor(sqrt(N));
        test1(testFileName, N, L, numTests);

        cout << "\nВычислительный эксперимент №2:\n\n";
        N = L = 10;
        test1(testFileName, N, L, numTests);
        N = L = 100;
        test1(testFileName, N, L, numTests);


        cout << "\nВычислительный эксперимент №3:\n\n";
        //N = rand() % 100;
        //L = rand() % 15;
        N = 200;
        L = 35;
        test2(testFileName, N, L, numTests);
    }
    
    return 0;
}
