#include "main.hpp"

double *generateRandBand(ofstream &file, int N, int L)
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

void MultiplyDiag(int N, int L, double **A, double *f, double *x, int k)
{
    double badK = pow(10, (-1) * k);
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

void test1(string testFileName, int N, int L, int numTests)
{
    double globDiff = 0.0, diff = 0.0, curDiff = 0.0;
    double *f = new double[N];
    double **A = new double *[N];
    for (int k = 0; k < numTests; k++)
    {
        ofstream fileWrite(testFileName + ".txt");
        double *sourceX = generateRandBand(fileWrite, N, L);
        if (N == 10)
        {
            cout << "Generated solution:->";
            printArray(sourceX, N);
            cout << '\n';
        }
        fileWrite.close();

        ifstream fileRead(testFileName + ".txt");
        fileRead >> N >> L;
        readMatrixFromFile(fileRead, N, L, f, A);

        double *newX = solveBottomBandMatrix(N, L, A, f);
        if (N == 10)
        {
            cout << "Solved solution:->";
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
    cout << "Matrix size = " << N << "x" << N << ", range = [-10; 10]: " << globDiff << endl;
}

void test2(string testFileName, int N, int L, int numTests)
{
    cout << "N = " << N << "  L = " << L << endl;

    double diff = 0.0, curDiff = 0.0;
    double *F = new double[N];
    double **A = new double *[N];
    double *globDiffers = new double[3];
    for (int i = 0; i < 3; i++)
        globDiffers[i] = 0.0;

    for (int j = 0; j < numTests; j++)
    {
        ofstream fileWrite(testFileName + ".txt");
        double *sourceX = generateRandBand(fileWrite, N, L);
        fileWrite.close();

        for (int k = 2; k < 8; k += 2)
        {
            ifstream fileRead(testFileName + ".txt");
            fileRead >> N >> L;
            readMatrixFromFile(fileRead, N, L, F, A);
            // Умножаем элементы на 10^-2к
            MultiplyDiag(N, L, A, F, sourceX, k);

            double *newX = solveBottomBandMatrix(N, L, A, F);

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
        cout << "Matrix size = " << N << "x" << N << ", range = [-10; 10], k = " << (i + 1) * 2 << ": " << globDiffers[i] / numTests << endl;
}


int main()
{
    setlocale(LC_ALL, "rus");
    srand(time(NULL));

    std::filesystem::create_directory(testDir);
    int n;
    cout << "1: Matrix from file\t2: Tests\t3: Exit\n->: ";
    cin >> n;

    if (n == 1)
    {
        cout << "Enter fileName: ";
        string fileName;
        cin >> fileName;
        fileName += ".txt";
        ifstream inputFile(fileName);
        if (!inputFile)
        {
            cerr << "Error openning file " << fileName << ".\n";
            return 1;
        }

        int N, L;
        inputFile >> N >> L;
        double *f = new double[N];
        double **A = new double *[N];

        readMatrixFromFile(inputFile, N, L, f, A);

        inputFile.close();

        printMatrixRounded(N, L, A, f);
        // Решение системы
        double *x = solveBottomBandMatrix(N, L, A, f);

        // Вывод результата
        cout << "Solution:" << endl;
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

        double *newX = solveBottomBandMatrix(N, L, newA, newF);

        double diff = 0.0, curDiff = 0.0;
        for (int i = 0; i < N; i++)
        {
            for (int i = 0; i < N; i++)
            {
                curDiff = abs(newX[i] - genX[i]);
                diff = max(diff, curDiff);
            }
        }
        cout << "Avarage difference for system = " << diff << endl;

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
        int numTests = 15;
        string testFileName = testDir + "testRand";
        cout << "\nTest 1:\n\n";
        // part1.1
        int N = 10;
        int L = N / 10; // 1/10
        test1(testFileName, N, L, numTests);

        L = floor(sqrt(N)); // 1/L
        test1(testFileName, N, L, numTests);

        // part 1.2
        N = 400;
        L = 40; // 1/10
        test1(testFileName, N, L, numTests);

        L = floor(sqrt(N)); // 1/L
        test1(testFileName, N, L, numTests);

        cout << "\nTest 2:\n\n";
        N = L = 10;
        test1(testFileName, N, L, numTests);
        N = L = 100;
        test1(testFileName, N, L, numTests);

        cout << "\nTest 3:\n\n";
        N = 100;
        L = 50;
        test2(testFileName, N, L, numTests);
    }
    if (n == 3)
    {
        exit(0);
    }
    return 0;
}
