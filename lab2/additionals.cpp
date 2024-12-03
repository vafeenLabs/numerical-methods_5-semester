// #include <iostream>
// #include <fstream>
// #include <cmath>
// #include <algorithm>

// using namespace std;

// void printRectMatrixWithoutF(int N, int L, double **A)
// {
//     for (int i = 0; i < N; i++)
//     {
//         for (int j = 0; j < L; j++)
//         {
//             cout << round(A[i][j] * 100.0) / 100.0 << '\t';
//         }
//         cout << endl;
//     }
//     for (int i = 0; i < 40; i++)
//         cout << '-';
//     cout << endl;
// }


// // Функция для решения СЛАУ с нижнеленточной матрицей
// double *solveBottomBandMatrixIncorrect(int N, int L, double **A, double *f)
// {
//     double mult;
//     double *cur_bottom = new double[L - 1];
//     // Прямой ход метода Гаусса с единственным делением
//     for (int i = 0; i < N; i++)
//     {
//         double diag = A[i][L - 1];
//         /*if (diag == 0.0)
//         {
//             cerr << "Error: diagonal element [" << i << "][" << i << "] equals 0.";
//             return nullptr;
//         }*/
//         f[i] /= diag;
//         A[i][L - 1] = 1.0;

//         for (int k = 1; k < L && i + k < N; k++) // сохраняем множители для операций со строками
//         {
//             cur_bottom[k - 1] = A[i + k][L - k - 1];
//         }

//         for (int j = 1; j < L && i + j < N; j++)
//         {
//             A[i + j][L - j - 1] /= diag;
//         }

//         for (int j = 1; j < L && i + j < N; j++) // вычитание текущей строки из следующих
//         {
//             mult = cur_bottom[j - 1];
//             for (int k = j; k < L && i + k < N; k++)
//             {
//                 A[i + k][L - (k - j) - 1] -= mult * A[i + k][L - k - 1];
//                 // printRectMatrixRounded(N, L, A, f);
//             }

//             f[i + j] -= f[i] * mult;

//             // printRectMatrixRounded(N, L, A, f);
//         }
//         printRectMatrixRounded(N, L, A, f);
//     }

//     // Обратный ход
//     double *x = new double[N];
//     for (int i = N - 1; i >= 0; i--)
//     {
//         x[i] = f[i];

//         for (int k = L - 2; k >= 0; k--)
//         {
//             x[i] -= A[i][k] * x[k - i];
//         }
//     }

//     delete[] cur_bottom;

//     return x;
// }


// void MultiplyDiag(int N, int L, double **A, double *f, double *x, int k)
// {
//     double badK = pow(10, (-2) * k);
//     for (int i = 0; i < N; i++)
//     {
//         A[i][L - 1] *= badK;
//     }

//     for (int i = N - 1; i >= 0; i--)
//     {
//         f[i] = 0.0;

//         for (int k = i; k < N && k - i < L; k++)
//             f[i] += A[k][L - (k - i) - 1] * x[k];

//         for (int j = 1; j < L && (i - j) >= 0; j++)
//             f[i] += A[i][L - j - 1] * x[i - j];
//     }
// }



// void testSolveFile(ifstream inputFile)
// {
//     int N, L;
//     inputFile >> N >> L;
//     double *f = new double[N];
//     double **A = new double *[N];

//     readMatrixFromFile(inputFile, N, L, f, A);
//     inputFile.close();

//     double *genX = new double[N];
//     cout << "Сгенерированные иксы:->";
//     for (int i = 0; i < N; i++)
//     {
//         genX[i] = (double)rand() / (double)RAND_MAX * 20.0 - 10.0;
//         cout << genX[i] << ' ';
//     }
//     cout << "\n";

//     for (int i = N - 1; i >= 0; i--)
//     {
//         f[i] = 0.0;

//         for (int k = i; k < N && k - i < L; k++)
//             f[i] += A[k][L - (k - i) - 1] * genX[k];

//         for (int j = 1; j < L; j++)
//             f[i] += A[i][L - j - 1] * genX[i - j];
//     }

//     double *newX = solveBottomMatrix(N, L, A, f);
//     cout << "Полученные иксы:\n->";
//     printArray(newX, N);
//     cout << "\n";
//     double diff = 0.0, curDiff = 0.0;
//     for (int i = 0; i < N; i++)
//     {
//         for (int i = 0; i < N; i++)
//         {
//             curDiff = abs(newX[i] - genX[i]);
//             diff = max(diff, curDiff);
//         }
//     }
//     cout << "Средняя относительная погрешность решения = " << diff << endl;
// }


// void printRectMatrix(int N, int L, double **A, double *f)
// {
//     for (int i = 0; i < N; i++)
//     {
//         for (int j = 0; j < L; j++)
//         {
//             cout << A[i][j] << '\t';
//         }
//         cout << "|\t" << f[i] << endl;
//     }
//     for (int i = 0; i < 40; i++)
//         cout << '-';
//     cout << endl;
// }