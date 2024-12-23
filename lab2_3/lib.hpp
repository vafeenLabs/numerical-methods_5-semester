
// // тут last - верхняя полулента, new - нижняя 

// double *solveBottomBandMatrixNew(int N, int L, double **A, double *f)
// {
//     double mult;
//     double *mult_set = new double[L - 1];

//     for (int i = N - 1; i >= 0; i--)
//     {
//         double diag = A[i][L - 1];

//         f[i] /= diag;
//         A[i][L - 1] = 1.0;

//         for (int k = 1; k < L; k++) // сохраняем множители для операций со строками
//         {
//             mult_set[k - 1] = A[i][L - k - 1];
//             A[i][L - k - 1] /= diag;
//         }

//         for (int j = 0; j < L - 1 && (i - j - 1) >= 0; j++) // вычитание текущей строки из следующих
//         {
//             mult = mult_set[j];
//             for (int k = 1; k < L && (L - k - j - 1) >= 0; k++)
//             {
//                 A[i - j - 1][L - k] -= mult * A[i][L - k - j - 1];
//                 // printMatrixRounded(N, L, A, f);
//             }

//             f[i - j - 1] -= f[i] * mult;

//             // printMatrixRounded(N, L, A, f);
//         }

//         // printMatrixRounded(N, L, A, f);
//     }

//     // Обратный ход
//     double *x = new double[N];
//     for (int i = 0; i < N; i++)
//     {
//         x[i] = f[i];

//         for (int k = 1; k < L && i - k >= 0; k++)
//         {
//             x[i] -= A[i][L - k - 1] * x[i - k];
//         }
//     }

//     delete[] mult_set;

//     return x;
// }

// // Функция для решения СЛАУ с нижнеленточной матрицей
// double *solveBottomBandMatrixLast(int N, int L, double **A, double *f)
// {
//     double mult;
//     double *cur_bottom = new double[L - 1];
//     // Прямой ход метода Гаусса с единственным делением
//     for (int i = 0; i < N; i++)
//     {
//         double diag = A[i][L - 1];
//         A[i][L - 1] = 1.0;
//         f[i] /= diag;

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
//                 // printMatrixRounded(N, L, A, f);
//             }
//             mult = cur_bottom[j - 1];
//             f[i + j] -= f[i] * mult;

//             // printRectMatrixRounded(N, L, A, f);
//         }
//     }

//     // Обратный ход
//     double *x = new double[N];
//     for (int i = N - 1; i >= 0; i--)
//     {
//         x[i] = f[i];

//         for (int k = i + 1; k < N && k - i < L; k++)
//         {
//             x[i] -= A[k][L - (k - i) - 1] * x[k];
//         }
//         // x[i] /= A[i][L - 1];
//     }

//     delete[] cur_bottom;

//     return x;
// }