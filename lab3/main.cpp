#include "main.hpp"
#include <iostream>
#include <cmath>
#include <random>
#include <ctime>

using namespace std;
// Генерация случайного нормализованного вектора
void generate_random_vector(long double* omega, int N)
{
    long double norm = 0.0;
    for (int i = 0; i < N; ++i)
    {
        omega[i] = (rand() / static_cast<long double>(RAND_MAX));
        norm += omega[i] * omega[i];
    }

    norm = sqrt(norm);
    for (int i = 0; i < N; ++i)
    {
        omega[i] /= norm;
    }
}

// Умножение матрицы на вектор
void mat_vec_mul(long double** A, long double* x, long double* result, int N)
{
    for (int i = 0; i < N; ++i)
    {
        result[i] = 0.0;
        for (int j = 0; j < N; ++j)
        {
            result[i] += A[i][j] * x[j];
        }
    }
}

// Генерация симметричной матрицы с заданными собственными значениями
void generate_symmetric_matrix(long double** A, long double** H, long double* eigenvalues, int N)
{
    long double* omega = new long double[N];
    generate_random_vector(omega, N);

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            H[i][j] = (i == j ? 1.0 : 0.0) - 2 * omega[i] * omega[j];
        }
    }

    long double** Amid = new long double* [N];
    for (int i = 0; i < N; ++i)
    {
        Amid[i] = new long double[N];
    }

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            Amid[i][j] = H[i][j] * eigenvalues[j];
        }
    }

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            A[i][j] = 0.0;
            for (int k = 0; k < N; ++k)
            {
                A[i][j] += Amid[i][k] * H[k][j];
            }
        }
    }

    for (int i = 0; i < N; ++i)
    {
        delete[] Amid[i];
    }
    delete[] Amid;
    delete[] omega;
}


void generate_orth_vector(long double* input_vec, int N, long double* out_vec)
{
    out_vec[0] = (-1.0) * input_vec[1];
    out_vec[1] = input_vec[0];
    for (int i = 2; i < N; i++)
        out_vec[i] = 0.0;
}

long double dot_product(long double* a, long double* b, int N)
{
    long double dot = 0.0;
    for (int i = 0; i < N; i++)
        dot += a[i] * b[i];
    return dot;
}

long double vec_length(long double* vec, int N)
{
    long double norm = 0.0;
    for (int i = 0; i < N; i++)
        norm += vec[i] * vec[i];

    return sqrt(norm);
}

long double angle_vectors(long double* a, long double* b, int N)
{
    long double dot = dot_product(a, b, N);
    long double len_a = vec_length(a, N);
    long double len_b = vec_length(b, N);

    long double cos = dot / (len_a * len_b);


    if (cos > 1.0) cos = 1.0;
    if (cos < -1.0) cos = -1.0;

    return fabs(acos(fabs(cos)));
}

// Метод обратных итераций с исчерпыванием
void straight_iteration_exhaust(long double** A, int N, long double lambda_1, long double* x_1, long double lambda_2, long double* x_2, long double& lambda_3, long double* x_3, long double epsilon, int M, int& K, long double& r, long double& avg_vec, long double lambda_true, long double* x_true, long double& avg_lambda)
{
    long double* v = new long double[N];
    long double** temp = new long double* [N];
    long double* prev_v = new long double[N];
    long double* prev_x_3 = new long double[N];
    for (int i = 0; i < N; i++)
    {
        temp[i] = new long double[N];
        prev_x_3[i] = x_2[i];
    }

    for (int i = 0; i < N; ++i)
    {
        generate_orth_vector(x_2, N, x_3);
        //x_3[i] = x_2[i];
    }

    //A(1) = A - lambda_n * x_n * x_n^t 
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            temp[i][j] = A[i][j] - lambda_1 * x_1[i] * x_1[j];
        }
    }

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            temp[i][j] = temp[i][j] - lambda_2 * x_2[i] * x_2[j];
        }
    }

    long double alpha_k = 0.0, prev_alpha = 0.0, vecDiff = 10.0, norm;
    K = 0;
    for (int k = 0; k < M; ++k)
    {
        ++K;
        // Нормализация вектора
        norm = 0.0;
        for (int i = 0; i < N; ++i)
        {
            norm += x_3[i] * x_3[i];
        }
        norm = sqrt(norm);
        for (int i = 0; i < N; ++i)
        {
            if (K >= 2)
                prev_v[i] = v[i];
            v[i] = x_3[i] / norm;
        }

        //x_k+1 = A(1) * v(k)
        mat_vec_mul(temp, v, x_3, N);


        // Вычисление alpha_k
        prev_alpha = alpha_k;
        alpha_k = 0.0;
        for (int i = 0; i < N; ++i)
        {
            alpha_k += v[i] * x_3[i];
        }
        if (k >= 2)
            vecDiff = angle_vectors(v, prev_v, N);

        //cout << fabs(alpha_k - prev_alpha) << "  " << vecDiff << endl;
        if (fabs(alpha_k - prev_alpha) < epsilon && vecDiff < epsilon && fabs(alpha_k - lambda_true) < epsilon)
        {
            break;
        }
    }
    avg_vec += vecDiff;
    avg_lambda += fabs(alpha_k - prev_alpha);
    //cout << vecDiff << endl;

    lambda_3 = alpha_k;
    for (int i = 0; i < N; i++)
        x_3[i] = v[i];


    // Вычисление меры точности r
    long double* Ax = new long double[N];
    long double* lambda_x = new long double[N];
    mat_vec_mul(A, x_3, Ax, N);
    for (int i = 0; i < N; ++i)
    {
        lambda_x[i] = lambda_3 * x_3[i];
    }
    r = 0.0;
    for (int i = 0; i < N; ++i)
    {
        r = max(r, fabs(Ax[i] - lambda_x[i]));
    }

    for (int i = 0; i < N; i++)
        delete[] temp[i];
    delete[] v;
    delete[] prev_v;
    delete[] temp;
    delete[] prev_x_3;
    delete[] Ax;
    delete[] lambda_x;
}

void sort_array_abs(long double* arr, int size)
{
    for (int i = 0; i < size - 1; ++i)
    {
        for (int j = 0; j < size - i - 1; ++j)
        {
            if (std::fabs(arr[j]) > std::fabs(arr[j + 1]))
            {
                // Меняем элементы местами
                long double temp = arr[j];
                arr[j] = arr[j + 1];
                arr[j + 1] = temp;
            }
        }
    }
}

// Функция для генерации массива случайных чисел
void generate_random_eigenvalues(long double* eigenvalues, int N, long double lower_bound, long double upper_bound)
{
    for (int i = 0; i < N; ++i)
    {
        eigenvalues[i] = lower_bound + (rand() / static_cast<long double>(RAND_MAX)) * (upper_bound - lower_bound);
    }
    sort_array_abs(eigenvalues, N);
}

void make_test(int N, int numTests, int M, long double lower_bound, long double upper_bound, long double epsilon)
{
    long double avg_diff_lambda = 0.0, avg_diff_vec = 0.0, avg_iterations = 0.0, avg_r = 0.0;
    for (int i = 0; i < numTests; i++)
    {
        long double* eigenvalues = new long double[N];
        generate_random_eigenvalues(eigenvalues, N, lower_bound, upper_bound);
        long double** A = new long double* [N];
        long double** H = new long double* [N];
        for (int i = 0; i < N; ++i)
        {
            A[i] = new long double[N];
            H[i] = new long double[N];
        }
        generate_symmetric_matrix(A, H, eigenvalues, N);
        long double lambda_last = eigenvalues[N - 1];
        long double lambda_prelast = eigenvalues[N - 2];
        long double lambda_true = eigenvalues[N - 3];
        long double* x_last = new long double[N];
        long double* x_prelast = new long double[N];
        long double* x_true = new long double[N];
        for (int i = 0; i < N; ++i)
        {
            x_last[i] = H[i][N - 1];
            x_prelast[i] = H[i][N - 2];
            x_true[i] = H[i][N - 3];
        }
        long double lambda_3; // что ищем
        long double* x_3 = new long double[N]; //что ищем
        int K;
        long double r;
        straight_iteration_exhaust(A, N, lambda_last, x_last, lambda_prelast, x_prelast, lambda_3, x_3, epsilon, M, K, r, avg_diff_vec, lambda_true, x_true, avg_diff_lambda);
        //avg_diff_lambda += fabs(lambda_3 - lambda_true);
        avg_r += r;
        avg_iterations += K;


        for (int i = 0; i < N; ++i)
        {
            delete[] A[i];
            delete[] H[i];
        }
        delete[] A;
        delete[] H;
        delete[] eigenvalues;
        delete[] x_last;
        delete[] x_prelast;
        delete[] x_true;
        delete[] x_3;
    }
    cout << N << 'x' << N << "  [" << lower_bound << "; " << upper_bound << "]  eps = " << epsilon << ":\nСред. оценка точности собств. значений = " << avg_diff_lambda / numTests << "\nСред. оценка точности точности собств. векторов = " << avg_diff_vec / numTests << "\nСред. мера точности r = " << avg_r / numTests << "\nСред. число операций = " << avg_iterations / numTests << "\n\n";

}

int main()
{
    setlocale(LC_ALL, "rus");
    srand(time(NULL));
    int numTests = 20;
    int M = 10000000; // макс кол-во итераций

    int var;
    cout << "1) Одиночный тест   2) Полный прогон" << endl;
    cin >> var;

    if (var == 1)
    {
        //#1
        int N = 30;
        long double lower_bound = -50.0, upper_bound = 50.0;
        long double* eigenvalues = new long double[N];
        generate_random_eigenvalues(eigenvalues, N, lower_bound, upper_bound);
        cout << "Набор собственных значений:" << endl;
        for (int i = 0; i < N; i++)
            cout << eigenvalues[i] << ' ';
        cout << endl;

        long double** A = new long double* [N];
        long double** H = new long double* [N];
        for (int i = 0; i < N; ++i)
        {
            A[i] = new long double[N];
            H[i] = new long double[N];
        }
        generate_symmetric_matrix(A, H, eigenvalues, N);


        long double lambda_last = eigenvalues[N - 1];
        long double lambda_prelast = eigenvalues[N - 2];
        long double lambda_true = eigenvalues[N - 3];
        long double* x_last = new long double[N];
        long double* x_prelast = new long double[N];
        long double* x_true = new long double[N];
        for (int i = 0; i < N; ++i)
        {
            x_last[i] = H[i][N - 1];
            x_prelast[i] = H[i][N - 2];
            x_true[i] = H[i][N - 3];
        }
        long double epsilon = 1e-5; // точность
        long double lambda_3; // что ищем
        long double* x_3 = new long double[N]; //что ищем
        int K;
        long double r;

        long double avg_diff_lambda = 0.0, avg_diff_vec = 0.0, avg_iterations = 0.0, avg_r = 0.0;
        straight_iteration_exhaust(A, N, lambda_last, x_last, lambda_prelast, x_prelast, lambda_3, x_3, epsilon, M, K, r, avg_diff_vec, lambda_true, x_true, avg_diff_lambda);
        //avg_diff_lambda += fabs(lambda_3 - lambda_true);
        avg_r += r;
        avg_iterations += K;

        cout << "Результаты метода прямых итераций с исчерпыванием:\n";
        cout << "Третье максимальное по модулю собственное значение: " << lambda_3 << " ~ " << eigenvalues[N - 3] << endl;
        cout << "Собственный вектор, соответствующий третьему максимальному собственному значению:\n";
        for (int i = 0; i < N; ++i)
        {
            cout << x_3[i] << " ";
        }
        cout << "\n";
        cout << endl;
        cout << "Число выполненных итераций: " << K << endl;
        cout << "Мера точности r: " << r << endl;
        cout << "Оценка точности собств. значений: " << avg_diff_lambda << endl;
        cout << "Оценка точности собств. векторов: " << avg_diff_vec;



        for (int i = 0; i < N; ++i)
        {
            delete[] A[i];
        }
        delete[] A;
        delete[] eigenvalues;
        delete[] x_last;
        delete[] x_prelast;

        return 0;
    }
    else if (var == 2)
    {
        //#1
        int N = 10;
        long double lower_bound = -2.0, upper_bound = 2.0;
        long double epsilon = 1e-05;
        make_test(N, numTests, M, lower_bound, upper_bound, epsilon);

        //#2
        epsilon = 1e-08;
        make_test(N, numTests, M, lower_bound, upper_bound, epsilon);

        //#3
        lower_bound = -50.0, upper_bound = 50.0;
        epsilon = 1e-05;
        make_test(N, numTests, M, lower_bound, upper_bound, epsilon);

        //#4
        epsilon = 1e-08;
        make_test(N, numTests, M, lower_bound, upper_bound, epsilon);

        //#5
        N = 30;
        lower_bound = -2.0, upper_bound = 2.0;
        epsilon = 1e-05;
        make_test(N, numTests, M, lower_bound, upper_bound, epsilon);

        //#6
        epsilon = 1e-08;
        make_test(N, numTests, M, lower_bound, upper_bound, epsilon);

        //#7
        lower_bound = -50.0, upper_bound = 50.0;
        epsilon = 1e-05;
        make_test(N, numTests, M, lower_bound, upper_bound, epsilon);

        //#8
        epsilon = 1e-08;
        make_test(N, numTests, M, lower_bound, upper_bound, epsilon);

        //#9
        N = 50;
        lower_bound = -2.0, upper_bound = 2.0;
        epsilon = 1e-05;
        make_test(N, numTests, M, lower_bound, upper_bound, epsilon);

        //#10
        epsilon = 1e-08;
        make_test(N, numTests, M, lower_bound, upper_bound, epsilon);

        //#11
        lower_bound = -50.0, upper_bound = 50.0;
        epsilon = 1e-05;
        make_test(N, numTests, M, lower_bound, upper_bound, epsilon);

        //#12
        epsilon = 1e-08;
        make_test(N, numTests, M, lower_bound, upper_bound, epsilon);

        return 0;
    }
}
