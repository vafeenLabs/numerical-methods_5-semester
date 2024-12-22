#include <iostream>
#include <cmath>
#include <random>
#include <ctime>

using namespace std;

// Сортировка массива по возрастанию модулей элементов
void sort_array_abs(long double *arr, int size)
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

// Генерация случайного нормализованного вектора
// Нормализация означает, что длина (норма) вектора равна 1
void generate_random_vector(long double *omega, int N)
{
    long double norm = 0.0; // Переменная для хранения длины вектора
    for (int i = 0; i < N; ++i)
    {
        // Генерация случайных чисел от 0 до 1
        omega[i] = (rand() / static_cast<long double>(RAND_MAX));
        // Вычисление квадрата нормы
        norm += omega[i] * omega[i];
    }

    // Вычисление длины (нормы) вектора
    norm = sqrt(norm);
    for (int i = 0; i < N; ++i)
    {
        // Нормализация каждого элемента вектора
        omega[i] /= norm;
    }

    // В этот момент вектор omega нормализован
}

// Умножение матрицы на вектор
// result = A * x
void mat_vec_mul(long double **A, long double *x, long double *result, int N)
{
    for (int i = 0; i < N; ++i)
    {
        result[i] = 0.0; // Инициализация текущего элемента результата
        for (int j = 0; j < N; ++j)
        {
            // Вычисление скалярного произведения строки матрицы на вектор
            result[i] += A[i][j] * x[j];
        }
    }
}

// Генерация вектора, ортогонального заданному входному вектору
// Выходной вектор находится в плоскости, заданной первыми двумя компонентами
void generate_orth_vector(long double *input_vec, int N, long double *out_vec)
{
    // Задание первых двух компонент ортогонального вектора
    out_vec[0] = (-1.0) * input_vec[1];
    out_vec[1] = input_vec[0];

    // Заполнение остальных компонент нулями
    for (int i = 2; i < N; i++)
        out_vec[i] = 0.0;
}

// Вычисление скалярного произведения двух векторов
long double dot_product(long double *a, long double *b, int N)
{
    long double dot = 0.0; // Переменная для накопления результата
    for (int i = 0; i < N; i++)
        dot += a[i] * b[i]; // Суммирование произведений соответствующих компонент
    return dot;
}

// Вычисление длины (нормы) вектора
long double vec_length(long double *vec, int N)
{
    long double norm = 0.0; // Переменная для накопления суммы квадратов компонент
    for (int i = 0; i < N; i++)
        norm += vec[i] * vec[i];

    // Возврат квадратного корня из суммы квадратов (длина вектора)
    return sqrt(norm);
}
