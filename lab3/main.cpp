#include "main.hpp"

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

// Генерация симметричной матрицы с заданными собственными значениями
void generate_symmetric_matrix(long double **A, long double **H, long double *eigenvalues, int N)
{
    long double *omega = new long double[N];
    generate_random_vector(omega, N);

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            H[i][j] = (i == j ? 1.0 : 0.0) - 2 * omega[i] * omega[j];
        }
    }

    cout << "\nМатрица H\n";
    printMatrix(H, N);

    long double **Amid = new long double *[N];
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

    cout << "\nМатрица A\n";
    printMatrix(A, N);

    for (int i = 0; i < N; ++i)
    {
        delete[] Amid[i];
    }
    delete[] Amid;
    delete[] omega;
}

long double angle_vectors(long double *a, long double *b, int N)
{
    long double dot = dot_product(a, b, N);
    long double len_a = vec_length(a, N);
    long double len_b = vec_length(b, N);

    long double cos = dot / (len_a * len_b);

    if (cos > 1.0)
        cos = 1.0;
    if (cos < -1.0)
        cos = -1.0;

    return fabs(acos(fabs(cos)));
}

// Метод обратных итераций с исчерпыванием
void straight_iteration_exhaust(long double **A, int N, long double lambda_1, long double *x_1, long double lambda_2, long double *x_2, long double &lambda_3, long double *x_3, long double epsilon, int M, int &K, long double &r, long double &avg_vec, long double lambda_true, long double *x_true, long double &avg_lambda)
{
    long double *v = new long double[N];        // Создание динамического массива v размером N для хранения текущего вектора
    long double **temp = new long double *[N];  // Создание массива указателей на строки временной матрицы temp
    long double *prev_v = new long double[N];   // Массив для хранения предыдущего вектора v
    long double *prev_x_3 = new long double[N]; // Массив для хранения предыдущего собственного вектора x_3

    for (int i = 0; i < N; i++) // Цикл для инициализации массивов
    {
        temp[i] = new long double[N]; // Выделение памяти под каждую строку временной матрицы temp
        prev_x_3[i] = x_2[i];         // Инициализация prev_x_3 значениями из x_2 перед началом итераций
    }

    for (int i = 0; i < N; ++i) // Генерация ортогонального вектора x_3 на основе x_2
    {
        generate_orth_vector(x_2, N, x_3); // Вызов функции для генерации ортогонального вектора
    }

    // Первый проход - вычитание x_1
    // Формирование матрицы A(1) = A - lambda_n * x_1 * x_1
    for (int i = 0; i < N; ++i) // Цикл по строкам матрицы temp
    {
        for (int j = 0; j < N; ++j) // Цикл по столбцам матрицы temp
        {
            temp[i][j] = A[i][j] - lambda_1 * x_1[i] * x_1[j]; // Заполнение матрицы temp вычитанием произведения lambda_1 и внешнего произведения x_1 с самим собой
        }
    }

    // Второй проход - вычитание x_2
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            temp[i][j] = temp[i][j] - lambda_2 * x_2[i] * x_2[j]; // Корректировка значений в temp с использованием второго собственного значения и вектора
        }
    }

    // PRINT(temp) // Закомментированная строка для вывода временной матрицы

    long double alpha_k = 0.0, prev_alpha = 0.0, vecDiff = 10.0, norm;
    K = 0;

    for (int k = 0; k < M; ++k) // Основной цикл итераций
    {
        ++K; // Увеличение счетчика итераций

        // вычисление нормы вектора
        norm = 0.0; // Инициализация переменной для вычисления нормы вектора
        for (int i = 0; i < N; ++i)
        {
            norm += x_3[i] * x_3[i]; // Вычисление суммы квадратов элементов вектора x_3 для нахождения нормы
        }
        norm = sqrt(norm);

        // нормализация вектора
        for (int i = 0; i < N; ++i)
        {
            if (K >= 2)
                prev_v[i] = v[i]; // Сохранение предыдущего значения вектора v при второй и последующих итерациях

            v[i] = x_3[i] / norm; // Нормализация вектора v путем деления каждого элемента на его норму
        }

        // Вычисление нового вектора: x_k+1 = A(1) * v(k)
        mat_vec_mul(temp, v, x_3, N); // Умножение матрицы temp на вектор v с сохранением результата в x_3

        prev_alpha = alpha_k; // Сохранение предыдущего значения alpha_k перед его обновлением
        alpha_k = 0.0;

        for (int i = 0; i < N; ++i)
        {
            alpha_k += v[i] * x_3[i]; // Вычисление alpha_k как скалярного произведения v и x_3
        }

        if (k >= 2)
            vecDiff = angle_vectors(v, prev_v, N); // Вычисление разности углов между текущим и предыдущим вектором

        // Проверка условий остановки: если погрешность меньше заданной точности epsilon,
        if (fabs(alpha_k - prev_alpha) < epsilon && vecDiff < epsilon && fabs(alpha_k - lambda_true) < epsilon)
        {
            break;
        }
    }

    avg_vec += vecDiff;                       // Обновление средней разности углов между векторами
    avg_lambda += fabs(alpha_k - prev_alpha); // Обновление средней разности собственных значений

    lambda_3 = alpha_k; // Сохранение найденного собственного значения
    for (int i = 0; i < N; i++)
        x_3[i] = v[i]; // Копирование нормализованного вектора v в результатирующий собственный вектор x_3

    // Вычисление меры точности r
    long double *Ax = new long double[N];       // Создание динамического массива Ax размером N для хранения результата умножения A на собственный вектор
    long double *lambda_x = new long double[N]; // Создание динамического массива lambda_x размером N для хранения произведения найденного собственного значения на собственный вектор

    mat_vec_mul(A, x_3, Ax, N); // Умножение матрицы A на собственный вектор x_3 с сохранением результата в Ax

    for (int i = 0; i < N; ++i)
    {
        lambda_x[i] = lambda_3 * x_3[i]; // Вычисление lambda_x как произведение найденного значения lambda_3 и собственного вектора x_3
    }

    r = 0.0;
    for (int i = 0; i < N; ++i)
    {
        r = max(r, fabs(Ax[i] - lambda_x[i])); // Обновление меры точности r как максимального значения разности между Ax и lambda_x
    }

    for (int i = 0; i < N; i++)
        delete[] temp[i]; // Освобождение памяти выделенной под строки временной матрицы temp

    delete[] v;        // Освобождение памяти для массива v
    delete[] prev_v;   // Освобождение памяти для массива prev_v
    delete[] temp;     // Освобождение памяти для массива указателей на строки временной матрицы temp
    delete[] prev_x_3; // Освобождение памяти для массива prev_x_3
    delete[] Ax;       // Освобождение памяти для массива Ax
    delete[] lambda_x; // Освобождение памяти для массива lambda_x
}

// Функция для генерации массива случайных чисел
void generate_random_eigenvalues(long double *eigenvalues, int N, long double from_number, long double to_number)
{
    for (int i = 0; i < N; ++i)
    {
        eigenvalues[i] = from_number + (rand() / static_cast<long double>(RAND_MAX)) * (to_number - from_number);
    }
    sort_array_abs(eigenvalues, N);
}

void make_test(int N, int numTests, int M, long double from_number, long double to_number, long double epsilon)
{
    long double avg_diff_lambda = 0.0, avg_diff_vec = 0.0, avg_iterations = 0.0, avg_r = 0.0;
    for (int i = 0; i < numTests; i++)
    {
        long double *eigenvalues = new long double[N];
        generate_random_eigenvalues(eigenvalues, N, from_number, to_number);
        long double **A = new long double *[N];
        long double **H = new long double *[N];
        for (int i = 0; i < N; ++i)
        {
            A[i] = new long double[N];
            H[i] = new long double[N];
        }
        generate_symmetric_matrix(A, H, eigenvalues, N);
        long double lambda_last = eigenvalues[N - 1];
        long double lambda_prelast = eigenvalues[N - 2];
        long double lambda_true = eigenvalues[N - 3];
        long double *x_last = new long double[N];
        long double *x_prelast = new long double[N];
        long double *x_true = new long double[N];
        for (int i = 0; i < N; ++i)
        {
            x_last[i] = H[i][N - 1];
            x_prelast[i] = H[i][N - 2];
            x_true[i] = H[i][N - 3];
        }
        long double lambda_3;                  // что ищем
        long double *x_3 = new long double[N]; // что ищем
        int K;
        long double r;
        straight_iteration_exhaust(A, N, lambda_last, x_last, lambda_prelast, x_prelast, lambda_3, x_3, epsilon, M, K, r, avg_diff_vec, lambda_true, x_true, avg_diff_lambda);
        // avg_diff_lambda += fabs(lambda_3 - lambda_true);
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
    cout << N << 'x' << N << "  [" << from_number << "; " << to_number << "]  eps = " << epsilon << ":\nСред. оценка точности собств. значений = " << avg_diff_lambda / numTests << "\nСред. оценка точности точности собств. векторов = " << avg_diff_vec / numTests << "\nСред. мера точности r = " << avg_r / numTests << "\nСред. число операций = " << avg_iterations / numTests << "\n\n";
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
        // #1
        int N = 10;
        long double from_number = -50.0, to_number = 50.0;
        long double *eigenvalues = new long double[N];
        generate_random_eigenvalues(eigenvalues, N, from_number, to_number);
        cout << "Набор собственных значений:" << endl;
        for (int i = 0; i < N; i++)
            cout << eigenvalues[i] << ' ';
        cout << endl;

        long double **A = new long double *[N];
        long double **H = new long double *[N];
        for (int i = 0; i < N; ++i)
        {
            A[i] = new long double[N];
            H[i] = new long double[N];
        }
        generate_symmetric_matrix(A, H, eigenvalues, N);

        long double lambda_last = eigenvalues[N - 1];
        long double lambda_prelast = eigenvalues[N - 2];
        long double lambda_true = eigenvalues[N - 3];
        long double *x_last = new long double[N];
        long double *x_prelast = new long double[N];
        long double *x_true = new long double[N];
        for (int i = 0; i < N; ++i)
        {
            x_last[i] = H[i][N - 1];
            x_prelast[i] = H[i][N - 2];
            x_true[i] = H[i][N - 3];
        }
        long double epsilon = 1e-8;            // точность
        long double lambda_3;                  // что ищем
        long double *x_3 = new long double[N]; // что ищем
        int K;
        long double r;

        long double avg_diff_lambda = 0.0, avg_diff_vec = 0.0, avg_iterations = 0.0, avg_r = 0.0;
        straight_iteration_exhaust(A, N, lambda_last, x_last, lambda_prelast, x_prelast, lambda_3, x_3, epsilon, M, K, r, avg_diff_vec, lambda_true, x_true, avg_diff_lambda);
        // avg_diff_lambda += fabs(lambda_3 - lambda_true);
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
        // #1
        int N = 10;
        long double from_number = -2.0, to_number = 2.0;
        long double epsilon = 1e-05;
        make_test(N, numTests, M, from_number, to_number, epsilon);

        // #2
        epsilon = 1e-08;
        make_test(N, numTests, M, from_number, to_number, epsilon);

        // #3
        from_number = -50.0, to_number = 50.0;
        epsilon = 1e-05;
        make_test(N, numTests, M, from_number, to_number, epsilon);

        // #4
        epsilon = 1e-08;
        make_test(N, numTests, M, from_number, to_number, epsilon);

        // #5
        N = 30;
        from_number = -2.0, to_number = 2.0;
        epsilon = 1e-05;
        make_test(N, numTests, M, from_number, to_number, epsilon);

        // #6
        epsilon = 1e-08;
        make_test(N, numTests, M, from_number, to_number, epsilon);

        // #7
        from_number = -50.0, to_number = 50.0;
        epsilon = 1e-05;
        make_test(N, numTests, M, from_number, to_number, epsilon);

        // #8
        epsilon = 1e-08;
        make_test(N, numTests, M, from_number, to_number, epsilon);

        // #9
        N = 50;
        from_number = -2.0, to_number = 2.0;
        epsilon = 1e-05;
        make_test(N, numTests, M, from_number, to_number, epsilon);

        // #10
        epsilon = 1e-08;
        make_test(N, numTests, M, from_number, to_number, epsilon);

        // #11
        from_number = -50.0, to_number = 50.0;
        epsilon = 1e-05;
        make_test(N, numTests, M, from_number, to_number, epsilon);

        // #12
        epsilon = 1e-08;
        make_test(N, numTests, M, from_number, to_number, epsilon);

        return 0;
    }
}
