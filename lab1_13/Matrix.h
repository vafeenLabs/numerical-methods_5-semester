#pragma once
#include<iostream>
#include <fstream>
#include <random>
/// 											
/// 											
/// 		*	*								*
/// 		*	*	*							*
/// 		*	*	*	*						*
/// 		*		*	*	*					*
/// 		*			*	*	*				*
/// 		*				*	*	*			*
/// 		*					*	*	*		*
/// 		*						*	*	*	*
/// 		*							*	*	*
/// 		*								*	*
/// 
///         p1                          c     b     a   pn
class Matrix
{
public:
    Matrix(size_t size, const std::string& fileName);
    Matrix(size_t size, double range_min, double range_max);

    void solve();
    double* getSolution();
    double getE1() { return E1; }
    double getE2() { return E2; }
    void getArrays(double a[], double b[], double c[], double p1[], double pn[], double rsh[]);
    bool isSolved() { return solved; }
    size_t getSize() { return SIZE; }

    double getDefaultValueMultipleSolutions() { return defaultValueMultipleSolutions; }
    void setDefaultValueMultipleSolutions(double defaultValueMultipleSolutions) { this->defaultValueMultipleSolutions = defaultValueMultipleSolutions; }

    double getQ() { return q; }
    void setQ(double _q) { this->q = _q; }
private:
    ~Matrix();
    void getRhsForSystemErrors();
    bool checkSolution();
    void printArray(double arr[], int size);
    void fillMatrix(const std::string& fileName);
    bool isCorrectToBeSolved();

    size_t SIZE;
    double defaultValueMultipleSolutions = -999;
    double q = 0.0000000001;
    double* a; double* b; double* c;
    double* p1; double* pn;
    double* rhs;
    double* x;

    bool solved;
    double* testA; double* testB; double* testC;
    double* testP1; double* testPn;
    double* testRhs;

    double* e1Rhs;
    double* e2Rhs;
    double* e2X;

    double* eSolution1X;
    double* eSolution2X;

    double E1 = 11;
    double E2 = 11;
    double E = 0.0000001;
};
double generateRandomNumber(double range_min, double range_max) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(range_min, range_max);
    return dis(gen);
}

void Matrix::solve()
{
    getRhsForSystemErrors();

    for (size_t j = 1; j < SIZE; j++)
    {
        p1[j - 1] = b[j - 1];
        p1[j] = c[j - 1];

        bool ab = false; // b находиться под a, может измениться при вычислениях
        bool ac = false;// с находиться под а, может измениться при вычислениях
        for (int i = j; i < SIZE; i++) {
            double multiplier = 0;
            if (p1[j - 1] == 0) // если элемент главной диагонали = 0, меняем его со следующей строкой матрицы
            {
                double tmp = p1[j - 1];
                p1[j - 1] = p1[j];
                p1[j] = tmp;

                tmp = pn[j - 1];
                pn[j - 1] = pn[j];
                pn[j] = tmp;

                tmp = rhs[j - 1];
                rhs[j - 1] = rhs[j];
                rhs[j] = tmp;
                tmp = e1Rhs[j - 1];
                e1Rhs[j - 1] = e1Rhs[j];
                e1Rhs[j] = tmp;
                tmp = e2Rhs[j - 1];
                e2Rhs[j - 1] = e2Rhs[j];
                e2Rhs[j] = tmp;

                tmp = a[j - 1];
                a[j - 1] = b[j];
                b[j] = tmp;

                if (j < SIZE - 2)
                {
                    x[j - 1] = a[j];
                    a[j] = 0;
                }
                continue;   // следующий элемент стал нулём, нет смысла вычислять множитель
            }

            multiplier = p1[i] / p1[j - 1]; // ищем множитель

            p1[i] = (-1) * (a[j - 1] * multiplier);

            if (!ab)
            {
                b[i] -= a[j - 1] * multiplier;
                ab = true;
            }
            else if (!ac)
            {
                c[i - 1] -= a[j - 1] * multiplier;
                ac = true;
            }
            pn[i] -= pn[j - 1] * multiplier;

            rhs[i] -= rhs[j - 1] * multiplier;
            e1Rhs[i] -= e1Rhs[j - 1] * multiplier;
            e2Rhs[i] -= e2Rhs[j - 1] * multiplier;

        }
    }

    p1[SIZE - 1] = pn[SIZE - 1];    // последние значения уже расчитаны
    a[SIZE - 2] = pn[SIZE - 2];

    if (isCorrectToBeSolved())
    {
        double* a1 = new double[SIZE - 1];
        double* a2 = new double[SIZE - 1];
        for (size_t i = 0; i < SIZE - 1; i++)
        {
            a1[i] = a[i];
            a2[i] = a[i];
        }

        x[SIZE - 1] = (p1[SIZE - 1] == 0) ? defaultValueMultipleSolutions : (rhs[SIZE - 1] / p1[SIZE - 1]);
        eSolution1X[SIZE - 1] = (p1[SIZE - 1] == 0) ? defaultValueMultipleSolutions : (e1Rhs[SIZE - 1] / p1[SIZE - 1]);
        eSolution2X[SIZE - 1] = (p1[SIZE - 1] == 0) ? defaultValueMultipleSolutions : (e2Rhs[SIZE - 1] / p1[SIZE - 1]);

        for (size_t i = 0; i < SIZE; i++)
        {
            rhs[i] -= pn[i] * x[SIZE - 1];
            e1Rhs[i] -= pn[i] * eSolution1X[SIZE - 1];
            e2Rhs[i] -= pn[i] * eSolution2X[SIZE - 1];
        }
        x[SIZE - 2] = (p1[SIZE - 2] == 0) ? defaultValueMultipleSolutions : (rhs[SIZE - 2] / p1[SIZE - 2]);
        eSolution1X[SIZE - 2] = (p1[SIZE - 2] == 0) ? defaultValueMultipleSolutions : (e1Rhs[SIZE - 2] / p1[SIZE - 2]);
        eSolution2X[SIZE - 2] = (p1[SIZE - 2] == 0) ? defaultValueMultipleSolutions : (e2Rhs[SIZE - 2] / p1[SIZE - 2]);
        for (int i = SIZE - 3; i >= 0; --i)
        {
            if (x[i] != 0)
            {
                rhs[i] -= x[i] * x[i + 2];
                e1Rhs[i] -= eSolution1X[i] * eSolution1X[i + 2];
                e2Rhs[i] -= eSolution2X[i] * eSolution2X[i + 2];
            }

            a[i] *= x[i + 1];
            x[i] = (p1[i] == 0) ? defaultValueMultipleSolutions : ((rhs[i] - a[i]) / p1[i]);

            a1[i] *= eSolution1X[i + 1];
            eSolution1X[i] = (p1[i] == 0) ? defaultValueMultipleSolutions : ((e1Rhs[i] - a1[i]) / p1[i]);
            a2[i] *= eSolution2X[i + 1];
            eSolution2X[i] = (p1[i] == 0) ? defaultValueMultipleSolutions : ((e2Rhs[i] - a2[i]) / p1[i]);

        }

        solved = true;
        solved = checkSolution();

        if (solved)
        {
            //printArray(eSolution2X, SIZE);

            double Er1 = 11;
            double Er2 = 11;
            for (size_t i = 0; i < SIZE; i++)
            {
                double er1 = (eSolution1X[i] - 1) < 0 ? (eSolution1X[i] * (-1) + 1) : (eSolution1X[i] - 1);
                double er2 = (eSolution2X[i] - e2X[i]) < 0 ? (eSolution2X[i] * (-1) + e2X[i]) : (eSolution2X[i] - e2X[i]);
                if (e2X[i] > q || (((-1) * e2X[i]) > q))
                {
                    if (e2X[i] < 0)
                        e2X[i] *= -1;
                    er2 /= e2X[i];
                }
                if (Er1 < er1 || Er1>10)
                {
                    Er1 = er1;
                }
                if (Er2 < er2 || Er2>10)
                {
                    Er2 = er2;
                }
            }
            E1 = Er1;
            E2 = Er2;
        }
    }
}

Matrix::Matrix(size_t size, const std::string& fileName)
{
    if (size > 0)
    {
        SIZE = size;
        a = new double[SIZE - 1];
        testA = new double[SIZE - 1];
        c = new double[SIZE - 1];
        testC = new double[SIZE - 1];
        b = new double[SIZE];
        testB = new double[SIZE];

        p1 = new double[SIZE];
        testP1 = new double[SIZE];
        pn = new double[SIZE];
        testPn = new double[SIZE];

        rhs = new double[SIZE];
        testRhs = new double[SIZE];
        x = new double[SIZE];

        e1Rhs = new double[SIZE];
        e2Rhs = new double[SIZE];
        e2X = new double[SIZE];

        eSolution1X = new double[SIZE];
        eSolution2X = new double[SIZE];

        solved = false;
        fillMatrix(fileName);
        // Вывод массивов на экран
        /*
        std::cout << "A: ";
        printArray(a, SIZE - 1);
        std::cout << "B: ";
        printArray(b, SIZE);
        std::cout << "C: ";
        printArray(c, SIZE - 1);
        std::cout << "P1: ";
        printArray(p1, SIZE);
        std::cout << "Pn: ";
        printArray(pn, SIZE);
        std::cout << "Rhs: ";
        printArray(rhs, SIZE);
        */
    }
    else
    {
        std::cout << "incorrect size for matrix!" << std::endl;
    }
}
Matrix::Matrix(size_t size, double range_min, double range_max)
{
    SIZE = size;
    a = new double[SIZE - 1];
    testA = new double[SIZE - 1];
    c = new double[SIZE - 1];
    testC = new double[SIZE - 1];
    b = new double[SIZE];
    testB = new double[SIZE];

    p1 = new double[SIZE];
    testP1 = new double[SIZE];
    pn = new double[SIZE];
    testPn = new double[SIZE];

    rhs = new double[SIZE];
    testRhs = new double[SIZE];
    x = new double[SIZE];
    e1Rhs = new double[SIZE];
    e2Rhs = new double[SIZE];
    e2X = new double[SIZE];

    eSolution1X = new double[SIZE];
    eSolution2X = new double[SIZE];

    solved = false;


    // Заполнение массивов
    for (size_t i = 0; i < SIZE; ++i) {
        x[i] = 0;
        double number;
        number = generateRandomNumber(range_min, range_max/*-0.1, 0.1*/);
        e2X[i] = number;
        //printArray(e2X, SIZE);
        //std::cout << e2X[i]<<' ';

        if (i < SIZE - 1)
        {
            number = number = generateRandomNumber(range_min, range_max);
            a[i] = number;
            number = number = generateRandomNumber(range_min, range_max);
            c[i] = number;
        }

        number = generateRandomNumber(range_min, range_max);
        p1[i] = number;

        number = generateRandomNumber(range_min, range_max);
        pn[i] = number;

        number = generateRandomNumber(range_min, range_max);
        rhs[i] = number;

        number = generateRandomNumber(range_min, range_max); // для b[i]
        if (i == 0 && number == 0)
            ++number;
        b[i] = number;
        if (b[i] == 0 && i > 0)
        {
            if (a[i - 1] == 0)
            {
                number = generateRandomNumber(range_min, range_max); // для b=0
                if (number == 0)
                    ++number;
                a[i - 1] = number; //  для a[i-1]
            }
            if (a[i] == 0)
            {
                number = generateRandomNumber(range_min, range_max); // для b=0
                if (number == 0)
                    ++number;
                a[i] = number; //  для a[i]
            }
        }
    }
    p1[0] = b[0];
    p1[1] = c[0];
    pn[SIZE - 1] = b[SIZE - 1];
    pn[SIZE - 2] = a[SIZE - 2];

    // Вывод массивов на экран
    /*
    std::cout << "NEWWWWWWWWWWWW: ";
    std::cout << "A: ";
    printArray(a, SIZE-1);
    std::cout << "B: ";
    printArray(b, SIZE);
    std::cout << "C: ";
    printArray(c, SIZE-1);
    std::cout << "P1: ";
    printArray(p1, SIZE);
    std::cout << "Pn: ";
    printArray(pn, SIZE);
    std::cout << "Rhs: ";
    printArray(rhs, SIZE);
    */

    for (size_t i = 0; i < SIZE; i++)
    {
        if (i < SIZE - 1)
        {
            testA[i] = a[i];
            testC[i] = c[i];
        }
        testB[i] = b[i];
        testP1[i] = p1[i];
        testPn[i] = pn[i];
        testRhs[i] = rhs[i];
    }

}
void Matrix::fillMatrix(const std::string& fileName)
{
    std::ifstream file(fileName);
    if (file.is_open())
    {
        double value; double min = 1, max = 1;
        for (int i = 0; i < SIZE; ++i)
        {
            x[i] = 0;
            for (int j = 0; j < SIZE; ++j)
            {
                double value;
                file >> value;
                if (value < min)
                {
                    min = value;
                }
                if (value > max)
                {
                    max = value;
                }

                if (i == j)
                {
                    b[i] = value;
                }
                else if ((j - i) == 1)
                {
                    a[i] = value;
                }
                else if ((i - j) == 1)
                {
                    c[i - 1] = value;
                }
                else if (j == 0)
                {
                    p1[i] = value;
                }
                else if (j == SIZE - 1)
                {
                    pn[i] = value;
                }
            }
        }
        p1[0] = b[0];
        p1[1] = c[0];
        pn[SIZE - 1] = b[SIZE - 1];
        pn[SIZE - 2] = a[SIZE - 2];

        for (size_t i = 0; i < SIZE; i++)
        {
            double value;
            file >> value;


            rhs[i] = value;
        }

        file.close();

        for (size_t i = 0; i < SIZE; i++)
        {
            if (i < SIZE - 1)
            {
                testA[i] = a[i];
                testC[i] = c[i];
            }
            testB[i] = b[i];
            testP1[i] = p1[i];
            testPn[i] = pn[i];
            testRhs[i] = rhs[i];
        }
        for (size_t i = 0; i < SIZE; i++)
            e2X[i] = generateRandomNumber(min, max);
    }
    else
    {
        std::cout << "Can`t open the file: " << fileName << std::endl;
    }

}

bool Matrix::checkSolution()
{
    if (solved)
    {
        double* testS = new double[SIZE];
        bool result = true;
        for (size_t i = 0; i < SIZE && result; i++)
        {
            testS[i] = testRhs[i];

            if (i == 0)
            {
                testS[i] -= (testP1[i] * x[0] + testPn[i] * x[SIZE - 1] + testA[i] * x[1]);
            }
            else if (i == SIZE - 1)
            {
                testS[i] -= (testP1[i] * x[0] + testPn[i] * x[SIZE - 1] + testC[i - 1] * x[SIZE - 2]);
            }
            else if (i == SIZE - 2)
            {
                testS[i] -= (testB[i] * x[i] + testP1[i] * x[0] + testPn[i] * x[SIZE - 1] + testC[i - 1] * x[SIZE - 3]);
            }
            else if (i == 1)
            {
                testS[i] -= (testB[i] * x[i] + testP1[i] * x[0] + testPn[i] * x[SIZE - 1] + testA[i] * x[2]);
            }
            else
            {
                testS[i] -= (testP1[i] * x[0] + testPn[i] * x[SIZE - 1] + testB[i] * x[i] + testA[i] * x[i + 1] + testC[i - 1] * x[i - 1]);
            }
            if (testS[i] < 0)
            {
                result = ((testS[i] * (-1)) <= E);
            }
            else
            {
                result = ((testS[i]) <= E);
            }
        }
        return result;
    }
    else
    {
        return false;
    }
}
bool Matrix::isCorrectToBeSolved()
{
    double* testRHS = new double[SIZE];


    if ((rhs[SIZE - 1] > E || rhs[SIZE - 1] < -E) && p1[SIZE - 1] == 0) // -1*1 + 0*x + 1*1 = 10
        return false;

    /*for (size_t i = 0; i < SIZE; i++)
        testRHS[i] = rhs[i] - pn[i] * x[SIZE - 1];


    if ((testRHS[SIZE - 2] > E || testRHS[SIZE - 2] < -E) && p1[SIZE - 2] == 0)
        return false;

    for (int i = SIZE - 3; i >= 0; --i)
    {
        if (x[i] != 0)
        {
            testRHS[i] = rhs[i] - (x[i] * x[i + 2]);
        }
        if (((testRHS[i] - a[i] * x[i + 1]) > E || (testRHS[i] - a[i] * x[i + 1]) < -E) && p1[i] == 0)
            return false;
    }*/

    return true;
}

void Matrix::getRhsForSystemErrors()
{
    e1Rhs[0] = testP1[0] + testA[0] + testPn[0];
    e1Rhs[1] = testP1[1] + testA[1] + testB[1] + testPn[1];

    e2Rhs[0] = testP1[0] * e2X[0] + testA[0] * e2X[1] + testPn[0] * e2X[SIZE - 1];
    e2Rhs[1] = testP1[1] * e2X[0] + testA[1] * e2X[2] + testB[1] * e2X[1] + testPn[1] * e2X[SIZE - 1];
    for (size_t i = 2; i < SIZE - 2; i++)
    {
        e1Rhs[i] = testP1[i] + testC[i - 1] + testA[i] + testB[i] + testPn[i];

        e2Rhs[i] = testP1[i] * e2X[0] + testC[i - 1] * e2X[i - 1] + testA[i] * e2X[i + 1] + testB[i] * e2X[i] + testPn[i] * e2X[SIZE - 1];
    }
    e1Rhs[SIZE - 2] = testP1[SIZE - 2] + testB[SIZE - 2] + testC[SIZE - 3] + testPn[SIZE - 2];
    e1Rhs[SIZE - 1] = testP1[SIZE - 1] + testC[SIZE - 2] + testPn[SIZE - 1];

    e2Rhs[SIZE - 2] = testP1[SIZE - 2] * e2X[0] + testB[SIZE - 2] * e2X[SIZE - 2] + testC[SIZE - 3] * e2X[SIZE - 3] + testPn[SIZE - 2] * e2X[SIZE - 1];
    e2Rhs[SIZE - 1] = testP1[SIZE - 1] * e2X[0] + testC[SIZE - 2] * e2X[SIZE - 2] + testPn[SIZE - 1] * e2X[SIZE - 1];

}

double* Matrix::getSolution()
{
    double* xCopy = new double[SIZE];
    std::memcpy(xCopy, x, SIZE * sizeof(double));
    return xCopy;
}
void Matrix::getArrays(double _a[], double _b[], double _c[], double _p1[], double _pn[], double _rsh[])
{
    for (size_t i = 0; i < SIZE; i++)
    {
        if (i < SIZE - 1)
        {
            _a[i] = testA[i];
            _c[i] = testC[i];
        }
        _b[i] = testB[i];
        _p1[i] = testP1[i];
        _pn[i] = testPn[i];
        _rsh[i] = testRhs[i];
    }
}

void Matrix::printArray(double arr[], int size)
{
    std::cout << std::endl;
    for (int j = 0; j < size; j++) {
        std::cout << '[' << j + 1 << ']' << '=' << arr[j] << std::endl;
    }
    std::cout << std::endl;
}
Matrix::~Matrix()
{
    delete[] a; a = nullptr;
    delete[] b; b = nullptr;
    delete[] c; c = nullptr;
    delete[] p1; p1 = nullptr;
    delete[] pn; pn = nullptr;
    delete[] rhs; rhs = nullptr;
    delete[] x; x = nullptr;

    delete[] testA; testA = nullptr;
    delete[] testB; testB = nullptr;
    delete[] testC; testC = nullptr;
    delete[] testP1; testP1 = nullptr;
    delete[] testPn; testPn = nullptr;
    delete[] testRhs; testRhs = nullptr;

    delete[] e1Rhs; e1Rhs = nullptr;
    delete[] e2Rhs; e2Rhs = nullptr;
    delete[] e2X; e2X = nullptr;

    delete[] eSolution1X; eSolution1X = nullptr;
    delete[] eSolution2X; eSolution2X = nullptr;
}