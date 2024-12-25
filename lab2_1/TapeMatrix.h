﻿#pragma once
#include "functions.h"
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <stdexcept>
#include <fstream>
#include <vector>


class TapeMatrix {
private:
    std::vector<std::vector<double>> matrix;
    std::vector<std::vector<double>> matrixCopy;

    std::vector<std::vector<double>> accuracyMatrixLU;
    int N; // Размер обычной матрицы NxN
    int L; // Половина ширины ленты

    bool solved = false;
    bool illConditionedMatrices = false;

    std::vector<double> x;
    std::vector<double> f;

    double q = 0.000000001;
    std::vector<double> accuracyX;
    std::vector<double> solutionForAccuracyX;
    std::vector<double> accuracyF;
    double meanRatioRelativeAccuracy = 0.0;

    std::vector<double> accuracyLUF;
    std::vector<double> solutionForAccuracyLUX;
    double meanRatioRelativeAccuracyIllConditionedMatrices = 0.0;

public:


    TapeMatrix(const std::string& filename, int n, int l) : N(n), L(l) {
        try {
            matrix.resize(N, std::vector<double>(2 * L - 1));
            matrix.reserve(N);
            matrixCopy.resize(N, std::vector<double>(2 * L - 1));
            matrixCopy.reserve(N);

            f.resize(N); f.reserve(N);
            x.resize(N); x.reserve(N);
            accuracyX.resize(N); accuracyX.reserve(N);
            accuracyF.resize(N); accuracyF.reserve(N);
            solutionForAccuracyX.resize(N); solutionForAccuracyX.reserve(N);
            std::ifstream file(filename);
            if (!file.is_open()) {
                throw std::runtime_error("Open file error.");
            }

            double value;
            for (int i = 0; i < N; ++i) {
                matrix[i].reserve(2 * L - 1);
                matrixCopy[i].reserve(2 * L - 1);
                accuracyX[i] = 1;
                for (int j = 0; j < N; ++j) {
                    file >> value;
                    int new_col = j - i + L - 1;
                    if (new_col >= 0 && new_col < 2 * L - 1)
                    {
                        matrix[i][new_col] = value;
                        matrixCopy[i][new_col] = value;
                    }
                }
            }
            matrix[N - 1][2 * L - 2] = 0;
            matrixCopy[N - 1][2 * L - 2] = 0;

            f.resize(N); f.reserve(N);
            for (size_t i = 0; i < N; i++)
            {
                file >> value;
                f[i] = (value);
            }

            getAccuracyF();
            file.close();
        }
        catch (const std::exception& e) {
            throw std::runtime_error( "Error: "+ std::string(e.what()));
        }
    }

    TapeMatrix(double minValue, double maxValue, int n, int l, bool illConditionedMatrices = false, int k = 0) : N(n), L(l) {
        try {
            matrix.resize(N, std::vector<double>(2 * L - 1));
            matrix.reserve(N);
            matrixCopy.resize(N, std::vector<double>(2 * L - 1));
            matrixCopy.reserve(N);
            this->illConditionedMatrices = illConditionedMatrices;

            if (illConditionedMatrices)
            {
                accuracyMatrixLU.resize(N, std::vector<double>(2 * L - 1));
                accuracyMatrixLU.reserve(N);
                accuracyLUF.resize(N); accuracyLUF.reserve(N);
                solutionForAccuracyLUX.resize(N); solutionForAccuracyLUX.reserve(N);
            }

            f.resize(N); f.reserve(N);
            x.resize(N); x.reserve(N);
            accuracyX.resize(N); accuracyX.reserve(N);
            accuracyF.resize(N); accuracyF.reserve(N);
            solutionForAccuracyX.resize(N); solutionForAccuracyX.reserve(N);
            int count = L - 1;
            for (size_t i = 0; i < N; i++) {
                matrix[i].reserve(2 * L - 1);
                matrixCopy[i].reserve(2 * L - 1);
                if (illConditionedMatrices)
                {
                    accuracyMatrixLU[i].reserve(2 * L - 1);
                }
                f[i] = generateRandomNumber(minValue, maxValue);
                accuracyX[i] = generateRandomNumber(minValue, maxValue);
                if (count < 2 * L - 1 && i < L) {
                    ++count;
                }
                else if (i > N - L) {
                    --count;
                }
                for (size_t j = 0; j < count; j++) {
                    double value = generateRandomNumber(minValue, maxValue);
                    if (i < L) {
                        if (2 * L - 1 - count + j == L - 1 && value == 0) {
                            value = 1;
                        }
                        matrix[i][2 * L - 1 - count + j] = value;
                        matrixCopy[i][2 * L - 1 - count + j] = value;
                        if (illConditionedMatrices)
                        {
                            double aValue = generateRandomNumber(minValue, maxValue);
                            accuracyMatrixLU[i][2 * L - 1 - count + j] = (aValue == 0)? 1: aValue;

                        }
                    }
                    else {
                        if (j == L - 1 && value == 0) {
                            value = 1;
                        }
                        matrix[i][j] = value;
                        matrixCopy[i][j] = value;
                        if (illConditionedMatrices)
                        {
                            double aValue = generateRandomNumber(minValue, maxValue);
                            accuracyMatrixLU[i][j] = (aValue == 0) ? 1 : aValue;
                        }
                    }
                }
            }
            matrix[N - 1][2 * L - 2] = 0;
            matrixCopy[N - 1][2 * L - 2] = 0;
            if (illConditionedMatrices)
            {
                accuracyMatrixLU[N - 1][2 * L - 2] = 0;
                for (size_t i = 0; i < N; i++)
                {
                    accuracyMatrixLU[i][L - 1] *= (std::pow(10, -k));
                }
                //PrintAccuracyLUMatrix();
                getAccuracyLUF();
            }
            getAccuracyF();
        }
        catch (const std::exception& e) {
            throw std::runtime_error( "Error: " + std::string(e.what()));
        }
    }

    TapeMatrix(const TapeMatrix& other) {
        matrix = other.matrix;
        matrixCopy = other.matrixCopy;
        N = other.N;
        L = other.L;
        solved = other.solved;
        x = other.x;
        f = other.f;
        q = other.q;
        accuracyX = other.accuracyX;
        solutionForAccuracyX = other.solutionForAccuracyX;
        accuracyF = other.accuracyF;
        meanRatioRelativeAccuracy = other.meanRatioRelativeAccuracy;
        
        illConditionedMatrices = other.illConditionedMatrices;
        if (other.illConditionedMatrices)
        {
            accuracyMatrixLU = other.accuracyMatrixLU;
            accuracyLUF = other.accuracyLUF;
            solutionForAccuracyLUX = other.solutionForAccuracyLUX;
            meanRatioRelativeAccuracyIllConditionedMatrices = other.meanRatioRelativeAccuracyIllConditionedMatrices;
        }
    }
    TapeMatrix& operator=(const TapeMatrix& other) {
        if (this == &other) {
            return *this;
        }

        matrix = other.matrix;
        matrixCopy = other.matrixCopy;
        N = other.N;
        L = other.L;
        solved = other.solved;
        x = other.x;
        f = other.f;
        q = other.q;
        accuracyX = other.accuracyX;
        solutionForAccuracyX = other.solutionForAccuracyX;
        accuracyF = other.accuracyF;
        meanRatioRelativeAccuracy = other.meanRatioRelativeAccuracy;

        illConditionedMatrices = other.illConditionedMatrices;
        if (other.illConditionedMatrices)
        {
            accuracyMatrixLU = other.accuracyMatrixLU;
            accuracyLUF = other.accuracyLUF;
            solutionForAccuracyLUX = other.solutionForAccuracyLUX;
            meanRatioRelativeAccuracyIllConditionedMatrices = other.meanRatioRelativeAccuracyIllConditionedMatrices;
        }

        return *this;
    }
    

    std::vector<std::vector<double>> getLUMatrixTape() const {
        return solved ? matrix : std::vector<std::vector<double>>();
    }
    std::vector<std::vector<double>> getaccuracyMatrixTapeLU() const {
        return illConditionedMatrices ? accuracyMatrixLU : std::vector<std::vector<double>>();
    }
    std::vector<std::vector<double>> getOriginalMatrixTape() const {
        return matrixCopy;
    }

    bool isSolved() { return solved; }

    int getN() const {
        return N;
    }
    int getL() const {
        return L;
    }
    double getQ() const {
        return q;
    }
    std::vector<double> getSolution() const {
        return solved ? x : std::vector<double>();
    }
    std::vector<double> getF() const {
        return f;
    }

    std::vector<double> getAccuracyX() const {
        return solved ? accuracyX : std::vector<double>();
    }
    
    std::vector<double> getSolutionForAccuracyX() const {
        return solved ? solutionForAccuracyX : std::vector<double>();
    }

    std::vector<double> getAccuracyF() const {
        return accuracyF;
    }
    std::vector<double> getAccuracyFLU() const {
        return f;
    }
    std::vector<double> getSolutionForAccuracyLUX() const {
        return illConditionedMatrices && solved? solutionForAccuracyLUX : std::vector<double>();
    }
    double getMeanRatioRelativeAccuracy() const {
        return meanRatioRelativeAccuracy;
    }
    double  getMeanRatioRelativeAccuracyIllConditionedMatrices() const {
        if (illConditionedMatrices)
        {
            return meanRatioRelativeAccuracyIllConditionedMatrices;
        }
        else {
            return -1;
        }
    }


    void solveSLAE() {

        if (getLUMatrix())
            getXSolution();

        if (checkSolution())
        {
            solved = true;
            getMeanRatioRelativeAccuracyBySolution();
        }
    }

private:
    bool getLUMatrix()
    {
        try
        {
            // Метод Халецкого для решения СЛАУ Ax=b
            if (matrix[0][L - 1] == 0)
            {
                return false;
            }
            for (size_t i = L; i < 2 * L - 1; i++)
            {
                matrix[0][i] = matrix[0][i] / matrix[0][L - 1];
            }
            for (int i = 1; i < N; i++)
            {
                //PrintMatrix();
                // col_v with w
                int newUpLine = i - 1;
                int newUpCol = L;
                int newLeftLine = i;
                int newLeftCol = L - 2;
                double sum = 0;

                int new_v = L - 2;
                for (int k = 0; k < L; k++)
                {
                    newUpLine = i - 1;
                    newUpCol = L;
                    newLeftLine = i + k;
                    newLeftCol = new_v;
                    sum = 0;
                    while ((newUpLine >= 0 && newUpLine < N) && newUpCol >= 0 && newUpCol <= 2 * L - 2 && newLeftLine >= 0 && newLeftLine < N && newLeftCol >= 0 && newLeftCol <= 2 * L - 2)
                    {
                        sum += (matrix[newUpLine][newUpCol] * matrix[newLeftLine][newLeftCol]);
                        newUpLine -= 1;
                        newUpCol += 1;
                        newLeftCol -= 1;
                    }
                    
                    if (i + k >= 0 && i + k < N && new_v + 1 >= 0 && new_v + 1 < 2 * L - 1)
                    {
                        matrix[i + k][new_v + 1] -= sum;
                    }
                    new_v -= 1;
                }

                //line_v
                new_v = L;
                for (size_t k = 0; k < L - 1; k++)
                {
                    newUpLine = i - 1;
                    newUpCol = new_v + 1;
                    newLeftLine = i;
                    newLeftCol = L - 2;
                    sum = 0;
                    while ((newUpLine >= 0 && newUpLine < N) && newUpCol >= 0 && newUpCol <= 2 * L - 2 && newLeftLine >= 0 && newLeftLine < N && newLeftCol >= 0 && newLeftCol <= 2 * L - 2)
                    {
                        sum += (matrix[newUpLine][newUpCol] * matrix[newLeftLine][newLeftCol]);
                        newUpLine -= 1;
                        newUpCol += 1;
                        newLeftCol -= 1;
                    }
                    if (matrix[i][L - 1] != 0)
                    {
                        matrix[i][new_v] = (matrix[i][new_v] - sum) / matrix[i][L - 1];
                    }
                    else
                    {
                        return false;
                    }
                    new_v += 1;
                }
            }
            return true;
        }
        catch (const std::exception& e)
        {
            throw std::runtime_error( "Error in getLUMatrix(): " + std::string(e.what()));
            return false;
        }
    }
    void getXSolution()
    {
        try
        {
            if (illConditionedMatrices)
            {
                // решение Ly=f
                std::vector<double> y;
                std::vector<double> accuracyY;
                std::vector<double> accuracyLUY;
                y.resize(N); y.reserve(N);
                accuracyY.resize(N); accuracyY.reserve(N);
                accuracyLUY.resize(N); accuracyLUY.reserve(N);
                for (int i = 0; i < N; ++i)
                {
                    double sum = 0;
                    double accuracySum = 0;
                    double accuracyLUSum = 0;
                    for (int j = 0; j < L - 1; j++)
                    {
                        if (i - j - 1 >= 0)
                        {
                            sum += y[i - j - 1] * matrix[i][L - j - 2];
                            accuracySum += accuracyY[i - j - 1] * matrix[i][L - j - 2];
                            accuracyLUSum += accuracyLUY[i - j - 1] * accuracyMatrixLU[i][L - j - 2];
                        }
                    }
                    y[i] = (f[i] - sum) / matrix[i][L - 1];
                    accuracyY[i] = (accuracyF[i] - accuracySum) / matrix[i][L - 1];
                    accuracyLUY[i] = (accuracyLUF[i] - accuracyLUSum) / accuracyMatrixLU[i][L - 1];
                    //std::cout << y[i]<<" , ";
                }
                // решение Ux=y

                for (int i = N - 1; i >= 0; --i)
                {
                    double sum = 0;
                    double accuracySum = 0;
                    double accuracyLUSum = 0;
                    for (int j = 0; j < L - 1; j++)
                    {
                        if (i + j + 1 < N)
                        {
                            sum += x[i + j + 1] * matrix[i][L + j];
                            accuracySum += solutionForAccuracyX[i + j + 1] * matrix[i][L + j];
                            accuracyLUSum += solutionForAccuracyLUX[i + j + 1] * accuracyMatrixLU[i][L + j];
                        }
                    }
                    x[i] = y[i] - sum;
                    solutionForAccuracyX[i] = accuracyY[i] - accuracySum;
                    solutionForAccuracyLUX[i] = accuracyLUY[i] - accuracyLUSum;
                }
                /*for (size_t i = 0; i < N; i++)
                {
                    std::cout << "XX" << i << ' ' << solutionForAccuracyLUX[i] << ' ';
                }*/
            }
            else
            {
                // решение Ly=f
                std::vector<double> y;
                std::vector<double> accuracyY;
                y.resize(N); y.reserve(N);
                accuracyY.resize(N); accuracyY.reserve(N);
                for (int i = 0; i < N; ++i)
                {
                    double sum = 0;
                    double accuracySum = 0;
                    for (int j = 0; j < L - 1; j++)
                    {
                        if (i - j - 1 >= 0)
                        {
                            sum += y[i - j - 1] * matrix[i][L - j - 2];
                            accuracySum += accuracyY[i - j - 1] * matrix[i][L - j - 2];
                        }
                    }
                    y[i] = (f[i] - sum) / matrix[i][L - 1];
                    accuracyY[i] = (accuracyF[i] - accuracySum) / matrix[i][L - 1];
                    //std::cout << y[i]<<" , ";
                }
                // решение Ux=y

                for (int i = N - 1; i >= 0; --i)
                {
                    double sum = 0;
                    double accuracySum = 0;
                    for (int j = 0; j < L - 1; j++)
                    {
                        if (i + j + 1 < N)
                        {
                            sum += x[i + j + 1] * matrix[i][L + j];
                            accuracySum += solutionForAccuracyX[i + j + 1] * matrix[i][L + j];
                        }
                    }
                    x[i] = y[i] - sum;
                    solutionForAccuracyX[i] = accuracyY[i] - accuracySum;
                }
            }
            
            
        }
        catch (const std::exception& e)
        {
            throw std::runtime_error( "Error in getXSolution(): " + std::string(e.what()));
        }
    }

    bool checkSolution()
    {
        try
        {
            bool check = true;
            double sum = 0;
            int count = L - 1;
            for (size_t i = 0; i < N && check; i++)
            {
                sum = 0;
                if (count < 2 * L - 1 && i < L)
                {
                    ++count;
                }
                else if (i > N - L && N != L)
                {
                    --count;
                }
                for (size_t j = 0; j < count && j < N; j++)
                {
                    if (i < L)
                    {
                        sum += x[j] * matrixCopy[i][2 * L - 1 - count + j];
                    }
                    else
                    {
                        sum += x[i - L + 1 + j] * matrixCopy[i][j];
                    }
                }
                if (f[i] - sum > q || f[i]-sum < q * (-1))
                {
                    check = false;
                }
            }
            return check;
        }
        catch (const std::exception& e)
        {
            throw std::runtime_error( "Error in checkSolution(): " + std::string(e.what()));
            return false;
        }
       
    }

    void getMeanRatioRelativeAccuracyBySolution()
    {
        double Er2 = 11;
        for (size_t i = 0; i < N; i++)
        {
            double er2 = (solutionForAccuracyX[i] - accuracyX[i]) < 0 ? (solutionForAccuracyX[i] * (-1) + accuracyX[i]) : (solutionForAccuracyX[i] - accuracyX[i]);
            if ((accuracyX[i] > q || (((-1) * accuracyX[i]) > q))&& accuracyX[i] != 0)
            {
                er2 /= accuracyX[i]<0? accuracyX[i]*(-1): accuracyX[i];
            }
            if (Er2 < er2 || Er2>10)
            {
                Er2 = er2;
            }
        }
        meanRatioRelativeAccuracy = Er2;
        if (illConditionedMatrices)
        {
            Er2 = 11;
            for (size_t i = 0; i < N; i++)
            {
                double erLU = (solutionForAccuracyLUX[i] - accuracyX[i]) < 0 ? (solutionForAccuracyLUX[i] * (-1) + accuracyX[i]) : (solutionForAccuracyLUX[i] - accuracyX[i]);
                if ((accuracyX[i] > q || (((-1) * accuracyX[i]) > q)) && accuracyX[i] != 0)
                {
                    erLU /= accuracyX[i] < 0 ? accuracyX[i] * (-1) : accuracyX[i];
                }
                if (Er2 < erLU || Er2>10)
                {
                    Er2 = erLU;
                }
            }
            meanRatioRelativeAccuracyIllConditionedMatrices = Er2;
        }
    }

    void getAccuracyF()
    {
        try
        {
            double sum = 0;
            int count = L - 1;
            for (size_t i = 0; i < N; i++)
            {
                sum = 0;
                if (count < 2 * L - 1 && i < L)
                {
                    ++count;
                }
                else if (i > N - L && N != L)
                {
                    --count;
                }
                for (size_t j = 0; j < count && j < N; j++)
                {
                    if (i < L)
                    {
                        sum += accuracyX[j] * matrixCopy[i][2 * L - 1 - count + j];
                    }
                    else
                    {
                        sum += accuracyX[i - L + 1 + j] * matrixCopy[i][j];
                    }
                }
                accuracyF[i] = sum;
            }
        }
        catch (const std::exception& e)
        {
            throw std::runtime_error( "Error in checkSolution(): " + std::string(e.what()));
        }

    }
    void getAccuracyLUF()
    {
        try
        {
            // решение Ux=y

            std::vector<double> y;
            std::vector<double> accuracyY;
            y.resize(N); y.reserve(N);
            accuracyY.resize(N); accuracyY.reserve(N);
            for (int i = N - 1; i >= 0; --i)
            {
                double accuracySum = 0;
                for (int j = 0; j < L - 1; j++)
                {
                    if (i + j + 1 < N)
                    {
                        accuracySum += accuracyX[i + j + 1] * accuracyMatrixLU[i][L + j];
                    }
                }
                accuracyY[i] = accuracyX[i] + accuracySum;
            }
            /*for (size_t i = 0; i < N; i++)
            {
                std::cout << "y" << i << ' ' << accuracyY[i]<<' ';
            }*/
            // решение Ly=f
            for (int i = 0; i < N; ++i)
            {
                double accuracySum = 0;
                for (int j = 0; j < L - 1; j++)
                {
                    if (i - j - 1 >= 0)
                    {
                        accuracySum += accuracyY[i - j - 1] * accuracyMatrixLU[i][L - j - 2];
                    }
                }
                accuracyLUF[i]= accuracyMatrixLU[i][L - 1]>0? accuracyY[i] * accuracyMatrixLU[i][L - 1] + accuracySum: (accuracyY[i] * accuracyMatrixLU[i][L - 1] + accuracySum);
            }
            /*for (size_t i = 0; i < N; i++)
            {
                std::cout << "f" << i << ' ' << accuracyLUF[i] << ' ';
            }
            for (size_t i = 0; i < N; i++)
            {
                std::cout << "x" << i << ' ' << accuracyX[i] << ' ';
            }*/
        }
        catch (const std::exception& e)
        {
            throw std::runtime_error( "Error in getXSolution(): " + std::string(e.what()));
        }
    }
};