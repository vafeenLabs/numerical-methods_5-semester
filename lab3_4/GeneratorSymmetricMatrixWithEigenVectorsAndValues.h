#pragma once


#include <iostream>
#include <vector>
#include <random>
//#include <algorithm>

class GeneratorSymmetricMatrixWithEigenVectorsAndValues {
private:
    double _size;
    std::vector<std::vector<double>> _eigenVectorsData;
    std::vector<std::vector<double>> _inverseEigenVectorsData;
    std::vector<double> _eigenValuesData;

    std::vector<std::vector<double>> _symmetricMatrix;
public:
    GeneratorSymmetricMatrixWithEigenVectorsAndValues(int size, double min, double max, double lambdaMin=0, double lambdaMax=0) :
        _size(size), 
        _eigenVectorsData(size, std::vector<double>(size)), 
        _eigenValuesData(size), 
        _inverseEigenVectorsData(size, std::vector<double>(size)),
        _symmetricMatrix(size, std::vector<double>(size))
    {
        lambdaMin = lambdaMin == 0 && lambdaMin != min ? min : lambdaMin;
        lambdaMax = lambdaMax == 0 && lambdaMax != max ? max : lambdaMax;
        std::random_device rd;
        std::mt19937 gen(rd());
        //std::uniform_int_distribution<> dis(min, max);
        //std::uniform_int_distribution<> lambdaDis(lambdaMin, lambdaMax);
        std::uniform_real_distribution<> dis(min, max);
        std::uniform_real_distribution<> lambdaDis(lambdaMin, lambdaMax);

        for (int i = 0; i < _size; ++i) {
            _eigenValuesData[i] = lambdaDis(gen);
            for (int j = 0; j < _size; ++j) {
                _eigenVectorsData[i][j] = dis(gen);
            }
        }
        
        for (int i = 0; i < _size; ++i)
        {
            for (int j = 0; j < _size; ++j)
            {
                _symmetricMatrix[i][j] = _eigenVectorsData[i][j] * _eigenValuesData[j];
            }
        }

        _inverseEigenVectorsData = inverseMatrix(_eigenVectorsData);

        _symmetricMatrix = multiply(_symmetricMatrix, _inverseEigenVectorsData, _size);
        
    }

    std::vector<std::vector<double>> getSymmetricMatrix()
    {
        return _symmetricMatrix;
    }
    std::vector<std::vector<double>> getEigenVectorsData()
    {
        return _eigenVectorsData;
    }
    std::vector<std::vector<double>> getInverseveEigenVectorsData()
    {
        return _inverseEigenVectorsData;
    }
    std::vector<double> getEigenValuesData()
    {
        return _eigenValuesData;
    }

private:
    std::vector<std::vector<double>> multiply(const std::vector<std::vector<double>>& matrix1, const std::vector<std::vector<double>>& matrix2,const int &size) {
        

        std::vector<std::vector<double>> result(size, std::vector<double>(size, 0.0));

        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                for (int k = 0; k < size; ++k) {
                    result[i][j] += matrix1[i][k] * matrix2[k][j];
                }
            }
        }

        return result;
    }


    std::vector<std::vector<double>> inverseMatrix(const std::vector<std::vector<double>>& A) {
        int n = A.size();

        //  метод Гаусса-Жордана
        //  [A | I] I - единичная матрица
        std::vector<std::vector<double>> augmentedMatrix(n, std::vector<double>(2 * n, 0.0));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                augmentedMatrix[i][j] = A[i][j];
            }
            augmentedMatrix[i][i + n] = 1.0;
        }

        for (int i = 0; i < n; ++i) {
            double pivot = augmentedMatrix[i][i];
            for (int j = 0; j < 2 * n; ++j) {
                augmentedMatrix[i][j] /= pivot;
            }
            for (int j = 0; j < n; ++j) {
                if (j != i) {
                    double factor = augmentedMatrix[j][i];
                    for (int k = 0; k < 2 * n; ++k) {
                        augmentedMatrix[j][k] -= factor * augmentedMatrix[i][k];
                    }
                }
            }
        }
        std::vector<std::vector<double>> inverse(n, std::vector<double>(n, 0.0));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                inverse[i][j] = augmentedMatrix[i][j + n];
            }
        }

        return inverse;
    }

    std::vector<double> normalizeVector(const std::vector<double>& vector) {
        double sum = 0.0;
        for (double element : vector) {
            sum += element * element;
        }

        double magnitude = std::sqrt(sum);

        std::vector<double> normalizedVector;
        for (double element : vector) {
            normalizedVector.push_back(element / magnitude);
        }

        return normalizedVector;
    }
};
