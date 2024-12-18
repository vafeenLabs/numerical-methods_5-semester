#pragma once

#include "TapeMatrix.h"

#include <iostream>
#include <vector>
#include <random>

class MethodInverseIterations{
private:
    // in
    
    int _size;
    std::vector<std::vector<double>> _symmetricMatrix;
    
    double _givenEigenVectorsE;
    double _givenEigenValuesE;
    int _maxIterationsNumber;
    
    // out

    double _firstMinEigenValue;
    std::vector<double> _eigenVectorByFirstMinEigenValue;
    int _IterationsNumber;
    double _resultedEigenVectorsE=0;
    double _resultedEigenValuesE=0;
    double r=10;
public:

    std::vector<double> getEigenVectorByFirstMinEigenValue() { return _eigenVectorByFirstMinEigenValue; }
    double getFirstMinEigenValue() { return _firstMinEigenValue; }
    double getIterationsNumber() { return _IterationsNumber; }
    double getR() { return r; }
    double getResultedEigenVectorsE() { return _resultedEigenVectorsE; }
    double getResultedEigenValuesE() { return _resultedEigenValuesE; }

    MethodInverseIterations(int size,
        std::vector<std::vector<double>> symmetricMatrix,
        double eigenVectorsE,
        double eigenValuesE,
        int maxIterationsNumber):
        _size(size),

        _givenEigenValuesE(eigenValuesE), 
        _givenEigenVectorsE(eigenVectorsE),
        
        _maxIterationsNumber(maxIterationsNumber),

        _eigenVectorByFirstMinEigenValue(size),

        _symmetricMatrix(size, std::vector<double>(size))
    {
        _symmetricMatrix = symmetricMatrix;
    }

    void Solve()
    {
        std::vector<double> x_rand(_size);
        for (size_t i = 0; i < _size; i++)
        {
            x_rand[i] = generateRandomNumber(-10, 10);
        }
        int k = 0;
        double q = 10;
        double qPrev = 0;
        double maxVecE = 10;
        TapeMatrix* system = new TapeMatrix(_symmetricMatrix, _size, _size);
        while ((( _resultedEigenValuesE > _givenEigenValuesE) || (maxVecE > _givenEigenVectorsE)) && (k <  _maxIterationsNumber))
        {
            std::vector<double> v = normalizeVector(x_rand);

            //std::cout << "Vk:\n";
            //printArr(v, _size);
            system->solveSLAE(v);
            if (system->isSolved())
            {
                x_rand = system->getSolution();
                //std::cout << "x k+1:\n";
                //printArr(x_rand, _size);
            
                qPrev = q;
                q = 0;
                for (int i = 0; i < _size; ++i)
                {
                    q += v[i] * x_rand[i];
                }
                _firstMinEigenValue = 1 / q;
                //std::cout << "2 q:\n"<< q<<'\n';
                //std::cout << "value:\n"<< _firstMinEigenValue << '\n';
                maxVecE = -10;
                for (size_t i = 0; i < _size && !_eigenVectorByFirstMinEigenValue.empty(); i++)
                {
                    double tmp = std::abs(std::abs(v[i]) - std::abs(_eigenVectorByFirstMinEigenValue[i]));
                    if (maxVecE < tmp)
                    {
                        maxVecE = tmp;
                    }
                }
                _resultedEigenVectorsE = maxVecE;
                _resultedEigenValuesE = std::abs(std::abs(q) - std::abs(qPrev));
                _eigenVectorByFirstMinEigenValue = v;
            }
            ++k;
        }
        _IterationsNumber = k;
        std::vector<double> _R(_size);
        for (size_t i = 0; i < _size; i++)
        {
            for (size_t j = 0; j < _size; j++)
            {
                _R[i] += _symmetricMatrix[i][j] * _eigenVectorByFirstMinEigenValue[j];
            }
            _R[i] -= _eigenVectorByFirstMinEigenValue[i] * _firstMinEigenValue;
            if (_R[i]<r)
            {
                r = _R[i];
            }
        }
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
