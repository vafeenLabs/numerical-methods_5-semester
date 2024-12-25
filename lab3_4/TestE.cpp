#pragma once
#include "TestE.h"
#include "GeneratorSymmetricMatrixWithEigenVectorsAndValues.h"
#include "MethodInverseIterations.h"

#include <iostream>
#include <iomanip>

void testOne()
{
    const int averageTests = 10;

    const int numTests = 12;
    const int sizes[] = { 10, 30, 50 };
    const int lambdaRange[] = { 2, 50 };
    const int eigenVectorsValuesERange[] = {-5, -8};

    std::cout << std::endl << std::string(150, '-') << std::endl;
    std::cout << std::setw(5) << "Test" 
              << std::setw(10) << "Size" 
              << std::setw(18) << "Lambda range" 
              << std::setw(10) << "E" 
              << std::setw(24) << "Average E val" 
              << std::setw(25) << "Average E vec" 
              << std::setw(25) << "Average measure r" 
              << std::setw(31) << "Average iterations number" 
              << std::endl;
    std::cout << std::string(150, '-') << std::endl;

    int counterE = -1;
    int counterLambdaRange = 0;
    int counterSizes = 0;
    int N = sizes[0];
    int lambda = lambdaRange[0];
    int E;

    int testNumber = 1;
    for (int counterSizes = 0; counterSizes < sizeof(sizes) / sizeof(sizes[0]); ++counterSizes) {
        N = sizes[counterSizes];

        for (int counterLambdaRange = 0; counterLambdaRange < sizeof(lambdaRange) / sizeof(lambdaRange[0]); ++counterLambdaRange) {
            lambda = lambdaRange[counterLambdaRange];

            for (int counterE = 0; counterE < sizeof(eigenVectorsValuesERange) / sizeof(eigenVectorsValuesERange[0]); ++counterE) {
                E = eigenVectorsValuesERange[counterE];

                double averageEVec = 0;
                double averageEVal = 0;
                double averageR = 0;
                int averageIterationsNumber = 0;
                for (size_t i = 0; i < averageTests; ++i)
                {

                    GeneratorSymmetricMatrixWithEigenVectorsAndValues* generator = new GeneratorSymmetricMatrixWithEigenVectorsAndValues(N, -100, 100, -lambda, lambda);
                    std::vector<std::vector<double>> symmetricMatrix = generator->getSymmetricMatrix();
                    
                    

                    MethodInverseIterations* finder = new MethodInverseIterations(N, symmetricMatrix, std::pow(10, E), std::pow(10, E), 1000);
                    finder->Solve();
                    double Evec = finder->getResultedEigenVectorsE();
                    double Eval = finder->getResultedEigenValuesE();
                    double r = finder->getR();
                    int IterationsNumber = finder->getIterationsNumber();

                    averageEVec += Evec;
                    averageEVal += Eval;
                    averageR += r;
                    averageIterationsNumber += IterationsNumber;

                    if (i == averageTests-1)
                    {
                        /*std::cout << "-------vector----------" << std::endl;
                        printArr(finder->getEigenVectorByFirstMinEigenValue(), N);
                        std::vector<std::vector<double>> eigenVectors = generator->getEigenVectorsData();
                        std::vector<std::vector<double>> IeigenVectors = generator->getInverseveEigenVectorsData();

                        std::cout << "-------eigenVectors----------" << std::endl;
                        printArr(eigenVectors, N, N);*/

                        /*
                        std::cout << "-------symmetricMatrix----------" << std::endl;
                        printArr(symmetricMatrix, N, N);
                        std::vector<double> eigenValues = generator->getEigenValuesData();
                        std::cout << "------eigenValues-----------" << std::endl;
                        printArr(eigenValues, N);
                        std::cout << "-----------------" << std::endl;
                        
                        std::cout << "-------min-value-------" << std::endl;
                        std::cout << finder->getFirstMinEigenValue() << std::endl;
                        std::cout << "----------r----------" << std::endl;
                        std::cout << std::scientific << finder->getR() << std::endl;
                        std::cout << "----------E-vec----------" << std::endl;
                        std::cout << std::scientific << finder->getResultedEigenVectorsE() << std::endl;
                        std::cout << "----------E-val---------" << std::endl;
                        std::cout << std::scientific << finder->getResultedEigenValuesE() << std::endl;
                        std::cout << "-----------------" << std::endl;
                        */
                        averageEVec /= averageTests;
                        averageEVal /= averageTests;
                        averageR /= averageTests;
                        averageIterationsNumber /= averageTests;
                        std::vector<double> firstEigenVector = finder->getEigenVectorByFirstMinEigenValue();
                        double minFirstEigenValue = finder->getFirstMinEigenValue();
                        std::cout << std::setw(5) << testNumber
                                  << std::setw(10) << N 
                                  << std::setw(12) << -lambda<<':'<< lambda 
                                  << std::scientific << std::setprecision(1) << std::setw(17) << std::pow(10, E)
                                  << std::scientific << std::setprecision(3) << std::setw(20) << averageEVal
                                  << std::scientific << std::setprecision(3) << std::setw(25) << averageEVec
                                  << std::scientific << std::setprecision(3) << std::setw(25) << averageR
                                  << std::setw(23) << averageIterationsNumber
                                  << std::endl;
                        ++testNumber;
                    }
                }
            }
        }
    }

    
}
