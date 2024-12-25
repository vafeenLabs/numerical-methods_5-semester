#pragma once
#include "TestE.h"

#include <iostream>
#include <iomanip>
#include "TapeMatrix.h"

void testOne()
{
    const int numTests = 2;
    const int sizes[numTests] = { 10, 100 };
    const int ratios[numTests] = { 10, 10 };
    const int rangeMin = -10;
    const int rangeMax = 10;

    std::cout << std::endl << std::string(70, '-') << std::endl;
    std::cout << std::setw(5) << "Test" << std::setw(10) << "Size" << std::setw(17) << "L/N ratio" << std::setw(27) << "Mean relative accuracy" << std::endl;
    std::cout << std::string(70, '-') << std::endl;

    for (int i = 0; i < numTests; ++i) {
        int N = sizes[i];
        int L = N / ratios[i];

        TapeMatrix* lenta = new TapeMatrix(rangeMin, rangeMax, N, L);
        lenta->solveSLAE();
        while (!lenta->isSolved())
        {
            lenta = new TapeMatrix(rangeMin, rangeMax, N, L);
        }
        double E = lenta->getMeanRatioRelativeAccuracy();
        std::cout << std::setw(5) << i + 1 << std::setw(10) << N << std::setw(13) << "1/"<<ratios[i] << std::scientific << std::setprecision(2) << std::setw(20) << E << std::endl;
        
    }
}

void testTwo()
{
    const int numTests = 2;
    const int countTests = 2;
    const int sizes[numTests] = { 10, 100 };
    const int rangeMin = -10;
    const int rangeMax = 10;

    std::cout << std::endl << std::string(70, '-') << std::endl;
    std::cout << std::setw(5) << "Test" << std::setw(10) << "Size" << std::setw(25) << "Mean relative accuracy" << std::endl;
    std::cout << std::string(70, '-') << std::endl;
    int testNumber = 1;
    for (int i = 0; i < numTests; ++i) {
        for (int j = 0; j < countTests; ++j, ++testNumber)
        {
            int N = sizes[i];
            int L = N;

            TapeMatrix* lenta = new TapeMatrix(rangeMin, rangeMax, N, L);
            lenta->solveSLAE();
            while (!lenta->isSolved())
            {
                lenta = new TapeMatrix(rangeMin, rangeMax, N, L);
            }
            double E = lenta->getMeanRatioRelativeAccuracy();
            std::cout << std::setw(5) << testNumber << std::setw(10) << N << std::setw(15) << std::scientific << std::setprecision(2) << E << std::endl;
        }
    }
}

void testThree()
{
    const int numK = 3;
    const int sizes[1] = { 10 };
    const int _k[numK] = { 2,3,4 };
    const int rangeMin = -10;
    const int rangeMax = 10;

    std::cout << std::endl << std::string(70, '-') << std::endl;
    std::cout << std::setw(5) << "Test" << std::setw(10) << "Size" << std::setw(17) << "Grade k" << std::setw(27) << "Mean relative accuracy" << std::endl;
    std::cout << std::string(70, '-') << std::endl;

    
        for (int j = 0; j < numK; ++j)
        {
            int N = sizes[0];
            int L = N;
            int k = _k[j];

            TapeMatrix* lenta = new TapeMatrix(rangeMin, rangeMax, N, L, true, k);
            lenta->solveSLAE();
            while (!lenta->isSolved())
            {
                lenta = new TapeMatrix(rangeMin, rangeMax, N, L, true, k);
            }
            double E = lenta->getMeanRatioRelativeAccuracyIllConditionedMatrices();
            std::cout << std::setw(5) << j+1 << std::setw(10) << N << std::setw(13) << k << std::scientific << std::setprecision(2) << std::setw(20) << E << std::endl;
        }
    
}