#pragma once
#include <iostream>
#include "Matrix.h"
#include <iomanip>



double roundError(double error) {
    // Находим порядок погрешности
    int power = std::floor(std::log10(std::abs(error)));
    // Округляем погрешность до 3 значащих цифр согласно условию
    double roundedError = std::round(error / std::pow(10, power - 2)) * std::pow(10, power - 2);
    return roundedError;
}
void printArray(double arr[], int size)
{
    std::cout << std::endl;
    for (int j = 0; j < size; j++) {
        std::cout << '[' << j + 1 << ']' << '=' << arr[j] << std::endl;
    }
    std::cout << std::endl;
}

struct matrixCase
{
    size_t size;
    double min;
    double max;

    matrixCase(size_t size, double min, double max)
    {
        this->size = size;
        this->max = max;
        this->min = min;
    }

    matrixCase(matrixCase& other)
    {
        size = other.size;
        min = other.min;
        max = other.max;
    }
};

void getTestResults(matrixCase* arr[], size_t size, size_t testNumberForAccuracy = 10)
{
    double MainСounter = 0;
    double failureСounter = 0;

    std::ofstream outputFile("results.txt");

    outputFile << std::endl << std::setfill('-') << std::setw(15) << "+"
        << std::setw(25) << "+"
        << std::setw(35) << "+"
        << std::setw(40) << "+"
        << std::setw(35) << "+" << std::setfill(' ') << std::endl;
    outputFile << std::endl << std::left << std::setw(15) << "Test #"
        << std::left << std::setw(25) << "Dimensionality system"
        << std::left << std::setw(35) << "Range of values elements matrix"
        << std::left << std::setw(40) << "Average relative error of the system"
        << std::left << std::setw(35) << "Average value of accuracy estimation" << std::endl;


    for (size_t i = 0; i < size; i++)
    {
        double sumE1 = 0;
        double sumE2 = 0;
        size_t j = 0;

        while (j < testNumberForAccuracy)
        {
            bool solved = false;

            while (!solved)
            {
                MainСounter++;
                Matrix* A = new Matrix(arr[i]->size, arr[i]->min, arr[i]->max);
                A->solve();
                solved = A->isSolved();

                if (solved)
                {
                    //double* solution = A->getSolution();
                    //printArray(solution, arr[i]->size);
                    sumE1 += A->getE1();
                    sumE2 += A->getE2();
                    //std::cout << "E1" << sumE2;
                    ++j;
                }
                else
                {
                    failureСounter++;
                }
            }
        }

        outputFile << std::endl << std::setfill('-') << std::setw(15) << "+"
            << std::setw(25) << "+"
            << std::setw(35) << "+"
            << std::setw(40) << "+"
            << std::setw(35) << "+" << std::setfill(' ') << std::endl;
        outputFile << std::endl << std::left << std::setw(15) << i + 1
            << std::left << std::setw(25) << arr[i]->size
            << std::left << "[" << (arr[i]->min) << ":" << (arr[i]->max) << std::setw(35) << "]"
            << std::left << std::setw(40) << roundError(sumE2 / testNumberForAccuracy)
            << std::left << std::setw(35) << roundError(sumE1 / testNumberForAccuracy) << std::endl;
    }


    outputFile << std::endl << std::setfill('-') << std::setw(15) << "+"
        << std::setw(25) << "+"
        << std::setw(35) << "+"
        << std::setw(40) << "+"
        << std::setw(35) << "+" << std::setfill(' ') << std::endl;

    outputFile << std::endl << std::left << std::setw(15) << "Failure rate"
        << std::left << std::round((failureСounter / MainСounter) * 100) << std::setw(25) << '%' << std::endl;

    outputFile.close();
}


/*
void getTestResults(matrixCase* arr[], size_t size, size_t testNumberForAccuracy = 10)
{
    double MainСounter = 0;
    double failureСounter = 0;
    for (size_t i = 0; i < size; i++)
    {
        double sumE1 = 0;
        double sumE2 = 0;
        size_t j = 0;
        while (j < testNumberForAccuracy)
        {
            bool solved = false;
            while (!solved)
            {
                MainСounter++;
                Matrix* A = new Matrix(arr[i]->size, arr[i]->min, arr[i]->max);
                A->solve();
                solved = A->isSolved();
                if (solved)
                {
                    double* solution = A->getSolution();
                    sumE1 += A->getE1();
                    sumE2 += A->getE2();
                    ++j;
                }
                else
                {
                    failureСounter++;
                }
            }
        }
        std::cout << std::endl << "Table row " << i + 1 << ':' << std::endl;
        std::cout << "Dimensionality system: " << arr[i]->size << std::endl;
        std::cout << "Range of values elements matrix : Min:" << arr[i]->min << " Max: " << arr[i]->max << std::endl;
        std::cout << "Average relative error of the system: " << roundError(sumE2 / testNumberForAccuracy) << std::endl;
        std::cout << "Average value of accuracy estimation: " << roundError(sumE1 / testNumberForAccuracy) << std::endl;
    }
    std::cout << std::endl << "Failure rate: " << std::round((failureСounter / MainСounter) * 100) << '%' << std::endl;
}
*/

