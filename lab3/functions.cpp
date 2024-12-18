#pragma once
#include "functions.h"
#include <iostream>
#include <random>
#include <vector>
#include <iomanip>
#include <stdexcept>

double generateRandomNumber(double min_val, double max_val) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(min_val, max_val);
    return dis(gen);
}
double roundError(double error) {
    // Находим порядок погрешности
    int power = std::floor(std::log10(std::abs(error)));
    // Округляем погрешность до 3 значащих цифр согласно условию
    double roundedError = std::round(error / std::pow(10, power - 2)) * std::pow(10, power - 2);
    return roundedError;
}


void printArr(const std::vector<std::vector<double>>& matrix, int rows, int cols) {
    try
    {
        std::cout << std::endl << std::string(1, ' ') << std::string(cols*10+cols, '-') << std::endl;
        for (int i = 0; i < rows; ++i) {
            std::cout  << '|';
            for (int j = 0; j < cols; ++j) {
                std::cout << std::fixed << std::setprecision(5) << std::setw(10) << matrix.at(i).at(j) << " ";
            }
            std::cout  << '|' << std::endl;
        }
        std::cout << std::string(1, ' ') << std::string(cols * 10+cols, '-') << std::endl;
    }
    catch (const std::exception& e) {
        throw std::runtime_error("Print matrix error: " + std::string(e.what()));
    }
}

void printArr(const std::vector<double>& vector, int size) {
    try
    {
        std::cout << std::endl << std::string(20, '-') << std::endl;
        for (int i = 0; i < size; ++i) {
            std::cout << std::fixed << std::setprecision(5) << '|'
                      << std::setw(2) << "[" << i << "] = "
                      << std::setw(10) << std::setfill(' ') << vector.at(i)
                      << std::setw(2) << '|' << std::endl;
        }
        std::cout << std::endl << std::string(20, '-') << std::endl;
    }
    catch (const std::exception& e) {
        throw std::runtime_error("Print array error: " + std::string(e.what()));
    }
    std::cout << std::endl;
}
