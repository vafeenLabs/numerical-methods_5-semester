#pragma once
#include "TapeMatrix.h"
#include "TestE.h"
#include "functions.h"

#include <Windows.h>

int main() {
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);

    
    /*
    int N = 10; // Размер обычной матрицы
    int L = 3; // Половина ширины ленты
    */

    //std::string filename = "matrix.txt";
    //TapeMatrix lenta(filename, N, L);

    /*
    TapeMatrix lenta(-10, 10, N, L);
    
    printArr(lenta.getOriginalMatrixTape(), lenta.getN(), 2*lenta.getL()-1);
    lenta.solveSLAE();
    if (lenta.isSolved())
    {
        printArr(lenta.getLUMatrixTape(), lenta.getN(), 2*lenta.getL()-1);
        printArr(lenta.getSolution(), lenta.getN());
        printArr(lenta.getAccuracyX(), lenta.getN());
        printArr(lenta.getSolutionForAccuracyX(), lenta.getN());
        double E = lenta.getMeanRatioRelativeAccuracy();
        std::cout << std::scientific << "E: " << std::setprecision(2) << E << std::endl;
    }

    TapeMatrix lentaTest3(-10, 10, N, L, true, 2);

    printArr(lentaTest3.getOriginalMatrixTape(), lentaTest3.getN(), 2 * lentaTest3.getL() - 1);
    printArr(lentaTest3.getaccuracyMatrixTapeLU(), lentaTest3.getN(), 2 * lentaTest3.getL() - 1);
    lentaTest3.solveSLAE();
    if (lentaTest3.isSolved())
    {
        //printArr(lentaTest3.getLUMatrixTape(), lentaTest3.getN(), 2 * lentaTest3.getL() - 1);
        //printArr(lentaTest3.getSolution(), lentaTest3.getN());
        printArr(lentaTest3.getSolutionForAccuracyLUX(), lentaTest3.getN());
        printArr(lentaTest3.getAccuracyX(), lentaTest3.getN());
        double E = lentaTest3.getMeanRatioRelativeAccuracyIllConditionedMatrices();
        std::cout << std::scientific << "E: " << std::setprecision(2) << E << std::endl;
    }
    */
    
    testOne();
    testTwo();
    testThree();

    return 0;
}
