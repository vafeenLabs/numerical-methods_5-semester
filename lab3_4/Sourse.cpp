#pragma once

//#include "GeneratorSymmetricMatrixWithEigenVectorsAndValues.h"
//#include "MethodInverseIterations.h"

#include "TestE.h"
#include <iostream>
#include <Windows.h>


int main()
{
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
    /*
    int size = 10;
    GeneratorSymmetricMatrixWithEigenVectorsAndValues* generator = new GeneratorSymmetricMatrixWithEigenVectorsAndValues(size, 1, 10);
    std::vector<std::vector<double>> symmetricMatrix = generator->getSymmetricMatrix();
    std::vector<std::vector<double>> eigenVectors = generator->getEigenVectorsData();
    std::vector<std::vector<double>> IeigenVectors = generator->getInverseveEigenVectorsData();
    std::vector<double> eigenValues = generator->getEigenValuesData();

    std::cout <<"-------symmetricMatrix----------" << std::endl;
    printArr(symmetricMatrix, size, size);
    std::cout <<"-------eigenVectors----------" << std::endl;
    printArr(eigenVectors, size, size);
    std::cout <<"------eigenValues-----------" << std::endl;
    printArr(eigenValues, size);
    std::cout <<"-----------------" << std::endl;
    printMatrix(IeigenVectors);
    std::cout <<"-----------------" << std::endl;
    
    MethodInverseIterations* finder = new MethodInverseIterations(size, symmetricMatrix, 0.0000000001, 0.0000000001, 100);
    
    //symmetricMatrix = {
    //    {1,3},
    //    {3,1}
    //};
    //MethodInverseIterations* finder = new MethodInverseIterations(size, symmetricMatrix, 0.0000001, 0.0000001, 100);
    

    finder->Solve();

    std::cout <<"-------vector----------" << std::endl;
    printArr(finder->getEigenVectorByFirstMinEigenValue(), size);
    std::cout << "-------min-value-------" << std::endl;
    std::cout << finder->getFirstMinEigenValue() << std::endl;
    std::cout <<"----------r----------" << std::endl;
    std::cout << std::scientific << finder->getR() << std::endl;
    std::cout <<"----------E-vec----------" << std::endl;
    std::cout << std::scientific<< finder->getResultedEigenVectorsE() << std::endl;
    std::cout <<"----------E-val---------" << std::endl;
    std::cout << std::scientific<< finder->getResultedEigenValuesE() << std::endl;
    std::cout <<"-----------------" << std::endl;
    */


    testOne();

    return 0;
}


