#include <iostream>
#include "Matrix.h"
#include"Utils.h"


int main() {

    /*
    */
    matrixCase* arr[] = {
        new matrixCase(10,-10,10),
        new matrixCase(10,-100,100),
        new matrixCase(10,-1000,1000),
        new matrixCase(100,-10,10),
        new matrixCase(100,-100,100),
        new matrixCase(100,-1000,1000),
        new matrixCase(1000,-10,10),
        new matrixCase(1000,-100,100),
        new matrixCase(1000,-1000,1000),
    };
    getTestResults(arr, sizeof(arr) / sizeof(arr[0]), 10);

    /*Matrix* A = new Matrix(10, "13.txt");
    A->solve();
    double* x = A->getSolution();
    double E1 = A->getE1();
    double E2 = A->getE2();
    printArray(x, 10);
    std::cout << E1 << ' ' << E2;*/


    return 0;
}
