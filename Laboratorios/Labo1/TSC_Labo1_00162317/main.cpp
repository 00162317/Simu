#include <iostream>
#include "math_tools.h"
#include "display_tools.h"
using namespace std;
int main() {

    Matrix MatrixA, MatrixB;
    zeroes(MatrixA, 3);
    //////////////////////////////////////
    MatrixA.at(0).at(0) = 2;
    MatrixA.at(0).at(1) = 2;
    MatrixA.at(0).at(2) = 3;
    //////////////////////////////////////
    MatrixA.at(1).at(0) = 4;
    MatrixA.at(1).at(1) = 5;
    MatrixA.at(1).at(2) = 6;
    //////////////////////////////////////
    MatrixA.at(2).at(0) = 7;
    MatrixA.at(2).at(1) = 8;
    MatrixA.at(2).at(2) = 9;
    //////////////////////////////////////
    showMatrix(MatrixInversa(MatrixA, MatrixB));


    return 0;
}

