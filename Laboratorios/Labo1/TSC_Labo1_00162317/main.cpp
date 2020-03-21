#include <iostream>
#include "math_tools.h"
#include "display_tools.h"
using namespace std;
int main() {

    Matrix MatrixA, MatrixB;
    zeroes(MatrixA, 3);
    //////////////////////////////////////
    MatrixA.at(0).at(0) = 5;
    MatrixA.at(0).at(1) = 8;
    MatrixA.at(0).at(2) = 2;
    //////////////////////////////////////
    MatrixA.at(1).at(0) = 9;
    MatrixA.at(1).at(1) = 1;
    MatrixA.at(1).at(2) = 7;
    //////////////////////////////////////
    MatrixA.at(2).at(0) = 2;
    MatrixA.at(2).at(1) = 6;
    MatrixA.at(2).at(2) = 1;
    //////////////////////////////////////
    showMatrix(MatrixInversa(MatrixA, MatrixB));


    return 0;
}

