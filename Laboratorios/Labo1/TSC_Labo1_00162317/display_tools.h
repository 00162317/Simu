//
// Created by Roberto Hernandez on 20/3/2020.
//

#include <iostream>
#include "math_tools.h"

#ifndef TSC_LABO1_00162317_DISPLAY_TOOLS_H
#define TSC_LABO1_00162317_DISPLAY_TOOLS_H

#endif //TSC_LABO1_00162317_DISPLAY_TOOLS_H

using namespace std;


void showVector(Vector b){
    cout << "[\t";
    for(int i=0;i<b.size();i++){
        cout << b.at(i) << "\t";
    }
    cout << "]\n";
}

void showMatrix(Matrix K){
    for(int i=0;i<K.at(0).size();i++){
        cout << "[\t";
        for(int j=0;j<K.size();j++){
            cout << K.at(i).at(j) << "\t";
        }
        cout << "]\n";
    }
}
