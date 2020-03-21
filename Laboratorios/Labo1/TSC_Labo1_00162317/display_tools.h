//
// Created by Roberto Hernandez on 20/3/2020.
//
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

void showMatrix(Matrix k){
    for (int i = 0; i < k.at(0).size() ; i++) {
        cout << "[\t";
        for (int j = 0; j < k.size() ; j++) {
            cout << k.at(i).at(j) << "\t";
        }
        cout << "]\n";
    }
}