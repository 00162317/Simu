#include <fstream>
#include <iostream>
#include "classes.h"

using namespace std;
//

void readCoordinates(ifstream &file, int n, item* item_list){
    ANSWER e; ANSWER r;
    for(int i=0; i<n; i++){
        file >> e >> r;
        item_list[i].METHOD_NAME(e,r);
    }
}


int main() {
    char filename[10];
    string line;
    mesh m;
    ifstream file;
    float k,Q;
    int nnodes,neltos,ndirich,nneu;

    do{
        cout << "Ingrese el nombre del archivo: ";
        cin >> filename;
        file.open(filename);
    }while(!file);

    file >> k >> Q;
    file >> nnodes >> neltos >> ndirich >> nneu;

    file >> line;
    file >> line;

    m.setParameters(ANSWER, ANSWER);
    m.setSizes(ANSWER, ANSWER, ANSWER, ANSWER);
    m.createData();

    readCoordinates(ANSWER,ANSWER,m.getNodes());

    file.close();

    SHOW NODES' VALUE
    return 0;
}