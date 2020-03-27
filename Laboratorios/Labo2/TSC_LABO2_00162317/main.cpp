#include <fstream>
#include <iostream>
#include "classes.h"

using namespace std;
//

void readCoordinates(ifstream &file, int n, item* item_list){
    int e; float r;
    for(int i=0; i<n; i++){
        file >> e >> r;
        item_list[i].setIntFloat(e,r);
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

    m.setParameters(k, Q);
    m.setSizes(nnodes, neltos, ndirich, nneu);
    m.createData();

    readCoordinates(file,nnodes,m.getNodes());

    file.close();

    for(int i = 0;i<10;++i){

        cout << m.getNode(i).getX() << endl;
    }
    return 0;
}