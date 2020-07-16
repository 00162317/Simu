#include <iostream>
#include "math_tools.h"
#include "classes.h"
#include "tools.h"
#include "sel.h"

int main(int argc, char *argv[])
{
    char filename[150];
    strcpy(filename,argv[1]);

    vector<Matrix> localKs;
    vector<Vector> localbs;
    vector<Vector> Ts;

    Matrix K;
    Vector b;

    cout << "IMPLEMENTACI"<<char(224)<<"N DEL M"<<char(144)<<"TODO DE LOS ELEMENTOS FINITOS\n"
         << "\t- NAVIER-STOKES" << "\t- 3 DIMENSIONES\n"
         << "\t- FUNCIONES DE FORMA LINEALES\n" << "\t- PESOS DE GALERKIN\n"
         << "*********************************************************************************\n\n";

    mesh m;
    leerMallayCondiciones(m,filename);
    cout << "Datos obtenidos correctamente\n********************\n";

    crearSistemasLocales(m,localKs,localbs);
    cout << "Sistemas locales creados correctamente\n********************\n";

    zeroes(K,4*m.getSize(NODES));
    zeroes(b,4*m.getSize(NODES));
    cout << "Comenzando ensamblaje\n********************\n";
    ensamblaje(m,localKs,localbs,K,b);

    cout << "Aplicando condiciones\n********************\n";
    applyDirichlet(m,K,b);

    showMatrix(K);

    showVector(b);

	cout << "Calculando valores\n********************\n";
    calculate(K,b,Ts,4*m.getSize(NODES)-m.getSize(DIRICHLET));

    cout << "Escribiendo resultados en archivo .post.res\n********************\n";
    writeResults(m,Ts,filename);

    delete[] &m;

    return 0;

}
