# Laboratorio 3
Técnicas de Simulación en Computadoras: Tercera práctica de laboratorio 

## Contenido

<p align="center"> Método de los elementos finitos para el problema de transferencia de calor
en una dimensión con funciones de forma lineales y con peso de Galerkin </p>

## Cambios en math_tools.h

```cpp
//Agregada la función para calcular la Inversa de una Matriz
void inverseMatrix(Matrix M, Matrix &Minv){
    Matrix Cof, Adj;
    float det = determinant(M);
    if(det == 0) 
        exit(EXIT_FAILURE);
    cofactors(M,Cof);
    transpose(Cof,Adj);
    productRealMatrix(1/det,Adj,Minv);
}
```

## Cambios en display_tools.h

```cpp
//Agregada una función para mostrar los vectores b
void showbs(vector<Vector> bs){
    for(int i=0;i<bs.size();i++){
        cout << "b del elemento "<< i+1 << ":\n";
        showVector(bs.at(i));
        cout << "*************************************\n";
    }
}

//Agregada una función para mostrar las matrices K
void showKs(vector<Matrix> Ks){
    for(int i=0;i<Ks.size();i++){
        cout << "K del elemento "<< i+1 << ":\n";
        showMatrix(Ks.at(i));
        cout << "*************************************\n";
    }
}
```

## Creando tools.h

```cpp
//Leer todos los datos del problem.msh
void leerMallayCondiciones(mesh &m){
    
    //Nombre del archivo, con longitud de 10 caracteres 
    char filename[10];
    
    ifstream file;
    float k,Q;
    int nnodes,neltos,ndirich,nneu;

   //Leer el nombre del archivo hasta que sea ingresado correctamente
  do{
     cout << "Ingrese el nombre del archivo que contiene los datos de la malla: ";
     cin >> filename;
     file.open(filename);
    }while(!file);

    //Leyendo los datos de Conductividad Térmica y de Calor
    file >> k >> Q;
    
    //Leyendo la cantidad de nodos, elementos, condiciones de dirichlet y Neumman
    file >> nnodes >> neltos >> ndirich >> nneu;

    //Almacenando los datos de Conductividad Térmica y de Calor en una malla
    m.setParameters(k,Q);
    
    //Estableciendo los tamaños de arreglos
    m.setSizes(nnodes,neltos,ndirich,nneu);
    
    //Inicializando la malla con los datos de los arreglos
    m.createData();

    //Leyendo los demas datos del problem.msh
    obtenerDatos(file,SINGLELINE,nnodes,INT_FLOAT,m.getNodes());
    obtenerDatos(file,DOUBLELINE,neltos,INT_INT_INT,m.getElements());
    obtenerDatos(file,DOUBLELINE,ndirich,INT_FLOAT,m.getDirichlet());
    obtenerDatos(file,DOUBLELINE,nneu,INT_FLOAT,m.getNeumann());

    file.close();
}
```

```cpp
//Leer datos del problem.msh segun el formato
void obtenerDatos(istream &file,int nlines,int n,int mode,item* item_list){
    string line;
    file >> line;
    
    //Omitir dos cadenas de caracteres
    if(nlines==DOUBLELINE) 
      file >> line;
   
   for(int i=0;i<n;i++){
        switch(mode){
        
        //Leer Coordenadas
        case INT_FLOAT:
            int e; float r;
            file >> e >> r;
            item_list[i].setIntFloat(e,r);
            break;
            
        //Leer elementos
        case INT_INT_INT:
            int e1,e2,e3;
            file >> e1 >> e2 >> e3;
            item_list[i].setIntIntInt(e1,e2,e3);
            break;
        }
    }
}
```

<hr>
<p align="center">Para servirles, <strong>Equipo de Instructores</strong> </p>

