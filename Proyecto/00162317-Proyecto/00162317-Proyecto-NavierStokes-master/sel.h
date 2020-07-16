void showMatrix(Matrix K){
    for(int i=0;i<K.at(0).size();i++){
        cout << "[\t";
        for(int j=0;j<K.size();j++){
            cout << K.at(i).at(j) << "\t";
        }
        cout << "]\n";
    }
}

void showKs(vector<Matrix> Ks){
    for(int i=0;i<Ks.size();i++){
        cout << "K del elemento "<< i+1 << ":\n";
        showMatrix(Ks.at(i));
        cout << "*************************************\n";
    }
}

void showVector(Vector b){
    cout << "[\t";
    for(int i=0;i<b.size();i++){
        cout << b.at(i) << "\t";
    }
    cout << "]\n";
}

void showbs(vector<Vector> bs){
    for(int i=0;i<bs.size();i++){
        cout << "b del elemento "<< i+1 << ":\n";
        showVector(bs.at(i));
        cout << "*************************************\n";
    }
}

//Mandamos un indice, elemento y la malla
node selectNode(int i, element e,mesh &m){
	node n;
	switch(i){
	    //Retorna el nodo de la malla que estemos obteniendo
		case 1: n = m.getNode(e.getNode1()-1); break;
		case 2: n = m.getNode(e.getNode2()-1); break;
		case 3: n = m.getNode(e.getNode3()-1); break;
        case 4: n = m.getNode(e.getNode4()-1); break;
	}
	return n;
}

float selectCoord(int c, node n){
	float v;
	switch(c){
	    //Aqui para obtener las coordenadas
		case EQUIS: v = n.getX(); break;
		case YE: v = n.getY(); break;
        case ZETA: v = n.getZ(); break;
	}
	return v;
}

//Mandamos el elemento, la coordenada, indice I, indice J y la malla
float calcularTenedor(element e, int coord, int i, int j,mesh &m){
    //Aqui calculamos la diferencia de por ejemplo: x2-x1,y2-y1
	node n1=selectNode(i,e,m),n2=selectNode(j,e,m);
	//Aqui retornamos esa diferencia
	return selectCoord(coord,n1) - selectCoord(coord,n2);
}

float operarRestaTenedor(float a, float b, float c, float d){
    return (a*b)-(c*d);
}

float calculateMagnitude(float v1, float v2){
    return sqrt(pow(v1,2)+pow(v2,2));
}

void ubicarSubMatriz(Matrix &K,int fi,int ff,int ci,int cf,Matrix M){
    int n = 0, m= 0;
    for(int i=fi;i<=ff;i++){
        for(int j=ci;j<=cf;j++){
            K.at(i).at(j) = M.at(n).at(m);
            m++;
        }
        n++; m = 0;
    }
}

float calculateLocalD(int i,mesh m){
    Matrix matriz;
    Vector row1, row2, row3;
    //Determinante
    element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);
    node n4 = m.getNode(e.getNode4()-1);

	row1.push_back(calcularTenedor(e,EQUIS,2,1,m));
    row1.push_back(calcularTenedor(e,YE,2,1,m));
    row1.push_back(calcularTenedor(e,ZETA,2,1,m));

	row2.push_back(calcularTenedor(e,EQUIS,3,1,m));
    row2.push_back(calcularTenedor(e,YE,3,1,m));
    row2.push_back(calcularTenedor(e,ZETA,3,1,m));

    row3.push_back(calcularTenedor(e,EQUIS,4,1,m));
    row3.push_back(calcularTenedor(e,YE,4,1,m));
    row3.push_back(calcularTenedor(e,ZETA,4,1,m));

	matriz.push_back(row1); matriz.push_back(row2); matriz.push_back(row3);

    return determinant(matriz);

}

float calculateLocalJ(int i,mesh m){
    Matrix matriz;
    //Filas
    Vector row1, row2, row3;
    //Jacobiano
    element e = m.getElement(i);

	row1.push_back(calcularTenedor(e,EQUIS,2,1,m));
    row1.push_back(calcularTenedor(e,EQUIS,3,1,m));
    row1.push_back(calcularTenedor(e,EQUIS,4,1,m));

    row2.push_back(calcularTenedor(e,YE,2,1,m));
    row2.push_back(calcularTenedor(e,YE,3,1,m));
    row2.push_back(calcularTenedor(e,YE,4,1,m));

    row3.push_back(calcularTenedor(e,ZETA,2,1,m));
    row3.push_back(calcularTenedor(e,ZETA,3,1,m));
    row3.push_back(calcularTenedor(e,ZETA,4,1,m));

	matriz.push_back(row1); matriz.push_back(row2); matriz.push_back(row3);

    return determinant(matriz);

}
    //ALPHA
void calculateLocalA(int i,Matrix &A,mesh m){
    //Video: 12:23
    //Matriz de 3x3
    /*Esta es una matriz ya obtenida en 2D, renombrada tambien como alfa
     * */
    zeroes(A,3);
    element e = m.getElement(i);

    /*|y3z4    y4z4    y4z3|
     *|x4z3    x2z4    x3z2|
     *|x3y4    x4y2    x2y3|
     * */
    //---------------------------------------------------------------------------PRIMERA FILA-------------------------------------------------------------------------------------------------
    //(y3-y1), (z4-z1), (y4-y1), (z3-z1) CALCULAR TENEDOR
    //((y31*z41)-(y41*z31)) RESTA TENEDOR
    A.at(0).at(0) = operarRestaTenedor(calcularTenedor(e,YE,3,1,m),calcularTenedor(e,ZETA,4,1,m),calcularTenedor(e,YE,4,1,m),calcularTenedor(e,ZETA,3,1,m));
    //(y4-y1), (z2-z1), (y2-y1), (z4-z1) CALCULAR TENEDOR
    //((y41*z21)-(y21*z41)) RESTA TENEDOR
    A.at(0).at(1) = operarRestaTenedor(calcularTenedor(e,YE,4,1,m),calcularTenedor(e,ZETA,2,1,m),calcularTenedor(e,YE,2,1,m),calcularTenedor(e,ZETA,4,1,m));
    //(y2-y1), (z3-z1), (y3-y1), (z31-21) CALCULAR TENEDOR
    //((y21*z31)-(y31*z21)) RESTA TENEDOR
    A.at(0).at(2) = operarRestaTenedor(calcularTenedor(e,YE,2,1,m),calcularTenedor(e,ZETA,3,1,m),calcularTenedor(e,YE,3,1,m),calcularTenedor(e,ZETA,2,1,m));
    //-------------------------------------------------------------------------SEGUNDA FILA---------------------------------------------------------------------------------------------------
    //((x41*z31)-(x31*z41))
    A.at(1).at(0) = operarRestaTenedor(calcularTenedor(e,EQUIS,4,1,m),calcularTenedor(e,ZETA,3,1,m),calcularTenedor(e,EQUIS,3,1,m),calcularTenedor(e,ZETA,4,1,m));
    //((y21*z41)-(x41*z21))
    A.at(1).at(1) = operarRestaTenedor(calcularTenedor(e,EQUIS,2,1,m),calcularTenedor(e,ZETA,4,1,m),calcularTenedor(e,EQUIS,4,1,m),calcularTenedor(e,ZETA,2,1,m));
    //((x31*z21)-(x21*z31))
    A.at(1).at(2) = operarRestaTenedor(calcularTenedor(e,EQUIS,3,1,m),calcularTenedor(e,ZETA,2,1,m),calcularTenedor(e,EQUIS,2,1,m),calcularTenedor(e,ZETA,3,1,m));
    //------------------------------------------------------------------------TERCERA FILA----------------------------------------------------------------------------------------------------
    //((x31*y41)-(x41*y31))
    A.at(2).at(0) = operarRestaTenedor(calcularTenedor(e,EQUIS,3,1,m),calcularTenedor(e,YE,4,1,m),calcularTenedor(e,EQUIS,4,1,m),calcularTenedor(e,YE,3,1,m));
    //((x41*y21)-(x21*y41))
    A.at(2).at(1) = operarRestaTenedor(calcularTenedor(e,EQUIS,4,1,m),calcularTenedor(e,YE,2,1,m),calcularTenedor(e,EQUIS,2,1,m),calcularTenedor(e,YE,4,1,m));
    //((x21*y31)-(x31*y21))
    A.at(2).at(2) = operarRestaTenedor(calcularTenedor(e,EQUIS,2,1,m),calcularTenedor(e,YE,3,1,m),calcularTenedor(e,EQUIS,3,1,m),calcularTenedor(e,YE,2,1,m));

}

//Matriz quemada
void calculateB(Matrix &B){
    zeroes(B,3,12);
    //Video: 11:35
    //Esta es una matriz de constantes
    /*Es un patron que se mantiene acarreado de 1D y 2D
     * Esta es Beta con dimensiones 3x12
     *
     *
     *      -1  1   0   0   -1  1   0   0   -1  1   0   0
     * */
    //---------------------------------------------------------------------------PRIMERA FILA-------------------------------------------------------------------------------------------------
    B.at(0).at(0) = -1; B.at(0).at(1) = 1; B.at(0).at(2) = 0; B.at(0).at(3) = 0;
    B.at(0).at(4) = -1; B.at(0).at(5) = 1; B.at(0).at(6) = 0; B.at(0).at(7) = 0;
    B.at(0).at(8) = -1; B.at(0).at(9) = 1; B.at(0).at(10) = 0; B.at(0).at(11) = 0;

    //---------------------------------------------------------------------------SEGUNDA FILA-------------------------------------------------------------------------------------------------
    B.at(1).at(0) = -1; B.at(1).at(1) = 0; B.at(1).at(2) = 1; B.at(1).at(3) = 0;
    B.at(1).at(4) = -1; B.at(1).at(5) = 0; B.at(1).at(6) = 1; B.at(1).at(7) = 0;
    B.at(1).at(8) = -1; B.at(1).at(9) = 0; B.at(1).at(10) = 1; B.at(1).at(11) = 0;

    //---------------------------------------------------------------------------TERCERA FILA-------------------------------------------------------------------------------------------------
    B.at(2).at(0) = -1; B.at(2).at(1) = 0; B.at(2).at(2) = 0; B.at(2).at(3) = 1;
    B.at(2).at(4) = -1; B.at(2).at(5) = 0; B.at(2).at(6) = 0; B.at(2).at(7) = 1;
    B.at(2).at(8) = -1; B.at(2).at(9) = 0; B.at(2).at(10) = 0; B.at(2).at(11) = 1;

}
//------------------------------------------------------------------COMPONENTES B,C,D,E------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------B

float componenteB1prueba(int i, mesh m){
    
    element e = m.getElement(i);

    node n1 = selectNode(1,e,m);

    float x21,x31,x41, x1;
    float y21,y31,y41, y1;
    float z21,z31,z41, z1;


    x1 = selectCoord(EQUIS,n1);

    y1 = selectCoord(YE,n1);
    
    z1 = selectCoord(ZETA,n1);
    
    //x2-x1
    x21 = calcularTenedor(e,EQUIS,2,1,m);
    //x3-x1
    x31 = calcularTenedor(e,EQUIS,3,1,m);
    //x4-x1
    x41 = calcularTenedor(e,EQUIS,4,1,m);
    
    //y2-y1
    y21 = calcularTenedor(e,YE,2,1,m);
    //y3-y1
    y31 = calcularTenedor(e,YE,3,1,m);
    //y4-y1
    y41 = calcularTenedor(e,YE,4,1,m);

    //z2-z1
    z21 = calcularTenedor(e,ZETA,2,1,m);
    //z3-z1
    z31 = calcularTenedor(e,ZETA,3,1,m);
    //z4-z1
    z41 = calcularTenedor(e,ZETA,4,1,m);
    //z21

    float expresion = ((x21*(3*z31+z41+6*z1)+x31*(3*z31+z41+6*z1)+2*x41*(z31+z41+3*z1)+6*x1*(2*z31+z41+5*z1)
        +y21*(3*z31+z41+6*z1)+y31*(3*z31+z41+6*z1)+2*(y41*(z31+z41+3*z1)+3*y1*(2*z31+z41+5*z1)))/720);

    return expresion;
}

//------------------------------------------------------------------------------------------------C
float componenteC1prueba(int i, mesh m ){
    
    element e = m.getElement(i);

    node n1 = selectNode(1,e,m);

    float x21,x31,x41, x1;
    float y21,y31,y41, y1;
    float z21,z31,z41, z1;


    x1 = selectCoord(EQUIS,n1);

    y1 = selectCoord(YE,n1);
    
    z1 = selectCoord(ZETA,n1);
    
    //x2-x1
    x21 = calcularTenedor(e,EQUIS,2,1,m);
    //x3-x1
    x31 = calcularTenedor(e,EQUIS,3,1,m);
    //x4-x1
    x41 = calcularTenedor(e,EQUIS,4,1,m);
    
    //y2-y1
    y21 = calcularTenedor(e,YE,2,1,m);
    //y3-y1
    y31 = calcularTenedor(e,YE,3,1,m);
    //y4-y1
    y41 = calcularTenedor(e,YE,4,1,m);

    //z2-z1
    z21 = calcularTenedor(e,ZETA,2,1,m);
    //z3-z1
    z31 = calcularTenedor(e,ZETA,3,1,m);
    //z4-z1
    z41 = calcularTenedor(e,ZETA,4,1,m);
    //z21

    float expresion = ((x21*(2*z21+z31+z41+6*z1)+x31*(z21+2*z31+z41+6*z1)+x41*(z21+z31+2*(z41+3*z1))
    +6*x1*(z21+z31+z41+5*z1)+y21*(2*z21+z31+z41+6*z1)+y31*(z21+2*z31+z41+6*z1)+y41*(z21+z31+2*(z41+3*z1))+6*z21*(z21+z31+z41+5*z1))/360);

    return expresion;
}

//------------------------------------------------------------------------------------------------D
float componenteD1prueba(int i, mesh m){
    
    element e = m.getElement(i);

    node n1 = selectNode(1,e,m);

    float x21,x31,x41, x1;
    float y21,y31,y41, y1;
    float z21,z31,z41, z1;

    x1 = selectCoord(EQUIS,n1);

    y1 = selectCoord(YE,n1);
    
    z1 = selectCoord(ZETA,n1);
    
    //x2-x1
    x21 = calcularTenedor(e,EQUIS,2,1,m);
    //x3-x1
    x31 = calcularTenedor(e,EQUIS,3,1,m);
    //x4-x1
    x41 = calcularTenedor(e,EQUIS,4,1,m);
    
    //y2-y1
    y21 = calcularTenedor(e,YE,2,1,m);
    //y3-y1
    y31 = calcularTenedor(e,YE,3,1,m);
    //y4-y1
    y41 = calcularTenedor(e,YE,4,1,m);

    //z2-z1
    z21 = calcularTenedor(e,ZETA,2,1,m);
    //z3-z1
    z31 = calcularTenedor(e,ZETA,3,1,m);
    //z4-z1
    z41 = calcularTenedor(e,ZETA,4,1,m);
    //z21

    float expresion = ((y21*(2*z21+z31+z41+6*z1)+y31*(z21+2*z31+z41+6*z1)+y41*(z21+z31+2*(z41+3*z1))
    +6*y1*(z21+z31+z41+5*z1))/720);

    return expresion;
}
//------------------------------------------------------------------------------------------------E
float componenteE1prueba(int i, mesh m){
    
    element e = m.getElement(i);

    node n1 = selectNode(1,e,m);

    float x21,x31,x41,x1;
    float y21,y31,y41,y1;
    float z21,z31,z41,z1;

    x1 = selectCoord(EQUIS,n1);

    y1 = selectCoord(YE,n1);
    
    z1 = selectCoord(ZETA,n1);
    
    //x2-x1
    x21 = calcularTenedor(e,EQUIS,2,1,m);
    //x3-x1
    x31 = calcularTenedor(e,EQUIS,3,1,m);
    //x4-x1
    x41 = calcularTenedor(e,EQUIS,4,1,m);
    
    //y2-y1
    y21 = calcularTenedor(e,YE,2,1,m);
    //y3-y1
    y31 = calcularTenedor(e,YE,3,1,m);
    //y4-y1
    y41 = calcularTenedor(e,YE,4,1,m);

    //z2-z1
    z21 = calcularTenedor(e,ZETA,2,1,m);
    //z3-z1
    z31 = calcularTenedor(e,ZETA,3,1,m);
    //z4-z1
    z41 = calcularTenedor(e,ZETA,4,1,m);
    //z21

    float expresion = ((x21*(2*z21+z31+z41+6*z1)+x31*(z21+2*z31+z41+6*z1)+x41*(z21+z31+2*z41+3*z1)
    +6*x1*(z21+z31+z41+5*z1))/240);

    return expresion;
}


void calculateGammaMatrix(Matrix& m){
    zeroes(m,12,3);
    //Video: 13:48
    /*Matriz 12x3
     * Respuesta a la integral triple de N transpuesta
     * dando como resultado un vector 4x1, que tiene 4 veces 1/24
     * Gamma
     * */
	float e = 1/24.0;
	m.at(0).at(0) = e; /**/m.at(0).at(1) = 0; /**/m.at(0).at(2) = 0;
	m.at(1).at(0) = e; /**/m.at(1).at(1) = 0; /**/m.at(1).at(2) = 0;
   	m.at(2).at(0) = e; /**/m.at(2).at(1) = 0; /**/m.at(2).at(2) = 0;
   	m.at(3).at(0) = e; /**/m.at(3).at(1) = 0; /**/m.at(3).at(2) = 0;
   	m.at(4).at(0) = 0; /**/m.at(4).at(1) = e; /**/m.at(4).at(2) = 0;
   	m.at(5).at(0) = 0; /**/m.at(5).at(1) = e; /**/m.at(5).at(2) = 0;
   	m.at(6).at(0) = 0; /**/m.at(6).at(1) = e; /**/m.at(6).at(2) = 0;
   	m.at(7).at(0) = 0; /**/m.at(7).at(1) = e; /**/m.at(7).at(2) = 0;
   	m.at(8).at(0) = 0; /**/m.at(8).at(1) = 0; /**/m.at(8).at(2) = e;
   	m.at(9).at(0) = 0; /**/m.at(9).at(1) = 0; /**/m.at(9).at(2) = e;
   	m.at(10).at(0) = 0; /**/m.at(10).at(1) = 0; /**/m.at(10).at(2) = e;
   	m.at(11).at(0) = 0; /**/m.at(11).at(1) = 0; /**/m.at(11).at(2) = e;

}

//Recibimos el elemento y la malla
Matrix createLocalK(int e,mesh &m){

    //Respuesta de mi forma debil
    //[k+m  L]
    //[T    0]

    //Preparacion de ingredientes
    //Jacobiano, determinante son lo mismo en forma de codigo
    float u_bar,nu,rho,Ve,J,D,P;
    Matrix alpha,beta,gamma,delta;
    Matrix K,Alpha,Beta,Alphat,Betat,BPrima;

    //Preparando alpha
    u_bar = m.getParameter(ADJECTIVE_VELOCITY);
    J = calculateLocalJ(e,m);//Jacobiano
    D = calculateLocalD(e,m);//Determinante

    if(D == 0){
        cout << "\n!---DETERMINANTE 0---!\n";
        exit(EXIT_FAILURE);
    }

    calculateLocalA(e,Alpha,m);//alpha
    calculateB(Beta);//Matriz quemada
    
   //-----------------------------------------------------------------------------------------------------COMPONENTE B
    //[k+m  L]
    //[T    0]
    //Matriz K
    float mP = componenteB1prueba(e,m);//Nuevo
    float realB = mP/(pow(D,2));//Nuevo
    transpose(Alpha,Alphat);
    transpose(Beta,Betat);
    //B = P/D**2 * Bt * alphaT * B * alpha
    productRealMatrix(realB, productMatrixMatrix(Betat, productMatrixMatrix(Alphat, productMatrixMatrix(Alpha, Beta, 3, 3, 12), 3, 3, 12), 12, 3, 12),delta);//Guardamos en beta

    //-----------------------------------------------------------------------------------------------------COMPONENTE C
    //[k+m  L]
    //[T    0]
    //Matriz M
    float Bprima = componenteC1prueba(e,m);
    float realC = Bprima/pow(D,2);
    transpose(Alpha,Alphat);
    transpose(Beta,Betat);
    //C = Bprima/D**2 * Bt * alphaT * B * alpha
    productRealMatrix(realC, productMatrixMatrix(Betat, productMatrixMatrix(Alphat, productMatrixMatrix(Alpha, Beta, 3, 3, 12), 3, 3, 12), 12, 3, 12),beta);

    //-----------------------------------------------------------------------------------------------------COMPONENTE D
    //[k+m  L]
    //[T    0]
    //Matriz L
    float alphaPrima = componenteD1prueba(e,m);
    float realD = alphaPrima/pow(D,2);
    transpose(Alpha,Alphat);
    transpose(Beta,Betat);
    //D = alphaPrima/D**2 * Bt * alphaT * B * alpha
    productRealMatrix(realD, productMatrixMatrix(Betat, productMatrixMatrix(Alphat, productMatrixMatrix(Alpha, Beta, 3, 3, 12), 3, 3, 12), 12, 3, 12),alpha);



    //-----------------------------------------------------------------------------------------------------COMPONENTE E
    //[k+m  L]
    //[T    0]
    //Matriz T
    float gama = componenteE1prueba(e,m);
    float realE = gama/pow(D,2);
    transpose(Alpha,Alphat);
    transpose(Beta,Betat);
    //E = gama/D**2 * Bt * alphaT * B * alpha
    productRealMatrix(realE, productMatrixMatrix(Betat, productMatrixMatrix(Alphat, productMatrixMatrix(Alpha, Beta, 3, 3, 12), 3, 3, 12), 12, 3, 12),gamma);

    

    //[k+m  L]
    //[T    0]
    //Colocando submatrices en K
    zeroes(K,16);
    ubicarSubMatriz(K,0,11,0,11,sumMatrix(alpha,beta,12,12));
    ubicarSubMatriz(K,0,11,12,15,gamma);
    ubicarSubMatriz(K,12,15,0,11,delta);
    return K;
}

Vector createLocalb(int e,mesh &m){
    Vector b0,b,f;
    Matrix e_star;

    float f_x = m.getParameter(EXTERNAL_FORCE_X);
    float f_y = m.getParameter(EXTERNAL_FORCE_Y);
    float f_z = m.getParameter(EXTERNAL_FORCE_Z);
    float J = calculateLocalJ(e,m);
    calculateGammaMatrix(e_star);
    zeroes(f,3);
    f.at(0) = f_x;
    f.at(1) = f_y;
    f.at(2) = f_z;

    zeroes(b0,12);
    productMatrixVector(e_star,f,b0);
    productRealVector(J,b0,b);
    b.push_back(0); b.push_back(0); b.push_back(0);

    return b;
}

void crearSistemasLocales(mesh &m,vector<Matrix> &localKs,vector<Vector> &localbs){
    for(int i=0;i<m.getSize(ELEMENTS);i++){
        localKs.push_back(createLocalK(i,m));
        localbs.push_back(createLocalb(i,m));
    }
}

void assemblyK(element e,Matrix localK,Matrix &K,int nnodes){
    int index1  = e.getNode1() - 1;
    int index2  = e.getNode2() - 1;
    int index3  = e.getNode3() - 1;
    int index4  = e.getNode4() - 1;

    int index5  = index1+nnodes;
    int index6  = index2+nnodes;
    int index7  = index3+nnodes;
    int index8  = index4+nnodes;

    int index9  = index1+2*nnodes;
    int index10 = index2+2*nnodes;
    int index11 = index3+2*nnodes;
    int index12 = index4+2*nnodes;
    int index13 = index1+3*nnodes;
    int index14 = index2+3*nnodes;
    int index15 = index3+3*nnodes;
    int index16 = index4+3*nnodes;


    K.at(index1).at(index1)  += localK.at(0).at(0);
    K.at(index1).at(index2)  += localK.at(0).at(1);
    K.at(index1).at(index3)  += localK.at(0).at(2);
    K.at(index1).at(index4)  += localK.at(0).at(3);
    K.at(index1).at(index5)  += localK.at(0).at(4);
    K.at(index1).at(index6)  += localK.at(0).at(5);
    K.at(index1).at(index7)  += localK.at(0).at(6);
    K.at(index1).at(index8)  += localK.at(0).at(7);
    K.at(index1).at(index9)  += localK.at(0).at(8);
    K.at(index1).at(index10) += localK.at(0).at(9);
    K.at(index1).at(index11) += localK.at(0).at(10);
    K.at(index1).at(index12) += localK.at(0).at(11);
    K.at(index1).at(index13) += localK.at(0).at(12);
    K.at(index1).at(index14) += localK.at(0).at(13);
    K.at(index1).at(index15) += localK.at(0).at(14);
    K.at(index1).at(index16) += localK.at(0).at(15);

    K.at(index2).at(index1)  += localK.at(1).at(0);
    K.at(index2).at(index2)  += localK.at(1).at(1);
    K.at(index2).at(index3)  += localK.at(1).at(2);
    K.at(index2).at(index4)  += localK.at(1).at(3);
    K.at(index2).at(index5)  += localK.at(1).at(4);
    K.at(index2).at(index6)  += localK.at(1).at(5);
    K.at(index2).at(index7)  += localK.at(1).at(6);
    K.at(index2).at(index8)  += localK.at(1).at(7);
    K.at(index2).at(index9)  += localK.at(1).at(8);
    K.at(index2).at(index10) += localK.at(1).at(9);
    K.at(index2).at(index11) += localK.at(1).at(10);
    K.at(index2).at(index12) += localK.at(1).at(11);
    K.at(index2).at(index13) += localK.at(1).at(12);
    K.at(index2).at(index14) += localK.at(1).at(13);
    K.at(index2).at(index15) += localK.at(1).at(14);
    K.at(index2).at(index16) += localK.at(1).at(15);

    K.at(index3).at(index1)  += localK.at(2).at(0);
    K.at(index3).at(index2)  += localK.at(2).at(1);
    K.at(index3).at(index3)  += localK.at(2).at(2);
    K.at(index3).at(index4)  += localK.at(2).at(3);
    K.at(index3).at(index5)  += localK.at(2).at(4);
    K.at(index3).at(index6)  += localK.at(2).at(5);
    K.at(index3).at(index7)  += localK.at(2).at(6);
    K.at(index3).at(index8)  += localK.at(2).at(7);
    K.at(index3).at(index9)  += localK.at(2).at(8);
    K.at(index3).at(index10) += localK.at(2).at(9);
    K.at(index3).at(index11) += localK.at(2).at(10);
    K.at(index3).at(index12) += localK.at(2).at(11);
    K.at(index3).at(index13) += localK.at(2).at(12);
    K.at(index3).at(index14) += localK.at(2).at(13);
    K.at(index3).at(index15) += localK.at(2).at(14);
    K.at(index3).at(index16) += localK.at(2).at(15);

    K.at(index4).at(index1)  += localK.at(3).at(0);
    K.at(index4).at(index2)  += localK.at(3).at(1);
    K.at(index4).at(index3)  += localK.at(3).at(2);
    K.at(index4).at(index4)  += localK.at(3).at(3);
    K.at(index4).at(index5)  += localK.at(3).at(4);
    K.at(index4).at(index6)  += localK.at(3).at(5);
    K.at(index4).at(index7)  += localK.at(3).at(6);
    K.at(index4).at(index8)  += localK.at(3).at(7);
    K.at(index4).at(index9)  += localK.at(3).at(8);
    K.at(index4).at(index10) += localK.at(3).at(9);
    K.at(index4).at(index11) += localK.at(3).at(10);
    K.at(index4).at(index12) += localK.at(3).at(11);
    K.at(index4).at(index13) += localK.at(3).at(12);
    K.at(index4).at(index14) += localK.at(3).at(13);
    K.at(index4).at(index15) += localK.at(3).at(14);
    K.at(index4).at(index16) += localK.at(3).at(15);

    K.at(index5).at(index1)  += localK.at(4).at(0);
    K.at(index5).at(index2)  += localK.at(4).at(1);
    K.at(index5).at(index3)  += localK.at(4).at(2);
    K.at(index5).at(index4)  += localK.at(4).at(3);
    K.at(index5).at(index5)  += localK.at(4).at(4);
    K.at(index5).at(index6)  += localK.at(4).at(5);
    K.at(index5).at(index7)  += localK.at(4).at(6);
    K.at(index5).at(index8)  += localK.at(4).at(7);
    K.at(index5).at(index9)  += localK.at(4).at(8);
    K.at(index5).at(index10) += localK.at(4).at(9);
    K.at(index5).at(index11) += localK.at(4).at(10);
    K.at(index5).at(index12) += localK.at(4).at(11);
    K.at(index5).at(index13) += localK.at(4).at(12);
    K.at(index5).at(index14) += localK.at(4).at(13);
    K.at(index5).at(index15) += localK.at(4).at(14);
    K.at(index5).at(index16) += localK.at(4).at(15);

    K.at(index6).at(index1)  += localK.at(5).at(0);
    K.at(index6).at(index2)  += localK.at(5).at(1);
    K.at(index6).at(index3)  += localK.at(5).at(2);
    K.at(index6).at(index4)  += localK.at(5).at(3);
    K.at(index6).at(index5)  += localK.at(5).at(4);
    K.at(index6).at(index6)  += localK.at(5).at(5);
    K.at(index6).at(index7)  += localK.at(5).at(6);
    K.at(index6).at(index8)  += localK.at(5).at(7);
    K.at(index6).at(index9)  += localK.at(5).at(8);
    K.at(index6).at(index10) += localK.at(5).at(9);
    K.at(index6).at(index11) += localK.at(5).at(10);
    K.at(index6).at(index12) += localK.at(5).at(11);
    K.at(index6).at(index13) += localK.at(5).at(12);
    K.at(index6).at(index14) += localK.at(5).at(13);
    K.at(index6).at(index15) += localK.at(5).at(14);
    K.at(index6).at(index16) += localK.at(5).at(15);

    K.at(index7).at(index1)  += localK.at(6).at(0);
    K.at(index7).at(index2)  += localK.at(6).at(1);
    K.at(index7).at(index3)  += localK.at(6).at(2);
    K.at(index7).at(index4)  += localK.at(6).at(3);
    K.at(index7).at(index5)  += localK.at(6).at(4);
    K.at(index7).at(index6)  += localK.at(6).at(5);
    K.at(index7).at(index7)  += localK.at(6).at(6);
    K.at(index7).at(index8)  += localK.at(6).at(7);
    K.at(index7).at(index9)  += localK.at(6).at(8);
    K.at(index7).at(index10) += localK.at(6).at(9);
    K.at(index7).at(index11) += localK.at(6).at(10);
    K.at(index7).at(index12) += localK.at(6).at(11);
    K.at(index7).at(index13) += localK.at(6).at(12);
    K.at(index7).at(index14) += localK.at(6).at(13);
    K.at(index7).at(index15) += localK.at(6).at(14);
    K.at(index7).at(index16) += localK.at(6).at(15);

    K.at(index8).at(index1)  += localK.at(7).at(0);
    K.at(index8).at(index2)  += localK.at(7).at(1);
    K.at(index8).at(index3)  += localK.at(7).at(2);
    K.at(index8).at(index4)  += localK.at(7).at(3);
    K.at(index8).at(index5)  += localK.at(7).at(4);
    K.at(index8).at(index6)  += localK.at(7).at(5);
    K.at(index8).at(index7)  += localK.at(7).at(6);
    K.at(index8).at(index8)  += localK.at(7).at(7);
    K.at(index8).at(index9)  += localK.at(7).at(8);
    K.at(index8).at(index10) += localK.at(7).at(9);
    K.at(index8).at(index11) += localK.at(7).at(10);
    K.at(index8).at(index12) += localK.at(7).at(11);
    K.at(index8).at(index13) += localK.at(7).at(12);
    K.at(index8).at(index14) += localK.at(7).at(13);
    K.at(index8).at(index15) += localK.at(7).at(14);
    K.at(index8).at(index16) += localK.at(7).at(15);

    K.at(index9).at(index1)  += localK.at(8).at(0);
    K.at(index9).at(index2)  += localK.at(8).at(1);
    K.at(index9).at(index3)  += localK.at(8).at(2);
    K.at(index9).at(index4)  += localK.at(8).at(3);
    K.at(index9).at(index5)  += localK.at(8).at(4);
    K.at(index9).at(index6)  += localK.at(8).at(5);
    K.at(index9).at(index7)  += localK.at(8).at(6);
    K.at(index9).at(index8)  += localK.at(8).at(7);
    K.at(index9).at(index9)  += localK.at(8).at(8);
    K.at(index9).at(index10) += localK.at(8).at(9);
    K.at(index9).at(index11) += localK.at(8).at(10);
    K.at(index9).at(index12) += localK.at(8).at(11);
    K.at(index9).at(index13) += localK.at(8).at(12);
    K.at(index9).at(index14) += localK.at(8).at(13);
    K.at(index9).at(index15) += localK.at(8).at(14);
    K.at(index9).at(index16) += localK.at(8).at(15);

    K.at(index10).at(index1)  += localK.at(9).at(0);
    K.at(index10).at(index2)  += localK.at(9).at(1);
    K.at(index10).at(index3)  += localK.at(9).at(2);
    K.at(index10).at(index4)  += localK.at(9).at(3);
    K.at(index10).at(index5)  += localK.at(9).at(4);
    K.at(index10).at(index6)  += localK.at(9).at(5);
    K.at(index10).at(index7)  += localK.at(9).at(6);
    K.at(index10).at(index8)  += localK.at(9).at(7);
    K.at(index10).at(index9)  += localK.at(9).at(8);
    K.at(index10).at(index10) += localK.at(9).at(9);
    K.at(index10).at(index11) += localK.at(9).at(10);
    K.at(index10).at(index12) += localK.at(9).at(11);
    K.at(index10).at(index13) += localK.at(9).at(12);
    K.at(index10).at(index14) += localK.at(9).at(13);
    K.at(index10).at(index15) += localK.at(9).at(14);
    K.at(index10).at(index16) += localK.at(9).at(15);

    K.at(index11).at(index1)  += localK.at(10).at(0);
    K.at(index11).at(index2)  += localK.at(10).at(1);
    K.at(index11).at(index3)  += localK.at(10).at(2);
    K.at(index11).at(index4)  += localK.at(10).at(3);
    K.at(index11).at(index5)  += localK.at(10).at(4);
    K.at(index11).at(index6)  += localK.at(10).at(5);
    K.at(index11).at(index7)  += localK.at(10).at(6);
    K.at(index11).at(index8)  += localK.at(10).at(7);
    K.at(index11).at(index9)  += localK.at(10).at(8);
    K.at(index11).at(index10) += localK.at(10).at(9);
    K.at(index11).at(index11) += localK.at(10).at(10);
    K.at(index11).at(index12) += localK.at(10).at(11);
    K.at(index11).at(index13) += localK.at(10).at(12);
    K.at(index11).at(index14) += localK.at(10).at(13);
    K.at(index11).at(index15) += localK.at(10).at(14);
    K.at(index11).at(index16) += localK.at(10).at(15);

    K.at(index12).at(index1)  += localK.at(11).at(0);
    K.at(index12).at(index2)  += localK.at(11).at(1);
    K.at(index12).at(index3)  += localK.at(11).at(2);
    K.at(index12).at(index4)  += localK.at(11).at(3);
    K.at(index12).at(index5)  += localK.at(11).at(4);
    K.at(index12).at(index6)  += localK.at(11).at(5);
    K.at(index12).at(index7)  += localK.at(11).at(6);
    K.at(index12).at(index8)  += localK.at(11).at(7);
    K.at(index12).at(index9)  += localK.at(11).at(8);
    K.at(index12).at(index10) += localK.at(11).at(9);
    K.at(index12).at(index11) += localK.at(11).at(10);
    K.at(index12).at(index12) += localK.at(11).at(11);
    K.at(index12).at(index13) += localK.at(11).at(12);
    K.at(index12).at(index14) += localK.at(11).at(13);
    K.at(index12).at(index15) += localK.at(11).at(14);
    K.at(index12).at(index16) += localK.at(11).at(15);

    K.at(index13).at(index1)  += localK.at(12).at(0);
    K.at(index13).at(index2)  += localK.at(12).at(1);
    K.at(index13).at(index3)  += localK.at(12).at(2);
    K.at(index13).at(index4)  += localK.at(12).at(3);
    K.at(index13).at(index5)  += localK.at(12).at(4);
    K.at(index13).at(index6)  += localK.at(12).at(5);
    K.at(index13).at(index7)  += localK.at(12).at(6);
    K.at(index13).at(index8)  += localK.at(12).at(7);
    K.at(index13).at(index9)  += localK.at(12).at(8);
    K.at(index13).at(index10) += localK.at(12).at(9);
    K.at(index13).at(index11) += localK.at(12).at(10);
    K.at(index13).at(index12) += localK.at(12).at(11);
    K.at(index13).at(index13) += localK.at(12).at(12);
    K.at(index13).at(index14) += localK.at(12).at(13);
    K.at(index13).at(index15) += localK.at(12).at(14);
    K.at(index13).at(index16) += localK.at(12).at(15);

    K.at(index14).at(index1)  += localK.at(13).at(0);
    K.at(index14).at(index2)  += localK.at(13).at(1);
    K.at(index14).at(index3)  += localK.at(13).at(2);
    K.at(index14).at(index4)  += localK.at(13).at(3);
    K.at(index14).at(index5)  += localK.at(13).at(4);
    K.at(index14).at(index6)  += localK.at(13).at(5);
    K.at(index14).at(index7)  += localK.at(13).at(6);
    K.at(index14).at(index8)  += localK.at(13).at(7);
    K.at(index14).at(index9)  += localK.at(13).at(8);
    K.at(index14).at(index10) += localK.at(13).at(9);
    K.at(index14).at(index11) += localK.at(13).at(10);
    K.at(index14).at(index12) += localK.at(13).at(11);
    K.at(index14).at(index13) += localK.at(13).at(12);
    K.at(index14).at(index14) += localK.at(13).at(13);
    K.at(index14).at(index15) += localK.at(13).at(14);
    K.at(index14).at(index16) += localK.at(13).at(15);

    K.at(index15).at(index1)  += localK.at(14).at(0);
    K.at(index15).at(index2)  += localK.at(14).at(1);
    K.at(index15).at(index3)  += localK.at(14).at(2);
    K.at(index15).at(index4)  += localK.at(14).at(3);
    K.at(index15).at(index5)  += localK.at(14).at(4);
    K.at(index15).at(index6)  += localK.at(14).at(5);
    K.at(index15).at(index7)  += localK.at(14).at(6);
    K.at(index15).at(index8)  += localK.at(14).at(7);
    K.at(index15).at(index9)  += localK.at(14).at(8);
    K.at(index15).at(index10) += localK.at(14).at(9);
    K.at(index15).at(index11) += localK.at(14).at(10);
    K.at(index15).at(index12) += localK.at(14).at(11);
    K.at(index15).at(index13) += localK.at(14).at(12);
    K.at(index15).at(index14) += localK.at(14).at(13);
    K.at(index15).at(index15) += localK.at(14).at(14);
    K.at(index15).at(index16) += localK.at(14).at(15);

    K.at(index16).at(index1)  += localK.at(15).at(0);
    K.at(index16).at(index2)  += localK.at(15).at(1);
    K.at(index16).at(index3)  += localK.at(15).at(2);
    K.at(index16).at(index4)  += localK.at(15).at(3);
    K.at(index16).at(index5)  += localK.at(15).at(4);
    K.at(index16).at(index6)  += localK.at(15).at(5);
    K.at(index16).at(index7)  += localK.at(15).at(6);
    K.at(index16).at(index8)  += localK.at(15).at(7);
    K.at(index16).at(index9)  += localK.at(15).at(8);
    K.at(index16).at(index10) += localK.at(15).at(9);
    K.at(index16).at(index11) += localK.at(15).at(10);
    K.at(index16).at(index12) += localK.at(15).at(11);
    K.at(index16).at(index13) += localK.at(15).at(12);
    K.at(index16).at(index14) += localK.at(15).at(13);
    K.at(index16).at(index15) += localK.at(15).at(14);
    K.at(index16).at(index16) += localK.at(15).at(15);

}

void assemblyb(element e,Vector localb,Vector &b,int nnodes){
    int index1  = e.getNode1() - 1;
    int index2  = e.getNode2() - 1;
    int index3  = e.getNode3() - 1;
    int index4  = e.getNode4() - 1;

    int index5  = index1+nnodes;
    int index6  = index2+nnodes;
    int index7  = index3+nnodes;
    int index8  = index4+nnodes;

    int index9  = index1+2*nnodes;
    int index10 = index2+2*nnodes;
    int index11 = index3+2*nnodes;
    int index12 = index4+2*nnodes;


    b.at(index1) += localb.at(0);
    b.at(index2) += localb.at(1);
    b.at(index3) += localb.at(2);
    b.at(index4) += localb.at(3);
    b.at(index5) += localb.at(4);
    b.at(index6) += localb.at(5);
    b.at(index7) += localb.at(6);
    b.at(index8) += localb.at(7);
    b.at(index9) += localb.at(8);
    b.at(index10) += localb.at(9);
    b.at(index11) += localb.at(10);
    b.at(index12) += localb.at(11);

}

void ensamblaje(mesh &m,vector<Matrix> &localKs,vector<Vector> &localbs,Matrix &K,Vector &b){
    int nnodes = m.getSize(NODES);
    for(int i=0;i<m.getSize(ELEMENTS);i++){
        element e = m.getElement(i);
        assemblyK(e,localKs.at(i),K,nnodes);
        assemblyb(e,localbs.at(i),b,nnodes);
    }
}

void applyDirichlet(mesh &m,Matrix &K,Vector &b){
    for(int i=0;i<m.getSize(DIRICHLET);i++){
        condition c = m.getCondition(i,DIRICHLET);
        int index = c.getNode1()-1;

        K.erase(K.begin()+index);
        b.erase(b.begin()+index);

        for(int row=0;row<K.size();row++){
            float cell = K.at(row).at(index);
            K.at(row).erase(K.at(row).begin()+index);
            b.at(row) += -1*c.getValue()*cell;
        }
    }
}

void calculate(Matrix &K, Vector &b, vector<Vector> &Ts,int size){
    Vector v1,v2,T,dT,Tnext;

    zeroes(v1,size);
    zeroes(v2,size);
    zeroes(T,size);
    zeroes(dT,size);

    float tf = 3;
    float delta_t = 0.5;
    float t = 0.5;
    Ts.push_back(T);

    do{
        //Tnext=T+dt*(b-K*T) //Forward Euler

        productMatrixVector(K,T,v1);
        productRealVector(-1,v1,v2);
        productRealVector(delta_t,sumVector(b,v2,size),dT);

        Tnext.clear();
        copyVector(sumVector(T,dT,size),Tnext);
        Ts.push_back(Tnext);

        T.clear();
        copyVector(Tnext,T);

        t+=delta_t;

    }while(t<=tf);
}

