Matrix createLocalK(int element,mesh &m){
    Matrix K;
    Vector row1, row2;

    float var = m.getParameter(THERMAL_CONDUCTIVITY),
            l = (m.getNode(element+1).getX() - m.getNode(element).getX());

    row1.push_back(0.5*(var)/l); row1.push_back(-0.5*(var)/l);
    row2.push_back(-0.5*(var)/l); row2.push_back(0.5*(var)/l);

    K.push_back(row1); K.push_back(row2);

    return K;
}

Vector createLocalb(int element,mesh &m){
    Vector b;
    float var = m.getParameter(HEAT_SOURCE),
            l = m.getNode(element+1).getX() - m.getNode(element).getX();

    b.push_back(((-1)*var*l)/2);
    b.push_back(((-1)*var*l)/2);

    return b;
}

void crearSistemasLocales(mesh &m,vector<Matrix> &localKs,vector<Vector> &localbs){
    for(int i=0;i<m.getSize(ELEMENTS);i++){
        localKs.push_back(createLocalK(i,m));
        localbs.push_back(createLocalb(i,m));
    }
}

void assemblyK(element e,Matrix localK,Matrix &K){
    int index1 = e.getNode1() - 1; //0 - 1
    int index2 = e.getNode2() - 1; //1 - 2

    K.at(index1).at(index1) += localK.at(0).at(0);
    K.at(index1).at(index2) += localK.at(0).at(1);
    K.at(index2).at(index1) += localK.at(1).at(0);
    K.at(index2).at(index2) += localK.at(1).at(1);
}

void assemblyb(element e,Vector localb,Vector &b){
    int index1 = e.getNode1() - 1;
    int index2 = e.getNode2() - 1;

    b.at(index1) += localb.at(0);
    b.at(index2) += localb.at(1);
}

void ensamblaje(mesh &m,vector<Matrix> &localKs,vector<Vector> &localbs,Matrix &K,Vector &b){
    for(int i=0;i<m.getSize(ELEMENTS);i++){
        element e = m.getElement(i);
        assemblyK(e,localKs.at(i),K);
        assemblyb(e,localbs.at(i),b);
    }
}

void applyNeumann(mesh &m,Vector &b){
    for(int i=0;i<m.getSize(NEUMANN);i++){
        condition c = m.getCondition(i,NEUMANN);
        b.at(c.getNode1()-1) += c.getValue();
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

void calculate(Matrix &K, Vector &b, Vector &T){
    Matrix Kinv;
    inverseMatrix(K,Kinv);
    productMatrixVector(Kinv,b,T);
}