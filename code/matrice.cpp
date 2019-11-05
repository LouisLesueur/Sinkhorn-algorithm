#pragma once
#include "matrice.h"

/*=======================================MATRICE================================================*/

/*------------------Opérateurs------------------*/

double Matrice::operator()(int i, int j) const {
    assert(0<=i && i<m && 0<=j && j<n);
    return tab[i*n+j];
}

double& Matrice::operator()(int i, int j) {
    assert(0<=i && i<m && 0<=j && j<n);
    return tab[i*n+j];
}

/*------------------Constructeurs------------------*/

Matrice::Matrice(int nligcol){
    compteur=new int(1);
    m=n=nligcol;
    tab=new double[n*m];
    cout<<"Nouvelle matrice carrée de taille "<<nligcol<<endl;
}

Matrice::Matrice(int nlig, int ncol){
    compteur=new int(1);
    m=nlig; n=ncol;
    tab=new double[m*n];
    cout<<"Nouvelle matrice de taille "<<nlig<<" "<<ncol<<endl;
}

Matrice::Matrice(const Matrice& A){
    compteur=A.compteur;
    *compteur+=1;
    m=A.m; n=A.n;
    tab=A.tab;
    cout<<"Matrice copiée ( shallow copy )"<<endl;
}

/*------------------Destructeur------------------*/

Matrice::~Matrice(){
    *compteur-=1;
    if(*compteur==0){
        delete compteur;
        delete [] tab;
        cout<<"Matrice détruite !"<<endl;
    }
}

/*=======================================================================================*/

Matrice operator*(const Matrice& A, const Matrice& B){
    assert(A.ncolonnes()==B.nlignes());
    int m=A.nlignes(); int n=B.ncolonnes();

    Matrice C(m,n);
    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            double S=0;
            for(int k=0; k<B.nlignes(); k++){
                S+=A(i,k)*B(k,j);
            }
            C(i,j)=S;
        }
    }
    return C;
}

ostream& operator<<( ostream& str, const Matrice& A) {
    for(int i=0; i<A.nlignes(); i++) {
        for(int j=0; j<A.ncolonnes(); j++) {
            str<<A(i,j)<<' ';
        }
        str<<endl;
    }
    return str;
}
