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

Matrice Matrice::operator*(double k){
	for(int i=0; i<n*m; i++)
		tab[i]*=k;
}

Matrice Matrice::operator/(double k){
	for(int i=0; i<n*m; i++)
		tab[i]/=k;
}

Matrice &Matrice::operator=(const Matrice &source){
	compteur=source.compteur;
	*compteur+=1;
	m=source.m; n=source.n;
	tab = source.tab;
	return *this;

}

/*------------------Constructeurs------------------*/

Matrice::Matrice(int nligcol){
    compteur=new int(1);
    m=n=nligcol;
    tab=new double[m*n];
}

Matrice::Matrice(int nlig, int ncol){
    compteur=new int(1);
    m=nlig; n=ncol;
    tab=new double[m*n];
}

Matrice::Matrice(const Matrice& A){
    compteur=A.compteur;
    *compteur+=1;
    m=A.m; n=A.n;
    tab=A.tab;
}

/*------------------Destructeur------------------*/

Matrice::~Matrice(){
    *compteur-=1;
    if(*compteur==0){
        delete compteur;
        delete [] tab;
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

Matrice transpose(const Matrice& A){
	Matrice B(A.ncolonnes(), A.nlignes());
	for(int i=0; i<A.nlignes(); i++){
		for(int j=0; j<A.ncolonnes(); j++)
			B(j,i)=A(i,j);
	}
	return B;
}

Matrice matrix_pow(const Matrice& A, double lambda){
	Matrice C(A.nlignes(), A.ncolonnes());
	for(int i=0; i<A.nlignes(); i++){
		for(int j=0; j<A.ncolonnes(); j++){
			C(i,j) = pow(A(i,j), lambda);
		}
	}
	return C;
}

Matrice product(const Matrice& A, const Matrice& B){
	Matrice C(A.nlignes(), A.ncolonnes());
		for(int i=0; i<A.nlignes(); i++){
			for(int j=0; j<A.ncolonnes(); j++){
				C(i,j) = A(i,j)*B(i,j);
			}
		}
	return C;
}
Matrice operator/(const Matrice& A, const Matrice& B){
	Matrice C(A.nlignes(), A.ncolonnes());
	for(int i=0; i<A.nlignes(); i++){
		for(int j=0; j<A.ncolonnes(); j++){
			C(i,j)=A(i,j)/B(i,j);
		}
	}
	return C;
}
