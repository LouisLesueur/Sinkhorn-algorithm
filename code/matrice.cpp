#include <cmath>
#include "matrice.h"
#include <iostream>
#include <cassert>
#include <time.h>
#include <cstdlib>

using namespace std;

void Matrice::Affiche(){
	for(int i=0;i<nlig;i++){
		for(int j=0;j<ncol;j++){
			cout << get(i,j);
			cout << " ";
		}
		cout << " " << endl;
	}	
}

void Matrice::set(double x, int i, int j){ // Accesseur
	tab[i*ncol+j] = x;	
}	

double Matrice::get(int i, int j) const{
	return tab[i*ncol+j];
}

double Matrice::operator()(int i, int j){
	return get(i,j);
}


Matrice::Matrice(int n, int m){
	cout << "Construction" << endl;
	compteur = new int(1);
	nlig = n;
	ncol = m;
	tab = new double[n*m];
}


Matrice::Matrice(const Matrice & A){
	cout << "Construction" << endl;
	++*compteur;
	nlig = A.nlig;
	ncol = A.ncol;
	tab = A.tab;
}

Matrice::~Matrice(){
	cout << "Destruction" << endl;
	--*compteur;
	if(*compteur==0){
		delete[] tab;
		delete compteur;
	}
}

Matrice Matrice::operator*(Matrice & M){
         assert(ncol == M.nlig);
         Matrice C(nlig,M.ncol);
         for(int i=0;i<nlig;i++){ // On parcourt toutes les lignes
                 for(int j=0;j<M.ncol;j++){ // Puis toutes les colonnes
                	double tmp = 0;
			for(int k=0;k<ncol;k++){
				tmp += get(i,k)*M.get(k,j); // Formule de la multiplication matricielle
			}
			C.set(tmp,i,j);
		 }

	}
	return C;
}

void Matrice::Alea(){
	srand(time(NULL));
	for(int i=0;i<nlig*ncol;i++){
		tab[i] = rand()%21 -10; // Nombres aléatoires entre [-100,100]
	}
}
