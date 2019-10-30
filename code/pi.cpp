#include "pi.h"

//===================================CLASSE PI=========================

//-------------------------CONSTRUCTEUR--------------------------------
pi::pi(int n){
	size = n;
	tab = new double[n*n];
	double sum = 0;

	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			double coef = Random(0,100);
			tab[i+n*j] = coef;
			sum += coef;
		}
	}
	for(int i=0; i<n*n; i++)
		tab[i] /= sum;
}

//----------------------------OPERATEURS----------------------------------

double pi::operator()(int i, int j) const{
	assert(i<size && j<size);
	return tab[i+size*j];
}

double& pi::operator()(int i, int j){
	assert(i<size && j<size);
	return tab[i+size*j];
}

//------------------------MARGINALES---------------------------------------

simplex pi::first_marginal(){
	double marg[size];
	for(int i=0; i<size; i++){
		double sum = 0;
		for(int j=0; j<size; j++)
			sum += tab[i+size*j];
		marg[i] = sum;
	}
	return simplex(marg, size);
}

simplex pi::second_marginal(){
	double marg[size];
	for(int j=0; j<size; j++){
		double sum = 0;
		for(int i=0; i<size; i++)
			sum += tab[i+size*j];
		marg[j] = sum;
	}
	return simplex(marg, size);
}

//====================================================================
