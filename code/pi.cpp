#include "pi.h"
#include "matrice.h"
#include <vector>
#include <cmath>
using namespace std;
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
	return simplex(marg, size, '1');
}

simplex pi::second_marginal(){
	double marg[size];
	for(int j=0; j<size; j++){
		double sum = 0;
		for(int i=0; i<size; i++)
			sum += tab[i+size*j];
		marg[j] = sum;
	}
	return simplex(marg, size, '2');
}

//====================================================================

void pi::plot(){
	ofstream monFlux("pi.csv");
	for(int i=0; i<size; i++){
		for(int j=0; j<size; j++)
			monFlux<<" "<<tab[i+size*j];
		monFlux<<endl;
	}
}


pi W(const simplex &s1, const simplex &s2, double eps, int n_iter){
	int m = s1.length();
	// Initialising
	Matrice ksi(m);
	for(int i=0; i<m; i++){
		for(int j=0; j<m; j++){
			double dis = pow((s1(i) - s2(j)), 2);
			ksi(i,j) = exp(-dis/eps);
		}
	}
	cout << "ksi=" << ksi << endl;
	Matrice p1(m, 1);
	Matrice p2(m, 1);
	Matrice v(m, 1);
	for(int i=0; i<m; i++){
		p1(i,0) = s1(i);
		p2(i,0) = s2(i);
		v(i,0) = 1;
	}
	Matrice u(m, 1);
	Matrice prod(ksi*v);
	for(int i=0; i<m; i++){
		u(i,0) = p1(i,0)/prod(i,0);
	}
	// Using recursion formula,
	for(int i=0; i<n_iter; i++){
		Matrice prod1(transpose(ksi)*u);
		for(int i=0; i<m; i++){
			v(i,0) = p2(i,0)/prod1(i,0);
		}
		Matrice prod2(ksi*v);
		for(int i=0; i<m; i++){
			u(i,0) = p1(i,0)/prod2(i,0);
		}
	}
	// Creating diag matrix
	Matrice diag1(m);
	Matrice diag2(m);
	for(int i=0; i<m; i++){
		for(int j=0; j<m; j++){
			if(i==j){
				diag1(i,j) = u(i,0);
				diag2(i,j) = v(i,0);
			}
			else{
				diag1(i,j) = 0;
				diag2(i,j) = 0;
			}
		}
	}
	Matrice final_prod(diag1*ksi*diag2);
	pi gamma(m);
	for(int i=0; i<m; i++){
		for(int j=0; j<m; j++){
			gamma(i,j) = final_prod(i,j);
		}
	}
	return gamma;
}
