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
			monFlux<<" "<<double(tab[i+size*j]);
		monFlux<<endl;
	}
}


Matrice div(const Matrice a, const Matrice b){
	int n = a.nlignes();
	int m = a.ncolonnes();
	Matrice c(n,m);
	for(int i=0; i<n; i++){
		for(int j=0; j<m; j++){
			if(b(i,j) == 0)
				cout << "Division par 0 !" << endl;
			c(i,j) = a(i,j)/b(i,j);
		}
	}
	return c;
}

pi W(const simplex &s1, const simplex &s2, double eps, int n_iter){
	int m = s1.length();
	// Initialising
	Matrice ksi(m);
	for(int i=0; i<m; i++){
		for(int j=0; j<m; j++){
			double x_i = ((double)i+1/2)/m;
			double x_j = ((double)j+1/2)/m;
			double dis = pow(x_i - x_j, 2);
			ksi(i,j) = exp(-dis/eps);
		}
	}
	Matrice p1(m, 1);
	Matrice p2(m, 1);
	Matrice v(m, 1);
	for(int i=0; i<m; i++){
		p1(i,0) = s1(i);
		p2(i,0) = s2(i);
		v(i,0) = 1;
	}
	// Using recursion formula,
	for(int i=0; i<n_iter; i++){
		v = div(p2, transpose(ksi)*div(p1, ksi*v));
	}

	// Creating diag matrix
	Matrice u(m, 1);
	u = div(p1, ksi*v);
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

simplex bar(const simplex & s1, const simplex & s2, double lambda, double eps, int n_iter){
	int m = s1.length();
	// Initialising
	Matrice ksi(m);
	for(int i=0; i<m; i++){
		for(int j=0; j<m; j++){
			double x_i = ((double)i+1/2)/m;
			double x_j = ((double)j+1/2)/m;
			double dis = pow(x_i - x_j, 2);
			ksi(i,j) = exp(-dis/eps);
		}
	}
	Matrice p1(m, 1);
	Matrice p2(m, 1);
	Matrice v1(m, 1);
	Matrice v2(m, 1);
	Matrice u1(m, 1);
	Matrice u2(m, 1);
	for(int i=0; i<m; i++){ // ukn and vkn
		p1(i, 0) = s1(i);
		p2(i, 0) = s2(i);
		v1(i, 0) = 1;
		v2(i, 0) = 1;
		u1(i, 0) = 1;
		u2(i, 0) = 1;
	}
	for(int i=0; i<n_iter; i++){
		// Computing vk(n+1)
		v1 = div(p1, transpose(ksi)*u1);
		v2 = div(p1, transpose(ksi)*u2);
		// Computing p(n+1)
		Matrice p = product(matrix_pow(product(u1, ksi*v1), lambda), matrix_pow(product(u2, ksi*v2), 1-lambda));
		// Computing uk(n+1)
		u1 = div(p, ksi*v1);
		u2 = div(p, ksi*v2);
	}
	// Computing the first projection matrix
	Matrice diag1(m);
	Matrice diag2(m);
	for(int i=0; i<m; i++){
		for(int j=0; j<m; j++){
			if(i==j){
				diag1(i,j) = u1(i,0);
				diag2(i,j) = v1(i,0);
			}
			else{
				diag1(i,j) = 0;
				diag2(i,j) = 0;
			}
		}
	}
	Matrice final_prod(diag1*ksi*diag2);
	pi gamma(m);
	for(int i=0; i<m; i++)
		for(int j=0; j<m; j++)
			gamma(i,j) = final_prod(i,j);
	return gamma.first_marginal();
}
