#include "pi.h"
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

pi W(simplex s1, simplex s2, double eps, int n_iter){
	int m = s1.length();
	// Initialising
	vector<double> p1;
	vector<double> p2;
	vector<double> ksi;
	vector<double> v;
	vector<double> u;
	for(int i=0; i<m; i++){
		p1.push_back(s1(i));
		p2.push_back(s2(i));
		double dis = pow((s1(i) - s2(i)), 2);
		ksi.push_back(exp(-dis/eps));
		v.push_back(1.);
		u.push_back(p1[i]/(ksi[i]*v[i]));
	}
	for(int i=0; i<n_iter; i++){
		v = p2*(1/scalar(ksi, u));
		u = p1/(ksi*v);	
	}
}

