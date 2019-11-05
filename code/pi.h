#pragma once
#include "simplex.h"
#include<vector>
using namespace std;
class pi{
	private:
		int size;
		double* tab;
	public:
		//Construire une matrice pi de taille n
		pi(int n);
		//Avoir la taille de la matrice
		int length() const{return size;};

		//Pour acceder aux éléments
		double operator()(int i, int j) const;
		double& operator()(int i, int j);

		//Pour accèder aux marginales
		simplex first_marginal();
		simplex second_marginal();

		~pi(){delete [] tab;}
};

vector<double> operator*(vector<double> p1, vector<double> p2){
	vector<double> q;
	for(int i=0; i<p1.size(); i++){
		q.push_back(p1[i]*p2[i]);
	}
	return q;
}


vector<double> operator/(vector<double> p1, vector<double> p2){
	vector<double> q;
	for(int i=0; i<p1.size(); i++){
		q.push_back(p1[i]/p2[i]);
	}
	return q;
}


vector<double> operator*(vector<double> p1, double lam){
	vector<double> q;
	for(int i=0; i<p1.size(); i++){
		q.push_back(p1[i]*lam);
	}
	return q;
}


double scalar(vector<double> ksi, vector<double> u){
	double sca = 0;
	for(int i=0; i<ksi.size(); i++){
		sca += ksi[i] * u[i];
	}
	return sca;
}
double W(double eps, pi p);
