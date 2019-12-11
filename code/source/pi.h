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

		void plot();

		//Pour accèder aux marginales
		simplex first_marginal();
		simplex second_marginal();

		~pi(){delete [] tab;}
};

pi W(const simplex &s1, const simplex &s2, double eps, int n_iter);
simplex bar(const simplex & s1, const simplex & s2, double lambda, double eps, int n_iter);
