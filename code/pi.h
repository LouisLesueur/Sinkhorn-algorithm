#pragma once
#include "simplex.h"

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


//À coder !
//On considère que le xi du papier = 1
double W(double eps, pi p);
