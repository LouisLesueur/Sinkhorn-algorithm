#pragma once

#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <cassert>

//À appeler une fois dans le main
void InitRandom();
//Pour générer un double aléatoire entre a et b
double Random(double a, double b);

class simplex{
	private:
		double* tab;
		int size;
	public:
		//Pour générer un simplexe de taille n
		simplex(int n);
		//Pour récupérer la taille du simplexe
		int length() const{return size;};

		//Pour accèder à un élément du simplexe
		double operator()(int i) const;
		double& operator()(int i);

		~simplex(){delete [] tab;}

};

//Pour afficher les éléments du simplexe
std::ostream& operator<<(std::ostream& str, const simplex& s);
