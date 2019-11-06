#pragma once

#include <iostream>
#include <fstream>
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
		char id;
	public:
		//constructeur vide, nécessaire pour pi
        simplex(){}
		//Pour générer un simplexe de taille n
		simplex(int n, char ID='0');
		//pour construire un simplexe à partir d'un tableau
		simplex(double values[], int n, char ID='0');
		//Pour récupérer la taille du simplexe
		int length() const{return size;};
		//To save in a csv file
		void plot();

		//Pour accèder à un élément du simplexe
		double operator()(int i) const;
		double& operator()(int i);

        ~simplex(){std::cout<<id<<std::endl; delete [] tab;}

};

//Pour afficher les éléments du simplexe
std::ostream& operator<<(std::ostream& str, const simplex& s);
