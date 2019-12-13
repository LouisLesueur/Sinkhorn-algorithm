#pragma once

#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include "../CImg/CImg.h"

using namespace cimg_library;

//À appeler une fois dans le main
void InitRandom();
//Pour générer un double aléatoire entre a et b
double Random(double a, double b);

class simplex{
	private:
		double* tab;
		int size;
		char id;
		int constant, width, height;
	public:
		//constructeur vide, nécessaire pour pi
        simplex(){}
		//Pour générer un simplexe de taille n
		simplex(int n, char ID='0');
		//pour construire un simplexe à partir d'un tableau
		simplex(double values[], int n, int WIDTH, int HEIGHT, char ID='0');
		//simplex à partir d'une image
		simplex(const char* const path, char ID='0');
		//Pour récupérer la taille du simplexe
		int length() const{return size;};
		//To save in a csv file
		void plot();

		//To export as an image
		void export_to_img();

		//Pour accèder à un élément du simplexe
		double operator()(int i) const;
		double& operator()(int i);

        ~simplex(){delete [] tab;}

};

//Pour afficher les éléments du simplexe
std::ostream& operator<<(std::ostream& str, const simplex& s);
