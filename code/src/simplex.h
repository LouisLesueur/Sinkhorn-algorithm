#pragma once

#include <iostream>
#include <fstream>
#include <ctime>
#include <cstring>
#include <cstdlib>
#include <dlib/gui_widgets.h>
#include <dlib/pixel.h>
#include <dlib/matrix.h>
#include <dlib/array2d.h>
#include <dlib/image_io.h>

using namespace dlib;
using namespace std;

//À appeler une fois dans le main
//void InitRandom();
//Pour générer un double aléatoire entre a et b
//double Random(double a, double b);

class simplex{
	private:
		matrix<double> tab;
		string name;
		int constant, width, height;
	public:
		//constructeur vide
		simplex(){}
		//constructeur de simplexe 'vide'
		simplex(matrix<double> VAL, int WIDTH, int HEIGHT, int CSTE, string path);
		//pour construire un simplexe à partir d'une image
		simplex(string path);
		
		//Pour récupérer la taille du simplexe
		int length() const{return width*height;};
		int w() const{return width;};
		int h() const{return height;};
		int cte() const{return constant;};
		matrix<double> val() const{return tab;};

		//To export as an image
		void export_to_img();

		//Pour accèder à un élément du simplexe
		double operator()(int i) const;
		double& operator()(int i);

};

//Pour afficher les éléments du simplexe
std::ostream& operator<<(std::ostream& str, const simplex& s);
