#include "simplex.h"


//======================CLASSE SIMPLEXE=======================

//---------------------CONSTRUCTEURS----------------------------
simplex::simplex(matrix<double> VAL, int WIDTH, int HEIGHT, int CSTE, string path){
	name=path;
	constant=CSTE;
	width=WIDTH;
	height=HEIGHT;
	tab = VAL;
}

simplex::simplex(string path){
	matrix<uint8> img;
	load_image(img, path);
	
	width = img.nc();
	height = img.nr();
	
	tab = matrix<double>(width*height, 1);
	int sum = 0;
	for(int i=0; i<width; i++){
		for(int j=0; j<height; j++){
			tab(i+width*j,0) = img(i,j);
			sum += img(i,j);
		}
	}
	constant = sum;
	for(int i=0; i<width*height; i++)
		tab(i,0)/=sum;
	name=path;
}

//--------------------OPERATEURS----------------------------
double simplex::operator()(int i) const{ return tab(i,0); }
double& simplex::operator()(int i){ return tab(i,0); }



void simplex::export_to_img(){
	matrix<uint8> img(height, width);
	for(int i=0; i<width; i++){
		for(int j=0; j<height; j++)
			img(i,j)=tab(i+width*j,0)*constant;
	}

	save_png(img, name);
}
