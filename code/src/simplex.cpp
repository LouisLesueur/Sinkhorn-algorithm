#include "simplex.h"


//======================SIMPLEX CLASS=======================

//---------------------CONSTRUCTORS----------------------------
simplex::simplex(matrix<float> VAL, int WIDTH, int HEIGHT, int CSTE, string path){
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

	int size = max(width, height);

	tab = ones_matrix<float>(size*size,1);
	int sum = 0;
	for(int i=0; i<width; i++){
		for(int j=0; j<height; j++){
			tab(i+width*j,0) = img(i,j);
			sum += img(i,j);
		}
	}
	constant = sum;

	tab=tab/sum;

	name=path;
}

//--------------------OPERATORS----------------------------
float simplex::operator()(int i) const{ return tab(i,0); }
float& simplex::operator()(int i){ return tab(i,0); }



void simplex::export_to_img(){
	matrix<uint8> img(height, width);
	for(int i=0; i<width; i++){
		for(int j=0; j<height; j++){
			float imag=tab(i+width*j,0)*constant;
			if(imag>255)
			  img(i,j)=255;
			else
			  img(i,j)=imag;
		}
	}

	save_png(img, name);
}
