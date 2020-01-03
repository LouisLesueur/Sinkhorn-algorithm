#include "simplex.h"


//======================CLASSE SIMPLEXE=======================

//---------------------CONSTRUCTEURS----------------------------
simplex::simplex(int n, char ID){
	id = ID;
	size = n;
	width = n;
	height = 1;
	tab = new double[n];
	/*double sum = 0;
	for(int i=0; i<n; i++){
		tab[i] = Random(0,100);
		sum += tab[i];
	}
	for(int i=0; i<n; i++)
		tab[i] /= sum;
	constant = sum;*/
}

simplex::simplex(double values[], int n, int WIDTH, int HEIGHT, char ID){
	id = ID;
	double sum = 0;
	for(int i=0; i<n; i++){
		sum += values[i];
	}
	constant = sum;
	size = n;
	width = WIDTH;
	height = HEIGHT;
	tab = new double[n];
	for(int i=0; i<n; i++){
		tab[i] = values[i]/sum;
	}
}

simplex::simplex(const char* const path, char ID){
	array2d<double> img;
	load_image(img, path);
	
	width = img.nc();
	height = img.nr();
	size = width*height;
	tab = new double[size];
	int sum = 0;
	for(int i=0; i<width; i++){
		for(int j=0; j<height; j++){
			tab[i+width*j] = img[i][j];
			sum += img[i][j];
		}
	}
	constant = sum;
	for(int i=0; i<size; i++)
		tab[i]/=sum;
	id=ID;
}

//--------------------OPERATEURS----------------------------
double simplex::operator()(int i) const{
	assert(i<size);
	return tab[i];
}

double& simplex::operator()(int i){
	assert(i<size);
	return tab[i];
}

std::ostream& operator<<(std::ostream& str, const simplex& s){
	str<<"[ ";
	for(int i=0; i<s.length(); i++)
		str << s(i) << " ";
	str<<"]";
	return str;
}

//======================================================

void simplex::plot(){
	std:: string name = "simplex";
	name += id;
	name += ".csv";
	std::ofstream monFlux(name.c_str());
	for(int i=0; i<size; i++)
		monFlux<<i<<" "<<tab[i]<<std::endl;
}


void simplex::export_to_img(){
	array2d<double> img(width, height);
	for(int i=0; i<width; i++){
		for(int j=0; j<height; j++)
			img[i][j]=tab[i+width*j]*constant;
	}
	std::string name = "img";
	name += id;
	name += ".png";
	const char* path = name.c_str();
	save_png(img, path);
}


//========================ALEA=======================
//void InitRandom(){srand((unsigned int)time(0));}
//double Random(double a, double b){return ( rand()/(double)RAND_MAX ) * (b-a) + a;}
