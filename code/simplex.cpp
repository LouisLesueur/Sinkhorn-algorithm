#include "simplex.h"


//======================CLASSE SIMPLEXE=======================

//---------------------CONSTRUCTEURS----------------------------
simplex::simplex(int n){
	size = n;
	tab = new double[n];
	double sum = 0;
	for(int i=0; i<n; i++){
		tab[i] = Random(0,100);
		sum += tab[i];
	}
	for(int i=0; i<n; i++)
		tab[i] /= sum;
}

simplex::simplex(double values[], int n){
	double sum = 0;
	size = n;
	tab = new double[n];
	for(int i=0; i<n; i++){
		sum += values[i];
		tab[i] = values[i];
	}
	std::cout<<sum<<std::endl;
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

//========================ALEA=======================
void InitRandom(){srand((unsigned int)time(0));}
double Random(double a, double b){return ( rand()/(double)RAND_MAX ) * (b-a) + a;}