#include "simplex.h"
#include "pi.h"
#include <cmath>
#include <iostream>
using namespace std;

double gaussian(double x, double sigma, double mu){
	return 1/(sigma*sqrt(2*M_PI)) * exp(-pow(x-mu,2)/(2*pow(sigma,2)));
}


int main()
{	
	InitRandom();	
	int range = 10; // On prend des valeurs de gaussiennes pour des x allant de -range à range
	double values1[2*range+1];
	double values2[2*range+1];
	double mu = 5;
	double sigma = 1;
	for(int i=0; i<2*range+1; i++){
		double x = (double) i-range;
		values1[i] = gaussian(x,sigma,mu);
		values2[i] = gaussian(x,sigma,-mu);
	}
	simplex s1(values1, 2*range + 1);
	simplex s2(values2, 2*range + 1);
	int n_iter = 100;
	double eps = 30/n_iter;
	pi gamma = W(s1, s2, eps, n_iter);

	bool show_marginals = true;
	if(show_marginals){
		cout << "Premier simplexe   = " << s1 << endl;
		cout << "Première marginale =" << gamma.first_marginal() << endl;
		cout << "Second simplexe    = " << s2 << endl;
		cout << "Seconde marginale  =" << gamma.second_marginal() << endl;
	}
	gamma.plot();

	return 0;
}
