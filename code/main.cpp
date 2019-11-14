#include "simplex.h"
#include "pi.h"
#include "matrice.h"
#include <cmath>
#include <iostream>
using namespace std;

double gaussian(double x, double sigma, double mu){
	return 1/(sigma*sqrt(2*M_PI)) * exp(-pow(x-mu,2)/(2*pow(sigma,2)));
}


int main(int argc, char *argv[])
{	
	InitRandom();
	const int N = 100; // Intervale [0,1] divisé en 100
	double values[N];
	for(int i=0; i<N; i++)
		values[i] = i;
	simplex s1(values, N, '1');
	simplex s2(values, N, '2');
	cout << s1 << endl;
	int n_iter = strtol(argv[1], nullptr, 0);;
	cout << "n_iter=" << n_iter << endl;
	double eps = strtol(argv[2], nullptr, 0);; // Pour des valeurs de argv[2]<1, eps=0 !
	cout << "eps=" << eps << endl;
	eps = 0.1;
      	pi gamma = W(s1, s2, eps, n_iter);

	bool show_marginals = false;
	if(show_marginals){
		cout << "Premier simplexe   = " << s1 << endl;
		cout << "Première marginale =" << gamma.first_marginal() << endl;
		cout << "Second simplexe    = " << s2 << endl;
		cout << "Seconde marginale  =" << gamma.second_marginal() << endl;
	}
	bool show_matrix = false;
	if(show_matrix){
		Matrice M(N);
		for(int i=0; i<N; i++){
			for(int j=0; j<N; j++){
				cout << i << " " << j << " " << gamma(i,j) << endl;
				M(i,j) = gamma(i,j);
			}
		}
		cout << "Matrice minimisante = " << M << endl;
	}
	s1.plot();
	s2.plot();
	gamma.plot();
	return 0;
}
