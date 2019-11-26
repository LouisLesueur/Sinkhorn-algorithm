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
	const int N = 200; // Intervale [0,1] divisé en N
	double values1[N];
	double values2[N];
	for(int i=0; i<N; i++){
		values1[i] = 1;
		values2[i] = 2* ((float)i) /N;
	}
	simplex s1(values1, N, '1');
	simplex s2(values2, N, '2');
	int n_iter = strtol(argv[1], nullptr, 0);;
	cout << "n_iter=" << n_iter << endl;
	//double eps = (double) strtol(argv[2], nullptr, 0);;
	double eps = 1/((double)3*n_iter);
	cout << "eps=" << eps << endl;
      	pi gamma = W(s1, s2, eps, n_iter);
	s1.plot();
	s2.plot();
	gamma.plot();
	return 0;
}
