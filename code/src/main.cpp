#include "simplex.h"
#include "pi.h"
#include <cmath>
#include <iostream>
#include <cstdlib>
using namespace std;

double gaussian(double x, double sigma, double mu){
	return 1/(sigma*sqrt(2*M_PI)) * exp(-pow(x-mu,2)/(2*pow(sigma,2)));
}


int main(int argc, char *argv[])
{	
	simplex square("../square.png");
	simplex circle("../circle.png");

	simplex test(square.val(), square.w(), square.h(), square.cte(), "lol.png");
	test.export_to_img();

	int n_iter = strtol(argv[1], nullptr, 0);
	double fact = strtold(argv[2], nullptr);
	double lambda = strtold(argv[3], nullptr);

	double eps=1/double(fact*square.length());
	simplex barycenter = bar(square, circle, lambda, eps, n_iter, "bary.png");
	barycenter.export_to_img();
	return 0;
}
