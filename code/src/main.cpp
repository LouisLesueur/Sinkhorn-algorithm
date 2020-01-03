#include "simplex.h"
#include "pi.h"
#include <cmath>
#include <iostream>
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

	int n_iter = strtol(argv[1], nullptr, 0);;
	double eps = 2/square.length();
	double lambda = 0.8;
	simplex barycenter = bar(square, circle, lambda, eps, n_iter, "bary.png");
	barycenter.export_to_img();
	return 0;
}
