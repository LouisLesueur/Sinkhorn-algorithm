#include "pi.h"

using namespace dlib;
using namespace std;

//====================================================================

simplex bar(const simplex & p1, const simplex & p2, double lambda, double eps, int n_iter, string name){
	double lamb1=lambda, lamb2=1-lambda;
	int m = p1.length();

	matrix<double> K(m,m);
	set_all_elements(K,0);

	// Initialising
	for(int i=0; i<m; i++){
		for(int j=0; j<m; j++){
			double x_i = ((double)i+1/2)/m;
			double x_j = ((double)j+1/2)/m;
			double dis = pow(x_i - x_j, 2);
			K(i,j) = exp(-dis/eps);
		}
	}

	matrix<double> a1, a2, p, b1, b2;
	a1 = ones_matrix<double>(m,1);
	a2 = ones_matrix<double>(m,1);

	for(int l=0; l<n_iter; l++){
		p = pointwise_multiply(pow(trans(K)*a1, lamb1), pow(trans(K)*a2, lamb2));
		b1 = pointwise_divide(p, trans(K)*a1);
		b2 = pointwise_divide(p, trans(K)*a2);
		a1 = pointwise_divide(p1.val(), K*b1);
		a2 = pointwise_divide(p2.val(), K*b2);
	}

	simplex out((lamb1*p1.cte() + lamb2*p2.cte())*p, p1.w(), p1.h(), 1, name);
	return out;
}

simplex bar2(const simplex & p1, const simplex & p2, double lambda, double eps, int n_iter, string name){
    double lamb1=lambda, lamb2=1-lambda;
    int m = p1.length();

    matrix<double> K(m,m);
    set_all_elements(K,0);

    // Initialising
    for(int i=0; i<m; i++){
        for(int j=0; j<m; j++){
            double x_i = ((double)i+1/2)/m;
            double x_j = ((double)j+1/2)/m;
            double dis = pow(x_i - x_j, 2);
            K(i,j) = exp(-dis/eps);
        }
    }

		matrix<double> u1, u2, v1, v2, p;
    u1 = ones_matrix<double>(m,1);
    u2 = ones_matrix<double>(m,1);
    v1 = ones_matrix<double>(m,1);
    v2 = ones_matrix<double>(m,1);

    for(int l=0; l<n_iter; l++){
        // Computing vk(n+1)
        v1 = pointwise_divide(p1.val(), trans(K)*u1);
        v2 = pointwise_divide(p2.val(), trans(K)*u2);
        // Computing p(n+1)
        p = pointwise_multiply(pow(pointwise_multiply(u1, K*v1), lamb1), pow(pointwise_multiply(u2, K*v2), lamb2));
        // Computing uk(n+1)
        u1 = pointwise_divide(p, K*v1);
        u2 = pointwise_divide(p, K*v2);
    }

    simplex out((lamb1*p1.cte() + lamb2*p2.cte())*p, p1.w(), p1.h(), 1, name);
    return out;
}
