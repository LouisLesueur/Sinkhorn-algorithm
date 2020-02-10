#include "barycenter.h"

using namespace dlib;
using namespace std;

//====================================================================
simplex bar(const matrix<float> &K, const simplex & p1, const simplex & p2, float lambda, float eps, int n_iter, string name){
    float lamb1=lambda, lamb2=1-lambda;
    int m = p1.length();

		matrix<float> u1, u2, v1, v2, p;
    u1 = ones_matrix<float>(m,1);
    u2 = ones_matrix<float>(m,1);
    v1 = ones_matrix<float>(m,1);
    v2 = ones_matrix<float>(m,1);


  int count = 0;

		#pragma omp parallel for
    for(int l=0; l<n_iter; l++){
        // Computing vk(n+1)
        v1 = pointwise_divide(p1.val(), K*u1);
        v2 = pointwise_divide(p2.val(), K*u2);
        // Computing p(n+1)
        p = pointwise_multiply(pow(pointwise_multiply(u1, K*v1), lamb1), pow(pointwise_multiply(u2, K*v2), lamb2));
        // Computing uk(n+1)
        u1 = pointwise_divide(p, K*v1);
        u2 = pointwise_divide(p, K*v2);
        count++;
        if(count<n_iter){
            cout << int((float)count/(float)n_iter * 100.0) << "% \r";
            cout.flush();
        }

    }

    simplex out(p, p1.w(), p1.h(), 255/max(p), name);
    return out;
}
