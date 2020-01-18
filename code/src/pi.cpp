#include "pi.h"

using namespace dlib;
using namespace std;

//====================================================================


matrix<float> gen_K(int m, float eps){
  matrix<float> K(m,m);

  // Initialising
  #pragma omp parallel for
  for(int i=0; i<m; i++){
      for(int j=0; j<m; j++){
          float x_i = ((float)i+1/2)/m;
          float x_j = ((float)j+1/2)/m;
          float dis = pow(x_i - x_j, 2);
          K(i,j) = -dis/eps;
      }
  }
  return exp(K);
}

simplex bar(const matrix<float> &K, const matrix<float> &tK, const simplex & p1, const simplex & p2, float lambda, float eps, int n_iter, string name){
    float lamb1=lambda, lamb2=1-lambda;
    int m = p1.length();

		matrix<float> u1, u2, v1, v2, p;
    #pragma omp parallel sections
    {

      #pragma omp section
      {
      u1 = ones_matrix<float>(m,1);
    }
      #pragma omp section
      {
      u2 = ones_matrix<float>(m,1);
    }
      #pragma omp section
      {
      v1 = ones_matrix<float>(m,1);
    }
      #pragma omp section
      {
      v2 = ones_matrix<float>(m,1);
    }
  }

		#pragma omp parallel for
    for(int l=0; l<n_iter; l++){
        // Computing vk(n+1)
        v1 = pointwise_divide(p1.val(), tK*u1);
        v2 = pointwise_divide(p2.val(), tK*u2);
        // Computing p(n+1)
        p = pointwise_multiply(pow(pointwise_multiply(u1, K*v1), lamb1), pow(pointwise_multiply(u2, K*v2), lamb2));
        // Computing uk(n+1)
        u1 = pointwise_divide(p, K*v1);
        u2 = pointwise_divide(p, K*v2);
    }



    simplex out((lamb1*p1.cte() + lamb2*p2.cte())*p, p1.w(), p1.h(), 1, name);
    return out;
}
