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
          float x_i = ((float)i)/m;
          float x_j = ((float)j)/m;
          float dis = pow(x_i - x_j, 2);
          K(i,j) = -dis/eps;
      }
  }
  return make_symmetric(exp(K));
}

matrix<float> gen_C(int m){
  matrix<float> C(m,m);
  // Initialising
  #pragma omp parallel for
  for(int i=0; i<m; i++){
      for(int j=0; j<m; j++){
          float x_i = ((float)i+1/2)/m;
          float x_j = ((float)j+1/2)/m;
          float dis = pow(x_i - x_j, 2);
          C(i,j) = dis;
      }
  }
  return make_symmetric(C);
}

matrix<float> lse(matrix<float> M){
    return log(sum_cols(exp(M)));
}

matrix<float> M(matrix<float> u, matrix<float> v, matrix<float> C, float eps, int m){
    matrix<float> H = ones_matrix<float>(m,1);
    return (-C + u*trans(H) + H*trans(v))/eps;
}

simplex bar_log(const matrix<float> &C, const simplex & p1, const simplex & p2, float lambda, float eps, int n_iter, string name){
    float lamb1=lambda, lamb2=1-lambda;
    int m = p1.length();

        matrix<float> u1, u2, v1, v2, Lp;
    #pragma omp parallel sections
    {

      #pragma omp section
      {
      u1 = zeros_matrix<float>(m,1);
    }
      #pragma omp section
      {
      u2 = zeros_matrix<float>(m,1);
    }
      #pragma omp section
      {
      v1 = zeros_matrix<float>(m,1);
    }
      #pragma omp section
      {
      v2 = zeros_matrix<float>(m,1);
    }
  }
    //#pragma omp parallel for
    for(int l=0; l<n_iter; l++){
        // Computing uk(n+1)
	cout << "commence LSE" << endl;
        matrix<float> LSE_1 = lse(M(u1, v1, C, eps, m));
        matrix<float> LSE_2 = lse(M(u2, v2, C, eps, m));
	cout << "commence u" << endl;
        u1 = eps*log(p1.val()) - eps*LSE_1 + u1;
        u2 = eps*log(p2.val()) - eps*LSE_2 + u2;
        // Computing Lp(n+1)
	cout << "commence Lp" << endl;
        Lp = lamb1 * LSE_1 + lamb2 * LSE_2;
        // Computing vk(n+1)
	cout << "commence v" << endl;
        v1 = eps*Lp - eps*LSE_1 + v1;
        v2 = eps*Lp - eps*LSE_2 + v2;
	cout << exp(Lp) << endl;
    }

    matrix<float> p = exp(Lp);
    p = (lamb1*p1.cte() + lamb2*p2.cte()) * p;


    simplex out(p, p1.w(), p1.h(), 1, name);
    return out;
}

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



    simplex out((lamb1*p1.cte() + lamb2*p2.cte())*p, p1.w(), p1.h(), 1, name);
    return out;
}
