#include "pi.h"

using namespace dlib;
using namespace std;

//====================================================================

matrix<float> gen_K(int m, float eps){
  matrix<float> K(m,m);
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


matrix<float> M(const matrix<float> &H, matrix<float> u, matrix<float> v, matrix<float> C, float eps, int m){
    return (-C + u*trans(H) + H*trans(v))/eps;
}

simplex bar_log(const matrix<float> &C, const simplex & p1, const simplex & p2, float lambda, float eps, int n_iter, string name){

    float lamb1=lambda, lamb2=1-lambda;
    int m = p1.length();
    const matrix<float> H = ones_matrix<float>(m,1);

    matrix<float> u1, u2, v1, v2, Lp, LSE_v1, LSE_v2;
    u1 = zeros_matrix<float>(m,1);
    u2 = zeros_matrix<float>(m,1);
    v1 = zeros_matrix<float>(m,1);
    v2 = zeros_matrix<float>(m,1);

    cout << 0 << "% \r";
    cout.flush();

    for(int l=0; l<n_iter; l++){


        // Computing uk(n+1)
        u1 += eps*(log(p1.val()) - lse(M(H, u1, v1, C, eps, m)));
        u2 += eps*(log(p2.val()) - lse(M(H, u2, v2, C, eps, m)));

        cout << int(0.25*((float)(l+1)/(float)n_iter) * 100.0) << "% \r (u updated)";
        cout.flush();

        LSE_v1 = lse(trans(M(H, u1, v1, C, eps, m)));
        LSE_v2 = lse(trans(M(H, u2, v2, C, eps, m)));

        cout << int(0.5*((float)l/(float)n_iter) * 100.0) << "% \r (LSEv computed)";
        cout.flush();

        Lp = lamb1 * LSE_v1 + lamb2 * LSE_v2 ;

        cout << int(0.75*((float)(l+1)/(float)n_iter) * 100.0) << "% \r (Lp computed)";
        cout.flush();

        // Computing vk(n+1)
        v1 += eps*(Lp - LSE_v1);
        v2 += eps*(Lp - LSE_v2);

        cout << int((float)(l+1)/(float)n_iter * 100.0) << "% \r (v updated)";
        cout.flush();

    }

    matrix<float> p = (lamb1*p1.cte() + lamb2*p2.cte()) * exp(Lp);

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
