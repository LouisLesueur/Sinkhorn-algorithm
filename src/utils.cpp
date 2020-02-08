#include "utils.h"

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

matrix<float> ave(float tau, const matrix<float> &u, const matrix<float> &u1){ return tau*u + (1-tau)*u1; }
matrix<float> lse(const matrix<float> &M){ return log(sum_cols(exp(M))); }
matrix<float> M(const matrix<float> &H, const matrix<float> &u, const matrix<float> &v, const matrix<float> &C, float eps, int m){ return (-C + u*trans(H) + H*trans(v))/eps; }
