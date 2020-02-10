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
