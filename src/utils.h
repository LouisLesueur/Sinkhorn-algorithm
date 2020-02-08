#pragma once
#include "simplex.h"
#include <dlib/matrix.h>

//To generate cost matrix
matrix<float> gen_K(int m, float eps);
matrix<float> gen_C(int m);


//Some functions for log_bar
matrix<float> ave(float tau, const matrix<float> &u, const matrix<float> &u1);
matrix<float> lse(const matrix<float> &M);
matrix<float> M(const matrix<float> &H, const matrix<float> &u, const matrix<float> &v, const matrix<float> &C, float eps, int m);
